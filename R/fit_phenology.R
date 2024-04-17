# -----------------------------------------------------------------------------
# phenology model: code to create nimbleFunction's describing the mean function
# of the Omeyer and Girondot phenology models
# -----------------------------------------------------------------------------

#' Returns a numbleFunction implementing the Omeyer et al. 2022 model of mean number of counts
#' as a function of days. By default this fits a 3-parameter symmetrical curve with no peak flattening
#' @export

phenology.omeyer <- nimble::nimbleFunction(
  run = function(t = double(0), alpha = double(0), s1 = double(0),
                 s2 = double(0,default=0), tf = double(0,default=0), tp = double(0)) {

    if(s2 == 0) s2 = s1 
    if(t < (tp - tf)) {
      mu <- alpha * exp(-((t - tp + tf) / s1)^2)
    } else {
      if(t < (tp + tf)){
        mu <- alpha
      } else {
        mu <- alpha * exp(-((t - tp - tf) / s2)^2)
      }
    }
    return(mu)
    returnType(double(0))
  }
)

#' Returns a numbleFunction implementing the phenology model of
#' Girondot et al. 2010 Estimating density of animals during migratory
#' waves: a new model applied to marine turtles at nesting sites.
#' @export

phenology.girondot <- nimble::nimbleFunction(
  run = function(t = double(0), B = double(0), E = double(0,default=0), Max = double(0),
                 MinE = double(0,default=0), MinB = double(0,default=0),tf = double(0,default=0), tp = double(0)) {
    
    if(E == 0) E = tp + (tp - B) 
    pi = 3.141593
    if(t < B){ mu <- MinB} else {
      if(t < (tp - tf)) {
        mu <- ((1+cos(pi*(tp - tf - t)/(tp - tf - B)))/2)*(Max - MinB)+MinB
      } else {
        if(t < (tp + tf)){
          mu <- Max
        } else {
          if(t < E){
            mu <- ((1+cos(pi*(t - (tp + tf))/(E - (tp + tf)))))/2*(Max - MinE)+MinE
          } else {
            mu <- MinE
          }
        }
      }
    }
    return(mu)
    returnType(double(0))
  }
)

# tibble(t = 1:365) %>%
# mutate(N = map_dbl(t,~phenology.girondot(.x, B = 20, tp = 110, Max = 400)))

# -----------------------------------------------------------------------
# MT_NimModel: Basic code for the NIMBLE model
# -----------------------------------------------------------------------

# By fixing some parameters to zero in this code (e.g. tf and s2) we could fit a range 
# of different possible shapes. In the priors section rather than stochastic, we make
# some determininstic nodes with values = 0, depending on which parameters are specified in the run
# function. 
# Generates `nimbleCode` for fitting the Omeyer et al. 2022 model in `Nimble`

MT_NimModel <- nimble::nimbleCode({
  ## likelihood

  for(r in 1:R){
    for(i in 1:N) {
      for(j in 1:nt[i]) {
        mu[i, j, r] <- phenology.omeyer(t[i, j], alpha[r,bch[i]], s1[r,bch[i]], s2[r,bch[i]], tf[r,bch[i]], tp[r,bch[i]])
      }
      Y[i,r] ~ dSnbinomNim(phi[r,bch[i]], mu[i, ,r], nt[i])
    }
  }
  ## priors

  for(r in 1:R){

    for(b in 1:B){
      alpha_mn[r,b] ~ dexp(1)
      alpha_rate[r,b] <- 1 / alpha_mn[r,b]
      alpha[r,b] ~ dexp(alpha_rate[r,b])

      tf[r,b] ~ dexp(0.2)
      phi[r,b] ~ T(dinvgamma(shape = 1, rate = 0.1), , 50)
      
      s_rate[r,b] ~ dunif(0.01, 10)
      s1[r,b] ~ dexp(s_rate[r,b])
      s2[r,b] ~ dexp(s_rate[r,b])

      #stp[r,b] ~ dunif(0, 10)
      #tp_ast[r,b] ~ dnorm(0, sd = 1)
      #tp[r,b] <- Pk + tp_ast[r,b] * stp[r,b]
      tp[r,b] ~ dunif(0,365)
    }

  }

})

# ------------------------------------------------------------------------------------------------------------------
# New code below to adapt the base model depending on different model parameterisations. Here
# model is either 'phenology.omeyer' or 'phenology.girondot'. 
# -------------------------------------------------------------------------------------------------------------------

# Priors are specified as a hierarchical named list which allows users to tweak elements

 phenology_priors = function(model,params=NULL){
         
          if(is.null(params)) params = formalArgs(model) %>% subset(.!='t')
          list(tf = 'tf ~ dexp(0.2)', 
             phi = 'phi ~ T(dinvgamma(shape = 1, rate = 0.1), , 50)',
              alpha = list(alpha_mn = 'alpha_mn ~ dexp(1)',
                           alpha = 'alpha ~ dexp(1/alpha_mn)'),
              tp = 'tp ~ dunif(0,365)',
              s1 = list(s_rate = 's_rate ~ dunif(0.01, 10)',
                        s1 = 's1 ~ dexp(s_rate)'),
             s2 = 's2 ~ dexp(s_rate)'
            )[params]           
 }
   

# Now a function to build the blocks

phenology_model = function(model,params = NULL,priors=NULL){

 if(is.null(params)) params = formalArgs(model) %>% subset(.!='t')
 if(is.null(priors)) priors = phenology_priors(model,params)

 lik = paste(paste0(model,'(t[i, j],'), paste(paste0(params,'=',params,'[r,bch[i]]'), collapse = ', '),')')
 lik  = c('for(r in 1:R){', 
            'for(i in 1:N) {', 
              'for(j in 1:nt[i]) {',
                 paste('mu[i, j, r] <-', lik),
              '}',
            'Y[i,r] ~ dSnbinomNim(phi[r,bch[i]], mu[i, ,r], nt[i])',
          '}',
          '}')
                
 priors = unlist(priors) 
 names = sub('.*\\.', '', names(priors))
 rep = paste0(names,'[r,b]') %>% setNames(paste0(names,'\\b'))
 priors = str_replace_all(priors,rep) %>%
          c('for(r in 1:R){','for(b in 1:B){',.,'}','}')
  
 c(lik,priors) %>%
 parse(text=.) %>%
 c(quote(`{`),.) %>%
 as.call()

}


# ----------------------------------------------------------------------------------
# MT_diff_matrix: Helper function for making a matrix containing the days observed in each sampling window
# ----------------------------------------------------------------------------------

#' @param t is a vector of survey days from an object created by `MT_prep`
#' @param w is a vector of sampling windows giving the number of days represented by each survey ub `t`

# ! This may need modifying in the case where w = 1 as single column matrices
# don't seem to be well handled by the custom distribution in NIMBLE

MT_diff_matrix = function(t,w){

  tdiffs =  purrr::map2(t, w, function(t,w) (t - w + 1):t)

  t <- matrix(NA, length(tdiffs), max(c(2,w)))     #Here it forces ncol to be minimum of 2 until issues can be resolved
  for(i in 1:length(tdiffs)) {
    t[i, 1:length(tdiffs[[i]])] <- tdiffs[[i]]
  }
  t
}


# ----------------------------------------------------------------------------------
# MT_make_model: Function for building the nimble model
# ----------------------------------------------------------------------------------

#' @param data is a dataframe created by MT_prep

#' @param response is a character vector giving the response

#' @param Pk is a single numeric giving an approximate peak day for counts which is used as a starting
#' point for model fits. Omeyer et al. 2022 found that using an approximate peak in priors resulted
#' in better covnvergence.

MT_make_model <- function(data, response = c('activities','nests'), null=FALSE){

  tmat = MT_diff_matrix(data$day,data$window)
  consts = list(t=tmat,nt = rowSums(!is.na(tmat)), N = nrow(tmat))
  consts$bch <- if(has_name(data,'beach')) as.numeric(data$beach) else rep(1,consts$N)
  consts$B = max(consts$bch)

  # Remove any response columns that are all NAs prior to fitting - unless fitting null model
  Y = as.matrix(dplyr::select(data,any_of(c('activities','nests'))))

  if(!null){
    Y = Y[, colSums(is.na(Y)) < nrow(Y), drop=F]
  } else Y[] <- NA

  consts$R = ncol(Y)
  consts$y.names = colnames(Y)

  nimble::nimbleModel(
    code = MT_NimModel,
    constants = consts,
    data = list(Y = Y),
    dimensions = list(mu = c(dim(consts$t),ncol(Y))),
    check = nimble::getNimbleOption("checkModel"),
    buildDerivs = nimble::getNimbleOption("buildModelDerivs")
  ) %>%
    nimble::compileNimble()

}


# ----------------------------------------------------------------------------------
# MT_initialize: Function for initializing the nimble model
# ----------------------------------------------------------------------------------

# Generates a valid set of initial parameters for the model by sampling from the prior
# distributions.

#' @param model is a model structure created by `MT_make_model`

#' @param smin,smax Single numerics giving the lower and upper bounds of a uniform distribution for
#' the s_rate parameter which are used for adjusting the default prior if needed. See Details.

#' @attemmpts Number of attempts to find a valid set of starting parameters before giving up.

#' In this function, if the prior on s_rate is kept at the default vague priors `runif(0.01,10)` then it can take
#' an impossibly long number of iterations to find a valid set of initial parameters, because s1 and s2
#' that are sufficiently high and close to real values are very unlikely to be sampled. You can make it fit by narrowing the
#' prior to something like `runif(0.01,0.2)`. Until some way can be found of quickly estimating a sensible
#' starting point for s_rate, options are provided for manually setting the limits of the uniform distribution. Some
#' trial and error may be needed to find values that initialize quickly.

MT_initialize <- function(model, tp.min, tp.max, smin, smax, attempts = 100) {

  valid <- 0
  consts = model$origData
  B = consts$B
  R = consts$R
  tdiffs = consts$t
  bch = consts$bch
  attempt = 1

  while(valid == 0 & attempt < attempts) {
    ## sample initial values from the priors

    alpha_mn <- matrix(rexp(B*R,1),nrow=R)
    alpha_rate <- 1 / alpha_mn
    alpha <- matrix(rexp(B*R,alpha_rate),nrow=R)

    s_rate <- matrix(runif(R*B, smin, smax),nrow=R)
    s2 <- s1 <- matrix(rexp(R*B,s_rate),nrow=R)

    tp <- matrix(runif(B*R, tp.min, tp.max),nrow=R)

    tf <- matrix(rexp(B*R,0.2),nrow=R)

    phi <- matrix(rep(51,R*B),nrow=R)
    while(any(phi > 50)){
      phi <-matrix(rexp(B*R,1),nrow=R)
    }

    mu = array(dim = c(dim(tdiffs),R))
    y.init <- model$Y

    for(r in 1:R){
      for(i in 1:nrow(tdiffs)) {
        for(j in 1:consts$nt[i]) {
          mu[i, j, r] <- phenology.omeyer(tdiffs[i, j], alpha[r,bch[i]], s1[r,bch[i]], s2[r,bch[i]], tf[r,bch[i]], tp[r,bch[i]])
        }
        y.init[i,r] <- rnbinom(1,phi,mu=mu[i,1,r])
      }
    }

    mu[is.na(mu)] <- 0
    y.init[!is.na(model$Y)] <- NA

    inits <- list(
      alpha_mn = alpha_mn,
      alpha_rate = alpha_rate,
      alpha = alpha,
      s_rate = s_rate,
      s1 = s1,
      s2 = s2,
      #stp = stp,
      #tp_ast = tp_ast,
      tp = tp,
      tf = tf,
      phi = phi,
      mu = mu,
      Y = y.init
    )
    ## check validity of initial values
    model$setInits(inits)
    valid <- ifelse(!is.finite(model$calculate()), 0, 1)
    attempt = attempt+1
  }
  if(valid == 0) stop(paste('Model failed to initialize after',attempts,'attempts. Try increasing attempts or adjusting smin and smax'))
  return(inits)
}

# initial values for peak could also be found automatically with
# loess(activities ~ day, data = test$data[[1]]) %>%
# predict(newdata = data.frame(day=1:365)) %>%
# which.max


# ----------------------------------------------------------------------------------------------------------------
# MT_run_model: Function for compiling and running the nimble model
# ----------------------------------------------------------------------------------------------------------------

#init.control = function() list(smin = 0.01,smax=10)

# ----------------------------------------------------------------------------------------------------------------
# MT_run_model: Function for compiling and running the nimble model
# ----------------------------------------------------------------------------------------------------------------

#' It should be possible to get this model to only fit certain parameters using the 'parameters.to.monitor' argument, with
#' other values fixed to zero (e.g. for flattening) or s1 = s2 for symmetrical.
#' Use fixed = setdiff(formalArgs(phenology.omeyer),params) %>% subset(.!='t') to find parameters to fix at 0 
                        
#' @param model A model created by MT_make_model

#' @param nchains,niter,nburnin,thin Single numerics giving the number of MCMC chains to run, the total number of iterations
#' per chain, the number of burnin samples to discard, and the sample thining rate, respectively

#' @param parameters.to.monitor Names of the model parameters to monitor. Defaults are generally ok.

#' @param init.attempts Number of attempts to make at finding a valid set of starting parameters for the model before
#' gracefully giving up.

MT_run_model = function(data,nchains=2, niter=40000, nburnin = 10000, thin = 1,
                        parameters.to.monitor = c("alpha", "s1", "s2", "tf", "tp", "phi"),
                        init.attempts = 100,
                        init.control = list(tp.min = 0, tp.max = 365, smin = 0.01,smax=1),null=FALSE) {


  assign('dSnbinomNim', dSnbinomNim, envir = .GlobalEnv)
  assign('rSnbinomNim', rSnbinomNim, envir = .GlobalEnv)

  nimble::nimbleOptions(verbose = FALSE)

  cat('#Compiling model\n')
  model = MT_make_model(data,null=null)
  if(null) data[TRUE,model$origData$y.names] <- NA

  cat('#Initialising model\n')
  #init.control = modifyList(init.control(),init.control)   # This line currently prevents MT_run_model in parallel - function init.control() can't be found)
  inits <- purrr::map(1:nchains, ~MT_initialize(model,attempts = init.attempts,
                                         smin = init.control$smin,
                                         smax=init.control$smax,
                                         tp.min = init.control$tp.min,
                                         tp.max = init.control$tp.max))

  cat('#Configuring model\n')
  config <- nimble::configureMCMC(model, monitors = parameters.to.monitor, thin = thin)
  built  <- nimble::buildMCMC(config)
  cbuilt <- nimble::compileNimble(built)

  ## run the model for multiple chains
  fit = nimble::runMCMC(cbuilt, nchains = nchains, niter = niter, nburnin = nburnin, inits = inits, thin = thin,
                progressBar = TRUE, samplesAsCodaMCMC = TRUE)

  class(fit) <- c('phenology_model',class(fit))
  attr(fit,'y.names') <- model$origData$y.names
  attr(fit,'data') <- dplyr::select(data,day,window,any_of(c('beach',model$origData$y.names)))
  attr(fit,'ref_date') <- attr(data,'ref_date')

  return(fit)

}


# ----------------------------------------------------------------------------------------------------------------
# MT_fit: Function for running a nimble model across multiple sites/seasons in a nested dataset produced by MT_prep
# ----------------------------------------------------------------------------------------------------------------

#' This function currently requires that you put some bounds on parameters for selecting valid initital values in order
#' to get it to initialise in any sensible amount of time. 
#'
#' @export

MT_fit <- function(data, nchains=2, niter=40000, nburnin = 10000, thin = 1,
                   params = c("alpha", "s1", "s2", "tf", "tp", "phi"),
                   init.attempts = 100,
                   init.control = list(tp.min = 0, tp.max = 365, smin = 0.01,smax=1),
                   ncores=1){

  if(is(data,'MTsim')) data <- MT_sim2df(data)
  if(!is(data,'phenology_df')) stop('data should be a phenology_df object created by as_phenology')

  if(ncores>1){

    cat(paste('#Fitting',nrow(data),'models in parallel'))
    future::plan(future::tweak(future::multisession,workers = ncores))
    data = mutate(data,fit = furrr::future_map(data,MT_run_model,nchains=nchains,niter=niter, nburnin = nburnin, thin = thin,
                                               parameters.to.monitor = params,init.attempts = init.attempts,init.control = init.control))
    future::plan(future::sequential)

  } else {

    data = dplyr::mutate(data,fit = purrr::map(data,MT_run_model,nchains=nchains,niter=niter, nburnin = nburnin, thin = thin,
                                        parameters.to.monitor = params,init.attempts = init.attempts,init.control = init.control))

  }

  return(data)

}
