
# ----------------------------------------------------------------------------------
# Function for sampling from the final model
# ----------------------------------------------------------------------------------

# Where this fails it has been because the default sep argument in
# spread draws "[, ]" doesn't work. Dissection of the source code
# has shown this to be the case and that replacing with sep = ', ' fixes it.

# Use pars argument with caution - not yet fully implemented!!!!!

MT_sample = function(model,n = 2000,pars){

  data = attr(model,'data')
  model %<>% tidybayes::recover_types(data)
  samples = tidybayes::spread_draws(model,alpha[Y,beach],s1[Y,beach],s2[Y,beach]
                         ,tf[Y,beach],tp[Y,beach],phi[Y,beach],ndraws=n,sep=', ') %>%
    mutate(.draw = 1:n) %>%
    dplyr::ungroup()
  samples$Y <- attr(model,'y.names')[samples$Y]

  if(!rlang::has_name(data,'beach')) {samples$beach <- NULL}

  samples = rename(samples,y.var=Y)

  # Subsetting of specific parameters could be done much more efficiently in spread draws above
  # but this works as quick fix. It doesn't currently respect different beaches or y.vars
  if(!missing(pars)) samples = samples[,pars,drop=TRUE]

  return(samples)

}

# ----------------------------------------------------------------------------------
# Function for predicting mean and values from fitted model
# ----------------------------------------------------------------------------------

# This currently matches column positions to argument names in meanFnNimm, but it won't work
# if meanFnNim has a different number of arguments e.g. a different functional form as
# we expand the code later.... we could this do this by collapsing the args to a list in a do.call maybe
# This is the better way of doing this using column-argument matching:
# MT_sample(model,samples) %>%
# mutate(.iter = row_number()) %>%
# expand_grid(t = days) %>%
# mutate(total = pmap_dbl(list(!!!rlang::parse_exprs(formalArgs('meanFnNim'))), meanFnNim))
# Will also enable different mean functions (models) to be passed.
#' @export

predict.phenology_model = function(model,days, CI=TRUE, full.posterior=FALSE, draws=2000){

  data = attr(model,'data')

  if(missing(days)) days = min(data$day):max(data$day)


 if(!full.posterior & !CI) {
    samples = summary(model) %>%
              dplyr::select(y.var:mean) %>%
              tidyr::spread(key=variable,value=mean)
    } else {
    samples = MT_sample(model,draws)
   }

  # Generate predictions from model
  samples =
    tidyr::expand_grid(samples,t = days) %>%
    mutate(mu = purrr::pmap_dbl(list(!!!rlang::parse_exprs(formalArgs('phenology.omeyer'))), phenology.omeyer)) %>%
    dplyr::rename(day = t) %>%
    mutate(N.pred = rnbinom(day,phi,mu=mu)) %>%
    append_observed(data) %>%
    mutate(date = attr(model,'ref_date') + day) %>%
    dplyr::select(any_of('.draw'),y.var,any_of('beach'),date,day,mu:N)

  # If no full posterior don't include daily simulated values in output
  if(!full.posterior) samples = dplyr::select(samples,-N.pred,-N)

  # Calculate CIs if required
  if(!full.posterior & CI){
    samples = daily_means(samples,by_site=T)
  }

  class(samples) <- c('MTpred',class(samples))
  return(samples)
}


#' @export

predict.phenology_df = function(object,days,ncores=1,...){

  if(!all(c('data','fit') %in% names(object))) stop ("object should be of class phenology_df with columns 'data' and 'fit'")

  if(missing(days)){
    days = bind_rows(object$data)$day
    days = min(days):max(days)
  }


  if(ncores>1){

    future::plan(future::tweak(multisession,workers = ncores))
    object = mutate(object,predict = furrr::future_map(fit,predict.phenology_model,days=days,...))
    plan(sequential)

  } else {

    object = mutate(object,predict = purrr::map(fit,predict.phenology_model,days=days,...))

  }

  return(object)

}


# ----------------------------------------------------------------------------------
# MT_append_observed: Replaces predictions with observed counts, where available
# ----------------------------------------------------------------------------------

append_observed = function(preds,data){

  t = MT_diff_matrix(data$day,data$window)

  obs.days = dplyr::select(data,any_of(c('beach','day',unique(preds$y.var)))) %>%
  tidyr::gather('y.var','N.obs',-any_of(c('beach','day')))

  obs.days =
    data.frame(t) %>%
    {if(has_name(data,'beach')) bind_cols(.,beach=data$beach) else .} %>%
    tidyr::gather('.col','day',-any_of('beach')) %>%
    dplyr::select(-.col) %>%
    tidyr::drop_na() %>%
    left_join(obs.days,by='day') %>%
    tidyr::replace_na(list(N.obs = 0))

  preds %>%
    left_join(obs.days,by=dplyr::intersect(c('day','beach','y.var'),names(obs.days))) %>%
    mutate(N = dplyr::coalesce(N.obs,N.pred))

}

# ---------------------------------------------------------------------------------------------------------------------------------------
# plot.MTpred: plot method for MTpred
# ---------------------------------------------------------------------------------------------------------------------------------------

plot.MTpred <- function(obj,by_site=T,colour = c(nests='orange',activities='blue'),shape=21){

  orig = na.omit(distinct(dplyr::select(obj,y.var,any_of('beach'),day,N.obs)))
  cat('Calculating daily means and CIs')
  mean.line = daily_means(obj,by_site=by_site)

  pl =
    ggplot(mean.line,aes(x = day,y = mu)) +
    geom_ribbon(aes(ymin = .lower,ymax = .upper,group = y.var),alpha=.3) +
    geom_line(aes(colour=y.var)) +
    scale_colour_manual(values = colour) +
    ylab('Count')

  if(by_site & length(orig)) pl = pl + geom_point(data=orig,aes(y = N.obs,colour=y.var),shape=shape)

  if(has_name(mean.line,'beach')) pl + facet_wrap(~beach) else pl

}




# ---------------------------------------------------------------------------------------------------------------------------------------
# summary.MTfit: summary method for MTfit
# ---------------------------------------------------------------------------------------------------------------------------------------

#' @param model A model fit using MT_fit

#' @returns A tibble containing the mean, median, sd, credible intervals, potential scale reduction factor and effective sample size for
#' parameters.
#' @export

summary.phenology_model = function(model){

  .s = suppressWarnings(tidybayes::spread_draws(model,alpha[Y,beach],
                                     phi[Y,beach],
                                     s1[Y,beach],
                                     s2[Y,beach],
                                     tf[Y,beach],
                                     tp[Y,beach])) %>%
    tidybayes::summarise_draws() %>%
    dplyr::select(-mad,-ess_tail) %>%
    rename(n_eff = ess_bulk)

  if(all(is.na(.s$beach))) .s$beach <- NULL

  .s$Y <- attr(model,'y.names')[.s$Y]
  rename(.s,y.var=Y)

}

#' @export

summary.phenology_df = function(object){

  mutate(object,summary = purrr::map(fit,summary.phenology_model))

}

# ----------------------------------------------------------------------------------------------------------
# season_merge: Helper function to combine predictions for all sites within seasons for summarizing
# ----------------------------------------------------------------------------------------------------------

MT_merge_sites = function(obj){

  if(!all(has_name(obj,c('season','beach')))) stop("obj should be of class MT_tbl with a 'season' and 'beach' column")

  condense = function(.x) list(setNames(.x,obj$beach) %>% bind_rows(.id = 'beach'))
  group_by(obj,season) %>%
    summarise(across(any_of(c('data','summary','predict')),condense),
              across(any_of('fit'),list))

}

# ----------------------------------------------------------------------------------------------------------
# daily_means: Returns the mean and credible interval for counts on each day/site for plotting seasonal trends
# ----------------------------------------------------------------------------------------------------------

daily_means = function(preds,by_site=T,interval = 0.95,what = c('counts','proportions')){

  what = match.arg(what)

  group_by(preds,y.var,date,day,.draw) %>%
    {if(by_site & has_name(preds,'beach')) group_by(.,beach,.add=T) else .} %>%
    summarise(mu = sum(mu),.groups = 'keep') %>%
    {if(what=='proportions') {
      ungroup(.,date,day) %>%
      mutate(mu = mu/sum(mu)) %>%
      group_by(date,day,.add=T)} else .} %>%
    ungroup(.draw) %>%
    tidybayes::mean_qi(mu,.width = interval) %>%
    dplyr::select(-(.width:.interval)) %>%
    ungroup()

}

# ==========================================================================================================
# get_phenology
# ==========================================================================================================

# This is a series of S3 Methods to retrieve different metrics of interest from the fitted
# phenology_model, prediction objects, phenology_df.
# It uses predictions from the model if available, or generates them internally if
# not using draws and days. Avoids storing full posterior predictions (millions of rows) in dataframe
# which can make very large objects. Only calls predict once for all metrics requested which speeds
# it up.

#'
#'
#' @export
#'

get_phenology = function(object,what, by_site = TRUE,CI = 0.95, full.posterior = FALSE, draws = 1000, days,
                         quantile = c(0.025,0.975)) {

  .v = c('daily.means','daily.proportions','annual.totals','quantiles')
  if(!all(what %in% .v)) stop (paste("'what' must be:", paste(.v,collapse=', ')))

  UseMethod("get_phenology")

}

#'
#'
#' @export
#'

get_phenology.MTpred = function(object,what, CI = 0.95, by_site = TRUE,
                                full.posterior = FALSE){

  # Change to getFromNamespace to make more robust
  .f = getFromNamespace(what,'TurtTools')
  .f(object,CI = CI, by_site = by_site, full.posterior = full.posterior)

}

#'
#'
#' @export
#'

get_phenology.phenology_model = function(object,..., draws = 1000, days){

  obj = predict(object,draws = draws,days = days, full.posterior = TRUE)
  get_phenology.MTpred(obj,...)

}

#'
#'
#' @export
#'

get_phenology.phenology_df = function(object,what,by_site = TRUE,...,draws = 1000, days){

  if(!has_name(object,'predict')) obj = predict(object,draws = draws,days = days, full.posterior=TRUE)

  if(by_site==FALSE & has_name(object,'beach')) obj = merge_sites(obj)

  for(i in what){
    obj = mutate(obj, !!i := purrr::map(predict,get_phenology,what=i,by_site=by_site,...))
  }

  dplyr::select(obj,any_of(c('season','beach',what)))

}

# ----------------------------------------------------------------------------------------------------------
# Series of functions used internally by get_phenology
# ----------------------------------------------------------------------------------------------------------

daily.means = function(preds,by_site=T,CI=0.95,...){

  group_by(preds,y.var,date,day,.draw) %>%
    {if(by_site & has_name(preds,'beach')) group_by(.,beach,.add=T) else .} %>%
    summarise(mu = sum(mu),.groups = 'keep') %>%
    ungroup(.draw) %>%
    tidybayes::mean_qi(mu,.width = CI) %>%
    dplyr::select(-(.width:.interval)) %>%
    ungroup()

}


daily.proportions = function(preds,by_site=T,CI=0.95,...){

  group_by(preds,y.var,date,day,.draw) %>%
    {if(by_site & has_name(preds,'beach')) group_by(.,beach,.add=T) else .} %>%
    summarise(mu = sum(mu),.groups = 'keep') %>%
    ungroup(.,date,day) %>%
    mutate(mu = mu/sum(mu)) %>%
    group_by(date,day,.add=T) %>%
    ungroup(.draw) %>%
    tidybayes::mean_qi(mu,.width = CI) %>%
    dplyr::select(-(.width:.interval)) %>%
    ungroup()

}


annual.totals = function(preds,by_site=T,CI=0.95,full.posterior=F,...){

  if(!inherits(preds,'MTpred')) stop ('predictions should be an MTpred object')

  preds =
    group_by(preds,y.var,.draw) %>%
    {if(by_site & 'beach' %in% names(.)) group_by(.,beach,.add=T) else .} %>%
    summarise(N = sum(N),.groups = 'keep')

  if(full.posterior) return(ungroup(preds))

  preds = ungroup(preds,.draw)
  tidybayes::mean_qi(preds,N,.width = CI) %>%
    inner_join(summarise(preds,log.sd = sd(log(N)))) %>%
    dplyr::select(-(.width:.interval)) %>%
    ungroup()

}

# Need to adjust this to cope with cases with or without site names.
annual.proportions = function(preds,by_site=T,CI=0.95,interval=0.95,full.posterior=F,...){

  group_by(preds,y.var,.draw) %>%

    {if(has_name(preds,'beach')){
        group_by(.,beach,.add=T) %>%
        summarise(N = sum(N),.groups = 'keep')
        ungroup(beach) %>%
        mutate(proportion = N/sum(N)) %>%
        group_by(beach,.add=T)
      } else {
        summarise(N = sum(N),.groups = 'keep')
      }
    } %>%

    {if(!full.posterior){
      ungroup(.,.draw) %>%
        tidybayes::mean_qi(proportion,.width = interval) %>%
        dplyr::select(-(.width:.interval))
    } else .} %>%

    ungroup()

}

# Check how this one works across multiple sites....might need to group by day too
quantiles = function(preds,by_site=T,CI=0.95,full.posterior=F,quantile = c(0.025,0.975),...){

  group_by(preds,y.var,.draw) %>%
    {if(by_site & has_name(preds,'beach')) group_by(.,beach,.add=T) else .} %>%
    mutate(s = cumsum(mu)/sum(mu)) %>%
    summarise(start = day[s>=quantile[1]][1],
              end = day[s>=quantile[2]][1],
              duration=end-start) %>%
    {if(!full.posterior){
      ungroup(.,.draw) %>%
        tidybayes::mean_qi(start,end,duration,.width = CI) %>%
        dplyr::select(-(.width:.interval))
    } else .} %>%
    ungroup()

}


# ----------------------------------------------------------------------------------------------------------
# Helper function to merge data within seasons for calculation of overall totals
# ----------------------------------------------------------------------------------------------------------

merge_sites = function(obj){

  if(!is(obj,'phenology_df')) stop("obj should be of class phenology_df")

  condense = function(.x) list(setNames(.x,obj$beach) %>% bind_rows(.id = 'beach'))
    dplyr::select(obj,season,any_of(c('beach','data','predict','summary'))) %>%
    group_by(season) %>%
    summarise(across(any_of(c('data','predict','summary')),condense))
}

