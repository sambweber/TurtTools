# -----------------------------------------------------------------------
# dSnbinomNim
# -----------------------------------------------------------------------

#' This is a custom probability distribution giving the density for the sum of mutiple
#' random variables with negative binomial distributions (e.g. in this case, numbers of activities over multiple days).
#' The function is a `Nimble` implementation of the function in package `HelpersMG`.

#' @param x,size,mu,nt,log Definitions as given in HelpersMG::dSnbinom

dSnbinomNim <- nimble::nimbleFunction(
  run = function(x = double(0), size = double(0), mu = double(1),
                 nt = integer(0), log = integer(0, default = 0)) {
    returnType(double(0))
    if(nt < 1) {
      print("’nt’ must be >= 1")
      return(NA)
    }
    if(nt > length(mu)) {
      print("’nt’ must be less than the length of ’mu’")
      return(NA)
    }
    mu1 <- numeric(nt)
    for(i in 1:nt) mu1[i] <- mu[i]
    if(any(mu1 < 0) | size <= 0) {
      print("dSnbinomNim must have positive mu and size")
      return(NA)
    }
    if(x < 0) {
      print("dSnbinomNim must have x >= 0")
      return(NA)
    }
    if(all(mu1 == 0) & x > 0) {
      if(log) {
        return(-Inf)
      } else {
        return(0) }
    }
    niter <- 1000
    tol <- 0.0000001
    p <- size / (size + mu1)
    if(length(mu1) == 1) {
      return(dnbinom(x, size, p[1], log = log))
    } else {
      q <- 1 - p
      p1 <- max(p)
      q1 <- 1 - p1
      qp1 <- (q * p1) / (q1 * p)
      R1 <- qp1
      for(i in 1:length(R1)) {
        R1[i] <- R1[i]^(-size)
      }
      R <- sum(log(R1))
      delta <- c(1, rep(NA, niter))
      xi <- rep(NA, niter)
      tqp1 <- 1 - ((q1 * p) / (q * p1))
      tqp2 <- tqp1
      xi[1] <- sum((size * tqp2))
      delta[2] <- xi[1] * delta[1]
      i <- 2
      val <- 2 * tol
      if(is.na(val) | is.nan(val)) return(NA)
      while(val >= tol & i <= niter) {
        for(j in 1:length(tqp1)) {
          tqp2[j] <- tqp1[j]^i
        }
        xi[i] <- sum((size * tqp2) / i)
        delta[i + 1] <- sum(c(1:i) * xi[1:i] * delta[i + 1 - c(1:i)]) / i
        val <- abs(xi[i] - xi[i - 1])
        if(is.na(val) | is.nan(val)) return(NA)
        i <- i + 1
      }
      delta1 <- delta[1:i]
      size1 <- length(mu1) * size + c(1:i) - 1
      p11 <- rep(p1, i)
      PrS <- sum(delta1 * dnbinom(x, size1, p11))
      if(log) {
        return(R + log(PrS))
      } else {
        return(exp(R) * PrS)
      }
    }
  })


# ----------------------------------------------------------------------------
# rSnbinomNim: A random number generator for the Snbinom distribution
# ----------------------------------------------------------------------------

# This assumes that nt (the sampling window in days) is always 1,
# as it always will be for predicting from unsampled days

rSnbinomNim <-nimble::nimbleFunction(
  run = function(n = integer(0), size = double(0), mu = double(1),
                 nt = integer(0)) {
    returnType(integer())
    prob <- size/(size+mu[1])
    return(rnbinom(1,size,prob))
  })
