## ------------------------------------------------------
## Example 2: change in a time series of distributions
## ------------------------------------------------------

## Using Bayesian online changepoint detection for trend change
## in the distance measure before 
## For now, we only consider the i.i.d noise case 

source("BOCPD_functions.R")

.libPaths(c("~/R/library", .libPaths()))
library(transport)
# library(DescTools)



## Functions to generate simulated time series (parameter change)
CIsim.cdf = function(seed=10, npre, npost, level, trend, ns=100) {
  # level and trend are vectors corresponding to (mu, sigma)
  # set trend to 0 if we want to keep some of the parameters fixed
  set.seed(seed)
  nt = npre + npost
  
  # simulate two series of samples
  mu = level[1]
  sigma = level[2]
  trend.interv.mu = trend[1]
  trend.interv.sigma = trend[2]
  
  mu.control = rep(mu, times=nt)
  mu.interv = c(rep(mu, times=npre), mu + trend.interv.mu * (1:npost))
  sigma.control = rep(sigma, times=nt)
  sigma.interv = c(rep(sigma, times=npre), sigma + trend.interv.sigma * (1:npost))
  
  sample.control = sapply(1:length(mu.control), FUN=function(x) { 
    rnorm(ns, mean=mu.control[x], sd=sigma.control[x]) })
  sample.interv = sapply(1:length(mu.interv), FUN=function(x) { 
    rnorm(ns, mean=mu.interv[x], sd=sigma.interv[x]) })
  
  return(list(sample.ctr=sample.control, sample.itv=sample.interv))
}

## The distance measures
# Wasserstein distance

WSdist = function(sample.ctr, sample.itv, bw=0.5, nd=50) {
  ns = ncol(sample.ctr)
  
  allWS = sapply(1:ns, FUN=function(x) {
    # Wasserstein between the empirical distributions (1d)
    y.ctr = sample.ctr[, x]
    y.itv = sample.itv[, x]
    WS = wasserstein1d(a=y.ctr, b=y.itv, p=1) 
    return(WS)
  })
  return(allWS)
}

# WSdist(sample.ctr = sample.control, sample.itv = sample.interv,
#        bw=0.3, nd=50)

## Area-under-curve

AUCdiff = function(sample.ctr, sample.itv, bw=0.5, nd=50) {
  ns = ncol(sample.ctr)
  
  allAUC = sapply(1:ns, FUN=function(x) {
    y.ctr = sample.ctr[, x]
    y.itv = sample.itv[, x]
    # we need common coordinate for the difference
    yrange = range(c(y.ctr, y.itv))
    d.ctr = density(y.ctr, bw=bw, n=nd, from=yrange[1], to=yrange[2])
    d.itv = density(y.itv, bw=bw, n=nd, from=yrange[1], to=yrange[2])
    d.dif = d.itv$y - d.ctr$y
    
    dAUC = AUC(x=d.ctr$x, y=d.dif, method="trapezoid", absolutearea=TRUE)
    return(dAUC)
  })
  return(allAUC)
}

# AUCdiff(sample.ctr = sample.control, sample.itv = sample.interv,
#         bw=0.3, nd=50)

## Parametric test (change in trend in parameters)
trend.stats.cdf = function(y, nt, monthly=TRUE) {
  if (monthly == TRUE) {
    mint = 12   # because we know the first one or two years are pre-intervention
    maxt = nt - 3
    seglen = mint:maxt
  } else {
    mint = 2   
    maxt = nt - 1
    seglen = mint:maxt
  }
  
  stats.vec = sapply(1:length(seglen), FUN=function(m){
    seg1 = 1:seglen[m]
    seg2 = (seglen[m]+1):nt
    Xmat = cbind(rep(1, nt), 
                 c(seg1, rep(seg1[length(seg1)], length(seg2))),
                 c(rep(0, length(seg1)), seg2-seg2[1]+1))
    
    lmfit = lm(y ~ -1 + Xmat)
    vec = c(0, -1, 1)
    sddiff = sqrt(t(vec) %*% summary(lmfit)$cov.unscaled %*% vec)
    stats = (lmfit$coef[2] - lmfit$coef[3]) / sddiff
  })
  
  return(c(Tmax = max(abs(stats.vec)), 
           cpt = which.max(abs(stats.vec))[1] + mint))
}


## ---------------------------
## Run the simulation study
## ---------------------------
## Simulation set up
npre = 24                        # 2 years pre-intervention
npost.vec = seq(24, 120, by=3)   # from 2 to 10 years after intervention

npost.max = npost.vec[length(npost.vec)]
ntt = npre + npost.max

mu = 10
sigma = 1
# trend increment, from 5% to 100% over 10 years
trend.increase.mu = round(mu * c(0.05, seq(0.1, 1, by=0.1)) / 120, 4)
# trend increment, from 50% to 200% over 10 years
trend.increase.sigma = round(sigma * seq(0.5, 2, by=0.5) / 120, 4)

# delay in effect, up to 20 time points
delay.set = 1:20


## Simulation to get the critical value
dist.measure = c("Wasserstein", "AUC")[1]   # choose a distance measure

alpha = 0.95   # significance level
bw = 0.3    # band width for density estimation
nd = 50     # number of division for density estimation

## Basic information and global parameters
alpha.ini <- c(0, 0.05)             # the trend model
sigma.ini <- (0.05)^2
params.ini.list <- list(list(alpha=alpha.ini, sigma=sigma.ini))

a.prior <- c(alpha.ini, 0.5/sigma.ini, 0.5/sigma.ini)     # Normal distribution
sig.prior <- c(100, 1)                                    # the scale of sd
prior.list <- list(list(a.prior=a.prior, sig.prior=sig.prior))

nmodel <- 1
pm <- 1      # prior model probability
ptr <- 0.02    # run-length distribution
global.params <- list(pm, ptr)

## likelihood
intfun.list <- list(integrate.fun5)

maxp = 50
np = 40
msl = 10

## Start the repetitions
simN = 1000
simResult = lapply(1:simN, FUN=function(s) {
  
  bocpd.list = sapply(1:length(trend.increase.mu), FUN=function(m) {
    
    # simulate time series with different scales of trend
    trend.interv.mu = trend.increase.mu[m]
    
    # if using delay in effect, otherwise, set delay to 0
    delay = sample(delay.set, size=1)
    # delay = 0
    npre.delay = npre + delay
    npost.delay = npost.max - delay
    sim.sample = CIsim.cdf(seed=s*length(trend.increase.mu)+m, 
                           npre=npre.delay, npost=npost.delay, 
                           level=c(mu, sigma), 
                           trend=c(trend.interv.mu, 0), 
                           ns=200)
    
    if (dist.measure == "Wasserstein") {
      diff.sample = WSdist(sample.ctr = sim.sample$sample.ctr,
                           sample.itv = sim.sample$sample.itv,
                           bw=bw, nd=nd)
    } else if (dist.measure == "AUC") {
      diff.sample = AUCdiff(sample.ctr = sim.sample$sample.ctr,
                            sample.itv = sim.sample$sample.itv,
                            bw=bw, nd=nd)
    }
    
    ## The model object
    # For the online algorithm (reads in data one by one)
    x = diff.sample
    hmm = NULL
    hmm$x <- as.matrix(x[1:(2*msl-1)], ncol=1)
    hmm$Nx <- nrow(hmm$x)
    hmm$pModels <- pm
    
    # hmm$transition <- cbind(rep(1-ptr, hmm$Nx), rep(ptr, hmm$Nx))   
    # column 1: stay, column 2: changepoint
    runl.fun <- runl.fun2
    Gxt0 <- sapply(msl:length(x), FUN=function(x) { runl.fun(x, ptr, msl)$Gx })
    Gxt1 <- sapply(msl:length(x)-1, FUN=function(x) { runl.fun(x, ptr, msl)$Gx })
    stay <- c(rep(1, msl-1), (1 - Gxt0) / (1 - Gxt1))
    change <- c(rep(0, msl-1), 1 - (1 - Gxt0) / (1 - Gxt1))
    hmm$transition <- cbind(stay[1:(2*msl-1)], change[1:(2*msl-1)])
    
    ## For the main recursion
    hmm$intPar <- list()    # the integrated log-likelihood
    hmm$intExtra <- list()  # additional calculation when adding a new cpt
    hmm$forward <- list()   # the log-posterior of the cpts
    hmm$ParsID <- list()    # indices of potential changepoints
    hmm$tmpCost <- list()   # tmplike in PELT code, the current loglike
    
    hmm$forward[[msl-1]] <- 1
    
    # from msl to 2*msl-1 (last cpt = 0)
    params.list <- params.ini.list
    
    # additional calculation for t=msl-1 (in order to get p(y_t | y_1:t))
    intPar.m <- integrate.par(hmm, cpt=msl-1, cps=0, nmodel=1, params.list,
                              prior.list, intfun.list)
    intPar <- unlist(intPar.m)
    hmm$intPar[[msl-1]] <- list(intPar)
    
    # the initial iterations before introducing the first changepoint
    for (i in msl:(2*msl-1)) {
      # all segments start at t=1
      intPar.m <- integrate.par(hmm, cpt=i, cps=0, nmodel=1, params.list,
                                prior.list, intfun.list)
      intPar <- unlist(intPar.m)
      # store the result
      hmm$intPar[[i]] <- list(matrix(intPar, nrow=nmodel, ncol=1))
      
      # calculate the forward quantity
      # don't need to calculate the I_t^{s, m} here
      hmm$forward[[i]] <- 1
      hmm$ParsID[[i]] <- i    # distance to the last cpt
      
      # calculate the probabilities used in the MAP estimation (sum over models)
      intPar.sum <- sapply(1:nmodel, FUN=function(m) {
        hmm$intPar[[i]][[m]][, 1] })
      hmm$tmpCost[[i]] <- logsumexp(matrix(intPar.sum + log(hmm$pModels), nrow=1)) + 
        log(runl.fun(i, ptr, msl)$gx) 
    }
    
    cpt.map = list(breaks=0, models=1)
    
    # from 2*msl to N
    # stop at first detection
    i = 2*msl
    while((length(cpt.map$breaks) == 1) & (i <= length(x))) {
      # read in new data
      hmm$x = as.matrix(x[1:i], ncol=1)
      hmm$Nx = nrow(hmm$x)
      hmm$transition <- cbind(stay[1:i], change[1:i])
      
      # number of changepoints
      l <- length(hmm$forward[[i-1]]) + 1
      pos <- hmm$ParsID[[i-1]]   
      # pos are the distances of the changepoints to i
      # i-pos are locations of the start of the segments
      
      # compute the integration I_t^{s, m} and update the parameters
      Integrate <- forward.cpt.int(hmm, cpt=i, nmodel, params.ini.list, 
                                   prior.list, intfun.list, msl)
      hmm <- Integrate$hmm
      hmm$intPar[[i]] <- Integrate$intPar.list
      hmm$intExtra[[i]] <- Integrate$intExtra.list
      
      # compute the forward quantities, the weights
      # and resample the changepoints
      hmm = forward.cpt.calc.resample(hmm, cpt=i, nmodel, maxp, np,
                                      prune=TRUE, delay=5, ptr, msl)
      
      # Run MAP and stop the algorithm if a changepoint is found
      cpt.map <- map.cpt(hmm, nmodel=1, global.params, runl.fun, msl)
      
      i = i + 1
      # if ((i %% round(length(x)/10)) == 0) { cat(".") }
    }
    
    cpt.est = cpt.map$breaks[2]
    time.est = i      # time when the first changepoint is declared
    c(cpt.est=cpt.est, time.est=time.est, delay=delay)
  })
  
  return(bocpd.list)
})

save(simResult, file="sim_result_distribution_change_BOCPD.RData")


