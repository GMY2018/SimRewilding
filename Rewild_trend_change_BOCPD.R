## ------------------------------------------------
## Example 1: trend change in time series data
## ------------------------------------------------

## Using Bayesian online changepoint detection (BOCPD) for detecting trend change
## Using functions for the standard BOCPD (without SMC or gradient ascent)
## For now, only consider the i.i.d noise case 

source("BOCPD_functions.R")


## Functions to generate simulated time series
CIsim = function(seed=10, npre, npost, level, trend, sigma) {
  set.seed(seed)
  nt = npre + npost
  
  trend.control = trend[1]
  trend.interv = trend[2]
  # with i.i.d. noises
  y.control = level + trend.control * 1:nt + rnorm(nt, mean=0, sd=1) * sigma
  y.interv = c(level + trend.control * 1:npre + rnorm(npre, mean=0, sd=1) * sigma,  
               level + trend.control * npre + trend.interv * 1:npost + rnorm(npost, mean=0, sd=1) * sigma)
  
  return(list(y.ctr=y.control, y.itv=y.interv))
}


## ---------------------------
## Run the simulation study
## ---------------------------
## Simulation set up, assuming a higher sampling frequency
npre = 24                         # e.g., 24 months (2 years) pre-intervention
npost.vec = seq(24, 120, by=12)   # e.g., from 2 to 10 years after intervention

npost.max = npost.vec[length(npost.vec)]
ntt = npre + npost.max

level = 10 
trend.control = 0.005
# trend increment, as percentage change to the mean level before intervention
trend.increase = round(level * c(0.05, seq(0.1, 1, by=0.1)) / 120, 4)  

sigma = 0.05     # noise level

# delay in taking effect
delay.set = 1:20

## Basic information and global parameters
alpha.ini <- c(0, 0.05)             # the trend model
sigma.ini <- (0.05)^2
params.ini.list <- list(list(alpha=alpha.ini, sigma=sigma.ini))

a.prior <- c(alpha.ini, 0.5/sigma.ini, 0.5/sigma.ini)     # Normal distribution
sig.prior <- c(100, 1)                                    # the scale of sd
prior.list <- list(list(a.prior=a.prior, sig.prior=sig.prior))

nmodel <- 1
pm <- 1        # prior model probability
ptr <- 0.02    # run-length distribution
global.params <- list(pm, ptr)

## likelihood
intfun.list <- list(integrate.fun5)

maxp = 50
np = 40
msl = 10

## Start the repetitions
Nsim = 1000
simResult = lapply(1:Nsim, FUN=function(s) {
  
   bocpd.list = sapply(1:length(trend.increase), FUN=function(m) {
    
    # simulate time series with different scales of trend
    trend.interv = trend.increase[m]
    
    delay = sample(delay.set, size=1)   # if there is a delay effect
    # delay = 0
    npre.delay = npre + delay
    npost.delay = npost.max - delay
    sim.ts = CIsim(seed=s*length(trend.increase)+m, 
                   npre=npre.delay, npost=npost.delay, level=level, 
                   trend=c(trend.control, trend.interv), sigma=sigma)
    
    ## The model object
    # For the online algorithm (reads in data one by one)
    x = sim.ts$y.itv - sim.ts$y.ctr
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
    return(c(cpt.est=cpt.est, time.est=time.est, delay=delay))
  })
  
  return(bocpd.list)
})

save(simResult, file="sim_result_trend_change_BOCPD.RData")


