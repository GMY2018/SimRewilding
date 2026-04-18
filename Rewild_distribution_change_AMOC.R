## ------------------------------------------------------
## Example 2: change in a time series of distributions
## ------------------------------------------------------

# In this case, we observe distributions each time, and we want to detect changes
# in the distributions (e.g., parameter changes, distribution type changes)

# npre: number of samples before intervention
# npost: number of samples after intervention
# tau: the changepoint
# mu: location parameter in the control and intervention time series
# sigma: scale parameter in the control and intervention time series
# trend.interv.mu: location parameter after intervention
# trend.interv.sigma: scale parameter after intervention

# The changepoint detection method involves the evaluation of a distance measure between
# the sample distributions from the control and intervention time series at each
# time point. We then look for the change in the distance measure to identify
# the time point when the distributions from the control and the intervention 
# time series start to differ.

library(transport)


## ----------------------
## Essential functions
## ----------------------
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
# library(transport)

# a function for this task
WSdist = function(sample.ctr, sample.itv, bw=0.5, nd=50) {
  ns = ncol(sample.ctr)
  
  allWS = sapply(1:ns, FUN=function(x) {
    # Wasserstein between two empirical distributions (1d)
    y.ctr = sample.ctr[, x]
    y.itv = sample.itv[, x]
    WS = wasserstein1d(a=y.ctr, b=y.itv, p=1) 
    return(WS)
  })
  
  return(allWS)
}

# WSdist(sample.ctr = sample.control, sample.itv = sample.interv,
#        bw=0.3, nd=50)

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
npost.vec = seq(24, 120, by=12)   # from 2 to 10 years after intervention

npost.max = npost.vec[length(npost.vec)]
ntt = npre + npost.max

mu = 10
sigma = 1
# trend increment, from 5% to 100% over 10 years
trend.increase.mu = round(mu * c(0.05, seq(0.1, 1, by=0.1)) / 120, 4)
# trend increment, from 50% to 200% over 10 years
trend.increase.sigma = round(sigma * seq(0.5, 2, by=0.5) / 120, 4)

# if we want to test the combinations of changes
# trend.combine = cbind(trend.increase.mu, rep(trend.increase.sigma, each=3)[-1])

# delay in effect
delay.set = 1:20

## Simulation to get the critical value
alpha = 0.95   # significance level
bw = 0.3    # band width for density estimation
nd = 50     # number of division for density estimation

Nsim = 1000
getCritical = sapply(npost.vec[c(1, 17, 33)], FUN=function(n) {
  npost = n
  
  Tmax.list = sapply(1:Nsim, FUN=function(i) {
    nt = npre + npost
    sim.sample = CIsim.cdf(seed=i, npre, npost, level=c(mu, sigma), 
                           trend=c(0, 0), ns=200)
    
    # calculate the Wasserstein distance measure
    diff.sample = WSdist(sample.ctr = sim.sample$sample.ctr,
                         sample.itv = sim.sample$sample.itv,
                         bw=bw, nd=nd)
    testStats = trend.stats.cdf(y=diff.sample, nt=nt, monthly=TRUE)
    return(testStats)
  })
  
  criticVal = quantile(Tmax.list[1, ], probs=0.95)
  return(criticVal)
})

save(getCritical, file="sim_result_distribution_change_AMOC.RData")


## Repeat a large number of times
simN = 1000
time.start = proc.time()
simResult.mu = lapply(1:simN, FUN=function(i) {
  
  teststat.list = lapply(1:length(trend.increase.mu), FUN=function(m) {
    
    trend.interv.mu = trend.increase.mu[m]
    
    delay = sample(delay.set, size=1)
    # delay = 0
    npre.delay = npre + delay
    npost.delay = npost.max - delay
    sim.sample = CIsim.cdf(seed=i*length(trend.increase.mu)+m, 
                           npre=npre.delay, npost=npost.delay, level=c(mu, sigma), 
                           trend=c(trend.interv.mu, 0), ns=200)
    
    teststat = sapply(npost.vec, FUN=function(n) {
      npost = n
      nt = npre + npost
      
      # calculate the distance measure
      diff.sample = WSdist(sample.ctr = sim.sample$sample.ctr[, 1:nt],
                           sample.itv = sim.sample$sample.itv[, 1:nt],
                           bw=bw, nd=nd)
      testStats = trend.stats.cdf(y=diff.sample, nt=nt, monthly=TRUE)
      return(testStats)
    })
    return(list(teststat=teststat, delay=delay))
  })
  return(teststat.list)
})
proc.time() - time.start

# save the result
save(simResult.mu, getCritical, file="sim_result_distribution_change_AMOC.RData")


## change in variance
time.start = proc.time()
simResult.sigma = lapply(1:simN, FUN=function(i) {
  
  teststat.list = lapply(1:length(trend.increase.sigma), FUN=function(m) {
    
    trend.interv.sigma = trend.increase.sigma[m]
    
    delay = sample(delay.set, size=1)
    # delay = 0
    npre.delay = npre + delay
    npost.delay = npost.max - delay
    sim.sample = CIsim.cdf(seed=i*length(trend.increase.sigma)+m, 
                           npre=npre.delay, npost=npost.delay, level=c(mu, sigma), 
                           trend=c(0, trend.interv.sigma), ns=200)
    
    teststat = sapply(npost.vec, FUN=function(n) {
      npost = n
      nt = npre + npost
      
      # calculate the distance measure
      diff.sample = WSdist(sample.ctr = sim.sample$sample.ctr[, 1:nt],
                           sample.itv = sim.sample$sample.itv[, 1:nt],
                           bw=bw, nd=nd)
      testStats = trend.stats.cdf(y=diff.sample, nt=nt, monthly=TRUE)
      return(testStats)
    })
    return(list(teststat=teststat, delay=delay))
  })
  return(teststat.list)
})
proc.time() - time.start


# save the result
save(simResult.mu, simResult.sigma, getCritical, 
     file="sim_result_distribution_change_AMOC.RData")


