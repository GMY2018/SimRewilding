## ------------------------------------------------------
## Example 2: change in a time series of distributions
## ------------------------------------------------------

## Assuming we only observe one time seires, with before and after prevention 
## data. Each time, we observe a sample (distributions)

## Changepoint detection strategy
## First estimate a before prevention density using all data from the 
## pre-intervention period, call it "baseline density".
## Then we evaluate the distance measure between the observed distribution at
## each time point and the "baseline density". 
## We look for change in the distance measure time series

library(transport)


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
WSdist2 = function(sample.itv, npre, bw=0.5, nd=50) {
  ns = ncol(sample.itv)
  y.ctr = sample.itv[, 1:npre]
  
  allWS = sapply(1:ns, FUN=function(x) {
    # compare the current distribution to the baseline distribution
    y.itv = sample.itv[, x]
    WS = wasserstein1d(a=y.ctr, b=y.itv, p=1) 
    return(WS)
  })
  
  return(allWS)
}

# WSdist2(sample.itv = sample.interv, npre=24, bw=0.3, nd=50)

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
npre = 24                         # 2 years pre-intervention
npost.vec = seq(24, 120, by=10)   # from 2 to 10 years after intervention

npost.max = npost.vec[length(npost.vec)]
ntt = npre + npost.max

mu = 10
sigma = 1
# trend increment, from 5% to 100% over 10 years
trend.increase.mu = round(mu * c(0.05, seq(0.1, 1, by=0.1)) / 120, 4)
# trend increment, from 50% to 200% over 10 years
trend.increase.sigma = round(sigma * seq(0.5, 2, by=0.5) / 120, 4)

# delay in effect
delay.set = 1:20

## Simulation to get the critical value
alpha = 0.95   # significance level
bw = 0.3    # band width for density estimation
nd = 50     # number of division for density estimation

Nsim = 1000
getCritical = sapply(npost.vec[c(8, 24)], FUN=function(n) {
  npost = n
  
  Tmax.list = sapply(1:Nsim, FUN=function(i) {
    nt = npre + npost
    sim.sample = CIsim.cdf(seed=i, npre, npost, level=c(mu, sigma), 
                           trend=c(0, 0), ns=200)
    
    # calculate the distance measure
    diff.sample = WSdist2(sample.itv = sim.sample$sample.itv,
                          npre=npre, bw=bw, nd=nd)
    
    testStats = trend.stats.cdf(y=diff.sample, nt=nt, monthly=TRUE)
    return(testStats)
  })
  
  criticVal = quantile(Tmax.list[1, ], probs=0.95)
  return(criticVal)
})


## Repeat a large number of times
simN = 1000
time.start = proc.time()
simResult.mu = lapply(1:simN, FUN=function(i) {
  
  teststat.list = lapply(1:length(trend.increase.mu), FUN=function(m) {
    
    trend.interv.mu = trend.increase.mu[m]
    
    npost.max = npost.vec[length(npost.vec)]
    ntt = npre + npost.max
    
    # if using delay in effect, otherwise, set delay to 0
    delay = sample(delay.set, size=1)
    # delay = 0
    npre.delay = npre + delay
    npost.delay = npost.max - delay
    sim.sample = CIsim.cdf(seed=i*length(trend.increase.mu)+m, 
                           npre=npre.delay, npost=npost.delay, 
                           level=c(mu, sigma), 
                           trend=c(trend.interv.mu, 0), 
                           ns=200)
    
    teststat = sapply(npost.vec, FUN=function(n) {
      npost = n
      nt = npre + npost
      
      # calculate the distance measure
      diff.sample = WSdist2(sample.itv = sim.sample$sample.itv[, 1:nt],
                            npre=npre, bw=bw, nd=nd)
      
      testStats = trend.stats.cdf(y=diff.sample, nt=nt, monthly=TRUE)
      return(testStats)
    })
    return(list(teststat=teststat, delay=delay))
  })
  return(teststat.list)
})
proc.time() - time.start

# save the result
save(simResult.mu, getCritical, file="sim_result_distribution_change_SingleTS.RData")

