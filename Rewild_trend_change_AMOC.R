## ------------------------------------------------
## Example 1: trend change in time series data
## ------------------------------------------------

# In this case, we observe the time series with trend change directly

# npre: number of samples before intervention
# npost: number of samples after intervention
# tau: the changepoint
# level: the mean level of the time series
# trend.control: the trend of the time series before intervention (often set to 0) 
# trend.interv: the trend of the intervention time series after intervention
# sigma: variance of the noise
# phi: AR(1) coefficient

# The changepoint detection is carried out on the difference time series
# y_t = y_interv - y_control
# y_t = level + noise   (for t <= tau)
# y_t = level + trend * (t-tau) + noise  (for t > tau)


## ----------------------
## Essential functions
## ----------------------
## Functions to generate simulated time series
CIsim = function(seed=10, npre, npost, level, trend, sigma) {
  set.seed(seed)
  nt = npre + npost
  
  trend.control = trend[1]
  trend.interv = trend[2]
  # with i.i.d. noises
  y.control = level + trend.control * 1:nt + rnorm(nt, mean=0, sd=1) * sigma
  y.interv = c(level + trend.control * 1:npre + rnorm(npre, mean=0, sd=1) * sigma,  
               level + trend.control * npre + trend.interv * 1:npost + 
               rnorm(npost, mean=0, sd=1) * sigma)
  
  return(list(y.ctr=y.control, y.itv=y.interv))
}

## Functions to calculate the test statistics of a trend difference
# beta_1 v.s. beta_0
trend.stats = function(y.ctr=NULL, y.itv, nt, monthly=TRUE) {
  # If we only have the intervention time series, set y.ctr to NULL
  
  if (monthly == TRUE) {
    mint = 24   # because we know the first one or two years are pre-intervention
    maxt = nt - 4
    seglen = mint:maxt
  } else {
    mint = 2
    maxt = nt - 1
    seglen = mint:maxt
  }
  
  # evaluate the test statistics for different choices of changpoint
  stats.vec = sapply(1:length(seglen), FUN=function(m){
    seg1 = 1:seglen[m]
    seg2 = (seglen[m]+1):nt
    Xmat = cbind(rep(1, nt), 
                 c(seg1, rep(seg1[length(seg1)], length(seg2))),
                 c(rep(0, length(seg1)), seg2-seg2[1]+1))
    
    if (!is.null(y.ctr)) {
      y.dif = y.itv - y.ctr
    } else { y.dif = y.itv }
    
    # fit simple trend model (linear regression)
    lmfit = lm(y.dif ~ -1 + Xmat)
    vec = c(0, -1, 1)
    sddiff = sqrt(t(vec) %*% summary(lmfit)$cov.unscaled %*% vec)
    stats = (lmfit$coef[2] - lmfit$coef[3]) / sddiff
  })
  
  # find the time point (tau) that maximises the test statistics
  # tau could be a potential changepoint, if it is greater that the critical value
  return(c(Tmax = max(abs(stats.vec)), 
           cpt = which.max(abs(stats.vec))[1] + mint))
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
# trend increment, as percentage change to the mean before intervention
trend.increase = round(level * c(0.05, seq(0.1, 1, by=0.1)) / 120, 4) 
sigma = 0.05            # noise level
# delay in effect, up to 20 time points
delay.set = 1:20


## (1) Simulation CIBA: look at the difference between intervention and control
## Simulation to get the critical values of the trend difference test
alpha = 0.95   # significance level
Nsim = 1000
getCritical.v1 = sapply(npost.vec[c(3, 7)], FUN=function(n) {
  npost = n
  Tmax.list = sapply(1:Nsim, FUN=function(i) {
    nt = npre + npost
    # the null hypothesis has trend 0
    sim.ts = CIsim(seed=i, npre, npost, level, 
                   trend=c(trend.control, trend.control), sigma)
    testStats = trend.stats(y.ctr=sim.ts$y.ctr, y.itv=sim.ts$y.itv, nt=nt,
                            monthly=TRUE)
    return(testStats)
  })
  criticVal = quantile(Tmax.list[1, ], probs=alpha)
  return(criticVal)
})

simN = 1000
time.start = proc.time()
simResult.v1 = lapply(1:simN, FUN=function(i) {
 
  teststat.list = lapply(1:length(trend.increase), FUN=function(m) {
    
    trend.interv = trend.control + trend.increase[m]
    
    delay = sample(delay.set, size=1)   # sample a delay period
    # delay = 0   # if there is no delay effect, set it to 0
    npre.delay = npre + delay
    npost.delay = npost.max - delay
    sim.ts = CIsim(seed=i*length(trend.increase) + m, 
                   npre=npre.delay, npost=npost.delay, level, 
                   trend=c(trend.control, trend.interv), sigma)
    
    # for each trend coefficient, simulate time series of increasing length
    teststat = sapply(npost.vec, FUN=function(n) {
      npost = n
      nt = npre + npost
      testStats = trend.stats(y.ctr=sim.ts$y.ctr[1:nt], 
                              y.itv=sim.ts$y.itv[1:nt], 
                              nt=nt, monthly=TRUE)
      return(testStats)
    })
    return(list(teststat=teststat, delay=delay))
  })
  return(teststat.list)
})
# computation time
proc.time() - time.start

# save the result
save(simResult.v1, getCritical.v1, file="sim_result_trend_change_AMOC.RData")


## (2) Simulation: BA, when there is only an intervention time series
getCritical.v2 = sapply(npost.vec[c(3, 7)], FUN=function(n) {
  npost = n
  Tmax.list = sapply(1:Nsim, FUN=function(i) {
    nt = npre + npost
    # the null hypothesis has trend 0
    sim.ts = CIsim(seed=i, npre, npost, level, 
                   trend=c(trend.control, trend.control), sigma)
    testStats = trend.stats(y.ctr=NULL, y.itv=sim.ts$y.itv, 
                            nt=nt, monthly=TRUE)
    return(testStats)
  })
  criticVal = quantile(Tmax.list[1, ], probs=alpha)
  return(criticVal)
})

time.start = proc.time()
simResult.v2 = lapply(1:simN, FUN=function(i) {
 
  teststat.list = lapply(1:length(trend.increase), FUN=function(m) {
    
    trend.interv = trend.control + trend.increase[m]
    
    delay = sample(delay.set, size=1)   # sample a delay period
    # delay = 0   # if there is no delay effect, set it to 0
    npre.delay = npre + delay
    npost.delay = npost.max - delay
    sim.ts = CIsim(seed=i*length(trend.increase) + m, 
                   npre=npre.delay, npost=npost.delay, level, 
                   trend=c(trend.control, trend.interv), sigma)
    
    # for each trend coefficient, simulate time series of increasing length
    teststat = sapply(npost.vec, FUN=function(n) {
      npost = n
      nt = npre + npost
      testStats = trend.stats(y.ctr=NULL, 
                              y.itv=sim.ts$y.itv[1:nt], 
                              nt=nt, monthly=TRUE)
      return(testStats)
    })
    return(list(teststat=teststat, delay=delay))
  })
  return(teststat.list)
})
# computation time
proc.time() - time.start

# save the result
save(simResult.v1, getCritical.v1, simResult.v2, getCritical.v2, 
     file="sim_result_trend_change_AMOC.RData")
