## ------------------------------------------------------
## Example 2: change in a time series of distributions
## ------------------------------------------------------

# Using the changepoint.forecast package for online detection
# Detecting changepoint as a change in the forecast errors
# We consider location changes of the distribution

library(transport)
library(changepoint.forecast)


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

mu = 10
sigma = 1
# trend increment, from 5% to 100% over 10 years
trend.increase.mu = round(mu * c(0.05, seq(0.1, 1, by=0.1)) / 120, 4)
# trend increment, from 20% to 200% over 10 years
trend.increase.sigma = round(sigma * c(0.2, seq(0.5, 2, by=0.5)) / 120, 4)

# delay in effect, up to 20 time points
delay.set = 1:20


## Simulation to get the critical value
alpha = 0.95   # significance level
bw = 0.3    # band width for density estimation
nd = 50     # number of division for density estimation


## Simulation: CIBA (control-intervention, before-after)
simN = 1000
simResult.mu1 = lapply(1:simN, FUN=function(s) {
  
  change.list = sapply(1:length(trend.increase.mu), FUN=function(m) {
    
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
    
    diff.sample = WSdist(sample.ctr = sim.sample$sample.ctr,
                         sample.itv = sim.sample$sample.itv,
                         bw=bw, nd=nd)
    
    ## Train the model on pre-data, and then predict the post-data
    y.ts = diff.sample
    
    Xmat = data.frame(x1=rep(1, ntt), x2=1:ntt)
    lmfit = lm(y.ts[1:npre] ~ -1 + x1 + x2, data=Xmat[1:npre, ])
    # h-step ahead prediction
    lm.predict = predict(lmfit, newdata=Xmat[-(1:npre), ])
    r.ts = c(lmfit$residuals, lm.predict - y.ts[-(1:npre)])
    
    # changepoint detection on residual time series
    cpt.fcst = cptForecast(errors=as.vector(r.ts), m=npre,
                           detector="PageCUSUM", forecastErrorType="Both")
    time.est = cpt.fcst@tau      # based on errors
    
    # retrieve the changepoint using AMOC for trend change
    if (!is.infinite(time.est)) {
      nt = npre + time.est
      testStats = trend.stats(y.ctr=NULL, y.itv=r.ts[1:nt], 
                              nt=nt, monthly=TRUE)
      cpt.est = testStats[2]
    } else { cpt.est = Inf }
    
    return(c(cpt.est=cpt.est, time.est=time.est, delay=delay))
  })
  
  return(change.list)
})


## Simulation: BA (before-after), intervention time series only
simResult.mu2 = lapply(1:simN, FUN=function(s) {
  
  change.list = lapply(1:length(trend.increase.mu), FUN=function(m) {
    
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
    
    diff.sample = WSdist2(sample.itv = sim.sample$sample.itv,
                          npre=npre, bw=bw, nd=nd)
    
    ## Train the model on pre-data, and then predict the post-data
    y.ts = diff.sample
    
    Xmat = data.frame(x1=rep(1, ntt), x2=1:ntt)
    lmfit = lm(y.ts[1:npre] ~ -1 + x1 + x2, data=Xmat[1:npre, ])
    # h-step ahead prediction
    lm.predict = predict(lmfit, newdata=Xmat[-(1:npre), ])
    r.ts = c(lmfit$residuals, lm.predict - y.ts[-(1:npre)])
    
    # changepoint detection on residual time series
    cpt.fcst = cptForecast(errors=as.vector(r.ts), m=npre,
                           detector="PageCUSUM", forecastErrorType="Both")
    time.est = cpt.fcst@tau      # based on errors
    
    # retrieve the changepoints
    if (!is.infinite(time.est)) {
      nt = npre + time.est
      testStats = trend.stats(y.ctr=NULL, y.itv=r.ts[1:nt], 
                              nt=nt, monthly=TRUE)
      cpt.est = testStats[2]
    } else {
      cpt.est = Inf
    }
    
    return(c(cpt.est=cpt.est, time.est=time.est, delay=delay))
  })
  
  return(change.list)
})

save(simResult.mu1, simResult.mu2, 
     file="sim_result_distribution_change_Forecast.RData")

