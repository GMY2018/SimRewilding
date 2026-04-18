## ------------------------------------------------
## Example 1: trend change in time series data
## ------------------------------------------------

# Using the changepoint.forecast package for online detection
# Detecting changepoint as a change in the forecast errors
# This can be applied to time series with auto-correlated noise

.libPaths(c("~/R/library", .libPaths()))
library(changepoint.forecast)


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

CIsim.ar = function(seed=10, npre, npost, level, trend, phi, sigma) {
  set.seed(seed)
  nt = npre + npost
  
  trend.control = trend[1]
  trend.interv = trend[2]
  # with auto-correlation
  y.control.ar = level + trend.control * 1:nt + arima.sim(model=list(ar=phi), n=nt, sd=sigma)
  y.interv.ar = c(level + trend.control * 1:npre + arima.sim(model=list(ar=phi), n=npre, sd=sigma),
                  level + trend.control * npre + trend.interv * 1:npost + 
                    arima.sim(model=list(ar=phi), n=npost, sd=sigma))
  
  return(list(y.ctr=y.control.ar, y.itv=y.interv.ar))
}


trend.stats = function(y.ctr=NULL, y.itv, nt, monthly=TRUE) {
  # If we only have the intervention time series, set y.ctr to NULL
  
  if (monthly == TRUE) {
    mint = 24  
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
## Simulation set up
npre = 24                          # 2 years pre-intervention
npost.vec = seq(24, 120, by=12)    # from 2 to 10 years after intervention

npost.max = npost.vec[length(npost.vec)]
ntt = npre + npost.max

level = 10 
trend.control = 0.005
# trend increment, from 1% to 50% over 10 years, or 5% to 100% in 10 years
# trend.increase = round(level * c(0.01, 0.05, seq(0.1, 0.5, by=0.05)) / 120, 4) 
trend.increase = round(level * c(0.05, seq(0.1, 1, by=0.1)) / 120, 4) 

sigma.vec = c(0.05, 0.1)            # noise levels
phi.vec = c(0.5, 0.8)               # auto-correlation strength

# delay in effect, up to 20 time points
delay.set = 1:20


## Simulation: trend with AR(1) noises
# using ARIMA(p, d, q) model (bad result...)
alpha = 0.95     # significance level
phi = phi.vec[2]
simN = 1000

## small noise
sigma = sigma.vec[1]
simResult.small.v2 = lapply(1:simN, FUN=function(s) {
  
  change.list = sapply(1:length(trend.increase), FUN=function(m) {
    
    # simulate time series with different scales of trend
    trend.interv = trend.increase[m]
    
    # if using delay in effect, otherwise, set delay to 0
    delay = sample(delay.set, size=1)
    # delay = 0
    npre.delay = npre + delay
    npost.delay = npost.max - delay
    sim.ts = CIsim(seed=s*length(trend.increase)+m, 
                   npre=npre.delay, npost=npost.delay, level=level, 
                   trend=c(trend.control, trend.interv), sigma=sigma)
    
    ## Train the model on pre-data, and then predict the post-data
    y.ts = sim.ts$y.itv  # - sim.ts$y.ctr
    
    # (1) with time as covariate
    Xmat = data.frame(x1=rep(1, ntt), x2=1:ntt)
    lmfit = lm(y.ts[1:npre] ~ -1 + x1 + x2, data=Xmat[1:npre, ])
    # h-step ahead prediction
    lm.predict = predict(lmfit, newdata=Xmat[-(1:npre), ])
    r.ts = c(lmfit$residuals, lm.predict - y.ts[-(1:npre)])
    
    # changepoint detection on residual time series
    cpt.fcst = cptForecast(errors=as.vector(r.ts), m=npre,
                           detector="PageCUSUM", forecastErrorType="Both")
    time.est = cpt.fcst@tau      # based on errors
    # time.est = cpt.fcst@tau2     # based on squared errors
    
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

# AR1 time series
simResult.ar.small = lapply(1:simN, FUN=function(s) {
  
  change.list = sapply(1:length(trend.increase), FUN=function(m) {
    
    # simulate time series with different scales of trend
    trend.interv = trend.increase[m]
    
    # if using delay in effect, otherwise, set delay to 0
    delay = sample(delay.set, size=1)
    # delay = 0
    npre.delay = npre + delay
    npost.delay = npost.max - delay
    sim.ts = CIsim.ar(seed=s*length(trend.increase)+m, 
                      npre=npre.delay, npost=npost.delay, level=level, 
                      trend=c(trend.control, trend.interv), 
                      phi, sigma=sigma)
    
    ## Train the model on pre-data, and then predict the post-data
    y.ts = sim.ts$y.itv    # - sim.ts$y.ctr
    
    # (1) with time as covariate
    Xmat = cbind(rep(1, ntt), 1:ntt)
    armafit = try(arima(y.ts[1:npre], xreg=Xmat[1:npre, ], order=c(1, 0, 0),
                        include.mean=FALSE), silent=TRUE)
    if (is.list(armafit)) {
      # h-step ahead prediction
      arima.predict = predict(armafit, newxreg=Xmat[-(1:npre), ],
                              n.ahead=npost.max)
      r.ts = c(armafit$residuals, arima.predict$pred - y.ts[-(1:npre)])
      
      # changepoint detection on residual time series
      cpt.fcst = cptForecast(errors=as.vector(r.ts), m=npre,
                             detector="PageCUSUM", forecastErrorType="Both")
      time.est = cpt.fcst@tau      # based on errors
      # time.est = cpt.fcst@tau2     # based on squared errors
      
      # retrieve changepoints
      if (!is.infinite(time.est)) {
        nt = npre + time.est
        testStats = trend.stats(y.ctr=NULL, y.itv=r.ts[1:nt], 
                                nt=nt, monthly=TRUE)
        cpt.est = testStats[2]
      } else { cpt.est = Inf }
    } else { 
      time.est = Inf
      cpt.est = Inf
    }
    
    return(c(cpt.est=cpt.est, time.est=time.est, delay=delay))
  })
  
  return(change.list)
})

save(simResult.small.v1, simResult.small.v2, simResult.ar.small, 
     file="sim_result_trend_change_Forecast.RData")


## large noise
sigma = sigma.vec[2]
simResult.large.v2 = lapply(1:simN, FUN=function(s) {
  
  change.list = sapply(1:length(trend.increase), FUN=function(m) {
    
    # simulate time series with different scales of trend
    trend.interv = trend.increase[m]
    
    # if using delay in effect, otherwise, set delay to 0
    delay = sample(delay.set, size=1)
    # delay = 0
    npre.delay = npre + delay
    npost.delay = npost.max - delay
    sim.ts = CIsim(seed=s*length(trend.increase)+m, 
                   npre=npre.delay, npost=npost.delay, level=level, 
                   trend=c(trend.control, trend.interv), sigma=sigma)
    
    ## Train the model on pre-data, and then predict the post-data
    y.ts = sim.ts$y.itv  # - sim.ts$y.ctr
    
    # (1) with time as covariate
    Xmat = data.frame(x1=rep(1, ntt), x2=1:ntt)
    lmfit = lm(y.ts[1:npre] ~ -1 + x1 + x2, data=Xmat[1:npre, ])
    # h-step ahead prediction
    lm.predict = predict(lmfit, newdata=Xmat[-(1:npre), ])
    r.ts = c(lmfit$residuals, lm.predict - y.ts[-(1:npre)])
    
    # changepoint detection on residual time series
    cpt.fcst = cptForecast(errors=as.vector(r.ts), m=npre,
                           detector="PageCUSUM", forecastErrorType="Both")
    time.est = cpt.fcst@tau      # based on errors
    # time.est = cpt.fcst@tau2     # based on squared errors
    
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

# AR1 time series
simResult.ar.large = lapply(1:simN, FUN=function(s) {
  
  change.list = sapply(1:length(trend.increase), FUN=function(m) {
    
    # simulate time series with different scales of trend
    trend.interv = trend.increase[m]
    
    # if using delay in effect, otherwise, set delay to 0
    delay = sample(delay.set, size=1)
    # delay = 0
    npre.delay = npre + delay
    npost.delay = npost.max - delay
    sim.ts = CIsim.ar(seed=s*length(trend.increase)+m, 
                      npre=npre.delay, npost=npost.delay, level=level, 
                      trend=c(trend.control, trend.interv), 
                      phi, sigma=sigma)
    
    ## Train the model on pre-data, and then predict the post-data
    y.ts = sim.ts$y.itv    # - sim.ts$y.ctr
    
    # (1) with time as covariate
    Xmat = cbind(rep(1, ntt), 1:ntt)
    armafit = try(arima(y.ts[1:npre], xreg=Xmat[1:npre, ], order=c(1, 0, 0),
                        include.mean=FALSE), silent=TRUE)
    if (is.list(armafit)) {
      # h-step ahead prediction
      arima.predict = predict(armafit, newxreg=Xmat[-(1:npre), ],
                              n.ahead=npost.max)
      r.ts = c(armafit$residuals, arima.predict$pred - y.ts[-(1:npre)])
      
      # changepoint detection on residual time series
      cpt.fcst = cptForecast(errors=as.vector(r.ts), m=npre,
                             detector="PageCUSUM", forecastErrorType="Both")
      time.est = cpt.fcst@tau      # based on errors
      # time.est = cpt.fcst@tau2     # based on squared errors
      
      # retrieve changepoints
      if (!is.infinite(time.est)) {
        nt = npre + time.est
        testStats = trend.stats(y.ctr=NULL, y.itv=r.ts[1:nt], 
                                nt=nt, monthly=TRUE)
        cpt.est = testStats[2]
      } else { cpt.est = Inf }
    } else { 
      time.est = Inf 
      cpt.est = Inf
    }
    
    return(c(cpt.est=cpt.est, time.est=time.est, delay=delay))
  })
  
  return(change.list)
})

save(simResult.small.v1, simResult.small.v2, simResult.ar.small, 
     simResult.large.v1, simResult.large.v2, simResult.ar.large, 
     file="sim_result_trend_change_Forecast.RData")


## --------------------------------------------------------------------
load("sim_result_trend_change_Forecast.RData")

# i.i.d data
Cpt.mat =  sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase), FUN=function(j) {
    simResult.small.v1[[i]][1, j]   # the estimated changepoint
  })
})

Time.mat =  sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase), FUN=function(j) {
    simResult.small.v1[[i]][2, j]   # the detection delay
  })
})

Delay.mat =  sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase), FUN=function(j) {
    simResult.small.v1[[i]][3, j]   # the delay in effect
  })
})

# auto-ccorrelated data
Cpt.mat =  sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase), FUN=function(j) {
    simResult.ar.small[[i]][1, j]   # the estimated changepoint
  })
})
Time.mat =  sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase), FUN=function(j) {
    simResult.ar.small[[i]][2, j]   # the detection delay
  })
})

Delay.mat =  sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase), FUN=function(j) {
    simResult.ar.small[[i]][3, j]   # the delay in effect
  })
})

Cpt.true = Delay.mat + npre
correctMat = abs(Cpt.mat - Cpt.true) <= 2
Detect.mat = apply(correctMat, MARGIN=1, FUN=sum, na.rm=TRUE) / simN

delayMat = Time.mat + npre - Cpt.true
Delay.median = apply(delayMat, MARGIN=1, FUN=median, na.rm=TRUE)
Delay.mean = apply(delayMat, MARGIN=1, FUN=function(x) {
  x.fin = x[!is.infinite(x)]   # remove the infinite values (i.e., no detection)
  mean(x.fin, na.rm=TRUE)
})
Delay.std = apply(delayMat, MARGIN=1, FUN=function(x) {
  x.fin = x[!is.infinite(x)]
  sd(x.fin, na.rm=TRUE)
})
NoDetection = apply(delayMat, MARGIN=1, FUN=function(x) {
  sum(is.infinite(x))
})

summary.mat = rbind(Detect.mat, Delay.median, Delay.mean, Delay.std, NoDetection)
colnames(summary.mat) = paste0("trend_", trend.increase)
rownames(summary.mat) = c("True_pos", "median_delay", "mean_delay", "std_delay",
                          "no_detection")
summary.mat

# the result for the i.i.d. case makes senes
# the result for the auto-correlated case doesn't...

