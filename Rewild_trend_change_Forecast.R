## ------------------------------------------------
## Example 1: trend change in time series data
## ------------------------------------------------

# Using the changepoint.forecast package for online detection
# Detecting changepoint as a change in the forecast errors
# This can be applied to time series with auto-correlated noise

# library(changepoint.forecast)    # if calling the package
source("Forecast_functions.R")     # if using the essential functions only


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
## Simulation set up, assuming a higher sampling frequency
npre = 24                          # e.g., 24 months (2 years) pre-intervention
npost.vec = seq(24, 120, by=12)    # e.g., from 2 to 10 years after intervention

npost.max = npost.vec[length(npost.vec)]
ntt = npre + npost.max

level = 10 
trend.control = 0.005
# trend increment, from 1% to 50% over 10 years, or 5% to 100% in 10 years
trend.increase = round(level * c(0.05, seq(0.1, 1, by=0.1)) / 120, 4) 

sigma = 0.05     # noise levels
phi = 0.8        # auto-correlation strength

# delay in effect, up to 20 time points
delay.set = 1:20


## (1) Simulation: BA (before-after) + i.i.d. noise time series
# alpha = 0.95     # significance level
simN = 1000
simResult = lapply(1:simN, FUN=function(s) {
  
  change.list = sapply(1:length(trend.increase), FUN=function(m) {
    
    # simulate time series with different scales of trend
    trend.interv = trend.control + trend.increase[m]
    
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
    
    # changepoint detection on error time series
    # Using the changepoint.forecast package
    # cpt.fcst = cptForecast(errors=as.vector(r.ts), m=npre,   
    #                        detector="PageCUSUM", forecastErrorType="Both")
    # time.est = cpt.fcst@tau    

    # Using the essential functions from "Forecast_functions.R"
    cpt.fcst = cptForecast(errors=as.vector(r.ts), m=npre, detector="PageCUSUM")
    time.est = cpt.fcst$tau
    
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


## (2) Simulation: BA (before-after) + AR(1) noise time series
simN = 1000
simResult.ar = lapply(1:simN, FUN=function(s) {
  
  change.list = sapply(1:length(trend.increase), FUN=function(m) {
    
    # simulate time series with different scales of trend
    trend.interv = trend.control + trend.increase[m]
    
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
      
      # changepoint detection on error time series
      # Using the changepoint.forecast package
      # cpt.fcst = cptForecast(errors=as.vector(r.ts), m=npre,   
      #                        detector="PageCUSUM", forecastErrorType="Both")
      # time.est = cpt.fcst@tau    

      # Using the essential functions from "Forecast_functions.R"
      cpt.fcst = cptForecast(errors=as.vector(r.ts), m=npre, detector="PageCUSUM")
      time.est = cpt.fcst$tau
      
      # retrieve the changepoint
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

save(simResult, simResult.ar, file="sim_result_trend_change_Forecast.RData")

