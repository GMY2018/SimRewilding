## ------------------------------------------------
## Example 1: trend change in time series data
## ------------------------------------------------

## CAN IGNOR THIS SCRIPT FOR NOW AS THERE ARE SOME PROBLEMS. 


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


## ---------------------------
## Run the simulation study
## ---------------------------
## Simulation set up
npre = 24                        # 2 years pre-intervention
npost.vec = seq(24, 120, by=3)   # from 2 to 10 years after intervention

npost.max = npost.vec[length(npost.vec)]
ntt = npre + npost.max

level = 10 
trend.control = 0.005
# trend increment, from 1% to 50% over 10 years
trend.increase = round(level * c(0.01, 0.05, seq(0.1, 0.5, by=0.05)) / 120, 4) 
# trend.increase = round(level * c(0.01, 0.05, seq(0.1, 1, by=0.1)) / 120, 4) 

sigma = 0.05            # noise level
phi = c(0.5, 0.8)[2]    # auto-correlation strength

# delay in taking effect
delay.set = 1:20

## Simulation: trend with AR(1) noises
# using ARIMA(p, d, q) model (also bad result...)
simN = 1000
simResult.ar = lapply(1:simN, FUN=function(s) {
  
  change.list = sapply(1:length(trend.increase), FUN=function(m) {
    
    # simulate time series with different scales of trend
    trend.interv = trend.increase[m]
    
    delay = sample(delay.set, size=1)
    # delay = 0
    npre.delay = npre + delay
    npost.delay = npost.max - delay
    sim.ts = CIsim.ar(seed=s*length(trend.increase)+m, 
                      npre=npre.delay, npost=npost.delay, level=level, 
                      trend=c(trend.control, trend.interv), 
                      phi, sigma=sigma)
    
    ## consider two types of model
    y.ts = sim.ts$y.itv
    
    # (1) with time as covariate
    # Xmat = cbind(rep(1, ntt), 1:ntt)
    # armafit = try(arima(y.ts[1:npre], xreg=Xmat[1:npre, ], order=c(1, 0, 0),
    #                     include.mean=FALSE), silent=TRUE)
    # if (is.list(armafit)) {
    #   # h-step ahead prediction
    #   arima.predict = predict(armafit, newxreg=Xmat[-(1:npre), ],
    #                           n.ahead=npost.max)
    #   r.ts = c(armafit$residuals, arima.predict$pred - y.ts[-(1:npre)])
    # 
    #   # changepoint detection on residual time series
    #   cpt.fcst = cptForecast(errors=as.vector(r.ts), m=npre,
    #                          detector="PageCUSUM", forecastErrorType="Both")
    #   cpt.est = cpt.fcst@tau     # based on errors
    #   # cpt.est = cpt.fcst@tau2    # based on squared errors
    # } else { cpt.est = NULL }
    
    # (2) without time covariate, using 1st difference
    armafit = try(arima(y.ts[1:npre], order=c(1, 1, 0), 
                        include.mean=TRUE), silent=TRUE)
    if (is.list(armafit)) {
      # h-step ahead prediction
      arima.predict = predict(armafit, n.ahead=npost.max)
      r.ts = c(armafit$residuals, arima.predict$pred - y.ts[-(1:npre)])
      
      # changepoint detection on residual time series
      cpt.fcst = cptForecast(errors=as.vector(r.ts), m=npre,
                             detector="PageCUSUM", forecastErrorType="Both")
      cpt.est = cpt.fcst@tau     # based on errors
      # cpt.est = cpt.fcst@tau2    # based on squared errors
    } else { cpt.est = NULL }

    return(list(cpt.est=cpt.est, delay=delay))
  })
  
  return(change.list)
})
# the returned value is the detection delay?

save(simResult.ar, file="sim_result_trend_change_Forecast.RData")
