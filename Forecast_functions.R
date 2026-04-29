## -----------------------------------------------------------
## Functions from "changepoint.forecast"
## -----------------------------------------------------------

## For use in the Python interactive interface

library(purrr)
library(tidyverse)
library(tibble)
library(dplyr)


## Ignore the class function, return a list

## The main function (wrapper of the CUSUM procedure)
cptForecast = function(errors,
                       m=ceiling(length(errors)/2),
                       detector='PageCUSUM',
                       gamma=0,
                       alpha=0.05){
  # Removed arguments (to simplify the code)
  # forecastErrorType='Raw',
  # critValue = 'Lookup',
  # samples=1000,
  # npts=500
  
  ## CUSUM for Raw and Squared errors (calculate both!)
  ans = cptSeqCUSUM(errors, m=m, detector=detector, gamma=gamma, sigma2=NULL,
                    alpha=alpha)
  ans2 = cptSeqCUSUM((errors-mean(errors[1:m]))^2, m=m, detector=detector, gamma=gamma,
                     sigma2=NULL, alpha=alpha)
  return(list(errors=errors, m=m, gamma=gamma, detector=detector,
              alpha=alpha, critValue=ans$critValue,
              errorsVar=ans$sigma2, errors2Var=ans2$sigma2,
              cusum=ans$cusum, cusum2=ans2$cusum,
              threshold=ans$threshold, threshold2=ans2$threshold,
              tau=ans$tau, tau2=ans2$tau))
  
}


## Functions for CUSUM (different test statistics)
# detector: c("CUSUM", "CUSUM1", "PageCUSUM", "PageCUSUM1")

cptSeqCUSUM = function(X,
                       m = ceiling(0.5*length(X)),
                       detector='PageCUSUM',
                       gamma=0,
                       sigma2=NULL,
                       alpha=0.05){
  
  N = length(X)
  n = N-m
  
  cusumValues = 0
  cusumA = 0
  cusumB = 0
  trainMean = mean(X[1:m])
  if(detector=='PageCUSUM'){
    for(i in 1:n){
      cusumA = max(cusumA + X[i+m] - trainMean, 0)
      cusumB = max(cusumB - X[i+m] + trainMean, 0)
      cusumValues[i+1] = max(cusumA, cusumB)
    }
  }else if(detector=='PageCUSUM1'){
    for(i in 1:n){
      cusumValues[i+1] = max(cusumValues[i] + X[m+i] - trainMean, 0)
    }
  }else if(detector=='CUSUM'){
    for(i in 1:n){
      cusumA = cusumA + X[i+m] - trainMean
      cusumValues[i+1] = abs(cusumA)
    }
  }else{
    for(i in 1:n){
      cusumValues[i+1] = cusumValues[i] + X[m+i] - trainMean
    }
  }
  
  cusumValues = cusumValues[-1]
  weightValues = purrr::map_dbl(1:n, ~weightFun(m=m, k=.x, gamma=gamma))
  # map_dbl, similar to lapply or sapply?
  
  # Only use the Lookup table for faster implementation
  critValue = critValLookup(gamma=gamma, alpha=alpha, detector=detector)
  
  if(is.null(sigma2)){
    sigma2 = stats::var(X[1:m])
  }
  thresholdValues = weightValues*critValue*sqrt(sigma2)
  thresholdExceeded = cusumValues > thresholdValues
  
  if(all(!thresholdExceeded)){
    tau = Inf
  }else{
    tau = min(which(thresholdExceeded==TRUE))
  }
  
  return(list(tau = tau,
              cusum = cusumValues,
              threshold = thresholdValues,
              sigma2 = sigma2,
              critValue = critValue))
}


## The weighting function
weightFun = function(m, k, gamma=0){
  if(!((gamma>=0)&&(gamma<0.5))){
    stop('gamma must be between 0 and 1/2')
  }
  return(sqrt(m)*(1 + k/m)*((k/(k+m))^gamma))
}


## Functions to read the critical values from the Lookup table
critValLookup = function(gamma=0, alpha=0.05, detector='PageCUSUM'){
  if(!(gamma %in% critValTable$Gamma)){
    stop(paste0(c('gamma value not stored in table. Either simulate a critical value using function simCritVal or try one of the following gamma values:\n',unique(critValTable$Gamma)), collapse=' '))
  }
  if(!(alpha %in% critValTable$Alpha)){
    stop(paste0(c('alpha value not stored in table. Either simulate a critical value using function simCritVal or try one of the following alpha values:\n',unique(critValTable$Alpha)), collapse=' '))
  }
  if(!(detector %in% critValTable$Detector)){
    stop(paste0(c('detector not in table. Either simulate a critical value using function simCritVal or try one of the following detectors:\n',unique(critValTable$Detector)), collapse=' '))
  }
  ans = dplyr::filter(critValTable, .data$Detector==detector, .data$Gamma==gamma, .data$Alpha==alpha)
  return(ans$CritVal)
}


## The look-up table (pre-run and saved locally)
critValTable = readRDS(file="CritValTable.rds")



## ---------------------------------------------------------------------------
## Test the reduced functions

# simulate some data
npre = 24                          # 2 years pre-intervention
npost.vec = seq(24, 120, by=12)    # from 2 to 10 years after intervention

npost.max = npost.vec[length(npost.vec)]
ntt = npre + npost.max

trend.control = 0
trend.increase = round(10 * c(0.05, seq(0.1, 1, by=0.1)) / 120, 4) 
trend.interv = trend.control + trend.increase[3]

# delay in effect, up to 20 time points
delay.set = 1:20
delay = sample(delay.set, size=1)

npre.delay = npre + delay
npost.delay = npost.max - delay
sim.ts = CIsim(seed=sample(1:50, size=1), 
               npre=npre.delay, npost=npost.delay, level=10, 
               trend=c(trend.control, trend.interv), sigma=0.1)

y.ts = sim.ts$y.itv  - sim.ts$y.ctr
Xmat = data.frame(x1=rep(1, ntt), x2=1:ntt)
lmfit = lm(y.ts[1:npre] ~ -1 + x1 + x2, data=Xmat[1:npre, ])
lm.predict = predict(lmfit, newdata=Xmat[-(1:npre), ])
r.ts = c(lmfit$residuals, lm.predict - y.ts[-(1:npre)])

# changepoint detection on residual time series
cpt.fcst = cptForecast(errors=as.vector(r.ts), m=npre, detector="PageCUSUM")
time.est = cpt.fcst$tau    



