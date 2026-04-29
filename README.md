# SimRewilding

Simulation study to demonstrate the use of changepoint detection methods on field sampling data from rewilding studies.

### This is a place for sharing the simulation code

We simulated time series with trend change and time series of distributions with distribution parameter change.

We use the at most one change (AMOC) method for detecting a single changepoint offline for both types of changes with higher and lower sampling frequencies respectively. 

We use the forecast error method [1] for online detection of a trend change (i.e., regression coefficient change) in higher frequency time seires with i.i.d. or auto-correlated noises. We use the Bayesian online changepoint detection (BOCPD) method [2] for online detection of trend change in time seires with i.i.d. noises when the sampling frequency is low. 


### References

[1]. Grundy, T., R.Killick, and I.Svetunkov. 2026. “Online Detection of Forecast Model Inadequacies Using Forecast Errors.” Journal of Time Series Analysis47, no. 3: 715–726. https://doi.org/10.1111/jtsa.12843.

[2]. Paul Fearnhead, Zhen Liu, On-Line Inference for Multiple Changepoint Problems, Journal of the Royal Statistical Society Series B: Statistical Methodology, Volume 69, Issue 4, September 2007, Pages 589–605, https://doi.org/10.1111/j.1467-9868.2007.00601.x
