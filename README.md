# SimRewilding

Simulation study to demonstrate the use of changepoint detection methods on field sampling data from rewilding studies.

### This is a place for sharing the simulation code

We simulated time series with trend change an time series with distribution change. 

We use the at most one change (AMOC) method for detecting one changepoint offline for both types of time series. 

We use Bayesian online changepoint detection (BOCPD) for online detection of trend change in time seires with i.i.d. noises, and we use changepoint.forecast package for online detection of trend change (based on forecast errors) in time seires with auto-correlated noises.
