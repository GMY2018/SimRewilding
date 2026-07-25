# SimRewilding

This repository contains the R code for running the simulation study to demonstrate the use of changepoint detection methods on field sampling data from rewilding studies. It accompanies the interacitive tool <b>TRACE</b> (https://github.com/M-DICE/trace) which helps ecologists answer the question: how quickly does a restoration intervention produce a measurable effect?

## The simulation study

We consider typical type of experimental design, that is, paired intervention and control time series, with observations from both before and after the management. <b>Two types of changes are considered, (1) the trend change and (2) the distribution change.</b> For distribution change, we consider a change in the distribution parameter. 

<b>For each type of change, we simulate time series of two different frequencies, with one being the 1/12 of the other, mimicing the difference between a monthly and an annual samplling frequencies</b>. Each time series has one changepoint, and it usually happens after the intervention point. The length of post-intervention period is set to 10 time points for the low frequency time series and 120 time points for the high frequency time series. For the high frequency time series, we consider both the independent and auto-correlated noises; whereas for the lower frequency time series, only the independent noises are considered. 

<b>We implement two types of changepoint detection approaches,</b> (a) the off-line analysis (or retrospective analysis) which detects at most one changepoint in the time series, and (b) the online analysis (or the streamline analysis) which scans through the time series from the beginning and stops once the first changepoint is detected. <b>We use the at most one change (AMOC) method for the off-line analysis and the forecast error-based (Forecast) method [1] and the Bayesian online changepoint detection (BOCPD) method [2] for online detection.</b> Specifically, we use the Forecast method to detect the change in high frequency time seires and the BOCPD method to detect the change in low frequency time seires . 

## What's in the repository

This repository contains all scripts required to run the above simulation. It is recommended that one download all the files from the repository, so that one will not miss any files sourced by the simulation scripts. 

<ul> 
  <li>The scripts named "Rewild_trend_change_XXXX.R" are for simulating the distribution change scenarios, where XXXX indicates the method used for changepoint detection.</li>
  <li>The scripts named "Rewild_distribution_change_XXXX.R" are for simulating the distribution change scenarios, where XXXX indicates the method used for changepoint detection.</li>
  <li>The scripts named "XXXX_functions.R" contain functions for implementing the two online chagnepoint detection methods.</li>
  <li>The script named "Simulation_figures_and_tables.R" is for producing summary tables and figures of the simulation result.</li>
  <li>The data file "CritValTable.rds" contains the table of critical values for the Forecast method.</li>
</ul>

## How to run the simulation

First <b>install R</b> (https://cran.r-project.org/) in one hasn't done so. A version from 4.3.1 onwards is recommended. 

Then install missing packages if required. In particular, one will also need to install R package <b>transport</b> (https://cran.r-project.org/web/packages/transport/index.html), which is used to calculate the Wasserstein distance for the distribution change scenarios. We pulled the essential functions from the R package <b>changepoint.forecast</b> (https://github.com/grundy95/changepoint.forecast) into the script "Forecast_functions". One can also install the package and call the functions from the package.

Finally, run the code on the simulation scripts (that is, those named "Rewild_XXXX_change_XXXX.R).

Remember to change the working directory to the folder containing the simulation scripts.

### References

[1]. Grundy, T., R.Killick, and I.Svetunkov. 2026. “Online Detection of Forecast Model Inadequacies Using Forecast Errors.” Journal of Time Series Analysis47, no. 3: 715–726. https://doi.org/10.1111/jtsa.12843.

[2]. Paul Fearnhead, Zhen Liu, On-Line Inference for Multiple Changepoint Problems, Journal of the Royal Statistical Society Series B: Statistical Methodology, Volume 69, Issue 4, September 2007, Pages 589–605, https://doi.org/10.1111/j.1467-9868.2007.00601.x
