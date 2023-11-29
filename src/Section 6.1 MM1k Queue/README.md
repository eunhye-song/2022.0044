This folder contains three R function codes and one R script (RunMM1kExpr.R) used to run the M/M/1/k queue example presented in the paper. The functions are described below:
* GG1k.LindleySteadyState.R: simulates one observation of time in system of a job in the G/G/1/k queue in (approximate) steady state. Both interarrival and service times are generated from the empirical distributions constructed from input data.
* MM1k.LindleySteadyState.Parametric.R: simulates one observation of time in system of a job in the M/M/1/k queue in (approximate) steady state. Both interarrival and service times are generated from exponential distributions. This function is needed for the Monte Carlo estimation of the true performance measure. 
* runMacro.shrinkage.MM1k.R: This function computes all bootstrap confidence intervals compared in the paper. 

