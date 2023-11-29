This folder contains three R function codes and one R script (RunMG1kExpr.R) used to run the M/M/1/k queue example presented in the paper. The functions are described below:
* GG1k.Lindley.R: simulates the average time in system of jobs in the G/G/1/k queue starting from the empty and idle state. Both interarrival and service times are generated from the empirical distributions constructed from input data.
* MM1k.Lindley.Parametric.R: simulates the average time in system of jobs in the M/M/1/k queue starting from the empty and idle state. Both interarrival and service times are generated from exponential distributions.
* runMacro.shrinkage.MG1k.R: This function computes all bootstrap confidence intervals compared in the paper. 
 
