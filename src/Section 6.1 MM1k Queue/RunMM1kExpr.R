# This script runs macro experiments for the M/M/1/k experiment in Song et al. (2023)
library(compiler)
library(resample)
library(pracma)
library(profvis)
library(boot)

source("runMacro.shrinkage.MM1k.R")
source("GG1k.LindleySteadyState.R")

# random number seed & number of macro replications
set.seed(1234)
macrorep = 10

# placeholder for the macro simulation results
result = NULL

# run conditions 
lambda = 0.7 # arrival rate 
mu = 1.0 # service rate 
k = 10 # system capacity of the M/M/1/k queue
alpha = 0.05 # nominal error rate for the confidence interval
n = rep(c(rep(100,3),rep(400,3),rep(1000,3),rep(4000,3))) # real-world sample sizes
B = 1000 # number of bootstrap sample size
R0 = ceiling(n^1.1) # number of replications to run at the original ecdf
R = ceiling(c(0.05,0.2,0.5)*n) # number of replications per bootstrap sample

# run twelve experiment settings specified above
for(i in 1:12){
  result = rbind(result, as.data.frame(runMacro.shrinkage.MM1k(lambda,mu,k,alpha,R0[i],R[i],macrorep,n[i],B)))
}

# run conditions 
lambda = 0.9 # arrival rate 
mu = 1.0 # service rate 
k = 10 # system capacity of the M/M/1/k queue
alpha = 0.05 # nominal error rate for the confidence interval
m = 1 # number of observation of the waiting time per replication
n = rep(c(rep(100,3),rep(400,3),rep(1000,3),rep(4000,3))) # real-world sample sizes
B = 1000 # number of bootstrap sample size
R0 = ceiling(n^1.1) # number of replications to run at the original ecdf
R = ceiling(c(0.05,0.2,0.5)*n) # number of replications per bootstrap sample

# run twelve experiment settings specified above
for(i in 1:12){
  result = rbind(result, as.data.frame(runMacro.shrinkage.MM1k(lambda,mu,k,alpha,R0[i],R[i],macrorep,n[i],B)))
}

# save the result as a .Rdata file
filename <- sprintf(paste("MM1k.Rdata",sep=""))
save(result,file = filename) 
