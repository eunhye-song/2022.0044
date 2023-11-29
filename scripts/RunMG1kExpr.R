# This script runs macro experiments for the M/G/1/k experiment in Song et al. (2023)
library(compiler)
library(resample)
library(pracma)
library(profvis)
library(boot)

source("runMacro.shrinkage.MG1k.R")
source("GG1k.Lindley.R") 
source("MM1k.Lindley.Parametric.R") 

# random number seed & number of macro replications
set.seed(1234)
macrorep = 10

# placeholder for the macro simulation results
result = NULL

# We test a mixture-of-betas service distribution
# Ghosh and Lam (2019): 0.3 * Beta(2, 6) + 0.7 * Beta(6, 2)
# service distribution parameters 
GLbetaMix_p = .3
GLalpha1a = 2
GLalpha2a = 6
GLalpha1b = 6
GLalpha2b = 2
bimod_GL_mean = 
  GLbetaMix_p*(GLalpha1a/(GLalpha1a+GLalpha2a)) + (1-GLbetaMix_p)*(GLalpha1b/(GLalpha1b+GLalpha2b))
list_GL = list(GLbetaMix_p,GLalpha1a,GLalpha2a,GLalpha1b,GLalpha2b) 

# arrival rate to result in rho = .7
lambda = .7*(1/bimod_GL_mean) 

# run conditions
k = 10 # system capacity for the M/G/1/k queue
alpha = 0.05 # nominal error rate for the confidence intervals
m = 20 # number of observed times in system per replication
n = rep(c(rep(100,2),rep(400,2),rep(1000,2),rep(4000,2))) # real-world sample sizes 
R0 = ceiling(n^1.1) # number of replications at the original ecdf of real-world data
R = ceiling(c(0.05,0.2)*n) # number of replications at each bootstrap sample
B = 1000 # bootstrap sample size

# the expected time in system estimated from 10^6 replications 
Expected_TIS = 1.17043 

# run eight experiment settings specified above
for(i in 1:8){
  result = rbind(result, as.data.frame(runMacro.shrinkage.MG1k(
    lambda,GLbetaMix_p,GLalpha1a,GLalpha2a,GLalpha1b,GLalpha2b,
    k,alpha,Expected_TIS,m,R0[i],R[i],macrorep,n[i],B)))
}


# save the result as a .Rdata file
filename <- sprintf(paste("MG1k.Rdata",sep=""))
save(result,file = filename) 
