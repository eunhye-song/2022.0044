# DirectIUJackson_4BROAR.R -- adds bootstrap lower CI for c
# Updated by RRB July 2023


# Purpose: direct resampling input uncertainty quantification
# Specialized for large-scale Jackson network-like (capacitated) example
# Routines general where possible, but treatment of service time
#  distributions and routing distributions treated separately, and specially

# Capacitated Jackson network simulated via Python code

# Uses R packages for latin hypercube sampling for MM doe (lhs) and 
#  spatial correlation with nugget (mlegp)
# ROAR Collab Version for R main program =============
library(compiler)
library(resample)
library(pracma)
#

library(reticulate) # This is the library you need to run python on R
library(geometry)   # Uses R package 'geometry' for cart2bary Cartesian to Barycentric
# conversion for routing probability vectors
# (bary2cart not used - this was easy to code directly in B2C.R)
library(boot)       # bootstrap lower CI for shrinkage c
library(tictoc)     # to assess computational effort to compute CIs

# Set to LEcuyer RNG
RNGkind(kind = "L'Ecuyer-CMRG")
# set different seed for each processor node
# AI <- Sys.getenv("SLURM_ARRAY_TASK_ID") # used for parallel runs on ROAR cluster
AI = "99" #TEMPORARY UNTIL RUN ON ROAR
NAI <- as.numeric(AI)
set.seed(NAI)


# Other functions for code clarity or network structure
source("BSgenerator2.R") # generate parameter values from resampled data
source("setPijk2.R")     # set values of network routing probabilities - setPijk2 is UNBALANCED
# setPijk is BALANCED
source("inputData.R")    # routine to generate 'original' sample of interarrival 
#  and routing data
#  replaces GG1ksim used by Wei Xie
source("getY.R")         # Jackson-like simulation for R0 run estimate of mean wait 


# Set key parameters:
#===================
#=================== 
#  input data characteristics (and compute from these either 
#  analytically or visa simulation the true expected value for coverage), 
#  ellipsoid %, sim run length, warmup and reps, CI significance level, 
#  number of bootstrap resamples, number of macro replications to check coverage

# overall timing
tic()

# Set input data characteristics (number of arrival/stage 2 nodes is numArv)
#===============================
numSamples = 100         # assume number of samples for each distribution the same 
#    use either 100 (with reps 5, 20, 50)
#    or 400 (with reps 20, 80, 200) (well behaved if large)
numArv = 4               # number of arrival nodes (assume exponential arrivals for now)
nomArate = .7            # set nominal arrival rate
rateArv = rep(nomArate,numArv) # set arrival distribution rates (exponential rates)
# move arrival rates away from balanced 
delRate = .5*(1:numArv)/sum(1:numArv)
rateArv = numArv*nomArate*delRate
# adjust to have mean of nomArate
rateArv = (nomArate-mean(rateArv)) + rateArv
# set initial mean arrival time for Y00 true Y calculation
mArvTime = 1/rateArv


deltaP = 1/(5*numArv)         # deltaP passed to setPijk2 it allows shifts in equally likely routing 
# deltaP = 1/5numArv gives shifts of +/-20% from equal


PijkDim = c(2,numArv,numArv)  # matrix of dimensions of: (# layers - 1), (numnodes layer 1,)    
#  (numnodes layer 2), where PijkRoute = routing probability from 
#  jth node in layer i to kth node in layer i+1
#  NOTE: when i = 2, k takes only three values, not numArv

# Set unbalanced Pijk parameters by calling 'setPijk2'
PijkRoute = array(0, dim = PijkDim)


# Set unbalanced routing to use balanced parameters, use setPijk instead of setPijk2
PijkRoute = setPijk2(PijkRoute,deltaP) 


# Set other key parameters (method, runlength and warmup set in runSim and getY)
#=========================
CImethod = "Direct"       # output file different for each: "Direct" or "Metamodel" 
nBootCI = 220             # number of bootstrap runs to compute CI - should be 1000 for 10, 220 for 4
mRuns = 10                # 1000 macro replications on 100 processors
                          # NOTE: 95% CI for nom. coverage of 90% at +/- 1% requires 3500 macroReps
                          # 95% CI for nom. coverage of 90% with 1000 reps is +/- 1.96%
                          # NOTE - these are conservative since less variance if nom. coverage near 95%
simReps = 10              # number of replications for sample mean, sd
                          # R0simReps is default number of initial replications R0
R0simReps = ceiling(numSamples^1.1)  
trueY = 0                 # must be determined by long simulation using bigSimReps
warmup = 100              # first customer waiting times deleted (set to 100 in real run)
runlength = 10            # number of customers per replication (excluding those deleted)
alp = 0.05                # setting for confidence interval percentage (95% if alp = 0.05)

# First, establish trueY with high replication simulation
# Get Y00 mean wait from original sample using 3,000,000 reps
#==================================================================
# output of getY: mean, variance(over reps),#reps
# trueY = 0 # 95%CI is  for runlength = 50
#trueY = 5.018 # SD .001218 CI "5.015","5.02" for runlength = 5
trueY = 5.040 # SD 0.001053 CI "5.037","5.042" for runlength = 10
#trueY = 5.151 # SD 0.0008215 CI "5.149","5.152" for runlength = 50
#trueY = 5.306 # SD 0.0020149 CI "5.303","5.311" for runlength = 400

#==================================================================

#print("initialization time:")
#toc()
# Begin macro loop:
#=================
#=================
# Initialize some vectors and arrays
#=================
macroResults = matrix(0,mRuns,22) # matrix for holding macro results of CIs, their size and coverage
runOutput = NULL
BootYbar = NULL
BootYvar = NULL

# mRun = 1 # use for testing before implementing macro loop
for(mRun in 1:mRuns){
  
  py_run_file("VBASim.py")
  py_run_file("JN2c.py")    # load the Jackson Network simulator (python)
 
# Generate original sample: (assume all sample sizes same)
#=========================
#=========================

# Sample based on the input distribution parameter values

Idata = inputData(numSamples,numArv,rateArv,PijkDim,PijkRoute)
                         # structure of "Idata"
                         # list(arrivals[numArv,numSamples],PijkHat[i,j,k])
                         # PijkHat converts RouteCount to sample estimate of probabilities

mArvTime = rowMeans(Idata[[1]])# the sample average interarrival times - not used
PijkHat = Idata[[2]]           # the sample estimates of routing probabilities - not used
#tic()
# Make initial simulations with sampled data
# reps and runlength R0simReps = ceiling(numSamples^1.1) runlength = 100
Y0results = getY(R0simReps,warmup, runlength, numArv, mArvTime, PijkHat) 
# output of getY: mean, variance(over reps),#reps
Y0 = Y0results[1]
#print("Y0 computation (158 reps rather than 5):")
#toc()
# Make Direct Bootstrap Runs:
#=========================== 
#===========================

# Resample mArvTime, PijkHat all B resamples in newSamples

newSamples = BSgenerator2(numSamples,numArv,PijkDim,nBootCI,Idata)
for(nBoot in 1:nBootCI){
  py_run_file("VBASim.py")
  py_run_file("JN2c.py")    
  # Extract bth bootstrap mArvTime PijkHat from newSamples
  mArvTime = newSamples[[1]][nBoot,]
  PijkHat = newSamples[[2]][nBoot,,,]
  # Make simulation run
  #tic()
  Yresults = getY(simReps,warmup, runlength, numArv, mArvTime, PijkHat)
  #toc()
  BootYbar[nBoot] = Yresults[1]
  BootYvar[nBoot] = Yresults[2]
  #print(nBoot)
}


# Basic Bootstrap CI
Lbasic = Y0 - quantile(BootYbar - Y0,1-alp/2)
Ubasic = Y0 - quantile(BootYbar - Y0, alp/2)
Coverbasic = as.integer((trueY<=Ubasic)&&(trueY>=Lbasic))

# Percentile Bootstrap CI
Lpercentile = quantile(BootYbar,alp/2) 
Upercentile = quantile(BootYbar,1-alp/2)
Coverpercentile = as.integer((trueY<=Upercentile)&&(trueY>=Lpercentile))

# Normal Bootstrap CI
Lnormal = Y0 - qnorm(1-alp/2)*sd(BootYbar) 
Unormal = Y0 + qnorm(1-alp/2)*sd(BootYbar) 
Covernormal = as.integer((trueY<=Unormal)&&(trueY>=Lnormal))

# Shrinkage basic CI
#Vhat = var(BootYbar)
#What = sum(BootYvar)/nBootCI
#chat = 1 - sqrt(max(0,(nBootCI/(nBootCI-1))*((Vhat-(What/simReps))/Vhat)))
# next two lines original version
#SSB = (nBootCI-1)*var(BootYbar)
#chat = 1-sqrt(max(0,1-sum(BootYvar)/simReps/SSB)) # here, 1 is used instead of B/(B-1)

#conservative chat based on bootstrap 95% lower CI
nBCI <<- nBootCI
Reps <<- simReps
Yvals = cbind(BootYbar,BootYvar)
c_val = function(Yvals,inds){
  SSB = (nBCI-1)*var(Yvals[inds,1])
  chat = 1-sqrt(max(0,1-sum(Yvals[inds,2])/Reps/SSB)) # here, 1 is used instead of B/(B-1)
}
boot.out = boot(Yvals,c_val,1000)
chat = boot.ci(boot.out,conf=.95,type="perc")$percent[4]
Y_shrinkage = chat*mean(BootYbar)+(1-chat)*BootYbar
Lshrinkage = Y0 - quantile(Y_shrinkage - Y0,1-alp/2)
Ushrinkage = Y0 - quantile(Y_shrinkage - Y0,alp/2)
Covershrinkage = as.integer((trueY<=Ushrinkage)&&(trueY>=Lshrinkage))

# Shrinkage percentile CI
Lshrinkage_p = quantile(Y_shrinkage,alp/2)
Ushrinkage_p = quantile(Y_shrinkage,1-alp/2)
Covershrinkage_p = as.integer((trueY<=Ushrinkage_p)&&(trueY>=Lshrinkage_p))

# shrinkage scale basic CI
Lshrinkagescale = Y0 - (1-chat)*quantile(BootYbar - Y0,1-alp/2)
Ushrinkagescale = Y0 - (1-chat)*quantile(BootYbar - Y0,alp/2)
Covershrinkagescale = as.integer((trueY<=Ushrinkagescale)&&(trueY>=Lshrinkagescale))

# shrinkage scale percentile CI
Lshrinkagescale_p = Y0 + (1-chat)*quantile(BootYbar - Y0,alp/2)
Ushrinkagescale_p = Y0 + (1-chat)*quantile(BootYbar - Y0,1-alp/2)
Covershrinkagescale_p = as.integer((trueY<=Ushrinkagescale_p)&&(trueY>=Lshrinkagescale_p))
 

# Add coverage result to matrix for macro result:
#===============================================
#===============================================
macroResults[mRun,] = matrix(c(Lbasic,Ubasic,Coverbasic,Lpercentile,Upercentile,Coverpercentile, 
                               Lnormal,Unormal,Covernormal,
                               chat,
                               Lshrinkage,Ushrinkage,Covershrinkage,Lshrinkage_p,Ushrinkage_p,Covershrinkage_p,
                               Lshrinkagescale,Ushrinkagescale,Covershrinkagescale,
                               Lshrinkagescale_p,Ushrinkagescale_p,Covershrinkagescale_p),1,22)

}
# End macroloop:
#==============
#==============
#==============
macroResults

toc()
# this write command for parallel processors
#write.table(format(macroResults, digits=4),paste0("./JacksonD4runs/",as.character(NAI),"/JacksonD4out.txt"),
#            append=TRUE,eol="\n",
#            quote=FALSE,
#            row.names=FALSE,col.names=FALSE)

write.table(format(macroResults, digits=4),paste0("./JacksonD4out.txt"),
             append=TRUE,eol="\n",
             quote=FALSE,
             row.names=FALSE,col.names=FALSE)




