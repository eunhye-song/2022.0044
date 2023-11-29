# MMIUJackson_10ROAR.R
# Updated by RRB July 2023

# Purpose: Metamodel-based input uncertainty quantification
# Specialized for large-scale Jackson network-like (capacitated) example
# Uses same strategy and R routines as Wei Xie code for Xie, Nelson, Barton
# Routines general where possible, but treatment of service time
#  distributions and routing distributions treated separately, and specially

# Uses R packages for latin hypercube sampling for MM doe (lhs) and 
#  spatial correlation with nugget (mlegp)
# ROAR Collab Version for R main program =============
library(compiler)
library(resample)
library(pracma)

# timing via tictoc
library(tictoc)

library(lhs)        # for Latin hypercube start for metamodel DOE
library(mlegp)      # for Gaussian process metamodel
library(reticulate) # This is the library you need to run python on R
library(geometry)   # Uses R package 'geometry' for cart2bary Cartesian to Barycentric
                      # conversion for routing probability vectors
                      # (bary2cart not used - this was easy to code directly in B2C.R)
# Key R functions used: sample; predict; solve; repmat (not R -- custom version)
# Also rexp, rmultinom 

# Set to LEcuyer RNG
RNGkind(kind = "L'Ecuyer-CMRG")
# set different seed for each processor node
# AI <- Sys.getenv("SLURM_ARRAY_TASK_ID") # used for parallel runs on ROAR cluster
AI = "99" #TEMPORARY UNTIL RUN ON ROAR
NAI <- as.numeric(AI)
set.seed(NAI)


# Substitute Wei Xie version of repmat??
source("repmat.R")

# Other functions from Wei Xie (ALL modified) for DOE on ellipsoid
source("BSgenerator2.R") # generate parameter values from resampled data
source("hyperE2.R")      # ellipsoid (cE=center, covE=covmatrix, aE=radius) 
                         # here, covmatrix modeled as either as full or diagonal
                         # covering 'ellipCov'*100 percent of bootstrap resampled parameters
source("radEv2.R")       # calculate radius of ellipsoid covering 'ellipCov'*100 percent of data
source("ePoints2.R")     # calculates evenly spaced design on hypersphere for large dimension
                         # (used by hyperE2 -- both modified versions of Wei Xie routines)
# Other functions added by Russell Barton for code clarity or network structure
# setPijk2 to be updated
source("setPijk2.R")     # set values of network routing probabilities - setPijk2 is UNBALANCED
                         # setPijk is BALANCED
source("inputData.R")    # routine to generate 'original' sample of interarrival 
                         #  and routing data
source("setDOE.R")       # DOE construction code from Wei Xie in her main program - but w/out Cov
source("getSampRow.R")   # add sampled parameter estimates to DOE
source("data2parms.R")   # convert bootstrapped interarrival means and Pijks to a
                         #  set of parameter vectors for the metamodel
source("parms2probs2.R") # convert probability parameter vectors from DOE to Pij form for sim. input
                         # old parms2probs put in Pijk form rather than Pij form
source("aB1.R")          # convert bootstrapped Pijks to vectors via all-but-1
source("IaB1.R")         # convert all-but-one probability parameter vector to Pijk form
source("pP.R")           # convert bootstrapped Pijks to vectors via pseudo probs
source("IpP.R")          # convert pseudo-probability parameter vector to Pijk form
source("B2C.R")          # convert bootstrapped (barycentric) Pijks to vectors via Cartesian coords
source("C2B.R")          # convert Cartesian probability parameter vector to (barycentric) Pijk form
source("getSimpVerts.R") # gets Cartesian coodinates in n-space of an n+1 simplex for B2C and C2B
source("mapDOEpts2.R")   # back-convert DOE parameters to arrival rates and Pij
                         #  for input to simulation runs (mapDOEpts mapped to Pijk rather than Pij)
source("runSim2c.R")     # Jackson-like simulation routine
                         #  replaces GG1ksim used by Wei Xie
source("getY.R")         # Jackson-like simulation for R0 run estimate of mean wait 
                         #  with original sample

# overall timing
tic()

# Set key parameters:
#===================
#=================== 
#  input data characteristics (and compute from these either 
#  analytically or visa simulation the true expected value for coverage), 
#  ellipsoid %, sim run length, warmup and reps, CI significance level, 
#  number of bootstrap resamples, number of macro replications to check coverage

# Set input data characteristics (number of arrival/stage 2 nodes is numArv)
#===============================
numSamples = 100         # assume number of samples for each distribution the same 
#    use either 100 (with reps 5, 20, 50)
#    or 400 (with reps 20, 80, 200) (well behaved if large)
numArv = 10               # number of arrival nodes (assume exponential arrivals for now)
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


PijkDim = c(2,numArv,numArv)          # matrix of dimensions of: (# layers - 1), (numnodes layer 1,)    
#  (numnodes layer 2), where PijkRoute = routing probability from 
#  jth node in layer i to kth node in layer i+1
#  NOTE: when i = 2, k takes only three values, not numArv

# Set unbalanced Pijk parameters by calling 'setPijk2'
PijkRoute = array(0, dim = PijkDim)


# Set unbalanced routing to use balanced parameters, use setPijk instead of setPijk2
PijkRoute = setPijk2(PijkRoute,deltaP) 


# Set other key parameters (method, runlength and warmup set in runSim and getY)
#=========================
CImethod = "Metamodel"    # output file different for each: "Direct" or "Metamodel" 
nBootCI = 1000            # number of bootstrap runs to compute CI - should be 1000
mRuns = 10                # 1000 macro replications on 100 processors
                          # NOTE: 95% CI for nom. coverage of 90% at +/- 1% requires 3500 macroReps
                          # 95% CI for nom. coverage of 90% with 1000 reps is +/- 1.96%
                          # NOTE - these are conservative since less variance if nom. coverage near 95%
simReps = 10              # number of replications for sample mean, sd
                          # R0simReps is default number of initial replications R0
R0simReps = ceiling(numSamples^1.1)  
trueY = 0                 # must be determined by long simulation using bigSimReps
warmup = 100              # first customer waiting times deleted (set to 100 in real run)
runlength = 10            # number of customers per replication (including those deleted)

alp = 0.05                # setting for confidence interval percentage (95% if alp = 0.05)
ellipCov = 0.99           # metamodel DOE ellipsoid sized to cover 'ellipCov' of 
                          # bootstrapped parameter vectors

#probMethod = 'allBut1'   # probMethod is method for handling routing probability parameters
probMethod = 'Cartesian'  # ('allBut1', 'pseudoP', 'Cartesian')
#probMethod = 'pseudoP'

# First, establish trueY with high replication simulation
#trueY = 3.365 # SD .00272 CI "3.360","3.370" for runlength = 5
trueY = 3.406 # SD .00245 CI "3.402","3.411" for runlength = 10
#trueY = 3.682 # SD .00214 CI "3.678","3.686" for runlength = 50
#trueY = 4.525 # SD .00183 CI "4.522","4.529" for runlength = 400

#==================================================================


# Begin macro loop:
#=================
#=================
# Initialize some vectors and arrays
#=================
macroResults = matrix(0,mRuns,6) # matrix for holding macro results of CIs, their size and coverage
runOutput = NULL
BootYbar = NULL


#mRuns = 2 # use for testing before implementing macro loop
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

#print("Y0sims")
#tic()
# Make initial simulations with sampled data
# reps and runlength R0simReps = ceiling(numSamples^1.1) runlength = 100
Y0results = getY(R0simReps,warmup, runlength, numArv, mArvTime, PijkHat) 
# output of getY: mean, variance(over reps),#reps
Y0 = Y0results[1]
#toc()
# Construct MM DoE:
#================= 
#=================

# Construct ellipsoid based on bootstrapping original sample
#===========================================================
#tic()
#print("ellipDOE")
ellipDOE = hyperE2(ellipCov,probMethod,numSamples,numArv,PijkDim,PijkRoute,Idata)
#toc()
# structure of "ellipDOE"
# list(nparms,cE=cent,covE=covP,aE=a)
# modified from We Xie code to allow covariance to be diagonal
nparms = ellipDOE[[1]]

# Construct nDOEpts design on ellipDOE
#=====================================
#nDOEpts set to 220 for m = 4 (22 parms) to 1000 for m = 10 (112 params)
nDOEpts = 999 # set to one less than desired since observed parameter estimate included 

tic()
print("setDOE")
MMDOEpts = setDOE(nDOEpts,ellipDOE) # set design in MM parameter space
toc()
# structure of MMDOEpts is an nDOEpts x nparms matrix

# Add sample DOE parameter vector as last of MMDOEpts
#======================================================
sampParms = getSampRow(numArv,mArvTime,PijkDim,PijkHat,probMethod)
MMDOEpts = rbind(MMDOEpts, sampParms)
nDOEpts1 = nDOEpts+1                # design augmented by 1 to include simulation of sampled parameters

# Run simulation at DOE points:
#=============================
#=============================

# Map DOEpts into arrival rates and routing probabilities for simulation (simDOEpts)
#==================================================================================
#  including fix any design MMDOEpts with negative means or infeasible probability vectors

simDOEpts = mapDOEpts2(numArv,PijkDim,PijkRoute,probMethod,nDOEpts1,MMDOEpts)
# structure of simDOEpts is a list of two arrays
# simDOEpts[1] is an nDOEpts1 x numArv matrix of arrival rates
# simDOEpts[2] is routing probabilities using ordinal node indices
# it is a nDOEpts x totnodes x totnodes matrix, totnodes = 2xnumArv+3
arrival_rates = simDOEpts[[1]]
routing_probs = simDOEpts[[2]]

# Make simulation runs
#=====================
tic()
print("runSim")
simResults = runSim2c(simReps,warmup, runlength, nDOEpts1, numArv, arrival_rates, routing_probs)
toc()
# runlength 
# structure of "simResults"
# an nDOEpts x 2 matrix, each row is for a different experimental condition
# first column is mean of Y over simReps replications (for each condition)
# second column is the sample variance of Y over simReps replications 

# simResults[,1]

# Fit metamodel to DOE points and simulation runs:
#================================================
#================================================
# Note: use simResults + MMDOEpts to fit metamodel, NOT simResults + simDOEpts

# Fit Gaussian process metamodel
tic()
print("mlegp")
SK2 = mlegp(MMDOEpts,simResults[,1],constantMean=1,
            nugget=simResults[,2]/simReps,nugget.known=1,
            verbose=0); # no reporting of intermediate iterations
toc()

# Construct MM-based bootstrap values to use to construct CIs:
#============================================================
#============================================================

# Resample data and compute MM parameters
#========================================
#tic()
#print("bsGenerator")
newSamples = BSgenerator2(numSamples,numArv,PijkDim,nBootCI,Idata); # (nBootCIxnparms)
#toc()
BSparms = data2parms(numArv,PijkDim,PijkRoute,probMethod,nBootCI,newSamples);
# structure of BSparms = number of parameters total in the new 
#  parameter vector, numVecs x BSparms[1] array
BSParmArray = BSparms[[2]]

# Make predictions at bootstrap points
#=====================================
#tic()
#print("MMpreds")
#Ypred = predict(SK2,BSParmArray,se.fit=TRUE) # prediction means and std err
Ypred = predict(SK2,BSParmArray) 
#Y_hat = Ypred$fit                            # prediction means if 
#S_hat = Ypred$se.fit  # standard error of predictions (not currently used) 
#toc()
# adjust the following code for metamodel output (Y_hat = BootYbar=============
#BootYbar = Y_hat
BootYbar = Ypred

# Basic Bootstrap CI
Lbasic = Y0 - quantile(BootYbar - Y0,1-alp/2)
Ubasic = Y0 - quantile(BootYbar - Y0, alp/2)
Coverbasic = as.integer((trueY<=Ubasic)&&(trueY>=Lbasic))

# Percentile Bootstrap CI
Lpercentile = quantile(BootYbar,alp/2) 
Upercentile = quantile(BootYbar,1-alp/2)
Coverpercentile = as.integer((trueY<=Upercentile)&&(trueY>=Lpercentile))


# Add coverage result to matrix for macro result:
#===============================================
#===============================================
macroResults[mRun,] = matrix(c(Lbasic,Ubasic,Coverbasic,Lpercentile,Upercentile,Coverpercentile),1,6)
print(mRun)
}
# End macroloop:
#==============
#==============
#==============

# overall timing
toc()

# This write statement for parallel processors
#write.table(format(macroResults, digits=4),paste0("./JacksonMM10runs/",as.character(NAI),"/JacksonMM10out.txt"),
#            append=TRUE,eol="\n",
#            quote=FALSE,
#            row.names=FALSE,col.names=FALSE)

write.table(format(macroResults, digits=4),paste0("./JacksonMM10out.txt"),
            append=TRUE,eol="\n",
            quote=FALSE,
            row.names=FALSE,col.names=FALSE)





