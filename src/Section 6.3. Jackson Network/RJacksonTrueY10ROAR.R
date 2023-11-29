# RJacksonTrueY4.R

# Purpose: Identify true mean waiting time (time in system less service time) for
# Metamodel-based input uncertainty quantification
#
# ROAR Collab Version for R main program =============
library(compiler)
library(resample)
library(pracma)
# Set to LEcuyer RNG
RNGkind(kind = "L'Ecuyer-CMRG")
# set different seed for each processor node
AI <- Sys.getenv("SLURM_ARRAY_TASK_ID")
#AI = 99 #TEMPORARY UNTIL RUN ON ROAR
NAI <- as.numeric(AI)
set.seed(NAI)


#=====================================================

library(reticulate) # This is the library you need to run python on R
#library(tictoc)     # Used in small test runs to see where large execution time

# Key R functions used: sample; predict; solve; repmat (not R -- custom version)
# Also rexp, rmultinom 

source("setPijk2.R")     # set values of network routing probabilities - setPijk2 is UNBALANCED
# setPijk is BALANCED
source("inputData.R")    # routine to generate 'original' sample of interarrival 
#  and routing data
source("runSim2.R")       # Jackson-like simulation routine
#  replaces GG1ksim used by Wei Xie
source("getY.R")         # Jackson-like simulation for Y0, Y00 and direct bootstrap evals

py_run_file("VBASim.py")
py_run_file("JN2c.py")    # load the Jackson Network simulator (python)
#  Do this only once within an experiment to avoid 
#  resetting the random number seeds


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
numArv = 10              # number of arrival nodes (assume exponential arrivals for now)
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
# True average wait determined by long simulation using bigSimReps
bigSimReps = 40000         # replications to establish trueY - 40000 per core
warmup = 100              # first customer waiting times deleted
runlength = 10           # number of customers recorded 
# Establish trueY with high replication simulation
# Get Y00 mean wait from original sample using large R0 replications and Y0 for CI:
#==================================================================
#tic()
#print("Y00")
Y00results = getY(bigSimReps,warmup, runlength, numArv, mArvTime, PijkHat) 
#toc()
# output of getY: mean, variance(over reps),#reps
trueY = Y00results[1]
varYhat = Y00results[2]

results = matrix(0,1,2)
results[1,] = matrix(c(trueY,varYhat),1,2)

write.table(format(results, digits=4),paste0("./JacksonD10runs/",as.character(NAI),"/JacksonD10out.txt"),
            append=TRUE,eol="\n",
            quote=FALSE,
            row.names=FALSE,col.names=FALSE)

