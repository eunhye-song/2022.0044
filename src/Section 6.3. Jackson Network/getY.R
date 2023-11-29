# getY.R

# getY(simReps, warmup, runlength, numArv, mArvTime, PijkHat)

# Purpose: compute output for Y0, Y00 and direct bootstrap Ys. 
#          Report variance in case t-based CIs later and for SSW reconstruction.

# Inputs:
#   numArv = m - number of customers simulated (excluding warmup customers, whose data is discarded)
#   mArvTime - m-vector of mean arrival times 
#   PijkHat - sample routing probabilities in 'ijk' format,  convert to Pij form,
#             = routing prob from node i to node j,i = 1,2m+3, j = 1,2m+3,
#             nonzero values for only some of i = 1, ..., 2m; j = m+1, ...2m+3

# Locals:
#   m = numArv - number of arrival points == stations in the first and second stages
#   warmup - number of completed initial customers whose data to discard      
#   runlength - total number of completed customers you will observe within the replication,
#   simReps - number of replications

# Output:
#   results - vector of wait time mean, variance, and numReps

getY = function(simReps, warmup, runlength, numArv, mArvTime, PijkHat){

# Define results matrix:
#======================
results = rep(0,3)

# Initialize local simulation parameters:
#=======================================
# warmup delete to remove some bias     
# runlength length of runs
# simReps number of replications
m = numArv
PijkDim = c(2,numArv,numArv)     # vector of dimensions of: (# layers - 1), (numnodes layer 1,)    
                                 #  (numnodes layer 2), where PijkHat = sample routing 
                                 #  probability from jth node in layer i to kth node in layer i+1


# initializing the run settings - DO NOT MODIFY
py$init_run_params(as.integer(m), as.integer(warmup), as.integer(runlength), as.integer(simReps)) 

# arate and Pij just defined here - values changed for every DOE point (DOEnum)
arate = rep(0,m) # a vector of arrival rates to the stations in the first stage
Pij = matrix(0,2*m+3,2*m+3) # Routing probability; (2m+3) by (2m+3) matrix

# srate and Qcapa fixed across all DOE points
srate = rep(1.0,2*m+2) # a vector of service rates at all 2m + 2 stations (not needed at exit node 2m+3)
Qcapa = rep(10,2*m+2) # queue capacity for all 2m + 2 stations

# set arrival_rates
arate = 1.0/mArvTime

# Convert PijkHat values to Pij values
    for (i in 1:2) {
      if(i == 1) kMax = PijkDim[i+1] else kMax = 3  
      for (j in 1:PijkDim[i+1]) {
        for (k in 1:kMax) {
          iprime = j + (i-1)*m
          jprime = k + i*m
          Pij[iprime,jprime] = PijkRoute[i,j,k]
        }
      }
    }
# Add exit probabilities for nodes 2m+1, 2m+2
    Pij[2*m+1,2*m+3] = 1.0
    Pij[2*m+2,2*m+3] = 1.0


  # initializing the input parameters - DO NOT MODIFY
  py$init_input_params(arate, srate, as.integer(Qcapa), Pij)

  # Run simulation replications
  py$run_replications()

  # retrieve output, which is 
  #   the average cumulative waiting time per customer per replication
  output = py$WaitTime
  results = c(mean(output),var(output),simReps)
  results

# Return Results:
#===============

return(results)
}