# runSim2c.R

# runSim(simReps,nDOEpts, numArv, arrival_rates, routing_probs)
                         # structure of "simResults"
                         # an nDOEpts x 2 matrix, each row is for a different experimental condition
                         # first column is mean of Y over simReps replications (for each condition)
                         # second column is the sample variance of Y over simReps replications 

# Purpose: given the length of the two stages (m), arrival rates and routing probabilities
#   simulate 'simReps' replications of the Jackson-like network and return the mean and sample
#   variance for each experimental condition 1, 2, ..., nDOEpts.

# Inputs:
#   simReps - number of replications for each condition
#   nDOEpts - number of conditions (values of interarrival vector and Pijk matrix)
#   numArv - number of customers simulated (excluding warmup customers, whose data is discarded)
#   arrival_rates - m-vector of arrival rates - for each of m input nodes
#   routing_probs - routing probabilities in 'ij' format, Pij = routing prob from node i to node j,
#                   i = 1,2m+3, j = 1,2m+3, nonzero values for some of i = 1, ..., 2m; j = m+1, ...2m+3

# Locals:
#   DOEnum - index of the next design point to be run
#   m = numArv - number of arrival points == stations in the first and second stages
#   warmup - number of completed initial customers whose data to discard 
#             if 200, the initial 200 wait time observations are discarded
#   runlength - total number of completed customers you will observe within the replication,
#                 so 400 total, if first 200 discarded and runlength = 200

# Output:
#   results - tally of means and variances of time across replications for each DOE point

# library(reticulate) # This is the library you need to run python on R -- moved to MAIN
# py_run_file("JN2.py") # load the Jackson Network simulator; Do this only once within an experiment to avoid resetting the random number seeds -- moved to MAIN

runSim2c = function(simReps,warmup,runlength,nDOEpts, numArv, arrival_rates, routing_probs){

# Define results matrix:
#======================
results = matrix(0,nDOEpts,2)

# Initialize local simulation parameters:
#=======================================
m = numArv

# arate and Pij just defined here - values changed for every DOE point (DOEnum)
arate = rep(0,m) # a vector of arrival rates to the stations in the first stage
Pij = matrix(0,2*m+3,2*m+3) # Routing probability; (2m+3) by (2m+3) matrix
# srate and Qcapa fixed across all DOE points
srate = rep(1.0,2*m+2) # a vector of service rates at all 2m + 2 stations (not needed at exit node 2m+3)
Qcapa = rep(10,2*m+2) # queue capacity for all 2m + 2 stations

# loop over experimental conditions:
#==================================

for(DOEnum in 1:nDOEpts){
  py_run_file("VBASim.py")
  py_run_file("JN2c.py")  
  # initializing the run settings - DO NOT MODIFY
  py$init_run_params(as.integer(m), as.integer(warmup), as.integer(runlength), as.integer(simReps)) 
  
  # set arrival_rates and Pij for this DOE point
  arate = arrival_rates[DOEnum,]
  Pij = routing_probs[DOEnum,,]

  # initializing the input parameters - DO NOT MODIFY
  py$init_input_params(arate, srate, as.integer(Qcapa), Pij)

  # Run simulation replications
  py$run_replications()

  # retrieve output, which is 
  #   the average cumulative waiting time per customer per replication
  output = py$WaitTime
  # metamodel response is output2: mean and variance across the replications
  output2 = c(mean(output),var(output))
  results[DOEnum,1:2] = output2

  # To recover utilization for all 2m+2 stations
  # util = py$Utilization
  # util[[2]] # returns vector of utilizations for all 2m+2 stations
}
# end of loop over experimental conditions:
#=========================================

# Return Results:
#===============

return(results)
}