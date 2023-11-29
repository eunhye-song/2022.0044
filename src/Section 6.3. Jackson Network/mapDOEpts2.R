# mapDOEpts2.R
# Purpose: take array of parameter vectors used by metamodel, convert to 
#  set of arrival rates and raw routing probabilities.
# Date: 10/8/19            Author: Russell Barton

# Inputs: 
#   numArv - number of arrival nodes
#   PijkDim -  vector of dimensions of: (# layers - 1), (numnodes layer 1,)    
              # (numnodes layer 2), where PijkRoute = routing probability from 
              #  jth node in layer i to kth node in layer i+1
#   PijkRoute - used to determine nonzero probabilities to expand the probability
#                 parameter vector back to full Pijk dimension
#   PijkHat - sample estimates for Pijk - used to fix any infeasible probability vectors
#   probMethod - how to back-transform to probabilities: 'allBut1','pseudoP','Cartesian') 
#   nDOEpts - number of parameter vectors to be converted for sim runs
#   DOEpts - the parameter vectors for each run in the DOE

# Locals:
#   arrival_means - nDOEpts x numArv array of mean interarrival times

# Outputs:
#   SimDOE[1] = arrival_rates array, nDOEpts x numArv
#   SimDOE[2] = nDOEpts x Pij (full set of routing_probs indexed just by node id) array

mapDOEpts2 = function(numArv,PijkDim,PijkRoute,probMethod,nDOEpts,MMDOEpts)
{

# Transform interarrival means to arrival rates for each of n=nDOEpts simulation runs
  arrival_means = MMDOEpts[,1:numArv]
# check for nonpositive interarrival means in DOE, adjust to smallest positive in MMDOE,
  if(min(arrival_means) <= 0) {
     warning("at least one negative mean interarrival time from DOE - fixed")
     arrival_means[arrival_means <= 0] = min(arrival_means[arrival_means > 0])
     # update MMDOEpts with corrected ones
     MMDOEpts[,1:numArv] = arrival_means
  }
  arrival_rates = 1/arrival_means

# Transform probability parameter vector to nDOEpts sets of routing probabilities Pijk 
  probParms = MMDOEpts[,(numArv+1):ncol(MMDOEpts)]
  temp = parms2probs2(nDOEpts,PijkDim,PijkRoute,probMethod,probParms)
  Pnij = temp[[1]]
# update probability parameters in MMDOEpts with corrections to negative values
  probParms = temp[[2]]
  MMDOEpts[,(numArv+1):ncol(MMDOEpts)] = probParms

# Assemble output
SimDOE = list(arrival_rates, Pnij)

  return(SimDOE);
}