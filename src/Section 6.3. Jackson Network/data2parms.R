# data2parms.R
# Purpose: take array of bootstrap resampled data, convert to arrival means 
#  and 'coded' (via probMethod) routing probabilities - rows correspond to resamples.
# Date: 7/15/19            Author: Russell Barton

# Inputs: 
#   numArv - number of arrival nodes
#   PijkDim -  vector of dimensions of: (# layers - 1), (numnodes layer 1,)    
              # (numnodes layer 2), where PijkRoute = routing probability from 
              #  jth node in layer i to kth node in layer i+1
#   PijkRoute - true routing probabilities, use in constructing reduced parm vector
#   numVecs - number of parameter vectors to generate (e.g. B0 or other)
#   samples - the bootstrapped (BSsamples or newSamples)data
#   probMethod - how to transform probabilities: 
#                'allBut1','pseudoP','Cartesian' 

# Locals:
#   probVec - vector of routing probabilities for a single node
#   probNotZero - elements allowed to be nonzero in that vector
#   probLength - adds lengths of probability vectors for each node for total
#   Parms - numVecs x totLength array of parameters for use in metamodeling 

# Outputs:
#   BSparms[1] = number of parameters total in the new parameter vector
#   BSparms[2] = numVecs x BSparms[1] array



data2parms = function(numArv,PijkDim,PijkRoute,probMethod,numVecs,samples)
{


# No transformation currently used for interarrival time means
  meanArv = samples[[1]]

# Transform routing probabilities to a probability vector
  PijkBoot = samples[[2]]

  for (n in 1:numVecs){
    Pijk = PijkBoot[n,,,]
    probLength = 0   # total length of probability parameter vector
    for (i in 1:2) {
      for (j in 1:PijkDim[i+1]) {
          # get probability vector for node j at level i
          probVec = Pijk[i,j,]
          probNotZero = rep(0, length(probVec))
          probNotZero[PijkRoute[i,j,] > 0] = 1
          # convert to parameters based on probMethod
          if (probMethod == "allBut1") {
             probParm = aB1(probVec, probNotZero)
          } else if (probMethod == "pseudoP") {
             probParm = pP(probVec, probNotZero)
          # otherwise use "Cartesian" transformation
          } else {
             probParm = B2C(probVec, probNotZero)
          }
          # keep track of total length of probability parameter vector
          probLength = probLength + length(probParm)

          # for each i,j concatenate its routing parameters to the overall vector
          if (i!=1 | j!=1) TprobParm = c(TprobParm, probParm) else TprobParm = probParm
      }
    }
# probLength needed to define Parms, so define Parms on first pass (n=1)
  if (n == 1) Parms = matrix(0,numVecs, (numArv+probLength));
  Parms[n,] = c(meanArv[n,],TprobParm)
  }


# Assemble output
  totLength = numArv + probLength
  BSparms = list(totLength,Parms)
  return(BSparms);
}