# parms2probs2.R
# Purpose: take array of bootstrap resampled and 'coded' (via probMethod) 
#  routing probabilities - rows corrspond to resamples.
#  NOTE: this version modified over parms2probs to give two indices, not three
#  in order to match "JN.py" structure
# Note: C2B, IaB1, IpP modified to pass back corrections to MMDOE probabilities when negative
# Date: 10/8/19            Author: Russell Barton

# Inputs: 
#   nDOEpts - number of Pijk sets - number of DOE runs
#   PijkDim -  vector of dimensions of: (# layers - 1), (numnodes layer 1,)    
              # (numnodes layer 2), where PijkRoute = routing probability from 
              #  jth node in layer i to kth node in layer i+1
#   PijkRoute - used to determine nonzero probabilities to expand the probability
#                 parameter vector back to full Pijk dimension
#   probMethod - how to transform probabilities: 
#                'allBut1','pseudoP','Cartesian' 
#   probParms - DOE xformed probability parameter matrix: nDOEpts x total # xformed parameters

# Locals:
#   probVec - vector of routing probabilities for a single node in Pijk form
#   probNotZero - elements allowed to be nonzero in that vector
#   probVecLength - length of probability vector for node i,j
#   probParm - metamodel's probability parameter vector for node i,j in nth run
#   pStart - where probParm starts in the overall metamodel probability vector
#   m = assume equal layer depth = PijkDim[2] = PijkDim[3] 

# Outputs:
#   Pnij - array of Pij values (i and j unique node numbers) for nth simulation run in DOE
#          structure of i and j: first phase nodes 1:m, second phase m+1:2m, 
#          third phase 2m+1:2m+2, exit 2m+3
#   Mapping from {i,j,k} to {i',j'}: i' = j + (i-1)*m; j' = k + i*m



parms2probs2 = function(nDOEpts,PijkDim,PijkRoute,probMethod,probParms)
{

  Pnijk = array(0,c(nDOEpts,PijkDim)) # set Pnijk entries to zero 

# compute layer depth
  m = PijkDim[2]

# set sizes for Pnijk and Pnij arrays
  Pnijk = array(0,c(nDOEpts,PijkDim)) # set Pnijk entries to zero
  Pnij = array(0,c(nDOEpts,(2*m)+3,(2*m)+3)) # set Pnij entries to zero 

  for (n in 1:nDOEpts){
    pStart = 1   # start index of probability parameter vector for Pij set,
    for (i in 1:2) {
      for (j in 1:PijkDim[i+1]) {
          # get probability vector for node j at level i
          probVec = PijkRoute[i,j,]          # temporary value for ProbVec - to be set later
          probNotZero = rep(0, length(probVec))
          probNotZero[PijkRoute[i,j,] > 0] = 1
          probVecLength = sum(probNotZero)   # total numer of nonzero prob entries for node i,j

          # convert to parameters based on probMethod
          if (probMethod == "allBut1") {
             probParm = probParms[n,pStart:(pStart+probVecLength-2)]
             temp = IaB1(probParm, probNotZero)
             probVec = temp[[1]]
             probParm = temp[[2]]
             # update ProbParms in case IaB1 did probability fix
             probParms[n,pStart:(pStart+probVecLength-2)]= probParm
             pStart = pStart+probVecLength-1 # update start for next Pij set in parm vector
             
          } else if (probMethod == "pseudoP") {
             probParm = probParms[n,pStart:(pStart+probVecLength-1)]
             temp = IpP(probParm, probNotZero)
             probVec = temp[[1]]
             probParm = temp[[2]]
             # update ProbParms in case pseudoP did probability fix
             probParms[n,pStart:(pStart+probVecLength-1)]= probParm
             pStart = pStart+probVecLength   # update start for next Pij set in parm vector

          # otherwise use "Cartesian" transformation
          } else {
             probParm = probParms[n,pStart:(pStart+probVecLength-2)]
             temp = C2B(probParm, probNotZero)
             probVec = temp[[1]]
             probParm = temp[[2]]
             # update ProbParms in case C2B did probability fix
             probParms[n,pStart:(pStart+probVecLength-2)]= probParm
             pStart = pStart+probVecLength-1 # update start for next Pij set in parm vector
          }
          Pnijk[n,i,j,] = probVec
      }
    }
  }
# Assemble output: put Pnijk in Pnij form

  
  for (n in 1:nDOEpts){
    for (i in 1:2) {

      # Careful! nonzero Pijk ONLY for k <= 3 when i = 2 - do not set higher than 3
      if(i == 1) kMax = PijkDim[i+1] else kMax = 3  

      for (j in 1:PijkDim[i+1]) {
        for (k in 1:kMax) {
          iprime = j + (i-1)*m
          jprime = k + i*m
          Pnij[n,iprime,jprime] = Pnijk[n,i,j,k]
        }
      }
    }
# Add exit probabilities for nodes 2m+1, 2m+2
    Pnij[n,2*m+1,2*m+3] = 1.0
    Pnij[n,2*m+2,2*m+3] = 1.0
  }


# Output is full probability vector 
  pList = list(Pnij,probParms)
  return(pList);
}