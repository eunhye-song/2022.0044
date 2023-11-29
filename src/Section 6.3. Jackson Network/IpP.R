# IpP.R

# Purpose: take pseudoprobability vector and produce routing probabilities 
#   of Pijk form
# Date: 10/8/19           Author: Russell Barton

# Inputs:
#   pPvec - input pseudoprobability vector
#   probNotZero - elements allowed to be nonzero in the full vector

# Locals:
#   NewpPvec - reduced dimension to replace reducedpPvec in case prob <0
#   Pvec - elements now nonnegative and sum to 1


# Outputs:
#   probVec - full dimensional probability vector probVec = Pijk[i,j,]

IpP = function(pPvec,probNotZero)
{
  # make sure no negative probabilities from ellipsoid DOE parameter
     NewpPvec = pPvec
     NewpPvec[pPvec < 0] = min(.01,pPvec[pPvec > 0])

  # Pvec is corrected pPvec, and sum to one
     Pvec = NewpPvec/sum(NewpPvec)

  # populate Pijk (i.e., probVec) with Pvec values
  probVec = probNotZero
  probVec[probNotZero==1] = Pvec

  pList = list(probVec,NewpPvec)
  return(pList);
}
