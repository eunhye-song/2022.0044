# IaB1.R

# Purpose: takes parameterized 'all-but-one' vector of probabilities, converts to 
#          Pijk form
#  used by mapDOEpts.R -- NOTE: DOE reducedPvec WILL be changed if has negative probs

# Date: 7/15/19            Author: Russell Barton

# Inputs:
#   reducedPvec - input probability parameters
#   probNotZero - elements allowed to be nonzero in the Pijk vector

# Locals:
#   augPvec   # reducedPvec with remaining probability added as last element

# Outputs:
#   probVec - full dimensional probability vector probVec = Pijk[i,j,]
#   NewreducedPvec - reduced dimension to replace reducedPvec in case prob <0

IaB1 = function(reducedPvec,probNotZero)

{
  # augment reducedPvec
  augPvec = c(reducedPvec, (1-sum(reducedPvec)))

  # make sure no negative probabilities from ellipsoid DOE parameter
  augPvec[augPvec < 0] = min(.01,augPvec[augPvec > 0])
  # make sure augPvec sums to one
  augPvec = augPvec/sum(augPvec)
  # UPDATE reducedPvec (based on augPvec in case some prob < 0) to pass back
  NewreducedPvec = augPvec[1:length(reducedPvec)]

  # populate Pijk (i.e., probVec) with augPvec values
  probVec = probNotZero
  probVec[probNotZero==1]=augPvec

  pList = list(probVec,NewreducedPvec)
  return(pList)
}

