# aB1.R

# Purpose: take 'generic' vector of routing probabilities and 
#   extract all-but-1 values from the elements not forced to be zero
# Date: 7/13/19           Author: Russell Barton

# Inputs:
#   probVec - input probability vector, but not all elements are 'active'
#   probNotZero - elements allowed to be nonzero in that vector

# Locals:
#    reducedPvec - holds probVec for return to data2parms function

# Outputs:
#   probability parameter vector suitable for metamodeling (full-dmensional variation)

aB1 = function(probVec,probNotZero)
{
  reducedPvec = probVec[probNotZero == 1]         # select elements allowed to be nonzero
  result = reducedPvec[1:(length(reducedPvec)-1)]      # remove last element (all but 1)

  return(result);
}