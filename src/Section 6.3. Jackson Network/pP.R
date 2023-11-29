# pP.R

# Purpose: take 'generic' vector of routing probabilities and 
#   produce values from elements not forced to zero - but now not forced to sum to 1
# Date: 7/11/19           Author: Russell Barton

# Inputs:
#   probvec - input probability vector, but not all elements are 'active'
#   probNotZero - elements allowed to be nonzero in that vector

# Locals:
#   reducedPvec - holds probVec elements not forced to zero
#   result - holds probVec for return to data2parms function
#   mult - random multiplier to move probability vector off hyperplane (with same p value)

# Outputs:
#   probability parameter vector suitable for metamodeling (full-dimensional variation)

pP = function(probVec,probNotZero)
{
  reducedPvec = probVec[probNotZero == 1]  # select elements allowed to be nonzero
  mult = runif(1,.5,1.5)
  result = mult*reducedPvec  # multiply all elements by 'mult' ~ U(.5, 1.5)

  return(result);
}
