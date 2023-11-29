# C2B.R

# Purpose: take a probability vector's representation in Cartesian coordinates 
#  on an n-1-dimensional simplex and generate n Barycentric coordinates relative
#  to a reference simplex (created by getSimpVerts)
#  used by mapDOEpts.R -- NOTE: DOE CPvec WILL be changed if has negative probs

# Date: 10/8/19            Author: Russell Barton

# Inputs:
#   CPvec - probability vector represented in Cartesian coordinates

# Locals:
# Cdim - length of CPvec = length(probVec) - 1
# SimpVert - matrix of simplex vertices in reduced space (Cdim+1) x Cdim
# CPvecMat - matrix form of probability vector needed by cart2bary

# Called Functions:
#   getSimpVerts(Cdim) - generates Cartesian coordinates of regular 
#     (Cdim+1) - dimension simplex centered at 0
#   cart2bary - from 'geometry' v0.4.1

# Locals:
#   augPvec   # CPvec converted to the nonzero probability elements

# Outputs:
#   probVec - full dimensional probability vector probVec = Pijk[i,j,]
#   NewCPvec - reduced dimension to replace reducedPvec in case prob <0

C2B = function(CPvec,probNotZero) 
{
  Cdim = length(CPvec)
  simpVert = getSimpVerts(Cdim)
  CPvecMat = t(as.matrix(CPvec))         # convert vector to matrix for cart2bary
  PvecMat = cart2bary(simpVert,CPvecMat)
  augPvec = as.vector(PvecMat)              # convert cart2bary matrix output to vector
  
  # make sure no negative probabilities from ellipsoid DOE parameter
     augPvec[augPvec < 0] = min(.01,augPvec[augPvec > 0])
  # make sure augPvec sums to one
  augPvec = augPvec/sum(augPvec)
  # create updated CPvec (based on augPvec in case some prob < 0) to pass back
  NewCPvec = B2C(augPvec,probNotZero)

  # populate Pijk (i.e., probVec) with augPvec values
  probVec = probNotZero
  probVec[probNotZero==1] = augPvec

  pList = list(probVec,NewCPvec)
  return(pList);

}
    
