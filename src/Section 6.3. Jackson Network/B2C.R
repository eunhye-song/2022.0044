# B2C.R

# Purpose: take 'generic' vector of routing probabilities and 
#   produce values from elements not forced to zero - as cartesian coordinates
#   representation on lower-dimension (by one) simplex, based on barycentric 
#   (p vector) representation
#   called by data2parms

# Inputs:
#   probVec - input probability vector, but not all elements are 'active'
#   probNotZero - elements allowed to be nonzero in that vector

# Locals:
#   reducedPvec - holds probVec elements not forced to zero
#   x - temporary holder of the Cartesian vector
#   Cdim - length of the Cartesian vector = length(Pvec) - 1
#   simpVert - matrix of simplex vertices in reduced space (Cdim+1) x Cdim

# Called Function:
#   getSimpVerts(Cdim) - generates Cartesian coordinates of regular 
#     (Cdim+1) - dimension simplex in positive orthant of the starting point 

# Outputs:
#   x - Cartesian coordinate representation of Pvec, initially = 0

B2C = function(probVec,probNotZero) 
{
  reducedPvec = probVec[probNotZero == 1]  # select elements allowed to be nonzero
  Cdim = length(reducedPvec)-1
  x = rep(0, Cdim)
  simpVert = getSimpVerts(Cdim)

for (i in 1:length(reducedPvec)) {  # <----- this line fixed 

    for (j in 1:Cdim) {
        x[j] = x[j] + reducedPvec[i]*simpVert[i,j]  
    }
}
return(x)
}
    
