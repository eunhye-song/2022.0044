# getSimpVerts.R

# Purpose: using simpDim as an argument, returns a simpDim x simpDim+1
#  matrix with each column a vertex of a regular simpDim dimensional simplex
#  in the positive orthant with 'lowermost' point at zero
# From https://www.scilab.org/sites/default/files/neldermead.pdf

getSimpVerts = function(simpDim) {
  simpDim1 = simpDim + 1
  p = (1/(simpDim*sqrt(2)))*(simpDim-1+sqrt(simpDim1))
  q = (1/(simpDim*sqrt(2)))*(sqrt(simpDim1)-1)
  X = matrix(0, simpDim1, simpDim)

  for (j in 1:simpDim) {
    X[1,j] = 0.0
  }
  for (i in 2:simpDim1) {
    for (j in 1:simpDim) {
        if (j == i-1) {
           X[i,j] = X[1,j] + p
        }
        else {
           X[i,j] = X[1,j] + q
        }
    }
  }
return(X)
}
    
