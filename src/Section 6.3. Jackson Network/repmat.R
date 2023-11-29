# equivalent of repmat function in matlab. 
# Purpose: repeat matrix X to generate a large matrix.

repmat = function(X,m,n){
  mx = dim(X)[1];
  nx = dim(X)[2];
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=TRUE);
}