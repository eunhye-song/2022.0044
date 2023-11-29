# radEv2.R
# Purpose: calculate radius of hyper-ellipse E which covers
#          ellipCov*100 percent of bootstrap points
# Code is identical to that of Wei Xie, but variable names changed by R. Barton
#  to match Jackson-like network setting.

# Inputs:
#   parmArray - matrix of resampled distribution parameter vectors
#       b indexes rows -- bootstrap resample index
#   covP - covariance matrix of parmArray
#   cent - mean (center) of parmArray
#   ellipCov - covering probability of ellipsoid

# Output: "a" - radius of covering ellipsoid (it is conservative to speed up running.)

radEv2 = function(parmArray,covP,cent,ellipCov)
{
  B = nrow(parmArray);                     # number of samples
  temp = sort(sqrt(diag((parmArray-repmat(cent,B,1))%*%solve(covP)%*%
    t(parmArray-repmat(cent,B,1)))));      # standardized distances to center
  a = temp[ceiling(B*ellipCov)];
  return(a);
}
