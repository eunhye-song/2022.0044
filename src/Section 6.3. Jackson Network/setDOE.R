# setDOE.R

# Purpose:
# Find nDOEpts space-filling on the ellipsoid capturing most bootstrap parameter values.
#  originally in main code of Wei Xie, pulled out for main code clarity.
# Date: 7/15/19            Author: Russell Barton

# Inputs:
#   nDOEpts - number of points in the design
#   ellipDOE - the ellipsoid determined to hold most bootstrap parameter vectors
#   ellipDOE = list(cE=cent,covE=covP,aE=a) the center, covariance, and size

# Outputs:
#   MMDOEpts - an nDOEpts x length(cE) set of parameter vector settings for
#                making simulation runs to fit the metamodel

setDOE = function(nDOEpts,ellipDOE)
{
  nparms = ellipDOE[[1]]
  ctr = ellipDOE[[2]]
  covP = ellipDOE[[3]] # NOTE: covP is diagonal for network model with exponential arr and serv
  a = ellipDOE[[4]]

  # generate k points based on Latin hypercube in unit hyper-sphere centered at 0:
  dp_4D = ePoints2(nDOEpts,nparms); # result is nDOEpts x nparms matrix
  C0_4D = chol(covP);
  DOE = t(repmat(t(ctr),1,nDOEpts)+a*C0_4D%*%t(dp_4D)); # Linear rescale to ellipsoid

  # remove points with non-positive means or invalid probs is done in mapDOEpts2

  return(DOE);
}
