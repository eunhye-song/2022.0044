# hyperE2.R
# Purpose: find the hyper-ellipsoid design space to ellipCov*100 
# percent of bootstrap samples of arrival and service times.
# Code similar to that of Wei Xie, but names changed by R. Barton
#  to match Jackson-like network setting, and added probabilities conversion.


# Inputs: 
#   ellipCov - covering percentage
#   probMethod - method for handling routing probability parameters
               # ('allBut1', 'pseudoP', 'Cartesian')
#   numSamples - assume number of samples for each input distribution the same
#   numArv - number of arrival distributions (use only one parameter for exponential)
#   PijkDim -  vector of dimensions of: (# layers - 1), (numnodes layer 1,)    
              # (numnodes layer 2), where PijkRoute = routing probability from 
              #  jth node in layer i to kth node in layer i+1
#   PijkRoute - true routing probabilities, use in constructing reduced parm vector
#   Idata - real world data
#   structure of "Idata"
#        list(arrivals[numSamples,numArv],PijkHat[i,j,k])

# Outputs: a list to define ellipsoid E
#   cE - center
#   covP - covariance matrix (may be modeled as diagonal fr large sample size, equal probs)
#   aE - radius

hyperE2 = function(ellipCov,probMethod,numSamples,numArv,PijkDim,PijkRoute,Idata)
{
  B0 = 721;         # B0 is the number of bootstrap samples initially generated, 
                    #  to fit initial ellipsoid 
  C = 707;          # C is the minimum number of points in a new sample that must fall 
                    #  in the ellpsoid. If less, then new points added to fit, ellipsoid 
                    #  recalculated, and a new bootstrap sample generated and tested. 


# generate arrivals and route probabilities
  BSsamples = BSgenerator2(numSamples,numArv,PijkDim,B0,Idata); 
              # structure of BSsamples = list(B0 x numArv means, B0 Pijk sets)
# convert arrivals and probabilities to parameter vector for MM
# BSparms structure: BSparms[1] = number of parameters, BSparms[2] = array of B0 x nparms
  BSparms = data2parms(numArv,PijkDim,PijkRoute,probMethod,B0,BSsamples);
  nparms = BSparms[[1]]
  parmArray = BSparms[[2]]

  sign = 1;  
  while(sign==1)    # generate more bootstrap points to do hypothesis test
  {
#    covP = diag(apply(parmArray,2,var));            # covariance structure modeled as diagonal
    covP = cov(parmArray)                           # general covariance structure
    cent = array(colMeans(parmArray),c(1,nparms));  # center (1xnparms)
    a = radEv2(parmArray,covP,cent,ellipCov);       # radius (1x1)
    
    # generate new bootstrap samples for test:
    newSamples = BSgenerator2(numSamples,numArv,PijkDim,B0,Idata); # (B0xnparms)
    BSparms = data2parms(numArv,PijkDim,PijkRoute,probMethod,B0,newSamples);
              # structure of BSparms = number of parameters total in the new 
              #  parameter vector, numVecs x BSparms[1] array
    newParmArray = BSparms[[2]]    
    # check number of samples in E:
    dist = newParmArray-repmat(cent,B0,1);
    temp = sum(diag(dist%*%solve(covP)%*%t(dist))<=a^2);
    if(temp>C) {                              # stop
      sign = 0;                               # accept this ellipsoid as design space
    }else {                                   # accumulate bootstrap samples to build E
      parmArray = rbind(parmArray,newParmArray);
    }
  }
  
  result = list(nparms,cE=cent,covE=covP,aE=a);  
  return(result);
}
