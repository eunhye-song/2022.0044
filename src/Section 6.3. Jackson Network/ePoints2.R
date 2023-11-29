# ePoints2.R
# Purpose: generate points evenly distributed in unit hyper-sphere
# Date : 11/10/2012              'ePoints' Author : Wei Xie
# Date: 7/16/2019   completely rewritten to allow for high-dimension generation
#                                'ePoints2' Author: Russell Barton

# Input: 
#   npoints - number of DOE points
#   ncoords - dimension of metamodel parameter space

# Locals:
# HypCube - an npoints latin hypercube design on [0,1]^ncoords
# nval(coord) - a normal random variate (quantile) based on a latin hypercube probability
# radius - a random radius chosen between zero and 1 with heavier
#          emphasis on values near 1 through(1/ncoords) power of U(0,1)
# nvalLength - length of the nval() vector for a particular point

# Called Functions:
#   qnorm(quantile, mean(default = 0), StdDev(default=1))
#   runif(nvals,lower,upper)

# Output:
#   X: points generated in unit hyper-sphere (kxd)

ePoints2 = function(npoints,ncoords)
{
  # the method employed by Wei Xie was presented in Sun and Farooq (2002). That method 
  # required integration of sin values raised to powers from 2 or 3 to on the order of 100
  # with impractical roundoff errors 

  # the method here uses the samples in the latin hypercube to drive univariate normal (0,1) generation,
  # then applies the method by Martin Roberts (Method #20) at
  # http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/ 

  nval = rep(0,ncoords)
  X = matrix(0,npoints,ncoords)
  HypCube = matrix(0,npoints,ncoords)
 
  HypCube = maximinLHS(npoints,ncoords,dup=5);  # (npointsxnccords) design matrix in [0,1]^d
                                                # dup is factor that determines the number of 
                                                # candidate points used in the search. A multiple 
                                                # of the number of remaining points than can be added.
  for(point in 1:npoints){
     for(coord in 1:ncoords){
        nval[coord] = qnorm(HypCube[point,coord]) # qnorm without extra args generates quantile for N(0,1)
     }
     radius = runif(1,0,1)^(1./ncoords)
     nvalLength = sqrt(sum(nval*nval))
     X[point,] = radius*nval/nvalLength
  }

  return(X);
}