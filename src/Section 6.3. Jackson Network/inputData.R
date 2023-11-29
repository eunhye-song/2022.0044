# inputData.R
# Purpose: given the distribution parameters, generates interarrival times
#           and routing counts
# Date : 7/6/2019            Author : Russell Barton
#        5/16/2023                    updated routing sample size to 10x

# Inputs:
#   numSamples - the sample size - assumed 10x for routing probs
#   numArv - number of arrival nodes (and interarrival distributions)
#   rateArv - true exponential arrival rates for numArv nodes
#   numNodes - array of number of nodes at first, second, third layers 
#   PijkRoute - true routing probabilities
#   PijkDim - vector of dimensions of: # layers - 1, numnodes stage 1, numnodes stage 2
#   PijkHat - routing probability from jth node in layer i to kth node in layer i+1

# Locals:
#   probVec - routing probabilities for a single node
# Outputs:
#   structure of "Idata" output
#     list(arrivals[numArv,numSamples],PijkHat[i,j,k])
#   arrivals - numArv x numSamples interarrival times
#   PijkHat - 2xnumArv x numArv empirical route probabilities
#            NOTE: route counts not needed to bootstrap resample
#                  only need empirical route probabilities.
inputData = function(numSamples,numArv,rateArv,PijkDim,PijkRoute)

{

# set sample size for routing probabilities, assumed 5x for routing probs
# typical setting numSamples = 100, numSamples2 = 1000
  numSamples2 = 10*numSamples

# set structure for interarrival times - all values initially set to zero
arrivals = matrix(0,numArv,numSamples)

# set structure for RouteCount (probabilities will be replaced by counts)
RouteCount = PijkRoute 

# generate interarrivals based on exponential distribution with rate parmArv[i,1]
for (i in 1:numArv) {
    arrivals[i,] = rexp(n = numSamples, rate = rateArv[i])
    }

# generate route counts based on multinomial - all routes out of j
for (i in 1:2) {
    for (j in 1:PijkDim[i+1]) {
        # get probability vector for node j at level i
        probVec = PijkRoute[i,j,]
        # generate multinomial data and fill RouteCount[i,j,]
        RouteCount[i,j, ] = rmultinom(1,numSamples2,probVec)
    }
}

# Convert routing counts into empirical probabilities
PijkHat = RouteCount/numSamples2

X = list(arrivals,PijkHat)

return(X)
}

