# BSgenerator2.R
# Purpose: generate B bootstrap resample distribution parameter vectors


# Inputs: 
#   B - number of bootstrap samples
#   Idata: list(arrivals([numSamples,numArv]),PijkHat[i,j,k])

# Intermediate:
#   Asamples: mean of bootstrapped interarrival times for each bootstrap resample

# Output:  'samples' - list with bootstrapped interarrival time means (Asamples) 
#          and bootstrapped routing probabilities (PijkBoot)
BSgenerator2 = function(numSamples,numArv,PijkDim,B,Idata)
{
  Asamples = matrix(0,B,numArv);  # B bootstrap resample moments (B x numArv)
  
  # bootstrap interarrival time data:
    # extract interarrival times
  arrivals = Idata[[1]];
  for (i in 1:numArv) {
    arrivalVec = arrivals[i,]
    temp1 = matrix(sample(arrivalVec,size=B*numSamples,replace=TRUE),B,numSamples)
    Asamples[,i] = rowMeans(temp1)
  }

  # bootstrap routing data:
    # extract routing probabilities
  PijkHat = Idata[[2]];

  # set dimension of PijkBoot, RouteCount
  PijkBoot = array(0, dim = c(B,PijkDim));
  RouteCount = PijkHat;

  # generate route counts based on multinomial - all routes out of j
  for (b in 1:B){
    for (i in 1:2) {
      for (j in 1:PijkDim[i+1]) {
        # get probability vector for node j at level i
        probVec = PijkHat[i,j,]
        # generate multinomial data and fill RouteCount[i,j,]
        RouteCount[i,j, ] = rmultinom(1,numSamples,probVec)
      }
    }
    # Convert routing counts into empirical probabilities
      PijkBoot[b,,,] = RouteCount/numSamples
  }

  # define samples
  samples = list(Asamples,PijkBoot)
  
  return(samples);
}
