MM1k.LindleySteadyState.Parametric <- function(lambda, mu, k){
  # simulates a single observation of the steady-state time in system for M/M/1/k queue
  # lambda: interarrival rate
  # mu: service rate
  # k: system capacity of the M/M/1/k queue
  
  # the initial number in system is sampled from the steady-state NIS distribution for M/M/1/k
  # setting up M/M/1/k steady-state NIS distribution
  rho = lambda/mu
  ssprob = 0:k
  ssprob = (1-rho)*rho^ssprob/(1-rho^(k+1))
  ssprob_rescaled = ssprob[1:(k+1)]/sum(ssprob[1:(k+1)]) 
  ssprob_rescaled = cumsum(ssprob_rescaled) 
  
  # we sample the initial NIS at the time of the first arrival (conditional on NIS>=1)
  initial_prob <- runif(1)
  initial_NIS <- min(which(ssprob_rescaled>=initial_prob))-1 # number in system before the first arrival    
  
  # sample all service times needed for this replication from the service time distribution
  Sr = rexp(1+initial_NIS,mu) 
  
  W = matrix(0,1+initial_NIS,1) # waiting time
  TIS = matrix(0,1+initial_NIS,1) # time in system
  D = matrix(0,1+initial_NIS,1) # departure time
  
  # variables to store time stamp
  clock = 0
  prevA = 0
  nextA = 0
  Adiff = 0
  
  if(initial_NIS>0)
  {
    D[1] <- Sr[1]
    if(initial_NIS == k) { 
      repeat{
        nextA = rexp(1,lambda)
        if(prevA + nextA > Sr[1])
        {
          break;
        }else{
          prevA = prevA + nextA
        }
      }
      Adiff = prevA+nextA 
    }
    if(initial_NIS>1)
    {
      for(x in 2:initial_NIS)
      {
        D[x] = D[x-1] + Sr[x]
      }  
    }
    W[initial_NIS+1] = max(D[initial_NIS] - Adiff,0)
    TIS[initial_NIS+1] = W[initial_NIS+1] + Sr[initial_NIS+1]
    D[initial_NIS+1] = TIS[initial_NIS+1]
    NIS = initial_NIS + 1
  }  else {
    # initially the system is empty
    W[1] = 0
    TIS[1] = Sr[1]            
    D[1] = clock + Sr[1] 
    NIS = 1
  } 
  return(TIS[initial_NIS+1])
}
