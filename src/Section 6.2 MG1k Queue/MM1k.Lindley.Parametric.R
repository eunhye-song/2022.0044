MM1k.Lindley.Parametric <- function(lambda, mu, k, r){
    # simulates the time in system starting empty and idle for the M/M/1/k queue
    # lambda: interarrival rate
    # mu: service rate
    # k: system capacity of the M/M/1/k queue
    # r: number of observed time in system for this replication

    # sample all service times needed for this replication from the service time distribution
    Sr = rexp(r,mu)
    
    bigr = 2*k + 2*r # extra arrivals to allow for blocking
    Avec = rexp(bigr,lambda) # pre-generated interarrival times
    
    W = matrix(0,r,1) # waiting time
    TIS = matrix(0,r,1) # time in system
    D = matrix(0,r,1) # departure time
    Arr = matrix(0,r,1) # arrival time
    
    prevAi = 0
    departed = 0;
    
    clock = Avec[1] 
    Arr[1,] = clock
    W[1,] = 0;
    TIS[1,] = Sr[1];
    D[1,] = clock + Sr[1];
    num = 1
    ii = 1
    for(i in 2:r)
    {
        repeat
        {
            ii = ii+1
            Ai = prevAi + Avec[ii] 
            clockPlusAi = clock + Ai
            j = sum(D[(i-num):(i-1),] > clockPlusAi)

            if(j < k)
            {
                break;   
            }
            else
            {
                prevAi = Ai;
            }
        }
        
        num = j + 1 
        clock = clockPlusAi
        Wi = W[(i-1),]+Sr[(i-1)]-Ai
        Wi = (Wi+abs(Wi))/2 
        W[i,] = Wi
        TIS[i,] = Wi + Sr[i]
        D[i,] = clock + TIS[i,]
        departed = departed + (j-1)
        prevAi = 0
        
        Arr[i,] = clock
    }
    
    list(data = TIS, Ybar = mean(TIS), S2 = var(TIS))
}
