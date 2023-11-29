runMacro.shrinkage.MM1k <- function(lambda,mu,k,alpha,R0,R,macrorep,n,B)
{
  # inputs
  # lambda: arrival rate
  # mu: service rate
  # k: system capacity
  # alpha: confidence interval nominal coverage rate is 1-alpha
  # R0: number of simulation replications made at the empirical cdf
  # R: number of simulation replications made at each bootstrap sample
  # macrorep: number of macro runs
  # n: number of real-world observations for each input model
  # B: number of bootstrap samples
  
  rho = lambda/mu;

  # True time-in-system; estimated from 10^8 macro runs of MM1k.LindleySteadyState.Parametric
  TIS = (rho==0.7)*3.090888 + (rho==0.9)*4.863077
  
  # estimates from the "real-world" sample in each macro replication
  lambdahat = matrix(0,macrorep,1)
  muhat = matrix(0,macrorep,1)

  #basic bootstrap resampling
  sim_BB_LB = matrix(0,macrorep,1)
  sim_BB_UB = matrix(0,macrorep,1)
  prob_correct_BB = rep(0,macrorep)
  
  #percentile bootstrap resampling
  sim_PB_LB = matrix(0,macrorep,1)
  sim_PB_UB = matrix(0,macrorep,1)
  prob_correct_PB = rep(0,macrorep)
  
  # z without R0 error term, near zero when R0 large
  sim_N_LB = matrix(0,macrorep,1)
  sim_N_UB = matrix(0,macrorep,1)
  prob_correct_N = rep(0,macrorep)    
  
  # sample shrinkage correction basic SSB
  sim_SSB_LB = matrix(0,macrorep,1)
  sim_SSB_UB = matrix(0,macrorep,1)
  prob_correct_SSB = rep(0,macrorep)
  
  # sample shrinkage correction percentile SSP
  sim_SSP_LB = matrix(0,macrorep,1)
  sim_SSP_UB = matrix(0,macrorep,1)
  prob_correct_SSP = rep(0,macrorep)
  
  #quantile basic bootstrap (shrinkage)
  sim_QSB_LB = matrix(0,macrorep,1)
  sim_QSB_UB = matrix(0,macrorep,1)
  prob_correct_QSB = rep(0,macrorep)
  
  #quantile percentile bootstrap (shrinkage)
  sim_QSP_LB = matrix(0,macrorep,1)
  sim_QSP_UB = matrix(0,macrorep,1)
  prob_correct_QSP = rep(0,macrorep)

  # c_lo estimated for shrinkage methods
  clo = matrix(0,macrorep,1)
  
  for(i in 1:macrorep)
  {
    arrival = rexp(n, lambda); # real-world data for arrival process
    service = rexp(n, mu); # real-world data for service process

    Yoriginal = matrix(0,R0,1)
    for(j in 1:R0)
    {
      Yoriginal[j] = GG1k.LindleySteadyState(arrival,service,k);
    }
    Yoriginal = sort(Yoriginal)
    Ybarbar = mean(Yoriginal)
    Soriginal = sd(Yoriginal)
    Voriginal = var(Yoriginal)
    
    a_bootstrap_indices = samp.bootstrap(n, B) # each column is a bootstrap sample
    s_bootstrap_indices = samp.bootstrap(n, B)
    
    boot_arrival = matrix(arrival[a_bootstrap_indices],n,B);
    boot_service = matrix(service[s_bootstrap_indices],n,B);
    
    Ybr = matrix(0,B,R)
    ebr = matrix(0,B,R)
    Yhat = matrix(0,B,1)
    S2 = matrix(0,B,1)
    
    for (b in 1:B)
    {
      # for bootstrap samples
      for (r in 1:R)
      {
        Ybr[b,r] =  GG1k.LindleySteadyState(boot_arrival[,b],boot_service[,b],k);
      }
      Yhat[b] = mean(Ybr[b,])
      S2[b] = var(Ybr[b,])
      ebr[b,] = Ybr[b,] - Yhat[b]
    }
    
    Yhatb = sort(Yhat) - Ybarbar # samples for basic bootstrap below
    Yhat = sort(Yhat) # Yhat assumed sorted below for quantile calculation via ceiling
    
    # Basic bootstrap method
    sim_BB_LB[i,] = Ybarbar - Yhatb[ceiling((1-alpha/2)*B)]
    sim_BB_UB[i,] = Ybarbar - Yhatb[ceiling((alpha/2)*B)]
    prob_correct_BB[i] = (sim_BB_LB[i]<TIS)&&(TIS<sim_BB_UB[i])
    
    # Normal CI method
    SSb = var(Yhat)*(B-1)
    SSw = sum(S2*(R-1))
    V2 = max((SSb/(B-1)) - (SSw/(B*R*(R-1))),0)
    ################ to avoid length-0 CI for the normal case
    if(V2==0) {V2=Voriginal/R0}
    ################
    sim_N_LB[i,] = Ybarbar-qnorm(1-alpha/2)*sqrt(V2); 
    sim_N_UB[i,] = Ybarbar+qnorm(1-alpha/2)*sqrt(V2);
    prob_correct_N[i] = (sim_N_LB[i]<TIS)&&(TIS<sim_N_UB[i])
    
    # Percentile bootstrap method
    sim_PB_LB[i,] = Yhat[ceiling((alpha/2)*B)]
    sim_PB_UB[i,] = Yhat[ceiling((1-alpha/2)*B)]
    prob_correct_PB[i] = (sim_PB_LB[i]<TIS)&&(TIS<sim_PB_UB[i])
    
    
    # SSB and SSP
    Yvals = cbind(Yhat, S2)
    nBCI <<- B
    Reps <<- R
    c_val = function(Yvals,inds){
      SSB = (nBCI-1)*var(Yvals[inds,1])
      chat = 1-sqrt(max(0,1-sum(Yvals[inds,2])/Reps/SSB)) # here, 1 is used instead of B/(B-1)
    }
    boot.out = boot(Yvals,c_val,1000)
    c = boot.ci(boot.out,conf=.95,type="perc")$percent[4]
    clo[i,] = c

    # basic
    Yhathat = c*mean(Yhat) + (1-c)*Yhat
    Yhathatb = Yhathat - Ybarbar # used for basic bootstrap
    sim_SSB_LB[i,] = Ybarbar - Yhathatb[ceiling((1-alpha/2)*B)]
    sim_SSB_UB[i,] = Ybarbar - Yhathatb[ceiling((alpha/2)*B)]
    prob_correct_SSB[i] = (sim_SSB_LB[i]<TIS)&&(TIS<sim_SSB_UB[i])
    
    # percentile
    sim_SSP_LB[i,] = Yhathat[ceiling((alpha/2)*B)]
    sim_SSP_UB[i,] = Yhathat[ceiling((1-alpha/2)*B)]
    prob_correct_SSP[i] = (sim_SSP_LB[i]<TIS)&&(TIS<sim_SSP_UB[i])
    
    # Quantile bootstrap shrinkage correction QSB, QSP
    sim_QSB_LB[i,] = Ybarbar - (1-c)*Yhatb[ceiling((1-alpha/2)*B)]
    sim_QSB_UB[i,] = Ybarbar - (1-c)*Yhatb[ceiling((alpha/2)*B)]
    prob_correct_QSB[i] = (sim_QSB_LB[i]<TIS)&&(TIS<sim_QSB_UB[i])
    sim_QSP_LB[i,] = Ybarbar + (1-c)*Yhatb[ceiling((alpha/2)*B)]
    sim_QSP_UB[i,] = Ybarbar + (1-c)*Yhatb[ceiling((1-alpha/2)*B)]
    prob_correct_QSP[i] = (sim_QSP_LB[i]<TIS)&&(TIS<sim_QSP_UB[i])
    
  }
  


  sim_BB_width = colMeans(sim_BB_UB-sim_BB_LB)
  sim_BB_sd = apply(sim_BB_UB-sim_BB_LB,2,sd)
  coverage_BB = mean(prob_correct_BB)
  sim_PB_width = colMeans(sim_PB_UB-sim_PB_LB)
  sim_PB_sd = apply(sim_PB_UB-sim_PB_LB,2,sd)
  coverage_PB = mean(prob_correct_PB)
  sim_N_width = colMeans(sim_N_UB-sim_N_LB)
  sim_N_sd = apply(sim_N_UB-sim_N_LB,2,sd)
  coverage_N = mean(prob_correct_N)
  sim_SSB_width = colMeans(sim_SSB_UB-sim_SSB_LB)
  sim_SSB_sd = apply(sim_SSB_UB-sim_SSB_LB,2,sd)
  coverage_SSB = mean(prob_correct_SSB)
  sim_SSP_width = colMeans(sim_SSP_UB-sim_SSP_LB)
  sim_SSP_sd = apply(sim_SSP_UB-sim_SSP_LB,2,sd)
  coverage_SSP = mean(prob_correct_SSP)
  sim_QSB_width = colMeans(sim_QSB_UB-sim_QSB_LB)
  sim_QSB_sd = apply(sim_QSB_UB-sim_QSB_LB,2,sd)
  coverage_QSB = mean(prob_correct_QSB)
  sim_QSP_width = colMeans(sim_QSP_UB-sim_QSP_LB)
  sim_QSP_sd = apply(sim_QSP_UB-sim_QSP_LB,2,sd)
  coverage_QSP = mean(prob_correct_QSP)
  clo_mean = mean(clo)
  clo_var = var(clo)
  

  list(lambda = lambda, alpha = alpha, k = k, R0 = R0, R = R, n = n,
       BB_cov = coverage_BB, sim_BB_width = sim_BB_width, sim_BB_sd = sim_BB_sd/sqrt(macrorep),
       PB_cov = coverage_PB, sim_PB_width = sim_PB_width, sim_PB_sd = sim_PB_sd/sqrt(macrorep),
       N_cov = coverage_N, sim_N_width = sim_N_width, sim_N_sd = sim_N_sd/sqrt(macrorep),
       SSB_cov = coverage_SSB, sim_SSB_width = sim_SSB_width, sim_SSB_sd = sim_SSB_sd/sqrt(macrorep),
       SSP_cov = coverage_SSP, sim_SSP_width = sim_SSP_width, sim_SSP_sd = sim_SSP_sd/sqrt(macrorep),
       QSB_cov = coverage_QSB, sim_QSB_width = sim_QSB_width, sim_QSB_sd = sim_QSB_sd/sqrt(macrorep),
       QSP_cov = coverage_QSP, sim_QSP_width = sim_QSP_width, sim_QSP_sd = sim_QSP_sd/sqrt(macrorep),
       clo_mean=clo_mean, clo_var = clo_var)
}
