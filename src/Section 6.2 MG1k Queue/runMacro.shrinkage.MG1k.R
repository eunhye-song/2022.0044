runMacro.shrinkage.MG1k <- function(lambda, betaMix_p, Balpha1a, Balpha2a, Balpha1b, Balpha2b,
                                         k,alpha,ETIS,m,R0,R,macrorep,n,B)
{
  # inputs
  # lambda: arrival rate
  # betaMix_p, Balpha1a, Balpha2a, Balpha1b, Balpha2b: parameters to set up the Beta mixture service time distribution
  # k: system capacity
  # alpha: confidence interval nominal coverage rate is 1-alpha
  # R0: number of simulation replications made at the empirical cdf
  # R: number of simulation replications made at each bootstrap sample
  # macrorep: number of macro runs
  # n: number of real-world observations for each input model
  # B: number of bootstrap samples
  # ETIS: the expected time in system
  
  lambdahat = matrix(0,macrorep,1)
  muhat = matrix(0,macrorep,1)
  m0 = 0 # only used for m_used calculation
  
  #basic bootstrap resampling
  sim_BB_LB = matrix(0,macrorep,1)
  sim_BB_UB = matrix(0,macrorep,1)
  prob_correct_BB = rep(0,macrorep)
  
  #percentile bootstrap resampling
  sim_PB_LB = matrix(0,macrorep,1)
  sim_PB_UB = matrix(0,macrorep,1)
  prob_correct_PB = rep(0,macrorep)
  
  #percentile parametric exponential-exponential bootstrap resampling
  sim_pPBe_LB = matrix(0,macrorep,1)
  sim_pPBe_UB = matrix(0,macrorep,1)
  prob_correct_pPBe = rep(0,macrorep)
  
  # z without R0 error term
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
  
  clo = matrix(0,macrorep,1)

  for(i in 1:macrorep)
  {
      arrival = rexp(n,lambda); # real-world data for arrival process
      
      mixIndex = as.numeric(runif(n,0,1) > betaMix_p) + 1
      bimodalBoth = matrix(c(rbeta(n,Balpha1a,Balpha2a),rbeta(n,Balpha1b,Balpha2b)),n,2)
      service = bimodalBoth[cbind(seq_len(nrow(bimodalBoth)), mixIndex)] # real-world data for service process
      
      Yoriginal = matrix(0,R0,1)
      for(j in 1:R0)
      {
          Yoriginal[j] = GG1k.Lindley(arrival,service,k,m)$Ybar;
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
              Ybr[b,r] =  GG1k.Lindley(boot_arrival[,b],boot_service[,b],k,m)$Ybar;
          }
          Yhat[b] = mean(Ybr[b,])
          S2[b] = var(Ybr[b,])
          ebr[b,] = Ybr[b,] - Yhat[b]
          
      }
      
      Yhatb = sort(Yhat) - Ybarbar 
      Yhat = sort(Yhat) 
      
      # Basic bootstrap method
      sim_BB_LB[i,] = Ybarbar - Yhatb[ceiling((1-alpha/2)*B)]
      sim_BB_UB[i,] = Ybarbar - Yhatb[ceiling((alpha/2)*B)]
      prob_correct_BB[i] = (sim_BB_LB[i]<ETIS)&&(ETIS<sim_BB_UB[i])
      
      # Normal CI
      SSb = var(Yhat)*(B-1)
      SSw = sum(S2*(R-1))
      V2 = max((SSb/(B-1)) - (SSw/(B*R*(R-1))),0)
      ################ to avoid length-0 CI for the normal case
      if(V2==0) {V2=Voriginal/R0}
      ################
      sim_N_LB[i,] = Ybarbar-qnorm(1-alpha/2)*sqrt(V2);
      sim_N_UB[i,] = Ybarbar+qnorm(1-alpha/2)*sqrt(V2);
      prob_correct_N[i] = (sim_N_LB[i]<ETIS)&&(ETIS<sim_N_UB[i])
      
      # Percentile bootstrap method
      sim_PB_LB[i,] = Yhat[ceiling((alpha/2)*B)]
      sim_PB_UB[i,] = Yhat[ceiling((1-alpha/2)*B)]
      prob_correct_PB[i] = (sim_PB_LB[i]<ETIS)&&(ETIS<sim_PB_UB[i])
      
      # With shrinkage correction, SSB and SSP
      Yvals = cbind(Yhat, S2)
      nBCI <<- B
      Reps <<- R
      c_val = function(Yvals,inds){
          SSB = (nBCI-1)*var(Yvals[inds,1])
          chat = 1-sqrt(max(0,1-sum(Yvals[inds,2])/Reps/SSB)) 
      }
      boot.out = boot(Yvals,c_val,1000)
      c = boot.ci(boot.out,conf=.95,type="perc")$percent[4]
      clo[i,] = c
      
      # basic
      Yhathat = c*mean(Yhat) + (1-c)*Yhat
      Yhathatb = Yhathat - Ybarbar 
      sim_SSB_LB[i,] = Ybarbar - Yhathatb[ceiling((1-alpha/2)*B)]
      sim_SSB_UB[i,] = Ybarbar - Yhathatb[ceiling((alpha/2)*B)]
      prob_correct_SSB[i] = (sim_SSB_LB[i]<ETIS)&&(ETIS<sim_SSB_UB[i])

      # percentile
      sim_SSP_LB[i,] = Yhathat[ceiling((alpha/2)*B)]
      sim_SSP_UB[i,] = Yhathat[ceiling((1-alpha/2)*B)]
      prob_correct_SSP[i] = (sim_SSP_LB[i]<ETIS)&&(ETIS<sim_SSP_UB[i])
      
      # Quantile bootstrap shrinkage correction QSB, QSP
      sim_QSB_LB[i,] = Ybarbar - (1-c)*Yhatb[ceiling((1-alpha/2)*B)]
      sim_QSB_UB[i,] = Ybarbar - (1-c)*Yhatb[ceiling((alpha/2)*B)]
      prob_correct_QSB[i] = (sim_QSB_LB[i]<ETIS)&&(ETIS<sim_QSB_UB[i])
      sim_QSP_LB[i,] = Ybarbar + (1-c)*Yhatb[ceiling((alpha/2)*B)]
      sim_QSP_UB[i,] = Ybarbar + (1-c)*Yhatb[ceiling((1-alpha/2)*B)]
      prob_correct_QSP[i] = (sim_QSP_LB[i]<ETIS)&&(ETIS<sim_QSP_UB[i])
      
      
      # Parametric percentile bootstrap - exponential service model
      Ybr_p = matrix(0,B,(2*R))
      Yp = matrix(0,B,1)

      # Set up the parametric bootstrap
      lambdaHat = 1/mean(arrival)
      muHat = 1/mean(service)
      boot_arrival = matrix(rexp((n*B),lambdaHat),n,B);
      boot_service = matrix(rexp((n*B),muHat),n,B);
      boot_lambdap = apply(boot_arrival,2,mean);
      boot_lambdap = 1/boot_lambdap; #convert mean interarrival time to rate
      boot_mup = apply(boot_service,2,mean);
      boot_mup = 1/boot_mup #convert mean service time to rate
      
      for (b in 1:B)
      {
          for (r in 1:(R))
          {
              Ybr_p[b,r] =  MM1k.Lindley.Parametric(boot_lambdap[b],boot_mup[b],k,m)$Ybar;
          }
          Yp[b] = mean(Ybr_p[b,])
      }
      Yp = sort(Yp)
      # Percentile interval calculation
      sim_pPBe_LB[i,] = Yp[ceiling((alpha/2)*B)]
      sim_pPBe_UB[i,] = Yp[ceiling((1-(alpha/2))*B)]
      prob_correct_pPBe[i] = (sim_pPBe_LB[i]<ETIS)&&(ETIS<sim_pPBe_UB[i])
      
  }

  sim_BB_width = colMeans(sim_BB_UB-sim_BB_LB)
  sim_BB_sd = apply(sim_BB_UB-sim_BB_LB,2,sd)
  coverage_BB = mean(prob_correct_BB)
  sim_PB_width = colMeans(sim_PB_UB-sim_PB_LB)
  sim_PB_sd = apply(sim_PB_UB-sim_PB_LB,2,sd)
  coverage_PB = mean(prob_correct_PB)
  sim_pPBe_width = colMeans(sim_pPBe_UB-sim_pPBe_LB)
  sim_pPBe_sd = apply(sim_pPBe_UB-sim_pPBe_LB,2,sd)
  coverage_pPBe = mean(prob_correct_pPBe)
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
  m_used = m-m0
  clo_mean = mean(clo)
  clo_var = var(clo)


  list(
       lambda = lambda, alpha = alpha, k = k, m_used = m_used, R0 = R0, R = R, n = n,
       pPBe_w = sim_pPBe_width, pPBe_sd = sim_pPBe_sd, pPBe_cov = coverage_pPBe,
       BB_w = sim_BB_width, BB_sd = sim_BB_sd, BB_cov = coverage_BB,
       PB_w = sim_PB_width, PB_sd = sim_PB_sd, PB_cov = coverage_PB,
       N_w = sim_N_width, N_sd = sim_N_sd, N_cov = coverage_N,
       SSB_w = sim_SSB_width, SSB_sd = sim_SSB_sd, SSB_cov = coverage_SSB,
       SSP_w = sim_SSP_width, SSP_sd = sim_SSP_sd, SSP_cov = coverage_SSP,
       QSB_w = sim_QSB_width, QSB_sd = sim_QSB_sd, QSB_cov = coverage_QSB,
       QSP_w = sim_QSP_width, QSP_sd = sim_QSP_sd, QSP_cov = coverage_QSP,
       clo_mean=clo_mean, clo_var = clo_var)
}
