This folder contains the script files used to generate the experiment data included in the paper. The same files are duplicated in 'src' folder with the rest of the source codes for convenience. More details can be found in 'src.' Brief descriptions follow:
* RunMM1kExpr.R: runs the M/M/1/k expeirments in Section 6.1 of the paper.
* RunMG1kExpr.R: runs the M/G/1/k expeirments in Section 6.2 of the paper.
* DirectIUJackson_10ROAR.R & DirectIUJackson_40ROAR.R: run the Jackson network experiments in Section 6.3 of the paper; these codes generate the direct bootstrap-based confidence intervals in comparison. The difference between the two scripts is the size of the Jackson network tested.
* MMIUJackson_10ROAR.R & MMIUJackson_40ROAR.R: Same as above, but applies the metamodeling technique in "A Bayesian framework for quantifying uncertainty in stochastic simulation (2014) Operations Research 62(6):1439–1452" by Xie W, Nelson BL, Barton RR to construct the confdience intervals. 
