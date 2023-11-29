There are four main R scripts to create the coverage results:

* DirectIUJackson_4ROAR.R:
  This code simulates capacitated Jackson network with 4 input nodes and computes confidence interval coverage for all methods except the metamodel-based input uncertainty intervals.

* DirectIUJackson_10ROAR.R:
  This code simulates capacitated Jackson network with 10 input nodes and computes confidence interval coverage for all methods except the metamodel-based input uncertainty intervals.

* MMIUJackson_4ROAR.R:
This code simulates capacitated Jackson network with 4 input nodes and computes confidence interval coverage for metamodel-based input uncertainty intervals.

* MMIUJackson_10ROAR.R: 
This code simulates capacitated Jackson network with 10 input nodes and computes confidence interval coverage for metamodel-based input uncertainty intervals.

Key parameters to set in these four routines:

* numSamples: sample size for arrival rate estimation (10x this for routing probability estimation).

* numArv: number of arrival nodes (m), set to 4 or 10 in the program files above, but can set to other values.

* nBootCI: number of bootstrap resamples used in computing a confidence interval, set too 220 for m = 4, 1000, for m = 10.

* mRuns: for computing coverage, multiple confidence intervals must be constructed. To use 1000, set mRuns = 1000 (very long run) or mRuns = 10 for each of 100 processors running the code in parallel. The latter for m = 10 metamodeling assessment took more than 43 hours (times 100 processors) on the Penn State ROAR computig cluster.

* runLength: number of samples from simulation run to include in the average cumulative waiting time.

* trueY: this is used to compute coverage and must be changed if you change runLength since after a 100-customer warmup the system is not in steady state.


All other routines in this folder are necessary to support the computations. Each has a brief header in the code describing ther function. Two Python files are required to do the network simulation:

* JN2c.py: This simulates the capacitated Jackson network.

* VBASim.py, BasicClassesMultiDraft.py, Basic_Classes.py: These contain a set of discrete-event simulation functions called by JN2c.py; This code is implemented by Barry L. Nelson. The details of the software can be found in his book, Foundations and Methods of Stochastic Simulation: A First Course, Springer (2013), and website: http://users.iems.northwestern.edu/~nelsonb/IEMS435/index1e.html
