# 2014 Ebola epidemic in Sierra Leone
Folder contains the materials to recreate the Applications part of the manuscript

Runs were done in 3 replicates and then combined after a burnin of 5% (I.e. the _out/Ebola_glm_mascot.log_ needs no burnin). The output is in the _out_ folder, including each of the replicates. If tracer cannot read them, running them through logCombiner fixes that.

To run the xmls, both the CoupledMCMC and the Mascot package have to be installed

`getConvergence.R` analyses the convergence of the 3 different runs using ESS values and scale reduction factors. The ESS values estimates by coda in this script can differ from the ESS values in tracer.

`getPredictors.R` plots the predictors.

`plotPrior.R` plots the prior and posterios on the sum of active predictors.
