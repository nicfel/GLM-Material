######################################################
######################################################
# Analyses inferred topology and node heights of
# species trees
######################################################
######################################################
library(ggplot2)
# needed to get the node heights
library(phytools)
# needed to read the trees
library(ape)
library(coda)
library(grid)
library(gridExtra)

# clear workspace
rm(list = ls())


col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)
col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)


# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)



# get the names of all SISCO first (of three) run log files
log <- list.files(path="./out", pattern="*.log", full.names = TRUE)

calc_stats = c("posterior", "prior", "treeLikelihood.CP1", "treeLikelihood.CP2", "treeLikelihood.CP3", "treeLikelihood.ig", "sum.migrationIndicator.", "sum.NeIndicator.", "migrationClock", "NeClock")

t = list()
for (i in seq(1,3)){
# for (i in seq(1,2)){
    
  print(i)
  # Make the filenames for all the three runs
  filename1 <- paste(log[i], sep="")
  
  # Read in the log files of the runs with fixed trees
  t_tmp = read.table(filename1, header=TRUE, sep="\t")
  print(length(t_tmp$Sample))
  # the 889 is required since the gelmam methods needs all chains to be equally long
  t[[i]] <- as.mcmc(t_tmp[seq(1,889),], start = 1, end = numeric(0), thin = 1)
  
  if (i==1){
    t.tot = t_tmp[-seq(1,ceiling(length(t_tmp$Sample)/20)), ]
  }else{
    t.tot = rbind(t.tot, t_tmp[-seq(1,ceiling(length(t_tmp$Sample)/20)), ])
  }

}

gelstats = gelman.diag(t, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=F)

gelstatsvals=gelstats$psrf
essvals = effectiveSize(t.tot)

print("minimal ess value for which the gelman stats is above 1.01")
print(min(essvals[which(gelstatsvals[,1]>1.05)]))
print("lowest ess value of any parameter")
print(min(essvals[seq(2,length(essvals))]))

plot(t.tot$posterior, type="line")


