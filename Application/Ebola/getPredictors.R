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
library("OutbreakTools")

# clear workspace
rm(list = ls())



# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# get the names of all files with 1 gene
log <- list.files(path="./dta", pattern="*.GLM.log", full.names = TRUE)

for (i in 1:length(log)){
  t <- read.table(log[i], header=TRUE, sep="\t")
  new.vals = data.frame(matrix(ncol=0,nrow=length(t$posterior)))
  for (j in seq(1, 46)){
    h_name <- paste("location.coefficientsTimesIndicators", j, sep="")
    new.vals[, h_name] <- t[,h_name]
  }
  if (i==1){
    indicators = new.vals
  }else{
    indicators = rbind(indicators, new.vals)
  }
}

predictor.mapping <- read.table("predictor_description.tsv", header=F, sep="\t")
  

ind.names = names(indicators)
for (i in seq(1, length(ind.names))){
  non_zero = indicators[,ind.names[i]]!=0
  new.vals = data.frame(name=predictor.mapping[i,2], description=predictor.mapping[i,1],
                        indicator=sum(non_zero)/length(indicators[,1]), coefficient=mean(indicators[non_zero, ind.names[i]]))
  if (i==1){
    predictors = new.vals
  }else{
    predictors = rbind(predictors, new.vals)
  }
  
}


# plot the indicators
indicator.val <- c()
for (i in seq(1, length(predictors$name))){
  indicator.val[i] = predictors[i,"indicator"]
}






sorted_indices = sort(indicator.val, decreasing = TRUE, index.return=TRUE)
bar.data = data.frame()
for (j in seq(length(sorted_indices[[2]]), 1)){
  i = sorted_indices[[2]][j]
  new.bar.data = data.frame(name = predictors[i,"name"], ind =  predictors[i,"indicator"], coeff=predictors[i,"coefficient"])
  if (sorted_indices[[1]][j] > 0.3){
    bar.data = rbind(bar.data, new.bar.data)
  }
}

plot.Sym.rate <- ggplot(data=bar.data, aes(x=name,y=ind)) +
  geom_bar(position="dodge",stat="identity") +
  coord_flip() +
  theme_bw() + theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     legend.position="none") + 
  scale_fill_manual(values = c("posterior" = "black", "prior" = "grey"))
plot(plot.Sym.rate)
ggsave(plot=plot.Sym.rate,"",width=2, height=2)




