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


post <- read.table("./out/Ebola_glm_mascot.log", header=TRUE, sep="\t")
prior <- read.table("./prior/Ebola_prior.log", header=TRUE, sep="\t")

Ne = data.frame(val=prior$sum.NeIndicator., method="prior")
Ne = rbind(Ne, data.frame(val=post$sum.migrationIndicator., method="posterior"))

p_Ne <- ggplot() +
  geom_histogram(data=prior,aes(x=sum.NeIndicator.,y=..density..,fill="Prior"), binwidth=1, colour="black",alpha=0.8) +
  geom_histogram(data=post,aes(x=sum.NeIndicator.,y=..density..,fill="Posterior"), binwidth=1, colour="black",alpha=0.8)+
  scale_fill_manual(name="",values=c("red", "blue")) + xlab("") + xlim(c(0,20)) +
  ggtitle("active Ne predictors")
  

plot(p_Ne)

p_mig <- ggplot() +
  geom_histogram(data=prior,aes(x=sum.migrationIndicator.,y=..density..,fill="Prior"), binwidth=1, colour="black",alpha=0.8) +
  geom_histogram(data=post,aes(x=sum.migrationIndicator.,y=..density..,fill="Posterior"), binwidth=1, colour="black",alpha=0.8) +
  scale_fill_manual(name="",values=c("red", "blue")) + xlab("")+ xlim(c(0,20)) +
  ggtitle("active migration predictors")
  
plot(p_mig)

ggsave(plot=p_Ne,"../../Figures/Ebola/Ne_prior.pdf",width=3.5, height=3.5)
ggsave(plot=p_mig,"../../Figures/Ebola/Migration_prior.pdf",width=3.5, height=3.5)


# p_all <- grid.arrange(p_ne, p_mig, layout_matrix=hlay)
# ggsave(plot= p_all ,paste("Ebola_predictors.pdf", sep=""),width=8, height=8)
