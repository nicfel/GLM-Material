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


t <- read.table("./out/Ebola_glm_mascot.log", header=TRUE, sep="\t")
t.labels = labels(t)
migration.predictors = data.frame(Sample=t$Sample)
ne.predictors = data.frame(Sample=t$Sample)

for (i in seq(1,length(t.labels[[2]]))){
  if (startsWith(t.labels[[2]][i], 'migrationGLMscaler')){
    name = strsplit(t.labels[[2]][i], split="\\.")[[1]][2]
    migration.predictors[,name] = t[, t.labels[[2]][i]]
  }
  if (startsWith(t.labels[[2]][i], 'NeGLMscaler')){
    name = strsplit(t.labels[[2]][i], split="\\.")[[1]][2]
    ne.predictors[,name] = t[, t.labels[[2]][i]]
  }
}


# make a dataframe with the support for incidence predictors
dat = data.frame(x=c("time to settlement with 100k", 
                     "gridded economic output", 
                     "population density", 
                     "population size", 
                     "mean temperature", 
                     "temperature seasonality", 
                     "mean precipitation", 
                     "precipitation seasonality", 
                     "cases happend 9 weeks earlier",
                     "cases happend 6 weeks earlier",
                     "cases happend 3 weeks earlier",
                     "cases happend 1 week earlier",
                     "unshifted case data",
                     "cases happend 1 week later",
                     "cases happend 3 weeks later",
                     "cases happend 6 weeks later",
                     "cases happend 9 weeks later"), y=c(

                                      length(which(ne.predictors[,"TT100k"]!=0))/length(ne.predictors$Sample),
                                      length(which(ne.predictors[,"Gecon"]!=0))/length(ne.predictors$Sample),
                                      length(which(ne.predictors[,"Pdens"]!=0))/length(ne.predictors$Sample),
                                      length(which(ne.predictors[,"Pop_Size"]!=0))/length(ne.predictors$Sample),
                                      length(which(ne.predictors[,"Temp"]!=0))/length(ne.predictors$Sample),
                                      length(which(ne.predictors[,"Tmpss"]!=0))/length(ne.predictors$Sample),
                                      length(which(ne.predictors[,"Prec"]!=0))/length(ne.predictors$Sample),
                                      length(which(ne.predictors[,"Precss"]!=0))/length(ne.predictors$Sample),                                                                                                                                                                        
                                      length(which(ne.predictors[,"incidence9weekearlier"]!=0))/length(ne.predictors$Sample),
                                      length(which(ne.predictors[,"incidence6weekearlier"]!=0))/length(ne.predictors$Sample),
                                      length(which(ne.predictors[,"incidence3weekearlier"]!=0))/length(ne.predictors$Sample),
                                      length(which(ne.predictors[,"incidence1weekearlier"]!=0))/length(ne.predictors$Sample),
                                      length(which(ne.predictors[,"incidence"]!=0))/length(ne.predictors$Sample),
                                      length(which(ne.predictors[,"incidence1weeklater"]!=0))/length(ne.predictors$Sample),
                                      length(which(ne.predictors[,"incidence3weeklater"]!=0))/length(ne.predictors$Sample),
                                      length(which(ne.predictors[,"incidence6weeklater"]!=0))/length(ne.predictors$Sample),
                                      length(which(ne.predictors[,"incidence9weeklater"]!=0))/length(ne.predictors$Sample)),
                 direction =c(
                   median(ne.predictors[which(ne.predictors[,"TT100k"]!=0),"TT100k"])>0,
                   median(ne.predictors[which(ne.predictors[,"Gecon"]!=0),"Gecon"])>0,
                   median(ne.predictors[which(ne.predictors[,"Pdens"]!=0),"Pdens"])>0,
                   median(ne.predictors[which(ne.predictors[,"Pop_Size"]!=0),"Pop_Size"])>0,
                   median(ne.predictors[which(ne.predictors[,"Temp"]!=0),"Temp"])>0,
                   median(ne.predictors[which(ne.predictors[,"Tmpss"]!=0),"Tmpss"])>0,
                   median(ne.predictors[which(ne.predictors[,"Prec"]!=0),"Prec"])>0,
                   median(ne.predictors[which(ne.predictors[,"Precss"]!=0),"Precss"])>0,                                                                                                                                                                        
                   median(ne.predictors[which(ne.predictors[,"incidence9weekearlier"]!=0),"incidence9weekearlier"])>0,
                   median(ne.predictors[which(ne.predictors[,"incidence6weekearlier"]!=0),"incidence6weekearlier"])>0,
                   median(ne.predictors[which(ne.predictors[,"incidence3weekearlier"]!=0),"incidence3weekearlier"])>0,
                   median(ne.predictors[which(ne.predictors[,"incidence1weekearlier"]!=0),"incidence1weekearlier"])>0,
                   median(ne.predictors[which(ne.predictors[,"incidence"]!=0),"incidence"])>0,
                   median(ne.predictors[which(ne.predictors[,"incidence1weeklater"]!=0),"incidence1weeklater"])>0,
                   median(ne.predictors[which(ne.predictors[,"incidence3weeklater"]!=0),"incidence3weeklater"])>0,
                   median(ne.predictors[which(ne.predictors[,"incidence6weeklater"]!=0),"incidence6weeklater"])>0,
                   median(ne.predictors[which(ne.predictors[,"incidence9weeklater"]!=0),"incidence9weeklater"])>0) )


p_ne <- ggplot(dat) + geom_bar(aes(x=x,y=y, fill=direction), stat="identity", colour="black") +
  xlab("") +
  geom_hline(yintercept = 0.7143, color="black", linetype="dashed") + 
  scale_fill_manual(values=c(col1,col3), guide=F)+
  ylab("inclusion probability") +
  ggtitle("Ne predictors") +
  scale_x_discrete(limits = c("time to settlement with 100k", 
                              "gridded economic output", 
                              "population density", 
                              "population size", 
                              "mean temperature", 
                              "temperature seasonality", 
                              "mean precipitation", 
                              "precipitation seasonality", 
                              "cases happend 9 weeks earlier",
                              "cases happend 6 weeks earlier",
                              "cases happend 3 weeks earlier",
                              "cases happend 1 week earlier",
                              "unshifted case data",
                              "cases happend 1 week later",
                              "cases happend 3 weeks later",
                              "cases happend 6 weeks later",
                              "cases happend 9 weeks later")) + coord_flip()

# plot(p_ne)
# ggsave(plot=p_ne,paste("Ebola_Ne_predictors.pdf", sep=""),width=4, height=2)




# make a dataframe with the support for incidence predictors
dat.mig = data.frame(x=c("great circle distance", 
                     "shared border", 
                     "time to settlement with 100k at origin", 
                     "time to settlement with 100k at destination", 
                     "gridded economic output at origin", 
                     "gridded economic output at destination", 
                     "population density at origin", 
                     "population density at destination", 
                     "population size at origin", 
                     "population size at destination", 
                     "mean temperature at origin", 
                     "mean temperature at destination", 
                     "temperature seasonality at origin", 
                     "temperature seasonality at destination", 
                     "mean precipitation at origin", 
                     "mean precipitation at destination", 
                     "precipitation seasonality at origin", 
                     "precipitation seasonality at destination", 
                     "Kailahun from",
                     "Kenema from",
                     "Kono from",
                     "Bombali from",
                     "Kambia from",
                     "Koinadugu from",
                     "PortLoko from",
                     "Tonkolili from",
                     "Bo from",
                     "Bonthe from",
                     "Moyamba from",
                     "Pujehun from",
                     "Western Rural from",
                     "Western Urban from",
                     "Kailahun to",
                     "Kenema to",
                     "Kono to",
                     "Bombali to",
                     "Kambia to",
                     "Koinadugu to",
                     "PortLoko to",
                     "Tonkolili to",
                     "Bo to",
                     "Bonthe to",
                     "Moyamba to",
                     "Pujehun to",
                     "Western Rural to",
                     "Western Urban to",
                     "samples at origin",
                     "samples at destination",
                     "ratio of weekly cases"), y=c(
                       
                       length(which(migration.predictors[,"greatCircleDistances"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"national_border_shared"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"originTT100k"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"destinationTT100k"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"originGecon"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"destinationGecon"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"originPdens"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"destinationPdens"]!=0))/length(migration.predictors$Sample),                                                                                                                                                                        
                       length(which(migration.predictors[,"originPop_Size"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"destinationPop_Size"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"originTemp"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"destinationTemp"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"originTmpss"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"destinationTmpss"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"originPrec"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"destinationPrec"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"originPrecss"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"destinationPrecss"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Kailahun_from"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Kenema_from"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Kono_from"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Bombali_from"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Kambia_from"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Koinadugu_from"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"PortLoko_from"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Tonkolili_from"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Bo_from"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Bonthe_from"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Moyamba_from"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Pujehun_from"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"WesternRural_from"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"WesternUrban_from"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Kailahun_to"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Kenema_to"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Kono_to"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Bombali_to"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Kambia_to"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Koinadugu_to"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"PortLoko_to"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Tonkolili_to"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Bo_to"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Bonthe_to"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Moyamba_to"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Pujehun_to"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"WesternRural_to"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"WesternUrban_to"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Samples_from"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"Samples_to"]!=0))/length(migration.predictors$Sample),
                       length(which(migration.predictors[,"incidenceRatio"]!=0))/length(migration.predictors$Sample)),
                     direction=c(
                       
                      median(migration.predictors[which(migration.predictors[,"greatCircleDistances"]!=0),"greatCircleDistances"])>0,
                      median(migration.predictors[which(migration.predictors[,"national_border_shared"]!=0),"national_border_shared"])>0,
                      median(migration.predictors[which(migration.predictors[,"originTT100k"]!=0),"originTT100k"])>0,
                      median(migration.predictors[which(migration.predictors[,"destinationTT100k"]!=0),"destinationTT100k"])>0,
                      median(migration.predictors[which(migration.predictors[,"originGecon"]!=0),"originGecon"])>0,
                      median(migration.predictors[which(migration.predictors[,"destinationGecon"]!=0),"destinationGecon"])>0,
                      median(migration.predictors[which(migration.predictors[,"originPdens"]!=0),"originPdens"])>0,
                      median(migration.predictors[which(migration.predictors[,"destinationPdens"]!=0),"destinationPdens"])>0,                                                                                                                                                                        
                      median(migration.predictors[which(migration.predictors[,"originPop_Size"]!=0),"originPop_Size"])>0,
                      median(migration.predictors[which(migration.predictors[,"destinationPop_Size"]!=0),"destinationPop_Size"])>0,
                      median(migration.predictors[which(migration.predictors[,"originTemp"]!=0),"originTemp"])>0,
                      median(migration.predictors[which(migration.predictors[,"destinationTemp"]!=0),"destinationTemp"])>0,
                      median(migration.predictors[which(migration.predictors[,"originTmpss"]!=0),"originTmpss"])>0,
                      median(migration.predictors[which(migration.predictors[,"destinationTmpss"]!=0),"destinationTmpss"])>0,
                      median(migration.predictors[which(migration.predictors[,"originPrec"]!=0),"originPrec"])>0,
                      median(migration.predictors[which(migration.predictors[,"destinationPrec"]!=0),"destinationPrec"])>0,
                      median(migration.predictors[which(migration.predictors[,"originPrecss"]!=0),"originPrecss"])>0,
                      median(migration.predictors[which(migration.predictors[,"destinationPrecss"]!=0),"destinationPrecss"])>0,
                      median(migration.predictors[which(migration.predictors[,"Kailahun_from"]!=0),"Kailahun_from"])>0,
                      median(migration.predictors[which(migration.predictors[,"Kenema_from"]!=0),"Kenema_from"])>0,
                      median(migration.predictors[which(migration.predictors[,"Kono_from"]!=0),"Kono_from"])>0,
                      median(migration.predictors[which(migration.predictors[,"Bombali_from"]!=0),"Bombali_from"])>0,
                      median(migration.predictors[which(migration.predictors[,"Kambia_from"]!=0),"Kambia_from"])>0,
                      median(migration.predictors[which(migration.predictors[,"Koinadugu_from"]!=0),"Koinadugu_from"])>0,
                      median(migration.predictors[which(migration.predictors[,"PortLoko_from"]!=0),"PortLoko_from"])>0,
                      median(migration.predictors[which(migration.predictors[,"Tonkolili_from"]!=0),"Tonkolili_from"])>0,
                      median(migration.predictors[which(migration.predictors[,"Bo_from"]!=0),"Bo_from"])>0,
                      median(migration.predictors[which(migration.predictors[,"Bonthe_from"]!=0),"Bonthe_from"])>0,
                      median(migration.predictors[which(migration.predictors[,"Moyamba_from"]!=0),"Moyamba_from"])>0,
                      median(migration.predictors[which(migration.predictors[,"Pujehun_from"]!=0),"Pujehun_from"])>0,
                      median(migration.predictors[which(migration.predictors[,"WesternRural_from"]!=0),"WesternRural_from"])>0,
                      median(migration.predictors[which(migration.predictors[,"WesternUrban_from"]!=0),"WesternUrban_from"])>0,
                      median(migration.predictors[which(migration.predictors[,"Kailahun_to"]!=0),"Kailahun_to"])>0,
                      median(migration.predictors[which(migration.predictors[,"Kenema_to"]!=0),"Kenema_to"])>0,
                      median(migration.predictors[which(migration.predictors[,"Kono_to"]!=0),"Kono_to"])>0,
                      median(migration.predictors[which(migration.predictors[,"Bombali_to"]!=0),"Bombali_to"])>0,
                      median(migration.predictors[which(migration.predictors[,"Kambia_to"]!=0),"Kambia_to"])>0,
                      median(migration.predictors[which(migration.predictors[,"Koinadugu_to"]!=0),"Koinadugu_to"])>0,
                      median(migration.predictors[which(migration.predictors[,"PortLoko_to"]!=0),"PortLoko_to"])>0,
                      median(migration.predictors[which(migration.predictors[,"Tonkolili_to"]!=0),"Tonkolili_to"])>0,
                      median(migration.predictors[which(migration.predictors[,"Bo_to"]!=0),"Bo_to"])>0,
                      median(migration.predictors[which(migration.predictors[,"Bonthe_to"]!=0),"Bonthe_to"])>0,
                      median(migration.predictors[which(migration.predictors[,"Moyamba_to"]!=0),"Moyamba_to"])>0,
                      median(migration.predictors[which(migration.predictors[,"Pujehun_to"]!=0),"Pujehun_to"])>0,
                      median(migration.predictors[which(migration.predictors[,"WesternRural_to"]!=0),"WesternRural_to"])>0,
                      median(migration.predictors[which(migration.predictors[,"WesternUrban_to"]!=0),"WesternUrban_to"])>0,
                      median(migration.predictors[which(migration.predictors[,"Samples_from"]!=0),"Samples_from"])>0,
                      median(migration.predictors[which(migration.predictors[,"Samples_to"]!=0),"Samples_to"])>0,
                      median(migration.predictors[which(migration.predictors[,"incidenceRatio"]!=0),"incidenceRatio"])>0))
p_mig <- ggplot(dat.mig) + geom_bar(aes(x=x,y=y,fill=direction), stat="identity", colour="black") +
  geom_hline(yintercept =  0.3947, color="black", linetype="dashed") + 
  xlab("") +
  scale_fill_manual(values=c(col1,col3), guide=F)+
  ylab("inclusion probability") +
  ggtitle("migration predictors") +
  scale_x_discrete(limits = rev(c("great circle distance", 
                              "shared border", 
                              "time to settlement with 100k at origin", 
                              "time to settlement with 100k at destination", 
                              "gridded economic output at origin", 
                              "gridded economic output at destination", 
                              "population density at origin", 
                              "population density at destination", 
                              "population size at origin", 
                              "population size at destination", 
                              "mean temperature at origin", 
                              "mean temperature at destination", 
                              "temperature seasonality at origin", 
                              "temperature seasonality at destination", 
                              "mean precipitation at origin", 
                              "mean precipitation at destination", 
                              "precipitation seasonality at origin", 
                              "precipitation seasonality at destination", 
                              "Kailahun from",
                              "Kenema from",
                              "Kono from",
                              "Bombali from",
                              "Kambia from",
                              "Koinadugu from",
                              "PortLoko from",
                              "Tonkolili from",
                              "Bo from",
                              "Bonthe from",
                              "Moyamba from",
                              "Pujehun from",
                              "Western Rural from",
                              "Western Urban from",
                              "Kailahun to",
                              "Kenema to",
                              "Kono to",
                              "Bombali to",
                              "Kambia to",
                              "Koinadugu to",
                              "PortLoko to",
                              "Tonkolili to",
                              "Bo to",
                              "Bonthe to",
                              "Moyamba to",
                              "Pujehun to",
                              "Western Rural to",
                              "Western Urban to",
                              "samples at origin",
                              "samples at destination",
                              "ratio of weekly cases"))) + coord_flip()

  coord_flip()

# plot(p_mig)
# ggsave(plot=p_mig,paste("Ebola_Migration_predictors.pdf", sep=""),width=4, height=8)
  
hlay <- rbind(c(1,2),
                c(1,2),
                c(NA,2),
                c(NA,2),
                c(NA,2),
                c(NA,2))
hlay <- rbind(c(1,1,1,1,2,2,2,2,2))
  
  
p_all <- grid.arrange(p_ne, p_mig, layout_matrix=hlay)
ggsave(plot= p_all ,paste("../../Figures/Ebola/Ebola_predictors.pdf", sep=""),width=8, height=8)
