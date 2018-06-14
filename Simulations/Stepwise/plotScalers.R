######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
# needed to calculate ESS values
library(coda)
library("methods")


# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# get the names of all SISCO first (of three) run log files
log <- list.files(path="./out", pattern="*.log", full.names = TRUE)

# use the matlab standard colors to plot
col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)
col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)


states <- 5

m_params <- read.table("mParameters.txt", header=TRUE, sep="\t", nrows = 0)
Ne_params <- read.table("NeParameters.txt", header=TRUE, sep="\t", nrows = 0)

Nenr <- 1
mnr <- 1

# Read In Data ---------------------------------
for (i in seq(1,length(log),1)){
  print(i)
  # Make the filenames for all the three runs
  filename1 <- paste(log[i], sep="")

  
  # Read in the SISCO *.logs
  t <- read.table(filename1, header=TRUE, sep="\t")
  t <- t[-seq(1,ceiling(length(t$m1)/10)), ]

  # calculate ess values
  ess <- effectiveSize(t)

    
  if (min(ess[2:3])<100){
    print("masco ESS value to low")
    print(sprintf("ESS value is %f for file %s",min(ess[2:6]),filename1))
  }else{
      dfname <- data.frame(filename = filename1)
      
      # get the number of the run
      tmp = strsplit(filename1, split="_")
      tmp = strsplit(tmp[[1]][2], split="S")
      
      runnumber = as.numeric(tmp[[1]][2])
      
      Ne_index = which(Ne_params$Nr==runnumber)
      m_index = which(m_params$Nr==runnumber)
      
      # loop over the Ne's
      for (i in seq(2,length(Ne_params))){
          name = names(Ne_params)[i]
          Ne.new <- data.frame(true=Ne_params[Ne_index,i], est=median(t[which(t[,name]!=0), name]), inc=length(which(t[,name]!=0))/length(t[,name]))
          if (Nenr==1){
            Ne <- Ne.new
            Nenr=Nenr+1
          }else{
            Ne <- rbind(Ne, Ne.new)
          }
      }
      # loop over the Ne's
      for (i in seq(2,length(m_params))){
        name = names(m_params)[i]
        m.new <- data.frame(true=m_params[m_index,i], est=median(t[which(t[,name]!=0), name]), inc=length(which(t[,name]!=0))/length(t[,name]))
        if (mnr==1){
          m <- m.new
          mnr=mnr+1
        }else{
          m <- rbind(m, m.new)
        }
      }
  }
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot the rate ratios
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


m_red1 = m[which(m$true!=0),]
m_red2 = m[which(m$true==0),]

p_mig <- ggplot()+
  geom_point(data=m_red1, aes(x=true, y=est))+
  geom_segment(data=m_red1, aes(x = -3.5, y = -3.5, xend = 3.5, yend = 3.5), color="red") +
  ylab("estimated") + xlab("true") + ggtitle("coefficients") + 
  theme(legend.position="none")
plot(p_mig)



p_mig_inc1 <- ggplot()+
  geom_point(data=m_red1, aes(x=abs(true), y=inc), size=1)+
  ylab("inclusion probability") + xlab("absolute value of true scaler") + ggtitle("Inclusion of active covariates") + 
  geom_hline(yintercept = 0.6818, color="black", linetype="dashed") + 
  theme(legend.position="none")
plot(p_mig_inc1)

p_mig_inc2 <- ggplot()+
  geom_violin(data=m_red2, aes(x=true, y=inc))+
  ylab("inclusion probability") + xlab("true") + ggtitle("Exclusion of inactive covariates") + 
  theme(legend.position="none")
plot(p_mig_inc2)


#geom_point(data=m, aes(x=true,y=est), size=0.001, alpha=0.1) +

Ne_red1 = Ne[which(Ne$true!=0),]
Ne_red2 = Ne[which(Ne$true==0),]


p_Ne <- ggplot()+
  geom_point(data=Ne_red1, aes(x=true, y=est))+
  geom_segment(data=Ne_red1, aes(x = -3.5, y = -3.5, xend = 3.5, yend = 3.5), color="red") +
  ylab("estimated") + xlab("true") + ggtitle("coefficients") + 
  theme(legend.position="none")
plot(p_Ne)


p_Ne_inc1 <- ggplot()+
  geom_point(data=Ne_red1, aes(x=abs(true), y=inc), size=1)+
  ylab("inclusion probability") + xlab("absolute value of true scaler") + ggtitle("Inclusion of active covariates") + 
  geom_hline(yintercept = 0.6818, color="black", linetype="dashed") + 
  theme(legend.position="none") 
plot(p_Ne_inc1)

p_Ne_inc2 <- ggplot()+
  geom_violin(data=Ne_red2, aes(x=true, y=inc))+
  ylab("inclusion probability") + xlab("true") + ggtitle("Exclusion of inactive covariates") + 
  theme(legend.position="none")
plot(p_Ne_inc2)


ggsave(plot=p_mig,"../../text/figures/Stepwise_migration.pdf",width=3.5, height=3.5)

ggsave(plot=p_mig_inc1,"../../text/figures/Stepwise_mig_inclusion.pdf",width=3.5, height=3.5)
ggsave(plot=p_mig_inc2,"../../text/figures/Stepwise_mig_exclusion.pdf",width=3.5, height=3.5)


ggsave(plot=p_Ne,"../../text/figures/Stepwise_ne.pdf",width=3.5, height=3.5)

ggsave(plot=p_Ne_inc1,"../../text/figures/Stepwise_Ne_inclusion.pdf",width=3.5, height=3.5)
ggsave(plot=p_Ne_inc2,"../../text/figures/Stepwise_Ne_exclusion.pdf",width=3.5, height=3.5)


print(sprintf("coverage migration = %f ", mean(cov.m$isIn)))
print(sprintf("coverage Ne = %f ",mean(cov.Ne$isIn)))
