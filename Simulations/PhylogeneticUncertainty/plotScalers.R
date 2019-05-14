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
library("colorblindr")


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
first = TRUE
firstCI = TRUE

# Read In Data ---------------------------------
for (i in seq(1,length(log),1)){
  print(i)
  # Make the filenames for all the three runs
  filename1 <- paste(log[i], sep="")
  
  
  # Read in the log files of the runs with fixed trees
  t <- read.table(filename1, header=TRUE, sep="\t")
  t <- t[-seq(1,ceiling(length(t$m1)/10)), ]
  # Read in the log files of the runs w/o fixed trees
  t.phylo <- read.table(gsub("/out", "/treeinfout", filename1), header=TRUE, sep="\t")
  t.phylo <- t.phylo[-seq(1,ceiling(length(t.phylo$m1)/10)), ]
  
  
  # calculate ess values
  ess <- effectiveSize(t)
  ess.phylo <- effectiveSize(t.phylo)
  
  
  
  if (min(ess[2:3])<100 || min(ess.phylo[2:3])<100){
    print("masco ESS value to low")
    print(sprintf("ESS value is %f for file %s",min(ess.phylo[2:3]),filename1))
  }else{
    dfname <- data.frame(filename = filename1)
    
    # get the number of the run
    tmp = strsplit(filename1, split="_")
    tmp = strsplit(tmp[[1]][2], split="S")
    
    runnumber = as.numeric(tmp[[1]][2])
    
    Ne_index = which(Ne_params$Nr==runnumber)
    m_index = which(m_params$Nr==runnumber)
    
    active.ne = 0;
    # loop over the Ne's
    for (i in seq(2,length(Ne_params))){
      name = names(Ne_params)[i]
      
      lower = quantile(t[which(t[,name]!=0), name],0.025)
      upper = quantile(t[which(t[,name]!=0), name],0.975)
      
      if (Ne_params[Ne_index,i] >=lower && Ne_params[Ne_index,i] <= upper){
        inint=1
      }else{
        inint=0
      }
      
      lower.phylo=quantile(t.phylo[which(t.phylo[,name]!=0), name],0.025)
      upper.phylo=quantile(t.phylo[which(t.phylo[,name]!=0), name],0.975)
      
      if (is.na(lower.phylo))lower.phylo=0.0
      if (is.na(upper.phylo))upper.phylo=0.0
      
      if (Ne_params[Ne_index,i] >=lower.phylo && Ne_params[Ne_index,i] <= upper.phylo){
        inint.phylo=1
      }else{
        inint.phylo=0
      }
      
      
      Ne.new <- data.frame(true=Ne_params[Ne_index,i], 
                           est=median(t[which(t[,name]!=0), name]), lower=lower, upper=upper, inc=length(which(t[,name]!=0))/length(t[,name]), inint=inint,
                           est.phylo=median(t.phylo[which(t.phylo[,name]!=0), name]), lower.phylo=lower.phylo, upper.phylo=upper.phylo, inc.phylo=length(which(t.phylo[,name]!=0))/length(t.phylo[,name]), inint.phylo=inint.phylo)
      
      if (Nenr==1){
        Ne <- Ne.new
        Nenr=Nenr+1
      }else{
        Ne <- rbind(Ne, Ne.new)
      }
      
      if (Ne_params[Ne_index,i]!=0){
        active.ne = active.ne+1;
      }
      
    }
    
    
    active.mig = 0;
    # loop over the Ne's
    for (i in seq(2,length(m_params))){
      name = names(m_params)[i]
      
      lower = quantile(t[which(t[,name]!=0), name],0.025)
      upper = quantile(t[which(t[,name]!=0), name],0.975)
      
      if (m_params[m_index,i] >=lower && m_params[m_index,i] <= upper){
        inint=1
      }else{
        inint=0
      }
      
      lower.phylo=quantile(t.phylo[which(t.phylo[,name]!=0), name],0.025)
      upper.phylo=quantile(t.phylo[which(t.phylo[,name]!=0), name],0.975)
      
      if (m_params[m_index,i] >=lower.phylo && m_params[m_index,i] <= upper.phylo){
        inint.phylo=1
      }else{
        inint.phylo=0
      }
      
      
      m.new <- data.frame(true=m_params[m_index,i], 
                          est=median(t[which(t[,name]!=0), name]), lower=lower, upper=upper, inc=length(which(t[,name]!=0))/length(t[,name]), inint=inint,
                          est.phylo=median(t.phylo[which(t.phylo[,name]!=0), name]), lower.phylo=lower.phylo, upper.phylo=upper.phylo, inc.phylo=length(which(t.phylo[,name]!=0))/length(t.phylo[,name]),inint.phylo=inint.phylo)
      if (mnr==1){
        m <- m.new
        mnr=mnr+1
      }else{
        m <- rbind(m, m.new)
      }
      
      if (m_params[m_index,i]!=0){
        active.mig = active.mig+1;
      }
    }
    
    
    confidence_vals = seq(0.05, 0.95, 0.05)
      
    for (j in seq(1,length(confidence_vals))){
      
      hpd.int = HPDinterval(as.mcmc(t), prob=1-confidence_vals[[j]])
      hpd.int.phylo = HPDinterval(as.mcmc(t.phylo), prob=1-confidence_vals[[j]])
      

      lower.mig = hpd.int["sum.migrationIndicator.","lower"]
      upper.mig = hpd.int["sum.migrationIndicator.","upper"]
      intsize.mig = length(which(t[,"sum.migrationIndicator."]>=lower.mig & t[,"sum.migrationIndicator."]<=upper.mig))/length(t$Sample)
      
      
      lower.mig.phylo = hpd.int.phylo["sum.migrationIndicator.","lower"]
      upper.mig.phylo = hpd.int.phylo["sum.migrationIndicator.","upper"]
      intsize.mig.phylo = length(which(t.phylo[,"sum.migrationIndicator."]>=lower.mig.phylo & t.phylo[,"sum.migrationIndicator."]<=upper.mig.phylo))/length(t.phylo$Sample)
      
      
      lower.ne = hpd.int["sum.NeIndicator.","lower"]
      upper.ne = hpd.int["sum.NeIndicator.","upper"]
      intsize.ne = length(which(t[,"sum.NeIndicator."]>=lower.ne & t[,"sum.NeIndicator."]<=upper.ne))/length(t$Sample)
      
      
      lower.ne.phylo = hpd.int.phylo["sum.NeIndicator.","lower"]
      upper.ne.phylo = hpd.int.phylo["sum.NeIndicator.","upper"]
      intsize.ne.phylo = length(which(t.phylo[,"sum.NeIndicator."]>=lower.ne.phylo & t.phylo[,"sum.NeIndicator."]<=upper.ne.phylo))/length(t.phylo$Sample)
      

  
      inint.mig=0
      inint.mig.phylo=0
      inint.ne=0
      inint.ne.phylo=0 
      
      if (active.mig >=lower.mig && active.mig <= upper.mig){inint.mig=1}
      if (active.mig >=lower.mig.phylo && active.mig <= upper.mig.phylo){inint.mig.phylo=1}
      if (active.ne >=lower.ne && active.ne <= upper.ne){inint.ne=1}
      if (active.ne >=lower.ne.phylo && active.ne <= upper.ne.phylo){inint.ne.phylo=1}
      
      new.sumactive <- data.frame(inint.mig=inint.mig,inint.mig.phylo=inint.mig.phylo,inint.ne=inint.ne,inint.ne.phylo=inint.ne.phylo,
                                  intsize.mig=intsize.mig, intsize.mig.phylo=intsize.mig.phylo, intsize.ne=intsize.ne,intsize.ne.phylo=intsize.ne.phylo,
                                  confidence_val=confidence_vals[[j]])
      
      if(firstCI){sumactive=new.sumactive; firstCI=F}
      else{sumactive=rbind(sumactive,new.sumactive)}
    }
    

  }
}

sumactive_tmp = sumactive[which(sumactive$confidence_val<0.25),]

p_conf <- ggplot(data = sumactive)+
  geom_smooth(aes(x=intsize.mig,y=inint.mig, color="sum active migration predictors fixed tree"), method="lm", formula=y~x)+
  geom_smooth(aes(x=intsize.mig.phylo,y=inint.mig.phylo, color="sum active migration predictors inferred tree"), method="lm", formula=y~x)+
  geom_smooth(aes(x=intsize.ne,y=inint.ne, color="sum active Ne predictors fixed tree"), method="lm", formula=y~x)+
  geom_smooth(aes(x=intsize.ne.phylo,y=inint.ne.phylo, color="sum active Ne predictors inferred tree"), method="lm", formula=y~x) +
  xlab("HPD interval size") +
  ylab("coverage") +
  theme_minimal()

plot(p_conf)


ggsave(plot=p_conf,"../../Figures/Phylo/Phylo_sumactive.pdf",width=7, height=3.5)


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

ggsave(plot=p_mig,"../../Figures/Phylo/Phylo_migration_fixed.pdf",width=3.5, height=3.5)


p_mig <- ggplot()+
  geom_point(data=m_red1, aes(x=true, y=est.phylo))+
  geom_segment(data=m_red1, aes(x = -3.5, y = -3.5, xend = 3.5, yend = 3.5), color="red") +
  ylab("estimated") + xlab("true") + ggtitle("coefficients") + 
  theme(legend.position="none")
plot(p_mig)

ggsave(plot=p_mig,"../../Figures/Phylo/Phylo_migration_inferred.pdf",width=3.5, height=3.5)



p_mig_inc1 <- ggplot()+
  geom_point(data=m_red1, aes(x=abs(true), y=inc), size=1)+
  ylab("inclusion probability") + xlab("absolute value of true scaler") + ggtitle("Inclusion of active covariates") + 
  geom_hline(yintercept = 0.4268, color="black", linetype="dashed") 
plot(p_mig_inc1)
ggsave(plot=p_mig_inc1,"../../Figures/Phylo/Phylo_mig_inclusion_fixed.pdf",width=3.5, height=3.5)


p_mig_inc1 <- ggplot()+
  geom_point(data=m_red1, aes(x=abs(true), y=inc.phylo), size=1)+
  ylab("inclusion probability") + xlab("absolute value of true scaler") + ggtitle("Inclusion of active covariates") + 
  geom_hline(yintercept = 0.4268, color="black", linetype="dashed") 
plot(p_mig_inc1)
ggsave(plot=p_mig_inc1,"../../Figures/Phylo/Phylo_mig_inclusion_inferred.pdf",width=3.5, height=3.5)





# p_mig_inc1 <- ggplot()+
#   geom_point(data=m_red1, aes(x=abs(true), y=inc), size=2)+
#   geom_point(data=m_red1, aes(x=abs(true), y=inc), size=1)+
#   ylab("inclusion probability") + xlab("absolute value of true scaler") + ggtitle("Inclusion of active covariates") + 
#   geom_hline(yintercept = 0.6383, color="black", linetype="dashed") + 
#   plot(p_mig_inc1)
# 
# 
# p_mig_inc2 <- ggplot()+
#   geom_violin(data=m_red2, aes(x="fixed tree", y=inc))+
#   geom_violin(data=m_red2, aes(x="inferred tree", y=inc.phylo))+
#   ylab("inclusion probability") + xlab("true") + ggtitle("Exclusion of inactive covariates") + 
#   theme_minimal()+
#   # theme(legend.position="none") + 
#   scale_color_OkabeIto()
# plot(p_mig_inc2)


#geom_point(data=m, aes(x=true,y=est), size=0.001, alpha=0.1) +

Ne_red1 = Ne[which(Ne$true!=0),]
Ne_red2 = Ne[which(Ne$true==0),]


p_Ne <- ggplot()+
  geom_point(data=Ne_red1, aes(x=true, y=est))+
  geom_segment(data=Ne_red1, aes(x = -3.5, y = -3.5, xend = 3.5, yend = 3.5), color="red") +
  ylab("estimated") + xlab("true") + ggtitle("coefficients") + 
  theme(legend.position="none")
plot(p_Ne)
ggsave(plot=p_Ne,"../../Figures/Phylo/Phylo_Ne_fixed.pdf",width=3.5, height=3.5)

p_Ne <- ggplot()+
  geom_point(data=Ne_red1, aes(x=true, y=est.phylo))+
  geom_segment(data=Ne_red1, aes(x = -3.5, y = -3.5, xend = 3.5, yend = 3.5), color="red") +
  ylab("estimated") + xlab("true") + ggtitle("coefficients") + 
  theme(legend.position="none")
plot(p_Ne)
ggsave(plot=p_Ne,"../../Figures/Phylo/Phylo_Ne_inferred.pdf",width=3.5, height=3.5)



p_Ne_inc1 <- ggplot()+
  geom_point(data=Ne_red1, aes(x=abs(true), y=inc), size=1)+
  ylab("inclusion probability") + xlab("absolute value of true scaler") + ggtitle("Inclusion of active covariates") + 
  geom_hline(yintercept = 0.4268, color="black", linetype="dashed") + 
  theme(legend.position="none") 
plot(p_Ne_inc1)
ggsave(plot=p_Ne_inc1,"../../Figures/Phylo/Phylo_Ne_inclusion_fixed.pdf",width=3.5, height=3.5)

p_Ne_inc1 <- ggplot()+
  geom_point(data=Ne_red1, aes(x=abs(true), y=inc.phylo), size=1)+
  ylab("inclusion probability") + xlab("absolute value of true scaler") + ggtitle("Inclusion of active covariates") + 
  geom_hline(yintercept = 0.4268, color="black", linetype="dashed") + 
  theme(legend.position="none") 
plot(p_Ne_inc1)
ggsave(plot=p_Ne_inc1,"../../Figures/Phylo/Phylo_Ne_inclusion_inferred.pdf",width=3.5, height=3.5)


# p_Ne_inc2 <- ggplot()+
#   geom_violin(data=Ne_red2, aes(x=true, y=inc))+
#   ylab("inclusion probability") + xlab("true") + ggtitle("Exclusion of inactive covariates") + 
#   theme(legend.position="none")
# plot(p_Ne_inc2)

p_m_inc <- ggplot()+
  geom_point(data=m_red1, aes(x=inc, y=inc.phylo, color=abs(true)),size=2)+
  geom_segment(data=m_red1, aes(x = 0, y = 0, xend =1, yend = 1), color="red") +
  ylab("inclusion probability inferred tree") +
  xlab("inclusion probability fixed tree") +
  ggtitle("migration rates")  +
  scale_colour_gradientn(values=c(0,0.5,1,1.5,2,2.5), breaks=c(0,0.5,1,1.5,2,2.5),  colours=c("white", "black", "black", "black", "black") , na.value="white")
plot(p_m_inc)
ggsave(plot=p_m_inc,"../../Figures/Phylo/Phylo_migration_inclusion.pdf",width=5, height=3.5)


p_Ne_inc <- ggplot()+
  geom_point(data=Ne_red1, aes(x=inc, y=inc.phylo, color=abs(true)),size=2)+
  geom_segment(data=Ne_red1, aes(x = 0, y = 0, xend =1, yend = 1), color="red") +
  ylab("inclusion probability inferred tree") +
  xlab("inclusion probability fixed tree") +
  ggtitle("effective population size") + 
  scale_colour_gradientn(values=c(0,0.5,1,1.5,2,2.5), breaks=c(0,0.5,1,1.5,2,2.5),  colours=c("white", "black", "black", "black", "black") , na.value="white")
plot(p_Ne_inc)
ggsave(plot=p_Ne_inc,"../../Figures/Phylo/Phylo_Ne_inclusion.pdf",width=5, height=3.5)


# ggsave(plot=p_mig,"../../Figures/Stepwise/Stepwise_migration.pdf",width=3.5, height=3.5)
# 
# ggsave(plot=p_mig_inc1,"../../Figures/Stepwise/Stepwise_mig_inclusion.pdf",width=3.5, height=3.5)
# ggsave(plot=p_mig_inc2,"../../Figures/Stepwise/Stepwise_mig_exclusion.pdf",width=3.5, height=3.5)
# 
# 
# ggsave(plot=p_Ne,"../../Figures/Stepwise/Stepwise_ne.pdf",width=3.5, height=3.5)
# 
# ggsave(plot=p_Ne_inc1,"../../Figures/Stepwise/Stepwise_Ne_inclusion.pdf",width=3.5, height=3.5)
# ggsave(plot=p_Ne_inc2,"../../Figures/Stepwise/Stepwise_Ne_exclusion.pdf",width=3.5, height=3.5)


# print(sprintf("coverage migration = %f ", mean(cov.m$isIn)))
# print(sprintf("coverage Ne = %f ",mean(cov.Ne$isIn)))


p_treecomp <- ggplot()+
  geom_point(data=m_red1, aes(x=est, y=est.phylo),size=2)+
  geom_errorbar(data=m_red1, aes(x=est,ymin=lower.phylo, ymax=upper.phylo)) +
  geom_errorbarh(data=m_red1, aes(y=est.phylo, xmin=lower,xmax=upper))+
  geom_segment(data=m_red1, aes(x = -3.5, y = -3.5, xend = 3.5, yend = 3.5), color="red") +
  ylab(paste("scaler inferred tree", "(coverage =", format(sum(m_red1$inint.phylo)/length(m_red1$inint)*100,digits=1), "%)")) +
  xlab(paste("fixed inferred tree", "(coverage =", format(sum(m_red1$inint)/length(m_red1$inint)*100,digits=1), "%)")) +
  ggtitle("migration rates") + 
  # theme_minimal()+
  theme(legend.position="none") 
plot(p_treecomp)
ggsave(plot=p_treecomp,"../../Figures/Phylo/Phylo_migration_treecomp.pdf",width=3.5, height=3.5)




# p_error <- ggplot()+
#   geom_point(data=m_red1, aes(x=true, y=est),size=1)+
#   geom_errorbar(data=m_red1, aes(x=true, ymin=lower, ymax=upper))+
#   geom_segment(data=m_red1, aes(x = -3.5, y = -3.5, xend = 3.5, yend = 3.5), color="red") +
#   # facet_grid(method~.)+
#   ylab("estimated") + xlab("true")+ 
#   ggtitle("migration rates") + 
#   theme(legend.position="none")
# plot(p_error)
# ggsave(plot=p_error,"../../Figures/Phylo/Phylo_migration_error.pdf",width=3.5, height=3.5)


p_Ne_treecomp <- ggplot()+
  geom_point(data=Ne_red1, aes(x=est, y=est.phylo),size=2)+
  geom_errorbar(data=Ne_red1, aes(x=est,ymin=lower.phylo, ymax=upper.phylo)) +
  geom_errorbarh(data=Ne_red1, aes(y=est.phylo, xmin=lower, xmax=upper))+
  geom_segment(data=Ne_red1, aes(x = -3.5, y = -3.5, xend = 3.5, yend = 3.5)) +
  ylab(paste("scaler inferred tree", "(coverage =", format(sum(Ne_red1$inint.phylo)/length(Ne_red1$inint)*100,digits=1), "%)")) +
  xlab(paste("fixed inferred tree", "(coverage =", format(sum(Ne_red1$inint)/length(Ne_red1$inint)*100,digits=1), "%)")) +
  ggtitle("effective population size") + 
  # theme_minimal()+
  theme(legend.position="none")
plot(p_Ne_treecomp)
ggsave(plot=p_Ne_treecomp,"../../Figures/Phylo/Phylo_Ne_treecomp.pdf",width=3.5, height=3.5)


# p_Ne_error <- ggplot()+
#   geom_point(data=Ne_red1, aes(x=true, y=est, color=method),size=0.1)+
#   geom_errorbar(data=Ne_red1, aes(x=true, ymin=lower, ymax=upper, color=method))+
#   geom_segment(data=Ne_red1, aes(x = -3.5, y = -3.5, xend = 3.5, yend = 3.5), color="red") +
#   # facet_grid(method~.)+
#   ylab("estimated") + xlab("true")+ 
#   ggtitle("effective population size") + 
#   theme(legend.position="none")
# plot(p_Ne_error)
# ggsave(plot=p_Ne_error,"../../Figures/Phylo/Phylo_Ne_error.pdf",width=3.5, height=3.5)
