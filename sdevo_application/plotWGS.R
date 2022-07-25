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
library(ape)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# get state clock runs
log <- list.files(path="./out", pattern=".*rep0\\.log", full.names = TRUE)

dat = data.frame()
for (i in seq(1, length(log))){
  t1 = read.table(log[[i]], header=T, sep="\t")
  t2 = read.table(gsub("rep0","rep1",log[[i]]), header=T, sep="\t")
  t3 = read.table(gsub("rep0","rep2",log[[i]]), header=T, sep="\t")
  
  t1 <- t1[-seq(1,ceiling(length(t1$Sample)/10)), ]
  t2 <- t2[-seq(1,ceiling(length(t2$Sample)/10)), ]
  t3 <- t3[-seq(1,ceiling(length(t3$Sample)/10)), ]
  
  t = rbind(t1,t2,t3)
  print(effectiveSize(as.mcmc(t$posterior)))
  print(length(t$posterior))
  if (effectiveSize(as.mcmc(t$posterior))>100){
    
    name = gsub("./out/", "", strsplit(log[[i]], split="_")[[1]][[1]])
    states = strsplit(log[[i]], split="_")[[1]][[3]]
    migration = strsplit(log[[i]], split="_")[[1]][[4]]
    clock = paste(strsplit(log[[i]], split="_")[[1]][[5]], "clock model")
    
    
    dat =rbind(dat, data.frame(x=t$birthRateSVCanonical.loc1/t$birthRateSVCanonical.loc0,
                               tumor = name,
                               states=states,
                               migration=migration,
                               clock=clock))
    }
  
}


p.rates.ratio <- ggplot(dat)+
  geom_density(aes(x=x, fill=clock, color=clock), alpha=0.6)+
  facet_wrap(migration~interaction(tumor, states), ncol=2)+
  theme_minimal() +
  scale_x_log10()

plot(p.rates.ratio)

ggsave(plot=p.rates.ratio, "ratios.pdf", width=6, height=6)

dsa
t1 = read.table("./out/T1_state_clocks_wgs_amplicon_300_sites_2_red.log", header=T, sep="\t")
t <- t1[-seq(1,ceiling(length(t1$Sample)/10)), ]

p.rates <- ggplot(t)+
  geom_density(aes(x=birthRateSVCanonical.loc0-deathRateSVCanonical.loc0+t$migrationRateSMCanonical.loc1_to_loc0, fill="birth rate in center"), alpha=0.6)+
  # geom_density(aes(x=birthRateSVCanonical.loc1, fill="birth rate on edge"), alpha=0.6) +
  geom_density(aes(x=birthRateSVCanonical.loc1-deathRateSVCanonical.loc1-t$migrationRateSMCanonical.loc1_to_loc0, fill="death rate in center"), alpha=0.6)+
  # geom_density(aes(x=deathRateSVCanonical.loc1, fill="death rate on edge"), alpha=0.6) +
  theme_minimal()
plot(p.rates)



p.rates.ratio <- ggplot(t)+
  geom_density(aes(x=birthRateSVCanonical.loc1/birthRateSVCanonical.loc0), fill="black", alpha=0.6)+
  theme_minimal()

plot(p.rates)
plot(p.rates.ratio)

dsa
ggsave(plot=p.rates, "../Figures/rates.pdf", width=9, height=9)



tr1 = read.table("./out/T1_state_clocks_wgs_amplicon_300_sites_2_red.HCCtumor.traj", header=T, sep="\t", colClasses=c("numeric", "character"))
tr <- tr1[-seq(1,ceiling(length(tr1$Sample)/2)), ]


trace = data.frame()
# for (i in seq(1,length(tr$Sample))){
for (i in seq(length(tr$Sample),length(tr$Sample))){
    
  tmp = strsplit(tr[i, "typedTraj.t.HCCtumor"], split=",")[[1]]
  new.trace = data.frame()
  for (j in seq(1,length(tmp), 50)){
    tmp2 = strsplit(tmp[[j]], split=":")[[1]]
    new.trace = rbind(new.trace, data.frame(time=as.numeric(tmp2[[1]]), center=as.numeric(tmp2[[6]]), edge=as.numeric(tmp2[[7]]), run=i))
  }
  tmp2 = strsplit(tmp[[length(tmp)]], split=":")[[1]]
  
  new.trace$time = new.trace$time-as.numeric(tmp2[[1]])
  trace = rbind(trace, new.trace)
}

ggplot(trace)+
  geom_line(aes(x=time, y=edge, group=run, color="edge"), alpha=1) +
  geom_line(aes(x=time, y=center, group=run, color="center"), alpha=1) +
  scale_y_log10()+
  theme_minimal()

