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



# prevents files hanging 
system("rm combined/*.*")
# get the names of all output files of the first replicate
log <- list.files(path="./out/", pattern=paste("*rep0.*", sep=""), full.names = TRUE)
for (i in seq(1, length(log))){
  # for (i in seq(1, 9)){
  
  print(log[i])
  if (grepl("trees", log[[i]]) || grepl("log", log[[i]])){
    in_command <- " -b 50 -log"
    for (j in seq(0,2)){
      in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), log[i]), sep="")
    }
    
    out_command = gsub("rep0_", "", log[i])
    out_command = gsub("out", "combined", out_command)
    
    combined_command = out_command
    combined_command = paste(" -o ", gsub("_rep0", "",combined_command), sep="")
    inc = gsub("-o ","",combined_command)
    
    # combine the trees
    if (grepl("node", log[[i]])){
      system(intern=T,ignore.stdout = TRUE, ignore.stderr = TRUE, paste("/Applications/BEAST\\ 2.6.7/bin/logcombiner", in_command, combined_command, sep=" "))
      system(intern=T, paste("/Applications/BEAST\\ 2.6.7/bin/treeannotator", inc, gsub(".trees",".tree",inc), sep=" "))
    }else if (grepl("log", log[[i]])){
      system(intern=T,ignore.stdout = TRUE, ignore.stderr = TRUE, paste("/Applications/BEAST\\ 2.6.7/bin/logcombiner", in_command, combined_command, sep=" "))
    }else{
      # system(intern=T,ignore.stdout = TRUE, ignore.stderr = TRUE, paste("/Applications/BEAST\\ 2.6.7/bin/logcombiner", in_command, combined_command, sep=" "))
    }
  }
}




# get state clock runs
log <- list.files(path="./combined", pattern=".*\\.log", full.names = TRUE)

dat = data.frame()
for (i in seq(1, length(log))){
  t = read.table(log[[i]], header=T, sep="\t")
  
  print(effectiveSize(as.mcmc(t$posterior)))
  print(length(t$posterior))
  if (effectiveSize(as.mcmc(t$posterior))>5){
    
    name = gsub("./tmp/", "", strsplit(log[[i]], split="_")[[1]][[1]])
    states = strsplit(log[[i]], split="_")[[1]][[3]]
    migration = strsplit(log[[i]], split="_")[[1]][[4]]
    clock = paste(strsplit(log[[i]], split="_")[[1]][[6]], "clock model")
    rep = strsplit(log[[i]], split="_")[[1]][[5]]
    
    
    dat =rbind(dat, data.frame(x=t$birthRateSVCanonical.loc1/t$birthRateSVCanonical.loc0,
                               tumor = name,
                               states=states,
                               migration=migration,
                               clock=clock,
                               repetition=rep))
    }
}
p.rates.ratio <- ggplot(dat)+
  geom_density(aes(x=x, fill=clock, color=clock, linetype=repetition), alpha=0.4)+
  facet_wrap(migration~interaction(tumor, states), ncol=4)+
  theme_minimal() +
  scale_x_log10()
plot(p.rates.ratio)

ggsave(plot=p.rates.ratio, "ratios.pdf", width=6, height=6)

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

