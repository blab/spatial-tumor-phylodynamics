#setup_nextstrain_metadata.R
####### SETUP #######
library(tidyverse)
library(phangorn)
library(ape)
library(ggtree)
library(treeio)
library(viridis)
library(tumortree)

colors <- tumortree::get_color_palette()

t1_li_edge_labels <- data.frame("punch" = c("t1z5", "t1l13", "t1z1", "t1f24", "t1z3", "t1f23", "t1f11", "t1l8", "t1f14", "t1f9", "t1l1", "t1l10", "t1l3", "t1f2", "t1l6", "t1f4"),
                                "edgeP" = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)) 

t2_li_edge_labels <- data.frame("punch" = c("t2z11", "t2z13", "t2z1", "t2f9", "t2f5", "t2f2", "t2f13", "t2z9", "t2z6"),
                                 "edgeP" = c(1, 1, 1, 0, 0, 0, 0, 0, 0)) 


#reated for the Augur pipeline
t1_metadata <- t1_li_edge_labels %>% 
    mutate(punch_label = toupper(gsub("t1", "", punch))) %>% 
    rename("strain" = punch, "edge" = edgeP) %>% 
    add_row("strain" = "Normal")

t2_metadata <- t2_li_edge_labels %>% 
    mutate(punch_label = toupper(gsub("t2", "", punch))) %>% 
    rename("strain" = punch, "edge" = edgeP) %>% 
    add_row("strain" = "Normal")

write_tsv(t1_metadata, file = "../li-application/nextstrain_analysis/data/t1_metadata.tsv")
write_tsv(t2_metadata, file = "../li-application/nextstrain_analysis/data/t2_metadata.tsv")

