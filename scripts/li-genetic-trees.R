#li-genetic-trees.R
####### SETUP #######
library(tidyverse)
library(phangorn)
library(ape)
library(ggtree)
library(treeio)
library(viridis)
library(tumortree)

t1_li_edge_labels <- data.frame("punch" = c("t1z5", "t1l13", "t1z1", "t1f24", "t1z3", "t1f23", "t1f11", "t1l8", "t1f14", "t1f9", "t1l1", "t1l10", "t1l3", "t1f2", "t1l6", "t1f4"),
                                "edgeP" = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)) #%>% 
    # dplyr::mutate("Punch" = toupper(Punch))


t2_li_edge_labels <- data.frame("punch" = c("t2z11", "t2z13", "t2z1", "t2f9", "t2f5", "t2f2", "t2f13", "t2z9", "t2z6"),
                                 "edgeP" = c(1, 1, 1, 0, 0, 0, 0, 0, 0))  #%>% 
    # dplyr::mutate("Punch" = toupper(Punch))

# t1_subclone_colors <- c("alpha" = "#F4CAC7" , "alpha1" = "#F4ADAE", "alpha2" = "#EFA6B9", "beta" = "#F6E3B1", "beta1" = "#F7C0A2", "gamma"= "#B3D9B0", "other" = "grey")
# t1_subclone_labels <- data.frame("punch" = c("t1z5", "t1l13", "t1z1", "t1f24", "t1z3", "t1f23", "t1f11", "t1l8", "t1f14", "t1f9", "t1l1", "t1l10", "t1l3", "t1f2", "t1l6", "t1f4"),
#                                  "subclone" = c("gamma", "beta", "beta", "beta", "beta1", "beta1", "alpha2", "alpha2", "alpha2", "alpha2", "alpha1", "alpha1", "alpha1", "alpha1", "alpha1", "alpha1"))
# 
# t2_subclone_colors <- c("theta" = "#68AEDE","delta2"="#A3CDAA", "delta1" = "#8BC394", "other" = "#CBCACB")
# t2_subclone_labels <- data.frame("punch" = c("t2z11", "t2z13", "t2z1", "t2f9", "t2f5", "t2f2", "t2f13", "t2z9", "t2z6"),
#                                  "subclone" = c("other", "theta", "theta", "delta2", "delta2", "delta1", "delta1", "delta1", "delta1"))
####### MAXIMUM PARSIMONY TREE ######

#Fasta alignment files compiled in process_li_wgs_data.R
#VAF cutoff is 0.05 for variants
fasta_file_t1 <- "../li-application/li_hcc_data/li_t1_wgs.fa"
fasta_file_t2 <- "../li-application/li_hcc_data/li_t2_wgs.fa"

set.seed(38172)
#read in data into phyDat format (phangorn/ape)
t1_data <- phangorn::read.phyDat(file = fasta_file_t1, format = "fasta", type = "DNA")
t2_data <- phangorn::read.phyDat(file = fasta_file_t2, format = "fasta", type = "DNA")

#Manual topology input from Li et al

##Nwk strings
manual_t1_wgs_tree <- ape::read.tree(text = "((t1z5:37290,(((t1f11:459.5,(t1l8:304,(t1f14:329,t1f9:224)100:118)99.6:39.5)100:1706,(t1l1:1216,(t1l10:678,(t1l3:422.5,(t1f2:221,(t1l6:371,t1f4:201)100:67)100:83.5)100:349)100:438)100:680)100:16807.5,(t1z1:48483,(((t1f23:8645,t1z3:18876)100:22360,t1f24:45272)100:418.5,t1l13:4214.5)100:103)100:696.5)100:2566),blood:2301)100;")

manual_t2_wgs_tree <- ape::read.tree(text = "(((((((t2z6, t2z9), t2f13), t2f2), (t2f5, t2f9)), (t2z1, t2z13)), t2z11), blood);")

#Reconstruct branch lengths
t1_treePars  <- acctran(manual_t1_wgs_tree, t1_data)
t2_treePars  <- acctran(manual_t2_wgs_tree, t2_data)


t1_anc.acctran <- ancestral.pars(t1_treePars , t1_data, "ACCTRAN")
t2_anc.acctran <- ancestral.pars(t2_treePars , t2_data, "ACCTRAN")

#Put in treeio tree object

t1_treePars_obj <- treeio::as.treedata(t1_treePars)
t2_treePars_obj <- treeio::as.treedata(t2_treePars)

#Add data 
##Tumor 1
t1_treePars_obj@data <- tibble("node" = 1:(length(t1_treePars_obj@phylo$tip.label) + length(t1_treePars_obj@phylo$node.label)), 
                               "punch" = c(t1_treePars_obj@phylo$tip.label, rep(NA, length(t1_treePars_obj@phylo$node.label)))) %>% 
    
    dplyr::left_join(., t1_li_edge_labels, by = "punch") %>% 
    mutate(punch_label = toupper(gsub("t1", "", punch)))


t1_treePars_obj@data$punch_label[t1_treePars_obj@data$punch_label == "BLOOD"] <- "Blood"
##Tumor 2
t2_treePars_obj@data <- tibble("node" = 1:(length(t2_treePars_obj@phylo$tip.label) + length(t2_treePars_obj@phylo$node.label)), 
                               "punch" = c(t2_treePars_obj@phylo$tip.label, rep(NA, length(t2_treePars_obj@phylo$node.label)))) %>% 
    dplyr::left_join(., t2_li_edge_labels, by = "punch") %>% 
    mutate(punch_label = toupper(gsub("t2", "", punch)))


t2_treePars_obj@data$punch_label[t2_treePars_obj@data$punch_label == "BLOOD"] <- "Blood"
##trees where blood is a leaf
colors <- get_color_palette()
t1_wgs_tree <- ggtree(t1_treePars_obj, aes(color = ifelse(edgeP==1, "edge","center"), label = punch_label)) +
    scale_color_manual(values = colors) +
    geom_tippoint(size =2) +geom_text(color = "black", nudge_x=7000, size =3) +
    theme(legend.position = "none")
t1_wgs_tree 
ggsave(filename = "../figures/t1_wgs_genetic_state_tree.png",
       t1_wgs_tree, height = 4, width = 3)

t2_wgs_tree <- ggtree(t2_treePars_obj, aes(color = ifelse(edgeP==1, "edge","center"), label = punch_label)) +
    scale_color_manual(values = colors) +
    geom_tippoint(size =2) +
    #ggtitle("Li et al T2") +
    geom_text(color = "black", nudge_x=7000, size = 3) +
    theme(legend.position = "none")
t2_wgs_tree 
ggsave(filename = "../figures/t2_wgs_genetic_state_tree.png",
       t2_wgs_tree, height = 4, width = 3)
