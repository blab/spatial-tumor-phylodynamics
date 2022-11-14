#li-genetic-trees.R
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

# t1_subclone_colors <- c("alpha" = "#F4CAC7" , "alpha1" = "#F4ADAE", "alpha2" = "#EFA6B9", "beta" = "#F6E3B1", "beta1" = "#F7C0A2", "gamma"= "#B3D9B0", "other" = "grey")
# t1_subclone_labels <- data.frame("punch" = c("t1z5", "t1l13", "t1z1", "t1f24", "t1z3", "t1f23", "t1f11", "t1l8", "t1f14", "t1f9", "t1l1", "t1l10", "t1l3", "t1f2", "t1l6", "t1f4"),
#                                  "subclone" = c("gamma", "beta", "beta", "beta", "beta1", "beta1", "alpha2", "alpha2", "alpha2", "alpha2", "alpha1", "alpha1", "alpha1", "alpha1", "alpha1", "alpha1"))
# 
# t2_subclone_colors <- c("theta" = "#68AEDE","delta2"="#A3CDAA", "delta1" = "#8BC394", "other" = "#CBCACB")
# t2_subclone_labels <- data.frame("punch" = c("t2z11", "t2z13", "t2z1", "t2f9", "t2f5", "t2f2", "t2f13", "t2z9", "t2z6"),
#                                  "subclone" = c("other", "theta", "theta", "delta2", "delta2", "delta1", "delta1", "delta1", "delta1"))
####### MAXIMUM PARSIMONY TREE ######

#Fasta alignment files compiled in process_li_wgs_data.R and write_ling_state_clocks.R
## Here these will only be used to convert divergence into number of genetic mutations
## VAF cutoff is 0.05 for variants
fasta_file_t1 <- "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/data/li_t1_wgs.fa"
fasta_file_t2 <- "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/data/li_t2_wgs.fa"
fasta_file_t3 <- "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/ling-application/data/ling-punch-sequences.fa"


##read in data into phyDat format (phangorn/ape)
t1_data <- phangorn::read.phyDat(file = fasta_file_t1, format = "fasta", type = "DNA")
t2_data <- phangorn::read.phyDat(file = fasta_file_t2, format = "fasta", type = "DNA")
t3_data <- phangorn::read.phyDat(file = fasta_file_t3, format = "fasta", type = "DNA")

t1_seq_length <- sum(attr(t1_data, "weight"))
t2_seq_length <- sum(attr(t2_data, "weight"))
t3_seq_length <- sum(attr(t3_data, "weight"))

######## MAXIMUM LIKELIHOOD TREE FROM AUGUR PIPELINE #################
max_likelihood_t1_wgs_tree <- ape::read.tree(file = "../li-application/nextstrain_analysis/li_t1_tree.nwk")
max_likelihood_t2_wgs_tree <- ape::read.tree(file = "../li-application/nextstrain_analysis/li_t2_tree.nwk")

t1_treePars_obj <- treeio::as.treedata(max_likelihood_t1_wgs_tree)
t2_treePars_obj <- treeio::as.treedata(max_likelihood_t2_wgs_tree)

##Tumor 1
t1_treePars_obj@data <- tibble("node" = 1:(length(t1_treePars_obj@phylo$tip.label) + length(t1_treePars_obj@phylo$node.label)), 
                               "punch" = c(t1_treePars_obj@phylo$tip.label, rep(NA, length(t1_treePars_obj@phylo$node.label)))) %>% 
    
    dplyr::left_join(., t1_li_edge_labels, by = "punch") %>% 
    mutate(punch_label = toupper(gsub("t1", "", punch)))


t1_treePars_obj@data$punch_label[t1_treePars_obj@data$punch_label == "BLOOD"] <- "Normal"

#Convert divergence to number of mutations for branch lengths by multiplying by sequence length
t1_treePars_obj@phylo$edge.length <- t1_treePars_obj@phylo$edge.length * t1_seq_length


##Tumor 2
t2_treePars_obj@data <- tibble("node" = 1:(length(t2_treePars_obj@phylo$tip.label) + length(t2_treePars_obj@phylo$node.label)), 
                               "punch" = c(t2_treePars_obj@phylo$tip.label, rep(NA, length(t2_treePars_obj@phylo$node.label)))) %>% 
    dplyr::left_join(., t2_li_edge_labels, by = "punch") %>% 
    mutate(punch_label = toupper(gsub("t2", "", punch)))


t2_treePars_obj@data$punch_label[t2_treePars_obj@data$punch_label == "BLOOD"] <- "Normal"

#Convert branch lengths in to mutations instead of divergence
t2_treePars_obj@phylo$edge.length <- t2_treePars_obj@phylo$edge.length * t2_seq_length

## TREE PLOTTING 

#### LI ET AT TUMOR 1
scale_bar_width <- 10000
scale_bar_start <- 70000
scale_bar_vert_pos <- 14

##trees where blood is a leaf
t1_wgs_tree <- ggtree(t1_treePars_obj, aes(color = ifelse(edgeP==1, "edge","center"), label = punch_label)) +
    scale_color_manual(values = colors) +
    geom_tippoint(size =2) +geom_text(color = "black", nudge_x=4000, size =4, hjust = "left") +
    theme(legend.position = "none") +
    
    #Scale bar
    geom_segment(aes(x=scale_bar_start, xend = scale_bar_start + scale_bar_width,
                     y = scale_bar_vert_pos, yend = scale_bar_vert_pos), color = "black") +
    geom_text(aes(x=scale_bar_start + scale_bar_width/2, y=scale_bar_vert_pos + 1,
                  label = paste(scale_bar_width, "mutations", sep = " ")),
              vjust = "bottom", hjust = "midde", color = "black")
t1_wgs_tree 

ggsave(filename = "../figures/t1_wgs_genetic_state_tree.png", t1_wgs_tree, height = 4, width = 3)

#t1_x_range <- layer_scales(t1_wgs_tree)$x$range$range
# ggsave(filename = "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/figures/t1_wgs_genetic_state_tree.png",
#        t1_wgs_tree, height = 4, width = 3)
#for scale bar
scale_bar_width <- 10000
scale_bar_start <- 70000
scale_bar_vert_pos <- 8

t2_wgs_tree <- ggtree(t2_treePars_obj, aes(color = ifelse(edgeP==1, "edge","center"), label = punch_label)) +
    scale_color_manual(values = colors) +
    geom_tippoint(size =2) +geom_text(color = "black", nudge_x=4000, size =4, hjust = "left") +
    theme(legend.position = "none") +
    geom_segment(aes(x=scale_bar_start, xend = scale_bar_start + scale_bar_width,
                     y = scale_bar_vert_pos, yend = scale_bar_vert_pos), color = "black") +
    geom_text(aes(x=scale_bar_start + scale_bar_width/2, y=scale_bar_vert_pos + 1,
                  label = paste(scale_bar_width, "mutations", sep = " ")),
              vjust = "bottom", hjust = "midde", color = "black")
t2_wgs_tree 
ggsave(filename = "../figures/t2_wgs_genetic_state_tree.png", t2_wgs_tree , height = 4, width = 3)
#Get sscale
#layer_scales(t2_wgs_tree)$x$range$range
#layer_scales(t2_wgs_tree)$y$range$range

##Tumor 3

t3_ling_edge_center_labels <- read_csv(file="../ling-application/data/edge-center-calls.csv") %>% 
    dplyr::select(Punch, edge) %>% 
    rename(punch=Punch)

max_likelihood_t3_wes_tree <- ape::read.tree(file = "../ling-application/nextstrain_analysis/ling_tree.nwk")

t3_treePars_obj <- treeio::as.treedata(max_likelihood_t3_wes_tree)

t3_treePars_obj@data <- tibble("node" = c(1:(length(t3_treePars_obj@phylo$tip.label) + length(t3_treePars_obj@phylo$node.label))), 
                               "punch" = c(t3_treePars_obj@phylo$tip.label, rep(NA, length(t3_treePars_obj@phylo$node.label)))) %>% 
    
    dplyr::left_join(., t3_ling_edge_center_labels, by = "punch")

#Convert divergence to number of mutations for branch lengths by multiplying by sequence length
t3_treePars_obj@phylo$edge.length <- t3_treePars_obj@phylo$edge.length * t3_seq_length

#rescale root -- will visualize this as broken axis
t3_treePars_obj@phylo$edge.length[2] <- 10

t3_treePars_obj@data$is_tumor_root <- "no"
t3_treePars_obj@data$is_tumor_root[t3_treePars_obj@data$node == 25] <- "yes"

#for scale bar
scale_bar_width <- 5
scale_bar_start <- 1
scale_bar_vert_pos <- 20

t3_wes_tree <- ggtree(t3_treePars_obj, aes(color = ifelse(edge==1, "edge","center"),
                                           label = punch,
                                           linetype = is_tumor_root)) +
    scale_color_manual(values = colors) +
    geom_tippoint(size =2) +geom_text(color = "black", nudge_x=0.5, size =8, hjust = "left") +
    theme(legend.position = "none") +
    
    #Scale bar
    geom_segment(aes(x=scale_bar_start, xend = scale_bar_start + scale_bar_width,
                     y = scale_bar_vert_pos, yend = scale_bar_vert_pos), color = "black") +
    geom_text(aes(x=scale_bar_start + scale_bar_width/2, y=scale_bar_vert_pos + 1,
                  label = paste(scale_bar_width, "mutations", sep = " ")),
              vjust = "bottom", hjust = "midde", color = "black") + xlim(c(0, 17))

t3_wes_tree 
ggsave(filename = "../figures/t3_wgs_genetic_state_tree.png", t3_wes_tree, height = 4, width = 3)
#useful for placing scale bars
layer_scales(t3_wes_tree )$x$range$range
#layer_scales(t3_wes_tree )$y$range$range


#t3_x_range <- layer_scales(t3_wgs_tree)$x$range$range

############# RESCALED TREES ##################

#Fit to t1 and t2, root for t3 is too long and should be condensed to visualized
# max_x_range <- max(c(t1_x_range, t2_x_range))
# 
# 
# #rescale root
# t3_rescaled <- t3_treePars_obj
# t3_rescaled@phylo$edge.length[2] <- 0.5
# 
# 
# #find normal-tumor branch (will make this dashed)
# normal_node <- t3_treePars_obj@data$node[which(t3_treePars_obj@data$punch == "Normal")]
# 
# normal_parent <- t3_treePars_obj@phylo$edge[which(t3_treePars_obj@phylo$edge[,2] == normal_node),1]
# 
# root_node <- t3_treePars_obj@phylo$edge[which(t3_treePars_obj@phylo$edge[,1] == normal_parent),2][2]
# 
# 
# #Label normal parent
# t3_rescaled@data$is_tumor_root <- "no"
# t3_rescaled@data$is_tumor_root[t3_rescaled@data$node == root_node] <- "yes"
# 
# t3_rescaled_tree <- t3_rescaled %>% 
#     ggtree(., aes(color = ifelse(edge==1, "edge","center"),
#                   label = punch, 
#                   linetype=is_tumor_root)) +
#     scale_color_manual(values = colors) +
#     geom_tippoint(size =2) +geom_text(color = "black", nudge_x=0.03, size =3) +
#     theme(legend.position = "none") + xlim(c(0, max_x_range)) +
#     scale_linetype_manual(values = c("yes" = "dashed", "no" = "solid"))
# t3_rescaled_tree 
# ggsave(filename = "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/figures/t3_wes_genetic_state_tree.png",
#        t3_rescaled_tree, height = 4, width = 3)
# 
# 
# #Now plot t1 and t2 on the same axes
# 
# t1_wgs_tree <- ggtree(t1_treePars_obj, aes(color = ifelse(edgeP==1, "edge","center"), label = punch_label)) +
#     scale_color_manual(values = colors) +
#     geom_tippoint(size =2) +geom_text(color = "black", nudge_x=0.03, size =3) +
#     theme(legend.position = "none") + xlim(c(0, max_x_range)) 
# t1_wgs_tree 
# 
# ggsave(filename = "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/figures/t1_wgs_genetic_state_tree.png",
#        t1_wgs_tree , height = 4, width = 3)
# 
# t2_wgs_tree <- ggtree(t2_treePars_obj, aes(color = ifelse(edgeP==1, "edge","center"), label = punch_label)) +
#     scale_color_manual(values = colors) +
#     geom_tippoint(size =2) +geom_text(color = "black", nudge_x=0.03, size =3) +
#     theme(legend.position = "none") + xlim(c(0, max_x_range)) 
# t2_wgs_tree 
# 
# ggsave(filename = "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/figures/t2_wgs_genetic_state_tree.png",
#        t2_wgs_tree , height = 4, width = 3)
# 
