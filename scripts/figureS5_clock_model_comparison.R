#Script to generate plots for Figure 5

###### SETUP ######
library(tidyverse)
library(ggtree)
library(treeio)
library(beastio)
library(tumortree)
library(HDInterval)
library(ggnewscale)
library(viridis)
library(phangorn)
library(stringr)

setwd("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/")
colors_loc <- get_color_palette(names=c("edge", "center"))
names(colors_loc) <- c("loc1", "loc0")
###### LOGS -- StateClocks #####
##Tumor 1
### Read in log file
T1_log <- beastio::readLog("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/logs/T1_state_clocks_amplicon_only_weighted_red.log")    
###Collect birthRate posteriors
T1_log_df <- as.data.frame(T1_log) %>% 
    dplyr::mutate("birthRate.loc1" = birthRateSVCanonical.loc1,
                  "birthRate.loc0" = birthRateSVCanonical.loc0) %>% 
    pivot_longer(cols = c(birthRate.loc1, birthRate.loc0), 
                 names_prefix = "birthRate.", names_to = "state", values_to = "birthRate") %>% 
    add_column("clock" = "state_dependent")

##Tumor 2
### Read in log file
T2_log <- beastio::readLog("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/logs/T2_state_clocks_amplicon_only_weighted_red.log")    
###Collect birthRate posteriors
T2_log_df <- as.data.frame(T2_log) %>% 
    dplyr::mutate("birthRate.loc1" = birthRateSVCanonical.loc1,
                  "birthRate.loc0" = birthRateSVCanonical.loc0) %>% 
    pivot_longer(cols = c(birthRate.loc1, birthRate.loc0), 
                 names_prefix = "birthRate.", names_to = "state", values_to = "birthRate") %>% 
    add_column("clock" = "state_dependent")

###### LOGS -- Strict clock #####
##Tumor 1
### Read in log file
T1_strict_clock_log <- beastio::readLog("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/logs/T1_state_clocks_amplicon_only_weighted_red_strict_clock.log")    
###Collect birthRate posteriors
T1_strict_clock_log_df <- as.data.frame(T1_strict_clock_log) %>% 
    dplyr::mutate("birthRate.loc1" = birthRateSVCanonical.loc1,
                  "birthRate.loc0" = birthRateSVCanonical.loc0) %>% 
    pivot_longer(cols = c(birthRate.loc1, birthRate.loc0), 
                 names_prefix = "birthRate.", names_to = "state", values_to = "birthRate") %>% 
    add_column("clock" = "strict")

##
##Tumor 2
### Read in log file
T2_strict_clock_log <- beastio::readLog("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/logs/T2_state_clocks_amplicon_only_weighted_red_strict_clock.log")    
###Collect birthRate posteriors
T2_strict_clock_log_df <- as.data.frame(T2_strict_clock_log) %>% 
    dplyr::mutate("birthRate.loc1" = birthRateSVCanonical.loc1,
                  "birthRate.loc0" = birthRateSVCanonical.loc0) %>% 
    pivot_longer(cols = c(birthRate.loc1, birthRate.loc0), 
                 names_prefix = "birthRate.", names_to = "state", values_to = "birthRate") %>% 
    add_column("clock" = "strict")

###### COMBINED LOG SUMMARIES FOR COMPARISON #####
## Combine state-dependent and strict clock models in same posteriors table
### Tumor 1
T1_df <- bind_rows(list(T1_log_df,T1_strict_clock_log_df))
### Summarize state-depedent birth rate differences
T1_summary_df <- T1_df %>% 
    tidyr::pivot_wider(names_from ="state", names_prefix = "birthRate_", values_from = "birthRate") %>% 
    dplyr::mutate("birthRateDiff" = birthRate_loc1 - birthRate_loc0) %>% 
    group_by(clock) %>% 
    dplyr::summarise("meanBirthRateDiff" = mean(birthRateDiff),
                  "birthRateDiff_hdi95_lower" = hdi(birthRateDiff,  credMass = 0.95)[1],
                  "birthRateDiff_hdi95_upper" = hdi(birthRateDiff,  credMass = 0.95)[2])
### Tumor 2
T2_df <- bind_rows(list(T2_log_df,T2_strict_clock_log_df))

### Summarize state-depedent birth rate differences for each model
T2_summary_df <- T2_df %>% 
    tidyr::pivot_wider(names_from ="state", names_prefix = "birthRate_", values_from = "birthRate") %>% 
    dplyr::mutate("birthRateDiff" = birthRate_loc1 - birthRate_loc0) %>% 
    group_by(clock) %>% 
    dplyr::summarise("meanBirthRateDiff" = mean(birthRateDiff),
                     "birthRateDiff_hdi95_lower" = hdi(birthRateDiff,  credMass = 0.95)[1],
                     "birthRateDiff_hdi95_upper" = hdi(birthRateDiff,  credMass = 0.95)[2])
###### PLOTTING #####

#T1 posteriors
t1_posteriors_plot <- ggplot(T1_df, aes(x = birthRate, fill = state)) + geom_density(alpha = 0.8) + theme_classic() +
    scale_fill_manual(values = colors_loc) + theme(legend.position = "none") +
    facet_wrap(~clock, nrow = 2) +
    theme(text = element_text(size = 20)) + xlab("Estimated birth rate")
t1_posteriors_plot
#T2 posteriors
t2_posteriors_plot <- ggplot(T2_df, aes(x = birthRate, fill = state)) + geom_density(alpha = 0.8) + theme_classic() +
    scale_fill_manual(values = colors_loc) + theme(legend.position = "none") +
    facet_wrap(~clock, nrow = 2) +
    theme(text = element_text(size = 20)) + xlab("Estimated birth rate")
t2_posteriors_plot 
ggsave(plot=t1_posteriors_plot,file ="manuscript/figures/t1_posteriors_model_comparison.png", height = 5, width = 4)
ggsave(plot=t2_posteriors_plot,file ="manuscript/figures/t2_posteriors_model_comparison.png", height = 5, width = 4)

#T1 birth rate differences for each model
t1_clock_model_summary <- ggplot(T1_summary_df, aes(x = clock, y = meanBirthRateDiff)) + geom_point(size = 4) + geom_errorbar(aes(ymin = birthRateDiff_hdi95_lower, ymax = birthRateDiff_hdi95_upper), width = 0, size = 2) +
    theme_classic() + coord_flip() +ylim(0,35) +ylab("") +xlab("") +
    scale_x_discrete(limits = rev) + theme(text = element_text(size = 25))
t1_clock_model_summary 
#T2 birth rate differences for each model
t2_clock_model_summary <- ggplot(T2_summary_df, aes(x = clock, y = meanBirthRateDiff)) + geom_point(size = 4) +
    geom_errorbar(aes(ymin = birthRateDiff_hdi95_lower, ymax = birthRateDiff_hdi95_upper), width = 0, size = 2) +
    theme_classic() + coord_flip() +ylim(0,35) +ylab("") +
    xlab("") +scale_x_discrete(limits = rev) + theme(text = element_text(size = 25))
t2_clock_model_summary 
ggsave(plot=t1_clock_model_summary,file ="manuscript/figures/t1_summary_model_comparison.png", height = 3, width = 7)
ggsave(plot=t2_clock_model_summary,file ="manuscript/figures/t2_summary_model_comparison.png", height = 3, width = 7)
#90% HPD interval T1
#coda::HPDinterval(as.mcmc(T1_log[,"birthRateSVCanonical.loc1"] /T1_log[,"birthRateSVCanonical.loc0"]), prob = .90)
#90% HPD interval
#coda::HPDinterval(as.mcmc(T2_log[,"birthRateSVCanonical.loc1"] /T2_log[,"birthRateSVCanonical.loc0"]), prob = .90)


### MCC TREE COMPARISON ####

#read in MCC trees 
## MCC trees for state-dependent clock
T1_mcc_tree <- read.beast("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/trees/T1_state_clocks_amplicon_only_weighted_red_mcc.tree")
T2_mcc_tree <- read.beast("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/trees/T2_state_clocks_amplicon_only_weighted_red_mcc.tree")

#MCC trees for strict clock
T1_strict_clock_mcc_tree <- read.beast("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/trees/T1_state_clocks_amplicon_only_weighted_red_strict_clock_mcc.tree")
T2_strict_clock_mcc_tree <- read.beast("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/trees/T2_state_clocks_amplicon_only_weighted_red_strict_clock_mcc.tree")

#MCC trees for strict clock fixed toplogy
T1_strict_clock_fixed_topology_mcc_tree <- read.beast("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/trees/T1_state_clocks_amplicon_only_weighted_red_strict_clock_fixed_topology_mcc.tree")
T2_strict_clock_fixed_topology_mcc_tree <- read.beast("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/trees/T2_state_clocks_amplicon_only_weighted_red_strict_clock_fixed_topology_mcc.tree")

ggtree(T1_strict_clock_fixed_topology_mcc_tree, aes(color = as.numeric(posterior)))
## tree plotting

#set sizes
tip_size <- 0.5
line_size <- 0.75

#set colors
colors_edge_center<- tumortree::get_color_palette(names = c("edge", "center"))
colors_loc <- tumortree::get_color_palette(names = c("edge", "center"))

names(colors_loc) <- c("loc1", "loc0")
support_cutoff <- 0.5
t1_nodes_to_collapse <- sort(as.numeric(T1_mcc_tree@data$node[which(as.numeric(T1_mcc_tree@data$posterior) < support_cutoff)]))
t2_nodes_to_collapse <- sort(as.numeric(T1_mcc_tree@data$node[as.numeric(T2_mcc_tree@data$posterior) < support_cutoff]))
#state-dependent clock
t1_treeplot <- ggtree(T1_mcc_tree, color = "darkgrey", size = line_size) +
    geom_tippoint(size = tip_size, aes(color = type), alpha = 0.8) +
    scale_color_manual(values = colors_loc) +
    theme(legend.position = "none") #+
#    geom_point2(aes(subset = as.numeric(posterior) < support_cutoff), shape = 8)
    
t1_treeplot

# test <- ggtree::collapse(t1_treeplot, node = c(184, 185))
# test + geom_point2(aes(subset=(node == 181)), size=5, shape=23, fill="steelblue")
t2_treeplot <- ggtree(T2_mcc_tree, color = "darkgrey", size = line_size) +
    geom_tippoint(size = tip_size, aes(color = type), alpha = 0.8) +
    scale_color_manual(values = colors_loc) +
    theme(legend.position = "none")

t2_treeplot

#strict clocks

t1_strict_clock_treeplot <- ggtree(T1_strict_clock_fixed_topology_mcc_tree, color = "darkgrey", size = line_size) +
    geom_tippoint(size = tip_size, aes(color = type), alpha = 0.8) +
    scale_color_manual(values = colors_loc) +
    theme(legend.position = "none") +
    scale_x_reverse()
t1_strict_clock_treeplot

t2_strict_clock_treeplot <- ggtree(T2_strict_clock_fixed_topology_mcc_tree, color = "darkgrey", size = line_size) +
    geom_tippoint(size = tip_size, aes(color = type), alpha = 0.8) +
    scale_color_manual(values = colors_loc) +
    theme(legend.position = "none") +
    scale_x_reverse()

t2_strict_clock_treeplot 


######### TANGLEGRAMS ###############

##followed this blog post on how to make tanglegrams from ggtrees
#arftrhmn.net/how-to-make-cophylogeny
#get data from both tree plots

#Tumor 1
d1 <- t1_treeplot$data
d2 <- t1_strict_clock_treeplot$data


#modify x of strict clock to put on same plot
d2$x <- max(d2$x) - d2$x + max(d1$x) +  max(d1$x)*0.3

#lines are going to be drawn just from tips
connections <- bind_rows(list(d1, d2)) %>% 
    filter(!is.na(label)) #find tips by those that have labels


#combine tree plot
pp1 <- t1_treeplot + geom_tree(data=d2, color = "darkgrey", size = line_size) + #base mcc tree + strict clock flipped
    geom_tippoint(data=d2, aes(color = type), alpha = 0.8, size = tip_size)  + #add tips to secon tree (strict clock)
    geom_line(data = connections , aes(x, y, group=label, color = type)) #add tangle
pp1
ggsave(plot=pp1,file ="manuscript/figures/t1_clocks_tanglegram.png", height = 5, width = 5)
#Tumor 2
d1 <- t2_treeplot$data
d2 <- t2_strict_clock_treeplot$data

#modify x of strict clock to put on same plot
d2$x <- max(d2$x) - d2$x + max(d1$x) +  max(d1$x)*0.3

#lines are going to be drawn just from tips
connections <- bind_rows(list(d1, d2)) %>% 
    filter(!is.na(label))

#combine tree plot
pp2 <- t2_treeplot + geom_tree(data=d2, color = "darkgrey", size = line_size) + #base mcc tree + strict clock flipped
    geom_tippoint(data=d2, aes(color = type), alpha = 0.8, size = tip_size)  + #add tips to secon tree (strict clock)
    geom_line(data = connections , aes(x, y, group=label, color = type)) #add tangles
pp2
ggsave(plot=pp2,file ="manuscript/figures/t2_clocks_tanglegram.png", height = 5, width = 5)

########## COMPARE ANCESTRAL NODE TIMINGS ###########

#tree <- T1_mcc_tree
#comparison_tree <- T1_strict_clock_mcc_tree
#get corresponding node based on clades (MRCA of subset of tips)

get_corresponding_node <- function(node, tree, comparison_tree) {
    #get daughter tips of node from original tree
    tips <- unlist(phangorn::Descendants(tree@phylo, node = node, type = "tips"))
    
    #get labels of tips
    labels <- tree@phylo$tip.label[tips]
    
    #get leaves of target comparison tree
    
    target_tree_tips <- which(comparison_tree@phylo$tip.label %in% labels)
    
    #get MRCA of target tree tips to get target node
    
    target_tree_node <- ape::getMRCA(comparison_tree@phylo, target_tree_tips )
    
    return(target_tree_node)
    
}
#get all internal nodes to compare for each tumor
t1_nodes <- seq((length(T1_mcc_tree@phylo$tip.label) + 1), length(T1_mcc_tree@phylo$tip.label) +T1_mcc_tree@phylo$Nnode)
t2_nodes <- seq((length(T2_mcc_tree@phylo$tip.label) + 1), length(T2_mcc_tree@phylo$tip.label) +T2_mcc_tree@phylo$Nnode)

#get target nodes for both tumors 

t1_strict_clock_nodes <- map_dbl(t1_nodes, function(n) get_corresponding_node(node = n, tree = T1_mcc_tree, comparison_tree = T1_strict_clock_mcc_tree))
t2_strict_clock_nodes <- map_dbl(t2_nodes, function(n) get_corresponding_node(node = n, tree = T2_mcc_tree, comparison_tree = T2_strict_clock_mcc_tree))

#make comparison dataset

#Tumor 1
t1_node_height_comparison_df <- T1_mcc_tree@data %>% 
    dplyr::filter(node %in% t1_nodes) 

t1_node_height_comparison_df$strict_clock_node <- t1_strict_clock_nodes[match(as.integer(t1_node_height_comparison_df$node), t1_nodes)]
t1_node_height_comparison_df$strict_clock_node_height <- as.numeric(T1_strict_clock_mcc_tree@data$height[match(t1_node_height_comparison_df$strict_clock_node, T1_strict_clock_mcc_tree@data$node)])


#Tumor 2
t2_node_height_comparison_df <- T2_mcc_tree@data %>% 
    dplyr::filter(node %in% t2_nodes) 

t2_node_height_comparison_df$strict_clock_node <- t2_strict_clock_nodes[match(as.integer(t2_node_height_comparison_df$node), t2_nodes)]
t2_node_height_comparison_df$strict_clock_node_height <- as.numeric(T2_strict_clock_mcc_tree@data$height[match(t2_node_height_comparison_df$strict_clock_node, T2_strict_clock_mcc_tree@data$node)])

##### HEIGHT COMPARISON PLOTS ######

joint_posterior_min <- min(t1_node_height_comparison_df$posterior, t2_node_height_comparison_df$posterior)
#close to 0, so just use 0-1 as the range
t1_node_height_comparison_plot <- ggplot(t1_node_height_comparison_df, aes(x=as.numeric(height), y=strict_clock_node_height, fill = type)) +
    geom_point(alpha =0.8, shape = 21, aes(size = as.numeric(posterior))) + theme_classic() +
    scale_fill_manual(values = colors_loc) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_y_continuous(trans='log10') +
    scale_x_continuous(trans='log10')+
    theme(text = element_text(size = 15)) +
    xlab("log(Node height with state-dependent clock)") +
    ylab("log(Node height with strict clock)") +
    scale_size(limits = c(0,1)) +
    theme(legend.position = c(0.25, 0.85)) +
    labs(size = "Node posterior") +
    guides(fill = "none") #+
    # theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='dashed'),
    #       legend.box.margin = margin(t = 2, l = 2))
    
t1_node_height_comparison_plot
ggsave(plot=t1_node_height_comparison_plot,file ="manuscript/figures/t1_clocks_height_comparison_plot.png", height = 7, width = 7)

t2_node_height_comparison_plot <- ggplot(t2_node_height_comparison_df, aes(x=as.numeric(height), y=strict_clock_node_height,
                                                                           fill= type)) +
    geom_point(alpha =0.8, shape = 21, aes(size = as.numeric(posterior))) + theme_classic() +
    scale_fill_manual(values = colors_loc) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
    scale_y_continuous(trans='log10') +
    scale_x_continuous(trans='log10') +
    theme(text = element_text(size = 15)) +
    xlab("log(Node height with state-dependent clock)") +
    ylab("log(Node height with strict clock)") +
    scale_size(limits = c(0,1)) +
    theme(legend.position = c(0.25, 0.85)) +
    labs(size = "Node posterior") +
    guides(fill = "none") #+
#    theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='dashed'),
#          legend.box.margin = margin(t = 2, l = 2))
t2_node_height_comparison_plot    
ggsave(plot=t2_node_height_comparison_plot,file ="manuscript/figures/t2_clocks_height_comparison_plot.png", height = 7, width = 7)


#### CLADE COMPARISON ######

#output of CladeSetComparator Beast2 tool
#/Applications/BEAST\ 2.6.2/bin/applauncher CladeSetComparator -tree1 trees/T2_state_clocks_amplicon_only_weighted_red_strict_clock.HCCtumor.trees -tree2 trees/T2_state_clocks_amplicon_only_weighted_red.HCCtumor.trees -out T2_clock_model_clade_comparison_out.txt
#/Applications/BEAST\ 2.6.2/bin/applauncher CladeSetComparator -tree1 trees/T1_state_clocks_amplicon_only_weighted_red_strict_clock.HCCtumor.trees -tree2 trees/T1_state_clocks_amplicon_only_weighted_red.HCCtumor.trees -out T1_clock_model_clade_comparison_out.txt

#Tumor 1
t1_clock_model_clade_comparison <- read_table("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/T1_clock_model_clade_comparison_out.txt", col_names = FALSE)
colnames(t1_clock_model_clade_comparison) <- c("clade", "state_dependent", "strict")
t1_clock_model_clade_comparison <- t1_clock_model_clade_comparison %>% 
    dplyr::mutate("deviation" = abs(state_dependent - strict))
#Tumor 2
t2_clock_model_clade_comparison <- read_table("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/T2_clock_model_clade_comparison_out.txt", col_names = FALSE)
colnames(t2_clock_model_clade_comparison) <- c("clade", "state_dependent", "strict")
t2_clock_model_clade_comparison <- t2_clock_model_clade_comparison %>% 
    dplyr::mutate("deviation" = abs(state_dependent - strict))

#find ranges for color 
max_deviation <- max(t1_clock_model_clade_comparison$deviation, t2_clock_model_clade_comparison$deviation)
#plotting

#Tumor 1 
t1_clock_model_clade_comparison_scatter <- t1_clock_model_clade_comparison %>% 
    dplyr::filter(state_dependent > 0.05) %>% 
    ggplot(., aes(x=state_dependent, y=strict, color = deviation)) + geom_point(size = 2, alpha =0.8) +
    theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_color_viridis(option = "magma", limits = c(0, max_deviation)) + theme(text= element_text(size = 15)) + xlab("Clade probability (state-dependent clock)") +
    ylab("Clade probability (strict clock)") + theme(legend.position = "none")
t1_clock_model_clade_comparison_scatter
#Tumor 2
t2_clock_model_clade_comparison_scatter <- ggplot(t2_clock_model_clade_comparison, aes(x=state_dependent, y=strict, color = deviation)) +
    geom_point(alpha =0.8, size =2) + theme_classic() + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_color_viridis(option = "magma", limits = c(0, max_deviation)) + theme(text= element_text(size = 15)) + xlab("Clade probability (state-dependent clock)") +
    ylab("Clade probability (strict clock)") + theme(legend.position = "none")
t2_clock_model_clade_comparison_scatter

#function to match clade to given MCC tree
get_clade_tips <- function(node, mcc_tree) {
    node_tips <- unlist(phangorn::Descendants(mcc_tree@phylo, node = node, type = "tips"))
    clade_tips <- mcc_tree@phylo$tip.label[node_tips]
    return(clade_tips)
}

#mcc_tree_clades_list <- t1_mcc_clades_list
match_clade <- function(tree_clade_vec, all_clades_list) {

    #find equal set
    f1 <- function(x,y) setequal(x,y)
    node_match <- which(mapply(f1, all_clades_list, list(tree_clade_vec)))

    return(node_match)
    if (length(node_match) == 0) {
    #     
    #     return(NA)
    #     
    } else {
    #     
    #     return(node_match)
    }

    
}

get_clade_string_vector <- function(clade_string) {
    
    tips <- unlist(str_split(clade_string, pattern = "[{},]" ))
    tips <- tips[tips != ""]
    
    return(tips)
}



#Tumor 1 -- state-dependent clock
##Get names of tips that define clades to match 
t1_mcc_clades_list <- purrr::map(t1_nodes, function(n) get_clade_tips(node=n, mcc_tree = T1_mcc_tree))
names(t1_mcc_clades_list) <- t1_nodes

##Extract names of tips that define clades from string (from CladeSetComparator)
t1_all_clades_list <- purrr::map(t1_clock_model_clade_comparison$clade, function(clade_string) get_clade_string_vector(clade_string))

##Match each node to corresponding clade in clade set
t1_mcc_matched_nodes <- purrr::map_dbl(t1_mcc_clades_list, function(c) match_clade(tree_clade_vec = c,
                                                                                   all_clades_list = t1_all_clades_list))
#Add matched info to tree

T1_mcc_tree@data$deviation <- 0

T1_mcc_tree@data$deviation[match(t1_nodes, T1_mcc_tree@data$node)] <- t1_clock_model_clade_comparison$deviation[t1_mcc_matched_nodes]

t1_mcc_tree_clade_divergence <- ggtree(T1_mcc_tree, aes(color = deviation)) +
    scale_color_viridis(option = "magma", limits = c(0,max_deviation)) + geom_point(aes(color = deviation, size = as.numeric(posterior))) +
    scale_size(range = c(0,2), limits = c(0, 1))
t1_mcc_tree_clade_divergence
#add heatmape with states

#make compatible data frame for heatmap with only type column and rownames as tip labels
T1_state_labels <- as.data.frame(T1_mcc_tree@data) %>% 
    dplyr::filter(as.numeric(node) <= length(T1_mcc_tree@phylo$tip.label))

rownames(T1_state_labels) <- T1_mcc_tree@phylo$tip.label[as.numeric(T1_state_labels$node)]
T1_state_labels <- T1_state_labels %>% 
    dplyr::select(type)


p2 <- t1_mcc_tree_clade_divergence + new_scale_fill()
t1_mcc_tree_clade_divergence_with_states <- gheatmap(p2, T1_state_labels, offset=0.01, width=.01,
         colnames = FALSE) +
    scale_fill_manual(values = colors_loc) + guides(fill = "none")

t1_mcc_tree_clade_divergence_with_states
#Tumor 2
##Get names of tips that define clades to match 
t2_mcc_clades_list <- purrr::map(t2_nodes, function(n) get_clade_tips(node=n, mcc_tree = T2_mcc_tree))
names(t2_mcc_clades_list) <- t2_nodes

##Extract names of tips that define clades from string (from CladeSetComparator)
t2_all_clades_list <- purrr::map(t2_clock_model_clade_comparison$clade, function(clade_string) get_clade_string_vector(clade_string))

##Match each node to corresponding clade in clade set
t2_mcc_matched_nodes <- purrr::map_dbl(t2_mcc_clades_list, function(c) match_clade(tree_clade_vec = c,
                                                                                   all_clades_list = t2_all_clades_list))
#Add matched info to tree

T2_mcc_tree@data$deviation <- 0

T2_mcc_tree@data$deviation[match(t2_nodes, T2_mcc_tree@data$node)] <- t2_clock_model_clade_comparison$deviation[t2_mcc_matched_nodes]

t2_mcc_tree_clade_divergence <- ggtree(T2_mcc_tree, color = "darkgrey") + scale_color_viridis(option = "magma", limits = c(0,max_deviation)) +
    geom_point(aes(color = deviation, size = deviation)) +
    scale_size(range = c(0,2), limits = c(0,max_deviation))
t2_mcc_tree_clade_divergence
#make compatible data frame for heatmap with only type column and rownames as tip labels
T2_state_labels <- as.data.frame(T2_mcc_tree@data) %>% 
    dplyr::filter(as.numeric(node) <= length(T2_mcc_tree@phylo$tip.label))

rownames(T2_state_labels) <- T2_mcc_tree@phylo$tip.label[as.numeric(T2_state_labels$node)]
T2_state_labels <- T2_state_labels %>% 
    dplyr::select(type)


p3 <- t2_mcc_tree_clade_divergence + new_scale_fill()
t2_mcc_tree_clade_divergence_with_states <- gheatmap(p3, T2_state_labels, offset=0.01, width=.01,
         colnames = FALSE) +
    scale_fill_manual(values = colors_loc) + guides(fill = "none")

t2_mcc_tree_clade_divergence_with_states  

#save plots 
ggsave(plot=t1_mcc_tree_clade_divergence_with_states,file ="manuscript/figures/t1_mcc_tree_clade_divergence_with_states.png", height = 7, width = 10)
ggsave(plot=t2_mcc_tree_clade_divergence_with_states,file ="manuscript/figures/t2_mcc_tree_clade_divergence_with_states.png", height = 7, width = 10)

ggsave(plot=t1_clock_model_clade_comparison_scatter,file ="manuscript/figures/t1_clock_model_clade_comparison_scatter.png", height = 5, width = 5)
ggsave(plot=t2_clock_model_clade_comparison_scatter,file ="manuscript/figures/t2_clock_model_clade_comparison_scatter.png", height = 5, width = 5)

############ COMPARE STRICT CLOCK WITH FIXED TOPLOGY ######
#Tumor 1

##Remove branch lengths and write newick tree
##Newick tree is used in fixed topology xml
T1_mcc_tree_topo_only <- T1_mcc_tree@phylo
T1_mcc_tree_topo_only$edge.length <- NA
T1_nwk <- write.tree(T1_mcc_tree_topo_only, digits = 0)
gsub(":NA", "", T1_nwk, fixed = TRUE)

#Tumor 2

##Remove branch lengths and write newick tree
##Newick tree is used in fixed topology xml
T2_mcc_tree_topo_only <- T2_mcc_tree@phylo
T2_mcc_tree_topo_only$edge.length <- NA
T2_nwk <- write.tree(T2_mcc_tree_topo_only, digits = 0)
gsub(":NA", "", T2_nwk, fixed = TRUE)



