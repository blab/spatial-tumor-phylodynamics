#figure_1.R
##script to generate visualizations for figure 1

### SETUP ####
library(tumortree)
library(tidyverse)
library(viridis)
library(ape)
library(ggtree)
library(cowplot)

figures_dir <- "../figures"

#color palette
sim_colors <- get_color_palette(c("boundary_driven", "unrestricted"))
#sim_colors <- c("boundary_driven" = "#e49a8b", "unrestricted" = "#3C3C3C")

##### REPRESENATIVE SIMULATIONS -- TUMOR GROWTH PATTERNS AND TREES (Figure 1 A-F) ######

###### TUMORS #####
# Read in representative boundary-driven tumor
## Can skip this and read in example cells

## All simulated cells

#record of local directory location

# all_cells_boundary_driven <- read_csv("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation/cells_death_rate_validation_pop_1000_dr_0.10.csv") %>%
#     normalize_locs

all_cells_boundary_driven <- read_csv("../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.050.csv") %>%
    normalize_locs

## Filter to alive cells at end of simulation (resected tumor)
alive_cells_boundary_driven <- all_cells_boundary_driven %>%
    filter_alive_cells


## read in example tumor for unrestricted cell growth

# Record of local directory location
# all_cells_unrestricted <- read_csv("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation/cells_pushing_pop_1000_dr_0.10.csv") %>%
#      normalize_locs

all_cells_unrestricted <- read_csv("../eden/simulation_data/cells_pushing_pop_10000_mu_1_dr_0.050.csv") %>%
    normalize_locs

alive_cells_unrestricted <- all_cells_unrestricted %>%
    filter_alive_cells

########## Add mean growth info #########

# Boundary-driven tumor
#add number of divisions and mean growth rates to alive cells
alive_cells_boundary_driven$n_divisions <- purrr::map_dbl(alive_cells_boundary_driven$index, function(index) get_n_divisions(index, all_cells_boundary_driven))

#all cells are collected at the same time point, so mean growth rate is just the number of divisions over simulation time
alive_cells_boundary_driven$mean_growth_rate <- alive_cells_boundary_driven$n_divisions  /  max(alive_cells_boundary_driven$deathdate)
alive_cells_boundary_driven$relative_mean_growth_rate <- alive_cells_boundary_driven$mean_growth_rate / mean(alive_cells_boundary_driven$mean_growth_rate)

# Unrestricted growth tumor

#calculate number divisions and growth rates
alive_cells_unrestricted$n_divisions <- purrr::map_dbl(alive_cells_unrestricted$index, function(index) get_n_divisions(index,all_cells_unrestricted))
alive_cells_unrestricted$mean_growth_rate <- alive_cells_unrestricted$n_divisions  / max(alive_cells_unrestricted$deathdate)

alive_cells_unrestricted$relative_mean_growth_rate <- alive_cells_unrestricted$mean_growth_rate / mean(alive_cells_unrestricted$mean_growth_rate)

# Get sampled cells from saved trees to stay consistent
# tree_boundary_driven <- readRDS(file = "../eden/simtrees/boundary_driven_timetree.rds")
# tree_molecular_boundary_driven <- readRDS(file = "../eden/simtrees/boundary_driven_moleculartree.rds")
# 

# tree_unrestricted <- readRDS(file = "../eden/simtrees/unrestricted_timetree.rds")
# tree_molecular_unrestricted <- readRDS(file = "../eden/simtrees/unrestricted_moleculartree.rds")


## Boundary-driven tumor

# alive_cells_boundary_driven$sampled <- alive_cells_boundary_driven$index %in% tree_boundary_driven@data$index
# sampled_cells_boundary_driven <- alive_cells_boundary_driven %>%
#     filter(sampled == TRUE)
# 
# ## Unrestricted tumor
# 
# alive_cells_unrestricted$sampled <- alive_cells_unrestricted$index %in% tree_unrestricted@data$index
# sampled_cells_unrestricted <- alive_cells_unrestricted %>%
#     filter(sampled == TRUE)


## Can skip to here if don't have simulation data generated
# alive_cells_boundary_driven <- read_csv("../eden/simulation_data/alive_cells_death_rate_validation_pop_1000_dr_0.10.csv")
# alive_cells_unrestricted <- read_csv("../eden/simulation_data/alive_cells_pushing_pop_1000_dr_0.10.csv")

###### TUMOR TREES ########

#To make trees from simulated cells
## Can skip ahead and load pre-constructed trees

## Sample from alive cells

set.seed(3812)
alive_cells_boundary_driven <- tumortree::sample_alive_cells(alive_cells_boundary_driven,
                                                               n=100,
                                                               diversified_sampling = TRUE)

sampled_cells_boundary_driven <- alive_cells_boundary_driven %>%
    dplyr::filter(sampled)

set.seed(3212)

alive_cells_unrestricted <- tumortree::sample_alive_cells(alive_cells_unrestricted, n=100, diversified_sampling = TRUE) 
sampled_cells_unrestricted <-  alive_cells_unrestricted %>%
    dplyr::filter(sampled == TRUE)

write_csv(alive_cells_unrestricted, "../eden/simulation_data/alive_cells_pushing_pop_10000_mu_1_dr_0.050.csv")
write_csv(alive_cells_boundary_driven, "../eden/simulation_data/alive_cells_death_rate_validation_pop_10000_mu_1_dr_0.050.csv")

# #
# # Get true simulated trees
# ## Unrestricted growth
#
# ## Time tree
# tree_boundary_driven <- convert_all_cells_to_tree_fast(all_cells = all_cells_boundary_driven,
#                                                     add_all_states = TRUE,
#                                                     branch_unit = "time")

#While reconstructing the full time tree and pruning is fairly fast, calculating the genetic branch lengths takes longer
## Therefore, it's faster to just re-calculate the branch lengths for the tree once it's pruned. 
library(stringr)
#pruned_tree <- tree_boundary_driven_pruned
time_to_genetic_bl <- function(pruned_tree) {
    
    #sequences for parent nodes
    
    parent_seqs <- purrr::map_chr(pruned_tree@data$sequence[match(pruned_tree@phylo$edge[,1], pruned_tree@data$node)],
                                  function(s) stringr::str_replace_all(s,"[^[:alnum:]]", ""))
    #sequences for daughter nodes
    daughter_seqs <- purrr::map_chr(pruned_tree@data$sequence[match(pruned_tree@phylo$edge[,2], pruned_tree@data$node)],
                                  function(s) stringr::str_replace_all(s,"[^[:alnum:]]", ""))
    #find number of character differences
    #bl <- purrr::map2_dbl(parent_seqs, daughter_seqs, function(a,b) adist(a,b))
    bl <- mapply(function(x,y) sum(x!=y),strsplit(parent_seqs,""),strsplit(daughter_seqs,""))
    genetic_tree <- pruned_tree
    
    genetic_tree@phylo$edge.length <- bl
    
    return(genetic_tree)
}
tree_boundary_driven <- readRDS("../eden/simtrees/cells_death_rate_validation_pop_10000_mu_1_dr_0.050.rds")

tree_molecular_boundary_driven <- time_to_genetic_bl(tree_boundary_driven)

saveRDS(tree_molecular_boundary_driven, file = "../eden/simtrees/boundary_driven_moleculartree_large_full.rds")
#Repeat prune to add root
tree_boundary_driven_pruned <- prune_simulated_tree(tree = tree_boundary_driven,
                                          sampled_cells_indices = sampled_cells_boundary_driven$index)

# ## Genetic tree

tree_molecular_boundary_driven_pruned <- time_to_genetic_bl(tree_boundary_driven_pruned)
# tree_molecular_boundary_driven <- convert_all_cells_to_tree_fast(all_cells = all_cells_boundary_driven,
#                                                               add_all_states = FALSE,
#                                                               branch_unit = "genetic")

 
# tree_molecular_boundary_driven_pruned <- prune_simulated_tree(tree = tree_molecular_boundary_driven,
#                                                     sampled_cells_indices = sampled_cells_boundary_driven$index)

# Unrestricted growth

## Time tree
# tree_unrestricted <- convert_all_cells_to_tree_fast(all_cells = all_cells_unrestricted,
#                                add_all_states = TRUE,
#                                branch_unit = "time")

tree_unrestricted <- readRDS("../eden/simtrees/cells_pushing_pop_10000_mu_1_dr_0.050.rds")
tree_unrestricted_pruned <- prune_simulated_tree(tree = tree_unrestricted,
                                          sampled_cells_indices = sampled_cells_unrestricted$index)



# ## Genetic tree
# tree_molecular_unrestricted <- convert_all_cells_to_tree_fast(all_cells = all_cells_unrestricted,
#                                                     add_all_states = TRUE,
#                                                     branch_unit = "genetic")
tree_molecular_unrestricted <- time_to_genetic_bl(tree_unrestricted)

saveRDS(tree_molecular_unrestricted, file = "../eden/simtrees/unrestricted_moleculartree_large_full.rds")

tree_molecular_unrestricted_pruned <- time_to_genetic_bl(tree_unrestricted_pruned)
# tree_molecular_unrestricted_pruned <- prune_simulated_tree(tree = tree_molecular_unrestricted,
#                                           sampled_cells_indices = sampled_cells_unrestricted$index)


### To save as RDS
saveRDS(tree_boundary_driven_pruned, file = "../eden/simtrees/boundary_driven_timetree_large.rds")
saveRDS(tree_molecular_boundary_driven_pruned, file = "../eden/simtrees/boundary_driven_moleculartree_large.rds")

saveRDS(tree_unrestricted_pruned, file = "../eden/simtrees/unrestricted_timetree_large.rds")
saveRDS(tree_molecular_unrestricted_pruned, file = "../eden/simtrees/unrestricted_moleculartree_large.rds")
# 
# tree_boundary_driven <- readRDS(file = "../eden/simtrees/boundary_driven_timetree.rds")
# tree_molecular_boundary_driven <- readRDS(file = "../eden/simtrees/boundary_driven_moleculartree.rds")
# 
# tree_unrestricted <- readRDS(file = "../eden/simtrees/unrestricted_timetree.rds")
# tree_molecular_unrestricted <- readRDS(file = "../eden/simtrees/unrestricted_moleculartree.rds")

####### Extract number of mutations for each cell ######
sampled_cells_boundary_driven$n_muts <- unlist(purrr::map(sampled_cells_boundary_driven$mutations, function(muts) str_count(muts, ",") + 1))
sampled_cells_unrestricted$n_muts <- unlist(purrr::map(sampled_cells_unrestricted$mutations, function(muts) str_count(muts, ",") + 1))

##### Normalize axes ########
#timetree_limits <- c(0, max(c(alive_cells_boundary_driven$deathdate, alive_cells_unrestricted$deathdate)))

#genetictree_limits <- c(0, max(sampled_cells_boundary_driven$n_muts, sampled_cells_unrestricted$n_muts))

#get growth rate limits for both tumors to standardize color map
growth_rate_limits = c(min(c(alive_cells_boundary_driven$mean_growth_rate, alive_cells_unrestricted$mean_growth_rate)),
                       max(c(alive_cells_boundary_driven$mean_growth_rate, alive_cells_unrestricted$mean_growth_rate)))

#Get max diverged sample for annotation layer
most_diverged_cell_boundary_driven <- sampled_cells_boundary_driven[which.max(sampled_cells_boundary_driven$n_muts),]
least_diverged_cell_boundary_driven <- sampled_cells_boundary_driven[which.min(sampled_cells_boundary_driven$n_muts),]
most_diverged_cell_unrestricted <- sampled_cells_unrestricted[which.max(sampled_cells_unrestricted$n_muts),]
least_diverged_cell_unrestricted <- sampled_cells_unrestricted[which.min(sampled_cells_unrestricted$n_muts),]

### SUBFIGURE tumor_sim_growth_rate_boundary_driven.png###

# To make mini x-y axes
min_x <- min(alive_cells_boundary_driven$norm_locx)
min_y <- min(alive_cells_boundary_driven$norm_locy)

#Figure 1A
a <- ggplot(alive_cells_boundary_driven, aes(x = norm_locx, y = norm_locy, color = mean_growth_rate)) +
    geom_point(size = 1, alpha = 1) +
    theme_void() +
    scale_color_viridis(limits = growth_rate_limits) +# +
    scale_fill_viridis(limits = growth_rate_limits) +
    #scale_color_manual(values = c("grey", sim_colors["boundary_driven"][[1]])) +
    labs(color = "Mean\nbirth\nrate") +
    guides(fill="none") +
    annotate("text",
            x = most_diverged_cell_boundary_driven$norm_locx, y = most_diverged_cell_boundary_driven$norm_locy, size = 6, color = "black", label = "X") +
    annotate( "text",
            x = least_diverged_cell_boundary_driven$norm_locx, y = least_diverged_cell_boundary_driven$norm_locy, size = 6, color = "black",label = "X" ) +
    annotate(x=min_x, xend=min_x + 8, y=min_y, yend=min_y, colour="black", lwd=0.5, geom="segment", arrow = arrow(ends = "last", length = unit(.2,"cm"), type = "closed")) +
    annotate(x=min_x, xend=min_x, y=min_y, yend=min_y + 8, colour="black", lwd=0.5, geom="segment", arrow = arrow(ends = "last", length = unit(.2,"cm"), type = "closed")) +
    annotate(x=min_x + 10, y=min_y, colour="black", geom="text", label = "X") +
    annotate(x=min_x, y=min_y + 11, colour="black", geom="text", label = "Y") +
    theme(text=element_text(size =15))
a

#Inset for Figure 2A
a_void <- ggplot(alive_cells_boundary_driven, aes(x = norm_locx, y = norm_locy, color = mean_growth_rate)) +
    geom_point(size = 1, alpha = 1) +
    theme_void() +
    scale_color_viridis(limits = growth_rate_limits) +# +
    scale_fill_viridis(limits = growth_rate_limits) +
    theme(legend.position = "none")
    
a_void 

##### SAVE FIGURES ######
ggsave(file = "tumor_sim_growth_rate_boundary_driven_void_large.png",
       plot = a_void , path = figures_dir, height = 4, width = 4)

ggsave(file = "tumor_sim_growth_rate_boundary_driven_large.png",
       plot = a, path = figures_dir, height = 4, width = 5)



### SUBFIGURE Figure 1B ###

#add mean growth rate to tree
tree_boundary_driven_pruned@data$mean_growth_rate <- sampled_cells_boundary_driven$mean_growth_rate[match(tree_boundary_driven_pruned@data$index,
                                                                                                          sampled_cells_boundary_driven$index)]

#mark cells of interest (least and most diverged)
tree_boundary_driven_pruned@data$mark <- purrr::map_chr(tree_boundary_driven_pruned@data$index,
                                                        function(i) ifelse(i == as.character(least_diverged_cell_boundary_driven$index), "*",
                                                                           ifelse(i == as.character(most_diverged_cell_boundary_driven$index), "*", NA)))

#Figure 1B
b <- ggtree(tree_boundary_driven_pruned, aes(color = mean_growth_rate)) + #geom_tippoint(color = sim_colors["boundary_driven"][[1]], alpha = 0.8, size = 2) +
    geom_tippoint(aes(color = mean_growth_rate), size = 2, alpha = 0.8) +
    theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_color_viridis(limits = growth_rate_limits) +

    geom_text(aes(label=mark), size = 6, nudge_x = 0, nudge_y = 0, color = "black") +
    geom_rootedge() +


    theme(axis.line.x = element_line(arrow = grid::arrow(length = unit(0.2, "cm"),
                                                         ends = "last",
                                                         type = "closed"))) +
    xlab("Time to target tumor size") +
    scale_y_reverse() + geom_rootedge(color = "#939393")
#new_lim <- ggplot_build(b)$layout$panel_scales_y[[1]]$range$range + c(-1, 0)
#b <- b + coord_cartesian(clip = "off", ylim = new_lim, xlim = timetree_limits)

b

b_unlabeled <- ggtree(tree_boundary_driven_pruned, aes(color = mean_growth_rate)) + #geom_tippoint(color = sim_colors["boundary_driven"][[1]], alpha = 0.8, size = 2) +
    geom_tippoint(aes(color = mean_growth_rate), size = 2, alpha = 0.8) +
    theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_color_viridis(limits = growth_rate_limits) +
    
    #geom_text(aes(label=mark), size = 6, nudge_x = 0, nudge_y = 0, color = "black") +
    geom_rootedge() +
    
    
    theme(axis.line.x = element_line(arrow = grid::arrow(length = unit(0.2, "cm"),
                                                         ends = "last",
                                                         type = "closed"))) +
    xlab("Time to target tumor size") +
    scale_y_reverse() + geom_rootedge(color = "#939393")
#new_lim <- ggplot_build(b)$layout$panel_scales_y[[1]]$range$range + c(-1, 0)
#b <- b + coord_cartesian(clip = "off", ylim = new_lim, xlim = timetree_limits)

b_unlabeled

ggsave(file = "time_tree_growth_rate_boundary_driven_large.png", plot = b, path = figures_dir, height = 6, width = 5)
ggsave(file = "time_tree_growth_rate_boundary_driven_unlabeled_large.png", plot = b_unlabeled, path = figures_dir, height = 6, width = 5)

#viewClade(b, MRCA(b,"cell_2528", "cell_1474"))

### END SUBFIGURE ###


### SUBFIGURE Figure 1C ###

#add mean growth rate to tree
tree_molecular_boundary_driven_pruned@data$mean_growth_rate <- sampled_cells_boundary_driven$mean_growth_rate[match(tree_molecular_boundary_driven_pruned@data$index, sampled_cells_boundary_driven$index)]

#mark cells of interest
tree_molecular_boundary_driven_pruned@data$mark <- purrr::map_chr(tree_molecular_boundary_driven_pruned@data$index, function(i) ifelse(i == as.character(least_diverged_cell_boundary_driven$index), "*", ifelse(i == as.character(most_diverged_cell_boundary_driven$index), "*", NA)))

#labeled version for cells of interest reference(not shown)
c <- ggtree(tree_molecular_boundary_driven_pruned, aes(color = mean_growth_rate)) +# geom_tippoint(color = sim_colors["boundary_driven"][[1]], alpha = 0.8, size = 2) +
    geom_tippoint(aes(color = mean_growth_rate), size = 2, alpha = 0.8) +theme(legend.position = "none") +
    coord_cartesian(clip = "off") + scale_color_viridis(limits = growth_rate_limits) +
    geom_text(aes(label=mark), size = 6, nudge_x = 0.6, nudge_y = -0.6, color = "black") +
    scale_y_reverse() + geom_rootedge(color = "#939393")
    #ggtitle("Genetic tree")
c

#unlabeled version for Figure 1C
c_unlabelled <- ggtree(tree_molecular_boundary_driven_pruned, aes(color = mean_growth_rate)) +# geom_tippoint(color = sim_colors["boundary_driven"][[1]], alpha = 0.8, size = 2) +
    geom_tippoint(aes(color = mean_growth_rate), size = 2, alpha = 0.8) +theme(legend.position = "none") +
    scale_color_viridis(limits = growth_rate_limits) +
    # theme(axis.line.x = element_line(arrow = grid::arrow(length = unit(0.2, "cm"),
    #                                                  ends = "last",
    #                                                  type = "closed"))) +
    #xlab("Genetic distance") =
    scale_y_reverse() +
    geom_rootedge(color = "#939393")

c_unlabelled

ggsave(file = "molecular_tree_growth_rate_boundary_driven_large.png", plot = c, path = figures_dir, height = 6, width = 5)
ggsave(file = "molecular_tree_growth_rate_boundary_driven_unlabeled_large.png", plot = c_unlabelled, path = figures_dir, height = 6, width = 5)


### SUBFIGURE 1E ###

#for mini axes
min_x <- min(alive_cells_unrestricted$norm_locx)
min_y <- min(alive_cells_unrestricted$norm_locy)

#labeled version for reference (not shown)
a <- ggplot(alive_cells_unrestricted, aes(x = norm_locx, y = norm_locy, color = mean_growth_rate)) +
    geom_point(size =1)  +
    theme_void() +
    #geom_point(data = sampled_cells_unrestricted, size = 3, alpha = 1, colour="black",pch=21, aes(fill = mean_growth_rate)) +

    scale_color_viridis(limits = growth_rate_limits) +
    scale_fill_viridis(limits = growth_rate_limits) +
    #scale_color_manual(values = c("grey", sim_colors["unrestricted"][[1]])) +
    #theme(legend.position = "none") +
    labs(color = "Mean\nbirth\nrate") +
    guides(fill="none") +
    guides(fill="none") +
    annotate("text",
             x = most_diverged_cell_unrestricted$norm_locx, y = most_diverged_cell_unrestricted$norm_locy, size = 6, color = "black", label = "X") +
    annotate( "text",
              x = least_diverged_cell_unrestricted$norm_locx, y = least_diverged_cell_unrestricted$norm_locy, size = 6, color = "black", label = "X" ) +
    annotate(x=min_x, xend=min_x + 8, y=min_y, yend=min_y, colour="black", lwd=0.5, geom="segment", arrow = arrow(ends = "last", length = unit(.2,"cm"), type = "closed")) +
    annotate(x=min_x, xend=min_x, y=min_y, yend=min_y + 8, colour="black", lwd=0.5, geom="segment", arrow = arrow(ends = "last", length = unit(.2,"cm"), type = "closed")) +
    annotate(x=min_x + 10, y=min_y, colour="black", geom="text", label = "X") +
    annotate(x=min_x, y=min_y + 11, colour="black", geom="text", label = "Y") +
    theme(text=element_text(size=15))



a

#unlabeled version for publication figure
a_unlabeled <- ggplot(alive_cells_unrestricted, aes(x = norm_locx, y = norm_locy, color = mean_growth_rate)) +
    geom_point(size =1) +
    theme_void() +
    #geom_point(data = sampled_cells_unrestricted, size = 3, alpha = 1, colour="black",pch=21, aes(fill = mean_growth_rate)) +
    
    scale_color_viridis(limits = growth_rate_limits) +
    theme(legend.position = "none")
    #labs(color = "Mean\ngrowth\nrate") +
    #guides(fill="none") 



a_unlabeled

ggsave(file = "tumor_sim_growth_rate_unrestricted_large.png", plot = a, path = figures_dir, height = 4, width = 5)
ggsave(file = "tumor_sim_growth_rate_unrestricted_unlabeled_large.png", plot = a, path = figures_dir, height = 4, width = 4)

### SUBFIGURE 1E ###

#add growth rate info to tree
tree_unrestricted_pruned@data$mean_growth_rate <- sampled_cells_unrestricted$mean_growth_rate[match(tree_unrestricted_pruned@data$index, sampled_cells_unrestricted$index)]

#mark cells of interest (least and most diverged)
tree_unrestricted_pruned@data$mark <- purrr::map_chr(tree_unrestricted_pruned@data$index, function(i) ifelse(i == as.character(least_diverged_cell_unrestricted$index), "*", ifelse(i == as.character(most_diverged_cell_unrestricted$index), "*", NA)))

#version with cells marked (for reference)
b <- ggtree(tree_unrestricted_pruned, aes(color = mean_growth_rate)) + #geom_tippoint(color = sim_colors["unrestricted"][[1]], alpha = 0.8, size = 2) +
    geom_tippoint(aes(color = mean_growth_rate), alpha = 0.8, size = 2) +
    theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_color_viridis(limits = growth_rate_limits)+ labs(x = "Time to target tumor size") +
    geom_text(aes(label=mark), size = 6, nudge_x = 0.4, nudge_y = -0.5, color = "black") +


    theme(axis.line.x = element_line(arrow = grid::arrow(length = unit(0.2, "cm"),
                                                         ends = "last",
                                                         type = "closed")))+
    coord_cartesian(clip = "off") +
    geom_rootedge(color = "#939393")


#new_lim <- ggplot_build(b)$layout$panel_scales_y[[1]]$range$range + c(-1, 0)
#b <- b + coord_cartesian(clip = "off", ylim = new_lim, xlim = timetree_limits)
b

#unlabeled version for Figure 1E
b_unlabeled <- ggtree(tree_unrestricted_pruned, aes(color = mean_growth_rate)) + #geom_tippoint(color = sim_colors["unrestricted"][[1]], alpha = 0.8, size = 2) +
    geom_tippoint(aes(color = mean_growth_rate), alpha = 0.8, size = 2) +
    theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_color_viridis(limits = growth_rate_limits)+
    geom_rootedge(color = "#939393")
b_unlabeled
ggsave(file = "time_tree_growth_rate_unrestricted.png", plot = b, path = figures_dir, height = 6, width = 5)
ggsave(file = "time_tree_growth_rate_unrestricted_unlabeled.png", plot = b_unlabeled, path = figures_dir, height = 6, width = 5)

### END SUBFIGURE ###


### SUBFIGURE 1F ###
# (unrestricted genetic tree)
# add growth rate info to tree

tree_molecular_unrestricted_pruned@data$mean_growth_rate <- sampled_cells_unrestricted$mean_growth_rate[match(tree_molecular_unrestricted_pruned@data$index, sampled_cells_unrestricted$index)]

#add growth rate info to tree
tree_molecular_unrestricted_pruned@data$mark <- purrr::map_chr(tree_molecular_unrestricted_pruned@data$index, function(i) ifelse(i == as.character(least_diverged_cell_unrestricted$index), "*", ifelse(i == as.character(most_diverged_cell_unrestricted$index), "*", NA)))

# labeled tree for reference (not shown)
c <- ggtree(tree_molecular_unrestricted_pruned, aes(color = mean_growth_rate)) + #geom_tippoint(color = sim_colors["unrestricted"][[1]], alpha = 0.8, size = 2) +
    geom_tippoint(aes(color = mean_growth_rate), alpha = 0.8, size = 2) +
    theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_color_viridis(limits = growth_rate_limits) + labs(x = "Genetic distance") +
    theme(axis.line.x = element_line(arrow = grid::arrow(length = unit(0.2, "cm"),
                                                         ends = "last",
                                                         type = "closed")))+ coord_cartesian(clip = "off") +
    geom_text(aes(label=mark), size = 6, nudge_x = 0.6, nudge_y = -0.6, color = "black") +
    geom_rootedge(color = "#939393")

#new_lim <- ggplot_build(c)$layout$panel_scales_y[[1]]$range$range + c(-1, 0)
#c <- c + coord_cartesian(clip = "off", ylim = new_lim, xlim = genetictree_limits)
c

# unlabeled tree for published figure
c_unlabeled <- ggtree(tree_molecular_unrestricted_pruned, aes(color = mean_growth_rate)) + #geom_tippoint(color = sim_colors["unrestricted"][[1]], alpha = 0.8, size = 2) +
    geom_tippoint(aes(color = mean_growth_rate), alpha = 0.8, size = 2) +
    theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_color_viridis(limits = growth_rate_limits) +
    geom_rootedge(color = "#939393")

c_unlabeled

ggsave(file = "molecular_tree_growth_rate_unrestricted_large.png", plot = c, path = figures_dir, height = 6, width = 5)
ggsave(file = "molecular_tree_growth_rate_unrestricted_unlabeled_large.png", plot = c_unlabeled, path = figures_dir, height = 6, width = 5)


##### FIGURES 1I + 1H Boundary versus Unrestricted growth variance plots #####

#tree_files <- paste("manuscript/analysis/simtrees/sampled_cells_death_rate_validation_pop_1000_dr_", dr_vec, "_diversified_sampling_timetree.rds", sep = "")
#tree_files <- list.files(path = "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/manuscript/analysis/simtrees", pattern = "*_diversified_sampling_timetree.rds", full.names = TRUE)

#tree_files <- list.files(path = "../eden/simtrees", pattern = "*_diversified_sampling_timetree.rds", full.names = TRUE)

# Move to  cluster script
#Local directory 
# sim_cells_dr_0.1_files <- paste("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation/",
#                                 list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation",
#                                            pattern = "*i_[0-9]+_dr_0.10.csv"), sep = "")

# sim_cells_dr_0.1_files <- list.files(path = "../eden/simulation_data",
#                                            pattern = "*i_[0-9]+_dr_0.10.csv", full.names = TRUE)
# get_tree_stats <- function(sim_cells_file) {
#     print(sim_cells_file)
#     all_cells <- read_csv(sim_cells_file)
#     alive_cells <- all_cells %>% 
#         filter_alive_cells %>% 
#         dplyr::mutate(n_muts = str_count(mutations, ",") + 1,
#                       seq_length = str_count(sequence, ",") + 1) 
#     
#     endpoint <- max(alive_cells$deathdate)
#     
#     alive_cells <- alive_cells %>% 
#         dplyr::mutate(clock_rate = n_muts/endpoint)
#     
#     # sampled_cells <- alive_cells %>%
#     #     sample_alive_cells(., n = 100, diversified_sampling = TRUE) %>%
#     #     filter(sampled)
#     # 
#     simtree <- convert_all_cells_to_tree_fast(all_cells)
#     simtree <- prune_simulated_tree(simtree, sampled_cells_indices = alive_cells$index)
#     #simtree <- prune_simulated_tree(simtree, sampled_cells_indices = sampled_cells$index)
# 
#     
#     i_extract <- regmatches(basename(sim_cells_file),
#                             gregexpr("(?<=i_)[[:digit:]]+", basename(sim_cells_file), perl = TRUE))[[1]]
#     
#     #compare variance in terminal branch lengths
#     terminal_branch_length_df <- data.frame("terminal_branch_length" = c(simtree@phylo$edge.length[which(simtree@phylo$edge[,2] <= length(simtree@phylo$tip.label))])) %>% 
#         dplyr::mutate(norm_terminal_branch_length = terminal_branch_length / endpoint) 
#     
#     norm_terminal_branch_length_variance <- var(terminal_branch_length_df$norm_terminal_branch_length) 
#     terminal_branch_length_variance <- var(terminal_branch_length_df$terminal_branch_length) 
#     clock_rate_variance <- var(alive_cells$clock_rate) 
#     #clock_rate_variance <- var(sampled_cells$clock_rate) 
# 
# 
# 
#     return(data.frame("terminal_branch_length_variance" = terminal_branch_length_variance,
#                       "norm_terminal_branch_length_variance" = norm_terminal_branch_length_variance,
#                       "clock_rate_variance" = clock_rate_variance, 
#                       "i" = i_extract,
#                       "model" = ifelse(grepl("pushing", sim_cells_file),
#                                        "unrestricted",
#                                        "boundary_driven")))
#     
# }
# 
# set.seed(3182)
# tree_terminal_branch_lengths <- purrr::map(sim_cells_dr_0.1_files, function(tf) get_tree_stats(tf)) %>%
#     bind_rows
# 
# #write to csv
# write_csv(tree_terminal_branch_lengths, file = "../eden/stats/terminal_bl_clock_rate_stats_full.csv")

#to skip computation time 
tree_terminal_branch_lengths <- read_tsv(file = "../eden/stats/terminal_bl_clock_rate_stats_full_large.tsv")

# values for legend

tbl_var_boundary_driven <- mean(tree_terminal_branch_lengths$norm_terminal_branch_length_variance[tree_terminal_branch_lengths$model == "boundary_driven"])
tbl_var_unrestricted <- mean(tree_terminal_branch_lengths$norm_terminal_branch_length_variance[tree_terminal_branch_lengths$model == "unrestricted"])

clock_rate_var_boundary_driven <- mean(tree_terminal_branch_lengths$clock_rate_variance[tree_terminal_branch_lengths$model == "boundary_driven"])
clock_rate_var_var_unrestricted <- mean(tree_terminal_branch_lengths$clock_rate_variance[tree_terminal_branch_lengths$model == "unrestricted"])


## Terminal branch length plot for example boundary-driven tumor (Figure 1H)
d <- tree_terminal_branch_lengths %>% 
    ggplot(., aes(y = norm_terminal_branch_length_variance, x = model, fill = model)) + 
    geom_violin() + 
    #ggridges::geom_density_ridges(data = tree_terminal_branch_lengths_boundary_driven, rel_min_height = 0.1, binwidth = 2) +
    scale_fill_manual(values = sim_colors) +
    theme_classic() +
    xlab("") + ylab("Variance(Terminal branch length)") +
    stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
                 colour = "white",
                 size = 0.25) +
    theme(legend.position = "none") +
    theme(text=element_text(size=15)) +
    scale_x_discrete(labels=c("boundary_driven" = "Boundary-\nDriven", "unrestricted" = "Unrestricted"))
    #ggtitle("variance_model_comparison_term_branchh")
d

## Terminal branch length plot for example unrestricted tumor (Figure 1I)
e <- tree_terminal_branch_lengths %>% 
    ggplot(., aes(y = clock_rate_variance, x = model, fill = model)) + 
    geom_violin() + 
    scale_fill_manual(values = sim_colors) +
    theme_classic() + ylab("Variance(Lineage clock rate)") + xlab("") +
    stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
                 colour = "white",
                 size = 0.25) +
    theme(legend.position = "none") +
    theme(text=element_text(size=15)) +
    scale_x_discrete(labels=c("boundary_driven" = "Boundary-\nDriven", "unrestricted" = "Unrestricted"))
    #ggtitle("Lineage clock rate")
e    
#e_d_comb <- grid.arrange(d, e, ncol=2)
#ggsave(filename = "../figures/variance_model_comparison_term_branch_clock_rate.png", e_d_comb, height = 5, width = 7)
ggsave(filename = "../figures/variance_model_comparison_term_branch_large.png", d, height = 5, width = 5)
ggsave(filename = "../figures/variance_model_comparison_clock_rate_large.png", e, height = 5, width = 5)

 
###### INSETS ######

### terminal branch length
tree_boundary_driven_alive <- prune_simulated_tree(tree_boundary_driven, alive_cells_boundary_driven$index)
terminal_branch_length_boundary_driven_df <- data.frame("terminal_branch_length" = c(tree_boundary_driven_alive@phylo$edge.length[which(tree_boundary_driven_alive@phylo$edge[,2] <= length(tree_boundary_driven_alive@phylo$tip.label))])) %>% 
    dplyr::mutate(norm_terminal_branch_length = terminal_branch_length / max(tree_boundary_driven_alive@data$deathdate))

variance_termnial_branch_length_boundary_driven <- var(terminal_branch_length_boundary_driven_df$norm_terminal_branch_length)


#unrestricted growth model
tree_unrestricted_alive <- prune_simulated_tree(tree_unrestricted, alive_cells_unrestricted$index)
terminal_branch_length_unrestricted_df <- data.frame("terminal_branch_length" = c(tree_unrestricted_alive@phylo$edge.length[which(tree_unrestricted_alive@phylo$edge[,2] <= length(tree_unrestricted_alive@phylo$tip.label))])) %>% 
    dplyr::mutate(norm_terminal_branch_length = terminal_branch_length / max(tree_unrestricted_alive@data$deathdate))

variance_termnial_branch_length_unrestricted <- var(terminal_branch_length_unrestricted_df$norm_terminal_branch_length)

#normalize axes

terminal_branch_max <- max(terminal_branch_length_unrestricted_df$norm_terminal_branch_length, terminal_branch_length_boundary_driven_df$norm_terminal_branch_length)
a <- ggplot(terminal_branch_length_boundary_driven_df, aes(x=norm_terminal_branch_length)) + geom_density(fill = sim_colors["boundary_driven"]) +
    theme_classic() + xlab("Terminal branch length") +
    # annotate("text",
    #          x = 0.75, y = 1.6, size = 6, color = sim_colors["boundary_driven"],
    #          label = paste0("Var(TBL in Fig 1B) = \n", format(variance_termnial_branch_length_boundary_driven, digits = 2))) +
    theme(axis.title.y=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(), 
          axis.line.y=element_blank()) +
    theme(text=element_text(size=15))
    #xlim(c(0, terminal_branch_max ))
a
b <- ggplot(terminal_branch_length_unrestricted_df, aes(x=norm_terminal_branch_length)) + geom_density(fill = sim_colors["unrestricted"]) +
    theme_classic() + xlab("Terminal branch length") +
    # annotate("text",
    #          x = 0.15, y = 7, size = 6, color = sim_colors["unrestricted"],
    #          label = paste0("Var(TBL in Fig 1E) = \n", format(variance_termnial_branch_length_unrestricted, digits = 2))) +
    theme(axis.title.y=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.line.y=element_blank()) +
    theme(text=element_text(size=15))
    #xlim(c(0, terminal_branch_max ))
a_b <- plot_grid(a, b, ncol = 1)
a_b
ggsave(filename = "../figures/terminal_branch_length_variance_densities_large.png", a_b, height = 4, width = 5)


###### Clock rate variances for individual examples #######

alive_cells_boundary_driven <- alive_cells_boundary_driven %>% 
    filter_alive_cells %>% 
    dplyr::mutate(n_muts = str_count(mutations, ",") + 1) %>% 
    dplyr::mutate(clock_rate = n_muts/max(alive_cells_boundary_driven$deathdate))

alive_cells_unrestricted <- alive_cells_unrestricted %>% 
    filter_alive_cells %>% 
    dplyr::mutate(n_muts = str_count(mutations, ",") + 1) %>% 
    dplyr::mutate(clock_rate = n_muts/max(alive_cells_unrestricted$deathdate))

variance_clock_rate_boundary_driven <- var(alive_cells_boundary_driven$clock_rate)
print(variance_clock_rate_boundary_driven)
variance_clock_rate_unrestricted <- var(alive_cells_unrestricted$clock_rate)
print(variance_clock_rate_unrestricted)
clock_rate_max <- max(alive_cells_boundary_driven$clock_rate, alive_cells_unrestricted$clock_rate)

a <- ggplot(alive_cells_boundary_driven, aes(x=clock_rate)) + geom_density(fill = sim_colors["boundary_driven"]) +
    theme_classic() + xlab("Clock rate") +
    # annotate("text",
    #          x = 0.4, y = 1.2, size = 6, color = sim_colors["boundary_driven"],
    #          label = paste0("Var(GD in Fig 1C) = \n", format(variance_clock_rate_boundary_driven, digits = 2))) +
    theme(axis.title.y=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(), 
          axis.line.y=element_blank()) +
    xlim(c(0, clock_rate_max )) +
    theme(text=element_text(size=15))
a
b <- ggplot(alive_cells_unrestricted, aes(x=clock_rate)) + geom_density(fill = sim_colors["unrestricted"]) +
    theme_classic() + xlab("Clock rate") +
    # annotate("text",
    #          x = 0.25, y = 1.2, size = 6,
    #          label = paste0("Var(GD in Fig 1F) = \n", format(variance_clock_rate_unrestricted, digits = 2)), 
    #          color = sim_colors["unrestricted"]) +
    theme(axis.title.y=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(), 
          axis.line.y=element_blank()) +
    xlim(c(0, clock_rate_max )) +
    theme(text=element_text(size = 15))
b
a_b <- plot_grid(a, b, ncol = 1)
a_b
ggsave(filename = "../figures/clock_rate_variance_densities_large.png", a_b, height = 4, width = 5)



##### Edge versus center labeled boundary driven tumor #####
### For inset in Figure 2A
alive_cells_boundary_driven <- mark_boundary(cells_to_mark = alive_cells_boundary_driven, alive_cells = alive_cells_boundary_driven)


edge_center_colors <- c("edge" = "#89352F", 
                        "center" = "#A2D2E2")
a_edge_center <- ggplot(alive_cells_boundary_driven, aes(x = norm_locx, y = norm_locy, color = ifelse(est_edge == 1, "edge", "center" ))) +
    geom_point(size = 1.5, alpha = 1) +
    #geom_point(data = sampled_cells_boundary_driven, size = 3, alpha = 1) +
    
    #geom_point(data = sampled_cells_boundary_driven, size = 3, alpha = 1, colour="black",pch=21, aes(fill = mean_growth_rate)) +
    #geom_point(data = sampled_cells_boundary_driven, aes(x = norm_locx, y = norm_locy, color = dist_from_edge), size = 2) +
    theme_void() +
    scale_color_manual(values = edge_center_colors) +
    theme(legend.position = "none")
a_edge_center
ggsave(file = "tumor_sim_edge_center_boundary_driven_unlabeled_large.png", plot = a_edge_center, path = figures_dir, height = 4, width = 4)

