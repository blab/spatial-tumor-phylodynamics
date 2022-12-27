#figureS5.R

## Comparison of 2D and 3D branching signals


##### SETUP #####
library(tidyverse)
library(coda)
library(qdapRegex)
library(ggtree)
library(treeio)

colors_edge_center<- tumortree::get_color_palette(names = c("edge", "center"))
colors_loc <- tumortree::get_color_palette(names = c("edge", "center"))
names(colors_loc) <- c("loc1", "loc0")

# Record of local directory
# mcc_tree_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/trees/3D_neut_bdg/diversified_100",
#                              pattern="mcc",
#                              full.names = TRUE)



#Read in stats generated in dimension_comparison_stats.R
terminal_branch_ratios_df_comb <- read_csv("../physicell/stats/terminal_branch_dim_comparison.csv")

##### Plotting branching ratio histogram #####
sim_colors <- c("#fba500", "#244e70")


dim_comparison_plot <- ggplot(terminal_branch_ratios_df_comb, aes(x=1/terminal_branch_length_ratio, fill = as.factor(sim))) + 
    geom_histogram(color = "black", bins = 20) + theme_classic() +
    xlab("Mean center / edge terminal branch length ratio") + scale_fill_manual(values = sim_colors) + theme(legend.position = "none") +
    theme(text=element_text(size = 20))

dim_comparison_plot


ggsave(file = "../figures/dim_comparison_plot.png", plot=dim_comparison_plot, height = 5, width = 5)

## Get representative plots for each dimension
#mcc_tree_files_2d_rep <- mcc_tree_files_2d[grepl("s63293", mcc_tree_files_2d)]
#mcc_tree_files_3d_rep <- mcc_tree_files[grepl("s73314", mcc_tree_files)]

mcc_tree_files_2d_rep <- "../physicell/trees/2D_neut_bdg/diversified_100/sampconfig_m0_w1_d0.1_t1_mg1_mm1_l2e+08_i7_s63293_mcc.tree"
mcc_tree_files_3d_rep <- "../physicell/trees/3D_neut_bdg/diversified_100/sampconfig_m0_w1_d0.1_t1_mg1_mm1_l2e+08_i8_s73314_mcc.tree"

#### 2D Example TREE ####
mcc_tree <- read.beast(mcc_tree_files_2d_rep)

anc_state_nodes <- mcc_tree@data %>% 
    #filter(node > length(T1_mcc_tree@phylo$tip.label)) %>% 
    mutate(type.prob = as.numeric(type.prob)) %>% 
    dplyr::mutate(edge = type.prob * as.integer(type == "loc1") + (1 - type.prob)*as.integer(type == "loc0")) %>% 
    dplyr::mutate(center = type.prob * as.integer(type == "loc0") + (1 - type.prob)*as.integer(type == "loc1")) %>% 
    dplyr::select(center, edge, type, node)

pies <- nodepie(anc_state_nodes, cols=1:2, alpha=0.8, color = colors_edge_center)

treeplot <- ggtree(mcc_tree, color = "darkgrey", size = 1) +
    #geom_point(size = 2) +
    scale_color_manual(values = colors_loc) +
    theme(legend.position = "none")

treeplot_pie <- ggtree::inset(treeplot, pies, width = 0.025, height = 0.05) + theme(legend.position = "none") +
    geom_nodelab(aes(label = ifelse(round(as.numeric(posterior),2) < 0.99,
                                    round(as.numeric(posterior), 2), "")),
                 nudge_x = -40, nudge_y = 1, size  = 3, hjust='right')

treeplot_pie 
ggsave(plot=treeplot_pie,
       file="../figures/2d_physicell_example_mcc_tree.png", height = 7, width = 7)

#plot violin plots
terminal_branch_length_df <- data.frame("state" = mcc_tree@data$type[match(1:100,mcc_tree@data$node)],
                                        "terminal_branch_length" = mcc_tree@phylo$edge.length[match(1:100,mcc_tree@phylo$edge[,2])])

summary_terminal_branch <- terminal_branch_length_df %>% 
    group_by(state) %>% 
    summarise(mean_terminal_branch_length = mean(terminal_branch_length), 
              "N" = n())

ratio_means <- summary_terminal_branch$mean_terminal_branch_length[summary_terminal_branch$state == "loc0"] /summary_terminal_branch$mean_terminal_branch_length[summary_terminal_branch$state == "loc1"]

terminal_branch_length_plot <-  terminal_branch_length_df %>% 
    ggplot(., aes(x=state,y=terminal_branch_length, fill = state))  +
    geom_violin(alpha=0.8) + theme_classic()  +
    scale_fill_manual(values = colors_loc) +
    stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
                 colour = "black",
                 size = 0.25) +
    xlab("") + ylab("Terminal branch length") +
    theme(legend.position = "none")  +
    annotate("text", x = 2, y = 10000, size = 4 ,
             label = paste0("frac('Mean center TBL', 'Mean edge TBL') == ", format(ratio_means, digits = 2)),
             color = "black", parse = TRUE) +
    theme(text=element_text(size=15)) +
    scale_x_discrete(labels = c('center','edge'))

terminal_branch_length_plot
ggsave(plot=terminal_branch_length_plot,
       file="../figures/2d_physicell_example_terminal_branch_lengths.png", height = 5, width = 5)


##### 3D EXAMPLE #####

mcc_tree_3d <- read.beast(mcc_tree_files_3d_rep)
anc_state_nodes_3d <- mcc_tree_3d@data %>% 
    #filter(node > length(T1_mcc_tree@phylo$tip.label)) %>% 
    mutate(type.prob = as.numeric(type.prob)) %>% 
    dplyr::mutate(edge = type.prob * as.integer(type == "loc1") + (1 - type.prob)*as.integer(type == "loc0")) %>% 
    dplyr::mutate(center = type.prob * as.integer(type == "loc0") + (1 - type.prob)*as.integer(type == "loc1")) %>% 
    dplyr::select(center, edge, type, node)

pies_3d <- nodepie(anc_state_nodes_3d, cols=1:2, alpha=0.8, color = colors_edge_center)

treeplot_3d <- ggtree(mcc_tree_3d, color = "darkgrey", size = 1) +
    #geom_point(size = 2) +
    scale_color_manual(values = colors_loc) +
    theme(legend.position = "none")

treeplot_pie_3d <- ggtree::inset(treeplot_3d, pies_3d, width = 0.025, height = 0.025) + theme(legend.position = "none") +
    geom_nodelab(aes(label = ifelse(round(as.numeric(posterior),2) < 0.99, round(as.numeric(posterior), 2), "")),
                 nudge_x = -40, nudge_y = 1, size  = 3, hjust='right')

ggsave(treeplot_pie_3d,
       file="../figures/3d_physicell_example_mcc_tree.png", height = 7, width = 7)
#plot violin plots
terminal_branch_length_df_3d <- data.frame("state" = mcc_tree_3d@data$type[match(1:100,mcc_tree_3d@data$node)],
                                        "terminal_branch_length" = mcc_tree_3d@phylo$edge.length[match(1:100,mcc_tree_3d@phylo$edge[,2])])

summary_terminal_branch_3d <- terminal_branch_length_df_3d %>% 
    group_by(state) %>% 
    summarise(mean_terminal_branch_length = mean(terminal_branch_length),
              "N" = n())

ratio_means_3d <- summary_terminal_branch_3d$mean_terminal_branch_length[summary_terminal_branch_3d$state == "loc0"] /summary_terminal_branch_3d$mean_terminal_branch_length[summary_terminal_branch_3d$state == "loc1"]


terminal_branch_length_plot_3d <-  terminal_branch_length_df_3d %>% 
    ggplot(., aes(x=state,y=terminal_branch_length, fill = state))  +
    geom_violin(alpha=0.8) + theme_classic()  +
    scale_fill_manual(values = colors_loc) +
    stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
                 colour = "black",
                 size = 0.25) +
    xlab("") + ylab("Terminal branch length") +
    theme(legend.position = "none")  +
    annotate("text", x = 2, y = 4500, size = 4 ,
             label = paste0("frac('Mean center TBL', 'Mean edge TBL') == ", format(ratio_means_3d, digits = 2)),
             color = "black", parse = TRUE) +
    theme(text=element_text(size=15)) +
    scale_x_discrete(labels = c('center','edge'))

terminal_branch_length_plot_3d
ggsave(plot=terminal_branch_length_plot_3d,
       file="../figures/3d_physicell_example_terminal_branch_lengths.png", height = 5, width = 5)


#
