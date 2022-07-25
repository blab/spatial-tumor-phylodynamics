#visualize_3d_trees.R

setwd("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/")
library(tidyverse)
library(coda)
library(qdapRegex)
colors_edge_center<- tumortree::get_color_palette(names = c("edge", "center"))
colors_loc <- tumortree::get_color_palette(names = c("edge", "center"))
names(colors_loc) <- c("loc1", "loc0")
mcc_tree_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/trees/3D_neut_bdg/diversified_100",
                             pattern="mcc", 
                             full.names = TRUE)
mcc_tree_files <- mcc_tree_files[grepl("d0.1", mcc_tree_files)]

terminal_branch_ratios_df_3d <- data.frame()
for (mcc_file in mcc_tree_files) {
    
    print(mcc_file)
    mcc_tree <- read.beast(file=mcc_file)
    
    log_file <- paste0("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/logs/3D_neut_bdg/diversified_100/",
                       gsub("_mcc.tree", "", basename(mcc_file)), ".log")
    
    log <- readLog(log_file)
    ess <- coda::effectiveSize(log)
    
    if (min(c(ess["birthRateCanonical.1"], ess["birthRateCanonical.0"])) < 200) {
        next
    }
    
    # dr_extract <- regmatches(basename(mcc_file),
    #                           gregexpr("(?<=_d).+[0-9]+", basename(mcc_file), perl = TRUE))[[1]]
    
    dr_extract <- qdapRegex::ex_between(basename(mcc_file), "_d", "_")[[1]]

    rep_extract <- regmatches(basename(mcc_file),
                              gregexpr("(?<=_s)[0-9]+", basename(mcc_file), perl = TRUE))[[1]]
    
    #plot tree
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
    
    treeplot_pie <- ggtree::inset(treeplot, pies, width = 0.06, height = 0.06) + theme(legend.position = "none") #+ ggtitle(paste0(tumor_extract, migration_model, rep_extract, sep = " "))
    
    #plot violin plots
    terminal_branch_length_df <- data.frame("state" = mcc_tree@data$type[match(1:100,mcc_tree@data$node)],
               "terminal_branch_length" = mcc_tree@phylo$edge.length[match(1:100,mcc_tree@phylo$edge[,2])])
    
    summary_terminal_branch <- terminal_branch_length_df %>% 
        group_by(state) %>% 
        summarise(mean_terminal_branch_length = mean(terminal_branch_length))
    
    ratio_means <- summary_terminal_branch$mean_terminal_branch_length[summary_terminal_branch$state == "loc1"] /summary_terminal_branch$mean_terminal_branch_length[summary_terminal_branch$state == "loc0"]
    
    terminal_branch_ratios_df_3d <- bind_rows(list(terminal_branch_ratios_df_3d, data.frame("dr" = dr_extract, "terminal_branch_length_ratio" = ratio_means )))
    terminal_branch_length_plot <- ggplot(terminal_branch_length_df, aes(x=state,y=terminal_branch_length, fill = state))  +
        geom_violin(alpha=0.8) + theme_classic()  +
        scale_fill_manual(values = colors_loc) +
        stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
                     colour = "black",
                     size = 0.25) +
        xlab("") + ylab("Terminal branch length") +
        theme(legend.position = "none")  +
        # annotate("text", x = 1, y = 12, size = 4,
        #          label = paste0("frac('Mean center TBL', 'Mean edge TBL') == ", format(mean_ratio_branch_length_unrestricted, digits = 2)), 
        #          color = "black", parse = TRUE) +
        theme(text=element_text(size=15))
    
    violin_fig_file <- paste0("manuscript/figures/3d_physicell_mcc_trees/3d_physicell_d", dr_extract, "_s",
                       rep_extract, "_",
                       "terminal_branch_plot.png")
    
    ggsave(plot=terminal_branch_length_plot,
           file=violin_fig_file, height = 5, width = 5)
    
    fig_file <- paste0("manuscript/figures/3d_physicell_mcc_trees/3d_physicell_d", dr_extract, "_s",
                       rep_extract, "_",
                       "mcc_tree.png")
    print(fig_file)
    ggsave(plot=treeplot_pie,
           file=fig_file, height = 9, width = 7)
}

mcc_tree_files_2d <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/trees/2D_neut_bdg/diversified_100",
                             pattern="mcc", 
                             full.names = TRUE)

mcc_tree_files_2d <- mcc_tree_files_2d[grepl("d0.1", mcc_tree_files_2d)]
mcc_tree_files_2d <- mcc_tree_files_2d[! grepl("n_[0-9]", mcc_tree_files_2d)]
terminal_branch_ratios_df_2d <- data.frame()
for (mcc_file in mcc_tree_files_2d ) {
    
    print(mcc_file)
    mcc_tree <- read.beast(file=mcc_file)
    
    log_file <- paste0("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/logs/2D_neut_bdg/diversified_100/",
                       gsub("_mcc.tree", "", basename(mcc_file)), ".log")
    
    log <- readLog(log_file)
    ess <- coda::effectiveSize(log)
    
    if (min(c(ess["birthRateCanonical.1"], ess["birthRateCanonical.0"])) < 200) {
        next
    }
    
    # dr_extract <- regmatches(basename(mcc_file),
    #                           gregexpr("(?<=_d).+[0-9]+", basename(mcc_file), perl = TRUE))[[1]]
    
    dr_extract <- qdapRegex::ex_between(basename(mcc_file), "_d", "_")[[1]]
    
    rep_extract <- regmatches(basename(mcc_file),
                              gregexpr("(?<=_s)[0-9]+", basename(mcc_file), perl = TRUE))[[1]]
    
    #plot tree
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
    
    treeplot_pie <- ggtree::inset(treeplot, pies, width = 0.06, height = 0.06) + theme(legend.position = "none") #+ ggtitle(paste0(tumor_extract, migration_model, rep_extract, sep = " "))
    
    #plot violin plots
    terminal_branch_length_df <- data.frame("state" = mcc_tree@data$type[match(1:100,mcc_tree@data$node)],
                                            "terminal_branch_length" = mcc_tree@phylo$edge.length[match(1:100,mcc_tree@phylo$edge[,2])])
    
    summary_terminal_branch <- terminal_branch_length_df %>% 
        group_by(state) %>% 
        summarise(mean_terminal_branch_length = mean(terminal_branch_length))
    
    ratio_means <- summary_terminal_branch$mean_terminal_branch_length[summary_terminal_branch$state == "loc1"] /summary_terminal_branch$mean_terminal_branch_length[summary_terminal_branch$state == "loc0"]
    
    terminal_branch_ratios_df_2d <- bind_rows(list(terminal_branch_ratios_df_2d, data.frame("dr" = dr_extract, "terminal_branch_length_ratio" = ratio_means )))
    terminal_branch_length_plot <- ggplot(terminal_branch_length_df, aes(x=state,y=terminal_branch_length, fill = state))  +
        geom_violin(alpha=0.8) + theme_classic()  +
        scale_fill_manual(values = colors_loc) +
        stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
                     colour = "black",
                     size = 0.25) +
        xlab("") + ylab("Terminal branch length") +
        theme(legend.position = "none")  +
        # annotate("text", x = 1, y = 12, size = 4,
        #          label = paste0("frac('Mean center TBL', 'Mean edge TBL') == ", format(mean_ratio_branch_length_unrestricted, digits = 2)), 
        #          color = "black", parse = TRUE) +
        theme(text=element_text(size=15))
    
    violin_fig_file <- paste0("manuscript/figures/2d_physicell_mcc_trees/2d_physicell_d", dr_extract, "_s",
                              rep_extract, "_",
                              "terminal_branch_plot.png")
    
    ggsave(plot=terminal_branch_length_plot,
           file=violin_fig_file, height = 5, width = 5)
    
    fig_file <- paste0("manuscript/figures/2d_physicell_mcc_trees/2d_physicell_d", dr_extract, "_s",
                       rep_extract, "_",
                       "mcc_tree.png")
    print(fig_file)
    ggsave(plot=treeplot_pie,
           file=fig_file, height = 9, width = 7)
}

##get represenative plots for each dimension

mcc_tree_files_2d_rep <- mcc_tree_files_2d[grepl("s63293", mcc_tree_files_2d)]
mcc_tree_files_3d_rep <- mcc_tree_files[grepl("s73314", mcc_tree_files)]

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

treeplot_pie <- ggtree::inset(treeplot, pies, width = 0.05, height = 0.05) + theme(legend.position = "none") #+ ggtitle(paste0(tumor_extract, migration_model, rep_extract, sep = " "))
ggsave(plot=treeplot_pie,
       file="manuscript/figures/2d_physicell_example_mcc_tree.png", height = 7, width = 7)
#plot violin plots
terminal_branch_length_df <- data.frame("state" = mcc_tree@data$type[match(1:100,mcc_tree@data$node)],
                                        "terminal_branch_length" = mcc_tree@phylo$edge.length[match(1:100,mcc_tree@phylo$edge[,2])])

summary_terminal_branch <- terminal_branch_length_df %>% 
    group_by(state) %>% 
    summarise(mean_terminal_branch_length = mean(terminal_branch_length))

ratio_means <- summary_terminal_branch$mean_terminal_branch_length[summary_terminal_branch$state == "loc1"] /summary_terminal_branch$mean_terminal_branch_length[summary_terminal_branch$state == "loc0"]

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
             label = paste0("frac('Mean center TBL', 'Mean edge TBL') == ", format(1/ratio_means, digits = 2)),
             color = "black", parse = TRUE) +
    theme(text=element_text(size=15)) +
    scale_x_discrete(labels = c('center','edge'))

terminal_branch_length_plot
ggsave(plot=terminal_branch_length_plot,
       file="manuscript/figures/2d_physicell_example_terminal_branch_lengths.png", height = 5, width = 5)




### summary plot of terminal branch length ratios #####

terminal_branch_ratios_df_3d <- terminal_branch_ratios_df_3d %>% 
    add_column("sim" = "3D")

terminal_branch_ratios_df_2d <- terminal_branch_ratios_df_2d %>% 
    add_column("sim" = "2D")

terminal_branch_ratios_df_comb <- bind_rows(list(terminal_branch_ratios_df_2d,terminal_branch_ratios_df_3d))

sim_colors <- c("#254E70", "#D0CF92")
dim_comparison_plot <- ggplot(terminal_branch_ratios_df_comb, aes(x=1/terminal_branch_length_ratio, fill = sim)) + 
    geom_histogram(color = "black", bins = 20) + theme_classic() + #scale_fill_brewer(type ="qual", palette = 7) +
    xlab("Mean center / edge terminal branch length ratio") + scale_fill_manual(values=sim_colors) + theme(legend.position = "none") +
    theme(text=element_text(size = 15))

dim_comparison_plot


ggsave(file = "manuscript/figures/dim_comparison_plot.png", plot=dim_comparison_plot, height = 5, width = 5)

##### 3D EXAMPLE #####

mcc_tree <- read.beast(mcc_tree_files_3d_rep)
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

treeplot_pie <- ggtree::inset(treeplot, pies, width = 0.05, height = 0.05) + theme(legend.position = "none") #+ ggtitle(paste0(tumor_extract, migration_model, rep_extract, sep = " "))
ggsave(plot=treeplot_pie,
       file="manuscript/figures/3d_physicell_example_mcc_tree.png", height = 7, width = 7)
#plot violin plots
terminal_branch_length_df <- data.frame("state" = mcc_tree@data$type[match(1:100,mcc_tree@data$node)],
                                        "terminal_branch_length" = mcc_tree@phylo$edge.length[match(1:100,mcc_tree@phylo$edge[,2])])

summary_terminal_branch <- terminal_branch_length_df %>% 
    group_by(state) %>% 
    summarise(mean_terminal_branch_length = mean(terminal_branch_length))

ratio_means <- summary_terminal_branch$mean_terminal_branch_length[summary_terminal_branch$state == "loc1"] /summary_terminal_branch$mean_terminal_branch_length[summary_terminal_branch$state == "loc0"]

terminal_branch_length_plot <-  terminal_branch_length_df %>% 
    ggplot(., aes(x=state,y=terminal_branch_length, fill = state))  +
    geom_violin(alpha=0.8) + theme_classic()  +
    scale_fill_manual(values = colors_loc) +
    stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
                 colour = "black",
                 size = 0.25) +
    xlab("") + ylab("Terminal branch length") +
    theme(legend.position = "none")  +
    annotate("text", x = 2, y = 4500, size = 4 ,
             label = paste0("frac('Mean center TBL', 'Mean edge TBL') == ", format(1/ratio_means, digits = 2)),
             color = "black", parse = TRUE) +
    theme(text=element_text(size=15)) +
    scale_x_discrete(labels = c('center','edge'))

terminal_branch_length_plot
ggsave(plot=terminal_branch_length_plot,
       file="manuscript/figures/3d_physicell_example_terminal_branch_lengths.png", height = 5, width = 5)



