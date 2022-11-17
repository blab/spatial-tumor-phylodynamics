#figureS7.R

library(tumortree)
library(tidyverse)
##### PLOT POSTERIORS ####

all_logs_birthRate_df <- read_tsv("../li-application/stats/sdevo_estimates_summary.tsv")
#for supplement compare clock models (both unidirectional + oldstates)

## Tumor 1
t1_wgs_posteriors_plot_oldstates_clock_comparison <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           states == "oldstates",
           tumor == "T1") %>% 
    mutate("category" = paste(subset, state)) %>% 
    ggplot(., aes(x=birthRate, group = category), color = "black") +
    geom_density(aes(fill = state), alpha=0.7) +
    facet_wrap(~clock_model, nrow = 1) +
    theme_classic() + scale_fill_manual(values=colors_loc) +
    theme(text=element_text(size=20))+
    xlab("Estimated birth rate") +
    theme(legend.position = "none") +ylab("") 

t1_wgs_posteriors_plot_oldstates_clock_comparison
ggsave(plot=t1_wgs_posteriors_plot_oldstates_clock_comparison,
       file ="../figures/t1_li_wgs_posteriors_oldstates_clock_comparison.png", height = 4, width = 6)

## Tumor 2
t2_wgs_posteriors_plot_oldstates_clock_comparison <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           states == "oldstates",
           tumor == "T2") %>% 
    mutate("category" = paste(subset, state)) %>% 
    ggplot(., aes(x=birthRate, group = category), color = "black") +
    geom_density(aes(fill = state), alpha=0.8) +
    facet_wrap(~clock_model, nrow = 1) +
    theme_classic() + scale_fill_manual(values=colors_loc) +
    theme(text=element_text(size=20))+
    xlab("Estimated birth rate") +
    theme(legend.position = "none") +ylab("") 

t2_wgs_posteriors_plot_oldstates_clock_comparison

ggsave(plot=t2_wgs_posteriors_plot_oldstates_clock_comparison,
       file ="../figures/t2_li_wgs_posteriors_oldstates_clock_comparison.png", height = 4, width = 6)


# Get ratios 

#Tumor 1
t1_wgs_clock_comparison_ratios <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           states == "oldstates",
           tumor == "T1") %>% 
    mutate(clock_model = factor(clock_model, levels = c("state-dependent", "strict"))) %>% 
    #mutate("category" = paste0(clock_model, subset)) %>% 
    group_by(clock_model, subset) %>% 
    summarise(mean = mean(birthRateRatio),
              hpd_lower = hdi(birthRateRatio, credMass = 0.95)[1],
              hpd_upper = hdi(birthRateRatio, credMass = 0.95)[2]) %>% 
    ggplot(., aes(x=subset, y=mean), color = "black") +
    geom_point(size =2) +
    facet_wrap(~clock_model, ncol = 2) +
    geom_errorbar(aes(ymin = hpd_lower, ymax= hpd_upper), width = 0, size = 1) +
    theme_classic() +
    
    theme(text=element_text(size=15))+
    ylab("") +
    xlab("") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    #coord_flip() +
    theme(strip.text = element_blank()) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())


t1_wgs_clock_comparison_ratios
ggsave(plot=t1_wgs_clock_comparison_ratios,
       file ="../figures/t1_li_wgs_posteriors_ratios_clock_comparison.png", height = 3, width = 3)


#Tumor 2
t2_wgs_clock_comparison_ratios <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           states == "oldstates",
           tumor == "T2") %>% 
    mutate(clock_model = factor(clock_model, levels = c("state-dependent", "strict"))) %>% 
    group_by(clock_model, subset) %>% 
    summarise(mean = mean(birthRateRatio),
              hpd_lower = hdi(birthRateRatio, credMass = 0.95)[1],
              hpd_upper = hdi(birthRateRatio, credMass = 0.95)[2]) %>% 
    ggplot(., aes(x=subset, y=mean), color = "black") +
    geom_point(size =2) +
    geom_errorbar(aes(ymin = hpd_lower, ymax= hpd_upper), width = 0, size = 1) +
    theme_classic() +
    
    theme(text=element_text(size=15))+
    ylab("") +
    xlab("") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    #coord_flip() +
    theme(strip.text = element_blank()) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank()) +
    facet_wrap(~clock_model, ncol = 2)


t2_wgs_clock_comparison_ratios
ggsave(plot=t2_wgs_clock_comparison_ratios,
       file ="../figures/t2_li_wgs_posteriors_ratios_clock_comparison.png", height = 3, width = 3)

########## COMPARE ANCESTRAL NODE TIMINGS ###########

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

T1_mcc_tree <- read.beast("../li-application/combined/T1_wgs_oristates_unidir_1_state.HCCtumor.typed.node.tree")
T2_mcc_tree <- read.beast("../li-application/combined/T2_wgs_oristates_unidir_1_state.HCCtumor.typed.node.tree")

T1_strict_clock_mcc_tree <- read.beast("../li-application/combined/T1_wgs_oristates_unidir_1_strict.HCCtumor.typed.node.tree")
T2_strict_clock_mcc_tree <- read.beast("../li-application/combined/T2_wgs_oristates_unidir_1_strict.HCCtumor.typed.node.tree")

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
    geom_point(size = 4, alpha =0.8, shape = 21) + theme_classic() +
    scale_fill_manual(values = colors_loc) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_y_continuous(trans='log10') +
    scale_x_continuous(trans='log10')+
    theme(text = element_text(size = 20)) +
    xlab("log(Node height with state-dependent clock)") +
    ylab("log(Node height with strict clock)") +
    theme(legend.position = c(0.25, 0.85)) +
    labs(size = "Node posterior") +
    guides(fill = "none") #+

t1_node_height_comparison_plot

ggsave(plot=t1_node_height_comparison_plot,file ="../figures/t1_clocks_height_comparison_plot.png", height = 5, width = 5)

t2_node_height_comparison_plot <- ggplot(t2_node_height_comparison_df, aes(x=as.numeric(height), y=strict_clock_node_height,
                                                                           fill= type)) +
    geom_point(size = 4, alpha =0.8, shape = 21) + theme_classic() +
    scale_fill_manual(values = colors_loc) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
    scale_y_continuous(trans='log10') +
    scale_x_continuous(trans='log10') +
    theme(text = element_text(size = 20)) +
    xlab("log(Node height with state-dependent clock)") +
    ylab("log(Node height with strict clock)") +
    scale_size(limits = c(0,1)) +
    theme(legend.position = c(0.25, 0.85)) +
    labs(size = "Node posterior") +
    guides(fill = "none") 

t2_node_height_comparison_plot    
ggsave(plot=t2_node_height_comparison_plot,file ="../figures/t2_clocks_height_comparison_plot.png", height = 5, width = 5)

