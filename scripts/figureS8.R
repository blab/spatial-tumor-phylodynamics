#figureS8.R


library(treeio)
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

# log_files <- list.files(path = "../li-application/out",
#                         pattern="mcc.tree",
#                         full.names = TRUE)
li_mcc_tree_files <- list.files(path = "../li-application/combined",
                                pattern="T[1-2]_wgs_oristates_unidir_1_state.HCCtumor.typed.node.tree", 
                                full.names = TRUE)
li_mcc_tree_files <- li_mcc_tree_files[! grepl(".trees", li_mcc_tree_files )]
ling_mcc_tree_files <- list.files(path = "../ling-application/out",
                                  pattern="hcc-wes_unidir_state_comb.HCCtumor.typed.node.mcc.tree", 
                                  full.names = TRUE)
mcc_tree_files <- c(li_mcc_tree_files, ling_mcc_tree_files)

T1_mcc_tree <- read.beast("../li-application/combined/T1_wgs_oristates_unidir_1_state.HCCtumor.typed.node.tree")
T2_mcc_tree <- read.beast("../li-application/combined/T2_wgs_oristates_unidir_1_state.HCCtumor.typed.node.tree")
T3_mcc_tree <- read.beast("../ling-application/out/hcc-wes_unidir_state_comb.HCCtumor.typed.node.mcc.tree")

T1_strict_clock_mcc_tree <- read.beast("../li-application/combined/T1_wgs_oristates_unidir_1_strict.HCCtumor.typed.node.tree")
T2_strict_clock_mcc_tree <- read.beast("../li-application/combined/T2_wgs_oristates_unidir_1_strict.HCCtumor.typed.node.tree")
T3_strict_clock_mcc_tree <- read.beast("../ling-application/out/hcc-wes_unidir_strict_clock_comb.HCCtumor.typed.node.mcc.tree")

#get all internal nodes to compare for each tumor
t1_nodes <- seq((length(T1_mcc_tree@phylo$tip.label) + 1), length(T1_mcc_tree@phylo$tip.label) +T1_mcc_tree@phylo$Nnode)
t2_nodes <- seq((length(T2_mcc_tree@phylo$tip.label) + 1), length(T2_mcc_tree@phylo$tip.label) +T2_mcc_tree@phylo$Nnode)
t3_nodes <- seq((length(T3_mcc_tree@phylo$tip.label) + 1), length(T3_mcc_tree@phylo$tip.label) +T3_mcc_tree@phylo$Nnode)

#get target nodes for both tumors 

t1_strict_clock_nodes <- map_dbl(t1_nodes, function(n) get_corresponding_node(node = n, tree = T1_mcc_tree, comparison_tree = T1_strict_clock_mcc_tree))
t2_strict_clock_nodes <- map_dbl(t2_nodes, function(n) get_corresponding_node(node = n, tree = T2_mcc_tree, comparison_tree = T2_strict_clock_mcc_tree))
t3_strict_clock_nodes <- map_dbl(t3_nodes, function(n) get_corresponding_node(node = n, tree = T3_mcc_tree, comparison_tree = T3_strict_clock_mcc_tree))

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


#Tumor 3
t3_node_height_comparison_df <- T3_mcc_tree@data %>% 
    dplyr::filter(node %in% t3_nodes) 

t3_node_height_comparison_df$strict_clock_node <- t3_strict_clock_nodes[match(as.integer(t3_node_height_comparison_df$node), t3_nodes)]
t3_node_height_comparison_df$strict_clock_node_height <- as.numeric(T3_strict_clock_mcc_tree@data$height[match(t3_node_height_comparison_df$strict_clock_node, T3_strict_clock_mcc_tree@data$node)])

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

t3_node_height_comparison_plot <- ggplot(t3_node_height_comparison_df, aes(x=as.numeric(height), y=strict_clock_node_height,
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

t3_node_height_comparison_plot    
ggsave(plot=t3_node_height_comparison_plot,
       file ="../figures/t3_clocks_height_comparison_plot.png", height = 5, width = 5)
