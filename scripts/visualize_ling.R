# Ling et al visualizations
library(raster)
library(tidyverse)
library(tumortree)
#perimeter:12.8cm
#area = 11.22 cm^2

state_colors <- get_color_palette(names = c("edge", "center"))
wes_punches <- c("A61", "A25", "A66",  "A5", "B4", "B5", "B6", "C3", "C31", "C2", "D25", "D29", 
                 "D62", "D16", "D54", "D63", "D58", "A58", "B33", "B9", "B45")
## from webDigitizer2 (tar file contains original project)
ling_boundary_points <- read.csv("../ling-application/spatial_measurements/boundary_points.csv")
ling_boundary_points <- bind_rows(ling_boundary_points, ling_boundary_points[1,])

ling_punches_points <- read.csv("../ling-application/spatial_measurements/punches_points.csv") %>% 
    dplyr::filter(Punch %in% wes_punches)

find_closest_dist_to_boundary <- function(x,y, boundary_points) {
    return(min(raster::pointDistance(p1 = c(x,y), p2 = boundary_points, lonlat = FALSE)))
}

ling_punches_points$dist_to_boundary <- map2_dbl(ling_punches_points$X, ling_punches_points$Y, function(x,y) find_closest_dist_to_boundary(x,y, boundary_points = ling_boundary_points))


edge_cutoff <- 0.35 #cm


ling_punches_points <- ling_punches_points %>%
    mutate("edge" = as.integer(dist_to_boundary < edge_cutoff))

scale_bar_length <- 0.5 #cm
scale_bar_position_x <- 3.6
scale_bar_position_y <- 0.5
scale_bar <- data.frame("X" = c(scale_bar_position_x, scale_bar_position_x + scale_bar_length),
                        "Y" = c(scale_bar_position_y, scale_bar_position_y))

ling_spatial_map <- ggplot(data = ling_boundary_points, aes(x = X, y = -Y)) +
    geom_polygon(size= 0.5, color =  "black", linetype = "dashed", fill = "#e6e6e6", alpha = 0.5) +

    geom_point(data = ling_punches_points, shape = 21, size = 6, color = "black", aes(fill = ifelse(edge == 1, "edge", "center"))) +
    #geom_text(data = ling_punches_points, aes(label = Punch), size = 2) +
    geom_label_repel(data = ling_punches_points,
                     label.size = 0,
                     fill = NA,
                     nudge_x = 0,
                     nudge_y = 0.07,
                     size = 20,
                     aes(label = Punch)) +
    theme_void() +
    scale_fill_manual(values = state_colors) +
    labs(color = "") + #ggtitle("Ling et al. HCC Tumor") +
    coord_fixed(ratio = 1) +
    geom_line(data = scale_bar) +
    theme(legend.position = "none") +
    geom_text(data = data.frame("X" = c(scale_bar_position_x + scale_bar_length/2), "Y" = c(scale_bar_position_y - 0.1)), label= c(paste0(scale_bar_length, "cm", sep = " ")), size = 4)

ggsave(file= "../figures/ling_spatial_map.png", ling_spatial_map, width = 5, height = 5)


## POSTERIOR ESTIMATES

process_logs <- function(log_file) {
    print(log_file)

    log <- readLog(log_file)
    log_df <- as.data.frame(log) %>%
        dplyr::mutate(birthRateDiff = birthRateSVCanonical.loc1 - birthRateSVCanonical.loc0,
                      birthRateRatio = birthRateSVCanonical.loc1 / birthRateSVCanonical.loc0)
    
    return(log_df)
}

all_logs_df <- process_logs("../ling-application/logs/hcc-wes_unidir_state_comb.log")

all_logs_birthRate_df <- all_logs_df %>% 
    pivot_longer(., cols = c("birthRateSVCanonical.loc0", "birthRateSVCanonical.loc1"),
                 names_to = "state", 
                 names_prefix = "birthRateSVCanonical.", 
                 values_to = "birthRate")

#for main figure plot only state-dependent clock + unidirectional migration

#Tumor 1
ling_posteriors <- all_logs_birthRate_df %>% 
    ggplot(., aes(x=birthRate), color = "black") +
    geom_density(aes(fill = state), alpha=0.8) +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + scale_fill_manual(values=colors_loc) +
    theme(text=element_text(size=25))+
    xlab("Estimated birth rate") +
    theme(legend.position = "none") +ylab("")
ling_posteriors

ling_posteriors_violin <- all_logs_birthRate_df %>% 
    ggplot(., aes(x=ifelse(state == "loc0", "center", "edge"), y=birthRate), color = "black") +
    geom_violin(aes(fill = state), alpha=0.8) +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + scale_fill_manual(values=colors_loc) +
    theme(text=element_text(size=40))+
    xlab("Estimated birth rate") +
    theme(legend.position = "none") +ylab("") +
    stat_summary(fun=median, geom="point", size=1, color="black")
ling_posteriors_violin



ggsave(plot=ling_posteriors,
       file ="../figures/ling_posteriors.png", height = 5, width = 5)

ggsave(plot=ling_posteriors_violin,
       file ="../figures/ling_posteriors_voilin.png", height = 5, width = 3)

ling_ratio_posteriors <- all_logs_birthRate_df %>%
    ggplot(., aes(x=birthRateRatio), color = "black",  fill = "black") +
    geom_density(alpha=0.8, fill = "black") +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() +
    theme(text=element_text(size=25))+
    xlab("Estimated birth rate ratio \n(edge / center)") +
    theme(legend.position = "none") +ylab("") +
    geom_vline(xintercept = 1, linetype="dashed") 

ling_ratio_posteriors_violin <- all_logs_birthRate_df %>%
    ggplot(., aes(x="", y=birthRateRatio), color = "black",  fill = "black") +
    geom_violin(alpha=0.8, fill = "black") +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() +
    theme(text=element_text(size=30))+
    xlab("") +
    theme(legend.position = "none") +ylab("") +
    geom_hline(yintercept = 1, linetype="dashed")+
    stat_summary(fun=median, geom="point", size=1, color="lightgrey")

ling_ratio_posteriors
ling_ratio_posteriors_violin

ggsave(plot=ling_ratio_posteriors,
       file ="../figures/ling_ratio_posteriors.png", height = 5, width = 5)
ggsave(plot=ling_ratio_posteriors_violin,
       file ="../figures/ling_ratio_posteriors_violin.png", height = 5, width = 1.5)

### get stats

all_logs_birthRate_df %>% 
    summarize(mean(birthRateRatio))

#### MCC tree

#To generate MCC run combined_typed_node_mcc_trees.sh in li-application/out directory

mcc_file <- "../ling-application/trees/hcc-wes_unidir_state_comb.HCCtumor.typed.node.mcc.tree"

mcc_tree <- read.beast(file=mcc_file)


#mcc_tree@phylo$tip.label <- toupper(mcc_tree@phylo$tip.label)
    
anc_state_nodes <- mcc_tree@data %>% 
        #filter(node > length(T1_mcc_tree@phylo$tip.label)) %>% 
        mutate(type.prob = as.numeric(type.prob)) %>% 
        dplyr::mutate(edge = type.prob * as.integer(type == "loc1") + (1 - type.prob)*as.integer(type == "loc0")) %>% 
        dplyr::mutate(center = type.prob * as.integer(type == "loc0") + (1 - type.prob)*as.integer(type == "loc1")) %>% 
        dplyr::select(center, edge, type, node)
    
pies <- nodepie(anc_state_nodes, cols=1:2, alpha=0.8, color = edge_center_colors)
    
treeplot <- ggtree(mcc_tree, color = "darkgrey", size = 2) +
        geom_nodelab(aes(label = round(as.numeric(posterior),2)),  hjust="right", nudge_x = -0.03, nudge_y = 0.35, size = 20) +
        scale_color_manual(values = colors_loc) +
        theme(legend.position = "none") + 
        geom_tiplab(offset = 0.02, size = 20) +
        coord_cartesian(clip="off")
   
treeplot_pie <- ggtree::inset(treeplot, pies, width = 0.06, height = 0.06) + theme(legend.position = "none") #+ ggtitle(paste0(tumor_extract, migration_model, rep_extract, sep = " "))


treeplot_pie 
ggsave(plot=treeplot_pie,
           file="../figures/ling-wes-mcc.png", height = 10, width = 8)

       
       