#figureS6.R

#scripts to generate alternate states

library(tidyverse)
library(ggtree)
library(treeio)
library(beastio)
library(tumortree)
library(smoothr)
library(coda)
library(HDInterval)
library(ggrepel)

edge_center_colors <- tumortree::get_color_palette(names = c("edge", "center"))
colors_loc <- tumortree::get_color_palette(names = c("edge", "center"))
names(colors_loc) <- c("loc1", "loc0")

###### SPATIAL MAPS ######
# These are modified into 3D format in Affinity Designer
#Published edge/center labels from Li et al (Table S8)

t1_li_edge_labels <- data.frame("Punch" = c("t1z5", "t1l13", "t1z1", "t1f24", "t1z3", "t1f23", "t1f11", "t1l8", "t1f14", "t1f9", "t1l1", "t1l10", "t1l3", "t1f2", "t1l6", "t1f4"),
                                "edgeP" = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)) %>% 
    dplyr::mutate("Punch" = toupper(Punch))


t2_li_edge_labels <- data.frame("Punch" = c("t2z11", "t2z13", "t2z1", "t2f9", "t2f5", "t2f2", "t2f13", "t2z9", "t2z6"),
                                "edgeP" = c(1, 1, 1, 0, 0, 0, 0, 0, 0))  %>% 
    dplyr::mutate("Punch" = toupper(Punch))


#Spatial data formatted in process_li_wgs_data.R
T1_punch_coordinates <- read_csv(file = "../li-application/li_hcc_data/T1_punch_coordinates.csv") %>% 
    dplyr::left_join(.,t1_li_edge_labels, by="Punch") %>% 
    dplyr::mutate("punch_label" = gsub("T1", "", Punch))

#Boundary points
T1_boundary_coordinates <- read_csv(file = "../li-application/li_hcc_data/T1_boundary_coordinates.csv") 

T1_slice_order <- c("T1F", "T1H","T1K","T1L","T1O","T1Q","T1S","T1V","T1Y","T1Z","T1AB")
T1_punch_coordinates$Slice <- factor(T1_punch_coordinates$Slice, levels = T1_slice_order)

# Punch coordinates frm 
T1_wgs_punch_coordinates <- T1_punch_coordinates %>% 
    dplyr::filter(! is.na(edgeP))

T2_punch_coordinates <- read_csv(file = "../li-application/li_hcc_data/T2_punch_coordinates.csv") %>% 
    dplyr::left_join(.,t2_li_edge_labels, by="Punch")  %>% 
    dplyr::mutate("punch_label" = gsub("T2", "", Punch))

T2_boundary_coordinates <- read_csv(file = "../li-application/li_hcc_data/T2_boundary_coordinates.csv")

T2_wgs_punch_coordinates <- T2_punch_coordinates %>% 
    dplyr::filter(! is.na(edgeP))

#### ALTERNATE STATES (10% to boundary)

t1_map_plot_alt_states <- T1_punch_coordinates %>% 
    
    #filter(! is.na(edgeP)) %>% 
    ggplot(., aes(X,-Y)) + geom_polygon(data = T1_boundary_coordinates,
                                        fill = "lightgrey",
                                        alpha = 0.5,
                                        color = "black",
                                        linetype = "dashed",
                                        aes(group = Slice)) +
    geom_point(size = 1, color = "darkgrey") +
    geom_point(data = T1_wgs_punch_coordinates, size = 4, shape = 21, color = "black",
               aes(fill = ifelse(edge, "edge", "center"))) +
    #scale_fill_manual(values = sim_colors ) +
    # labs(fill = "") +
    #facet_wrap(~factor(Slice, levels = T1_slices), scales = "free") +
    theme_void() + theme(legend.position = "none") +
    scale_fill_manual(values = edge_center_colors) +
    geom_label_repel(data=T1_wgs_punch_coordinates,
                     label.size = 0,
                     fill = NA,
                     nudge_x = 0,
                     nudge_y = 2,
                     size = 5,
                     aes(label = punch_label))
t1_map_plot_alt_states 
ggsave(plot=t1_map_plot_alt_states,
       file ="../figures/t1_li_wgs_punch_map_alt_states.png", height = 5, width = 5*t1_coord_ratio)
t2_map_plot_alt_states <- T2_punch_coordinates %>% 
    
    #filter(! is.na(edgeP)) %>% 
    ggplot(., aes(X,-Y)) + geom_polygon(data = T2_boundary_coordinates,
                                        fill = "lightgrey",
                                        alpha = 0.5,
                                        color = "black",
                                        linetype = "dashed",
                                        aes(group = Slice)) +
    geom_point(size = 1, color = "darkgrey") +
    geom_point(data = T2_wgs_punch_coordinates, size = 4, shape = 21, color = "black",
               aes(fill = ifelse(edge, "edge", "center"))) +
    theme_void() + theme(legend.position = "none") +
    scale_fill_manual(values = edge_center_colors) +
    geom_label_repel(data=T2_wgs_punch_coordinates,
                     label.size = 0,
                     fill = NA,
                     nudge_x = 0,
                     nudge_y = 2,
                     size = 5,
                     aes(label = punch_label))
t2_map_plot_alt_states
ggsave(plot=t2_map_plot_alt_states,
       file ="../figures/t2_li_wgs_punch_map_alt_states.png", height = 5, width = 5*t2_coord_ratio)
##### PLOT POSTERIORS ####

# Get summary from sdevo_li_estimates_summary.R
all_logs_birthRate_df <- read_tsv("../li-application/stats/sdevo_estimates_summary.tsv")
# Also for supplement -- compare published states versus 10% diameter cutoff
## Tumor 1
t1_wgs_posteriors_plot_alt_states <- all_logs_birthRate_df %>% 
    dplyr::filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           states == "newstates",
           tumor == "T1") %>% 
    mutate("category" = paste0(state, subset)) %>% 
    ggplot(., aes(x=birthRate, group=category), color = "black") +
    geom_density(aes(fill = state, group = category), alpha=0.5) +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + scale_fill_manual(values=colors_loc) +
    theme(text=element_text(size=20))+
    xlab("Estimated birth rate") +
    theme(legend.position = "none") +ylab("") 

t1_wgs_posteriors_plot_alt_states
ggsave(plot=t1_wgs_posteriors_plot_alt_states ,file ="../figures/t1_li_wgs_posteriors_newstates.png", height = 5, width = 5)

# Tumor 2
t2_wgs_posteriors_plot_alt_states <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           states == "newstates",
           tumor == "T2") %>% 
    mutate("category" = paste0(state, subset)) %>% 
    ggplot(., aes(x=birthRate), color = "black") +
    geom_density(aes(fill = state, group = category), alpha=0.8) +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + scale_fill_manual(values=colors_loc) +
    theme(text=element_text(size=20))+
    xlab("Estimated birth rate") +
    theme(legend.position = "none") +ylab("") 
t2_wgs_posteriors_plot_alt_states
ggsave(plot=t2_wgs_posteriors_plot_alt_states ,file ="../figures/t2_li_wgs_posteriors_newstates.png", height = 5, width = 5)

#Ratios posteriors
## Tumor 1
t1_wgs_ratio_posteriors_plot_newstates_state_clock <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           states == "newstates",
           tumor == "T1") %>% 
    ggplot(., aes(x=birthRateRatio), color = "black",  fill = "black") +
    geom_density(alpha=0.6, fill = "black", aes(group = subset)) +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + 
    theme(text=element_text(size=20))+
    xlab("Estimated birth rate ratio (edge / center)") +
    theme(legend.position = "none") +ylab("") +
    geom_vline(xintercept = 1, linetype="dashed")

t1_wgs_ratio_posteriors_plot_newstates_state_clock
ggsave(plot=t1_wgs_ratio_posteriors_plot_newstates_state_clock,
       file ="../figures/t1_li_wgs_ratio_posteriors_newstates_stateclock.png", height = 5, width = 5)

# Tumor 2
t2_wgs_ratio_posteriors_plot_newstates_state_clock <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           states == "newstates",
           tumor == "T2") %>% 
    ggplot(., aes(x=birthRateRatio), color = "black",  fill = "black") +
    geom_density(alpha=0.6, fill = "black", aes(group = subset)) +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + 
    theme(text=element_text(size=20))+
    xlab("Estimated birth rate ratio (edge / center)") +
    theme(legend.position = "none") +ylab("") +
    geom_vline(xintercept = 1, linetype="dashed")

t2_wgs_ratio_posteriors_plot_newstates_state_clock
ggsave(plot=t2_wgs_ratio_posteriors_plot_newstates_state_clock,
       file ="../figures/t2_li_wgs_ratio_posteriors_newstates_stateclock.png", height = 5, width = 5)

####### MCC TREES ###############

####### MCC TREES ###############

#To generate MCC run combined_typed_node_mcc_trees.sh in li-application/out directory

mcc_tree_files <- list.files(path = "../li-application/logs",
                             pattern="T[1-2]_wgs_newstates_unidir_state_comb.HCCtumor_mcc.tree", 
                             full.names = TRUE)

#only do unidirectional

for (mcc_file in mcc_tree_files) {
    mcc_tree <- read.beast(file=mcc_file)
    migration_model <- ifelse(grepl("unidir", basename(mcc_file)), "unidirectional", "bidirectional")
    tumor_extract <- regmatches(basename(mcc_file), gregexpr("T[12]", basename(mcc_file), perl = TRUE))[[1]]
    states_extract <- ifelse(grepl("newstates", basename(mcc_file)), "newstates", "oldstates")
    clock_extract <- ifelse(grepl("strict", basename(mcc_file)), "strict", "state-dependent")
    
    mcc_tree@phylo$tip.label <- toupper(gsub("t[12]", "", mcc_tree@phylo$tip.label))
    
    anc_state_nodes <- mcc_tree@data %>% 
        #filter(node > length(T1_mcc_tree@phylo$tip.label)) %>% 
        mutate(type.prob = as.numeric(type.prob)) %>% 
        dplyr::mutate(edge = type.prob * as.integer(type == "loc1") + (1 - type.prob)*as.integer(type == "loc0")) %>% 
        dplyr::mutate(center = type.prob * as.integer(type == "loc0") + (1 - type.prob)*as.integer(type == "loc1")) %>% 
        dplyr::select(center, edge, type, node)
    
    pies <- nodepie(anc_state_nodes, cols=1:2, alpha=0.8, color = edge_center_colors)
    
    treeplot <- ggtree(mcc_tree, color = "darkgrey", size = 2) +
        #geom_point(size = 2) +
        scale_color_manual(values = colors_loc) +
        theme(legend.position = "none") + 
        geom_tiplab(offset = 0.05, size = 5) +
        coord_cartesian(clip="off") 
    
    treeplot_pie <- ggtree::inset(treeplot, pies, width = 0.05, height = 0.05) + theme(legend.position = "none") +
        geom_nodelab(aes(label = ifelse(as.numeric(posterior) <= 1, round(as.numeric(posterior), 2), "")),
                     nudge_x = -0.06, nudge_y = 0.2, size  = 5, hjust='right')
    
    fig_file <- paste0("../figures/", tumor_extract,  "_",
                       migration_model, "_",
                       clock_extract,  "_",
                       states_extract, "_",
                       
                       "mcc_tree.png")
    print(fig_file)
    ggsave(plot=treeplot_pie,
           file=fig_file, height = 10, width = 8)
}

#For legend

all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           states == "newstates") %>% 
    group_by(tumor) %>% 
    summarise(mean(birthRateRatio))

