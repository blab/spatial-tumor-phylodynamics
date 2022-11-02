# figure5.R

## Script to generate subfigures for Figure 5

###### SETUP ######
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

###### PUBLISHED STATES #######

t1_map_plot <- T1_punch_coordinates %>% 
    
    #filter(! is.na(edgeP)) %>% 
    ggplot(., aes(X,-Y)) + geom_polygon(data = T1_boundary_coordinates,
                                        fill = "lightgrey",
                                        alpha = 0.5,
                                        color = "black",
                                        linetype = "dashed",
                                        aes(group = Slice)) +
    #geom_point(size = 1, color = "darkgrey") +
    geom_point(data = T1_wgs_punch_coordinates, size = 4,
               shape = 21, color = "black", 
               aes(fill = ifelse(edgeP == 1, "edge", "center"))) +
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
t1_map_plot
t1_coord_ratio <- (max(T1_punch_coordinates$X) - min(T1_punch_coordinates$X)) / (max(T1_punch_coordinates$Y) - min(T1_punch_coordinates$Y))
ggsave(plot=t1_map_plot,file ="../figures/t1_li_wgs_punch_map_published_states.png", height = 5, width = 5*t1_coord_ratio )

t2_map_plot <- T2_punch_coordinates %>% 
    
    #filter(! is.na(edgeP)) %>% 
    ggplot(., aes(X,-Y)) + geom_polygon(data = T2_boundary_coordinates,
                                        fill = "lightgrey",
                                        alpha = 0.5,
                                        color = "black",
                                        linetype = "dashed",
                                        aes(group = Slice)) +
    #geom_point(size = 1, color = "darkgrey") +
    geom_point(data = T2_wgs_punch_coordinates, size = 4, shape = 21, color = "black",  aes(fill = ifelse(edgeP == 1, "edge", "center"))) +
    theme_void() + theme(legend.position = "none") +
    scale_fill_manual(values = edge_center_colors) +
    geom_label_repel(data=T2_wgs_punch_coordinates,
                     label.size = 0,
                     fill = NA,
                     nudge_x = 0,
                     nudge_y = 2,
                     size = 5,
                     aes(label = punch_label))
t2_map_plot

t2_coord_ratio <- (max(T2_punch_coordinates$X) - min(T2_punch_coordinates$X)) / (max(T2_punch_coordinates$Y) - min(T2_punch_coordinates$Y))
ggsave(plot=t2_map_plot,file ="../figures/t2_li_wgs_punch_map_published_states.png", height = 5, width = 5*t2_coord_ratio)

##### PLOT POSTERIORS ####

#Record of local directory
# log_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/curr_runs/out",
#                         pattern=".log", 
#                         full.names = TRUE)

log_files <- list.files(path = "../li-application/out",
                        pattern=".log",
                        full.names = TRUE)


#temporarily get rid of last corrupted log file
log_files <- log_files[! grepl("T1_wgs_newstates_bidir_state_rep2.log", log_files)]

#logs <- purrr::map(log_files, readLog)

process_logs <- function(log_file) {
    print(log_file)
    rep_extract <- regmatches(basename(log_file),
                              gregexpr("(?<=rep)[0-9]", basename(log_file), perl = TRUE))[[1]]
    migration_model <- ifelse(grepl("unidir", basename(log_file)), "unidirectional", "bidirectional")
    tumor_extract <- regmatches(basename(log_file), gregexpr("T[12]", basename(log_file), perl = TRUE))[[1]]
    states_extract <- ifelse(grepl("newstates", basename(log_file)), "newstates", "oldstates")
    clock_extract <- ifelse(grepl("strict", basename(log_file)), "strict", "state-dependent")
    log <- readLog(log_file)
    log_df <- as.data.frame(log) %>%
        add_column("migration_model" = migration_model,
                   "rep" = rep_extract,
                   "tumor" = tumor_extract,
                   "states" = states_extract,
                   "clock_model" = clock_extract) %>%
        dplyr::mutate(birthRateDiff = birthRateSVCanonical.loc1 - birthRateSVCanonical.loc0,
                      birthRateRatio = birthRateSVCanonical.loc1 / birthRateSVCanonical.loc0)
    
    return(log_df)
}

all_logs_df <- purrr::map(log_files, process_logs) %>% 
    bind_rows

all_logs_birthRate_df <- all_logs_df %>% 
    pivot_longer(., cols = c("birthRateSVCanonical.loc0", "birthRateSVCanonical.loc1"),
                 names_to = "state", 
                 names_prefix = "birthRateSVCanonical.", 
                 values_to = "birthRate")

#for main figure plot only state-dependent clock + unidirectional migration

#Tumor 1
t1_wgs_posteriors_plot_oldstates_state_clock <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           states == "oldstates",
           tumor == "T1") %>% 
    ggplot(., aes(x=birthRate), color = "black") +
    geom_density(aes(fill = state), alpha=0.8) +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + scale_fill_manual(values=colors_loc) +
    theme(text=element_text(size=25))+
    xlab("Estimated birth rate") +
    theme(legend.position = "none") +ylab("")

t1_wgs_violin_posteriors_plot_oldstates_state_clock <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           states == "oldstates",
           tumor == "T1") %>% 
    ggplot(., aes(x = ifelse(state == "loc0", "center", "edge"), y=birthRate), color = "black") +
    geom_violin(aes(fill = state), alpha=0.8) +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + scale_fill_manual(values=colors_loc) +
    theme(text=element_text(size=40))+
    xlab("Estimated birth rate") +
    theme(legend.position = "none") +ylab("") +
    stat_summary(fun=median, geom="point", size=1, color="black")

t1_wgs_violin_posteriors_plot_oldstates_state_clock
t1_wgs_posteriors_plot_oldstates_state_clock 

ggsave(plot=t1_wgs_posteriors_plot_oldstates_state_clock,
       file ="../figures/t1_li_wgs_posteriors_oldstates_stateclock.png", height = 5, width = 5)

ggsave(plot=t1_wgs_violin_posteriors_plot_oldstates_state_clock,
       file ="../figures/t1_li_wgs_posteriors_oldstates_stateclock_violin.png", height = 5, width = 3)

#Tumor 2
t2_wgs_posteriors_plot_oldstates_state_clock <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           states == "oldstates",
           tumor == "T2") %>% 
    ggplot(., aes(x=birthRate), color = "black") +
    geom_density(aes(fill = state), alpha=0.8) +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + scale_fill_manual(values=colors_loc) +
    theme(text=element_text(size=25))+
    xlab("Estimated birth rate") +
    theme(legend.position = "none") +ylab("") 

t2_wgs_violin_posteriors_plot_oldstates_state_clock <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           states == "oldstates",
           tumor == "T2") %>% 
    ggplot(., aes(x = ifelse(state == "loc0", "center", "edge"), y=birthRate), color = "black") +
    geom_violin(aes(fill = state), alpha=0.8) +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + scale_fill_manual(values=colors_loc) +
    theme(text=element_text(size=40))+
    xlab("Estimated birth rate") +
    theme(legend.position = "none") +ylab("") +
    stat_summary(fun=median, geom="point", size=1, color="black")

t2_wgs_posteriors_plot_oldstates_state_clock  

ggsave(plot=t2_wgs_posteriors_plot_oldstates_state_clock ,
       file ="../figures/t2_li_wgs_posteriors_oldstates_stateclock.png", height = 5, width = 3)

ggsave(plot=t2_wgs_violin_posteriors_plot_oldstates_state_clock ,
       file ="../figures/t2_li_wgs_posteriors_oldstates_stateclock_violin.png", height = 5, width = 3)

#also plot estimated differences
## Tumor 1
t1_wgs_ratio_posteriors_plot_oldstates_state_clock <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           states == "oldstates",
           tumor == "T1") %>% 
    ggplot(., aes(x=birthRateRatio), color = "black",  fill = "black") +
    geom_density(alpha=0.8, fill = "black") +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + 
    theme(text=element_text(size=30))+
    xlab("Estimated birth rate ratio (edge / center)") +
    theme(legend.position = "none") +ylab("") +
    geom_vline(xintercept = 1, linetype="dashed")

t1_wgs_ratio_posteriors_plot_oldstates_state_clock_violin <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           states == "oldstates",
           tumor == "T1") %>% 
    ggplot(., aes(x="", y=birthRateRatio), color = "black",  fill = "black") +
    geom_violin(alpha=0.8, fill = "black") +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + 
    theme(text=element_text(size=40))+
    xlab("Estimated birth rate ratio (edge / center)") +
    theme(legend.position = "none") +ylab("") +
    geom_hline(yintercept = 1, linetype="dashed") +
    stat_summary(fun=median, geom="point", size=1, color="lightgrey")

t1_wgs_ratio_posteriors_plot_oldstates_state_clock 
t1_wgs_ratio_posteriors_plot_oldstates_state_clock_violin
ggsave(plot=t1_wgs_ratio_posteriors_plot_oldstates_state_clock ,
       file ="../figures/t1_li_wgs_ratio_posteriors_oldstates_stateclock.png", height = 5, width = 5)

ggsave(plot=t1_wgs_ratio_posteriors_plot_oldstates_state_clock_violin ,
       file ="../figures/t1_li_wgs_ratio_posteriors_oldstates_stateclock_violin.png", height = 5, width = 1.5)

## Tumor 2
t2_wgs_ratio_posteriors_plot_oldstates_state_clock <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           states == "oldstates",
           tumor == "T2") %>% 
    ggplot(., aes(x=birthRateRatio), color = "black",  fill = "black") +
    geom_density(aes(), alpha=0.8, fill = "black") +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + 
    theme(text=element_text(size=30))+
    xlab("Estimated birth rate ratio (edge / center)") +
    theme(legend.position = "none") +ylab("") +
    geom_vline(xintercept = 1, linetype="dashed")


t2_wgs_ratio_posteriors_plot_oldstates_state_clock_violin <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           states == "oldstates",
           tumor == "T2") %>% 
    ggplot(., aes(x="", y=birthRateRatio), color = "black",  fill = "black") +
    geom_violin(alpha=0.8, fill = "black") +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + 
    theme(text=element_text(size=40))+
    xlab("Estimated birth rate ratio (edge / center)") +
    theme(legend.position = "none") +ylab("") +
    geom_hline(yintercept = 1, linetype="dashed") +
    stat_summary(fun=median, geom="point", size=1, color="lightgrey")

t2_wgs_ratio_posteriors_plot_oldstates_state_clock 

ggsave(plot=t2_wgs_ratio_posteriors_plot_oldstates_state_clock ,
       file ="../figures/t2_li_wgs_ratio_posteriors_oldstates_stateclock.png", height = 5, width = 5)
ggsave(plot=t2_wgs_ratio_posteriors_plot_oldstates_state_clock_violin ,
       file ="../figures/t2_li_wgs_ratio_posteriors_oldstates_stateclock_violin.png", height = 5, width = 1.5)

#get HDI intervals to put in text
t1_estimate_summaries <- all_logs_birthRate_df %>% 
    dplyr::filter(tumor == "T1", 
                  migration_model == "unidirectional",
                  clock_model == "state-dependent",
                  states == "oldstates") %>% 
    group_by(state) %>% 
    summarise(birthRate_HDI90_lower = hdi(birthRate, credMass = 0.9)[1], 
              birthRate_HDI90_upper = hdi(birthRate, credMass = 0.9)[2],
              birthRate_HDI50_lower = hdi(birthRate, credMass = 0.5)[1], 
              birthRate_HDI50_upper = hdi(birthRate, credMass = 0.5)[2],
              birthRateRatio_HDI90_lower = hdi(birthRateRatio, credMass = 0.9)[1], 
              birthRateRatio_HDI90_upper = hdi(birthRateRatio, credMass = 0.9)[2],
              birthRateRatio_HDI50_lower = hdi(birthRateRatio, credMass = 0.5)[1], 
              birthRateRatio_HDI50_upper = hdi(birthRateRatio, credMass = 0.5)[2],
              birthRateRatio_mean = mean(birthRateRatio))

t1_estimate_summaries %>% 
    select(contains("birthRateRatio"))
write_csv(t1_estimate_summaries, file = "../li-application/stats/t1_wgs_birth_rate_estimate_summary.csv")

#get HDI intervals to put in text
t2_estimate_summaries <- all_logs_birthRate_df %>% 
    dplyr::filter(tumor == "T2", 
                  migration_model == "unidirectional",
                  clock_model == "state-dependent",
                  states == "oldstates") %>% 
    group_by(state) %>% 
    summarise(birthRate_HDI90_lower = hdi(birthRate, credMass = 0.9)[1], 
              birthRate_HDI90_upper = hdi(birthRate, credMass = 0.9)[2],
              birthRate_HDI50_lower = hdi(birthRate, credMass = 0.5)[1], 
              birthRate_HDI50_upper = hdi(birthRate, credMass = 0.5)[2],
              birthRateRatio_HDI90_lower = hdi(birthRateRatio, credMass = 0.9)[1], 
              birthRateRatio_HDI90_upper = hdi(birthRateRatio, credMass = 0.9)[2],
              birthRateRatio_HDI50_lower = hdi(birthRateRatio, credMass = 0.5)[1], 
              birthRateRatio_HDI50_upper = hdi(birthRateRatio, credMass = 0.5)[2],
              birthRateRatio_mean = mean(birthRateRatio))

t2_estimate_summaries %>% 
    select(contains("birthRateRatio"))
write_csv(t2_estimate_summaries, file = "../li-application/stats/t2_wgs_birth_rate_estimate_summary.csv")

####### MCC TREES ###############

#To generate MCC run combined_typed_node_mcc_trees.sh in li-application/out directory

mcc_tree_files <- list.files(path = "../li-application/out",
                             pattern="oristates_unidir_state_comb.HCCtumor_mcc.tree", 
                             full.names = TRUE)

#only do unidirectional
# For testing
#mcc_file <- "../li-application/out/T2_wgs_oristates_unidir_state_comb.HCCtumor_mcc.tree"
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
        geom_tiplab(offset = 0.05, size = 20) +
        coord_cartesian(clip="off")
    
    treeplot_pie <- ggtree::inset(treeplot, pies, width = 0.06, height = 0.06) + theme(legend.position = "none") #+ ggtitle(paste0(tumor_extract, migration_model, rep_extract, sep = " "))
    treeplot_pie <- treeplot_pie + geom_nodelab(aes(label = ifelse(as.numeric(posterior) <= 1, round(as.numeric(posterior), 2), "")),
                                                nudge_x = -0.05, nudge_y = 0.2, size  = 20, hjust='right')
    #treeplot_pie
    fig_file <- paste0("../figures/", tumor_extract,  "_",
                       migration_model, "_",
                       clock_extract,  "_",
                       states_extract, "_",
                       
                       "mcc_tree.png")
    print(fig_file)
    ggsave(plot=treeplot_pie,
           file=fig_file, height = 10, width = 8)
}


