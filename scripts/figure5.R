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
library(rjson)


edge_center_colors <- tumortree::get_color_palette(names = c("edge", "center"))
colors_loc <- tumortree::get_color_palette(names = c("edge", "center"))
names(colors_loc) <- c("loc1", "loc0")

###### SPATIAL MAPS ######

#Published edge/center labels from Li et al (Table S8)
t1_li_edge_labels <- data.frame("punch" = c("t1z5", "t1l13", "t1z1", "t1f24", "t1z3", "t1f23", "t1f11", "t1l8", "t1f14", "t1f9", "t1l1", "t1l10", "t1l3", "t1f2", "t1l6", "t1f4"),
                                "edgeP" = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)) %>% 
    dplyr::mutate("Punch" = toupper(punch))


t2_li_edge_labels <- data.frame("punch" = c("t2z11", "t2z13", "t2z1", "t2f9", "t2f5", "t2f2", "t2f13", "t2z9", "t2z6"),
                                "edgeP" = c(1, 1, 1, 0, 0, 0, 0, 0, 0))  %>% 
    dplyr::mutate("Punch" = toupper(punch))


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

###### PLOT MAXIMUM LIKELIHOOD TREES ######
#Fasta alignment files compiled in process_li_wgs_data.R and write_ling_state_clocks.R
## Generate MLE trees with run_nextstrain_divergence_tree.sh
## Here these will only be used to convert divergence into number of genetic mutations
## VAF cutoff is 0.05 for variants
fasta_file_t1 <- "../li-application/nextstrain_analysis/data/li_t1_wgs.fa"
fasta_file_t2 <- "../li-application/nextstrain_analysis/data/li_t2_wgs.fa"

##read in data into phyDat format (phangorn/ape)
t1_data <- phangorn::read.phyDat(file = fasta_file_t1, format = "fasta", type = "DNA")
t2_data <- phangorn::read.phyDat(file = fasta_file_t2, format = "fasta", type = "DNA")
#t3_data <- phangorn::read.phyDat(file = fasta_file_t3, format = "fasta", type = "DNA")

t1_seq_length <- sum(attr(t1_data, "weight"))
t2_seq_length <- sum(attr(t2_data, "weight"))
#t3_seq_length <- sum(attr(t3_data, "weight"))

######## MAXIMUM LIKELIHOOD TREE FROM FASTTREE + AUGUR PIPELINE #################

#Read in trees
max_likelihood_t1_wgs_tree <- ape::read.tree(file = "../li-application/nextstrain_analysis/li_t1_tree.nwk")
max_likelihood_t2_wgs_tree <- ape::read.tree(file = "../li-application/nextstrain_analysis/li_t2_tree.nwk")

#read in branch info
t1_branch_info <- fromJSON(file = "../li-application/nextstrain_analysis/t1_branch_lengths.json")
t2_branch_info <- fromJSON(file = "../li-application/nextstrain_analysis/t2_branch_lengths.json")

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

t1_treePars_obj@data$node_label <- c(t1_treePars_obj@phylo$tip.label, t1_treePars_obj@phylo$node.label)
t2_treePars_obj@data$node_label <- c(t2_treePars_obj@phylo$tip.label, t2_treePars_obj@phylo$node.label)

t1_node_confidence <- unlist(t1_branch_info$nodes)[grepl(".confidence", names(unlist(t1_branch_info$nodes)))]
names(t1_node_confidence) <- gsub(".confidence", "", names(t1_node_confidence))

t2_node_confidence <- unlist(t2_branch_info$nodes)[grepl(".confidence", names(unlist(t2_branch_info$nodes)))]
names(t2_node_confidence) <- gsub(".confidence", "", names(t2_node_confidence))

t1_treePars_obj@data$confidence <- NA
t2_treePars_obj@data$confidence <- NA

t1_treePars_obj@data$confidence[match(names(t1_node_confidence), t1_treePars_obj@data$node_label)] <- round(t1_node_confidence,2)
t2_treePars_obj@data$confidence[match(names(t2_node_confidence), t2_treePars_obj@data$node_label)] <- round(t2_node_confidence,2)

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
              vjust = "bottom", hjust = "midde", color = "black") +
    geom_nodelab(aes(label = ifelse(confidence < 1, confidence, "")),
                 hjust = "right",
                 vjust = "bottom", 
                 size = 3.5,
                 nudge_x = -1000, nudge_y = 0.2,
                 color = "black")

x_max_limit <- layer_scales(t1_wgs_tree)$x$get_limits()
new_x_limit <- 1.1*x_max_limit[2]
t1_wgs_tree  <- t1_wgs_tree + xlim(0, new_x_limit)
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
              vjust = "bottom", hjust = "midde", color = "black") +
    geom_nodelab(aes(label = ifelse(confidence < 1, confidence, "")),
                 hjust = "right",
                 vjust = "bottom", 
                 color = "black",
                 size = 3.5,
                 nudge_x = -1000, nudge_y = 0.2)
x_max_limit <- layer_scales(t2_wgs_tree)$x$get_limits()
new_x_limit <- 1.1*x_max_limit[2]
t2_wgs_tree  <- t2_wgs_tree + xlim(0, new_x_limit)

t2_wgs_tree 
ggsave(filename = "../figures/t2_wgs_genetic_state_tree.png", t2_wgs_tree , height = 4, width = 3)
#Get sscale
#layer_scales(t2_wgs_tree)$x$range$range
#layer_scales(t2_wgs_tree)$y$range$range

##### PLOT POSTERIORS ####

#Record of local directory
# log_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/curr_runs/out",
#                         pattern=".log", 
#                         full.names = TRUE)

li_log_files <- list.files(path = "../li-application/out",
                        pattern=".log",
                        full.names = TRUE)

#get reid of inidivudal chains
li_log_files <- li_log_files[! grepl("chain1", li_log_files)]
li_log_files <- li_log_files[! grepl("T1red", li_log_files)]
li_log_files <- li_log_files[! grepl("bidir", li_log_files)]

log_files <- c(li_log_files)

#logs <- purrr::map(log_files, readLog)

process_logs <- function(log_file) {
    print(log_file)
    rep_extract <- regmatches(basename(log_file),
                              gregexpr("(?<=rep)[0-9]", basename(log_file), perl = TRUE))[[1]]
    subset_extract <- regmatches(basename(log_file),
                              gregexpr("(?<=unidir_)[0-9]", basename(log_file), perl = TRUE))[[1]]
    if (length(subset_extract) == 0) {
        subset_extract <- "1"
    }
    
    migration_model <- ifelse(grepl("unidir", basename(log_file)), "unidirectional", "bidirectional")
    tumor_extract <- regmatches(basename(log_file), gregexpr("T[12]", basename(log_file), perl = TRUE))[[1]]
    if (length(tumor_extract) == 0) {
        tumor_extract <- "T3"
    }
    states_extract <- ifelse(grepl("newstates", basename(log_file)), "newstates", "oldstates")
    clock_extract <- ifelse(grepl("strict", basename(log_file)), "strict", "state-dependent")
    log <- readLog(log_file, burnin = 0.2)
    ess <- coda::effectiveSize(log)
    minESSbirth <- min(ess["birthRateSVCanonical.loc1"], ess["birthRateSVCanonical.loc0"])
    log_df <- as.data.frame(log) %>%
        add_column("migration_model" = migration_model,
                   "tumor" = tumor_extract,
                   "rep" = rep_extract,
                   "subset" = subset_extract, 
                   "states" = states_extract,
                   "clock_model" = clock_extract, 
                   "minESS" = minESSbirth) %>%
        dplyr::mutate(birthRateDiff = birthRateSVCanonical.loc1 - birthRateSVCanonical.loc0,
                      birthRateRatio = birthRateSVCanonical.loc1 / birthRateSVCanonical.loc0)
    
    return(log_df)
}

all_logs_df <- purrr::map(log_files, process_logs) %>% 
    bind_rows

all_logs_birthRate_df <- all_logs_df %>% 
#    dplyr::filter(minESS > 50) %>% 
    pivot_longer(., cols = c("birthRateSVCanonical.loc0", "birthRateSVCanonical.loc1"),
                 names_to = "state", 
                 names_prefix = "birthRateSVCanonical.", 
                 values_to = "birthRate")

#for main figure plot only state-dependent clock + unidirectional migration

#Tumor 1
t1_wgs_violin_posteriors_plot_oldstates_state_clock_summary <- all_logs_birthRate_df %>% 
    dplyr::filter(migration_model == "unidirectional",
                  clock_model == "state-dependent",
                  states == "oldstates",
                  tumor == "T1") %>% 
    group_by(subset, state) %>% 
    summarize(birthRate_hdi95_lower = hdi(birthRate, credMass = 0.95)[1],
              birthRate_hdi95_upper = hdi(birthRate, credMass = 0.95)[2],
              birthRate_mean = mean(birthRate))

t1_wgs_violin_posteriors_plot_oldstates_state_clock <- all_logs_birthRate_df %>% 
    dplyr::filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           states == "oldstates",
           tumor == "T1") %>% 
    #ggplot(., aes(x = ifelse(state == "loc0", "center", "edge"), y=birthRate), color = "black") +
    ggplot(., aes(x = subset, y=birthRate, fill=state), color = "black") +
    geom_violin(aes(fill = state), alpha=0.8) +
    geom_boxplot(width=.1, outlier.alpha = 0) +
    #geom_point(data=t1_wgs_violin_posteriors_plot_oldstates_state_clock_summary, aes(x= subset, y=birthRate_mean), color = "black") +
    #geom_errorbar(data=t1_wgs_violin_posteriors_plot_oldstates_state_clock_summary,
    #              aes(x= subset, ymin=birthRate_hdi95_lower, y=birthRate_mean, ymax=birthRate_hdi95_upper), color = "black", width =0) +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + scale_fill_manual(values=colors_loc) +
    theme(text=element_text(size=15))+
    ylab("Estimated birth rate") +
    theme(legend.position = "none") +xlab("") +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x = element_blank()) +
    facet_wrap(~state) +
    theme( strip.background = element_blank(),
           strip.text = element_blank())
    
    #stat_summary(aes(group =state) fun=median, geom="point", size=1, color="black")

t1_wgs_violin_posteriors_plot_oldstates_state_clock
#t1_wgs_posteriors_plot_oldstates_state_clock 

# ggsave(plot=t1_wgs_posteriors_plot_oldstates_state_clock,
#        file ="../figures/t1_li_wgs_posteriors_oldstates_stateclock.png", height = 5, width = 5)

ggsave(plot=t1_wgs_violin_posteriors_plot_oldstates_state_clock,
       file ="../figures/t1_li_wgs_posteriors_oldstates_stateclock_violin.png", height = 5, width = 3)

#Tumor 2

t2_wgs_violin_posteriors_plot_oldstates_state_clock_summary <- all_logs_birthRate_df %>% 
    dplyr::filter(migration_model == "unidirectional",
                  clock_model == "state-dependent",
                  states == "oldstates",
                  tumor == "T2") %>% 
    group_by(subset, state) %>% 
    summarize(birthRate_hdi95_lower = hdi(birthRate)[1],
              birthRate_hdi95_upper = hdi(birthRate)[2],
              birthRate_mean = mean(birthRate))

t2_wgs_violin_posteriors_plot_oldstates_state_clock <- all_logs_birthRate_df %>% 
    dplyr::filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           states == "oldstates",
           tumor == "T2") %>% 
    #ggplot(., aes(x = ifelse(state == "loc0", "center", "edge"), y=birthRate), color = "black") +
    ggplot(., aes(x = subset, y=birthRate, fill=state), color = "black") +
    geom_violin(aes(fill = state), alpha=0.8) +
    geom_boxplot(width=.1, outlier.alpha = 0) +
    # geom_point(data=t2_wgs_violin_posteriors_plot_oldstates_state_clock_summary, aes(x= subset, y=birthRate_mean), color = "black") +
    # geom_errorbar(data=t2_wgs_violin_posteriors_plot_oldstates_state_clock_summary,
    #               aes(x= subset, ymin=birthRate_hdi95_lower, y=birthRate_mean, ymax=birthRate_hdi95_upper), color = "black", width =0) +
    # 
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + scale_fill_manual(values=colors_loc) +
    theme(text=element_text(size=15))+
    ylab("Estimated birth rate") +
    theme(legend.position = "none") +xlab("") +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x = element_blank()) +
    facet_wrap(~state) +
    theme( strip.background = element_blank(),
           strip.text = element_blank())
#t2_wgs_posteriors_plot_oldstates_state_clock  
t2_wgs_violin_posteriors_plot_oldstates_state_clock


ggsave(plot=t2_wgs_violin_posteriors_plot_oldstates_state_clock ,
       file ="../figures/t2_li_wgs_posteriors_oldstates_stateclock_violin.png", height = 5, width = 3)


########## Plot edge/center birth ratios ###############
## Tumor 1

t1_wgs_ratio_posteriors_plot_oldstates_state_clock_violin_summary <- all_logs_birthRate_df %>% 
    dplyr::filter(migration_model == "unidirectional",
                  clock_model == "state-dependent",
                  states == "oldstates",
                  tumor == "T1") %>% 
    group_by(subset, state) %>% 
    summarize(birthRateRatio_hdi95_lower = hdi(birthRateRatio)[1],
              birthRateRatio_hdi95_upper = hdi(birthRateRatio)[2],
              birthRateRatio_mean = mean(birthRateRatio))

t1_wgs_ratio_posteriors_plot_oldstates_state_clock_violin <- all_logs_birthRate_df %>% 
    dplyr::filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           states == "oldstates",
           tumor == "T1") %>% 
    ggplot(., aes(x=subset, y=birthRateRatio), color = "black",  fill = "black") +
    geom_violin(alpha=0.8, fill = "black") +
    geom_boxplot(width=.1, outlier.alpha = 0, color = "grey") +
    # geom_point(data=t1_wgs_ratio_posteriors_plot_oldstates_state_clock_violin_summary, aes(x= subset, y=birthRateRatio_mean), color = "grey") +
    # geom_errorbar(data=t1_wgs_ratio_posteriors_plot_oldstates_state_clock_violin_summary,
    #               aes(x= subset, ymin=birthRateRatio_hdi95_lower, y=birthRateRatio_mean, ymax=birthRateRatio_hdi95_upper), color = "grey", width =0) +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + 
    theme(text=element_text(size=15))+
    ylab("Estimated birth rate ratio \n (edge / center)") +
    theme(legend.position = "none") +xlab("") +
    geom_hline(yintercept = 1, linetype="dashed") +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x = element_blank())


t1_wgs_ratio_posteriors_plot_oldstates_state_clock_violin

ggsave(plot=t1_wgs_ratio_posteriors_plot_oldstates_state_clock_violin ,
       file ="../figures/t1_li_wgs_ratio_posteriors_oldstates_stateclock_violin.png", height = 5, width = 1.5)

## Tumor 2

t2_wgs_ratio_posteriors_plot_oldstates_state_clock_violin_summary <- all_logs_birthRate_df %>% 
    dplyr::filter(migration_model == "unidirectional",
                  clock_model == "state-dependent",
                  states == "oldstates",
                  tumor == "T2") %>% 
    group_by(subset, state) %>% 
    summarize(birthRateRatio_hdi95_lower = hdi(birthRateRatio)[1],
              birthRateRatio_hdi95_upper = hdi(birthRateRatio)[2],
              birthRateRatio_mean = mean(birthRateRatio))

t2_wgs_ratio_posteriors_plot_oldstates_state_clock_violin <- all_logs_birthRate_df %>% 
    dplyr::filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           states == "oldstates",
           tumor == "T2") %>% 
    ggplot(., aes(x=subset, y=birthRateRatio), color = "black",  fill = "black") +
    geom_violin(alpha=0.8, fill = "black") +
    geom_boxplot(width=.1, outlier.alpha = 0, color = "grey") +
    # geom_point(data=t2_wgs_ratio_posteriors_plot_oldstates_state_clock_violin_summary, aes(x= subset, y=birthRateRatio_mean), color = "grey") +
    # geom_errorbar(data=t2_wgs_ratio_posteriors_plot_oldstates_state_clock_violin_summary,
    #               aes(x= subset, ymin=birthRateRatio_hdi95_lower, y=birthRateRatio_mean, ymax=birthRateRatio_hdi95_upper), color = "grey", width =0) +
    # 
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + 
    theme(text=element_text(size=15))+
    ylab("Estimated birth rate ratio \n(edge / center)") +
    theme(legend.position = "none") +xlab("") +
    geom_hline(yintercept = 1, linetype="dashed") +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x = element_blank())

t2_wgs_ratio_posteriors_plot_oldstates_state_clock_violin
ggsave(plot=t2_wgs_ratio_posteriors_plot_oldstates_state_clock_violin ,
       file ="../figures/t2_li_wgs_ratio_posteriors_stateclock_violin.png", height = 5, width = 1.5)


#get HDI intervals to put in text
t1_estimate_summaries <- all_logs_birthRate_df %>% 
    dplyr::filter(tumor == "T1", 
                  migration_model == "unidirectional",
                  clock_model == "state-dependent",
                  states == "oldstates") %>% 
    group_by(state) %>% 
    summarise(birthRate_HDI95_lower = hdi(birthRate, credMass = 0.95)[1], 
              birthRate_HDI95_upper = hdi(birthRate, credMass = 0.95)[2],
              birthRate_HDI50_lower = hdi(birthRate, credMass = 0.5)[1], 
              birthRate_HDI50_upper = hdi(birthRate, credMass = 0.5)[2],
              birthRateRatio_HDI95_lower = hdi(birthRateRatio, credMass = 0.95)[1], 
              birthRateRatio_HDI95_upper = hdi(birthRateRatio, credMass = 0.95)[2],
              birthRateRatio_HDI50_lower = hdi(birthRateRatio, credMass = 0.5)[1], 
              birthRateRatio_HDI50_upper = hdi(birthRateRatio, credMass = 0.5)[2],
              birthRateRatio_mean = mean(birthRateRatio))

t1_estimate_summaries %>% 
    ungroup() %>% 
    #group_by(subset) %>% 
    dplyr::select(contains("birthRateRatio")) %>% 
    distinct()
write_csv(t1_estimate_summaries, file = "../li-application/stats/t1_wgs_birth_rate_estimate_summary.csv")

#get HDI intervals to put in text
t2_estimate_summaries <- all_logs_birthRate_df %>% 
    dplyr::filter(tumor == "T2", 
                  migration_model == "unidirectional",
                  clock_model == "state-dependent",
                  states == "oldstates") %>% 
    group_by(state) %>% 
    summarise(birthRate_HDI95_lower = hdi(birthRate, credMass = 0.95)[1], 
              birthRate_HDI95_upper = hdi(birthRate, credMass = 0.95)[2],
              birthRate_HDI50_lower = hdi(birthRate, credMass = 0.5)[1], 
              birthRate_HDI50_upper = hdi(birthRate, credMass = 0.5)[2],
              birthRateRatio_HDI95_lower = hdi(birthRateRatio, credMass = 0.95)[1], 
              birthRateRatio_HDI95_upper = hdi(birthRateRatio, credMass = 0.95)[2],
              birthRateRatio_HDI50_lower = hdi(birthRateRatio, credMass = 0.5)[1], 
              birthRateRatio_HDI50_upper = hdi(birthRateRatio, credMass = 0.5)[2],
              birthRateRatio_mean = mean(birthRateRatio))

t2_estimate_summaries %>% 
    ungroup() %>% 
#    group_by(subset) %>% 
    dplyr::select(contains("birthRateRatio")) %>% 
    distinct()
write_csv(t2_estimate_summaries, file = "../li-application/stats/t2_wgs_birth_rate_estimate_summary.csv")


####### MCC TREES ###############

#To generate MCC run combined_typed_node_mcc_trees.sh in li-application/out directory

li_mcc_tree_files <- list.files(path = "../li-application/combined",
                             pattern="T[1-2]_wgs_oristates_unidir_1_state.HCCtumor.typed.node.tree", 
                             full.names = TRUE)
li_mcc_tree_files <- li_mcc_tree_files[! grepl(".trees", li_mcc_tree_files )]

for (mcc_file in li_mcc_tree_files) {
    mcc_tree <- read.beast(file=mcc_file)
    migration_model <- ifelse(grepl("unidir", basename(mcc_file)), "unidirectional", "bidirectional")
    tumor_extract <- regmatches(basename(mcc_file), gregexpr("T[12]", basename(mcc_file), perl = TRUE))[[1]]
    if(length(tumor_extract) == 0) {
        tumor_extract <- "T3"
    }
    states_extract <- ifelse(grepl("newstates", basename(mcc_file)), "newstates", "oldstates")
    clock_extract <- ifelse(grepl("strict", basename(mcc_file)), "strict", "state-dependent")
    
    mcc_tree@phylo$tip.label <- toupper(gsub("t[12]", "", mcc_tree@phylo$tip.label))
    hmax <- max(as.numeric(mcc_tree@data$height))
    
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
        geom_tiplab(offset = 0.04*hmax, size = 8) +
        coord_cartesian(clip="off")
    
    x_max_limit <- layer_scales(treeplot)$x$get_limits()
    new_x_limit <- 1.1*x_max_limit[2]
    
    treeplot <- treeplot + xlim(c(0,new_x_limit))
    treeplot_pie <- ggtree::inset(treeplot, pies, width = 0.06, height = 0.06) + theme(legend.position = "none") #+ ggtitle(paste0(tumor_extract, migration_model, rep_extract, sep = " "))
    treeplot_pie <- treeplot_pie + geom_nodelab(aes(label = ifelse(round(as.numeric(posterior),2) < 1, round(as.numeric(posterior), 2), "")),
                                                nudge_x = -0.03*hmax, nudge_y = 0.05*hmax, size  = 8, hjust='right', vjust="bottom")
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

t1_estimate_summaries %>% 
    ungroup() %>% 
    dplyr::select(contains("birthRateRatio"))

t2_estimate_summaries %>% 
    dplyr::select(contains("birthRateRatio"))

