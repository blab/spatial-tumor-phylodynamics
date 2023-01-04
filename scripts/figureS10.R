#figureS6.R

#scripts to generate figures for Li et al with T1L13 removed

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

##### PLOT POSTERIORS ####

#Record of local directory
# log_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/curr_runs/out",
#                         pattern=".log", 
#                         full.names = TRUE)

log_files <- list.files(path = "../li-application/logs",
                        pattern="T1red_wgs_oristates_unidir_[1-3]_state_rep[0-2].log",
                        full.names = TRUE)

log_files <- log_files[! grepl("chain1",log_files)] 

#log_files <- log_files[! grepl("chain1", log_files)]
process_logs <- function(log_file) {
    print(log_file)
    migration_model <- ifelse(grepl("unidir", basename(log_file)), "unidirectional", "bidirectional")
    tumor_extract <- regmatches(basename(log_file), gregexpr("T[12]", basename(log_file), perl = TRUE))[[1]]
    states_extract <- ifelse(grepl("newstates", basename(log_file)), "newstates", "oldstates")
    clock_extract <- ifelse(grepl("strict", basename(log_file)), "strict", "state-dependent")
    subset_extract <- regmatches(basename(log_file),
                                 gregexpr("(?<=unidir_)[0-9]", basename(log_file), perl = TRUE))[[1]]
    if (length(subset_extract) == 0) {
        subset_extract <- "1"
    }
    
    log <- readLog(log_file, burnin = 0.2)
    log_df <- as.data.frame(log) %>%
        add_column("migration_model" = migration_model,
                   "tumor" = tumor_extract,
                   "states" = states_extract,
                   "subset" = subset_extract,
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

# Also for supplement -- compare published states versus 10% diameter cutoff
## Tumor 1
t1_wgs_posteriors_plot <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           tumor == "T1") %>% 
    mutate("category" = paste0(state, subset)) %>% 
    ggplot(., aes(x=birthRate), color = "black") +
    geom_density(aes(fill = state, group = category), alpha=0.6) +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + scale_fill_manual(values=colors_loc) +
    theme(text=element_text(size=20))+
    xlab("Estimated birth rate") +
    theme(legend.position = "none") +ylab("") 
t1_wgs_posteriors_plot
ggsave(plot=t1_wgs_posteriors_plot ,file ="../figures/t1_li_wgs_posteriors_rt1l13.png", height = 5, width = 5)


## Tumor 1
t1_wgs_ratio_posteriors_plot<- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           clock_model == "state-dependent",
           tumor == "T1") %>% 
    ggplot(., aes(x=birthRateRatio), color = "black",  fill = "black") +
    geom_density(alpha=0.6, fill = "black", aes(group = subset)) +
    #facet_grid(cols = vars(tumor)) +
    theme_classic() + 
    theme(text=element_text(size=20))+
    xlab("Estimated birth rate ratio (edge / center)") +
    theme(legend.position = "none") +ylab("") +
    geom_vline(xintercept = 1, linetype="dashed")

print(mean(all_logs_birthRate_df$birthRateRatio))
hdi(all_logs_birthRate_df$birthRateRatio, credMass = 0.95)
t1_wgs_ratio_posteriors_plot
ggsave(plot=t1_wgs_ratio_posteriors_plot,
       file ="../figures/t1_li_wgs_ratio_posteriors_oristates_stateclock_rt1l13.png", height = 5, width = 5)


####### MCC TREES ###############

#To generate MCC run combined_typed_node_mcc_trees.sh in li-application/out directory

mcc_file <- "../li-application/logs/T1red_wgs_oristates_unidir_state_comb.HCCtumor_mcc.tree"

#only do unidirectional


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
        geom_tiplab(offset = 0.05, size = 6) +
        coord_cartesian(clip="off")
    
    treeplot_pie <- treeplot_pie <- ggtree::inset(treeplot, pies, width = 0.06, height = 0.06) + theme(legend.position = "none") +
        geom_nodelab(aes(label = ifelse(as.numeric(posterior) <= 1, round(as.numeric(posterior), 2), "")),
                     nudge_x = -0.05, nudge_y = 0.2, size  = 6, hjust='right')
    
    fig_file <- paste0("../figures/", tumor_extract,  "_",
                       migration_model, "_",
                       clock_extract,  "_",
                       states_extract, "_",
                       
                       "rt1l13_mcc_tree.png")
    print(fig_file)
    ggsave(plot=treeplot_pie,
           file=fig_file, height = 10, width = 8)


