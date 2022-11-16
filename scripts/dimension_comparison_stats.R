#dimension_comparison_stats.R
library(tidyverse)

#To run this need to run SDevo on physicell 2D neutral + 3D neutral w/ boundary-driven growth
## Requires log files and MCC trees

##### 3D branching signals #####
mcc_tree_files <- list.files(path = "../physicell/trees/3D_neut_bdg/diversified_100",
                             pattern="mcc",
                             full.names = TRUE)

# Filter to desired runs (d0.1, 100 samples and estimated state-dependent death rates)
mcc_tree_files <- mcc_tree_files[grepl("d0.1", mcc_tree_files)]
mcc_tree_files <- mcc_tree_files[! grepl("n_[0-9]", mcc_tree_files)]
mcc_tree_files <- mcc_tree_files[! grepl("single_death", mcc_tree_files)]

# Dataframe to keep track of branching stats
terminal_branch_ratios_df_3d <- data.frame()


for (mcc_file in mcc_tree_files) {
    
    print(mcc_file)
    mcc_tree <- read.beast(file=mcc_file)
    
    # log_file <- paste0("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/logs/3D_neut_bdg/diversified_100/",
    #                    gsub("_mcc.tree", "", basename(mcc_file)), ".log")
    
    log_file <- paste0("../physicell/logs/3D_neut_bdg/diversified_100/",
                       gsub("_mcc.tree", "", basename(mcc_file)), ".log")
    
    log <- readLog(log_file)
    ess <- coda::effectiveSize(log)
    
    # Only use trees that meet ESS threshold
    if (min(c(ess["birthRateCanonical.1"], ess["birthRateCanonical.0"])) < 200) {
        next
    }
    
    
    dr_extract <- qdapRegex::ex_between(basename(mcc_file), "_d", "_")[[1]]
    
    rep_extract <- regmatches(basename(mcc_file),
                              gregexpr("(?<=_s)[0-9]+", basename(mcc_file), perl = TRUE))[[1]]
    
    
    terminal_branch_length_df <- data.frame("state" = mcc_tree@data$type[match(1:100,mcc_tree@data$node)],
                                            "terminal_branch_length" = mcc_tree@phylo$edge.length[match(1:100,mcc_tree@phylo$edge[,2])])
    
    summary_terminal_branch <- terminal_branch_length_df %>% 
        group_by(state) %>% 
        summarise(mean_terminal_branch_length = mean(terminal_branch_length))
    
    ratio_means <- summary_terminal_branch$mean_terminal_branch_length[summary_terminal_branch$state == "loc1"] /summary_terminal_branch$mean_terminal_branch_length[summary_terminal_branch$state == "loc0"]
    
    terminal_branch_ratios_df_3d <- bind_rows(list(terminal_branch_ratios_df_3d, data.frame("dr" = dr_extract,
                                                                                            "terminal_branch_length_ratio" = ratio_means )))
    
}

##### 2D branching signals #####
# Record of local directory
# mcc_tree_files_2d <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/trees/2D_neut_bdg/diversified_100",
#                              pattern="mcc",
#                              full.names = TRUE)

mcc_tree_files_2d <- list.files(path = "../physicell/trees/2D_neut_bdg/diversified_100",
                                pattern="mcc", 
                                full.names = TRUE)

# Filter tree files to d0.1
mcc_tree_files_2d <- mcc_tree_files_2d[grepl("d0.1", mcc_tree_files_2d)]
mcc_tree_files_2d <- mcc_tree_files_2d[! grepl("n_[0-9]", mcc_tree_files_2d)]
mcc_tree_files_2d <- mcc_tree_files_2d[! grepl("single_death", mcc_tree_files_2d)]


#Create data.frame to store terminal branch ratios for each tree

terminal_branch_ratios_df_2d <- data.frame()

for (mcc_file in mcc_tree_files_2d ) {
    
    print(mcc_file)
    mcc_tree <- read.beast(file=mcc_file)
    
    # log_file <- paste0("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/logs/2D_neut_bdg/diversified_100/",
    #                    gsub("_mcc.tree", "", basename(mcc_file)), ".log")
    
    log_file <- paste0("..physicell/logs/2D_neut_bdg/diversified_100/",
                       gsub("_mcc.tree", "", basename(mcc_file)), ".log")
    
    
    
    log <- readLog(log_file)
    ess <- coda::effectiveSize(log)
    
    if (min(c(ess["birthRateCanonical.1"], ess["birthRateCanonical.0"])) < 200) {
        next
    }
    
    
    
    dr_extract <- qdapRegex::ex_between(basename(mcc_file), "_d", "_")[[1]]
    
    rep_extract <- regmatches(basename(mcc_file),
                              gregexpr("(?<=_s)[0-9]+", basename(mcc_file), perl = TRUE))[[1]]
    
    
    terminal_branch_length_df <- data.frame("state" = mcc_tree@data$type[match(1:100,mcc_tree@data$node)],
                                            "terminal_branch_length" = mcc_tree@phylo$edge.length[match(1:100,mcc_tree@phylo$edge[,2])])
    
    summary_terminal_branch <- terminal_branch_length_df %>% 
        group_by(state) %>% 
        summarise(mean_terminal_branch_length = mean(terminal_branch_length))
    
    ratio_means <- summary_terminal_branch$mean_terminal_branch_length[summary_terminal_branch$state == "loc1"] /summary_terminal_branch$mean_terminal_branch_length[summary_terminal_branch$state == "loc0"]
    
    terminal_branch_ratios_df_2d <- bind_rows(list(terminal_branch_ratios_df_2d, data.frame("dr" = dr_extract, "terminal_branch_length_ratio" = ratio_means )))
    
}

### summary plot of terminal branch length ratios #####

terminal_branch_ratios_df_3d <- terminal_branch_ratios_df_3d %>% 
    add_column("sim" = "3D")

terminal_branch_ratios_df_2d <- terminal_branch_ratios_df_2d %>% 
    add_column("sim" = "2D")

terminal_branch_ratios_df_comb <- bind_rows(list(terminal_branch_ratios_df_2d,terminal_branch_ratios_df_3d))

# Save branching stats
write_csv(terminal_branch_ratios_df_comb, "../physicell/stats/terminal_branch_dim_comparison.csv")
