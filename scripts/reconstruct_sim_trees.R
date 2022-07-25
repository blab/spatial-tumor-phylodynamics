#reconstruct_sim_trees.R
## Script to reconstruct simulation trees from cells, included state info

library(tumortree)
library(tidyverse)

get_edge_time <- function(leaf, tree_full) {
    mean_edge_time <- weighted.mean(tree_full@data$state[c(leaf, Ancestors(tree_full@phylo, node = leaf, type ="all"))], tree_full@data$branch_length[c(leaf, Ancestors(tree_full@phylo, node = leaf, type ="all"))], na.rm = TRUE)
    
    
    return(mean_edge_time)
}

get_tree_fast <- function(sim_cells_file,
                          overwrite = FALSE,
                          cell_locations_df_file = NULL) {
    print(sim_cells_file)
    
    # Get simulation info from filename
    dr_extract <- regmatches(basename(sim_cells_file),
                             gregexpr("(?<=dr_)[[:digit:]]+.[[:digit:]]+", basename(sim_cells_file), perl = TRUE))[[1]]
    
    # Make tree filename
    tree_file_name <- paste0("../analysis/simtrees/", gsub("csv", "rds", basename(sim_cells_file)))
    
    if(file.exists(tree_file_name) & (! overwrite)) {
        
        #tree_full <- readRDS(tree_file_name)
        
    } else {
        all_cells <- read_csv(sim_cells_file) %>%
            normalize_locs
        
        alive_cells <- all_cells %>%
            filter_alive_cells
        
        tree_full <- convert_all_cells_to_tree_fast(all_cells = all_cells,
                                                    add_all_states = TRUE,
                                                    sampled_cells_indices = all_cells$index, 
                                                    branch_unit = "time",
                                                    cell_locations_df_file = cell_locations_df_file)
        #tree_alive <- prune_simulated_tree(tree_full, sampled_cells_indices = alive_cells$index, add_all_states = TRUE, all_cells = all_cells)
        tree_full@data$frac_time_on_edge <- NA
        tree_full@data$frac_time_on_edge[1:length(tree_full@phylo$tip.label)] <- purrr::map_dbl(1:length(tree_full@phylo$tip.label), function(l) get_edge_time(leaf = l, tree_full = tree_full))
        
        tree_full@data <- tree_full@data %>% 
            dplyr::mutate(n_muts = str_count(mutations, ",") + 1) %>% 
            dplyr::mutate(clock_rate = n_muts/max(node.depth.edgelength(tree_full@phylo))[1])
        
        
        saveRDS(tree_full, file = tree_file_name)
    }
    
}

#To make boundary-driven trees
sim_cells_files <- list.files(path="../simulation_data",
                                      pattern = "cells_death_rate_validation_pop_1000_dr_[0-9].[0-9]+.csv", include.dirs = TRUE, full.names = TRUE)

purrr::map(sim_cells_files, function(cells) get_tree_fast(sim_cells_file = cells,
                                                                   overwrite = TRUE))

#To make unrestricted trees
pushing_sim_cells_files <- list.files(path="../simulation_data",
                                      pattern = "cells_pushing_pop_1000_dr_[0-9].[0-9]+.csv", include.dirs = TRUE, full.names = TRUE)

pushing_sim_cells_locs_files <- list.files(path="../simulation_data",
                                           pattern = "cells_pushing_pop_1000_dr_[0-9].[0-9]+_locs.csv", include.dirs = TRUE, full.names = TRUE)

purrr::map2(pushing_sim_cells_files, pushing_sim_cells_locs_files,
            function(cells, locs) get_tree_fast(sim_cells_file = cells,
                                                overwrite = TRUE,
                                                cell_locations_df_file = locs))
