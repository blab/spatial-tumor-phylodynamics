#!/usr/bin/env Rscript
#Script to calculate true edge and center birth rates from simullations

library(tumortree)
library(dplyr)

#Record of local directory
# sim_cells_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation",
#                               pattern = "cells_death_rate_validation_pop_1000_dr_0.[0-9][0-9].csv",
#                               full.names = TRUE)

sim_cells_files <- list.files(path = "../../simulated_data",
                                     pattern = "cells_death_rate_validation_pop_10000_dr_0.[0-9][0-9][0-9].csv",
                              full.names = TRUE)

extract_rates_dr_validation <- function(all_cells_file, time_interval = 2/24, cutoff_frac = 0.1) {
    
    print(all_cells_file)
    
    all_cells <- read.csv(all_cells_file)
    
    
    dr_extract <- regmatches(basename(all_cells_file),
                             gregexpr("[[:digit:]]+\\.[[:digit:]]*(?=.csv)", basename(all_cells_file), perl = TRUE))[[1]]
    #remove early fraction of time to calculate parameters
    cutoff = max(all_cells$deathdate) * cutoff_frac
    
    
    time_intervals <- seq(cutoff, max(all_cells$deathdate) - time_interval, by = time_interval)
    
    rates_over_time <- purrr::map(time_intervals, function(time) get_compartment_growth_rates_at_t(all_cells = all_cells, time = time, time_interval = time_interval)) %>%
        bind_rows 
    

    mean_death_rate = weighted.mean(rates_over_time$death_rate, rates_over_time$N, na.rm = TRUE)
    mean_edge_growth_rate <- weighted.mean(rates_over_time$growth_rate[rates_over_time$state == "edge" & rates_over_time$time > cutoff], 
                                 rates_over_time$N[rates_over_time$state == "edge" & rates_over_time$time > cutoff], na.rm = TRUE)
    mean_center_growth_rate <- weighted.mean(rates_over_time$growth_rate[rates_over_time$state == "center" & rates_over_time$time > cutoff], rates_over_time$N[rates_over_time$state == "center" & rates_over_time$time > cutoff], na.rm = TRUE)
    
    mean_edge_birth_rate <- weighted.mean(rates_over_time$birth_rate[rates_over_time$state == "edge" & rates_over_time$time > cutoff],
                                 rates_over_time$N[rates_over_time$state == "edge" & rates_over_time$time > cutoff], na.rm = TRUE)
    mean_center_birth_rate <- weighted.mean(rates_over_time$birth_rate[rates_over_time$state == "center" & rates_over_time$time > cutoff], rates_over_time$N[rates_over_time$state == "center" & rates_over_time$time > cutoff], na.rm = TRUE)

    return(data.frame("mean_edge_birth_rate" = mean_edge_birth_rate,
                      "mean_center_birth_rate" = mean_center_birth_rate,
                      "mean_death_rate" = mean_death_rate,
                      "mean_edge_growth_rate" = mean_edge_growth_rate,
                      "mean_center_growth_rate" = mean_center_growth_rate,
                      "dr" = dr_extract))
    
    #return(rates_over_time)
}


growth_rates <- purrr::map(sim_cells_files, function(all_cells_file) extract_rates_dr_validation(all_cells_file, cutoff_frac = 0.25))


growth_rates_df <- growth_rates  %>%
    bind_rows

write.csv(growth_rates_df,
          file = "../eden/stats/validation_growth_and_death_rates_weighted_large.csv")
