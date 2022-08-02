#script to analyze relationship between number of free cells at birth and fitness (number of progeny)

###########
## SETUP ##
###########


library(tumortree)
library(tidyverse)
library(colortools)
library(shades)
library(scales)

model_colors <- c("boundary_driven" = "#e49a8b", "unrestricted" = "#3C3C3C")

#Local directories
# validation_sim_files_boundary_driven <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation",
#                                                    pattern = "cells_death_rate_validation_pop_1000_dr_0.[0-6]5_i_[0-9].csv",
#                                                    full.names = TRUE)
# 
# validation_sim_files_unrestricted <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation",
#                                                 pattern = "cells_pushing_pop_1000_dr_0.[0-6]5_i_[0-9].csv",
#                                                 full.names = TRUE)

validation_sim_files_boundary_driven <- list.files(path = "../eden/simulation_data",
                                                   pattern = "cells_death_rate_validation_pop_1000_dr_0.[0-6]5_i_[0-9].csv", 
                                                   full.names = TRUE)

validation_sim_files_unrestricted <- list.files(path = "../eden/simulation_data",
                                                pattern = "cells_pushing_pop_1000_dr_0.[0-6]5_i_[0-9].csv", 
                                                full.names = TRUE)

#function to mark if cell has progeny or not
mark_replication <- function(sim_file){
    
    #read in simulate cells
    cells_df <- read_csv(sim_file)
    
    #extract simulation info from filename
    
    
    dr_extract <- regmatches(basename(sim_file),
                             gregexpr("(?<=dr_)[[:digit:]]+.[[:digit:]]", basename(sim_file), perl = TRUE))[[1]]
    # i_extract <- regmatches(basename(sim_file),
    #                         gregexpr("(?<=i_)[[:digit:]]", basename(sim_file), perl = TRUE))[[1]]
    
    unrestricted <- any(grepl("cells_pushing", sim_file, fixed = TRUE))
    
    cells_df <- cells_df %>%
        dplyr::mutate("replicates" = index %in% cells_df$parent_index, 
                      "model" = ifelse(unrestricted, "unrestricted", "boundary_driven")) %>% 
        tibble::add_column("dr" = as.numeric(dr_extract)) %>% 
        #dplyr::mutate("has_progeny" <- as.numeric(replicates)) %>% 
        group_by(dr, model, cells_free) %>% 
        summarize(p_progeny = mean(replicates))
    return(cells_df)
}

cells_boundary_driven <- purrr::map(validation_sim_files_boundary_driven, mark_replication) %>% 
    bind_rows

cells_list_unrestricted <- purrr::map(validation_sim_files_unrestricted, mark_replication) %>% 
    bind_rows

prob_replicate_df <- bind_rows(list(cells_boundary_driven, cells_list_unrestricted))


prob_replicate_df <- prob_replicate_df %>%
    
    group_by(dr, cells_free, model) %>%
    summarize(sd_p_progeny = sd(p_progeny), mean_p_progeny = mean(p_progeny)) %>%
    filter(cells_free != 8)


write_csv(prob_replicate_df, file = "../eden/stats/prob_replicate_v_cells_free.csv")

#To skip computation
prob_replicate_df <- read_csv(file = "../eden/stats/prob_replicate_v_cells_free.csv")
#model_colors <- tumortree::get_color_palette(c("boundary_driven", "unrestricted"))
boundary_driven_ramp_palatte <- colorRampPalette(c(model_colors["boundary_driven"], "white"))


unrestricted_ramp_palatte <- colorRampPalette(c(model_colors["unrestricted"], "white"))

boundary_driven_fitness_plot <- prob_replicate_df %>% 
    filter(cells_free != 8, model == "boundary_driven") %>% 
    ggplot(., aes(x = cells_free, y = mean_p_progeny,
                  color = as.factor(dr/2)) )+
    geom_point(size = 3) +
    geom_line() +
    theme_classic() +
    geom_errorbar(aes(ymin=mean_p_progeny-sd_p_progeny,ymax=mean_p_progeny+sd_p_progeny), width = 0.1) +
    scale_color_manual(values = boundary_driven_ramp_palatte(9)) + 
    labs(color = "death rate") +
    xlab("Number adjacent cells free at cell birth") +
    ylab("P(progeny)") +
    theme(text = element_text(size = 20)) +
    labs(color = expression(alpha))
boundary_driven_fitness_plot
ggsave(file = "../figures/death_rate_versus_fitness_boundary_driven.png", plot = boundary_driven_fitness_plot, height = 5, width = 6)

unrestricted_fitness_plot <- prob_replicate_df %>% 
    filter(cells_free != 8, model == "unrestricted") %>% 
    ggplot(., aes(x = cells_free, y = mean_p_progeny,
                  color = as.factor(dr/2))) +
    geom_point(size = 3) +
    geom_line() +
    theme_classic() +
    geom_errorbar(aes(ymin=mean_p_progeny-sd_p_progeny,ymax=mean_p_progeny+sd_p_progeny), width = 0.1) +
    labs(color = "death rate") +
    scale_color_manual(values = unrestricted_ramp_palatte(9)) + 
    xlab("Number adjacent cells free at cell birth") +
    ylab("P(progeny)") +
    ylim(0,0.8) +
    theme(text = element_text(size = 20)) +
    labs(color = expression(alpha))
unrestricted_fitness_plot
ggsave(file = "../figures/death_rate_versus_fitness_unrestricted.png", plot = unrestricted_fitness_plot, height = 5, width = 6)
