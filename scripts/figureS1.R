#script to analyze relationship between number of free cells at birth and fitness (number of progeny)

###########
## SETUP ##
###########
setwd("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation")

library(tumortree)
library(tidyverse)
library(colortools)
library(shades)
library(scales)

validation_sim_files_boundary_driven <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation",
                                                   pattern = "cells_death_rate_validation_pop_1000_dr_0.[0-6]5_i_[0-9].csv", 
                                                   full.names = TRUE)

validation_sim_files_unrestricted <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation",
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

#put cells from each simulation in a list
# cells_list_reps <- list()
# #read in data
# for (i in reps) {
#   for (dr in death_rates) {
#     cells_list_reps[[as.character(paste0(dr, "_", i))]] <- read.csv(paste0("outputs/raw_simulation_results/validation/cells_death_rate_",
#                                                            dr,
#                                                            "_"
#                                                            ,i,
#                                                            ".csv")) %>% 
#       add_column("rep" = i, "free_diffusion" = FALSE)
#   }
# }

# cells_list_free_diff <- list()
# 
# for (i in reps) {
# 
#     cells_list_free_diff[[as.character(i)]] <- read.csv(paste0("outputs/raw_simulation_results/cells_free_diffusion_",
#                                                                            i,
#                                                                            ".csv")) %>% 
#       add_column("rep" = i, "alpha" = 0, "free_diffusion" = TRUE)
# }

#cells_list_reps <- purrr::map(cells_list_reps, mark_replication)


# cell_reps_df <-  cells_list_reps %>% 
#   bind_rows
# cell_reps_df$has_progeny <- as.numeric(cell_reps_df$replicates)



prob_replicate_df <- prob_replicate_df %>%
    
    group_by(dr, cells_free, model) %>%
    summarize(sd_p_progeny = sd(p_progeny), mean_p_progeny = mean(p_progeny)) %>%
    filter(cells_free != 8)


model_colors <- tumortree::get_color_palette(c("boundary_driven", "unrestricted"))
boundary_driven_ramp_palatte <- colorRampPalette(c(model_colors["boundary_driven"], "white"))


unrestricted_ramp_palatte <- colorRampPalette(c(model_colors["unrestricted"], "white"))

# analogous(model_colors[1])
# analogous(model_colors[2])
# show_col(model_colors)
# show_col(model_colors)
# show_col(saturation(model_colors[1], 0.3)[[1]])
# show_col(saturation(model_colors[2], 0.3)[[1]])

boundary_driven_fitness_plot <- prob_replicate_df %>% 
    filter(cells_free != 8, model == "boundary_driven") %>% 
    ggplot(., aes(x = cells_free, y = mean_p_progeny,
                  color = as.factor(dr)) )+
    geom_point(size = 3) +
    theme_classic() +
    geom_errorbar(aes(ymin=mean_p_progeny-sd_p_progeny,ymax=mean_p_progeny+sd_p_progeny), width = 0.1) +
    scale_color_manual(values = boundary_driven_ramp_palatte(9)) + 
    labs(color = "death rate") +
    xlab("Number adjacent cells free at cell birth") +
    ylab("P(progeny)") +
    
    # ggtitle("Spatial constraints on cell fitness at various death rates") +
    ylim(c(0,0.8)) +
    theme(text = element_text(size = 20))
boundary_driven_fitness_plot
ggsave(file = "../figures/death_rate_versus_fitness_boundary_driven.png", plot = boundary_driven_fitness_plot, height = 5, width = 6)

unrestricted_fitness_plot <- prob_replicate_df %>% 
    filter(cells_free != 8, model == "unrestricted") %>% 
    ggplot(., aes(x = cells_free, y = mean_p_progeny,
                  color = as.factor(dr))) +
    geom_point(size = 3) +
    theme_classic() +
    geom_errorbar(aes(ymin=mean_p_progeny-sd_p_progeny,ymax=mean_p_progeny+sd_p_progeny), width = 0.1) +
    labs(color = "death rate") +
    scale_color_manual(values = unrestricted_ramp_palatte(9)) + 
    xlab("Number adjacent cells free at cell birth") +
    ylab("P(progeny)") +
    #ggtitle("Spatial constraints on cell fitness at various death rates") +
    ylim(0,0.8) +
    theme(text = element_text(size = 20))
unrestricted_fitness_plot
ggsave(file = "../figures/death_rate_versus_fitness_unrestricted.png", plot = unrestricted_fitness_plot, height = 5, width = 6)
