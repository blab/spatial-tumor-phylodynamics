#figure4.R

##script to generate visualizations for figure 4

library(tumortree)
library(tidyverse)
library(cowplot)
library(viridis)
library(spatstat)
library(beastio)
library(HDInterval)
library(coda)

figures_dir <- "../figures"


sim_colors <- get_color_palette(c("boundary_driven", "unrestricted"))

colors_edge_center <- get_color_palette(c("edge", "center"))


sim_rates <- read_csv("../eden/stats/sim_validation_rates.csv")

#From extract_validation_sims_rates.R
weighted_sim_rates <- read_csv("../eden/stats/validation_growth_and_death_rates_weighted2.csv") %>% 
    dplyr::mutate(true_birth_rate_diff_weighted = mean_edge_birth_rate - mean_center_birth_rate) %>% 
    dplyr::select(dr, true_birth_rate_diff_weighted)

#### Figure 4A (random sampling eden simulation) #####

calculate_growth_rate_posterior_random_sampling <- function(mcmc_log, log_file) {
    
    print(log_file)
    if (is.null(log_file)) {
        
        return(NA)
    }
    #calculate effective sample size
    
    ess <- try(coda::effectiveSize(mcmc_log))
    if (! is.numeric(ess)) {
        ess <- 0
    }
    
    
    
    
    #extract simulation info from filename
    dr_extract <- regmatches(basename(log_file),
                             gregexpr("(?<=dr_)[[:digit:]]+.[[:digit:]]+", basename(log_file), perl = TRUE))[[1]]
    
    length(dr_extract)
    
    n_extract <- regmatches(basename(log_file),
                            gregexpr("(?<=n_)[[:digit:]]+", basename(log_file), perl = TRUE))[[1]]
    
    #single_death <- any(grepl("single_death", log_file, fixed = TRUE))
    
    
    
    #estimates from mcmc logs
    birthRate_loc0_posteriors <- mcmc_log[,"birthRateCanonical.0"]
    birthRate_loc1_posteriors <- mcmc_log[,"birthRateCanonical.1"]
    
    length(birthRate_loc0_posteriors)
    
    birthRateDiff_posteriors <- birthRate_loc1_posteriors - birthRate_loc0_posteriors
    birthRateRatio_posteriors <- (birthRate_loc1_posteriors/birthRate_loc0_posteriors)
    

    loc0_posteriors <- mcmc_log[,"birthRateCanonical.0"] -  mcmc_log[,"deathRateCanonical.0"]
    loc1_posteriors <- mcmc_log[,"birthRateCanonical.1"] -  mcmc_log[,"deathRateCanonical.1"]

    
    
    growthRateDiff_posteriors <- loc1_posteriors - loc0_posteriors
    growthRateRatio_posteriors <- (loc1_posteriors/loc0_posteriors)
    
    
    posterior_df <-  data.frame("mean_growth_rate_diff" = mean(growthRateDiff_posteriors),
                                "median_growth_rate_diff" = median(growthRateDiff_posteriors),
                                "mean_growth_rate_ratio" = mean(growthRateRatio_posteriors),
                                "median_growth_rate_ratio" = median(growthRateRatio_posteriors),
                                "mean_birth_rate_diff" = mean(birthRateDiff_posteriors),
                                "median_birth_rate_diff" = median(birthRateDiff_posteriors),
                                "mean_birth_rate_ratio" = mean(birthRateRatio_posteriors),
                                "median_birth_rate_ratio" = median(birthRateRatio_posteriors),
                                "growthRate_hdi95_lower" = hdi(growthRateDiff_posteriors,  credMass = 0.95)[1],
                                "growthRate_hdi95_upper" = hdi(growthRateDiff_posteriors,  credMass = 0.95)[2],
                                "growthRate_hdi85_lower" = hdi(growthRateDiff_posteriors,  credMass = 0.85)[1],
                                "growthRate_hdi85_upper" = hdi(growthRateDiff_posteriors,  credMass = 0.85)[2],
                                "growthRate_hdi75_lower" = hdi(growthRateDiff_posteriors, credMass = 0.75)[1],
                                "growthRate_hdi75_upper" = hdi(growthRateDiff_posteriors,  credMass = 0.75)[2],
                                "birthRate_hdi95_lower" = hdi(birthRateDiff_posteriors,  credMass = 0.95)[1],
                                "birthRate_hdi95_upper" = hdi(birthRateDiff_posteriors,  credMass = 0.95)[2],
                                "birthRate_hdi85_lower" = hdi(birthRateDiff_posteriors,  credMass = 0.85)[1],
                                "birthRate_hdi85_upper" = hdi(birthRateDiff_posteriors,  credMass = 0.85)[2],
                                "birthRate_hdi75_lower" = hdi(birthRateDiff_posteriors, credMass = 0.75)[1],
                                "birthRate_hdi75_upper" = hdi(birthRateDiff_posteriors,  credMass = 0.75)[2],
                                "dr" = dr_extract,
                                "minBirthRateESS" = min(ess["birthRateCanonical.1"], ess["birthRateCanonical.0"]),
                                "n" = as.numeric(n_extract))
    
    
    
    return(posterior_df)
    
    
}

## Record of local location
# log_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs2",
#                         pattern ="*estimate_dr_random_sampling.log",
#                         full.names = TRUE,
#                         include.dirs = TRUE)


log_files <- list.files(path = "../eden/logs",
                        pattern ="*estimate_dr_random_sampling.log",
                        full.names = TRUE,
                        include.dirs = TRUE)

logs <- purrr::map(log_files, readLog)

#saveRDS(logs, file = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs/random_sampling_mcmc_logs.rds")
#logs <- readRDS(file = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs/random_sampling_mcmc_logs.rds")


growth_rate_posteriors <- purrr::map2(logs, log_files, 
                                      function(l, f) calculate_growth_rate_posterior_random_sampling(mcmc_log = l,
                                                                                                     log_file = f)) %>% 
    bind_rows

write_csv(growth_rate_posteriors, file="../eden/stats/random_sampling_growth_rate_posteriors.csv")

# To skip to visualization without recalculating posterior means + intervals
growth_rate_posteriors <- read_csv(file="../eden/stats/random_sampling_growth_rate_posteriors.csv")

birth_rate_random_sampling_fig <- growth_rate_posteriors %>% 
    filter(minBirthRateESS > 200) %>% 
    mutate(dr = as.numeric(dr)) %>%
    left_join(., weighted_sim_rates, by = "dr") %>% 
    #mutate("true_birth_rate_diff" = (mean_edge_growth_rate + mean_edge_death_rate) - (mean_center_growth_rate + mean_center_death_rate)) %>% 
    ggplot(aes(x=true_birth_rate_diff_weighted, y=mean_birth_rate_diff), color = "black") +
    geom_point(color = "black") +
    theme_classic() +
    geom_abline(slope=1, intercept = 0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(ymin = birthRate_hdi95_lower, ymax = birthRate_hdi95_upper), color = "black", alpha = 0.8) +
    xlab("True birth rate difference (edge - center)") +
    ylab("Estimated birth rate difference (edge - center)") +
    theme(text = element_text(size = 15)) #+
#facet_wrap(~n)
birth_rate_random_sampling_fig
ggsave(file = "../figures/birth_rate_random_sampling_fig.png", birth_rate_random_sampling_fig, height = 5, width = 5)    

#example random sampling diagram

#local directory
# example_cells <- read_csv("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation/cells_death_rate_validation_pop_1000_dr_0.27.csv") %>% 
#     tumortree::filter_alive_cells()
example_cells <- read_csv("../eden/simulation_data/cells_death_rate_validation_pop_1000_dr_0.27.csv") %>%
    tumortree::filter_alive_cells()
example_cells <- mark_boundary(cells_to_mark = example_cells, alive_cells = example_cells)

#Write example tumor to CSV
write_csv(example_cells, "../eden/simulation_data/alive_cells_death_rate_validation_pop_1000_dr_0.27.csv")

#To skip to visualization
example_cells <- read_csv(file = "../eden/simulation_data/alive_cells_death_rate_validation_pop_1000_dr_0.27.csv")

#Note this is repeat code from figure S4, but used to generate inset
set.seed(392)
example_cells_diversified <- example_cells %>% 
    tumortree::sample_alive_cells(., n =100, diversified_sampling = TRUE) %>% 
    dplyr::filter(sampled)

example_cells_random <- example_cells %>% 
    tumortree::sample_alive_cells(., n =100, diversified_sampling = FALSE) %>% 
    dplyr::filter(sampled)

colors_edge_center <- get_color_palette(names = c("edge", "center"))

example_tumor_inset_diversified <- ggplot(example_cells, aes(locx, locy, color = ifelse(est_edge == 1, "edge", "center"))) +
    geom_point() + theme_void() +
    #geom_point(data = example_cells_diversified, shape = 21, aes(fill = ifelse(est_edge == 1, "edge", "center")), color = "black", size =3)+
    geom_point(data = example_cells_diversified, shape = 1, fill = "black",  color = "black", alpha = 1, size = 3)+
    scale_fill_manual(values = colors_edge_center) +
    scale_color_manual(values = colors_edge_center) +
    theme(legend.position = "none") 

example_tumor_inset_diversified
ggsave(file = "../figures/diversified_sampling_example_tumor.png", example_tumor_inset_diversified, height = 5, width = 5)    


example_tumor_inset_random <- ggplot(example_cells, aes(locx, locy, color = ifelse(est_edge == 1, "edge", "center"))) +
    geom_point() + theme_void() +
    #geom_point(data = example_cells_diversified, shape = 21, aes(fill = ifelse(est_edge == 1, "edge", "center")), color = "black", size =3)+
    geom_point(data = example_cells_random, shape = 1, fill = "black",  color = "black", alpha = 1, size = 3)+
    scale_fill_manual(values = colors_edge_center) +
    scale_color_manual(values = colors_edge_center) +
    theme(legend.position = "none") 

example_tumor_inset_random
ggsave(file = "../figures/random_sampling_example_tumor.png", example_tumor_inset_random, height = 5, width = 5)    

#### Bulk sampling (Figure 4B) #####

#Function specifically for bulk sequencing, very similar to previous function, so could clean this up later
calculate_growth_rate_posterior_bulk <- function(mcmc_log, log_file) {
    
    print(log_file)
    if (is.null(log_file)) {
        
        return(NA)
    }
    #calculate effective sample size
    
    ess <- try(coda::effectiveSize(mcmc_log))
    if (! is.numeric(ess)) {
        ess <- 0
    }
    
    
    
    
    #extract simulation info from filename
    dr_extract <- regmatches(basename(log_file),
                             gregexpr("(?<=dr_)[[:digit:]]+.[[:digit:]]+", basename(log_file), perl = TRUE))[[1]]
    
    length(dr_extract)
    
    # n_extract <- regmatches(basename(log_file),
    #                         gregexpr("(?<=n_)[[:digit:]]+", basename(log_file), perl = TRUE))[[1]]
    
    single_death <- any(grepl("single_death", log_file, fixed = TRUE))
    
    
    
    #estimates from mcmc logs
    birthRate_loc0_posteriors <- mcmc_log[,"birthRateCanonical.0"]
    birthRate_loc1_posteriors <- mcmc_log[,"birthRateCanonical.1"]
    
    length(birthRate_loc0_posteriors)
    
    birthRateDiff_posteriors <- birthRate_loc1_posteriors - birthRate_loc0_posteriors
    birthRateRatio_posteriors <- (birthRate_loc1_posteriors/birthRate_loc0_posteriors)
    
    if (single_death) {
        
        loc0_posteriors <- mcmc_log[,"birthRateCanonical.0"] -  mcmc_log[,"deathRateCanonical"]
        loc1_posteriors <- mcmc_log[,"birthRateCanonical.1"] -  mcmc_log[,"deathRateCanonical"]
        
    } else {
        
        loc0_posteriors <- mcmc_log[,"birthRateCanonical.0"] -  mcmc_log[,"deathRateCanonical.0"]
        loc1_posteriors <- mcmc_log[,"birthRateCanonical.1"] -  mcmc_log[,"deathRateCanonical.1"]
    }
    
    
    
    growthRateDiff_posteriors <- loc1_posteriors - loc0_posteriors
    growthRateRatio_posteriors <- (loc1_posteriors/loc0_posteriors)
    
    
    posterior_df <-  data.frame("mean_growth_rate_diff" = mean(growthRateDiff_posteriors),
                                "median_growth_rate_diff" = median(growthRateDiff_posteriors),
                                "mean_growth_rate_ratio" = mean(growthRateRatio_posteriors),
                                "median_growth_rate_ratio" = median(growthRateRatio_posteriors),
                                "mean_birth_rate_diff" = mean(birthRateDiff_posteriors),
                                "median_birth_rate_diff" = median(birthRateDiff_posteriors),
                                "mean_birth_rate_ratio" = mean(birthRateRatio_posteriors),
                                "median_birth_rate_ratio" = median(birthRateRatio_posteriors),
                                "growthRate_hdi95_lower" = hdi(growthRateDiff_posteriors,  credMass = 0.95)[1],
                                "growthRate_hdi95_upper" = hdi(growthRateDiff_posteriors,  credMass = 0.95)[2],
                                "growthRate_hdi85_lower" = hdi(growthRateDiff_posteriors,  credMass = 0.85)[1],
                                "growthRate_hdi85_upper" = hdi(growthRateDiff_posteriors,  credMass = 0.85)[2],
                                "growthRate_hdi75_lower" = hdi(growthRateDiff_posteriors, credMass = 0.75)[1],
                                "growthRate_hdi75_upper" = hdi(growthRateDiff_posteriors,  credMass = 0.75)[2],
                                "birthRate_hdi95_lower" = hdi(birthRateDiff_posteriors,  credMass = 0.95)[1],
                                "birthRate_hdi95_upper" = hdi(birthRateDiff_posteriors,  credMass = 0.95)[2],
                                "birthRate_hdi85_lower" = hdi(birthRateDiff_posteriors,  credMass = 0.85)[1],
                                "birthRate_hdi85_upper" = hdi(birthRateDiff_posteriors,  credMass = 0.85)[2],
                                "birthRate_hdi75_lower" = hdi(birthRateDiff_posteriors, credMass = 0.75)[1],
                                "birthRate_hdi75_upper" = hdi(birthRateDiff_posteriors,  credMass = 0.75)[2],
                                "dr" = dr_extract,
                                "minBirthRateESS" = min(ess["birthRateCanonical.1"], ess["birthRateCanonical.0"]),
                                "single_death" = single_death)
    
    
    
    return(posterior_df)
    
    
}


#Record of local directories
bulk_log_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/bulk/logs",
                             pattern ="*[0-9].log",
                             full.names = TRUE,
                             include.dirs = TRUE)
# 

bulk_log_files <- list.files(path = "../eden/logs/bulk",
                             pattern ="*[0-9].log",
                             full.names = TRUE,
                             include.dirs = TRUE)


logs_bulk <- purrr::map(bulk_log_files, readLog)

#saveRDS(logs_bulk, file = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/bulk/logs/all_mcmc_logs.rds")

#logs_bulk <- readRDS("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/bulk/logs/all_mcmc_logs.rds")

growth_rate_posteriors_bulk <- purrr::map2(logs_bulk, bulk_log_files, 
                                           function(l, f) calculate_growth_rate_posterior_bulk(mcmc_log = l,
                                                                                               log_file = f)) %>% 
    bind_rows %>% 
    add_column("clock_model" = "state_dependent")


write_csv(growth_rate_posteriors_bulk, file="../eden/stats/bulk_sampling_posteriors.csv")

#To skip computation
growth_rate_posteriors_bulk <- read_csv(file="../eden/stats/bulk_sampling_posteriors.csv")

    
bulk_true_versus_estimated_plot <- growth_rate_posteriors_bulk %>% 
    
    filter(minBirthRateESS > 200) %>% 
    mutate(dr = as.numeric(dr)) %>%
    left_join(.,weighted_sim_rates, by = "dr") %>% 
    #mutate("true_birth_rate_diff" = (mean_edge_growth_rate + mean_edge_death_rate) - (mean_center_growth_rate + mean_center_death_rate)) %>% 
    ggplot(aes(x=true_birth_rate_diff_weighted, y=mean_birth_rate_diff)) +
    geom_point() +
    theme_classic() +
    geom_abline(slope=1, intercept = 0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(ymin = birthRate_hdi95_lower, ymax = birthRate_hdi95_upper)) +
    xlab("True birth rate difference (edge - center)") +
    ylab("Estimated birth rate difference (edge - center)") +
    theme(text = element_text(size = 15)) 
bulk_true_versus_estimated_plot

# x_range_bulk_plot <- diff(ggplot_build(bulk_true_versus_estimated_plot)$layout$panel_scales_x[[1]]$range$range)
# y_range_bulk_plot <- diff(ggplot_build(bulk_true_versus_estimated_plot)$layout$panel_scales_y[[1]]$range$range)
# mean_range_bulk_plot <- mean(x_range_bulk_plot, y_range_bulk_plot)
ggsave("../figures/bulk_true_versus_estimated.png", bulk_true_versus_estimated_plot, height = 5 , width = 5)

###### Example bulk sampled tumor (Figure 4B inset)

#Sampled in bulk_sample_simulated_tumor.R

#generated in bulk_sample_simulated_tumor.R
alive_cells <-read_csv(file = "../eden/simulation_data/example_bulk_sequencing_alive_cells.csv")

punched_cells <- alive_cells %>% 
    filter(! is.na(punch_id))

g <- ggplot(alive_cells, aes(x = locx, y = locy, color = ifelse(est_edge == 1, "edge", "center") )) +
    geom_point(size = 2) +
    scale_color_manual(values = punch_colors)+
    geom_point(data=punched_cells, size = 2, color = "black") +
    theme_void() +
    theme(legend.position = "none")
g
ggsave("../figures/bulk_exampled_tumor_sampled_punches.png", g, height = 5, width = 5)

###### Physicell Visualizations ########

# Using posterior  stats from summarize_physicell_state_clocks_runs.R

############# PHYSICELL 2D NEUTRAL BOUNDARY-DRIVEN GROWTH #####
growth_rate_posteriors_2D_net_bdg_n100 <- read_csv("../physicell/stats/growth_rate_posteriors_2D_neut_bdg_n100.csv")
d_baseplot <- growth_rate_posteriors_2D_net_bdg_n100 %>% 
    dplyr::filter(minBirthRateESS > 200) %>% 
    ggplot(., aes(x = true_birth_rate_diff, y = mean_birth_rate_diff), color = "black") +
    theme_classic() + geom_errorbar(aes(ymin = birthRate_hdi95_lower, ymax = birthRate_hdi95_upper), width = 0, alpha = 0.5, size = 1) +
    geom_point(aes( y = mean_birth_rate_diff), alpha = 0.8, size = 2) + scale_color_viridis(option = "mako", discrete = TRUE, direction = -1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") + xlab("True birth rate difference (edge - center)") +
    ylab("True birth rate difference (edge - center)") + theme(text=element_text(size=15))
    #ggtitle("Physicell 2D Neutral Boundary-Driven Growth") +    theme(text=element_text(size=15))

d_baseplot

ggsave(plot=d_baseplot , file ="../figures/physicell_2D_neutral_bdg.png", height = 5, width = 6)

############# PHYSICELL 3D NEUTRAL BOUNDARY-DRIVEN GROWTH #####

growth_rate_posteriors_3D_net_bdg_n100 <- read_csv("../physicell/stats/growth_rate_posteriors_3D_neut_bdg_n100.csv")
e_baseplot <- growth_rate_posteriors_3D_net_bdg_n100 %>% 
    dplyr::filter(minBirthRateESS > 200, mean_birth_rate_diff > -0.001) %>% #this exludes run with poor convergence
    ggplot(., aes(x = true_birth_rate_diff, y = mean_birth_rate_diff), color = "black") +
    theme_classic() + geom_errorbar(aes(ymin = birthRate_hdi95_lower, ymax = birthRate_hdi95_upper), width = 0, alpha = 0.5, size = 1) +
    geom_point(aes( y = mean_birth_rate_diff), alpha = 0.8, size = 2) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") + xlab("True birth rate difference (edge - center)") +
    theme(text=element_text(size=15)) +
    ylab("True birth rate difference (edge - center)")

e_baseplot
ggsave(plot=e_baseplot, file ="../figures/physicell_3D_neutral_bdg.png", height = 5, width = 5)

############# PHYSICELL 2D SELECTION BOUNDARY-DRIVEN GROWTH #####
growth_rate_posteriors_2D_sel_bdg_n100 <- read_csv("../physicell/stats/growth_rate_posteriors_2D_sel_bdg_n100.csv")

f_baseplot <- growth_rate_posteriors_2D_sel_bdg_n100 %>% 
    dplyr::filter(minBirthRateESS > 200) %>% 
    ggplot(., aes(x = true_birth_rate_diff, y = mean_birth_rate_diff), color = "black") +
    theme_classic() + geom_errorbar(aes(ymin = birthRate_hdi95_lower, ymax = birthRate_hdi95_upper), width = 0, alpha = 0.5, size = 1) +
    geom_point(aes( y = mean_birth_rate_diff), alpha = 0.8, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") + xlab("True birth rate difference (edge - center)") +
    ylab("True birth rate difference (edge - center)") +
    theme(text = element_text(size = 15))

f_baseplot

ggsave(plot=f_baseplot, file ="../figures/physicell_2D_sel_bdg.png", height = 5, width = 5)


