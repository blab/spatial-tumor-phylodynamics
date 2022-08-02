#figureS4.R

##Scripts to make figure S4

library(tidyverse)
library(tumortree)
library(HDInterval)
library(beastio)
#example random sampling diagram
## Local directories
# example_cells <- read_csv("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation/cells_death_rate_validation_pop_500_dr_0.27.csv") %>%
#     tumortree::filter_alive_cells()

example_cells <- read_csv("../eden/simulation_data/cells_death_rate_validation_pop_1000_dr_0.27.csv") %>%
    tumortree::filter_alive_cells()

example_cells <- mark_boundary(cells_to_mark = example_cells, alive_cells = example_cells)

#Write example tumor to CSV
write_csv(example_cells, "../eden/simulation_data/alive_cells_death_rate_validation_pop_1000_dr_0.27.csv")

#To skip to visualization
example_cells <- read_csv(file = "../eden/simulation_data/alive_cells_death_rate_validation_pop_1000_dr_0.27.csv")

set.seed(392)
example_cells_diversified <- example_cells %>% 
    tumortree::sample_alive_cells(., n =50, diversified_sampling = TRUE) %>% 
    dplyr::filter(sampled)

example_cells_random <- example_cells %>% 
    tumortree::sample_alive_cells(., n =50, diversified_sampling = FALSE) %>% 
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

######### Sampling and clock model comparison (posteriors) ###########

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
    
    sampling_extract <- ifelse(grepl("random_sampling", basename(log_file)), "random", "diversified")
    
    #single_death <- any(grepl("single_death", log_file, fixed = TRUE))
    
    
    
    #estimates from mcmc logs
    birthRate_loc0_posteriors <- mcmc_log[,"birthRateCanonical.0"]
    birthRate_loc1_posteriors <- mcmc_log[,"birthRateCanonical.1"]
    
    length(birthRate_loc0_posteriors)
    
    birthRateDiff_posteriors <- birthRate_loc1_posteriors - birthRate_loc0_posteriors
    birthRateRatio_posteriors <- (birthRate_loc1_posteriors/birthRate_loc0_posteriors)
    
    # if (single_death) {
    #     
    #     loc0_posteriors <- mcmc_log[,"birthRateCanonical.0"] -  mcmc_log[,"deathRateCanonical"]
    #     loc1_posteriors <- mcmc_log[,"birthRateCanonical.1"] -  mcmc_log[,"deathRateCanonical"]
    #     
    # } else {
    
    loc0_posteriors <- mcmc_log[,"birthRateCanonical.0"] -  mcmc_log[,"deathRateCanonical.0"]
    loc1_posteriors <- mcmc_log[,"birthRateCanonical.1"] -  mcmc_log[,"deathRateCanonical.1"]
    #}
    
    
    
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
                                "n" = as.numeric(n_extract),
                                "sampling" = sampling_extract)
    
    
    
    return(posterior_df)
    
    
}


# log_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs",
#                         pattern ="n_100_state_clock_estimate_dr.log",
#                         full.names = TRUE,
#                         include.dirs = TRUE)
# 
# log_files2 <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs2",
#                         pattern ="n_100_state_clock_estimate_dr.log",
#                         full.names = TRUE,
#                         include.dirs = TRUE)
# 
# strict_clock_log_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs",
#                                      pattern ="n_100_state_clock_estimate_dr_strict_clock.log",
#                                      full.names = TRUE,
#                                      include.dirs = TRUE)
# 
# strict_clock_log_files2 <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs2",
#                                      pattern ="n_100_state_clock_estimate_dr_strict_clock.log",
#                                      full.names = TRUE)
# 
# log_files_random_sampling_log_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs2",
#                                                   pattern ="*n_100_state_clock_estimate_dr_random_sampling.log",
#                                                   full.names = TRUE,
#                                                   include.dirs = TRUE)
# 
# strict_clock_random_sampling_log_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs2",
#                                                      pattern ="*n_100_state_clock_estimate_dr_random_sampling_strict_clock.log",
#                                                      full.names = TRUE,
#                                                      include.dirs = TRUE)

log_files <- list.files(path = "../eden/logs",
                        pattern ="n_100_state_clock_estimate_dr.log",
                        full.names = TRUE,
                        include.dirs = TRUE)
strict_clock_log_files <- list.files(path = "../eden/logs",
                                     pattern ="n_100_state_clock_estimate_dr_strict_clock.log",
                                     full.names = TRUE,
                                     include.dirs = TRUE)

log_files_random_sampling_log_files <- list.files(path = "../eden/logs",
                        pattern ="n_100_state_clock_estimate_dr_random_sampling.log",
                        full.names = TRUE,
                        include.dirs = TRUE)

strict_clock_random_sampling_log_files <- list.files(path = "../eden/logs",
                                     pattern ="n_100_state_clock_estimate_dr_random_sampling_strict_clock.log",
                                     full.names = TRUE,
                                     include.dirs = TRUE)

logs <- purrr::map(log_files, readLog)
logs2 <- purrr::map(log_files2, readLog)
strict_clock_logs <- purrr::map(strict_clock_log_files, readLog)
strict_clock_logs2 <- purrr::map(strict_clock_log_files2, readLog)
random_sampling_logs <- purrr::map(log_files_random_sampling_log_files, readLog)
strict_clock_random_sampling_logs <- purrr::map(strict_clock_random_sampling_log_files, readLog)

## Local records of logs
# logs <- readRDS(file="/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs/all_sampled_sizes_mcmc_logs.rds")
# strict_clock_logs <- readRDS(file="/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs/strict_clock_all_sampled_sizes_mcmc_logs.rds")
# strict_clock_random_sampling_logs <- readRDS(file="/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs/strict_clock_random_sampling_mcmc_logs.rds")
# random_sampling_logs <- readRDS(file="/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs/random_sampling_mcmc_logs.rds")

growth_rate_posteriors_df <- purrr::map2(logs, log_files, 
                                      function(l, f) calculate_growth_rate_posterior_random_sampling(mcmc_log = l,
                                                                                                     log_file = f)) %>% 
    bind_rows

# growth_rate_posteriors_d2f <- purrr::map2(logs2, log_files2, 
#                                          function(l, f) calculate_growth_rate_posterior_random_sampling(mcmc_log = l,
#                                                                                                         log_file = f)) %>% 
#     bind_rows
# 
# growth_rate_posteriors_df_comb <- bind_rows(list(growth_rate_posteriors_df, growth_rate_posteriors_d2f)) %>% 
#     filter(minBirthRateESS > 200) %>% 
#     arrange(-minBirthRateESS) %>% 
#     distinct(., dr, .keep_all = TRUE)

strict_clock_growth_rate_posteriors_df <- purrr::map2(strict_clock_logs, strict_clock_log_files, 
                                                   function(l, f) calculate_growth_rate_posterior_random_sampling(mcmc_log = l,
                                                                                                                  log_file = f)) %>% 
    bind_rows

# strict_clock_growth_rate_posteriors_df2 <- purrr::map2(strict_clock_logs2, strict_clock_log_files2, 
#                                                       function(l, f) calculate_growth_rate_posterior_random_sampling(mcmc_log = l,
#                                                                                                                      log_file = f)) %>% 
#     bind_rows
# 
# strict_clock_growth_rate_posteriors_df_comb <- bind_rows(list(strict_clock_growth_rate_posteriors_df, strict_clock_growth_rate_posteriors_df2 )) %>% 
#     filter(minBirthRateESS > 200) %>% 
#     arrange(-minBirthRateESS) %>% 
#     distinct(., dr, .keep_all = TRUE)

random_sampling_posteriors_df <- purrr::map2(random_sampling_logs, log_files_random_sampling_log_files, 
                                                      function(l, f) calculate_growth_rate_posterior_random_sampling(mcmc_log = l,
                                                                                                                     log_file = f)) %>% 
    bind_rows

strict_clock_random_sampling_posteriors_df <- purrr::map2(strict_clock_random_sampling_logs, strict_clock_random_sampling_log_files, 
                                                      function(l, f) calculate_growth_rate_posterior_random_sampling(mcmc_log = l,
                                                                                                                     log_file = f)) %>% 
    bind_rows

#true simulated rates
#From extract_validation_sims_rates.R
weighted_sim_rates <- read_csv("../eden/stats/validation_growth_and_death_rates_weighted.csv") %>% 
    dplyr::mutate(true_birth_rate_diff_weighted = mean_edge_birth_rate - mean_center_birth_rate) %>% 
    dplyr::select(dr, true_birth_rate_diff_weighted)

#add true growth rates 
growth_rate_posteriors_df <- growth_rate_posteriors_df %>% 
    dplyr::mutate(dr = as.numeric(dr)) %>%
    dplyr::left_join(., weighted_sim_rates, by = "dr") %>% 
    add_column("clock_model" = "state_dependent")

strict_clock_growth_rate_posteriors_df <- strict_clock_growth_rate_posteriors_df %>% 
    dplyr::mutate(dr = as.numeric(dr)) %>%
    dplyr::left_join(., weighted_sim_rates, by = "dr") %>% 
    add_column("clock_model" = "strict")

growth_rate_posteriors_random_sampling_df <- random_sampling_posteriors_df %>% 
    dplyr::mutate(dr = as.numeric(dr)) %>%
    dplyr::left_join(., weighted_sim_rates, by = "dr") %>% 
    add_column("clock_model" = "state_dependent")

strict_clock_growth_rate_posteriors_random_sampling_df <- strict_clock_random_sampling_posteriors_df %>% 
    dplyr::mutate(dr = as.numeric(dr)) %>%
    dplyr::left_join(., weighted_sim_rates, by = "dr") %>% 
    add_column("clock_model" = "strict")

# random_sampling_growth_rate_posterior_df <- random_sampling_growth_rate_posterior_df %>% 
#     dplyr::mutate(dr = as.numeric(dr)) %>%
#     dplyr::left_join(., weighted_sim_rates, by = "dr") %>% 
#     add_column("sampling" = "random") %>% 
#     add_column("clock_model" = "strict")

#write.csv(growth_rate_posteriors_df, "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/outputs/simulation_analysis/birth_death_rate_posteriors_estimate_dr.csv")
#write.csv(strict_clock_growth_rate_posteriors_df, "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/outputs/simulation_analysis/birth_death_rate_posteriors_estimate_dr_strict_clock.csv")

# write.csv(growth_rate_posteriors_df, "../eden/stats/birth_death_rate_posteriors_estimate_dr.csv")
# write.csv(strict_clock_growth_rate_posteriors_df, "../eden/stats/birth_death_rate_posteriors_estimate_dr_strict_clock.csv")
# write.csv(clock_comparison_growth_rate_posteriors_df, file = "../eden/stats/clock_model_sampling_comparison_posteriors.csv")
# # To skip computation
# growth_rate_posteriors_df <- read_csv("../eden/stats/birth_death_rate_posteriors_estimate_dr.csv")
# strict_clock_growth_rate_posteriors_df <- read_csv("../eden/stats/birth_death_rate_posteriors_estimate_dr_strict_clock.csv")
# random_sampling_growth_rate_posterior_df <- read_csv("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/outputs/simulation_analysis/random_sampling_posteriors.csv")
# 
# clock_comparison_growth_rate_posteriors_df <- bind_rows(list(growth_rate_posteriors_df,
#                                                              strict_clock_growth_rate_posteriors_df,
#                                                              random_sampling_growth_rate_posterior_df)) %>%
#     dplyr::mutate(clock_model = factor(clock_model, levels = c('strict', 'state_dependent'))) %>%
#     dplyr::arrange(clock_model)

#Write all results to local directory (too big for github)
clock_comparison_growth_rate_posteriors_df <- bind_rows(list(growth_rate_posteriors_df,
                                                             strict_clock_growth_rate_posteriors_df,
                                                             growth_rate_posteriors_random_sampling_df,
                                                             strict_clock_growth_rate_posteriors_random_sampling_df)) %>% 
    dplyr::mutate(clock_model = factor(clock_model, levels = c('strict', 'state_dependent'))) %>% 
    dplyr::arrange(clock_model)

# write.csv(clock_comparison_growth_rate_posteriors_df,
#           file ="/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/outputs/simulation_analysis/clock_model_sampling_comparison_posteriors.csv")
# 
clock_comparison_growth_rate_posteriors_df_subset <- clock_comparison_growth_rate_posteriors_df %>% 
    dplyr::filter(n==100, minBirthRateESS > 200)

write.csv(clock_comparison_growth_rate_posteriors_df_subset, file = "../eden/stats/clock_model_sampling_comparison_posteriors.csv")

#To skip computation 
clock_comparison_growth_rate_posteriors_df_subset <- read_csv("../eden/stats/clock_model_sampling_comparison_posteriors.csv")

clock_colors <- c("state_dependent" = "#3A5A40", "strict" = "#ECA966")

facet_labels <- c("strict" = "Strict",
                  "state_dependent" = "State-dependent", 
                  "diversified" = "Diversified", 
                  "random" = "Random")

facet_labeller <- function(variable,value){
    return(facet_labels[value])
}

clock_rate_comparison_posteriors_plot <- clock_comparison_growth_rate_posteriors_df_subset %>% 
    ggplot(., aes(x=true_birth_rate_diff_weighted, y=mean_birth_rate_diff, color = clock_model)) +
    geom_point(aes(color = clock_model)) +
    theme_bw() +
    facet_grid(rows = vars(sampling), cols = vars(clock_model), labeller = facet_labeller) +
    geom_abline(slope=1, intercept = 0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(ymin = birthRate_hdi95_lower, ymax = birthRate_hdi95_upper, color=clock_model), alpha = 0.6, size = 1) +
    xlab("True birth rate difference (edge - center)") +
    ylab("Estimated birth rate difference (edge - center)") +
    scale_color_manual(values = clock_colors) +
    theme(legend.position = "none") +
    theme(text=element_text(size = 15))
clock_rate_comparison_posteriors_plot
ggsave(file = "../figures/clock_rate_comparison_birth_rate_posteriors.png",
       clock_rate_comparison_posteriors_plot, height = 6, width = 6)


