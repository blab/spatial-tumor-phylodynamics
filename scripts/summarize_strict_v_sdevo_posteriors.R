#summarize_strict_v_sdevo_posteriors.R
## Generate results table used in figure 5

library(tidyverse)

calculate_growth_rate_posterior <- function(mcmc_log, log_file) {
    
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


# MCMC logs files
# Need log files from running all validation xmls
# Not included in file
# State-dependent clock (SDevo)

# Local dirs
# log_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs2",
#                         pattern ="*estimate_dr.log",
#                         full.names = TRUE,
#                         include.dirs = TRUE)
# 
# # Strict clock
# strict_clock_log_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs2",
#                         pattern ="*estimate_dr_strict_clock.log",
#                         full.names = TRUE,
#                         include.dirs = TRUE)


# State-dependent clock (SDevo)
log_files <- list.files(path = "../eden/logs",
                        pattern ="*estimate_dr.log",
                        full.names = TRUE,
                        include.dirs = TRUE)

# Strict clock
strict_clock_log_files <- list.files(path = "../eden/logs",
                                     pattern ="*estimate_dr_strict_clock.log",
                                     full.names = TRUE,
                                     include.dirs = TRUE)

# Read in log files
logs <- purrr::map(log_files, readLog)
strict_clock_logs <- purrr::map(strict_clock_log_files, readLog)

#true simulated rates
sim_rates <- read_csv("../eden/stats/sim_validation_rates.csv")

#From extract_validation_sims_rates.R
weighted_sim_rates <- read_csv("../eden/stats/validation_growth_and_death_rates_weighted.csv") %>% 
    dplyr::mutate(true_birth_rate_diff_weighted = mean_edge_birth_rate - mean_center_birth_rate) %>% 
    dplyr::select(dr, true_birth_rate_diff_weighted)

#add true growth rates 
growth_rate_posteriors_df <- growth_rate_posteriors %>% 
    dplyr::mutate(dr = as.numeric(dr)) %>%
    dplyr::left_join(., sim_rates, by = "dr") %>% 
    dplyr::mutate("true_birth_rate_diff" = (mean_edge_growth_rate + mean_edge_death_rate) - (mean_center_growth_rate + mean_center_death_rate)) %>% 
    add_column("clock_model" = "state_dependent")

strict_clock_growth_rate_posteriors_df <- strict_clock_growth_rate_posteriors %>% 
    dplyr::mutate(dr = as.numeric(dr)) %>%
    dplyr::left_join(., sim_rates, by = "dr") %>% 
    dplyr::mutate("true_birth_rate_diff" = (mean_edge_growth_rate + mean_edge_death_rate) - (mean_center_growth_rate + mean_center_death_rate)) %>% 
    add_column("clock_model" = "strict")

#add weighted sim rates
growth_rate_posteriors_df <- growth_rate_posteriors_df %>% 
    left_join(., weighted_sim_rates, by = "dr")

strict_clock_growth_rate_posteriors_df <- strict_clock_growth_rate_posteriors_df %>% 
    left_join(., weighted_sim_rates, by = "dr")

# Write results to CSV
write.csv(growth_rate_posteriors_df, "../eden/stats/birth_death_rate_posteriors_estimate_dr.csv")
write.csv(strict_clock_growth_rate_posteriors_df, "../eden/stats/birth_death_rate_posteriors_estimate_dr_strict_clock.csv")
