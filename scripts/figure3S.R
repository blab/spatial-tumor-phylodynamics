##visualize_eden_strict_clock_runs.R

library(tumortree)
library(tidyverse)
library(HDInterval)
library(coda)
library(beastio)

setwd("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/")

calculate_growth_rate_posterior_strict_clock <- function(mcmc_log, log_file) {
    
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
                                "n" = as.numeric(n_extract))
    
    
    
    return(posterior_df)
    
    
}

# logs_bulk_single_death <- read_mcmc_logs(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/bulk/logs",
#                                    pattern = "[0-9]_single_death.log", burnin = 0.2, different_lengths = TRUE)

# logs_bulk <- read_mcmc_logs(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/bulk/logs",
#                                          pattern = "[0-9].log", burnin = 0.2, different_lengths = TRUE)

log_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs2",
                             pattern ="*strict_clock.log",
                             full.names = TRUE,
                             include.dirs = TRUE)

logs <- purrr::map(log_files, readLog)

saveRDS(logs, file = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs/strict_clock_mcmc_logs.rds")

growth_rate_posteriors <- purrr::map2(logs, log_files, 
                                           function(l, f) calculate_growth_rate_posterior_strict_clock(mcmc_log = l,
                                                                                               log_file = f)) %>% 
    bind_rows

sim_rates <- read_csv("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/outputs/simulation_analysis/sim_validation_rates.csv")

birth_rate_strict_clock_fig <- growth_rate_posteriors %>% 
    filter(minBirthRateESS > 100) %>% 
    mutate(dr = as.numeric(dr)) %>%
    left_join(., sim_rates, by = "dr") %>% 
    mutate("true_birth_rate_diff" = (mean_edge_growth_rate + mean_edge_death_rate) - (mean_center_growth_rate + mean_center_death_rate)) %>% 
    ggplot(aes(x=true_birth_rate_diff, y=mean_birth_rate_diff), color = "black") +
    geom_point(color = "black") +
    theme_classic() +
    geom_abline(slope=1, intercept = 0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(ymin = birthRate_hdi95_lower, ymax = birthRate_hdi95_upper), color = "black", alpha = 0.8) +
    xlab("True birth rate difference (edge - center)") +
    ylab("Estimated birth rate difference (edge - center)") +
    theme(text = element_text(size = 15)) +
    facet_wrap(~n)
birth_rate_strict_clock_fig
ggsave(file = "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/manuscript/figures/birth_rate_strict_clock_fig.png", birth_rate_strict_clock_fig, height = 5, width = 7)    
