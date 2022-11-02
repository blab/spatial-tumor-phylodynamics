#extended_physicell_viz.R
## This script only summarizes extended biological physicell simulations

library(tidyverse)

#Get means and 95% HPD intervals
## These are calculated from MCMC logs by write_physicell_posteriors_summary.Rscript for each file and combined into posterior_summary.tsv

posterior_summary <- read_tsv("../physicell/stats/posteriors/posterior_summary.tsv") %>% 
    dplyr::mutate("base" = gsub("_diversified_m1_n100", "", run))

#Get true birth rates
## True birth rates (labeled as growth rates in original files) are calculated from weighted means of edge and center birth rates over time
### get_true_rates_physicell.Rscript inputs each XML file and maps to the simulation growth rates recorded in "updated_growth_rate_differences" in the simulation data

true_birth_rates_df <- read_tsv("../physicell/stats/true_birth_rates/true_birth_rates_extended_physicell.tsv") %>% 
    rename("run" = simulation)

#Convert logical to string categories
true_birth_rates_df$edgeAssociated <- ifelse(true_birth_rates_df$edgeAssociated, "edge", "center")

#Pivot
true_birth_rates_df <- true_birth_rates_df %>% 
    pivot_wider(id_cols = run, names_from = edgeAssociated, values_from = mean_growth_rate) %>% 
    dplyr::mutate("true_birth_diff" = edge-center)


posterior_summary <- posterior_summary %>% 
    left_join(.,true_birth_rates_df, by = "run")
posterior_summary %>% 
    filter(minBirthRateESS > 200) %>% 

    ggplot(., aes(x=true_birth_diff, y = mean_birth_rate_diff,  color = simulation)) +
    geom_point() +
    geom_errorbar(aes(ymin = birthRate_hdi95_lower, ymax = birthRate_hdi95_upper), width = 0) +
    geom_abline(slope=1, intercept = 0, linetype = "dashed") +
    theme_classic() + facet_wrap(~simulation)


sse_summary <- posterior_summary %>% 
    group_by(n, simulation) %>% 
    summarise("mse" = mean((mean_birth_rate_diff - true_birth_diff)^2, na.rm = TRUE),
              "se_mse" = sd((mean_birth_rate_diff - true_birth_diff)^2, na.rm = TRUE) / sqrt(n()))


sse_summary %>% 
    ggplot(., aes(x = simulation, y = mse)) +
    geom_point(aes(color=simulation), size = 0.75) +
    theme_classic() + xlab("Sample size") + ylab("MSE of edge - center birth rate difference") +
    geom_errorbar(aes(ymin = mse - se_mse, ymax = mse + se_mse, color = simulation)) +
    theme(text = element_text(size = 15)) +
    #scale_x_continuous(breaks = pretty(clock_comparison_sample_size_results_summary$n, n = 10)) + 
    #scale_color_manual(values = clock_colors) + geom_smooth(aes(color=clock_model), linetype="dashed", se = FALSE, size = 0.5)+
    theme(legend.position = "none")

