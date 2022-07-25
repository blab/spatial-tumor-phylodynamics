#figure4.R

##script to generate visualizations for figure 4

library(tumortree)
library(tidyverse)
library(cowplot)
library(viridis)
library(spatstat)

figures_dir <- "../figures"


sim_colors <- get_color_palette(c("boundary_driven", "unrestricted"))

colors_edge_center <- get_color_palette(c("edge", "center"))


sim_rates <- read_csv("outputs/simulation_analysis/sim_validation_rates.csv")

## SUBFIGURE growth_rate_effect_estimated_versus_simulated_StDEvo.png
state_dependent_diversified_mcmc_logs <- read_mcmc_logs(path = "outputs/beast_analysis/state_dependent_clock_model/validation/logs",
                                                        pattern = "*state_dependent_diversified_sampling\\.log")
state_dependent_diversified_death_rates <- purrr::map_dbl(names(state_dependent_diversified_mcmc_logs), extract_death_rate)

growth_rate_state_clock_diversified_df <- purrr::map2(state_dependent_diversified_mcmc_logs,
                                                state_dependent_diversified_death_rates,
                                                function(x, dr) calc_state_dependent_growthRateDiff_means(mcmc_obj = x,
                                                                                                          dr = dr,
                                                                                                          fixed_dr = 0.1)) %>%
    bind_rows %>%
    left_join(sim_rates, by = "dr") %>%
    dplyr::mutate("mean_growth_rate_diff" = mean_edge_growth_rate - mean_center_growth_rate)


a <- ggplot(growth_rate_state_clock_diversified_df  , aes(x = mean_growth_rate_diff, y = mean), color = "grey") +
    geom_errorbar(alpha = 0.7, aes(ymin = hdi95_lower, ymax = hdi95_upper),  size=1, width=0) +
    geom_errorbar(alpha = 0.4, aes(ymin = hdi85_lower, ymax = hdi85_upper),  size=2, width=0) +
    geom_errorbar(alpha = 0.4, aes(ymin = hdi75_lower, ymax = hdi75_upper),  size=3, width=0) +
    geom_point(size = 2) +
    #geom_hline(yintercept = 0, linetype = "dashed") +
    #geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0) +
    theme_classic() + labs("color" = "") +
    ylab("Estimated state-dependent effect on growth rate") + xlab("True growth rate difference")
a






## SUBFIGURE ##
## read in example tumor
all_cells_boundary_driven <- read_csv("outputs/raw_simulation_results/validation/cells_death_rate_validation_pop_1000_dr_0.12.csv") %>%
    normalize_locs

alive_cells_boundary_driven <- all_cells_boundary_drive %>%
    filter_alive_cells
set.seed(81721)
sampled_cells_boundary_driven_diversified <- sample_alive_cells(alive_cells = alive_cells_boundary_driven,
                                                                n = 100, diversified_sampling = TRUE)
a_sub <- ggplot(sampled_cells_boundary_driven_diversified, aes(x = locx, y = locy)) +
    geom_point(aes(color = sampled, size = as.integer(sampled))) +
    theme_void() + scale_color_manual(values = c("grey", "#809070")) +
    scale_size(range = c(0.5,1)) + theme(legend.position = "none")
a_sub



## SUBFIGURE growth_rate_effect_estimated_versus_simulated_StDEvo_random.png
state_clock_random_mcmc_logs <- read_mcmc_logs(path = "outputs/beast_analysis/state_dependent_clock_model/validation/logs",
                                               pattern = "*state_dependent\\.log",
                                               burnin = 0.5,
                                               different_lengths = TRUE)

state_clock_random_death_rates <- purrr::map_dbl(names(state_clock_random_mcmc_logs), extract_death_rate)


growth_rate_state_clock_random_df <- purrr::map2(state_clock_random_mcmc_logs,
                                                 state_clock_random_death_rates,
                                           function(x, dr) calc_state_dependent_growthRateDiff_means(mcmc_obj = x,
                                                                                                     dr = dr,
                                                                                                     fixed_dr = 0.1)) %>%
    bind_rows %>%
    left_join(sim_rates, by = "dr") %>%
    dplyr::mutate("mean_growth_rate_diff" = mean_edge_growth_rate - mean_center_growth_rate)

b <- ggplot(growth_rate_state_clock_random_df  , aes(x = mean_growth_rate_diff, y = mean), color = "grey") +
    geom_errorbar(alpha = 0.7, aes(ymin = hdi95_lower, ymax = hdi95_upper),  size=1, width=0) +
    geom_errorbar(alpha = 0.4, aes(ymin = hdi85_lower, ymax = hdi85_upper),  size=2, width=0) +
    geom_errorbar(alpha = 0.4, aes(ymin = hdi75_lower, ymax = hdi75_upper),  size=3, width=0) +
    geom_point(size = 2) +
    #geom_hline(yintercept = 0, linetype = "dashed") +
    #geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0) +
    theme_classic() + labs("color" = "") +
    ylab("Estimated state-dependent effect on growth rate") + xlab("True growth rate difference")
b

set.seed(3333)
sampled_cells_boundary_driven_random <- sample_alive_cells(alive_cells = alive_cells_boundary_driven,
                                                           n = 100, diversified_sampling = FALSE)

b_sub <- ggplot(sampled_cells_boundary_driven_random, aes(x = locx, y = locy)) +
    geom_point(aes(color = sampled, size = as.integer(sampled))) +
    theme_void() + scale_color_manual(values = c("grey", "#809070")) +
    scale_size(range = c(0.5,1)) + theme(legend.position = "none")
b_sub
## SUBFIGURE ##


## INSETS WITH COMMOON AXES


global_ymin <- min(c(growth_rate_state_clock_diversified_df$hdi95_lower,growth_rate_state_clock_random_df$hdi95_upper))
global_ymax <- max(c(growth_rate_state_clock_diversified_df$hdi95_lower,growth_rate_state_clock_random_df$hdi95_upper))

global_xmin <- min(c(growth_rate_state_clock_diversified_df$mean_growth_rate_diff))
global_xmax <- max(c(growth_rate_state_clock_diversified_df$mean_growth_rate_diff))



a <- a + coord_cartesian(ylim = c(global_ymin, global_ymax), xlim = c(global_xmin, global_xmax ))
a_sub <- a_sub + border(linetype = "dashed", color = "darkgrey", size = 0.8)
a.with.inset <-
    ggdraw() +
    draw_plot(a)  +
    draw_plot(a_sub, x = 0.25, y = .65, width = .3, height = .3)
a.with.inset

## SUBFIGURE growth_rate_effect_estimated_versus_simulated_StDEvo_fixed_diversified.png ##
ggsave(file = "growth_rate_effect_estimated_versus_simulated_state_clocks_diversified.png", plot = a.with.inset, path = figures_dir, height = 7, width = 7)

b <- b + coord_cartesian(ylim = c(global_ymin, global_ymax), xlim = c(global_xmin, global_xmax ))
b_sub <- b_sub + border(linetype = "dashed", color = "darkgrey", size = 0.8)
b.with.inset <-
    ggdraw() +
    draw_plot(b)  +
    draw_plot(b_sub, x = 0.25, y = .65, width = .3, height = .3)
b.with.inset

ggsave(file = "growth_rate_effect_estimated_versus_simulated_state_clocks_random.png", plot = b.with.inset, path = figures_dir, height = 7, width = 7)


