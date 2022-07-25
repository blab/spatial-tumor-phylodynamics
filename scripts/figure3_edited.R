#figure3.R

##script to generate visualizations for figure 3

library(tumortree)
library(tidyverse)
library(cowplot)
library(viridis)
library(spatstat)
library(magick)

setwd("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/")
figures_dir <- "manuscript/figures"


sim_colors <- get_color_palette(c("boundary_driven", "unrestricted"))

#colors_edge_center <- get_color_palette(c("edge", "center"))
colors_loc <- c("edge" = "#89352F", 
                "center" = "#A2D2E2")
colors_edge_center <- colors_loc
tree_boundary_driven <- readRDS(file = "manuscript/analysis/simtrees/boundary_driven_timetree2.rds")
tree_molecular_boundary_driven <- readRDS(file = "manuscript/analysis/simtrees/boundary_driven_moleculartree2.rds")

sim_rates <- read_csv("outputs/simulation_analysis/sim_validation_rates.csv")

## SUBFIGURE growth_rate_effect_estimated_versus_simulated_StDEvo_fixed_diversified.png ##
fixed_clock_diversified_mcmc_logs <- read_mcmc_logs(path = "outputs/beast_analysis/state_dependent_clock_model/validation/logs",
                                                        pattern = "*matched_diversified_fixed\\.log")
fixed_clock_diversified_death_rates <- purrr::map_dbl(names(fixed_clock_diversified_mcmc_logs), extract_death_rate)


growth_rate_fixed_diversified_df <- purrr::map2(fixed_clock_diversified_mcmc_logs,
                                                fixed_clock_diversified_death_rates,
                      function(x, dr) calc_state_dependent_growthRateDiff_means(mcmc_obj = x,
                                                            dr = dr,
                                                            fixed_dr = 0.1)) %>%
    bind_rows %>%
    left_join(sim_rates, by = "dr") %>%
    dplyr::mutate("mean_growth_rate_diff" = mean_edge_growth_rate - mean_center_growth_rate)

b <- ggplot(growth_rate_fixed_diversified_df  , aes(x = mean_growth_rate_diff, y = mean), color = "grey") +
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

## SUBFIGURE ##
## read in example tumor
all_cells_boundary_driven <- read_csv("outputs/raw_simulation_results/validation/cells_death_rate_validation_pop_1000_dr_0.10.csv") %>%
    normalize_locs

alive_cells_boundary_driven <- all_cells_boundary_drive %>%
    filter_alive_cells
set.seed(81721)
sampled_cells_boundary_driven_diversified <- sample_alive_cells(alive_cells = alive_cells_boundary_driven,
                                                    n = 100, diversified_sampling = TRUE)
b_sub <- ggplot(sampled_cells_boundary_driven_diversified, aes(x = locx, y = locy)) +
    geom_point(aes(color = sampled, size = as.integer(sampled))) +
    theme_void() + scale_color_manual(values = c("grey", "#809070")) +
    scale_size(range = c(1,2)) + theme(legend.position = "none")
b_sub



## SUBFIGURE growth_rate_effect_estimated_versus_simulated_StDEvo_random.png
fixed_clock_random_mcmc_logs <- read_mcmc_logs(path = "outputs/beast_analysis/state_dependent_clock_model/validation/logs",
                                                    pattern = "*matched_fixed\\.log")
fixed_clock_random_death_rates <- purrr::map_dbl(names(fixed_clock_random_mcmc_logs), extract_death_rate)


growth_rate_fixed_random_df <- purrr::map2(fixed_clock_random_mcmc_logs,
                                           fixed_clock_random_death_rates,
                                                function(x, dr) calc_state_dependent_growthRateDiff_means(mcmc_obj = x,
                                                                                                          dr = dr,
                                                                                                          fixed_dr = 0.1)) %>%
    bind_rows %>%
    left_join(sim_rates, by = "dr") %>%
    dplyr::mutate("mean_growth_rate_diff" = mean_edge_growth_rate - mean_center_growth_rate)

c <- ggplot(growth_rate_fixed_random_df  , aes(x = mean_growth_rate_diff, y = mean), color = "grey") +
    geom_errorbar(alpha = 0.7, aes(ymin = hdi95_lower, ymax = hdi95_upper),  size=1, width=0) +
    geom_errorbar(alpha = 0.4, aes(ymin = hdi85_lower, ymax = hdi85_upper),  size=2, width=0) +
    geom_errorbar(alpha = 0.4, aes(ymin = hdi75_lower, ymax = hdi75_upper),  size=3, width=0) +
    geom_point(size = 2) +
    #geom_hline(yintercept = 0, linetype = "dashed") +
    #geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0) +
    theme_classic() + labs("color" = "") +
    ylab("Estimated state-dependent effect on growth rate") + xlab("True growth rate difference")
c

set.seed(3333)
sampled_cells_boundary_driven_random <- sample_alive_cells(alive_cells = alive_cells_boundary_driven,
                                                    n = 100, diversified_sampling = FALSE)

c_sub <- ggplot(sampled_cells_boundary_driven_random, aes(x = locx, y = locy)) +
    geom_point(aes(color = sampled, size = as.integer(sampled))) +
    theme_void() + scale_color_manual(values = c("grey", "#809070")) +
    scale_size(range = c(1,2)) + theme(legend.position = "none")
c_sub
## SUBFIGURE ##


## INSETS WITH COMMOON AXES


global_ymin <- min(c(growth_rate_fixed_diversified_df$hdi95_lower,growth_rate_fixed_random_df$hdi95_upper))
global_ymax <- max(c(growth_rate_fixed_diversified_df$hdi95_lower,growth_rate_fixed_random_df$hdi95_upper)) + 0.4

global_xmin <- min(c(growth_rate_fixed_diversified_df$mean_growth_rate_diff))
global_xmax <- max(c(growth_rate_fixed_diversified_df$mean_growth_rate_diff))



b <- b + coord_cartesian(ylim = c(global_ymin, global_ymax), xlim = c(global_xmin, global_xmax ))
b_sub <- b_sub + border(linetype = "dashed", color = "darkgrey", size = 0.8)
b.with.inset <-
    ggdraw() +
    draw_plot(b)  +
    draw_plot(b_sub, x = 0.65, y = .65, width = .3, height = .3)
b.with.inset

## SUBFIGURE growth_rate_effect_estimated_versus_simulated_StDEvo_fixed_diversified.png ##
ggsave(file = "growth_rate_effect_estimated_versus_simulated_fixed_diversified.png", plot = b.with.inset, path = figures_dir, height = 7, width = 7)

c <- c + coord_cartesian(ylim = c(global_ymin, global_ymax), xlim = c(global_xmin, global_xmax ))
c_sub <- c_sub + border(linetype = "dashed", color = "darkgrey", size = 0.8)
c.with.inset <-
    ggdraw() +
    draw_plot(c)  +
    draw_plot(c_sub, x = 0.65, y = .65, width = .3, height = .3)
c.with.inset

ggsave(file = "growth_rate_effect_estimated_versus_simulated_fixed_random.png", plot = c.with.inset, path = figures_dir, height = 7, width = 7)

############# PHYSICELL 2D NEUTRAL BOUNDARY-DRIVEN GROWTH #####
growth_rate_posteriors_2D_net_bdg_n100 <- read_csv("~/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/manuscript/analysis/beast_results/growth_rate_posteriors_2D_neut_bdg_n100.csv")
d_baseplot <- growth_rate_posteriors_2D_net_bdg_n100 %>% 
    dplyr::filter(minBirthRateESS > 150, dr != 0) %>% 
    ggplot(., aes(x = true_birth_rate_diff, y = mean_birth_rate_diff), color = "black") +
    theme_classic() + geom_errorbar(aes(ymin = birthRate_hdi95_lower, ymax = birthRate_hdi95_upper), width = 0, alpha = 0.5, size = 1) +
    geom_point(aes( y = mean_birth_rate_diff), alpha = 0.8, size = 2) + scale_color_viridis(option = "mako", discrete = TRUE, direction = -1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") + xlab("True growth rate difference (edge - center)") +
    ylab("True growth rate difference (edge - center)") +
    ggtitle("Physicell 2D Neutral Boundary-Driven Growth") #+
    #labs(color = "P(cell death)") 

#labs(color = "P(cell death)") 
d_baseplot
d_lowdeath <- ggdraw() +
    draw_image("manuscript/figures/sampconfig_m0_w1_d0.1_t1_mg1_mm1_l2e+08_i2_s42585.png") +
    theme_void()

d_middeath <- ggdraw() +
    draw_image("manuscript/figures/sampconfig_m0_w1_d0.3_t1_mg1_mm1_l2e+08_i3_s31628.png") +
    theme_void()

d_highdeath <- ggdraw() +
    draw_image("manuscript/figures/sampconfig_m0_w1_d0.7_t1_mg1_mm1_l2e+08_i3_s73619.png") +
    theme_void()

tumor_axis <- plot_grid(
    d_highdeath, d_middeath,d_lowdeath, nrow = 1)

d_full <- plot_grid(d_baseplot, tumor_axis, ncol = 1, rel_heights = c(1,0.2), rel_widths = c(1, 0.6))
d_full


ggsave(plot=d_full, file ="manuscript/figures/physicell_2D_neutral_bdg.png", height = 5, width = 4.5)

############# PHYSICELL 3D NEUTRAL BOUNDARY-DRIVEN GROWTH #####
growth_rate_posteriors_3D_net_bdg_n100 <- read_csv("~/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/manuscript/analysis/beast_results/growth_rate_posteriors_3D_neut_bdg_n100_single_death.csv")
e_baseplot <- growth_rate_posteriors_3D_net_bdg_n100 %>% 
    dplyr::filter(minBirthRateESS > 150) %>% 
    ggplot(., aes(x = true_growth_rate_diff, y = mean_growth_rate_diff), color = "black") +
    theme_classic() + geom_errorbar(aes(ymin = growthRate_hdi95_lower, ymax = growthRate_hdi95_upper), width = 0, alpha = 0.5, size = 1) +
    geom_point(aes( y = mean_growth_rate_diff), alpha = 0.8, size = 2) + scale_color_viridis(option = "mako", discrete = TRUE, direction = -1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") + xlab("True growth rate difference (edge - center)") +
    ylab("True growth rate difference (edge - center)") + ggtitle("Physicell 3D Neutral Boundary-Driven Growth") +
    ylim(-0.001, 0.004) ##TODO: take this away when we figure out what's going with sims
#labs(color = "P(cell death)") 
e_baseplot
e_lowdeath <- ggdraw() +
    draw_image("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/pressure_thumbnails_2/3D_neut_bdg/pressure_plots_recol/sampconfig_m0_w1_d0.1_t1_mg1_mm1_l2e+08_i7_s91196.png") +
    theme_void()

e_middeath <- ggdraw() +
    draw_image("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/pressure_thumbnails_2/3D_neut_bdg/pressure_plots_recol/sampconfig_m0_w1_d0.4_t1_mg1_mm1_l2e+08_i2_s8484.png") +
    theme_void()

e_highdeath <- ggdraw() +
    draw_image("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/pressure_thumbnails_2/3D_neut_bdg/pressure_plots_recol/sampconfig_m0_w1_d0.7_t1_mg1_mm1_l2e+08_i7_s28098.png") +
    theme_void()



e_tumor_axis <- plot_grid(
    e_highdeath, e_middeath,e_lowdeath, nrow = 1)

e_full <- plot_grid(e_baseplot, e_tumor_axis, ncol = 1, rel_heights = c(1,0.2), rel_widths = c(1, 0.6))
e_full


ggsave(plot=e_full, file ="manuscript/figures/physicell_3D_neutral_bdg.png", height = 5, width = 4.5)

############# PHYSICELL 2D SELECTION BOUNDARY-DRIVEN GROWTH #####
growth_rate_posteriors_2D_sel_bdg_n100 <- read_csv("~/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/manuscript/analysis/beast_results/growth_rate_posteriors_2D_sel_bdg_n100_single_death.csv")
f_baseplot <- growth_rate_posteriors_2D_sel_bdg_n100 %>% 
    dplyr::filter(minBirthRateESS > 150) %>% 
    ggplot(., aes(x = true_growth_rate_diff, y = mean_growth_rate_diff), color = "black") +
    theme_classic() + geom_errorbar(aes(ymin = growthRate_hdi95_lower, ymax = growthRate_hdi95_upper), width = 0, alpha = 0.5, size = 1) +
    geom_point(aes( y = mean_growth_rate_diff), alpha = 0.8, size = 2) + scale_color_viridis(option = "mako", discrete = TRUE, direction = -1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") + xlab("True growth rate difference (edge - center)") +
    ylab("True growth rate difference (edge - center)") + ggtitle("Physicell 2D Selection Boundary-Driven Growth") #+

f_baseplot

sel_bdg_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/logs/2D_sel_bdg_highermu/diversified_100", 
                        pattern = "single_death.log")

f_highest_diff_sim <- paste0("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/pressure_thumbnails_2/2D_sel_bdg_highermu/pressure_plots_recol/",
                           gsub("_diversified_m1_n100_single_death", "", gsub(".log", ".png", sel_bdg_files[which.max(growth_rate_posteriors_2D_sel_bdg_n100$mean_growth_rate_diff)])), sep = "")

f_lowest_diff_sim <- paste0("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/pressure_thumbnails_2/2D_sel_bdg_highermu/pressure_plots_recol/",
                          gsub("_diversified_m1_n100_single_death", "", gsub(".log", ".png", sel_bdg_files[which.min(growth_rate_posteriors_2D_sel_bdg_n100$mean_growth_rate_diff)])), sep = "")

f_median_sim <- sort(growth_rate_posteriors_2D_sel_bdg_n100$mean_growth_rate_diff)[round(length(growth_rate_posteriors_2D_sel_bdg_n100$mean_growth_rate_diff)/2)]

f_middle_diff_sim <- paste0("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/pressure_thumbnails_2/2D_sel_bdg_highermu/pressure_plots_recol/",
                          gsub("_diversified_m1_n100_single_death", "", gsub(".log", ".png", sel_bdg_files[which(growth_rate_posteriors_2D_sel_bdg_n100$mean_growth_rate_diff == f_median_sim)])), sep = "")


f_low <- ggdraw() +
    draw_image(f_lowest_diff_sim) +
    theme_void()

f_mid<- ggdraw() +
    draw_image(f_middle_diff_sim) +
    theme_void()

f_high <- ggdraw() +
    draw_image(f_highest_diff_sim) +
    theme_void()

f_high


f_tumor_axis <- plot_grid(
    f_low, f_mid,f_high, nrow = 1)

f_full <- plot_grid(f_baseplot, f_tumor_axis, ncol = 1, rel_heights = c(1,0.2), rel_widths = c(1, 0.6))
f_full


ggsave(plot=f_full, file ="manuscript/figures/physicell_2D_sel_bdg.png", height = 5, width = 4.5)

############# PHYSICELL 2D SELECTION UNRESTRICTED GROWTH #####
growth_rate_posteriors_2D_sel_n100 <- read_csv("~/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/manuscript/analysis/beast_results/growth_rate_posteriors_2D_sel_n100_single_death.csv")
g_baseplot <- growth_rate_posteriors_2D_sel_n100 %>% 
    dplyr::filter(minBirthRateESS > 150) %>% 
    ggplot(., aes(x = mean_growth_rate_diff), color = "black") +
    theme_classic() + geom_histogram(fill = "grey", color = "black") +
    xlab("Estimated growth rate difference") + ylab("Simulation count") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggtitle("2D Selection Unrestricted Growth")

g_baseplot

# growth_rate_posteriors_2D_sel_n100[which.max(growth_rate_posteriors_2D_sel_n100$mean_growth_rate_diff),]$i
# growth_rate_posteriors_2D_sel_n100[which.max(growth_rate_posteriors_2D_sel_n100$mean_growth_rate_diff),]$dr
sel_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/logs/2D_sel_highermu/diversified_100", 
           pattern = "single_death.log")

highest_diff_sim <- paste0("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/pressure_thumbnails_2/2D_sel_highermu/pressure_plots_recol/",
       gsub("_diversified_m1_n100_single_death", "", gsub(".log", ".png", sel_files[which.max(growth_rate_posteriors_2D_sel_n100$mean_growth_rate_diff)])), sep = "")

lowest_diff_sim <- paste0("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/pressure_thumbnails_2/2D_sel_highermu/pressure_plots_recol/",
                           gsub("_diversified_m1_n100_single_death", "", gsub(".log", ".png", sel_files[which.min(growth_rate_posteriors_2D_sel_n100$mean_growth_rate_diff)])), sep = "")

median_sim <- sort(growth_rate_posteriors_2D_sel_n100$mean_growth_rate_diff)[round(length(growth_rate_posteriors_2D_sel_n100$mean_growth_rate_diff)/2)]
middle_diff_sim <- paste0("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/pressure_thumbnails_2/2D_sel_highermu/pressure_plots_recol/",
                          gsub("_diversified_m1_n100_single_death", "", gsub(".log", ".png", sel_files[which(growth_rate_posteriors_2D_sel_n100$mean_growth_rate_diff == median_sim)])), sep = "")


g_highest_diff_sim <- ggdraw() +
    draw_image(highest_diff_sim) +
    theme_void()

g_lowest_diff_sim <- ggdraw() +
    draw_image(lowest_diff_sim) +
    theme_void()

g_middle_diff_sim <- ggdraw() +
    draw_image(middle_diff_sim ) +
    theme_void()

g_tumor_axis <- plot_grid(
    g_lowest_diff_sim, g_middle_diff_sim, g_highest_diff_sim, nrow = 1)

g_full <- plot_grid(g_baseplot, g_tumor_axis, ncol = 1, rel_heights = c(1,0.2), rel_widths = c(1, 0.6))
g_full

# g.with.inset <-
#     ggdraw() +
#     draw_plot(g_baseplot)  +
#     draw_plot(g_highest_diff_sim, x = 0.8, y = 0.3, width = 0.15, height = 0.15) +
#     draw_plot(g_lowest_diff_sim, x = 0.1, y = 0.3, width = 0.15, height = 0.15)
# g.with.inset
# d.with.inset

ggsave(plot=g_full, file ="manuscript/figures/physicell_2D_selection_hist.png", height = 5, width = 4.5)
############# EDEN BULK  #####
true_eden_rates <- read_csv("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/validation_growth_and_death_rates_weighted.csv") %>% 
    dplyr::rename(true_mean_edge_growth_rate = mean_edge_growth_rate,
                  true_mean_center_growth_rate = mean_center_growth_rate, 
                  true_mean_death_rate = mean_death_rate) %>% 
    dplyr::mutate("true_mean_growth_rate_diff" = true_mean_edge_growth_rate - true_mean_center_growth_rate)

growth_rate_posteriors_bulk_eden <- read_csv("~/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/manuscript/analysis/beast_results/growth_rate_posteriors_bulk_single_death.csv") %>% 
    left_join(., true_eden_rates, by = "dr")

h_baseplot <- growth_rate_posteriors_bulk_eden%>% 
    dplyr::filter(minBirthRateESS > 200) %>% 
    ggplot(., aes(x = true_mean_growth_rate_diff, y = mean_growth_rate_diff), color = "white") +
    theme_classic() + geom_errorbar(color = "white", aes(ymin = growthRate_hdi95_lower, ymax = growthRate_hdi95_upper), width = 0, alpha = 0.5, size = 1) +
    geom_point(color = "white", aes( y = mean_growth_rate_diff), alpha = 0.8, size = 2) + scale_color_viridis(option = "mako", discrete = TRUE, direction = -1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") + xlab("True growth rate difference (edge - center)") +
    ylab("True growth rate difference (edge - center)") + #ggtitle("Bulk sampling") + 
    ylim(0, 0.4)

h_baseplot

# example_bulk_cells <- read_csv("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/manuscript/analysis/example_bulk_sequencing_alive_cells.csv")
# bulk_tumor_inset <- ggplot(example_bulk_cell, aes(x = locx, y = locy, color = ifelse(punched & punch_on_edge, "edge_punch",
#                                                                 ifelse(punched, "center_punch", "not_sampled")))) +
#     geom_point(size = size) +
#     theme_void() +
#     scale_color_manual(values = c("not_sampled" = "darkgrey",
#                                   "edge_punch" = scales::muted('red'),
#                                   "center_punch" = "black")) +
#     theme(legend.position = "none")

#bulk_tumor_inset <- bulk_tumor_inset + border(linetype = "dashed", color = "darkgrey", size = 0.8)

bulk_tumor_inset <- ggdraw() +
    draw_image("manuscript/figures/example_bulk_sampling.png") +
    theme_void() +
    theme(panel.border = element_rect(colour = "darkgrey", fill=NA, size=0.5, linetype = "dashed"))

h.with.inset <-
    ggdraw() +
    draw_plot(h_baseplot)  +
    draw_plot(bulk_tumor_inset, x = 0.25, y = .65, width = .25, height = .25)
h.with.inset
ggsave(plot=h.with.inset, file ="manuscript/figures/bulk_true_versus_estimated_with_inset.png", height = 5, width = 4.5)


################ EXAMPLE POSTERIOR ##################

example_mcmc_log_file <- "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs2/death_rate_validation_pop_1000_dr_0.29_n_100_state_clock_estimate_dr.log"
example_mcmc_log <- as.data.frame(readLog(example_mcmc_log_file))
example_birth_rates_diff_posteriors <- example_mcmc_log %>% 
    dplyr::mutate(birthRateDiff = birthRateCanonical.1 -  birthRateCanonical.0) %>% 
    ggplot(aes(x = birthRateDiff)) + geom_density(fill = "darkgrey", alpha = 0.8) +
    theme_classic() + xlab("Estimated birth rate difference") + ylab("") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y=element_blank()) +
    theme(text = element_text(size = 20))

example_birth_rates_posteriors <- example_mcmc_log %>% 
    tidyr::pivot_longer(cols = c("birthRateCanonical.1", "birthRateCanonical.0"), 
                        names_prefix = "birthRateCanonical.",
                        values_to = "birthRate", 
                        names_to = "state") %>% 
    
    ggplot(aes(x = birthRate, fill = ifelse(state == 1, "edge", "center"))) + geom_density(alpha = 0.8, aes(fill = ifelse(state == 1, "edge", "center"))) +
    theme_classic() + xlab("Estimated birth rate") + ylab("") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y=element_blank()) +
    scale_fill_manual(values = colors_edge_center) +
    theme(legend.position = "none") +
    theme(text = element_text(size = 20))


 
ggsave(plot=example_birth_rates_diff_posteriors, file ="manuscript/figures/example_birth_rates_diff_posteriors.png", height = 5, width = 4.5)

ggsave(plot=example_birth_rates_posteriors, file ="manuscript/figures/example_birth_rates_posteriors.png", height = 5, width = 4.5)

############## EXAMPLE ANCESTRAL STATE RECONSTRUCTION TREE #####

#mcc_tree <- read.beast("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs2/death_rate_validation_pop_1000_dr_0.31_n_50_state_clock_estimate_single_dr.typed.node.mcc.tree")
mcc_tree <- read.beast("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs2/death_rate_validation_pop_1000_dr_0.29_n_100_state_clock_estimate_dr.typed.node.trees.mcc")

#colors_loc <- c("edge" = "#75ACBA", "center" = "#A28F7B")
colors_loc <- colors_edge_center
names(colors_loc) <- c("loc1", "loc0")

anc_state_nodes <- mcc_tree@data %>% 
    mutate(type.prob = as.numeric(type.prob)) %>% 
    dplyr::mutate(edge = type.prob * as.integer(type == "loc1") + (1 - type.prob)*as.integer(type == "loc0")) %>% 
    dplyr::mutate(center = type.prob * as.integer(type == "loc0") + (1 - type.prob)*as.integer(type == "loc1")) %>% 
    dplyr::select(center, edge, type, node)

#colors_edge_center <- c("edge" = "#75ACBA", "center" = "#A28F7B")
pies <- nodepie(anc_state_nodes, cols=1:2, alpha=0.8, color = colors_edge_center)
g <- ggtree(mcc_tree, color = "darkgrey") +
    #geom_point(size = 2) +
    scale_color_manual(values = colors_loc) +
    theme(legend.position = "none")
g_pie <- ggtree::inset(g, pies, width = 0.066, height = 0.066) + theme(legend.position = "none")

# ggtree(mcc_tree, aes(color = type)) +
#     #geom_point(size = 2) +
#     scale_color_manual(values = colors_loc) +
#     theme(legend.position = "none")
# g2 <- ggtree(mcc_tree, color = "darkgrey") +
#     #geom_point(size = 2) +
#     scale_color_manual(values = colors_loc) +
#     theme(legend.position = "none")
# g2_pie <- ggtree::inset(g2, pies, width = 0.07, height = 0.07) + theme(legend.position = "none")
# g2_pie_toy <- viewClade(g2_pie, MRCA(g,"cell5338loc1", "cell5378loc1"))

#ggsave(plot=g2_pie_toy, file ="manuscript/figures/toy_ancestral_state_recon_tree.png", height = 5, width = 5)

ggsave(plot=g_pie, file ="manuscript/figures/example_ancestral_state_recon_tree.png", height = 5, width = 5)
sim_colors <- c("boundary_driven" = "#B6465F", "unrestricted" = "#3C3C3C")

# data=data.frame("x" = rexp(1000, rate = 0.1))
# 
# ggplot(data, aes(x)) + geom_density(fill = "#B6465F") + theme_classic()  +
#     theme(axis.title.y=element_blank(),
#           axis.text.y=element_blank(),
#           axis.ticks.y=element_blank(),
#           axis.line.y=element_blank()) +
#     xlab("Terminal branch length") 

# data2=data.frame("x" = rexp(1000, rate = 10))
# ggplot(data2, aes(x)) + geom_density(fill = "#3C3C3C") + theme_classic()  +
#     theme(axis.title.y=element_blank(),
#           axis.text.y=element_blank(),
#           axis.ticks.y=element_blank(),
#           axis.line.y=element_blank()) +
#     xlab("Terminal branch length") 

# example_birth_rates_diff_posteriors <- example_mcmc_log %>% 
#     dplyr::mutate(birthRateDiff = birthRateCanonical.1 -  birthRateCanonical.0) %>% 
#     ggplot(aes(x = birthRateDiff)) + geom_density(fill = "darkgrey", alpha = 0.8) +
#     theme_classic() + xlab("Estimated birth rate difference") + ylab("") +
#     theme(axis.title.y=element_blank(),
#           axis.text.y=element_blank(),
#           axis.ticks.y=element_blank(),
#           axis.line.y=element_blank())
## example tumor

example_all_cells <- read_csv("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation/cells_death_rate_validation_pop_1000_dr_0.29.csv")

### Ancestral state reconstruction ######

sim_tree <- readRDS("~/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/manuscript/analysis/simtrees/cells_death_rate_validation_pop_1000_dr_0.29.rds")

example_alive_cells <- example_all_cells %>% 
    filter_alive_cells()
example_alive_cells <- example_alive_cells %>% 
    mark_boundary(alive_cells = example_alive_cells)

example_alive_cells$sampled <- example_alive_cells$index %in% as.numeric(gsub("loc[01]", "", gsub("cell","",  mcc_tree@phylo$tip.label)))

example_sampled_cells <- example_alive_cells %>% 
    filter(sampled)
example_tumor_inset <- ggplot(example_alive_cells, aes(locx, locy, color = ifelse(est_edge == 1, "edge", "center"))) +
    geom_point() + theme_void() +
    geom_point(data = example_sampled_cells, shape = 21, aes(fill = ifelse(est_edge == 1, "edge", "center")), color = "black", size =2)+
    scale_fill_manual(values = colors_edge_center) +
    scale_color_manual(values = colors_edge_center) +
    theme(legend.position = "none") 
example_tumor_inset
ggsave(plot=example_tumor_inset, file ="manuscript/figures/example_tumor_edge_center_sampled.png", height = 5, width = 5)

match_beast_to_sim_tree <- function(sim_tree, beast_tree) {
    

    #get leaves for both trees
    beast_tips <- purrr::map_chr(beast_tree@phylo$tip.label,
                                 function(s) str_extract(string = s, pattern = "[0-9]+"))
    sim_tips <- purrr::map_chr(sim_tree@phylo$tip.label,
                               function(s) str_extract(string = s, pattern = "[0-9]+"))
    
    #find internal nodes and children
    sim_internal_nodes <- length(sim_tips)+1:Nnode(sim_tree@phylo)
    node_child_list <- purrr::map(sim_internal_nodes, function(n) Descendants(sim_tree@phylo, node = n, type = "tips"))
    
    #find corresponding beast nodes
    beast_node_vec <- c()
    for (i in 1:length(node_child_list)) {
        
        child_vec <- sim_tips[unlist(node_child_list[[i]])]
        
        child_nodes <- which(beast_tips %in% as.character(child_vec))
        
        parent_node <- phytools::findMRCA(beast_tree@phylo, tips=child_nodes)
        
        
        
        beast_node_vec <- c(beast_node_vec, parent_node)
    }
    
    
    n <- length(sim_tree@phylo$tip.label)
    m <- n + sim_tree@phylo$Nnode
    
    node_conversion_df <- data.frame("sim_node" = sim_internal_nodes,
                                     "beast_node" = beast_node_vec,
                                     "true_state" = sim_tree@data$state[sim_internal_nodes], 
                                     "frac_time_on_edge" = sim_tree@data$frac_time_on_edge[sim_internal_nodes],
                                     "predicted_state" = beast_tree@data$type[match(beast_node_vec, beast_tree@data$node)],
                                     "posterior_prob" = beast_tree@data$type.prob[match(beast_node_vec, beast_tree@data$node)]) %>% 
        mutate("sim_height" = node.height(sim_tree@phylo)[(n+1):m],
               "beast_height" =  node.height(beast_tree@phylo)[(n+1):m]) %>% 
        dplyr::mutate("is.edge" = as.integer(predicted_state =="loc1")) %>% 
        dplyr::mutate("edge_prob" = (1-is.edge)*(1-as.numeric(posterior_prob)) + is.edge*as.numeric(posterior_prob))

    return(node_conversion_df)
    
}

sim_tree <- prune_simulated_tree(sim_tree, sampled_cells = example_sampled_cells$index)
node_conversion_df <- match_beast_to_sim_tree(sim_tree = sim_tree, beast_tree = mcc_tree)
logistic_model <- glm(true_state ~ edge_prob, node_conversion_df, family = binomial)
model_fit <- data.frame("edge_prob" = seq(0, 1, by = 0.01))
model_fit$true_state <- predict(logistic_model, newdata=model_fit, type="response")

# Plot the modeled probability values

true_versus_post_prob <- ggplot(node_conversion_df, aes(edge_prob, true_state, color = ifelse(true_state == 1, "edge", "center"))) +
    geom_jitter(width = 0, height = 0.1) + theme_classic() + scale_color_manual(values = colors_edge_center) +
    theme(legend.position = "none") + xlab("Posterior probability edge") + ylab("True state") + xlim(0, 1)+
    geom_line(data=model_fit,color = "black") +
    theme(text = element_text(size = 20))

true_versus_post_prob
ggsave(plot=true_versus_post_prob, file ="manuscript/figures/example_true_versus_post_prob_ancestor_recon.png", height = 5, width = 5)

