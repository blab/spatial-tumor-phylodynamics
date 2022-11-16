#viz_extended_physicell_R2.R

library(tidyverse)
######################## 

##For reviews
# 2D_neut_bdg_usecarryingcapacity + greyed out 2D_neut_bdg
# 2D_sele_bdg_deleterious + greyed out 2D_sel_bdg

#2D neutral boundary-driven
growth_rate_posteriors_2D_net_bdg_n100 <- read_csv("../physicell/stats/growth_rate_posteriors_2D_neut_bdg_n100.csv") %>% 
    dplyr::rename("true_birth_diff" = true_birth_rate_diff)


growth_rate_posteriors_2D_neut_bdg_carrying <- posterior_summary %>% 
    dplyr::filter(simulation == "2D_neut_bdg_usecarryingcapacity", run != "sampconfig_m0_w1_d0.2_t1_mg1_mm1_l2e+08_i7_s95791_diversified_m1_n100")

baseplot_2d_neut_bdg <- growth_rate_posteriors_2D_net_bdg_n100 %>% 
    dplyr::filter(minBirthRateESS > 200) %>% 
    ggplot(., aes(x = true_birth_diff,
                  y = mean_birth_rate_diff), color = "darkgrey") +
    theme_classic() + geom_errorbar(aes(ymin = birthRate_hdi95_lower,
                                        ymax = birthRate_hdi95_upper),
                                    width = 0,
                                    alpha = 0.5,
                                    size = 1,  color = "darkgrey") +
    geom_point(aes(y = mean_birth_rate_diff), alpha = 0.8, size = 2, color = "darkgrey") +
    scale_color_viridis(option = "mako", discrete = TRUE, direction = -1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") + xlab("True birth rate difference (edge - center)") +
    ylab("Estimated birth rate difference (edge - center)") + theme(text=element_text(size=40))

plot_carryingcapacity <- baseplot_2d_neut_bdg + geom_errorbar(data = growth_rate_posteriors_2D_neut_bdg_carrying,
                                                              aes(ymin = birthRate_hdi95_lower,
                                                                  ymax = birthRate_hdi95_upper),
                                                              width = 0,
                                                              alpha = 0.5,
                                                              size = 1,  color = "darkred") +
    geom_point(data = growth_rate_posteriors_2D_neut_bdg_carrying,
               aes(y = mean_birth_rate_diff), alpha = 0.8, size = 2, color = "darkred")

plot_carryingcapacity
ggsave(plot_carryingcapacity,
       file = "../figures/physicell_2D_neut_carryingcapcity.png")

# 2D_sel_bdg_deleterious + grayed out 2D_sel_bdg

growth_rate_posteriors_2D_sel_bdg_deleterious <- posterior_summary %>% 
    filter(simulation == "2D_sel_bdg_deleterious", minBirthRateESS > 200) 


baseplot_2d_sel_bdg <- posteriors_2D_sel_bdg %>% 
    dplyr::filter(minBirthRateESS > 200) %>% 
    ggplot(., aes(x = true_birth_diff,
                  y = mean_birth_rate_diff), color = "darkgrey") +
    theme_classic() + geom_errorbar(aes(ymin = birthRate_hdi95_lower,
                                        ymax = birthRate_hdi95_upper),
                                    width = 0,
                                    alpha = 0.5,
                                    size = 1,  color = "darkgrey") +
    geom_point(aes(y = mean_birth_rate_diff), alpha = 0.8, size = 2, color = "darkgrey") +
    scale_color_viridis(option = "mako", discrete = TRUE, direction = -1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") + xlab("True birth rate difference (edge - center)") +
    ylab("Estimated birth rate difference (edge - center)") + theme(text=element_text(size=40))

plot_deleterious <- baseplot_2d_sel_bdg +
    geom_errorbar(data = growth_rate_posteriors_2D_sel_bdg_deleterious,
                  aes(ymin = birthRate_hdi95_lower,
                      ymax = birthRate_hdi95_upper),
                  width = 0,
                  alpha = 0.5,
                  size = 1,  color = "darkred") +
    geom_point(data = growth_rate_posteriors_2D_sel_bdg_deleterious,
               aes(y = mean_birth_rate_diff), alpha = 0.8, size = 2, color = "darkred")

plot_deleterious

ggsave(plot_deleterious,
       file = "../figures/physicell_2D_sel_deleterious.png", height = 5, width = 5)
