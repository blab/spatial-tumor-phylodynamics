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

#check overall convergence
# posterior_summary %>% 
#     filter(minBirthRateESS > 200) %>% 
# 
#     ggplot(., aes(x=true_birth_diff, y = mean_birth_rate_diff,  color = simulation)) +
#     geom_point() +
#     geom_errorbar(aes(ymin = birthRate_hdi95_lower, ymax = birthRate_hdi95_upper), width = 0) +
#     geom_abline(slope=1, intercept = 0, linetype = "dashed") +
#     theme_classic() + facet_wrap(~simulation)

#save the replacement
posteriors_2D_sel_bdg <- posterior_summary %>% 
    filter(minBirthRateESS > 200, simulation == "2D_sel_bdg")

write_csv(posteriors_2D_sel_bdg, file = "../physicell/stats/growth_rate_posteriors_2D_sel_bdg_n100_updated.csv")


#Also write everything
write_csv(posterior_summary, 
          file = "../physicell/stats/growth_rate_posteriors_extended_physicell.csv")

######################## 

##For reviews
# 2D_neut_bdg_usecarryingcapacity + greyed out 2D_neut_bdg
# 2D_sele_bdg_deleterious + greyed out 2D_sel_bdg

#2D neutral boundary-driven
growth_rate_posteriors_2D_net_bdg_n100 <- read_csv("../physicell/stats/growth_rate_posteriors_2D_neut_bdg_n100.csv") %>% 
    rename("true_birth_diff" = true_birth_rate_diff)


growth_rate_posteriors_2D_neut_bdg_carrying <- posterior_summary %>% 
    filter(simulation == "2D_neut_bdg_usecarryingcapacity", run != "sampconfig_m0_w1_d0.2_t1_mg1_mm1_l2e+08_i7_s95791_diversified_m1_n100")

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
# For supplement:
#     2D_neut_bdg_motility + greyed out 2D_neut_bdg

#Add per motolity 
library(stringr)
growth_rate_posteriors_2D_neut_bdg_motility <- posterior_summary %>% 
    dplyr::filter(simulation == "2D_neut_bdg_motility",
           minBirthRateESS > 200) %>% 
    mutate("motility" = as.numeric(str_extract(base, "(?<=_m).+(?=_w1)")))

plot_2D_neut_bdg_motility <-  baseplot_2d_neut_bdg +
    geom_errorbar(data = growth_rate_posteriors_2D_neut_bdg_motility,
                  aes(ymin = birthRate_hdi95_lower,
                      ymax = birthRate_hdi95_upper, 
                      color = motility),
                  width = 0,
                  alpha = 0.5,
                  size = 1) +
    geom_point(data = growth_rate_posteriors_2D_neut_bdg_motility,
               aes(y = mean_birth_rate_diff, color = motility), alpha = 0.8, size = 2) +
    #scale_color_gradient(low = "darkgrey", high = "darkred")
    scale_color_viridis(option="inferno", begin=0.2, end = 0.8) 


plot_2D_neut_bdg_motility 
plot_2D_neut_bdg_motility <- plot_2D_neut_bdg_motility + theme(legend.position = "none")
ggsave(plot_2D_neut_bdg_motility, file = "../figures/physicell_2D_neut_bdg_motility.png")
# 3D_sel_bgd + greyed out 2D_sel_bdg


growth_rate_posteriors_3D_sel_bdg <- posterior_summary %>% 
    dplyr::filter(simulation == "3D_sel",
                  minBirthRateESS > 200)

plot_3D_sel_bdg <-  baseplot_2d_sel_bdg +
    geom_errorbar(data = growth_rate_posteriors_3D_sel_bdg,
                  aes(ymin = birthRate_hdi95_lower,
                      ymax = birthRate_hdi95_upper),
                  width = 0,
                  alpha = 0.5,
                  size = 1,  color = "darkred") +
    geom_point(data = growth_rate_posteriors_3D_sel_bdg,
               aes(y = mean_birth_rate_diff), alpha = 0.8, size = 2, color = "darkred")
plot_3D_sel_bdg
ggsave(plot_3D_sel_bdg , file = "../figures/physicell_3D_sel.png")


#sigmoid growth 
# exclude d=0.8
growth_rate_posteriors_2D_neut_bdg_sigmoid<- posterior_summary %>% 
    dplyr::filter(simulation == "2D_neut_bdg_sigmoid",
                  minBirthRateESS > 200)

plot_2D_neut_bdg_sigmoid <-  baseplot_2d_neut_bdg +
    geom_errorbar(data = growth_rate_posteriors_2D_neut_bdg_sigmoid,
                  aes(ymin = birthRate_hdi95_lower,
                      ymax = birthRate_hdi95_upper),
                  width = 0,
                  alpha = 0.5,
                  size = 1,  color = "darkred") +
    geom_point(data = growth_rate_posteriors_2D_neut_bdg_sigmoid,
               aes(y = mean_birth_rate_diff), alpha = 0.8, size = 2, color = "darkred")
plot_2D_neut_bdg_sigmoid 

ggsave(plot_2D_neut_bdg_sigmoid,
       file = "../figures/physicell_2D_neut_bdg_sigmoid.png")


#make toy division prob curves
dat <- data.frame("x" = 1:6, "y" = c(1,1,2, 2, 2, 2))
p = ggplot(dat, aes(x = x, y = y)) + 
    geom_step(size = 2, color = "darkgrey")  +
    # theme_classic() +
    # theme(axis.text.x=element_blank(), #remove x axis labels
    #       axis.ticks.x=element_blank(), #remove x axis ticks
    #       axis.text.y=element_blank(),  #remove y axis labels
    #       axis.ticks.y=element_blank()  #remove y axis ticks
    #) +xlab("") + ylab("Division probability") +
    theme(text = element_text(size = 50)) + theme_void()   

ggsave(p,
       file = "../figures/toy_binary_pressure_division.png", height =3, width = 5)

logist <- function(x){
    #y = exp(x) / (1 + exp(x))
    y = 1/(1 + exp(-5 * (x - 1)))
}

p1 <- ggplot(data = data.frame(x = c(-2, 4)), aes(x))+ 
    stat_function(fun = logist, size = 2, color = "darkgrey") + theme_void()

p1

ggsave(p1,
       file = "../figures/toy_sigmoid_pressure_division.png", height =3, width = 5)
