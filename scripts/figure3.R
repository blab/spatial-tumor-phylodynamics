#figure3.R
#Script to generate plot for Figure 3

library(tumortree)
library(tidyverse)
library(HDInterval)
library(coda)
library(beastio)
library(treeio)
library(ggtree)
library(phangorn)

colors_edge_center <- get_color_palette(names = c("edge", "center"))

################ EXAMPLE POSTERIORS ##################

#example_mcmc_log_file <- "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs2/death_rate_validation_pop_1000_dr_0.29_n_100_state_clock_estimate_dr.log"
example_mcmc_log_file <- "../eden/logs/death_rate_validation_pop_1000_dr_0.29_n_100_state_clock_estimate_dr.log"

example_mcmc_log <- as.data.frame(readLog(example_mcmc_log_file))

# Get true simulated rates
#From extract_validation_sims_rates.R
true_example_diff_df <- read_csv("../eden/stats/validation_growth_and_death_rates_weighted.csv") %>% 
    dplyr::mutate(true_birth_rate_diff_weighted = mean_edge_birth_rate - mean_center_birth_rate) %>% 
    dplyr::select(dr, true_birth_rate_diff_weighted) %>% 
    dplyr::filter(dr == "0.29")

#Get HPD intervals

example_birth_rates_diff_posteriors_summary <- example_mcmc_log %>% 
    dplyr::mutate(birthRateDiff = birthRateCanonical.1 -  birthRateCanonical.0) %>% 
    dplyr::summarise(mean = mean(birthRateDiff),
                  hpd_90_lower = hdi(birthRateDiff,
                                     credMass =0.9)[1],
                  hpd_90_upper = hdi(birthRateDiff,
                                     credMass =0.9)[2])
# Example edge - center difference posterior distibution
example_birth_rates_diff_posteriors <- example_mcmc_log %>% 
    dplyr::mutate(birthRateDiff = birthRateCanonical.1 -  birthRateCanonical.0) %>% 
    ggplot(aes(x = birthRateDiff)) + geom_density(fill = "darkgrey", alpha = 0.8) +
    theme_classic() + xlab("Estimated birth rate difference") + ylab("") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y=element_blank()) +
    theme(text = element_text(size = 20)) +
    geom_vline(xintercept = true_example_diff_df$true_birth_rate_diff_weighted[1], linetype = "dashed")

# Example edge - center posterior distibutions
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


example_birth_rates_diff_posteriors
example_birth_rates_posteriors

ggsave(plot=example_birth_rates_diff_posteriors, file ="../figures/example_birth_rates_diff_posteriors.png", height = 5, width = 4.5)

ggsave(plot=example_birth_rates_posteriors, file ="../figures/example_birth_rates_posteriors.png", height = 5, width = 4.5)

############## EXAMPLE ANCESTRAL STATE RECONSTRUCTION TREE #####

#mcc_tree <- read.beast("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs2/death_rate_validation_pop_1000_dr_0.29_n_100_state_clock_estimate_dr.typed.node.trees.mcc")
mcc_tree <- read.beast("../eden/trees/death_rate_validation_pop_1000_dr_0.29_n_100_state_clock_estimate_dr.typed.node.trees.mcc")


colors_loc <- colors_edge_center
names(colors_loc) <- c("loc1", "loc0")

anc_state_nodes <- mcc_tree@data %>% 
    mutate(type.prob = as.numeric(type.prob)) %>% 
    dplyr::mutate(edge = type.prob * as.integer(type == "loc1") + (1 - type.prob)*as.integer(type == "loc0")) %>% 
    dplyr::mutate(center = type.prob * as.integer(type == "loc0") + (1 - type.prob)*as.integer(type == "loc1")) %>% 
    dplyr::select(center, edge, type, node)


pies <- nodepie(anc_state_nodes, cols=1:2, alpha=0.8, color = colors_edge_center)
g <- ggtree(mcc_tree, color = "darkgrey") +
    #geom_point(size = 2) +
    scale_color_manual(values = colors_loc) +
    theme(legend.position = "none")
g_pie <- ggtree::inset(g, pies, width = 0.066, height = 0.066) + theme(legend.position = "none")

## To visualize tip labels
# ggtree(mcc_tree, aes(color = type)) +
#     #geom_point(size = 2) +
#     scale_color_manual(values = colors_loc) +
#     theme(legend.position = "none") + geom_tiplab()

## Toy tree for model schematic
pies2 <- nodepie(anc_state_nodes, cols=1:2, alpha=1, color = colors_edge_center)

g2 <- ggtree(mcc_tree, color = "darkgrey", size = 2) +
    geom_point(size = 2) +
    scale_color_manual(values = colors_loc) +
    theme(legend.position = "none") 
g2_pie <- ggtree::inset(g2, pies2, width = 0.08, height = 0.08) + theme(legend.position = "none")

g2_pie_toy <- viewClade(g2_pie, MRCA(g2,"cell4628loc1", "cell4232loc0"))



ggsave(plot=g2_pie_toy,
       file ="../figures/toy_ancestral_state_recon_tree.png", height = 5, width = 5)

ggsave(plot=g_pie, file ="../figures/example_ancestral_state_recon_tree.png", height = 5, width = 5)


#example_all_cells <- read_csv("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation/cells_death_rate_validation_pop_1000_dr_0.29.csv")
example_all_cells <- read_csv("../simulation_data/cells_death_rate_validation_pop_1000_dr_0.29.csv")

### Ancestral state reconstruction ######

sim_tree <- readRDS("../eden/simtrees/cells_death_rate_validation_pop_1000_dr_0.29.rds")

example_alive_cells <- example_all_cells %>% 
    filter_alive_cells()

example_alive_cells <- example_alive_cells %>% 
    mark_boundary(alive_cells = example_alive_cells)

example_alive_cells$sampled <- example_alive_cells$index %in% as.numeric(gsub("loc[01]", "", gsub("cell","",  mcc_tree@phylo$tip.label)))

write_csv(example_alive_cells, "../simulation_data/alive_cells_death_rate_validation_pop_1000_dr_0.29.csv")

# Can skip to here
example_alive_cells <- read_csv(file="../simulation_data/alive_cells_death_rate_validation_pop_1000_dr_0.29.csv")

example_sampled_cells <- example_alive_cells %>% 
    filter(sampled)

example_tumor_inset <- ggplot(example_alive_cells, aes(locx, locy, color = ifelse(est_edge == 1, "edge", "center"))) +
    geom_point() + theme_void() +
    geom_point(data = example_sampled_cells, shape = 21, aes(fill = ifelse(est_edge == 1, "edge", "center")), color = "black", size =2)+
    scale_fill_manual(values = colors_edge_center) +
    scale_color_manual(values = colors_edge_center) +
    theme(legend.position = "none") 

example_tumor_inset
ggsave(plot=example_tumor_inset, file ="../figures/example_tumor_edge_center_sampled.png", height = 5, width = 5)

# Figure 3D
# Example ancestral state reconstruction calibration plots
# Function to match nodes between simulated and true tree
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
ggsave(plot=true_versus_post_prob, file ="../figures/example_true_versus_post_prob_ancestor_recon.png", height = 5, width = 5)


# Figures 3 E + F 
## Strict versus state-dependent + true versus estimated plots 

#TODO: Move to separate analysis file
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

# Local records
#saveRDS(logs, file = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs/all_sampled_sizes_mcmc_logs.rds")
#saveRDS(strict_clock_logs, file = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs/strict_clock_all_sampled_sizes_mcmc_logs.rds")

# logs <- readRDS(file="/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs/all_sampled_sizes_mcmc_logs.rds")
# strict_clock_logs <- readRDS(file="/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/logs/strict_clock_all_sampled_sizes_mcmc_logs.rds")

growth_rate_posteriors <- purrr::map2(logs, log_files, 
                                      function(l, f) calculate_growth_rate_posterior_random_sampling(mcmc_log = l,
                                                                                                     log_file = f)) %>% 
    bind_rows

strict_clock_growth_rate_posteriors <- purrr::map2(strict_clock_logs, strict_clock_log_files, 
                                      function(l, f) calculate_growth_rate_posterior_random_sampling(mcmc_log = l,
                                                                                                     log_file = f)) %>% 
    bind_rows


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

# To skip above computation
growth_rate_posteriors_df <- read_csv("../eden/stats/birth_death_rate_posteriors_estimate_dr.csv")
strict_clock_growth_rate_posteriors_df <- read_csv("../eden/stats/birth_death_rate_posteriors_estimate_dr_strict_clock.csv")  


# True versus estimated plots (stat-dependent and strict comparisons,  Figure 3E)
clock_comparison_growth_rate_posteriors_df <- bind_rows(list(growth_rate_posteriors_df, strict_clock_growth_rate_posteriors_df)) %>% 
    dplyr::mutate(clock_model = factor(clock_model, levels = c('strict', 'state_dependent'))) %>% 
    dplyr::arrange(clock_model)

clock_colors <- c("state_dependent" = "#3A5A40", "strict" = "#ECA966")

clock_rate_comparison_posteriors_plot <- clock_comparison_growth_rate_posteriors_df %>% 
    dplyr::filter(n==50, minBirthRateESS > 100) %>% 
    ggplot(., aes(x=true_birth_rate_diff_weighted, y=mean_birth_rate_diff, color = clock_model)) +
    geom_point(aes(color = clock_model)) +
    theme_classic() +
    geom_abline(slope=1, intercept = 0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(ymin = birthRate_hdi95_lower, ymax = birthRate_hdi95_upper, color=clock_model), alpha = 0.6, size = 1) +
    xlab("True birth rate difference (edge - center)") +
    ylab("Estimated birth rate difference (edge - center)") +
    scale_color_manual(values = clock_colors) +
    theme(legend.position = "none") +
    theme(text=element_text(size=15))
clock_rate_comparison_posteriors_plot

ggsave(file = "../figures/clock_rate_comparison_birth_rate_posteriors.png",
       clock_rate_comparison_posteriors_plot, height = 5, width = 5)

# MSE plots (Figure 3F)
clock_comparison_sample_size_results_summary <- clock_comparison_growth_rate_posteriors_df %>% 
    dplyr::filter(minBirthRateESS > 200) %>% 
    group_by(n, clock_model) %>% 
    summarise("mse" = mean((mean_birth_rate_diff - true_birth_rate_diff_weighted)^2, na.rm = TRUE),
              "se_mse" = sd((mean_birth_rate_diff - true_birth_rate_diff_weighted)^2, na.rm = TRUE) / sqrt(n()))

clock_comparison_sample_size_mse  <- clock_comparison_sample_size_results_summary %>% 
    filter(as.numeric(n) <= 100) %>% 
    ggplot(., aes(x = n, y = mse)) +
    geom_point(aes(color=clock_model), size = 0.75) +
    theme_classic() + xlab("Sample size") + ylab("MSE of edge - center birth rate difference") +
    geom_errorbar(aes(ymin = mse - se_mse, ymax = mse + se_mse, color = clock_model)) +
    theme(text = element_text(size = 15)) +
    scale_x_continuous(breaks = pretty(clock_comparison_sample_size_results_summary$n, n = 10)) + 
    scale_color_manual(values = clock_colors) + geom_smooth(aes(color=clock_model), linetype="dashed", se = FALSE, size = 0.5)+
    theme(legend.position = "none")
clock_comparison_sample_size_mse
ggsave(file = "../figures/clock_comparison_sample_size_mse.png",
       clock_comparison_sample_size_mse, height = 5, width = 5)
