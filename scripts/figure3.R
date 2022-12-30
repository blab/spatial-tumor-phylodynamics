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

example_mcmc_log_file <- "../eden/logs/death_rate_validation_pop_1000_dr_0.29_n_100_state_clock_estimate_dr.log"

example_mcmc_log <- as.data.frame(readLog(example_mcmc_log_file))

# Get true simulated rates
## Calculated empirically by `extract_validation_sims_rates.R`
true_example_diff_df <- read_csv("../eden/stats/validation_growth_and_death_rates_weighted.csv") %>% 
    dplyr::mutate(true_birth_rate_diff_weighted = mean_edge_birth_rate - mean_center_birth_rate) %>% 
    dplyr::select(dr, true_birth_rate_diff_weighted) %>% 
    dplyr::filter(dr == "0.29")

#Get HPD intervals

example_birth_rates_diff_posteriors_summary <- example_mcmc_log %>% 
    dplyr::mutate(birthRateDiff = birthRateCanonical.2 -  birthRateCanonical.1) %>% 
    dplyr::summarise(mean = mean(birthRateDiff),
                  hpd_90_lower = hdi(birthRateDiff,
                                     credMass =0.9)[1],
                  hpd_90_upper = hdi(birthRateDiff,
                                     credMass =0.9)[2])
# Example edge - center difference posterior distibution
example_birth_rates_diff_posteriors <- example_mcmc_log %>% 
    dplyr::mutate(birthRateDiff = birthRateCanonical.2 -  birthRateCanonical.1) %>% 
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
    tidyr::pivot_longer(cols = c("birthRateCanonical.2", "birthRateCanonical.1"), 
                        names_prefix = "birthRateCanonical.",
                        values_to = "birthRate", 
                        names_to = "state") %>% 
    
    ggplot(aes(x = birthRate, fill = ifelse(state == 2, "edge", "center"))) +
    geom_density(alpha = 0.8, aes(fill = ifelse(state == 2, "edge", "center"))) +
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

# MCC tree is reconstructed from  .typed.node.trees file by BEAST2 TreeAnnotator
mcc_tree <- read.beast("../eden/trees/death_rate_validation_pop_1000_dr_0.29_n_100_state_clock_estimate_dr.typed.node.tree.mcc")


colors_loc <- colors_edge_center
names(colors_loc) <- c("loc1", "loc0")

anc_state_nodes <- mcc_tree@data %>% 
    mutate(type.prob = as.numeric(type.prob)) %>% 
    dplyr::mutate(edge = type.prob * as.integer(type == "loc1") +
                      (1 - type.prob)*as.integer(type == "loc0")) %>% 
    dplyr::mutate(center = type.prob * as.integer(type == "loc0") +
                      (1 - type.prob)*as.integer(type == "loc1")) %>% 
    dplyr::select(edge, center, type, node)


pies <- nodepie(anc_state_nodes,
                cols=1:2, 
                alpha=0.8,
                color = colors_edge_center)


g <- ggtree(mcc_tree, color = "darkgrey") +
    #geom_point(size = 2) +
    scale_color_manual(values = colors_loc) +
    theme(legend.position = "none")
    
g_pie <- ggtree::inset(g, pies, width = 0.03, height = 0.03) +
    theme(legend.position = "none")

g_pie <- g_pie + 
    geom_nodelab(aes(label = ifelse(round(as.numeric(posterior),2) < 0.99,
                                    round(as.numeric(posterior),
                                          2), "")),
                 hjust="right",
                 vjust="bottom",
                 nudge_x = -0.9,
                 nudge_y = 0.4, size = 3)

g_pie
ggsave(plot=g_pie,
       file ="../figures/example_ancestral_state_recon_tree.png",
       height = 5,
       width = 5)

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

g2_pie_toy

ggsave(plot=g2_pie_toy,
       file ="../figures/toy_ancestral_state_recon_tree.png", height = 5, width = 5)

#example_all_cells <- read_csv("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation/cells_death_rate_validation_pop_1000_dr_0.29.csv")
# example_alive_cells <- example_all_cells %>% 
#     filter_alive_cells()
# example_alive_cells <- example_alive_cells %>% 
#     mark_boundary(alive_cells = example_alive_cells)
# 
# example_alive_cells$sampled <- example_alive_cells$index %in% as.numeric(gsub("loc[01]", "", gsub("cell","",  mcc_tree@phylo$tip.label)))
#write_csv(example_alive_cells, "../simulation_data/alive_cells_death_rate_validation_pop_1000_dr_0.29.csv")

### Ancestral state reconstruction ######

# All simulated trees are reconstructed by reconstruct_sim_trees.R
sim_tree <- readRDS("../eden/simtrees/cells_death_rate_validation_pop_1000_dr_0.29.rds")

# Can skip to here
example_alive_cells <- read_csv(file="../eden/simulation_data/alive_cells_death_rate_validation_pop_1000_dr_0.29.csv")

example_sampled_cells <- example_alive_cells %>% 
    dplyr::filter(sampled)

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

#Get true differences again
true_diff_df <- read_csv("../eden/stats/validation_growth_and_death_rates_weighted.csv") %>% 
    dplyr::mutate(true_birth_rate_diff_weighted=mean_edge_birth_rate-mean_center_birth_rate)

# Posterior summary statistics for all sdevo analyses of simulated tumors are 
# generated by run_beast_to_ess.sh and compile_posterior_summaries.sh
all_posteriors_summary_df <- read_tsv("../eden/stats/posteriors/posterior_summary_all.tsv") %>% 
    left_join(., true_diff_df, by = "dr")

#strict_clock_growth_rate_posteriors_df <- read.csv("../eden/stats/birth_death_rate_posteriors_estimate_dr_strict_clock_updated.csv")

# True versus estimated plots (stat-dependent and strict comparisons,  Figure 3E)
clock_comparison_growth_rate_posteriors_df <- all_posteriors_summary_df %>% 
    dplyr::mutate(clock_model = factor(clock_model, levels = c('strict', 'state_dependent'))) %>% 
    dplyr::arrange(clock_model)

clock_colors <- c("state_dependent" = "#3A5A40", "strict" = "#ECA966")

clock_rate_comparison_posteriors_plot <- clock_comparison_growth_rate_posteriors_df %>% 
    dplyr::filter(n==50, minBirthRateESS > 200) %>% 
    ggplot(., aes(x=true_birth_rate_diff_weighted, y=mean_birth_rate_diff, color = clock_model)) +
    geom_point(aes(color = clock_model)) +
    theme_classic() +
    geom_abline(slope=1, intercept = 0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(ymin = birthRate_hdi95_lower, ymax = birthRate_hdi95_upper, color=clock_model),
                  alpha = 0.6, size = 1) +
    xlab("True birth rate difference (edge - center)") +
    ylab("Estimated birth rate difference (edge - center)") +
    scale_color_manual(values = clock_colors) +
    theme(legend.position = "none") +
    theme(text=element_text(size=15))
clock_rate_comparison_posteriors_plot

ggsave(file = "../figures/clock_rate_comparison_birth_rate_posteriors.png",
       clock_rate_comparison_posteriors_plot, height = 5, width = 5)

#strict only for presentations
clock_rate_comparison_posteriors_plot_strict_only <- clock_comparison_growth_rate_posteriors_df %>% 
    dplyr::filter(n==50, minBirthRateESS > 200) %>% 
    ggplot(., aes(x=true_birth_rate_diff_weighted, y=mean_birth_rate_diff, color = clock_model, alpha = clock_model)) +
    geom_point(aes(color = clock_model, alpha = clock_model)) +
    theme_classic() +
    geom_abline(slope=1, intercept = 0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(ymin = birthRate_hdi95_lower, ymax = birthRate_hdi95_upper, color=clock_model, alpha = clock_model), size = 1) +
    xlab("True birth rate difference (edge - center)") +
    ylab("Estimated birth rate difference (edge - center)") +
    scale_color_manual(values = clock_colors) +
    theme(legend.position = "none") +
    theme(text=element_text(size=15)) +
    scale_alpha_manual(values = c("state_dependent"=0, "strict" =0.6))
clock_rate_comparison_posteriors_plot_strict_only

ggsave(file = "../figures/clock_rate_comparison_birth_rate_posteriors_strict_only.png",
       clock_rate_comparison_posteriors_plot_strict_only, height = 5, width = 5)

# MSE plots (Figure 3F)
## Get only runs for MSE plot with same sample sizes

#function to subsampled each sample size category to the same number of simulations
subsample_results <- function(sample_size,
                              priority_drs,
                              simulation_number = 17, 
                              clock_comparison_growth_rate_posteriors_df) {
    
    posteriors_N_state_clock <- clock_comparison_growth_rate_posteriors_df %>% 
        dplyr::filter(minBirthRateESS > 200) %>% 
        dplyr::filter(n == sample_size,
                      clock_model == "state_dependent", 
                      dr %in% priority_drs)
    
    if (nrow(posteriors_N_state_clock) < simulation_number) {
        
        posteriors_N_state_clock_extra <- clock_comparison_growth_rate_posteriors_df %>% 
            dplyr::filter(minBirthRateESS > 200) %>% 
            dplyr::filter(n == sample_size,
                          clock_model == "state_dependent", 
                          ! dr %in% priority_drs)
        
        posteriors_N_state_clock_extra <- posteriors_N_state_clock_extra[sample(1:nrow(posteriors_N_state_clock_extra), 
                                              (simulation_number-nrow(posteriors_N_state_clock)), 
                                              replace = FALSE),]
        
        posteriors_N_state_clock <- bind_rows(posteriors_N_state_clock, posteriors_N_state_clock_extra)
        
        
    } else {
        
        posteriors_N_state_clock <- posteriors_N_state_clock[sample(1:nrow(posteriors_N_state_clock), 
                                                                          simulation_number, 
                                                                          replace = FALSE),]
    }
    
    posteriors_N_strict <- clock_comparison_growth_rate_posteriors_df %>% 
        dplyr::filter(minBirthRateESS > 200) %>% 
        dplyr::filter(n == sample_size,
                      clock_model == "strict", 
                      dr %in% priority_drs)
    
    if (nrow(posteriors_N_strict) < simulation_number) {
        
        posteriors_N_strict_extra <- clock_comparison_growth_rate_posteriors_df %>% 
            dplyr::filter(minBirthRateESS > 200) %>% 
            dplyr::filter(n == sample_size,
                          clock_model == "strict", 
                          ! dr %in% priority_drs)
        
        posteriors_N_strict_extra <- posteriors_N_strict_extra[sample(1:nrow(posteriors_N_strict_extra), 
                                                                                (simulation_number-nrow(posteriors_N_strict)), 
                                                                                replace = FALSE),]
        
        posteriors_N_strict <- bind_rows(posteriors_N_strict, posteriors_N_strict_extra)
        
        
    } else {
        
        posteriors_N_strict <- posteriors_N_strict[sample(1:nrow(posteriors_N_strict), 
                                                                    simulation_number, 
                                                                    replace = FALSE),]
    }
    
    posteriors_N <- bind_rows(posteriors_N_state_clock, posteriors_N_strict)
    
    return(posteriors_N)
    
    
}
set.seed(9811)
priority_drs = unique(clock_comparison_growth_rate_posteriors_df$dr[clock_comparison_growth_rate_posteriors_df$n == 5])
standardized_sim_number_df <- purrr::map(unique(clock_comparison_growth_rate_posteriors_df$n),
                                         function(ss) subsample_results(sample_size = ss,
                                                                        simulation_number = 17,
                                                                        priority_drs = priority_drs,
                                                                        clock_comparison_growth_rate_posteriors_df = clock_comparison_growth_rate_posteriors_df)) %>% 
    
    bind_rows()
clock_comparison_sample_size_results_summary <- clock_comparison_growth_rate_posteriors_df  %>% 
    
    dplyr::filter(minBirthRateESS > 200) %>% 
    group_by(n, clock_model) %>% 
    summarise("mse" = mean((mean_birth_rate_diff - true_birth_rate_diff_weighted)^2, na.rm = TRUE),
              "se_mse" = sd((mean_birth_rate_diff - true_birth_rate_diff_weighted)^2, na.rm = TRUE) / sqrt(n()),
              "N" = n())

#For caption stats
sum(clock_comparison_sample_size_results_summary$N)
max(clock_comparison_sample_size_results_summary$N)
min(clock_comparison_sample_size_results_summary$N)
mean(clock_comparison_sample_size_results_summary$N)

ggplot(clock_comparison_sample_size_results_summary, aes(x=n, y=N)) + geom_point()
clock_comparison_sample_size_mse  <- clock_comparison_sample_size_results_summary %>% 
    dplyr::filter(as.numeric(n) <= 100) %>% 
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
