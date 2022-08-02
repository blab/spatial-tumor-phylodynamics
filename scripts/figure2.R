#figure2.R

##script to generate visualizations for figure 2
### Note that Figure 2A insets are generated from Figure 1 script 


##### SETUP #####
library(tumortree)
library(tidyverse)
library(viridis)
library(ggtree)
library(phangorn)
library(ape)

figures_dir <- "../figures"

#sim_colors <- tumortree::get_color_palette(c("boundary_driven", "unrestricted"))
sim_colors <- c("boundary_driven" = "#e49a8b", "unrestricted" = "#3C3C3C")
colors_edge_center <- tumortree::get_color_palette(c("edge", "center"))



# Function to get lineage time spent on edge
get_edge_time <- function(leaf, tree_full) {
    mean_edge_time <- weighted.mean(tree_full@data$state[c(leaf, phangorn::Ancestors(tree_full@phylo, node = leaf, type ="all"))],
                                    tree_full@data$branch_length[c(leaf, phangorn::Ancestors(tree_full@phylo, node = leaf, type ="all"))], na.rm = TRUE)
    return(mean_edge_time)
}

# Read in example simulated tumors to compare
## BOUNDARY_DRIVEN
all_cells_boundary_driven <- read_csv("../eden/simulation_data/cells_death_rate_validation_pop_1000_dr_0.10.csv") %>%
    normalize_locs


## filter to alive cells
alive_cells_boundary_driven <- all_cells_boundary_driven %>%
    filter_alive_cells

alive_cells_boundary_driven <- mark_boundary(alive_cells_boundary_driven, alive_cells_boundary_driven)

## UNRESTRICTED 
all_cells_unrestricted <- read_csv("../eden/simulation_data/cells_pushing_pop_1000_dr_0.10.csv") %>%
    normalize_locs

alive_cells_unrestricted <- all_cells_unrestricted %>%
    filter_alive_cells

alive_cells_unrestricted <- mark_boundary(alive_cells_unrestricted, alive_cells_unrestricted)

# Read in subsampled boundary-driven trees
tree_boundary_driven <- readRDS(file = "../eden/simtrees/boundary_driven_timetree.rds")
tree_molecular_boundary_driven <- readRDS(file = "../eden/simtrees/boundary_driven_moleculartree.rds")


# Get sampled cells from saved trees
sampled_cells_boundary_driven <- alive_cells_boundary_driven
sampled_cells_boundary_driven$sampled <- alive_cells_boundary_driven$index %in%
    tree_boundary_driven@data$index

sampled_cells_boundary_driven <- sampled_cells_boundary_driven %>%
    dplyr::filter(sampled == TRUE)

### Saved trees for unrestricted growth

tree_unrestricted <- readRDS(file = "../eden/simtrees/unrestricted_timetree.rds")
tree_molecular_unrestricted <- readRDS(file = "../eden/simtrees/unrestricted_moleculartree.rds")


## Also created sampled cells data frame from saved tree for unrestricted growth simulation
sampled_cells_unrestricted <- alive_cells_unrestricted
sampled_cells_unrestricted$sampled <- alive_cells_unrestricted$index %in% tree_unrestricted@data$index
sampled_cells_unrestricted <- sampled_cells_unrestricted %>%
    filter(sampled == TRUE)

## mark boundary to label edge and center
# sampled_cells_boundary_driven <- mark_boundary(sampled_cells_boundary_driven, alive_cells_boundary_driven)
# nearest_points <- spatstat.geom::nndist(X= sampled_cells_boundary_driven$locx, Y=sampled_cells_boundary_driven$locy, by = as.logical(sampled_cells_boundary_driven$est_edge))
# sampled_cells_boundary_driven$dist_from_edge <- (1 - sampled_cells_boundary_driven$est_edge) * nearest_points[,which(colnames(nearest_points) == "TRUE")]

# Make full tree with state labels added to each node

## Time tree
tree_boundary_driven_full <- convert_all_cells_to_tree_fast(all_cells = all_cells_boundary_driven,
                                            add_all_states = TRUE,
                                            branch_unit = "time",
                                            cell_locations_df_file = NULL)


tree_boundary_driven_full <- prune_simulated_tree(tree_boundary_driven_full, sampled_cells_indices = alive_cells_boundary_driven$index)
## Genetic (molecular) tree
tree_molecular_boundary_driven_full <- convert_all_cells_to_tree_fast(all_cells = all_cells_boundary_driven,
                                                                   add_all_states = FALSE,
                                                                   branch_unit = "genetic",
                                                                   cell_locations_df_file = NULL)

tree_molecular_boundary_driven_full <- prune_simulated_tree(tree_molecular_boundary_driven_full,
                                                            sampled_cells_indices = alive_cells_boundary_driven$index, 
                                                            branch_unit = "genetic")

tree_molecular_boundary_driven_full@phylo$edge.length[is.na(tree_molecular_boundary_driven_full@phylo$edge.length)] <- 0

## Add lineage time on edge for each tip
tree_boundary_driven_full@data$frac_time_on_edge <- NA
tree_boundary_driven_full@data$frac_time_on_edge[1:length(tree_boundary_driven_full@phylo$tip.label)] <- purrr::map_dbl(1:length(tree_boundary_driven_full@phylo$tip.label), function(l) get_edge_time(leaf = l, tree_full = tree_boundary_driven_full))


## Add clock rate info
tree_boundary_driven_full@data <- tree_boundary_driven_full@data %>% 
    dplyr::mutate(n_muts = str_count(mutations, ",") + 1) %>% 
    dplyr::mutate(clock_rate = n_muts/max(tree_boundary_driven_full@data$deathdate))

## Visualizations of lineage time spend on edge in extant tumor
tree_boundary_driven_full@data %>% 
    filter(deathdate == max(tree_boundary_driven_full@data$deathdate, na.rm = TRUE)) %>% 
    ggplot(aes(x=locx, y=locy, color=frac_time_on_edge)) + geom_point(size =4) + theme_void() + scale_color_viridis()

#### Unrestricted growth -- make full tree w/ state labels

tree_unrestricted_full <- convert_all_cells_to_tree_fast(all_cells = all_cells_unrestricted,
                                                            add_all_states = TRUE,
                                                            cell_locations_df_file =  "../eden/simulation_data/cells_pushing_pop_1000_dr_0.10_locs.csv",
                                                            branch_unit = "time")

tree_unrestricted_full <- prune_simulated_tree(tree_unrestricted_full, sampled_cells_indices = alive_cells_unrestricted$index)


tree_molecular_unrestricted_full <- convert_all_cells_to_tree_fast(all_cells = all_cells_unrestricted,
                                                         add_all_states = FALSE,
                                                         cell_locations_df_file = NULL,
                                                         branch_unit = "genetic")

tree_molecular_unrestricted_full <- prune_simulated_tree(tree_molecular_unrestricted_full,
                                                         sampled_cells_indices = alive_cells_unrestricted$index, 
                                                         branch_unit = "genetic")

tree_molecular_unrestricted_full@phylo$edge.length[is.na(tree_molecular_unrestricted_full@phylo$edge.length)] <- 0

# Add fraction time on edge per tip
tree_unrestricted_full@data$frac_time_on_edge <- NA
tree_unrestricted_full@data$frac_time_on_edge[1:length(tree_unrestricted_full@phylo$tip.label)] <- purrr::map_dbl(1:length(tree_unrestricted_full@phylo$tip.label), function(l) get_edge_time(leaf = l, tree_full = tree_unrestricted_full))


# Add clock rate info
tree_unrestricted_full@data <- tree_unrestricted_full@data %>% 
    dplyr::mutate(n_muts = str_count(mutations, ",") + 1) %>% 
    dplyr::mutate(clock_rate = n_muts/max(tree_unrestricted_full@data$deathdate))

# Transfer info to genetic tree
tree_molecular_unrestricted_full@data <- tree_unrestricted_full@data
tree_molecular_boundary_driven_full@data <- tree_boundary_driven_full@data


# Prune tree to sampled cells

## Boundary-driven
tree_boundary_driven <- prune_simulated_tree(tree_boundary_driven_full,
                                             sampled_cells_indices = sampled_cells_boundary_driven$index)
tree_unrestricted <- prune_simulated_tree(tree_unrestricted_full,
                                          sampled_cells_indices = sampled_cells_unrestricted$index)

## Unrestricted
tree_molecular_unrestricted <- prune_simulated_tree(tree_molecular_unrestricted_full,
                                                    sampled_cells_indices = sampled_cells_unrestricted$index)

tree_molecular_boundary_driven <- prune_simulated_tree(tree_molecular_boundary_driven_full,
                                                       sampled_cells_indices = sampled_cells_boundary_driven$index)

# Save trees
# saveRDS(tree_unrestricted_full, file = "../eden/simtrees/unrestricted_timetree_full.rds")
# saveRDS(tree_boundary_driven_full, file = "../eden/simtrees/boundary_driven_timetree_full.rds")
# saveRDS(tree_molecular_boundary_driven_full, file = "../eden/simtrees/boundary_driven_molecular_full.rds")
# saveRDS(tree_molecular_unrestricted_full, file = "../eden/simtrees/unrestricted_molecular_full.rds")

# To skip above steps load trees directly 

#tree_boundary_driven_full <- readRDS(file = "../eden/simtrees/boundary_driven_timetree_full.rds")
#tree_unrestricted_full <- readRDS(file = "../eden/simtrees/unrestricted_timetree_full.rds")
#tree_molecular_boundary_driven_full <- readRDS(file = "../eden/simtrees/boundary_driven_molecular_full.rds")
#tree_molecular_unrestricted_full <- readRDS(file = "../eden/simtrees/unrestricted_molecular_full.rds")

##### FIGURES 2A and 2D #######

# Function to get distance from edge versus number of divisions across simulations at same death rate
get_sim_edge_dist_versus_divisions <- function(all_cells_file, cell_locations_df_file = NULL, time_cutoff=0) {
    
    print(all_cells_file)
    
    # Get all simulated cells
    all_cells <- read_csv(all_cells_file)
    
    #extract simulation information from fillle name
    dr_extract <- regmatches(basename(all_cells_file),
                             gregexpr("(?<=dr_)[[:digit:]]+.[[:digit:]]+", basename(all_cells_file), perl = TRUE))[[1]]
    i_extract <- regmatches(basename(all_cells_file),
                             gregexpr("(?<=i_)[[:digit:]]+", basename(all_cells_file), perl = TRUE))[[1]]
    model_extract <- ifelse(grepl("cells_death_rate_validation", all_cells_file), "boundary_driven",
                            ifelse(grep("cells_pushing", all_cells_file), "unrestricted", NA))
    #get time at end of simulation
    endpoint <- max(all_cells$deathdate)
    
    #get vector of timespoints to check 
    timepoints <- seq(2/24, endpoint - 6/24, by = 2/24)
    
    #
    timepoints <- timepoints[timepoints > (time_cutoff*max(timepoints))]

    edge_dist_versus_divisions_df <- purrr::map(timepoints , function(t) get_edge_dist_versus_divisions_at_t(all_cells,
                                                                                                             time = t,
                                                                                                             cell_locations_df_file = cell_locations_df_file)) %>%
        bind_rows() %>% 
        add_column("i" = i_extract, "dr" = dr_extract, "model" = model_extract)
    
    return(edge_dist_versus_divisions_df)
    
}


all_cells_files_dr_0.10 <- list.files(path="../eden/simulation_data",
                                      pattern="cells_death_rate_validation_pop_1000_dr_0.10_i_[0-9]+.csv",
                                      full.names = TRUE)

all_cells_files_unrestricted_dr_0.10 <- list.files(path="../eden/simulation_data",
                                      pattern="cells_pushing_pop2_1000_dr_0.10_i_[0-9]+.csv",
                                      full.names = TRUE)

all_cells_files_unrestricted_dr_0.10_locs <- list.files(path="../eden/simulation_data",
                                                   pattern="cells_pushing_pop2_1000_dr_0.10_i_[0-9]+_locs.csv",
                                                   full.names = TRUE)

#get distance to edge versus birth rate calculations for dr 0.10 + boundary-driven growth
all_sims_edge_dist_versus_divisions_dr_0.10_df <- purrr::map(all_cells_files_dr_0.10, function(cells) get_sim_edge_dist_versus_divisions(all_cells_file = cells,
                                                                                                                                         time_cutoff = 0)) %>%
    bind_rows()

# get distance to edge versus birth rate calculations for dr 0.10 + unrestricted
all_sims_edge_dist_versus_divisions_unrestricted_dr_0.10_df <- purrr::map2(all_cells_files_unrestricted_dr_0.10,all_cells_files_unrestricted_dr_0.10_locs,
                                                                           function(cells,locs) get_sim_edge_dist_versus_divisions(all_cells_file = cells,
                                                                                                                                   cell_locations_df_file = locs,
                                                                                                                                   time_cutoff = 0)) %>% 
    bind_rows()


# Save calculations to CSV file

## These files are too large for github
# write.csv(all_sims_edge_dist_versus_divisions_dr_0.10_df,
#           file="../eden/stats/edge_dist_versus_growth_rate_sim_data_dr_0.10.csv")
# 
# 
# write.csv(all_sims_edge_dist_versus_divisions_unrestricted_dr_0.10_df,
#           file="../eden/stats/edge_dist_versus_growth_rate_sim_data_unrestricted_dr_0.10.csv")


# To skip calculations and get output directly
# all_sims_edge_dist_versus_divisions_dr_0.10_df <- read_csv(file = "..eden/stats/edge_dist_versus_growth_rate_sim_data_dr_0.10.csv")
# ll_sims_edge_dist_versus_divisions_unrestricted_dr_0.10_df <- read_csv(file = "../eden/stats/edge_dist_versus_growth_rate_sim_data_unrestricted_dr_0.10.csv")

#Calculate mean and se of birth rate across binned disistances from the tumor edge
time_step <- 2/24
all_dr_0.10_sims_boundary_driven_summary <- all_sims_edge_dist_versus_divisions_dr_0.10_df %>% 
    
    #first normalize distances accross simulations
    dplyr::group_by(i,time) %>% 
    dplyr::mutate(norm_dist_from_edge = dist_from_edge/max(dist_from_edge)) %>% 

    ungroup %>% 
    
    #then bin distances to get binned birth rate
    dplyr::mutate(norm_dist_from_edge_bin = cut(norm_dist_from_edge, breaks=13)) %>% 
    dplyr::group_by(norm_dist_from_edge_bin, i) %>% 
    dplyr::summarise("birth_rate" = sum(divided)/n()/time_step) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(norm_dist_from_edge_bin) %>% 
    
    #get mean and se accross simulations
    dplyr::summarise("birth_rate_mean" = mean(birth_rate), "birth_rate_se" = sqrt(var(birth_rate) / n())) 



#to make histogram extract bin ranges and get midpoint for plotting
all_dr_0.10_sims_boundary_driven_summary$binRange <- str_extract(all_dr_0.10_sims_boundary_driven_summary$norm_dist_from_edge_bin, "[0-9].*[0-9]+")
all_dr_0.10_sims_boundary_driven_summary <- all_dr_0.10_sims_boundary_driven_summary %>% 
    tidyr::separate(binRange, into=c("binStart", "binEnd"), sep = ",") %>% 
    dplyr::mutate(binStart = as.numeric(binStart), binEnd = as.numeric(binEnd)) %>% 
    dplyr::mutate(binMid = (binEnd + binStart)/2)

#get bin width to fit bars together
binWidth <- mean(all_dr_0.10_sims_boundary_driven_summary$binEnd - all_dr_0.10_sims_boundary_driven_summary$binStart, na.rm = TRUE)


#Write summary stats for quick plotting
write.csv(all_dr_0.10_sims_boundary_driven_summary,
          file="../eden/stats/edge_dist_versus_growth_rate_sim_data_dr_0.10_summary.csv")

all_dr_0.10_sims_boundary_driven_summary <- read.csv("../eden/stats/edge_dist_versus_growth_rate_sim_data_dr_0.10_summary.csv")
#plotting
all_dr_0.10_sims_boundary_driven_hist <- ggplot(all_dr_0.10_sims_boundary_driven_summary,
                                                alpha = 0.5,
                                                aes(x = binMid, y = birth_rate_mean
                               ),
           color = "black") +
    
    geom_bar(stat="identity", fill = sim_colors["boundary_driven"], color = "black", width=binWidth*0.99, alpha = 0.9) +
    geom_point() +
    geom_errorbar(aes(ymin=birth_rate_mean - birth_rate_se, ymax=birth_rate_mean + birth_rate_se), width = binWidth*0.2) +
    xlab("Fraction distance from tumor edge") +
    ylab("Birth rate") +
    theme_classic() +
    theme(text = element_text(size = 15))

all_dr_0.10_sims_boundary_driven_hist   

##histogram for unrestricted growth

all_dr_0.10_sims_unrestricted_summary <- all_sims_edge_dist_versus_divisions_unrestricted_dr_0.10_df %>% 
    
    #first normalize distances across simulations
    dplyr::group_by(i,time) %>% 
    dplyr::mutate(norm_dist_from_edge = dist_from_edge/max(dist_from_edge)) %>% 
    ungroup %>% 
    
    #then bin distances to get binned birth rate
    dplyr::mutate(norm_dist_from_edge_bin = cut(norm_dist_from_edge, breaks=13)) %>% 
    dplyr::group_by(norm_dist_from_edge_bin, i) %>% 
    dplyr::summarise("birth_rate" = sum(divided)/n()/time_step) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(norm_dist_from_edge_bin) %>% 
    
    #get mean and se accross simulations
    dplyr::summarise("birth_rate_mean" = mean(birth_rate), "birth_rate_se" = sqrt(var(birth_rate) / n())) 

#to make histogram extract bin ranges and get midpoint for plotting
all_dr_0.10_sims_unrestricted_summary$binRange <- str_extract(all_dr_0.10_sims_unrestricted_summary$norm_dist_from_edge_bin, "[0-9].*[0-9]+")
all_dr_0.10_sims_unrestricted_summary <- all_dr_0.10_sims_unrestricted_summary %>% 
    tidyr::separate(binRange, into=c("binStart", "binEnd"), sep = ",") %>% 
    dplyr::mutate(binStart = as.numeric(binStart), binEnd = as.numeric(binEnd)) %>% 
    dplyr::mutate(binMid = (binEnd + binStart)/2)

#get bin width to fit bars together
binWidth <- mean(all_dr_0.10_sims_unrestricted_summary$binEnd - all_dr_0.10_sims_unrestricted_summary$binStart, na.rm = TRUE)

write.csv(all_dr_0.10_sims_unrestricted_summary,
          file="../eden/stats/edge_dist_versus_growth_rate_sim_data_unrestricted_dr_0.10_summary.csv")

all_dr_0.10_sims_unrestricted_summary <- read.csv("../eden/stats/edge_dist_versus_growth_rate_sim_data_unrestricted_dr_0.10_summary.csv")
#plotting
all_dr_0.10_sims_unrestricted_hist <- ggplot(all_dr_0.10_sims_unrestricted_summary, alpha = 0.5,
                                             aes(x = binMid, y = birth_rate_mean),color = "black") +
    
    geom_bar(stat="identity", fill = sim_colors["unrestricted"], color = "black", width=binWidth*0.99, alpha = 0.9) +
    geom_errorbar(aes(ymin=birth_rate_mean - birth_rate_se,
                      ymax=birth_rate_mean + birth_rate_se), width = binWidth*0.2) +
    xlab("Fraction distance from tumor edge") +
    ylab("Birth rate") +
    theme_classic() +
    theme(text = element_text(size = 15))

all_dr_0.10_sims_unrestricted_hist

#Normalize Y-axes
max_y <- max(c(all_dr_0.10_sims_unrestricted_summary$birth_rate_mean + all_dr_0.10_sims_unrestricted_summary$birth_rate_se,
                all_dr_0.10_sims_boundary_driven_summary$birth_rate_mean + all_dr_0.10_sims_boundary_driven_summary$birth_rate_se))

all_dr_0.10_sims_unrestricted_hist <- all_dr_0.10_sims_unrestricted_hist + coord_cartesian(ylim = c(0, max_y))
all_dr_0.10_sims_boundary_driven_hist <- all_dr_0.10_sims_boundary_driven_hist + coord_cartesian(ylim = c(0, max_y))

ggsave(file = "edge_dist_versus_birth_rate_hist_boundary_driven_all_dr_0.10.png", plot = all_dr_0.10_sims_boundary_driven_hist, path = figures_dir, height = 5, width = 5)
ggsave(file = "edge_dist_versus_birth_rate_hist_unrestricted_all_dr_0.10.png", plot = all_dr_0.10_sims_unrestricted_hist, path = figures_dir, height = 5, width = 5)

### END SUBFIGURE ###

### Figure 2C BOUNDARY-DRIVEN INSET ###
# Inset genetic tree clade colored by time spend on edge and center
c <- ggtree(tree_molecular_boundary_driven, aes(color = frac_time_on_edge)) +
    geom_tippoint(aes(color = frac_time_on_edge), size = 3, alpha = 1) +
    scale_color_gradient(low = colors_edge_center["center"],  high = colors_edge_center["edge"])  +
    coord_cartesian(clip = "off") + theme(legend.position = "none")

c + geom_tiplab()
c
c <- viewClade(c, MRCA(c,"cell_2620", "cell_1408"))
c
### Figure 2F INSET UNRESTRICTED ###
tree_molecular_unrestricted@data <- tree_unrestricted@data
d <- ggtree(tree_molecular_unrestricted, aes(color = frac_time_on_edge)) +
    geom_tippoint(aes(color = frac_time_on_edge), size = 3, alpha = 1) +
    scale_color_gradient(low = colors_edge_center["center"],  high = colors_edge_center["edge"])  +
    coord_cartesian(clip = "off") +
    theme(legend.position = "none") #+ geom_tiplab(color = "black")
d + geom_tiplab()

d <- viewClade(d, MRCA(d,"cell_1919", "cell_880"))
d

## BOUNDARY-DRIVEN Growth Figure 1C

# Normalize axe
ymax <- max(c(tree_boundary_driven_full@data$clock_rate, tree_unrestricted@data$clock_rate), na.rm = TRUE)
ymin <- min(c(tree_boundary_driven_full@data$clock_rate, tree_unrestricted@data$clock_rate), na.rm = TRUE)
#tree_boundary_driven@data$clock_rate <- sampled_cells_boundary_driven$clock_rate[match(tree_boundary_driven@data$index, sampled_cells_boundary_driven$index)]
alive_boundary_driven_data <- tree_boundary_driven_full@data %>% 
    filter(deathdate == max(tree_boundary_driven@data$deathdate))

model <- lm(clock_rate~frac_time_on_edge, data=alive_boundary_driven_data)


boundary_driven_r_squared <- summary(model)$r.squared
e <- ggplot(alive_boundary_driven_data, aes(x = frac_time_on_edge,  y = clock_rate, color = frac_time_on_edge)) +
    geom_point(size = 2,
               alpha = 1) +
    theme_classic() +
    geom_smooth(method = "lm", linetype = "dashed", color = "black", size = 0.5, fill = "grey") +
    theme(legend.position = "none") +
    scale_color_gradient(low = colors_edge_center["center"],  high = colors_edge_center["edge"]) +
    xlab("Fraction lineage time on edge") +
    ylab("Mean clock rate") +
    # theme(axis.title.x=element_blank(),
    #       axis.text.x=element_blank()) +
    annotate("text", x = 0.8, y = 0.35, size = 4,
             label = paste0("R ** {2} == ", format(boundary_driven_r_squared , digits = 2)), 
             color = "black", parse = TRUE) +
    theme(text=element_text(size = 15)) +
    ylim(c(ymin,ymax))

e

## UNRESTRICTED_GROWTH Figure 1F


alive_unrestricted_data <- tree_unrestricted_full@data %>% 
    filter(deathdate == max(tree_unrestricted@data$deathdate))

model <- lm(clock_rate~frac_time_on_edge, data=alive_unrestricted_data )

unrestricted_r_squared <- summary(model)$r.squared
f <- ggplot(alive_unrestricted_data, aes(x = frac_time_on_edge,  y = clock_rate, color = frac_time_on_edge)) +
    geom_point(size = 2,
               alpha = 1) +
    theme_classic() +
    geom_smooth(method = "lm", linetype = "dashed", color = "black", size = 0.5, fill = "grey") +
    theme(legend.position = "none") +
    scale_color_gradient(low = colors_edge_center["center"],  high = colors_edge_center["edge"]) +
    xlab("Fraction lineage time on edge") +
    ylab("Mean clock rate") +
    annotate("text", x = 0.55, y = 0.35, size = 4,
             label = paste0("R ** {2} == ", format(unrestricted_r_squared , digits = 2)), 
             color = "black", parse = TRUE) +
    theme(text=element_text(size = 15)) +
    ylim(c(ymin,ymax))



f
#
# ###COMBINED FIGURE ###
#
# e <- e + coord_cartesian(ylim = c(0.005, 0.045))
# f <- f + coord_cartesian(ylim = c(0.02, 0.1))

#c <- c + border(linetype = "dashed", color = "darkgrey", size = 0.8)
#d <- d + border(linetype = "dashed", color = "darkgrey", size = 0.8)
#
g <- ggdraw() +
    draw_plot(e)  +
    draw_plot(c, x = 0.2, y = 0.7, width = 0.25, height = 0.3 )
g

h <- ggdraw() +
    draw_plot(f)  +
    draw_plot(d, x = 0.2, y = .73, width = .25, height = .3)
h
g
ggsave(file = "mean_clock_rate_vs_frac_lineage_time_on_edge_with_tree_boundary_driven.png", plot = g, path = figures_dir, height = 5, width = 5)
ggsave(file = "mean_clock_rate_vs_frac_lineage_time_on_edge_with_tree_unrestricted.png", plot = h, path = figures_dir, height = 5, width = 5)


#Terminal branch length figure


#compare variance in terminal branch lengths
tree_height_boundary_driven <- max(tree_boundary_driven@data$deathdate)
tree_height_unrestricted <- max(tree_unrestricted@data$deathdate)
terminal_branch_length_df <- data.frame("norm_terminal_branch_length" = c(tree_boundary_driven_full@phylo$edge.length[which(tree_boundary_driven_full@phylo$edge[,2] <= length(tree_boundary_driven_full@phylo$tip.label))]/tree_height_boundary_driven,
                                                                     tree_unrestricted_full@phylo$edge.length[which(tree_unrestricted_full@phylo$edge[,2] <= length(tree_unrestricted_full@phylo$tip.label))]/tree_height_unrestricted),
                                        "state" = c(tree_boundary_driven_full@data$state[1:length(tree_boundary_driven_full@phylo$tip.label)],
                                                    tree_unrestricted_full@data$state[1:length(tree_unrestricted_full@phylo$tip.label)]),
                                        "model" = c(rep("boundary-driven", length(tree_boundary_driven_full@phylo$tip.label)), rep("unrestricted",length(tree_unrestricted_full@phylo$tip.label))),
                                        "frac_time_on_edge" = c(tree_boundary_driven_full@data$frac_time_on_edge[1:length(tree_boundary_driven_full@phylo$tip.label)],
                                                                tree_unrestricted_full@data$frac_time_on_edge[1:length(tree_unrestricted_full@phylo$tip.label)]))

## BOUNDARY-DRIVEN GROWTH ##
mean_ratio_branch_length <- mean(terminal_branch_length_df$norm_terminal_branch_length[terminal_branch_length_df$state == 0 & terminal_branch_length_df$model == "boundary-driven"]) / mean(terminal_branch_length_df$norm_terminal_branch_length[terminal_branch_length_df$state == 1  & terminal_branch_length_df$model == "boundary-driven"])
i <- terminal_branch_length_df %>%
    filter(model == "boundary-driven") %>%
    ggplot(., aes(x = ifelse(state == 1, "edge", "center"), y = norm_terminal_branch_length, fill = ifelse(state == 1, "edge", "center"))) +
    geom_violin(position="dodge", alpha = 1) + theme_classic() +
    #ggplot(., aes(x = ifelse(state == 1, "edge", "center"), y = terminal_branch_length)) +
    #geom_violin(position="dodge", alpha = 0.8, color = "black") + theme_classic() +
    # geom_jitter(size =2, aes(x = ifelse(state == 1, "edge", "center"),
    #                 y = terminal_branch_length, color = frac_time_on_edge),
    #             width = 0.2, height = 0) +
    scale_fill_manual(values = colors_edge_center) +
    stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
                 colour = "black",
                 size = 0.25) +
    xlab("") + ylab("Terminal branch length") +
    theme(text=element_text(size = 15)) +
    theme(legend.position = "none") +
    # theme(axis.title.x=element_blank(),
    # #       axis.text.x=element_blank())   +
    annotate("text", x = 1, y = 0.9, size = 4,
             label = paste0("frac('Mean center TBL', 'Mean edge TBL') == ", format(mean_ratio_branch_length, digits = 2)),
                                                    color = "black", parse = TRUE) +
    ylim(c(0,1)) 

i


## UNRESTRICTED ##
mean_ratio_branch_length_unrestricted <- mean(terminal_branch_length_df$norm_terminal_branch_length[terminal_branch_length_df$state == 0 & terminal_branch_length_df$model == "unrestricted"]) / mean(terminal_branch_length_df$norm_terminal_branch_length[terminal_branch_length_df$state == 1  & terminal_branch_length_df$model == "unrestricted"])

j <- terminal_branch_length_df %>%
    filter(model == "unrestricted") %>%
    ggplot(., aes(x = ifelse(state == 1, "edge", "center"), y = norm_terminal_branch_length, fill = ifelse(state == 1, "edge", "center"))) +
    geom_violin(position="dodge", alpha = 1) + theme_classic() +
    scale_fill_manual(values = colors_edge_center) +
    stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
                 colour = "black",
                 size = 0.25) +
    xlab("") + ylab("Terminal branch length") +
    theme(text=element_text(size = 15)) +
    theme(legend.position = "none")  +
    annotate("text", x = 1, y = 0.9, size = 4,
             label = paste0("frac('Mean center TBL', 'Mean edge TBL') == ", format(mean_ratio_branch_length_unrestricted, digits = 2)), 
             color = "black", parse = TRUE) +
    ylim(c(0,1))

j

### BOUNDARY-DRIVEN ###
i_sub <- ggtree(tree_boundary_driven, aes(color = ifelse(state == 1 & node <= 100, "edge", ifelse(node <= 100, "center", NA)))) +
    geom_tippoint( aes(color = ifelse(state == 1, "edge","center")), size = 3, alpha = 1) +
    scale_color_manual(values = colors_edge_center)  +
    coord_cartesian(clip = "off") +
    theme(legend.position = "none")

i_sub <- viewClade(i_sub, MRCA(c,"cell_2620", "cell_1408"))
i_sub
### UNRESTRICTED ###

j_sub <- ggtree(tree_unrestricted, aes(color = ifelse(state == 1 & node <= 100, "edge", ifelse(node <= 100, "center", NA)))) +
    geom_tippoint( aes(color = ifelse(state == 1, "edge","center")), size = 3, alpha = 1) +
    scale_color_manual(values = colors_edge_center)  +
    coord_cartesian(clip = "off") +
    theme(legend.position = "none")

j_sub <- viewClade(j_sub, MRCA(d,"cell_1919", "cell_880"))

j_sub

# ###COMBINED FIGURE.png ###
library(cowplot)
#
i.with.inset <- ggdraw() +
    draw_plot(i)  +
    draw_plot(i_sub, x = 0.7, y = .7, width = .25, height = .3)
i.with.inset

j.with.inset <- ggdraw() +
    draw_plot(j)  +
    draw_plot(j_sub, x = 0.7, y = .7, width = .25, height = .3)
j.with.inset

ggsave(file = "terminal_branch_length_violin_with_tree_boundary_driven.png", plot = i.with.inset, path = figures_dir, height = 5, width = 5)
ggsave(file = "terminal_branch_length_violin_with_tree_unrestricted.png", plot = j.with.inset, path = figures_dir, height = 5, width = 5)


# Function to get tree stats broken into edge and center states from trees
get_corr_edge_time_clock_rate_with_tree <- function(sim_tree_file) {
    print(sim_tree_file)
    #tree_file_name <- paste0(gsub("csv", "rds", basename(sim_cells_tree_file)))
    dr_extract <- regmatches(basename(sim_tree_file),
                             gregexpr("(?<=dr_)[[:digit:]]+.[[:digit:]]+", basename(sim_tree_file), perl = TRUE))[[1]]
    
    tree_full <- readRDS(sim_tree_file)
    alive_tree_data <- tree_full@data %>% 
        filter(deathdate == max(tree_full@data$deathdate, na.rm = TRUE))
    
    model <- lm(clock_rate~frac_time_on_edge, data=alive_tree_data)
    
    if (grepl("cells_death_rate_validation", sim_tree_file)) {
        model_extract = "boundary_driven"
    } else {
        model_extract = "unrestricted"
    }
    
    alive_tree <- prune_simulated_tree(tree_full, sampled_cells_indices = alive_tree_data$index)
    
    edge_terminal_branch_length <- alive_tree@data$cell_name[which(alive_tree@data$state == 1)]
    leaf_states <- alive_tree@data$state[match(alive_tree@phylo$tip.label, alive_tree@data$cell_name)]
    leaf_terminal_branches <- alive_tree@phylo$edge.length[match(1:length(alive_tree@phylo$tip.label), alive_tree@phylo$edge[,2])]
    

    #don't need to normalize because ratio stays the same
    edge_mean_terminal_branch <- mean(leaf_terminal_branches[leaf_states == 1])
    center_mean_terminal_branch <- mean(leaf_terminal_branches[leaf_states == 0])
    return(data.frame("r_squared" = summary(model)$r.squared, 
                      "dr" = dr_extract,
                      "model" = model_extract,
                      "terminal_branch_ratio" = center_mean_terminal_branch / edge_mean_terminal_branch))
}


#### Get stats for range of simulations (boundary-driven and unrestricted)

# Get saved simulated tree files to load 
## These are generated by reconstruct_sim_trees.R
sim_trees_files <- list.files(path="../eden/simtrees",
                              pattern = "[0-9][0-9].rds", include.dirs = TRUE, full.names = TRUE)

sim_trees_files <- sim_trees_files[(! grepl("_i_", sim_trees_files))]


summary_stats <- purrr::map(sim_trees_files, function(file) get_corr_edge_time_clock_rate_with_tree(file)) %>% 
    bind_rows()

# Save to CSV
write.csv(summary_stats, file="../eden/stats/model_summary_stats.csv")

# To get these stats without having trees generated
# summary_stats <- read_csv(file="../eden/stats/model_summary_stats.csv")
# Extract values to highlight from example tumors
highlighted_sims <- summary_stats  %>% 
    dplyr::filter(dr == "0.10") 


k <- ggplot(summary_stats, aes(x = as.numeric(dr)/2, y = terminal_branch_ratio, color = model)) +
    geom_point( size = 2) + theme_classic() + scale_color_manual(values = sim_colors,
                                                                             labels = c("Boundary-driven", "Unrestricted")) +
    labs(x = "Cell turnover - P(cell death)", y = "Mean center / edge terminal branch lengths") +
    theme(legend.position = c(0.8, 0.85),
          #legend.background = element_rect(linetype = "dashed", size = 0.5, colour = "darkgrey"),
          legend.title = element_blank(),
          text=element_text(size = 15)) +
    
    guides(color = guide_legend(override.aes = list(size = 3) ) ) +
    geom_point(data = highlighted_sims, color = "darkgrey", size = 5, alpha = 0.5) 
k

ggsave(file = "center_over_edge_terminal_branch_length.png", plot = k, path = figures_dir, height = 5, width = 5)

#edge_time_vs_clock_rate_r_squared_df <- read_csv(file = "../analysis/stats/edge_time_vs_clock_rate_r_squared2.csv", col_names = FALSE)
#colnames(edge_time_vs_clock_rate_r_squared_df) <- c("r_squared", "dr", "model")
l <- ggplot(summary_stats, aes(x=as.numeric(dr)/2, y=r_squared, color = model)) +
    geom_point(size = 2) +
    geom_point(data = highlighted_sims, color = "darkgrey", size = 5, alpha = 0.5) +
    theme_classic() +
    scale_color_manual(values = sim_colors, labels = c("Boundary-driven", "Unrestricted")) +
    labs(color = "", x = "Cell turnover - P(cell death)", y = bquote(''~R^2*' clock rate and fraction lineage time on edge')) +
    theme(legend.position = c(0.8, 0.85),
          #legend.background = element_rect(linetype = "dashed", size = 0.5, colour = "darkgrey"),
          legend.title = element_blank(),
          text = element_text(size=15)) +
    guides(color = guide_legend(override.aes = list(size = 3) ) )
l
ggsave(file = "edge_time_vs_clock_rate_r_squared.png", plot = l, path = figures_dir, height = 5, width = 5)



