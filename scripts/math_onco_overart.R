# Figures for Math Onco Coverart
library(tumortree)
library(tidyverse)
library(ggtree)

setwd("~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics")
all_cells <- read.csv("~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/simulation_data/cells_death_rate_validation_pop_15000_dr_0.10.csv")

endpoint <- max(all_cells$deathdate)

all_cells$n_divisions <- purrr::map_dbl(all_cells$index, function(i) get_n_divisions(index=i, all_cells=all_cells))
#write.csv(all_cells, "~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/simulation_data/cells_death_rate_validation_pop_15000_dr_0.10_with_division.csv")

alive_cells <- all_cells %>% 
    dplyr::filter(deathdate == endpoint)

#alive_cells <- tumortree::mark_boundary(cells_to_mark = alive_cells, alive_cells)

#alive_cells$n_divisions <- purrr::map_dbl(alive_cells$index, function(i) get_n_divisions(index=i, all_cells=all_cells))

set.seed(44411)
sampled_cells <- alive_cells %>% 
    tumortree::sample_alive_cells(., n = 100) %>% 
    dplyr::filter(sampled)



get_cell_paths <- function(index, all_cells) {
    
    ancestral_cells <- collect_all_ancestors(index = index, all_cells = all_cells)
    
    df <- data.frame("locx" = all_cells$locx[match(ancestral_cells, all_cells$index)],
                     "locy" = all_cells$locy[match(ancestral_cells, all_cells$index)],
                     "deathdate" = all_cells$deathdate[match(ancestral_cells, all_cells$index)],
                     "birthdate" = all_cells$birthdate[match(ancestral_cells, all_cells$index)],
                     "index" = ancestral_cells,
                     "sampled_index" = index)
    
    return(df)

}
sampled_ancestral_paths <- purrr::map(sampled_cells$index, function(i) get_cell_paths(index = i, all_cells = all_cells)) %>% 
    bind_rows()


tumor_sun <- ggplot(alive_cells, aes(x=locx, y=locy, color= n_divisions)) +
    geom_point(size =2) + theme_void() + scale_color_gradient(low = "#15948b", high = "#ffff00", limits = c(0, max(alive_cells$n_divisions))) +
    #theme(plot.background = element_rect(fill = "black")) +
    theme(legend.position = "none")+ geom_path(data = sampled_ancestral_paths, aes(group=sampled_index), color = "white")   
tumor_sun 

ggsave("~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/figures/tumor_sun.png", tumor_sun, height = 5, width = 5)

tumor_sun_ylims <- layer_scales(tumor_sun)$y$get_limits()
tumor_sun_xlims <- layer_scales(tumor_sun)$x$get_limits()
# genetic_tree <- convert_all_cells_to_tree_fast(all_cells,
#                                                add_all_states =  FALSE,
#                                                branch_unit = "genetic")

# time_tree <- convert_all_cells_to_tree_fast(all_cells,
#                                                add_all_states =  FALSE,
#                                                branch_unit = "time",
#                                             sampled_cells_indices=sampled_cells$index)

genetic_tree <- convert_all_cells_to_tree_fast(all_cells,
                                            add_all_states =  FALSE,
                                            branch_unit = "divisions",
                                            sampled_cells_indices=sampled_cells$index)

pruned_genetic_tree <- prune_simulated_tree(genetic_tree, sampled_cells_indices=sampled_cells$index,
                                             branch_unit = "divisions",
                                            all_cells = all_cells)
#"#55a791"
tree_plot <-ggtree(pruned_genetic_tree, layout = "ellipse",
                    aes(color = divisions + pruned_genetic_tree@phylo$root.edge, size = -4*divisions)) +
    geom_rootedge(color = "#6aae8b", size = 3) +
    scale_color_gradient(low = "#55a791", high = "#ffff00", limits = c(0, max(pruned_genetic_tree@data$divisions + pruned_genetic_tree@phylo$root.edge, na.rm = TRUE))) +
    geom_tippoint(size = 3)  +
    theme_void() +
    theme(plot.background = element_rect(fill = "black")) + theme(legend.position = "none")+ coord_cartesian(clip = "off") +
    scale_size_continuous(range = c(0.1, 2)) + coord_flip()

tree_plot 

ggsave("~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/figures/artsy_tree.png", tree_plot , height = 5, width = 12)


library(ggpubr)
library(cowplot)
##### ANIMATIONS
set.seed(2123)
endpoints <- seq(30, max(all_cells$deathdate))

for (e in endpoints) {


    alive_cells_step <- all_cells %>% 
        dplyr::filter(deathdate >= e & birthdate < e)



    #alive_cells_step$n_divisions <- purrr::map_dbl(alive_cells_step$index, function(i) get_n_divisions(index=i, all_cells=all_cells))

    sampled_ancestral_paths_step <- sampled_ancestral_paths %>% 
        filter(birthdate < e)

    #upper_color <- lighten("#ecf42f", amount = 0.1)
    upper_color <- lighten("#ffff00", amount = 0.1)
    #upper_color <- "#ffff00"
    tumor_sun_step <- ggplot(alive_cells_step, aes(x=locx, y=locy, color= n_divisions)) +
        geom_point(size=0.5) + theme_void() +
        scale_color_gradient(low = "#15948b", high = upper_color,
                             limits = c(0, max(alive_cells$n_divisions))) +
        theme(plot.background = element_rect(fill = "black")) +
        theme(legend.position = "none")+ geom_path(data = sampled_ancestral_paths_step, aes(group=sampled_index), color = "white", size = 0.25) +
        xlim(tumor_sun_xlims) + ylim(tumor_sun_ylims)
    tumor_sun_step 
    
    # tumor_sun_step2 <- ggplot(alive_cells_step, aes(x=locx, y=locy, color= n_divisions)) +
    #     geom_point(size=1) + theme_void() +
    #     scale_color_gradient(low = "#15948b", high = "#ffff00",
    #                          limits = c(0, max(alive_cells$n_divisions))) +
    #     #theme(plot.background = element_rect(fill = "black")) +
    #     theme(legend.position = "none")+ geom_path(data = sampled_ancestral_paths_step, aes(group=sampled_index), color = "white", size = 0.25) +
    #     xlim(tumor_sun_xlims) + ylim(tumor_sun_ylims)
    # tumor_sun_step2 
    #ggsave(paste0("~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/figures/tumor_sun", e, ".png"), tumor_sun_step, height = 5, width = 5)

    


    trimmed_tree <- pruned_genetic_tree 
    
    trimmed_tree@data$divisions[which(trimmed_tree@data$birthdate > e)] <- NA
    
    
    growing_tree <- pruned_genetic_tree
    # 
    # growing_tree@data$divisions[which(growing_tree@data$birthdate > e+1)] <- NA
    # step_tips <- growing_tree@data$index[which(growing_tree@data$birthdate > e+1)]
    
    step_tips <- sampled_ancestral_paths %>%
        filter(birthdate < e) %>%
        group_by(sampled_index) %>%
        summarise(birthdate= max(birthdate, rm.na = TRUE)) %>%
        left_join(., sampled_ancestral_paths_step)
    
    excluded_tips <- sampled_ancestral_paths %>%
        filter((birthdate < e) & (! index %in% step_tips$index))
    
    #potential_children <- match(step_tips$index, growing_tree@data$parent_index)
    potential_children <- growing_tree@data$index[which(growing_tree@data$parent_index %in% step_tips$index)]
    true_tips <- step_tips$index[which(! step_tips$index %in% growing_tree@data$parent_index)]
    
    children <- c(potential_children, true_tips)
    children <- children[which(! children %in% excluded_tips$index)]
    #children <- children[which(! children %in% growing_tree@data$parent_index[which(growing_tree@data$index %in% children)])]
    #true_tips <- step_tips$index[which(is.na(potential_children))]
    #children_with_parents <-  growing_tree@data$index[potential_children[which(! is.na(potential_children))]]
    
    selected_indices <- children
    
    #growing_tree@data$divisions[! which(growing_tree@data$index %in% selected_indices)] <- NA

    
    tree_plot_step <- ggtree(trimmed_tree,layout = "ellipse",
                        aes(color = divisions + pruned_genetic_tree@phylo$root.edge,
                            size = -divisions)) +
        geom_rootedge(color = "#6aae8b", size = 1) +
        scale_color_gradient(low = "#55a791", high = upper_color, na.value = "black",
                             limits = c(0, max(pruned_genetic_tree@data$divisions + pruned_genetic_tree@phylo$root.edge, na.rm = TRUE))) +
        #geom_tippoint(size = 1)  +
        geom_point2(data=growing_tree, aes(subset = index %in% selected_indices, color = divisions), size = 1) +
        theme_void() +
        theme(plot.background = element_rect(fill = "black")) + theme(legend.position = "none")+ coord_cartesian(clip = "off") +
        #scale_size_continuous(range = c(0.01, 0.7), limits = c(-max(all_cells$n_divisions),0)) +
        scale_size_continuous(range = c(0.01, 0.7)) +
        coord_flip()
    tree_plot_step 
    
    # tree_plot_ste2p <- ggtree(trimmed_tree,layout = "ellipse",
    #                          aes(color = divisions + pruned_genetic_tree@phylo$root.edge,
    #                              size = -divisions)) +
    #     geom_rootedge(color = "#6aae8b", size = 1) +
    #     scale_color_gradient(low = "#55a791", high = "#ffff00", na.value = "black",
    #                          limits = c(0, max(pruned_genetic_tree@data$divisions + pruned_genetic_tree@phylo$root.edge, na.rm = TRUE))) +
    #     #geom_tippoint(size = 1)  +
    #     geom_point2(data=growing_tree, aes(subset = index %in% selected_indices, color = divisions), size = 1) +
    #     theme_void() +
    #     #theme(plot.background = element_rect(fill = "black")) +
    #     theme(legend.position = "none")+ coord_cartesian(clip = "off") +
    #     scale_size_continuous(range = c(0.03, 0.65), limits = c(-max(all_cells$n_divisions),0)) + coord_flip()
    # 
    # 
    # tree_plot_step2

    #ggsave(paste0("~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/figures/artsy_tree", e, ".png"), tree_plot_step , height = 5, width = 12)
    
    w = 20.81*2
    h = 13.90*2
    
    comb_plot <- ggdraw(xlim = c(0,w), ylim = c(0,h)) +
        draw_plot(tree_plot_step, x = 0.05*w, y = 0, width = 0.9*w, height = 0.5*h) +
        draw_plot(tumor_sun_step, x = 0.5*w - 0.25*h, y = 0.5*h, width = 0.5*h, height = 0.5*h) +
        theme(plot.background = element_rect(fill = "black"))

    comb_plot 
    # comb_plot2 <- ggdraw(xlim = c(0,w), ylim = c(0,h)) +
    #     draw_plot(tree_plot_step2, x = 0.05*w, y = 0, width = 0.9*w, height = 0.5*h) +
    #     draw_plot(tumor_sun_step, x = 0.5*w - 0.45*h, y = 0.53*h, width = 0.45*h, height = 0.7*h) +
    #     theme(plot.background = element_rect(fill = "black"))
    # comb_plot2 
    ggsave(paste0("~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/figures/comb_artsy_tree_tumor_sun", e, ".png"), comb_plot,
           height = 1390, width = 2081,
           units = "px")
    
#     ggsave(paste0("~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/figures/comb_artsy_tree_tumor_sun2", e, ".png"), comb_plot2,
#            height = 1390, width = 2081,
#            units = "px")
}


 
