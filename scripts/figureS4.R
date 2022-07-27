#figureS4.R

##Scripts to make figure S4

#example random sampling diagram
example_cells <- read_csv("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation/cells_death_rate_validation_pop_1000_dr_0.27.csv") %>% 
    tumortree::filter_alive_cells()
example_cells <- mark_boundary(cells_to_mark = example_cells, alive_cells = example_cells)

#Write example tumor to CSV
write_csv(example_cells, "../eden/simulation_data/alive_cells_death_rate_validation_pop_1000_dr_0.27.csv")

#To skip to visualization
example_cells <- read_csv(file = "../eden/simulation_data/alive_cells_death_rate_validation_pop_1000_dr_0.27.csv")

set.seed(392)
example_cells_diversified <- example_cells %>% 
    tumortree::sample_alive_cells(., n =100, diversified_sampling = TRUE) %>% 
    dplyr::filter(sampled)

example_cells_random <- example_cells %>% 
    tumortree::sample_alive_cells(., n =100, diversified_sampling = FALSE) %>% 
    dplyr::filter(sampled)

colors_edge_center <- get_color_palette(names = c("edge", "center"))

example_tumor_inset_diversified <- ggplot(example_cells, aes(locx, locy, color = ifelse(est_edge == 1, "edge", "center"))) +
    geom_point() + theme_void() +
    #geom_point(data = example_cells_diversified, shape = 21, aes(fill = ifelse(est_edge == 1, "edge", "center")), color = "black", size =3)+
    geom_point(data = example_cells_diversified, shape = 1, fill = "black",  color = "black", alpha = 1, size = 3)+
    scale_fill_manual(values = colors_edge_center) +
    scale_color_manual(values = colors_edge_center) +
    theme(legend.position = "none") 

example_tumor_inset_diversified
ggsave(file = "../figures/diversified_sampling_example_tumor.png", example_tumor_inset_diversified, height = 5, width = 5)    


example_tumor_inset_random <- ggplot(example_cells, aes(locx, locy, color = ifelse(est_edge == 1, "edge", "center"))) +
    geom_point() + theme_void() +
    #geom_point(data = example_cells_diversified, shape = 21, aes(fill = ifelse(est_edge == 1, "edge", "center")), color = "black", size =3)+
    geom_point(data = example_cells_random, shape = 1, fill = "black",  color = "black", alpha = 1, size = 3)+
    scale_fill_manual(values = colors_edge_center) +
    scale_color_manual(values = colors_edge_center) +
    theme(legend.position = "none") 

example_tumor_inset_random
ggsave(file = "../figures/random_sampling_example_tumor.png", example_tumor_inset_random, height = 5, width = 5)    
