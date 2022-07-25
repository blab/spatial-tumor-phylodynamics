# Write_state_cllocks_xml_from_physicell.R

## Script to convert physicell simulations into state clocks xml

library(tumortree)
library(tidyverse)
library(rdist)


sample_physicell_diversified <- function(all_ids, tri, n_sampled_cells = 100) {
    
    #filter for single tumor to sample 
    trial_ids <- all_ids %>% 
        dplyr::filter(trial == tri)
    
    #randomize order of cell ids
    trial_ids <- trial_ids[sample(1:nrow(trial_ids), nrow(trial_ids), 
                                      replace = FALSE), ]
    
    #pairwise distance matrix between all alive cells in tumor
    dist_mat <- rdist::pdist(as.matrix(data.frame(x = trial_ids$x, 
                                                  y = trial_ids$y), ncol = 2))
    
    #find sampled set of cells that maximized distance between sampled cells
    fps <- rdist::farthest_point_sampling(dist_mat, k = n_sampled_cells)
    
    #label cells as sampled or unsampled 
    trial_ids$sampled <- FALSE
    trial_ids$sampled[fps] <- TRUE
    
    #filter to only sampled cells 
    sampled_trial_ids <- trial_ids %>% 
        dplyr::filter(sampled)
    
    return(sampled_trial_ids)
    
    
}
#set.seed(3123)

#apply diversified sampling to all tumors in dataset and combine into one dataset
# all_sampled_trials_df <- purrr::map(all_trials, function(tri) sample_physicell_diversified(all_ids = all_ids, tri = tri, n_sampled_cells = 100)) %>% 
#     bind_rows

#write_tsv(all_sampled_trials_df, file = "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/simulation_data/sampled_endstates.tsv")
#set.seed(381)
#all_cells <- read_tsv(data_files_2D_sel_bdg_diversified[[1]])


get_generated_sequence <- function(mut_row, all_sampled_muts, reference_bps, alternate_bps) {
    
    sequence <- c()
    for (pos in 1:length(all_sampled_muts)) {
        
        if(all_sampled_muts[pos] %in% mut_row) {
            
            sequence <- c(sequence, alternate_bps[pos])
            
        } else {
            sequence <- c(sequence, reference_bps[pos])
        }
    }
   
    return(paste(sequence, collapse =""))
}
draw_alt_bp <- function(ref, nucleotides = c("A", "G", "C", "T")) {
    
    seq_changes <- nucleotides[-which(nucleotides == ref)]
    alt_bp <- sample(seq_changes, size = 1)
    
    return(alt_bp)
}

write_state_clocks_xml_from_physicell_data <- function(all_cells, template_file, xml_file,
                                                       edge_cutoff=20, N=5000, center_edge = NULL, subsample = NULL) {
    
    sampled_cells <- tumortree::filter_alive_cells(all_cells) 
    
    if (! is.null(subsample)) {
        
        sampled_cells <- sampled_cells[sample(1:nrow(sampled_cells), size = subsample, replace = FALSE),]
    }
    #get sequence data from JC model (mutations are just recorded as integers)
    nucleotides <- c("A", "G", "C", "T")

    sampled_muts <- sampled_cells$mutations
    sampled_muts_list <- purrr::map(sampled_muts, 
                                function(mut_row) as.integer(tumortree::extract_mutations(mut_row)))

    all_sampled_muts <- sort(unique(unlist(sampled_muts_list)))
    reference_bps <- sample(nucleotides, size = length(all_sampled_muts), replace = TRUE)
    alternate_bps <- purrr::map_chr(reference_bps, draw_alt_bp)


    cell_sequences <- purrr::map_chr(sampled_muts_list, function(mut_row) get_generated_sequence(mut_row = mut_row, 
                                                                           all_sampled_muts = all_sampled_muts, 
                                                                           reference_bps = reference_bps, 
                                                                          alternate_bps = alternate_bps))

    # sampled_sim_tree <- tumortree::convert_all_cells_to_tree_fast(all_cells = all_cells, sampled_cells_indices = sampled_cells$index)
    sampled_tree_length <- all_cells %>%
        dplyr::mutate(br = deathdate - birthdate) %>%
        dplyr::summarize(tree_length = sum(br))

    
    #use this information to calculate clock rate

    
    clock_rate_mean <- length(unique(all_sampled_muts)) / sampled_tree_length /  length(all_sampled_muts)


    clock_rate_string <- format(clock_rate_mean, scientific = FALSE, digits = 5)
    #calculate clock rate

    clock_rate_mean <- mean(purrr::map_dbl(sampled_muts_list, length) / length(all_sampled_muts)) / max(all_cells$deathdate)
    clock_rate_string <- format(clock_rate_mean, scientific = FALSE, digits = 4)
    
    #origin string
    
    origin_string <- as.character(max(sampled_cells$deathdate))
    

    #get edge/center states

    sampled_cell_states <- as.integer(sampled_cells$distToHull < edge_cutoff)
    sampled_cells$cell_name <- paste0("cell", sampled_cells$index, 
                                  "loc", sampled_cell_states, sep = "")
    
    #rho string
    
    if (is.null(center_edge)) {
        
        rho_string <- format(nrow(sampled_cells) / N, digit = 2)
        
    } else {
        
        
        rho_string <- paste0(format(max(0, sum(1 - sampled_cell_states) / center_edge$center, na.rm = TRUE), digit = 3), 
                             " ",
                             format(max(0, sum(sampled_cell_states) / center_edge$edge, na.rm = TRUE), digit = 3),
                             sep = "")
        
        print(rho_string)
        
    }
    
    
    
    clock_rate_state_df <- data.frame("clock_rate" = purrr::map_dbl(sampled_muts_list, length) / length(all_sampled_muts) / max(all_cells$deathdate),
               "state" = sampled_cell_states) 
    
    # state_colors <- get_color_palette(names = c("edge", "center"))
    # g <- ggplot(clock_rate_state_df, aes(x=as.factor(state), y=clock_rate)) +
    #     geom_violin(aes(fill = ifelse(state == 1, "edge", "center"))) + theme_classic() + theme(legend.position = "none") +
    #     scale_fill_manual(values = state_colors) + xlab("State") + ylab("Clock rate")
    # 
    #print(g)
    
    sampled_cells$state <- sampled_cell_states
    #ggplot(sampled_cells, aes(x=locx, y=locy, color = state)) + geom_point() + theme_void() + theme(legend.position = "none")
    
    alignment_string <- ""
    for (i in 1:nrow(sampled_cells)) {
        alignment_string <- paste0(alignment_string, paste0("\t", 
                                                        "<sequence id=\"seq_", sampled_cells$cell_name[i], 
                                                        "\" spec=\"Sequence\" taxon=\"", sampled_cells$cell_name[i], 
                                                        "\" totalcount=\"4\" value=\"", cell_sequences[i], 
                                                        "\"/>", sep = ""), sep = "\n")
    }
    taxon_traits <- c()
    for (i in 1:nrow(sampled_cells)) {
        taxon_traits <- c(taxon_traits, paste0(sampled_cells$cell_name[i], 
                                           "=loc", sampled_cell_states[i], sep = ""))
    }
    taxon_traits_string <- paste(taxon_traits, collapse = ",")

    ##write xml from template
    con = file(template_file, "r")
    first_line = TRUE

    while (length(x <- readLines(con, n = 1)) > 0) {
    
        if (grepl("insert_sequence", x, fixed = TRUE)) {
        
            write(alignment_string, file = xml_file, append = !first_line)
        
        } else if (grepl("insert_trait", x, fixed = TRUE)) {
        
            write(gsub("insert_trait", taxon_traits_string, x), 
              file = xml_file, append = !first_line)
        
        } else if (grepl("insert_clock_rate", x, fixed = TRUE)) {
        
            write(gsub("insert_clock_rate", clock_rate_string, x), 
                  file = xml_file, append = !first_line)
        
        } else if (grepl("insert_origin", x, fixed = TRUE)) {
        
            write(gsub("insert_origin", origin_string, x), 
                  file = xml_file, append = !first_line)
        
        } else if (grepl("insert_rho", x, fixed = TRUE)) {
        
            write(gsub("insert_rho", rho_string, x), 
                  file = xml_file, append = !first_line)
        
        } else {
        
            write(x, file = xml_file, append = !first_line)
        }
        first_line = FALSE
    }
    close(con)
}

# write_state_clocks_xml_from_physicell_data(all_cells = all_cells, 
#                                            template_file = "outputs/beast_analysis/state_dependent_clock_model/validation/physicell/xml_files/state_dependent_canonical_estimate_dr_template.xml",
#                                            xml_file = "outputs/beast_analysis/state_dependent_clock_model/validation/physicell/xml_files/sampconfig_m0_w1.1_d0.8_t1_mg0.99_mm1_l2e+06_i1_s78461_state_dependent_canonical_estimate_dr.xml",
#                                            N = 5000, 
#                                            edge_cutoff = 20)


#wrapper function
write_state_clocks_xml_from_physicell_data_file <- function(data_file, 
                                                            template_file,
                                                            xml_file = NULL,
                                                            N = 5000, 
                                                            edge_cutoff = 20, 
                                                            sample_size_vary = FALSE, 
                                                            overwrite = FALSE, 
                                                            edge_center_proportions_df = NULL, 
                                                            center_edge = NULL, 
                                                            subsample = NULL) {
    if (is.null(xml_file)) {
        
        if (is.null(subsample)) {
            
            xml_file <- paste0(dirname(data_file), "/xml_files/", gsub(".tsv", "", basename(data_file)), ".xml", sep = "")
            
        } else {
            
            xml_file <- paste0(dirname(data_file), "/xml_files/", gsub(".tsv", "", basename(data_file)), "_n_", subsample, ".xml", sep = "")
        }
        
    }
    
    if (file.exists(xml_file) & (! overwrite)) {
        
        message("file exists, set overwrite=True if desired")
    } else {
        all_cells <- read_tsv(data_file)
    
        write_state_clocks_xml_from_physicell_data(all_cells = all_cells, 
                                               template_file = template_file,
                                               xml_file = xml_file,
                                               N = N, 
                                               edge_cutoff = edge_cutoff, 
                                               center_edge = center_edge, 
                                               subsample = subsample)
    }
    
}

template_file <- "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/state_clocks_template.xml"




data_files_power <- paste0("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/simulated_data/2D_neut_bdg/varyN_d0.3/to_beast_format/",
                     list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/simulated_data/2D_neut_bdg/varyN_d0.3/to_beast_format",
                                pattern="sampconfig.*"))

edge_center_proportions_2D_neut_bdg_file <- "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/simulated_data/2D_neut_bdg/edge_v_center_counts_2D_neut_bdg.tsv"
edge_center_proportions_2D_neut_bdg_df <- read_tsv(edge_center_proportions_2D_neut_bdg_file)

set.seed(7717)
for (file in data_files_power) {
    
    center_edge <- edge_center_proportions_2D_neut_bdg_df[gsub(".tsv", "", edge_center_proportions_2D_neut_bdg_df$file) == gsub("_diversified\\s*(.*?)\\.tsv", "", basename(file), perl = TRUE), c("center", "edge")]
    write_state_clocks_xml_from_physicell_data_file(data_file = file, 
                                                    template_file = template_file, 
                                                    overwrite = TRUE, 
                                                    N = 10000, 
                                                    center_edge = center_edge)
    
}

data_files_2D_neut_bdg_diversified <- paste0("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/simulated_data/2D_neut_bdg/diversified_100/to_beast_format_corrected_hulls/",
                                list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/simulated_data/2D_neut_bdg/diversified_100/to_beast_format_corrected_hulls",
                                           pattern="sampconfig.*"))



set.seed(3713)
for (file in data_files_2D_neut_bdg_diversified) {
    
    print(file)
    center_edge <- edge_center_proportions_2D_neut_bdg_df[gsub(".tsv", "",
                                                               edge_center_proportions_2D_neut_bdg_df$file) == gsub(".tsv", "", gsub("_diversified\\s*(.*?)\\.tsv", "", basename(file), perl = TRUE)),
                                                          c("center", "edge")]
    
    write_state_clocks_xml_from_physicell_data_file(data_file = file, 
                                                    template_file = template_file, 
                                                    overwrite = FALSE,
                                                    N = 10000, 
                                                    center_edge = center_edge)
    
}

#subsampled diversified 
subsampled_ns <- c(25,  50, 75)
for (file in data_files_2D_neut_bdg_diversified) {
    
    
    print(file)
    center_edge <- edge_center_proportions_2D_neut_bdg_df[gsub(".tsv", "",
                                                               edge_center_proportions_2D_neut_bdg_df$file) == gsub(".tsv", "", gsub("_diversified\\s*(.*?)\\.tsv", "", basename(file), perl = TRUE)),
                                                          c("center", "edge")]
    
    for (n_subsampled in subsampled_ns) {
        write_state_clocks_xml_from_physicell_data_file(data_file = file, 
                                                    template_file = template_file, 
                                                    overwrite = FALSE,
                                                    N = 10000, 
                                                    center_edge = center_edge, 
                                                    subsample = n_subsampled)
    }
    
}


data_files_2D_neut_bdg_vary_mu_diversified <- paste0("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/simulated_data/2D_neut_bdg/varyMu_d0.3/to_beast_format_corrected_hull/",
                                             list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/simulated_data/2D_neut_bdg/varyMu_d0.3/to_beast_format_corrected_hullls",
                                                        pattern="sampconfig.*"))



set.seed(4221)
for (file in data_files_2D_neut_bdg_vary_mu_diversified) {
    
    print(file)
    center_edge <- edge_center_proportions_2D_neut_bdg_df[gsub(".tsv", "",
                                                               edge_center_proportions_2D_neut_bdg_df$file) == gsub(".tsv", "", gsub("_diversified\\s*(.*?)\\.tsv", "", basename(file), perl = TRUE)),
                                                          c("center", "edge")]
    
    write_state_clocks_xml_from_physicell_data_file(data_file = file, 
                                                    template_file = template_file, 
                                                    overwrite = FALSE,
                                                    N = 10000, 
                                                    center_edge = center_edge)
    
}

#subsampled

set.seed(3312)
for (file in data_files_2D_neut_bdg_vary_mu_diversified) {
    
    print(file)
    center_edge <- edge_center_proportions_2D_neut_bdg_df[gsub(".tsv", "",
                                                               edge_center_proportions_2D_neut_bdg_df$file) == gsub(".tsv", "", gsub("_diversified\\s*(.*?)\\.tsv", "", basename(file), perl = TRUE)),
                                                          c("center", "edge")]
    
    for (n_subsampled in subsampled_ns) {
        write_state_clocks_xml_from_physicell_data_file(data_file = file, 
                                                        template_file = template_file, 
                                                        overwrite = FALSE,
                                                        N = 10000, 
                                                        center_edge = center_edge, 
                                                        subsample = n_subsampled)
    }
    
}


data_files_2D_sel_diversified <- paste0("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/simulated_data/2D_sel_highermu/diversified_100/to_beast_format_corrected_hulls/",
                                            list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/simulated_data/2D_sel_highermu/diversified_100/to_beast_format_corrected_hulls",
                                                       pattern="sampconfig.*"))

edge_center_proportions_2D_sel_df <-read_tsv("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/simulated_data/2D_sel_highermu/edge_v_center_counts_2D_sel_highermu.tsv")

set.seed(3713)
for (file in data_files_2D_sel_diversified) {
    
    center_edge <- edge_center_proportions_2D_sel_df[gsub(".tsv", "", edge_center_proportions_2D_sel_df$file) == gsub(".tsv", "", gsub("_diversified\\s*(.*?)\\.tsv", "", basename(file), perl = TRUE)), c("center", "edge")]
    
    write_state_clocks_xml_from_physicell_data_file(data_file = file, 
                                                    template_file = template_file, 
                                                    overwrite = FALSE, 
                                                    center_edge = center_edge, 
                                                    N = 10000)
    
}

set.seed(3121)
for (file in data_files_2D_sel_diversified) {
    
    center_edge <- edge_center_proportions_2D_sel_df[gsub(".tsv", "", edge_center_proportions_2D_sel_df$file) == gsub(".tsv", "", gsub("_diversified\\s*(.*?)\\.tsv", "", basename(file), perl = TRUE)), c("center", "edge")]
    
    for (n_subsampled in subsampled_ns) {
        write_state_clocks_xml_from_physicell_data_file(data_file = file, 
                                                        template_file = template_file, 
                                                        overwrite = FALSE,
                                                        N = 10000, 
                                                        center_edge = center_edge, 
                                                        subsample = n_subsampled)
    }
    
}

data_files_2D_sel_bdg_diversified <- paste0("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/simulated_data/2D_sel_bdg_highermu/diversified_100/to_beast_format_corrected_hulls/",
                                              list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/simulated_data/2D_sel_bdg_highermu/diversified_100/to_beast_format_corrected_hulls",
                                                         pattern="sampconfig.*"))
edge_center_proportions_2D_sel_bdg_df <-read_tsv("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/simulated_data/2D_sel_bdg_highermu/edge_v_center_counts_2D_sel_bdg_highermu.tsv")

set.seed(323)
for (file in data_files_2D_sel_bdg_diversified) {
    
    center_edge <- edge_center_proportions_2D_sel_bdg_df[gsub(".tsv", "", edge_center_proportions_2D_sel_bdg_df$file) == gsub(".tsv", "", gsub("_diversified\\s*(.*?)\\.tsv", "", basename(file), perl = TRUE)), c("center", "edge") ]
    
    write_state_clocks_xml_from_physicell_data_file(data_file = file, 
                                                    template_file = template_file, 
                                                    center_edge = center_edge,
                                                    overwrite = FALSE)
    
}

set.seed(3323)
for (file in data_files_2D_sel_bdg_diversified) {
    
    center_edge <- edge_center_proportions_2D_sel_bdg_df[gsub(".tsv", "", edge_center_proportions_2D_sel_bdg_df$file) == gsub(".tsv", "", gsub("_diversified\\s*(.*?)\\.tsv", "", basename(file), perl = TRUE)), c("center", "edge") ]
    
    for (n_subsampled in subsampled_ns) {
        write_state_clocks_xml_from_physicell_data_file(data_file = file, 
                                                        template_file = template_file, 
                                                        overwrite = FALSE,
                                                        N = 10000, 
                                                        center_edge = center_edge, 
                                                        subsample = n_subsampled)
    }
    
}

data_files_3D_neut_bdg_diversified <- paste0("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/simulated_data/3D_neut_bdg/diversified_100/to_beast_format_corrected_hulls/",
                                            list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/simulated_data/3D_neut_bdg/diversified_100/to_beast_format_corrected_hulls",
                                                       pattern="sampconfig.*"))

edge_center_proportions_3D_neut_bdg_df <-read_tsv("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/physicell/simulated_data/3D_neut_bdg/edge_v_center_counts_3D_neut_bdg.tsv")

set.seed(2331)
for (file in data_files_3D_neut_bdg_diversified) {
    
    center_edge <- edge_center_proportions_3D_neut_bdg_df[gsub(".tsv", "", edge_center_proportions_3D_neut_bdg_df$file) == gsub(".tsv", "", gsub("_diversified\\s*(.*?)\\.tsv", "", basename(file), perl = TRUE)), c("center", "edge")]
    
    write_state_clocks_xml_from_physicell_data_file(data_file = file, 
                                                    template_file = template_file, 
                                                    center_edge = center_edge, 
                                                    overwrite = FALSE, 
                                                    N = 15000)
    
}

set.seed(23211)
for (file in data_files_3D_neut_bdg_diversified) {
    
    center_edge <- edge_center_proportions_3D_neut_bdg_df[gsub(".tsv", "", edge_center_proportions_3D_neut_bdg_df$file) == gsub(".tsv", "", gsub("_diversified\\s*(.*?)\\.tsv", "", basename(file), perl = TRUE)), c("center", "edge")]
    
    for (n_subsampled in subsampled_ns) {
        write_state_clocks_xml_from_physicell_data_file(data_file = file, 
                                                        template_file = template_file, 
                                                        overwrite = FALSE,
                                                        N = 10000, 
                                                        center_edge = center_edge, 
                                                        subsample = n_subsampled)
    }
    
}






#
# for (file in data_files) {
#     
#     xml_file <- paste0("outputs/beast_analysis/state_dependent_clock_model/validation/physicell/xml_files/", gsub(".tsv", "", basename(file)), ".xml")
#     
#     write_state_clocks_xml_from_physicell_data_file(data_file = file, 
#                                                     template_file = "outputs/beast_analysis/state_dependent_clock_model/validation/physicell/xml_files/state_dependent_canonical_estimate_dr_template.xml",
#                                                     xml_file = xml_file)
#     
# }




# all_cells <- read_tsv(data_files_power[1])
# 
# t <- 1:max(all_cells$deathdate)
# b <- 0.0008
# N_final <- 10000
# n <- exp(b*t)
# ggplot(data = data.frame("t" = t, "n" = n), aes(x=t, y=n)) + geom_line() + theme_classic() + geom_hline(yintercept = N_final, color = "red", linetype = "dashed")
