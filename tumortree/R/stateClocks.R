#' Write BEAST StateClocks XML
#'
#' Insert sequence and trait data from simulated tumor cells into template xml to run BEAST StateClocks
#' Template must include "insert_sequence" to indicate spot for alignment data and "insert_trait" to indicate location for cell states.
#'
#' @param sim_cells_file string File for csv of all simulated cells
#' @param sampled_cell_file string File for csv of only sampled cells. If doesn't exist, cells will be sampled and newly sampled cells will be written to this file.
#' @param xml_file string File to write new XML file
#' @param template_file string Filename of XML template that includes insert_sequence and insert_trait placeholders
#' @param rewrite_sampled_cells_csv logical FALSE. Rewrite sampled cells file with edge/center info
#' @param n_samples integer default 100. Number of cells to sample if sampling.
#' @param diversified logical default FALSE. If sampling new cells should be diversified where distance is maximized between sampled cells.
#' @importFrom purrr map_chr
#' @importFrom magrittr %>%
#' @export
#'
#'
#'
write_state_clocks_xml <-function(sim_cells_file, sampled_cells_file, xml_file, template_file, 
          rewrite_sampled_cells_csv = FALSE, n_samples = 100, diversified = FALSE, 
          resample = FALSE, alive_cells_only = FALSE) {
    print("Starting")
    
    #Check if alive cells already saved
    alive_cells_file <- gsub("cells", "alive_cells", sim_cells_file)
    
    if (file.exists(alive_cells_file) & alive_cells_only) {

	    print("Found alive cells")
	    alive_cells <- read.csv(alive_cells_file)

    } else {
	
	print("Filtering alive cells")
    	sim_cells <- read.csv(sim_cells_file) %>% normalize_locs
    	endpoint <- max(sim_cells$deathdate)
    	alive_cells <- sim_cells %>% dplyr::filter(deathdate == endpoint)
    	assert(all(alive_cells$deathdate == endpoint))
    	print(paste0("Number of alive cells: ", nrow(alive_cells), sep = ""))
    }

    if (file.exists(sampled_cells_file) & (!resample)) {
      
	    sampled_cells <- read.csv(sampled_cells_file) %>% normalize_locs
    
    } else {
    
	message("Sampling")
        sampled_cells <- sample_alive_cells(alive_cells = alive_cells, 
                                            n = n_samples, diversified_sampling = diversified) %>% 
            normalize_locs %>% filter(sampled == TRUE)

    	write_csv(sampled_cells, file = sampled_cells_file)
    }
    #sampled_sim_tree <- tumortree::convert_all_cells_to_tree_fast(all_cells = sim_cells, 
#                                                                  sampled_cells_indices = sampled_cells$index)
#    sampled_tree_length <- sum(sampled_sim_tree@phylo$edge.length)
#    sim_tree <- tumortree::convert_all_cells_to_tree_fast(all_cells = sim_cells)
#    tree_length <- sum(sim_tree@phylo$edge.length)

    print("Collapsing sequences")
    sampled_cells$sequence_collapsed <- purrr::map_chr(sampled_cells$sequence,function(seq) extract_sequences(seq))
    
    print("Calculating tree length")

    if (exists("sim_cells")) {
    	all_sim_tree_length <- sim_cells %>% dplyr::mutate(br = deathdate - 
                                                           birthdate) %>% dplyr::summarize(tree_length = sum(br))
    	all_muts <- as.integer(unique(unlist(purrr::map(sim_cells$mutations, 
                                                    function(mut_row) extract_mutations(mut_row)))))
#    all_alive_muts <- as.integer(unique(unlist(purrr::map(alive_cells$mutations, 
#                                                          function(mut_row) extract_mutations(mut_row)))))
#    all_sampled_muts <- as.integer(unique(unlist(purrr::map(sampled_cells$mutations, 
#                                                            function(mut_row) extract_mutations(mut_row)))))
    	clock_rate_mean <- length(unique(all_muts))/all_sim_tree_length/length(all_muts)

    } else {
	
	muts_per_lineage <- (nchar(alive_cells$mutations) - 1)/2
        print(muts_per_lineage)
        seq_length <- nchar(sampled_cells$sequence_collapsed[1])
	print(seq_length)
     	clock_rate_per_lineage <- muts_per_lineage / alive_cells$deathdate / seq_length
	clock_rate_mean <- mean(clock_rate_per_lineage)
    }

    clock_rate_string <- format(clock_rate_mean, scientific = FALSE,
                                digits = 5)
    print(paste0("Clock rate: ", clock_rate_string, sep = ""))

    #All sequences should be the same length
    assert(all(nchar(sampled_cells$sequence_collapse) == nchar(sampled_cells$sequence_collapse[1])))
     
    print("Labeling edge states")
    if ((!"edge_adjacent" %in% colnames(sampled_cells)) | (!"most_extreme" %in% 
                                                           colnames(sampled_cells))) {
        sampled_cells$most_extreme <- purrr::map2_dbl(sampled_cells$locx, 
                                                      sampled_cells$locy, function(x_loc, y_loc) check_if_cell_most_extreme(x_loc, 
                                                                                                                            y_loc, alive_cells = alive_cells))
        alive_cells$most_extreme <- purrr::map2_dbl(alive_cells$locx, 
                                                    alive_cells$locy, function(x_loc, y_loc) check_if_cell_most_extreme(x_loc, 
                                                                                                                        y_loc, alive_cells = alive_cells))
        alive_cells$edge_adjacent <- purrr::map2_dbl(alive_cells$locx, 
                                                     alive_cells$locy, function(x_loc, y_loc) check_boundary_adjacent(x_loc, 
                                                                                                                      y_loc, alive_cells = alive_cells))
        sampled_cells$edge_adjacent <- purrr::map2_dbl(sampled_cells$locx, 
                                                       sampled_cells$locy, function(x_loc, y_loc) check_boundary_adjacent(x_loc, 
                                                                                                                          y_loc, alive_cells = alive_cells))
    }
    if (rewrite_sampled_cells_csv) {
        write_csv(sampled_cells, sampled_cells_file)
    }


    sampled_cell_states <- as.integer(sampled_cells$most_extreme | 
                                          sampled_cells$edge_adjacent)
    sampled_cells$cell_name <- paste0("cell", sampled_cells$index, 
                                      "loc", sampled_cell_states, sep = "")
    print("Writing alignment string")
    alignment_string <- ""
    for (i in 1:nrow(sampled_cells)) {
        alignment_string <- paste0(alignment_string, paste0("\t", 
                                                            "<sequence id=\"seq_", sampled_cells$cell_name[i], 
                                                            "\" spec=\"Sequence\" taxon=\"", sampled_cells$cell_name[i], 
                                                            "\" totalcount=\"4\" value=\"", sampled_cells$sequence_collapsed[i], 
                                                            "\"/>", sep = ""), sep = "\n")
    }
    taxon_traits <- c()
    for (i in 1:nrow(sampled_cells)) {
        taxon_traits <- c(taxon_traits, paste0(sampled_cells$cell_name[i], 
                                               "=loc", sampled_cell_states[i], sep = ""))
    }
    taxon_traits_string <- paste(taxon_traits, collapse = ",")
    rho_edge <- round(sum(sampled_cell_states)/sum(alive_cells$most_extreme | 
                                                       alive_cells$edge_adjacent), 5)
    rho_center <- round(sum(!sampled_cell_states)/sum(!(alive_cells$most_extreme | 
                                                            alive_cells$edge_adjacent)), 5)
    rho_string <- paste0(rho_center, " ", rho_edge, sep = "")

    print("Writing XML")
    con = file(template_file, "r")
    first_line = TRUE
    while (length(x <- readLines(con, n = 1)) > 0) {
        if (grepl("insert_sequence", x, fixed = TRUE)) {
            write(alignment_string, file = xml_file, append = !first_line)
        }
        else if (grepl("insert_trait", x, fixed = TRUE)) {
            write(gsub("insert_trait", taxon_traits_string, x), 
                  file = xml_file, append = !first_line)
        }
        else if (grepl("insert_rho", x, fixed = TRUE)) {
            write(gsub("insert_rho", rho_string, x), file = xml_file, 
                  append = !first_line)
        }
        else if (grepl("insert_clock_rate", x, fixed = TRUE)) {
            write(gsub("insert_clock_rate", clock_rate_string, 
                       x), file = xml_file, append = !first_line)
        }
        else {
            write(x, file = xml_file, append = !first_line)
        }
        first_line = FALSE
    }
    close(con)
}
