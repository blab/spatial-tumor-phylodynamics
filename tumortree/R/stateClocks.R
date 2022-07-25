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
write_state_clocks_xml <- function(sim_cells_file,
                                   sampled_cells_file,
                                   xml_file,
                                   template_file,
                                   rewrite_sampled_cells_csv = FALSE,
                                   n_samples = 100,
                                   diversified = FALSE,
                                   resample = FALSE) {

    sim_cells <- read.csv(sim_cells_file) %>%
        normalize_locs


    endpoint <- max(sim_cells$deathdate)

    alive_cells <- sim_cells %>%
        dplyr::filter(deathdate == endpoint)

    if (file.exists(sampled_cells_file) & (! resample)) {

        sampled_cells <- read.csv(sampled_cells_file) %>%
            normalize_locs

    } else {

        message("Sampling")

        sampled_cells <- sample_alive_cells(alive_cells = alive_cells,
                                                       n = n_samples,
                                                       diversified_sampling = diversified) %>%
            normalize_locs %>%
            filter(sampled == TRUE)


        write_csv(sampled_cells, file = sampled_cells_file)
    }


    ## find unique sampled mutations to filter sequence data
    ## (only include segregating sites)
    if (! "sequence_collapsed" %in% colnames(sampled_cells)) {
        sampled_muts <- sampled_cells$mutations

        all_sampled_muts <- as.integer(unique(unlist(purrr::map(sampled_muts, function(mut_row) extract_mutations(mut_row)))))
        ## extract sequence data


        sampled_cells$sequence_collapsed <- purrr::map_chr(sampled_cells$sequence, function(seq) extract_sequences(seq, sampled_muts = all_sampled_muts))
    }
    ####################
    ## Find edge states

    if ((! "edge_adjacent" %in% colnames(sampled_cells)) | (! "most_extreme" %in% colnames(sampled_cells)) ) {
        sampled_cells$most_extreme <- purrr::map2_dbl(sampled_cells$locx, sampled_cells$locy,
                                                      function(x_loc, y_loc) check_if_cell_most_extreme(x_loc, y_loc, alive_cells = alive_cells))

        alive_cells$most_extreme <- purrr::map2_dbl(alive_cells$locx, alive_cells$locy,
                                                    function(x_loc, y_loc) check_if_cell_most_extreme(x_loc, y_loc, alive_cells = alive_cells))

        sampled_cells$edge_adjacent <- purrr::map2_dbl(sampled_cells$locx, sampled_cells$locy,
                                                       function(x_loc, y_loc) check_boundary_adjacent(x_loc, y_loc, alive_cells = alive_cells))
    }

    if(rewrite_sampled_cells_csv) {

        write_csv(sampled_cells, sampled_cells_file)
    }

    sampled_cell_states <- as.integer(sampled_cells$most_extreme | sampled_cells$edge_adjacent)

    ## Add edge state to name

    sampled_cells$cell_name <- paste0("cell", sampled_cells$index, "loc", sampled_cell_states, sep = "")

    alignment_string <- ""

    for (i in 1:nrow(sampled_cells)) {
        alignment_string <- paste0(alignment_string, paste0('\t', '<sequence id="seq_',
                                                            sampled_cells$cell_name[i],
                                                            '" spec="Sequence" taxon="',
                                                            sampled_cells$cell_name[i],
                                                            '" totalcount="4" value="',
                                                            sampled_cells$sequence_collapsed[i],
                                                            '"/>', sep = ""), sep = "\n")
    }

    taxon_traits <- c()

    for (i in 1:nrow(sampled_cells)) {

        taxon_traits <- c(taxon_traits, paste0(sampled_cells$cell_name[i], "=loc", sampled_cell_states[i], sep = ""))
    }

    taxon_traits_string <- paste(taxon_traits, collapse = ",")


    con = file(template_file, "r")

    first_line = TRUE

    while(length(x <- readLines(con, n = 1)) > 0) {



        if (grepl("insert_sequence", x, fixed = TRUE)) { #alignment insertion
            write(alignment_string, file=xml_file, append=!first_line)

        } else if (grepl("insert_trait", x, fixed = TRUE)) { #taxon insertion
            write(gsub("insert_trait", taxon_traits_string , x), file=xml_file, append = !first_line)

        }  else {
            write(x, file=xml_file, append=!first_line)
        }

        first_line = FALSE


    }
    close(con)


}
