#' Extract mutations
#'
#' Extract binary mutations from cells data frame 
#' @param mut_row integer
#' @importFrom stringr str_replace_all
#' 
#' @return vector
#' @export
extract_mutations <- function(mut_row) {
    split_1 <- strsplit(as.character(mut_row), split = ',')[[1]]
    split_2 <- stringr::str_replace_all(split_1, "[[:punct:]]", "")
    split_2 <- stringr::str_replace_all(split_2, " s", "")
    return(split_2)
}

#' Extract sequences
#'
#' Extract sequences from seq column 
#' @param seq string
#' @param sampled_muts vector, if NULL then sequence will not be filtere to observed variation
#' @param type string
#' @importFrom stringr str_replace_all
#' 
#' @return string
#' @export
extract_sequences <- function(seq, sampled_muts = NULL) {
    split_1 <- strsplit(as.character(seq), split = ',')[[1]]
    split_2 <- stringr::str_replace_all(split_1, "[[:punct:]]", "")
    if(is.null(sampled_muts)) {
        split_2 <- stringr::str_replace_all(split_2, " ", "")
    } else {
        split_2 <- stringr::str_replace_all(split_2, " ", "")[sampled_muts]
    }
    return(paste(split_2, collapse=""))
}

#' Convert JC sim to nexus
#'
#' Write simulated tumor to nexus file. 
#' @param sim_cells_file string
#' @param nexus_file string
#' @param sampled_cells_file string
#' @param n_sampled integer
#' @param new_sampled_cells_file string
#' @param type string
#' @param diversified_sampling  FALSE, if TRUE will select from categories based on weight
#' @param edge_weight weight for random diversified sampling
#' @param fixed_edge_num number of edge states to sample
#' @importFrom dplyr filter
#' @importFrom treeio as.treedata
#' @importFrom ape read.tree write.nexus
#' @importFrom stringr str_replace_all
#' @importFrom purrr map
#' 
#' @return numeric
#' @export
convert_sim_to_nex <- function(sim_cells_file,
                               nexus_file,
                               sampled_cells_file = NULL,
                               n_sampled = 100,
                               new_sampled_cells_file = NULL,
                               type = c("tree", "topology", "alignment", "states"),
                               diversified_sampling = FALSE,
                               edge_weight = 0.5,
                               fixed_edge_num = NULL) {
    
    
    ## RAW SIMULATION DATA ##
    
    sim_cells <- read.csv(file=sim_cells_file) %>% 
        normalize_locs %>% 
        dplyr::arrange(index)
    
    
    ## SAMPLING ##
    
    
    endpoint <- max(sim_cells$deathdate)
    
    alive_cells <- sim_cells %>% 
        dplyr::filter(deathdate == endpoint)
    
    if (is.null(sampled_cells_file)) { #if sample cells not provided
        
        message("Sampled cells not provided, sampling from alive cells...")
        sampled_cells <- sample_alive_cells(alive_cells = alive_cells,  n=n_sampled, diversified_sampling = diversified_sampling, edge_weight = edge_weight, fixed_edge_num = fixed_edge_num) %>% 
            dplyr::filter(sampled == TRUE)
        
        if (is.null(new_sampled_cells_file)) {
            warning("Location for newly sampled cells not provided and will not be saved")
            
        } else {
            write.csv(sampled_cells, new_sampled_cells_file)
        }
        
    } else if (file.exists(sampled_cells_file)) {
        
        sampled_cells <- read.csv(sampled_cells_file)
        
    } else {
        
        stop("Sampled cells file does not exist")
    }
    
    
    if (type == "tree") {
        message ("Generating tree nexus file...")
        

        ## FIND TIME TREE FROM SAMPLED CELLS ##

        
        
        #tree for fixed topology
        sim_nwk_list <- convert_nodes_to_string(sampled_cells$index,
                                                sim_cells,
                                                branch_unit = "time")
        
        tree <- treeio::as.treedata(ape::read.tree(text = sim_nwk_list$tree.text))
        
        phylo <- tree@phylo
        
        
        #add "cell" to avoid repeat between tip label and node numbering
        phylo$tip.label <- paste0("cell_", phylo$tip.label)
        phylo$node.label <- paste0("cell_", phylo$node.label)
        
        

        ## CONVERT TREE TO NEXUS FORMAT ##

        ape::write.nexus(phylo, file = nexus_file, translate = FALSE)
        
    } else if (type == "topology") {
        message ("Generating tree topology nexus file...")
        
        
        ## FIND TREE FROM SAMPLED CELLS ##
        

        #tree for fixed topology
        sim_nwk_list <- convert_nodes_to_string(sampled_cells$index,
                                                sim_cells,
                                                branch_unit = "none")
        

        tree <- as.treedata(read.tree(text = sim_nwk_list$tree.text))
        
        phylo <- tree@phylo
        
        
        #add "cell" to avoid repeat between tip label and node numbering
        phylo$tip.label <- paste0("cell_", phylo$tip.label)
        phylo$node.label <- paste0("cell_", phylo$node.label)
        
        
        ape::write.nexus(phylo, file = nexus_file, translate = FALSE)
        
    } else if (type == "states") {
        message ("Generating states nexus file...")
        
        ## FIND LEAF STATES ##
        
        #check surrounding number of cells for each sampled cell
        #this will be used for taxa states
        
        sampled_cells <- mark_boundary(sampled_cells, alive_cells)
        sampled_cell_states <- sampled_cells$est_edge
        
        ## CONVERT STATES TO NEXUS ##
        
        #translate into nexus format
        
        #initialize
        write(paste0("#NEXUS",
                     "\n",
                     "\nBEGIN DATA;",
                     "\n",
                     paste0("\nDIMENSIONS NTAX=",nrow(sampled_cells), " NCHAR=1;"),
                     '\nFORMAT datatype=Standard symbols="01" missing=? GAP=-;',
                     "\n",
                     "\nMATRIX",
                     "\n"),
              file = nexus_file,
              append=FALSE)
        
        
        #then append other lines
        for (i in 1:nrow(sampled_cells)) {
            write(paste0("cell_",sampled_cells$index[i],"\t", sampled_cell_states[i]),
                  file = nexus_file,
                  append=TRUE)
            
        }
        
        write(paste0(";","\nEND;"),
              file = nexus_file,
              append=TRUE)
        
    } else if ( type == "alignment") {
        
        message ("Generating alignment nexus file...")
        
        ## WRITE JC SEQUENCE DATA TO NEXUS ##
        
        #sampled_muts <- sampled_cells$mutations
        
        #all_sampled_muts <- sort(as.integer(unique(unlist(purrr::map(sampled_muts, function(mut_row) extract_mutations(mut_row))))))
        ## extract sequence data
        all_sampled_muts <- NULL
        
        sampled_cells$sequence_collapsed <- purrr::map_chr(sampled_cells$sequence, function(seq) extract_sequences(seq, sampled_muts = all_sampled_muts))
        
        
        ## WRITE JC SEQUENCE DATA TO NEXUS ##
        
        
        
        len_seq <- nchar(sampled_cells$sequence_collapsed[1])
        
        #initialize
        write(paste0("#NEXUS",
                     "\n",
                     "\nBEGIN TAXA;",
                     "\n",
                     paste0("\tDIMENSIONS NTAX=",nrow(sampled_cells), ";"),
                     '\n', 
                     '\t', 'TAXLABELS'),
              file = nexus_file,
              append=FALSE)
        
        
        #then append other lines
        
        for (i in 1:nrow(sampled_cells)) {
            write(paste0("\t", "cell_",sampled_cells$index[i]),
                  file = nexus_file,
                  append=TRUE)
            
        }
        
        write(paste0(";", "\n", "END;", "\n", "BEGIN CHARACTERS;",
                     "\n",
                     "\t",paste0("DIMENSIONS NCHAR =", len_seq,";"),
                     "\n", "\t", "FORMAT DATATYPE=DNA MISSING=? GAP=-;",
                     "\n", "\t", "MATRIX"),
              file = nexus_file,
              append=TRUE)
        
        for (i in 1:nrow(sampled_cells)) {
            write(paste0("\t", "cell_",sampled_cells$index[i],"\t", sampled_cells$sequence_collapsed[i]),
                  file = nexus_file,
                  append=TRUE)
            
        }
        
        write(paste0(";","\nEND;"),
              file = nexus_file,
              append=TRUE)
        
        
    } else {
        stop("Type not provided. Provide one of the following output types: tree, topology, states, alignment")
    }
    
}
