#' Get punched cells
#'
#' Collect punch in radius surrounding sampled cell
#' 
#' @param sampled_cell_index integer sampled cell that is center of punch
#' @param alive_cells data.frame all alive cells in simulated tumor
#' @param punch_radius_n_cells integer number of cells from index cells to punch, value0 will return single cels
#' @return data.frame
#' @export
get_punched_cells <- function(sampled_cell_index, alive_cells, punch_radius_n_cells) {

  locx <- alive_cells[alive_cells$index == sampled_cell_index, "locx" ]
  locy <- alive_cells[alive_cells$index == sampled_cell_index, "locy" ]

  #get coordinates of boundaries of punch biopsy
  max_x = locx + punch_radius_n_cells
  max_y = locy + punch_radius_n_cells
  min_x = locx - punch_radius_n_cells
  min_y = locy - punch_radius_n_cells

  #collect all alive cells within radius

  punched_cells <- alive_cells[(alive_cells$locx <= max_x) & 
                               (alive_cells$locx >= min_x) &
                               (alive_cells$locy <= max_y) &
                               (alive_cells$locy >= min_y), "index"]

  return(punched_cells)
}




#' Plot punch biopsies
#'
#' Visualize punch biopsies on simulated tumor
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom dplyr full_join
#' @importFrom purrr map map2
#' @importFrom scales muted
#' 
#' @param alive_cells data.frame all alive cells in simulated tumor
#' @param color_by_edge logical if TRUE, color by edge and center classification
#' @param color_by_punch_id logical if TRUE, color by punch id
#' @param size numeric size scale for scatter plot points
#' @return ggplot
#' 
#' @export
plot_punch_biopsies <- function(alive_cells, punched_cell_indices_list,
                                color_by_edge = FALSE,
                                color_by_punch_id = FALSE,
                                size = 2) {
  
  punched_cell_indices <- unlist(punched_cell_indices_list)
  alive_cells$punched <- alive_cells$index %in% punched_cell_indices
  
  if (color_by_edge & color_by_punch_id) {
    
    error("Can only color by edge or punch id, not both")
  }
  
  if (color_by_punch_id) {
    if (! "punch_id" %in% colnames(alive_cells) ) {
      punch_sizes <- purrr::map(punched_cell_indices_list, length)
      punch_id_vec <- unlist(purrr::map2(1:length(punched_cell_indices_list), punch_sizes, function(x,y) rep(x,y)))
      punches_df <- data.frame("index" = punched_cell_indices, "punch_id" = punch_id_vec)
      alive_cells <- dplyr::full_join(alive_cells, punches_df)
    }
    
    g <- ggplot(alive_cells, aes(x = locx, y = locy, color = as.factor(punch_id))) +
      geom_point(size = size) +
      theme_void() +
      theme(legend.position = "none")
    
  } else if (color_by_edge) {
    
    if (! "punch_on_edge" %in% colnames(alive_cells)) {
      warning("Must provide punch_on_edge column to color by edge")
      
      g <- ggplot(alive_cells, aes(x = locx, y = locy, color = as.factor(punched))) +
        geom_point(size = size) +
        theme_void() +
        scale_color_manual(values = c("darkgrey", muted('red'))) +
        theme(legend.position = "none")
    }
    
    g <- ggplot(alive_cells, aes(x = locx, y = locy, color = ifelse(punched & punch_on_edge, "edge_punch",
                                                                    ifelse(punched, "center_punch", "not_sampled")))) +
      geom_point(size = size) +
      theme_void() +
      scale_color_manual(values = c("not_sampled" = "darkgrey",
                                    "edge_punch" = scales::muted('red'),
                                    "center_punch" = "black")) +
      theme(legend.position = "none")
    
  } else {
    
    
    g <- ggplot(alive_cells, aes(x = locx, y = locy, color = as.factor(punched))) +
      geom_point(size = size) +
      theme_void() +
      scale_color_manual(values = c("darkgrey", muted('red'))) +
      theme(legend.position = "none")
  }
  return(g)
}


#' Shift cells
#'
#' Utility function for bulk sampling tumor, will shift punch over if there are not enough cells sampled. This is mostly useful for sampling near edge. 
#' 
#' @param proposed_sampled_cell_index integer, original punch index to shift
#' @param alive_cells data.frame of alive cells in simulated tumor
#' @return integer, index cell of new proposed center of punch
#' 
#' @export 
shift_punch <- function(proposed_sampled_cell_index, alive_cells) {
  
  
  sampled_locx <- alive_cells[alive_cells$index == proposed_sampled_cell_index,"locx"]
  sampled_locy <- alive_cells[alive_cells$index == proposed_sampled_cell_index,"locy"]
  
  
  check_x_vec <- c(sampled_locx + 1, sampled_locx, sampled_locx - 1, sampled_locx - 1, sampled_locx + 1, sampled_locx - 1, sampled_locx, sampled_locx + 1)
  check_y_vec <- c(sampled_locy + 1, sampled_locy + 1, sampled_locy + 1, sampled_locy, sampled_locy, sampled_locy - 1, sampled_locy - 1, sampled_locy - 1)
  
  adjacent_cells <- c()
  for (i in 1:length(check_x_vec)) {
    
    check_x <- check_x_vec[i]
    check_y <- check_y_vec[i]
    adjacent_cells <- c(adjacent_cells, any(alive_cells$locx == check_x & alive_cells$locy == check_y))
    
  }
  
  #choose direction to another free cell
  if (!any(adjacent_cells)) {
    return(NA)
  } else {
    direction <- which(adjacent_cells)[sample(1:sum(adjacent_cells), size = 1)]
    
    
    
    new_sampled_cell_locx <- check_x_vec[direction]
    new_sampled_cell_locy <- check_y_vec[direction]
    
    proposed_sampled_cell_index <- alive_cells[alive_cells$locx == new_sampled_cell_locx & alive_cells$locy == new_sampled_cell_locy, "index"]
    
    return(proposed_sampled_cell_index)
  }
}



#' Bulk sample tumor
#'
#' Select non-overlapping punch biopsies to sample from simulated tumor
#' 
#' @importFrom testit assert
#' @importFrom rlist list.append
#' @param alive_cells data.frame of alive cells in simulated tumor
#' @param n_samples integer number of punches to sample
#' @param punch_radius_n_cells integer radius of punches in number of cells
#' @param acceptable_overlap number of cells punches can overlap by
#' @param max_sampling_attempts number of sampling attempts before aborted
#' @param max_shifting_attempts number of punch shifts attempted before punch is re-drawn
#' @param max_punch_size_deviation accepted devation from expected punch size, due to punch overlap or edge
#' @return list of punched cells

bulk_sample_tumor <- function(alive_cells,
                              n_samples = 25,
                              punch_radius_n_cells = 2,
                              acceptable_overlap = 1,
                              max_sampling_attempts = 100, 
                              max_shifting_attempts = 10,
                              max_punch_size_deviation = 5) {
  
  message(paste0("Selecting random bulk punch biopsies... \n", "Tumor population size: ", nrow(alive_cells), " cells \n", "Desired number of punches: ",
                 n_samples, "\n", "Punch size: ", (1+2*punch_radius_n_cells)^2, " cells \n"))
  
  testit::assert("Expected number of cells sampled <= tumor population", n_samples * (1+2*punch_radius_n_cells^2) <= nrow(alive_cells))
  
  #keep track of all sampled cells so that cells will not be sampled x2
  curr_sampled_cells <- c()
  
  #keep list of cells within punch biopsies
  curr_punch_list <-list()
  
  #to ensure sampling from available cells
  curr_unsampled_cells <- alive_cells$index
  
  #for loop breaking if too many sampling attempts
  stop = FALSE
  
  for (i in 1:n_samples) {
    
    print(paste0("Punch ", i))
    
    #keep track of number of attempts to limit loop length
    #(controlled by max_number_of_sampling_attempts)
    
    number_of_sampling_attempts <- 1
    
    # initial punch attempt
    
    
    
    proposed_sampled_cell_index <- sample(curr_unsampled_cells, size = 1)
    
    #cell that leads to edge cut off
    
    
    #proposed_sampled_cell_index <- "71367"
    
    proposed_punched_cells <- get_punched_cells(proposed_sampled_cell_index,
                                                alive_cells = alive_cells,
                                                punch_radius_n_cells = punch_radius_n_cells)
    
    
    # sometimes if a punch is on the edge, less cells than expected will be sampled
    # we don't want to bias the edge versus sampling by excluding these, so instead will try to shift over to 
    # fit within the margins of the tumor
    
    number_of_shifting_attempts <- 1
    
    while (length(proposed_punched_cells) < ((2*punch_radius_n_cells+1)^2 - max_punch_size_deviation)) { #if punch smaller than expected
      
      #if shifting has been unsuccessful
      
      if (number_of_shifting_attempts >= max_shifting_attempts) {
        
        message("Unable to find good shift, resampling")
        proposed_sampled_cell_index <- sample(curr_unsampled_cells, size = 1)
        
        break
      }
      
      #
      message("Not enough cells in punch, attempting to shift..,")
      
      proposed_sampled_cell_index <- shift_punch(proposed_sampled_cell_index, alive_cells)
      
      # give up if no adjacent cells to go to
      if (is.na(proposed_sampled_cell_index )) {
        
        message("Unable to find good shift, resampling")
        
        proposed_sampled_cell_index <- sample(curr_unsampled_cells, size = 1)
        proposed_punched_cells <- get_punched_cells(proposed_sampled_cell_index,
                                                    alive_cells = alive_cells,
                                                    punch_radius_n_cells = punch_radius_n_cells)
        break
      }
      
      proposed_punched_cells <- get_punched_cells(proposed_sampled_cell_index,
                                                  alive_cells = alive_cells,
                                                  punch_radius_n_cells = punch_radius_n_cells)
      
      
      number_of_shifting_attempts <- number_of_shifting_attempts + 1
      
      
    }
    
    #check if overlapping with existing punch or if punch is too small
    while(sum(proposed_punched_cells %in% curr_sampled_cells) > acceptable_overlap | length(proposed_punched_cells) < ((2*punch_radius_n_cells+1)^2 - max_punch_size_deviation)) {
      
      message("Too much overlap with previous punch, resampling....")
      
      
      # propose punch biopsy
      
      proposed_sampled_cell_index <- sample(curr_unsampled_cells, size = 1)
      proposed_punched_cells <- get_punched_cells(proposed_sampled_cell_index,
                                                  alive_cells = alive_cells,
                                                  punch_radius_n_cells = punch_radius_n_cells)
      
      
      # if punch is on the edge try to shift over to the center
      number_of_shifting_attempts <- 1
      
      while (length(proposed_punched_cells) < ((2*punch_radius_n_cells+1)^2 - max_punch_size_deviation)) { #if punch smaller than expected
        
        #if shifting has been unsuccessful
        
        if (number_of_shifting_attempts >= max_shifting_attempts) {
          message("Unable to find good shift, resampling")
          proposed_sampled_cell_index <- sample(curr_unsampled_cells, size = 1)
          break
        }
        
        #
        message("Not enough cells in punch, attempting to shift..,")
        
        
        proposed_sampled_cell_index <- shift_punch(proposed_sampled_cell_index, alive_cells)
        message(paste0("Shifted to ", proposed_sampled_cell_index))
        
        # give up if no adjacent cells to go to
        if (is.na(proposed_sampled_cell_index)) {
          message("Unable to find good shift, resampling")
          proposed_sampled_cell_index <- sample(curr_unsampled_cells, size = 1)
          proposed_punched_cells <- get_punched_cells(proposed_sampled_cell_index,
                                                      alive_cells = alive_cells,
                                                      punch_radius_n_cells = punch_radius_n_cells)
          break
        }
        
        proposed_punched_cells <- get_punched_cells(proposed_sampled_cell_index,
                                                    alive_cells = alive_cells,
                                                    punch_radius_n_cells = punch_radius_n_cells)
        
        
        number_of_shifting_attempts <- number_of_shifting_attempts + 1
        
        
      }
      
      
      #break loop if attempted sampling too many times
      
      if (number_of_sampling_attempts > max_sampling_attempts) {
        
        warning("Too many sampling attempts: try lower number of desired samples, sample a larger tumor, or increase acceptable punch overlap.")
        stop = TRUE
        
        break
      
        } else {
        number_of_sampling_attempts  <- number_of_sampling_attempts + 1
      }
    }
    
    if (stop) {
      message("Stopping...")
      
      break
    }
    
    
    
    # final punched cells minus overlap with previously sampled punches (avoid double sampling)
    
    
    
    punched_cells <- proposed_punched_cells[which(!(proposed_punched_cells %in% curr_sampled_cells))]
    
    
    # keep track of sampled / unsampled cells and final punch list
    curr_punch_list <- rlist::list.append(curr_punch_list, punched_cells)
    curr_sampled_cells <- c(curr_sampled_cells, punched_cells)
    curr_unsampled_cells <- curr_unsampled_cells[-which(curr_unsampled_cells %in% punched_cells)]
    
  }
  
  
  message ("Done sampling!")
  return(curr_punch_list)
  
  
  
}

#' Extract punch allele freqs
#'
#' Find true allele frequencies of sampled frequencies in punch
#' 
#' @param id integer id of punch
#' @param mut_presence_absence_df data.frame of binary variant presense absence status for each variant
#' @param sampled_cells data.frame of sampled cells
#' 
#' @export


extract_punch_allele_freqs <- function(id, mut_presence_absence_df, sampled_cells) {
  # note that columns "punch_id" must exist for sampled cells
  
  if("punch_id" %in% colnames(sampled_cells)) {
    punch_muts <- mut_presence_absence_df[which(sampled_cells$punch_id == id),]
    
    if (is.null(nrow(punch_muts))) {
      punch_allele_freqs <- punch_muts
    } else {
      punch_allele_freqs <- matrix(colMeans(punch_muts), nrow = 1, byrow = TRUE)
    }
    return(punch_allele_freqs)
  } else {
    error("Column 'punch_id' must exist")
  }
  
}

#' Calc punch allele freqs
#'
#' Find true allele frequencies of sampled frequencies in punch
#' 
#' @param id integer id of punch
#' @param mut_presence_absence_df data.frame of binary variant presense absence status for each variant
#' @param sampled_cells data.frame of sampled cells
#' 
#' @export

calc_punch_allele_freqs <- function(punched_cell_indices_list, alive_cells) {
  
  #add punch information
  punch_sizes <- purrr::map(punched_cell_indices_list, length)
  punch_id_vec <- unlist(purrr::map2(1:length(punched_cell_indices_list), punch_sizes, function(x,y) rep(x,y)))
  punches_df <- data.frame("index" = unlist(punched_cell_indices_list), "punch_id" = punch_id_vec)
  alive_cells <- full_join(alive_cells, punches_df)
  
  #add punched column if not already
  if(! 'punched' %in% colnames(alive_cells)) {
    alive_cells$punched <- alive_cells$index %in% unlist(punched_cell_indices_list)
  }
  sampled_cells <- alive_cells %>% 
    dplyr::filter(punched == TRUE)
  
  #get full presence-absence matrix for all sampled cells (within all punches)
  mut_presence_absence_df <- define_mut_presence_absence(sampled_cells)
  
  punch_allele_freqs_list <- purrr::map(sort(unique(sampled_cells$punch_id)), function(id) extract_punch_allele_freqs(id = id,
                                                                                                                      mut_presence_absence_df = mut_presence_absence_df,
                                                                                                                      sampled_cells = sampled_cells))
  #names(punch_allele_freqs_list) = as.character(1:length(punch_allele_freqs_list))
  punch_allele_freqs <- do.call(rbind,punch_allele_freqs_list)
  
  
  
  return(punch_allele_freqs)
}


#based off parameters described in https://github.com/kchkhaidze/CHESS.cpp/blob/master/r_package/R/distribution_functions.R
draw_sequences_from_allele_freq <- function(true_allele_freq, min_alt_reads = 2, min_coverage = 10, depth = 100) {
  
  snv.depth <- rpois(1, depth) #coverage on specific SNV is poisson distributed
  alt_reads <- rbinom(1, snv.depth, true_allele_freq) #draw reads with probably of true SNV allele frequency
  
  
  #check if passes QC thresholds
  
  ref_reads <-  snv.depth - alt_reads #get ref reads
  
  if (snv.depth < min_coverage) { #if doesn't meet coverage threshold, call SNV as unknown
    return(NA)
  } else if (alt_reads < min_alt_reads) { # if only alternative reads are less than threshold then call as reference variant
    return(0)
  } else {
    return(alt_reads/snv.depth) #otherwise returned observed VAF
  }
  
}


generate_punch_presense_absense_matrix <- function(sequenced_punch_allele_freqs, vaf_threshold = 0) {
  
  
  presense_absense_matrix <- sequenced_punch_allele_freqs
  
  #enforce minimum VAF and call everything below as 0
  presense_absense_matrix[presense_absense_matrix <= vaf_threshold] = 0
  
  #otherwise 
  presense_absense_matrix[presense_absense_matrix  > vaf_threshold] = 1
  
  #filter to variants that have been called
  # i <- (colSums(presense_absense_matrix, na.rm=TRUE) != 0)
  
  #presense_absense_matrix <- presense_absense_matrix[,i]
  
  return(presense_absense_matrix)
}






# table of "mutations" with reference and alternative alleles

get_sampled_reference_seq <- function(all_cells, all_sampled_muts) {
  
  origin_cell_sequence <- all_cells[all_cells$index == 1,"sequence"]
  sampled_reference_sequence <- extract_sequences(origin_cell_sequence, sampled_muts = all_sampled_muts)
  
  return(sampled_reference_sequence)
}

extract_sampled_alt_allele <- function(sampled_cells, mut_id, all_sampled_muts = NA) {
  
  #get list with vector of mutation ids contained in each sequence
  sampled_muts_list <- purrr::map(sampled_cells$mutations, function(mut_row) extract_mutations(mut_row))
  
  rep_seq_contains_mut_id <- which(unlist(purrr::map(sampled_muts_list, function(muts_vec) mut_id %in% muts_vec)))[1]
  
  if (all(is.na(all_sampled_muts))) {
    all_sampled_muts <- sort(as.integer(unique(unlist(purrr::map(sampled_muts_list, function(mut_row) extract_mutations(mut_row))))))
  }
  
  if (! ("sequence_collapsed" %in% colnames(sampled_cells))) {
    
    
    sampled_cells$sequence_collapsed <- purrr::map_chr(sampled_cells$sequence, function(seq) extract_sequences(seq, sampled_muts = all_sampled_muts)) 
  }
  
  alt_allele <- substring(sampled_cells[rep_seq_contains_mut_id, "sequence_collapsed"],which(mut_id == all_sampled_muts), which(mut_id == all_sampled_muts))
  return(alt_allele)
  
}

make_mutations_table <- function(sampled_cells, all_cells) {
  
  sampled_muts <- sampled_cells$mutations
  
  all_sampled_muts <- sort(as.integer(unique(unlist(purrr::map(sampled_muts, function(mut_row) extract_mutations(mut_row))))))
  ## extract sequence data
  
  
  sampled_cells$sequence_collapsed <- purrr::map_chr(sampled_cells$sequence, function(seq) extract_sequences(seq, sampled_muts = all_sampled_muts)) 
  
  sampled_reference_sequence <- strsplit(get_sampled_reference_seq(all_cells, all_sampled_muts = all_sampled_muts), "")
  
  mut_alleles_df <- data.frame("mut_id" = all_sampled_muts, "ref" = sampled_reference_sequence[[1]])
  
  mut_alleles_df$alt <- purrr::map_chr(mut_alleles_df$mut_id, function(mut_id) extract_sampled_alt_allele(sampled_cells = sampled_cells,
                                                                                                          mut_id = mut_id, 
                                                                                                          all_sampled_muts = all_sampled_muts))
  
  return(mut_alleles_df)
  
}

choose_allele <- function(mut_num, punch_num, sequenced_punch_presense_absense_matrix, mut_alleles_df) {
  
  mut_call <- sequenced_punch_presense_absense_matrix[punch_num, mut_num]
  
  if (is.na(mut_call)) {
    
    allele <- NA
    
  } else if (mut_call == 1) {
    
    allele = mut_alleles_df$alt[mut_num]
    
  } else if (mut_call == 0) {
    allele = mut_alleles_df$ref[mut_num]
  } else {
    error("Must input valid presense absense matrix")
  }
  
  return(allele)
}

#function to get alignment sequence for each punch

extract_punch_sequences <- function(sequenced_punch_presense_absense_matrix, mut_alleles_df) {
  
  punch_sequence_vec <- c()
  
  for (punch in 1:nrow(sequenced_punch_presense_absense_matrix)) {
    punch_seq <- purrr::map_chr(1:ncol(sequenced_punch_presense_absense_matrix), function(mut_num) choose_allele(mut_num = mut_num,
                                                                                                                 punch_num = punch,
                                                                                                                 sequenced_punch_presense_absense_matrix = sequenced_punch_presense_absense_matrix,
                                                                                                                 mut_alleles_df = mut_alleles_df) )
    
    punch_seq <- paste(punch_seq, sep="", collapse = "")
    punch_sequence_vec <- c(punch_sequence_vec, punch_seq)
  }
  
  return(punch_sequence_vec)
  
}

#get actual simulated sequences

#wrapper function to "sequence" sampled biopsies

sequence_punch_biopsies <- function(punch_indices, alive_cells,  mut_alleles_df, min_alt_reads = 2, min_coverage = 10, depth = 100, vaf_threshold = 0) {
  
  true_punch_allele_freqs <- calc_punch_allele_freqs(punched_cell_indices_list = punch_indices, alive_cells = alive_cells)
  
  sequenced_punch_allele_freqs <- apply(true_punch_allele_freqs, 1:2, draw_sequences_from_allele_freq, min_alt_reads = min_alt_reads, min_coverage = min_coverage, depth = depth)
  
  sequenced_punch_presense_absense_matrix <- generate_punch_presense_absense_matrix(sequenced_punch_allele_freqs, vaf_threshold = vaf_threshold)
  
  punch_sequences_vec <- extract_punch_sequences(sequenced_punch_presense_absense_matrix = sequenced_punch_presense_absense_matrix,
                                                 mut_alleles_df = mut_alleles_df)
  
  return(punch_sequences_vec)
}

#wrapper function to sample and sequence

bulk_sample_and_sequence <- function(all_cells, punch_radius_n_cells = 2, n_samples = 20, min_alt_reads = 2, min_coverage = 10, depth = 100, vaf_threshold = 0, seed = NA) {
  
  if (! is.na(seed)) {
    set.seed(seed)
  }
  
  alive_cells <- all_cells %>% 
    dplyr::filter(deathdate == max(deathdate))
  
  message("Sampling...")
  punch_indices <- bulk_sample_tumor(alive_cells, punch_radius_n_cells = punch_radius_n_cells, n_samples = n_samples)
  
  sampled_cells <- alive_cells %>% 
    dplyr::filter(index %in% unlist(punch_indices))
  
  mut_alleles_df <- make_mutations_table(sampled_cells = sampled_cells, all_cells = all_cells)
  
  message("Sequencing...")
  punch_sequences_vec <- sequence_punch_biopsies(punch_indices,
                                                 alive_cells = alive_cells,
                                                 mut_alleles_df = mut_alleles_df,
                                                 min_alt_reads = min_alt_reads,
                                                 min_coverage = min_coverage,
                                                 depth = depth,
                                                 vaf_threshold = vaf_threshold)
  
  return(list("indices" = punch_indices, "sequences" = punch_sequences_vec))
  
}


calc_cell_distance_from_edge <- function(cell_index, alive_cells) {
  
  if (! "est_edge" %in% colnames(alive_cells)) {
    
    warning("Edge cells were not marked, marking now...")
    alive_cells <-  mark_boundary(cells_to_mark = alive_cells, alive_cells = alive_cells)
  }
  
  #get coordinates of boundary cells
  edge_x_locs <- alive_cells[alive_cells$est_edge == 1, "locx"]
  edge_y_locs <- alive_cells[alive_cells$est_edge == 1, "locy"]
  
  cell_x_loc <- alive_cells[alive_cells$index == cell_index, "locx"]
  cell_y_loc <- alive_cells[alive_cells$index == cell_index, "locy"]
  
  dist_vec <- sqrt((edge_x_locs -  cell_x_loc)^2 + (edge_y_locs -  cell_y_loc)^2)
  
  #return minimum distance from boundary cell
  return(min(dist_vec))
  
}
# mark edge vs center
mark_edge_punches <- function(punched_cell_indices_list, alive_cells, fraction_edge_cutoff = 0) {
  
  if (! "punch_id" %in% colnames(alive_cells)) {
    punch_sizes <- purrr::map(punched_cell_indices_list, length)
    punch_id_vec <- unlist(purrr::map2(1:length(punched_cell_indices_list), punch_sizes, function(x,y) rep(x,y)))
    punches_df <- data.frame("index" = unlist(punched_cell_indices_list), "punch_id" = punch_id_vec)
    alive_cells <- full_join(alive_cells, punches_df)
  }
  
  #first mark boundary single cells 
  if (! "est_edge" %in% colnames(alive_cells)) {
    
    alive_cells <-  mark_boundary(cells_to_mark = alive_cells, alive_cells = alive_cells)
  }
  
  sampled_cells <- alive_cells[alive_cells$index %in% unlist(punched_cell_indices_list),]
  
  sampled_cells$distance_from_edge <-  purrr::map_dbl(sampled_cells$index, function(cell_index) calc_cell_distance_from_edge(cell_index = cell_index,
                                                                                                                             alive_cells = alive_cells))
  
  #calculate total span of tumor 
  #estimate by averages x and y span
  
  est_tumor_span <- mean((max(alive_cells$locx) - min(alive_cells$locx)),(max(alive_cells$locy) - min(alive_cells$locy)) )
  
  #define distance from boundary to be defined as edge 
  dist_cutoff <- fraction_edge_cutoff * est_tumor_span
  
  
  punches_edge_assigment_df <- sampled_cells %>% 
    group_by(punch_id) %>% 
    summarize("punch_on_edge" = any(distance_from_edge <= dist_cutoff))
  
  
  return(punches_edge_assigment_df)
  
}


add_missing_punch_info_to_cells <- function(alive_cells, punch_info = NULL, punches_edge_assigment_df = NULL) {
  
  
  
  if (! (("punch_id" %in% colnames(alive_cells)) & ( "punched" %in% colnames(alive_cells)))) {
    if (is.null(punch_info)) {
      
      warning("Missing punch info: punch info not added")
      
    } else {
      punched_cell_indices_list <- punch_info$indices
      punch_sizes <- purrr::map(punched_cell_indices_list, length)
      punch_id_vec <- unlist(purrr::map2(1:length(punched_cell_indices_list), punch_sizes, function(x,y) rep(x,y)))
      punches_df <- data.frame("index" = unlist(punched_cell_indices_list), "punch_id" = punch_id_vec)
      alive_cells <- full_join(alive_cells, punches_df)
      message("Punch info added")
    }
  }
  
  if (! "punch_on_edge" %in% colnames(alive_cells)) {
    alive_cells <- full_join(alive_cells, punches_edge_assigment_df)
    message("Punch edge status added")
    
  }
  
  
  return(alive_cells)
}


generate_bulk_seq_nexus_files <- function(punch_info, punches_edge_assigment_df, nexus_file_basename ) {
  
  n_punches <- length(punch_info$sequences)
  #states nexus file
  
  nexus_file <- paste0(nexus_file_basename, ".nex", sep = "")
  nexus_state_file <- paste0(nexus_file_basename, "_edge_state.nex", sep = "")
  
  message ("Generating states nexus file...")
  
  
  ## CONVERT STATES TO NEXUS ##
  
  #translate into nexus format
  
  #initialize
  write(paste0("#NEXUS",
               "\n",
               "\nBEGIN DATA;",
               "\n",
               paste0("\nDIMENSIONS NTAX=",n_punches, " NCHAR=1;"),
               '\nFORMAT datatype=Standard symbols="01" missing=? GAP=-;',
               "\n",
               "\nMATRIX",
               "\n"),
        file = nexus_state_file,
        append=FALSE)
  
  
  #then append other lines
  for (i in 1:n_punches) {
    write(paste0("punch_",i,"\t", as.integer(punches_edge_assigment_df$punch_on_edge[i])),
          file = nexus_state_file,
          append=TRUE)
    
  }
  
  write(paste0(";","\nEND;"),
        file = nexus_state_file,
        append=TRUE)
  
  message ("Generating alignment nexus file...")
  
  
  ## WRITE JC SEQUENCE DATA TO NEXUS ##
  
  
  len_seq <- nchar(punch_info$sequences[1])
  
  #initialize
  write(paste0("#NEXUS",
               "\n",
               "\nBEGIN TAXA;",
               "\n",
               paste0("\tDIMENSIONS NTAX=",n_punches, ";"),
               '\n', 
               '\t', 'TAXLABELS'),
        file = nexus_file,
        append=FALSE)
  
  
  #then append other lines
  
  for (i in 1:n_punches) {
    write(paste0("\t", "punch_",i),
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
  
  for (i in 1:n_punches) {
    write(paste0("\t", "punch_",i,"\t", punch_info$sequences[i]),
          file = nexus_file,
          append=TRUE)
    
  }
  
  write(paste0(";","\nEND;"),
        file = nexus_file,
        append=TRUE)
  
}

#calculate expected tree length based on pop gen theory
#this can be used to normalize a prior

calc_expected_total_tree_length <- function(n_tips, root_age) {
  
  sum_i <- 0
  for (i in 1:(n_tips - 1)) {
    sum_i <- sum_i + 1
  }
  
  E_total_tree_length <- root_age * sum_i / (1- 1/n_tips)
  return(E_total_tree_length)
}
write_rb_header_file <- function(rb_header_file, nexus_file_basename, simulation_name, alive_cells, n_tips, pop) {
  
  root_age <- max(alive_cells$deathdate)
  E_tree_length <- calc_expected_total_tree_length(n_tips, root_age)
  rate_pr <- E_tree_length / 10
  
  # write model info
  write('#########################################################
## RevBayes script to combine time tree                
## estimation and BiSSE model for differences in       
## diversification rates of edge and center            
## on simulated tumor data 
## Based on: 
## https://revbayes.github.io/tutorials/clocks/        
## https://revbayes.github.io/tutorials/sse/bisse.html
## ', rb_header_file)
  write(paste0("## Script generated: ", Sys.time()), rb_header_file, append = TRUE)
  write('#########################################################
                  ', rb_header_file, append = TRUE)
  
  # write simulation-specific parameters
  write("## input parameters", rb_header_file, append = TRUE)
  write(paste0('SIMULATION= "', simulation_name, '"', sep = ""), rb_header_file, append = TRUE)
  write(paste0('ROOT_AGE= ', round(root_age,2), sep = ""), rb_header_file, append = TRUE)
  write(paste0('N= ', pop, sep = ""), rb_header_file, append = TRUE)
  #write(paste0('rate_pr = ', rate_pr, sep = ""), rb_header_file, append = TRUE)
  
  
}


bisse_bulk_sample_and_sequence <- function(all_cells, simulation_name, rb_header_file, nexus_file_basename, punch_radius_n_cells = 2, n_samples = 20, min_alt_reads = 2, min_coverage = 10, depth = 100, vaf_threshold = 0, seed = NA) {
  
  # sampling, sequence
  alive_cells <- all_cells %>%
    dplyr::filter(deathdate == max(deathdate))
  
  punch_info <- bulk_sample_and_sequence(all_cells, punch_radius_n_cells = punch_radius_n_cells, n_samples = n_samples, min_alt_reads = min_alt_reads, min_coverage = min_coverage, depth = depth, vaf_threshold = vaf_threshold, seed = seed)
  punches_edge_assigment_df <- mark_edge_punches(punch_info$indices, alive_cells = alive_cells)
  
  alive_cells <- add_missing_punch_info_to_cells(alive_cells, punch_info, punches_edge_assigment_df )
  
  a <- plot_punch_biopsies(alive_cells, punch_info$indices, color_by_edge = TRUE, size = 2)
  
  print(a)
  
  b <- plot_punch_biopsies(alive_cells, punch_info$indices, color_by_punch_id = TRUE, size = 2)
  
  print(b)
  
  generate_bulk_seq_nexus_files(punch_info = punch_info,
                                punches_edge_assigment_df = punches_edge_assigment_df,
                                nexus_file_basename = nexus_file_basename)
  
  n_tips <- length(punch_info$sequences)
  pop <- nrow(alive_cells)
  write_rb_header_file(rb_header_file = rb_header_file, nexus_file_basename = nexus_file_basename, simulation_name = simulation_name, alive_cells = alive_cells, n_tips = n_tips, pop = pop)
}
