
#' Get color palette
#'
#' Get color palette to use to plot edge versus center, boundary-driven versus unrestricted growth.
#'
#' @param names vector of names to include for colors. For example "edge", "center"
#'
#' @return vector
#' @export
get_color_palette <- function(names = c("edge", "center", "boundary_driven", "unrestricted")){

  all_colors <- c("edge" = "#A2D2E2", "center" = "#89352F", "boundary_driven" = "#849324", "unrestricted" = "#FFB30F")
  return(all_colors[names])

}

#' Normalize locations
#'
#' Normalize x-y coordinates in of simulated cells.
#'
#' @param cells data.frame
#'
#' @return data.frame
#' @export
normalize_locs <- function(cells){
  center_x <- min(cells$locx) + (max(cells$locx) - min(cells$locx))/2
  center_y <- min(cells$locy) + (max(cells$locy) - min(cells$locy))/2
  norm_cells <- cells %>%
    mutate("norm_locx" = locx - center_x, "norm_locy" = locy - center_y)
  return(norm_cells)
}

#' Filter alive cells
#'
#' Filter cells to those alive at the end of the simulation.
#'
#' @param cells_df data.frame
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#'
#' @return data.frame
#' @export
filter_alive_cells <- function(cells_df) {
  endpoint <- max(cells_df$deathdate)

  alive_cells <- cells_df %>%
    dplyr::filter(deathdate == endpoint)

  return(alive_cells)
}

#' Diversified sampling
#'
#' Samples alive cells from simulated tumor with maximized distance between points
#'
#' @param alive_cells data.frame Alive cells in tumor
#' @importFrom rdist farthest_point_sampling
#' @importFrom magrittr %>%
#'
#' @return data.frame
#' @export
diversified_sampling <- function(alive_cells, n_sampled_cells) {

  #randomize alive cells order to prevent bias in sampling

  alive_cells <- alive_cells[sample(1:nrow(alive_cells), nrow(alive_cells), replace = FALSE),]

  dist_mat <- rdist::pdist(as.matrix(data.frame("x"=alive_cells$locx, "y"=alive_cells$locy), ncol = 2))
  fps <- rdist::farthest_point_sampling(dist_mat)

  #mark sampled cells
  alive_cells$sampled <- FALSE

  alive_cells$sampled[fps[1:n_sampled_cells]] <- TRUE

  #re-sort alive_cells by index
  alive_cells <- alive_cells %>%
    arrange(index)

  return(alive_cells)


}

#' Sample alive cells
#'
#' Randomly sample alive cells.
#'
#' @param alive_cells data.frame
#' @param n integer
#' @param diversified_sampling logical for diversified sampling
#' @param edge_weight probability of edge sample for diversified sampling -- parameter for binomial distribution
#' @param fixed_edge_num number of edge cells to sample, if not provided number of edge samples are randomly drawn
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#'
#' @return data.frame
#' @export
sample_alive_cells <- function(alive_cells, n, diversified_sampling= FALSE, weighted_sampling = FALSE, edge_weight = 0.5,  fixed_edge_num = NULL){

  if (diversified_sampling & (! weighted_sampling ) & is.null(fixed_edge_num)) {

      sampled_alive_cells <- diversified_sampling(alive_cells, n_sampled_cells = n)

  } else {

    if (weighted_sampling | (!(is.null(fixed_edge_num)))) { # if sampling is  weighted towards the edge

      #mark boundary cells
      alive_cells <- mark_boundary(alive_cells, alive_cells)

      #separate edge versus center cells
      edge_cells <- alive_cells[alive_cells$est_edge == 1,]
      center_cells <- alive_cells[alive_cells$est_edge == 0,]


      if(!(is.null(fixed_edge_num))) { #if fixed_edge_num is povided, then pool choices are not random

        message(paste0("Fixing number of edge samples to ", fixed_edge_num))

        pool_choices <- c(rep(1, fixed_edge_num), rep(0, n - fixed_edge_num))

      } else {

        #for each draw first choose from edge or center cell pools to sample from
        pool_choices <- rbinom(n = min(n, nrow(edge_cells)), size = 1, prob = edge_weight)
      }

      #draw edge and center pools
      edge_choices <- sample(edge_cells$index, size = min(n, nrow(edge_cells)))
      center_choices <- sample(center_cells$index, size = min(n, nrow(edge_cells)))

      sampled_indexes <- pool_choices * edge_choices + (1- pool_choices) * center_choices

    } else if ((! diversified_sampling) & (! weighted_sampling ) & is.null(fixed_edge_num)) {

      sampled_indexes <- sample(alive_cells$index,
                                          size = min(n, nrow(alive_cells)),
                                          replace = FALSE)



    } else {

      error("Conflicting sampling instructions,  choose one")
    }

  sampled_alive_cells <- alive_cells %>%
    dplyr::mutate("sampled" = index %in% sampled_indexes)
  }

  return(sampled_alive_cells)
}


#' Process and sample cells
#'
#' Wrapper function to filter and sample alive cells at the end of the simulation.
#'
#' @param alive_cells data.frame
#' @param n integer
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#'
#' @return data.frame
#' @export
process_and_sample_cells <- function(cells,n) {
  cells <- cells %>%
    filter_alive_cells %>%
    sample_alive_cells(.,n = n) %>%
    dplyr::filter(sampled == TRUE)

  return(cells)
}


#' Calculate pairwise genetic distance
#'
#' Calculate molecular branch lengths when generating tree given child and ancestor indices.
#'
#' @param child integer
#' @param ancestor integer
#' @param all_muts data.frame
#'
#' @return numeric
#' @export
calc_pairwise_genetic_distance <- function(child, ancestor, all_muts) {
  child_muts <- all_muts[child,] #get mutations in child node
  ancestor_muts <- all_muts[ancestor,] #get mutations in ancestor node
  mean_pairwise_differences <- sum(child_muts != ancestor_muts) #get pairwise dist/total muts
  return(mean_pairwise_differences)
}

#' Isolate sampled mutations
#'
#' Generate presence-absence matrix based on sampled mutations.
#'
#' @param mutations_column integer
#'
#' @return vector
#' @export
isolate_muts <- function(mutations_column) {

  #for each row split and filter to unique mutations in that cell
  i = 1
  mutations_vec <- c()
  while (i <= length(mutations_column)){
    row_muts <- na.omit(
      as.numeric(
        unlist(
          strsplit(as.character(mutations_column[i]),"[^0-9]")
        )
      )
    )

    mutations_vec <- c(mutations_vec, row_muts)
    i <- i + 1
  }

  return(sort(unique(mutations_vec)))
}

#' Get row mutations
#'
#' Extract vector of 0-1 mutations from string.
#'
#' @param mutations_row integer
#'
#' @return vector
#' @export
get_row_muts <- function(mutations_row) {

  row_muts <- na.omit(
    as.numeric(
      unlist(
        strsplit(as.character(mutations_row),"[^0-9]")
      )
    )
  )
  return(row_muts)

}

#' Get row mutations by i
#'
#' Alternate function to extract mutations by column.
#'
#' @param i integer
#' @param mutations_column integer
#'
#' @return vector
#' @export
get_row_muts_by_i <- function(i, mutations_column) {
  row_muts <- na.omit(
    as.numeric(
      unlist(
        strsplit(as.character(mutations_column[i]),"[^0-9]")
      )
    )
  )
  return(row_muts)
}

#' Compare mutations
#'
#' Compare all mutations with row.
#'
#' @param row_muts vector
#' @param all_muts data.frame
#'
#' @return numeric
#' @export
compare_muts <- function(row_muts, all_muts) {
  return(as.numeric(all_muts %in% row_muts))
}

#' Define mutations presence-absence
#'
#' Create presence-absence data.frame of samples vs sampled mutations.
#'
#' @param sampled_cells data.frame
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom future plan multisession
#'
#' @return data.frame
#' @export
define_mut_presence_absence <- function(sampled_cells) {
  #find all mutations in sample
  all_muts <- isolate_muts(mutations_column = sampled_cells$mutations)


  rowise_muts <- purrr::map(sampled_cells$mutations, get_row_muts)

  mut_presence_absence <- purrr::map(rowise_muts, ~compare_muts(.,all_muts)) %>%
    do.call(rbind,.)

  return(mut_presence_absence)
}

#' Collect all ancestors
#'
#' get all ancestors traversing up tree of indexed node.
#'
#' @param index integer
#' @param all_cells data.frame
#'
#' @return vector
#' @export
collect_all_ancestors <- function(index, all_cells) {
  ancestors = c()
  while (index != 1) {
    new_index <- all_cells[all_cells$index == index,"parent_index"][[1]]
    ancestors <- c(ancestors, new_index)
    index = new_index
  }
  return(ancestors)
}

#' Find MRCA
#'
#' Find MRCA between two cells given indexes.
#'
#' @param index_1 integer
#' @param index_2 integer
#' @param all_cells data.frame
#'
#' @return integer
#' @export
find_MRCA <- function(index_1, index_2, all_cells) {
  ancestors_1 <- collect_all_ancestors(index_1, all_cells)
  ancestors_2 <- collect_all_ancestors(index_2, all_cells)
  mrca <- ancestors_1[min(which(ancestors_1 %in% ancestors_2))]
  return(mrca)
}

#' Convert nodes to string
#'
#' Convert sampled simulated cells into newick formatted string and node list.
#'
#' @param orphan_cells vector of sampled leaves
#' @param all_cells data.frame of all cells in simulation
#' @param branch_unit Unit to calcualte branch length. Options: "time", "molecular", "generations", "none".
#' @importFrom dplyr arrange select
#' @importFrom magrittr %>%
#' @importFrom purrr map2_dbl
#' @importFrom tcltk tkProgressBar
#' @importFrom future plan multisession
#'
#' @return list
#' @export
convert_nodes_to_string <- function(orphan_cells,
                                    all_cells,
                                    branch_unit = c("time", "molecular", "generations", "none")){

  #should be already sorted by index, but just to make sure
  all_cells <- all_cells %>% dplyr::arrange(index)

  #make presence-absence df
  all_muts <- define_mut_presence_absence(all_cells) %>%
    as.data.frame %>%
    tibble::add_column("index" = all_cells$index) %>%
    dplyr::arrange(index) %>%
    dplyr::select(-index)

  #start out with each leaf as its own subtree
  subtrees <- as.character(orphan_cells)

  #keep track of all nodes in tree
  node_vec <- orphan_cells
  n <- length(orphan_cells)

  starting_orphans <- length(orphan_cells)
  # create progress bar
  #pb <- tkProgressBar(title = "progress bar", min = 0,
   #                   max = starting_orphans, width = 300)

  pb <- txtProgressBar(min = 0, max = starting_orphans, style = 3, label = "Making tree")
  while(length(subtrees) > 1) { #will connect subtrees by MCRA until reach single tree

    #make all pariwise combinations of cells without mrca
    all_combs <- as.data.frame(t(combn(orphan_cells, 2)))


    #find ancestors of candidate cells
    common_ancestors <- purrr::map2_dbl(all_combs$V1, all_combs$V2, function(index_1, index_2) find_MRCA(index_1 = index_1,
                                                                                                         index_2 = index_2,
                                                                                                         all_cells = all_cells))
    common_ancestors_data <- all_cells[common_ancestors,]

    #find youngest ancestor and two children
    youngest_ancestor <- common_ancestors_data$index[common_ancestors_data$birthdate == max(common_ancestors_data$birthdate)]

    #in case more than 2 cells have same ancestor,
    #choose one to start, we will come back to the other
    youngest_ancestor <- youngest_ancestor[1]




    #keep track of children descended from that ancestor

    children <- unique(unlist(all_combs[which(common_ancestors == youngest_ancestor),]))


    branch_lengths <- c()
    for (child in children) {


      c <- tail(which(node_vec == child), n = 1)

      if (branch_unit == "molecular") {

        bl <- calc_pairwise_genetic_distance(child, youngest_ancestor, all_muts = all_muts)

        branch_lengths <- c(branch_lengths, as.character(round(bl, 5)))

      } else if(branch_unit == "generations") {

        ancestors <- collect_all_ancestors(index = child, all_cells)
        gens <- which(ancestors == youngest_ancestor)
        branch_lengths <- c(branch_lengths, as.character(gens))

      } else {

        t_a <- all_cells[youngest_ancestor, "deathdate"]
        t_c <- all_cells[child, "deathdate"] #use deathdate

        branch_lengths <- c(branch_lengths, as.character(round(t_c-t_a,4)))

      }


    }


    #make new subtree from combining children and MCRA
    if (branch_unit == "none") { #don't include branch lengths
      subtrees <- c(subtrees, paste0("(",paste(subtrees[which(orphan_cells %in% children)], collapse = ","),")",as.character(youngest_ancestor)))

    } else {

      subtrees <- c(subtrees, paste0("(",paste(paste(subtrees[which(orphan_cells %in% children)],branch_lengths,sep=":"), collapse = ","),")",as.character(youngest_ancestor)))
    }



    subtrees <- subtrees[-which(orphan_cells %in% children)] #remove child trees after used

    orphan_cells <- c(orphan_cells, youngest_ancestor) #add ancestor as orphan
    node_vec <- c(node_vec, youngest_ancestor)
    orphan_cells <- orphan_cells[-which(orphan_cells %in% children)] #remove children

    prog <- starting_orphans - length(orphan_cells)
    #setTkProgressBar(pb, prog, label=paste("Tree building ", round(prog/starting_orphans*100, 0),
    #                                     "% done"))

    setTxtProgressBar(pb, prog)

  }

  #setTkProgressBar(pb, pb, label="Tree building 100% done")

  close(pb)
  return(list("tree.text"= paste0(tail(subtrees, n = 1),";"), "node_vec" = node_vec))
}

#' Convert nodes to string fast
#'
#' Convert sampled simulated cells into newick formatted string and node list but faster.
#'
#' @param orphan_cells vector of sampled leaves
#' @param all_cells data.frame of all cells in simulation
#' @param branch_unit Unit to calcualte branch length. Options: "time", "molecular", "generations", "none".
#' @importFrom dplyr arrange select
#' @importFrom magrittr %>%
#' @importFrom purrr map2_dbl
#'
#' @return list
#' @export
convert_nodes_to_string_fast <- function(orphan_cells,
                                    all_cells,
                                    branch_unit = c("time", "molecular", "generations", "none")){

  #should be already sorted by index, but just to make sure
  all_cells <- all_cells %>% dplyr::arrange(index)

  #make presence-absence df
  all_muts <- define_mut_presence_absence(all_cells) %>%
    as.data.frame %>%
    tibble::add_column("index" = all_cells$index) %>%
    dplyr::arrange(index) %>%
    dplyr::select(-index)

  #start out with each leaf as its own subtree
  subtrees <- as.character(orphan_cells)

  #keep track of all nodes in tree
  node_vec <- orphan_cells
  n <- length(orphan_cells)

  #make all pariwise combinations of cells without mrca
  all_combs <- as.data.frame(t(combn(orphan_cells, 2)))

  #find ancestors of candidate cells
  common_ancestors <- purrr::map2_dbl(all_combs$V1, all_combs$V2, function(index_1, index_2) find_MRCA(index_1 = index_1,
                                                                                                       index_2 = index_2,
                                                                                                       all_cells = all_cells))


  while(length(subtrees) > 1) { #will connect subtrees by MCRA until reach single tree

    #make all pariwise combinations of cells without mrca
    #all_combs <- as.data.frame(t(combn(orphan_cells, 2)))


    #find ancestors of candidate cells
    #common_ancestors <- purrr::map2_dbl(all_combs$V1, all_combs$V2, function(index_1, index_2) find_MRCA(index_1 = index_1,
    #                                                                                                     index_2 = index_2,
    #                                                                                                     all_cells = all_cells))
    common_ancestors_data <- all_cells[common_ancestors,]

    #find youngest ancestor and two children
    youngest_ancestor <- common_ancestors_data$index[common_ancestors_data$birthdate == max(common_ancestors_data$birthdate)]

    #in case more than 2 cells have same ancestor,
    #choose one to start, we will come back to the other
    youngest_ancestor <- youngest_ancestor[1]




    #keep track of children descended from that ancestor

    children <- unique(unlist(all_combs[which(common_ancestors == youngest_ancestor),]))


    branch_lengths <- c()
    for (child in children) {

      c <- tail(which(node_vec == child), n = 1)

      if (branch_unit == "molecular") {

        bl <- calc_pairwise_genetic_distance(child, youngest_ancestor, all_muts = all_muts)

        branch_lengths <- c(branch_lengths, as.character(round(bl, 5)))

      } else if(branch_unit == "generations") {
        ancestors <- collect_all_ancestors(index = child, all_cells)
        gens <- which(ancestors == youngest_ancestor)
        branch_lengths <- c(branch_lengths, as.character(gens))

      }else {
        t_a <- all_cells[youngest_ancestor, "deathdate"]
        t_c <- all_cells[child, "deathdate"] #use deathdate

        branch_lengths <- c(branch_lengths, as.character(round(t_c-t_a,4)))

      }


    }


    #make new subtree from combining children and MCRA
    if (branch_unit == "none") { #don't include branch lengths
      subtrees <- c(subtrees, paste0("(",paste(subtrees[which(orphan_cells %in% children)], collapse = ","),")",as.character(youngest_ancestor)))

    } else {
      subtrees <- c(subtrees, paste0("(",paste(paste(subtrees[which(orphan_cells %in% children)],branch_lengths,sep=":"), collapse = ","),")",as.character(youngest_ancestor)))
    }
    subtrees <- subtrees[-which(orphan_cells %in% children)] #remove child trees after used

    orphan_cells <- c(orphan_cells, youngest_ancestor) #add ancestor as orphan
    node_vec <- c(node_vec, youngest_ancestor)

    #adjust all data to only include orphans
    orphans_to_remove <- orphan_cells[orphan_cells %in% children]

    all_combs <- all_combs[-which(all_combs$V1 %in% orphans_to_remove | all_combs$V2 %in% orphans_to_remove)]

    orphan_cells <- orphan_cells[-which(orphan_cells %in% children)] #remove children

  }
  return(list("tree.text"= paste0(tail(subtrees, n = 1),";"), "node_vec" = node_vec))
}

#' Convert newick to treedata
#'
#' Convert tree in newick formatted string to ggtree object with location data.
#' For easy plotting and use with R phylogenetic tools
#'
#' @param nwk_list list of newick formated tree string and vector of indices in tree. Output of convert_nodes_to_string.
#' @param all_cells data.frame of all cells in the simulation
#' @param sampled_cells data.frame of sampled alive cells
#' @param include_edge_states logical to specify if edge states will be included in treedata. Default FALSE
#' @importFrom dplyr arrange select
#' @importFrom ggtree read.tree
#' @importFrom magrittr %>%
#' @importFrom purrr map_dbl map
#' @importFrom tibble as_tibble add_column
#' @importFrom tidyr unite
#' @importFrom treeio as.treedata
#'
#' @return treedata
#' @export
convert_nwk_to_treedata <- function(nwk_list,
                                    all_cells,
                                    sampled_cells,
                                    include_edge_states = FALSE) {

  # normalize x y coordinate if not done

  if (!("norm_locx" %in% colnames(sampled_cells))) {

    sampled_cells <- sampled_cells %>%
      normalize_locs
  }

  # read newick string into tree and convert into correct class

  tree <- treeio::as.treedata(ggtree::read.tree(text = nwk_list$tree.text))

  #just in case
  all_cells <- all_cells %>%
    dplyr::arrange(index)

  tree_data <- data.frame("X_coord" = all_cells[nwk_list$node_vec,"norm_locx"],
                          "Y_coord" = all_cells[nwk_list$node_vec,"norm_locy"],
                          "index" = as.character(nwk_list$node_vec))

  tree_data <- tree_data[match(c(tree@phylo$tip.label,tree@phylo$node.label), tree_data$index),]
  tree_data$node <- 1:nrow(tree_data)

  tree@data <- as_tibble(tree_data)

  # to include edge states and fraction of time spent on edge in treedata
  if (include_edge_states) {

    # if already found fraction of time on edge
    if ("frac_time_on_edge" %in% colnames(sampled_cells)) {

      tree@data$frac_time_on_edge <- sampled_cells[match(sampled_cells$index, tree@data$index),"frac_time_on_edge"]

    # otherwise calculate
    } else {

      print("Calculating lineage edge time")
      tree@data$frac_time_on_edge <- purrr::map_dbl(tree@data$node,
                                                  function(node) calc_segment_edge_time(node,
                                                                                        tree = tree,
                                                                                        sim_cells = all_cells))
    }

    # if already found edge state at cell death
    if ("death_edge_state" %in% colnames(sampled_cells)) {

      print("Finding node states")
      tree@data$death_edge_state <- sampled_cells[match(sampled_cells$index, tree@data$index),"death_edge_state"]

    # otherwise calculate
    } else {

      tree@data$death_edge_state <- purrr::map_dbl(tree@data$index, function(index) find_cell_state_at_death(index = index, all_cells = all_cells))
    }
  }

  #get muts again to include in metadata as allele info
  all_muts <- define_mut_presence_absence(all_cells) %>%
    as.data.frame %>%
    tibble::add_column("index" = all_cells$index) %>%
    dplyr::arrange(index) %>%
    dplyr::select(-index)



  sampled_muts <- unlist(purrr::map(1:nrow(sampled_cells), function(i) get_row_muts_by_i(i,sampled_cells$mutations))) %>%
    unique %>%
    sort

  sampled_mut_cols <- c(sampled_muts)

  muts_subset <- all_muts[,sampled_mut_cols] %>%
    tidyr::unite("mut", 1:length(sampled_mut_cols), sep = "")

  tree@data$allele <- muts_subset[tree_data$index,]
  tree@treetext <- nwk_list$tree.text



  return(tree)
}

#' Convert simulation to tree
#'
#' Convert data.frame of simulated cell to g to ggtree object with location data.
#'
#' @param all_cells data.frame
#' @param sampled_cells data.frame
#' @param branch_unit string
#' @param include_edge_states logical
#'
#' @return tree
#' @export
convert_sim_to_tree <- function(sampled_cells, all_cells, branch_unit= "time", include_edge_states = FALSE) {
  sim_nwk_list <- convert_nodes_to_string(sampled_cells$index,
                                          all_cells = all_cells,
                                          branch_unit = branch_unit)
  tree <- convert_nwk_to_treedata(sim_nwk_list, all_cells = all_cells, sampled_cells = sampled_cells,
                                  include_edge_states = include_edge_states)
  return(tree)
}

#' Label node types
#'
#' Label leaves and nodes to easily separate
#'
#' @param tree treedata
#'
#' @return treedata
#' @export
label_node_types <- function(tree) {
  tree@data$type <- ifelse(as.integer(tree@data$node) <= length(tree@phylo$tip.label), "leaf", "node")
  return(tree)
}

#' Calculate frequency spectrums
#'
#' Calculate frequency spectrums from all simulations with given death rate
#'
#' @param muts data.frame
#' @param death_rate numeric
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#'
#' @return data.frame
#' @export
calc_freq_spectrums <- function(muts, death_rate) {
  freqs <- sort(colMeans(muts), decreasing = TRUE)
  df <- data.frame("freq" = freqs,
                   "death_rate" = death_rate,
                   "mut_num_above" = 1:length(freqs))
  max_mut_num_above <- max(df$mut_num_above)
  df <- df %>%
    dplyr::mutate("max_mut_num_above" = max_mut_num_above)
  return(df)
}

#' Plot spatial phylogeny
#'
#' Plot tree on x-y coordinates.
#'
#' @param tree treedata
#' @param spatial_file string
#' @importFrom ape as.phylo
#' @importFrom dplyr arrange filter select
#' @importFrom phytools phylomorphospace
#' @importFrom treeio get.data
#'
#' @return NULL
#' @export
plot_space_phylo <- function(tree,
                             spatial_file = "spatial.png") {

  n <- length(tree@phylo$tip.label)
  m <- tree@phylo$Nnode

  #matrix with trait info for leaves
  X <- treeio::get.data(tree) %>%
    dplyr::arrange(node) %>%
    dplyr::filter(node %in% 1:n) %>%
    dplyr::select(X_coord, Y_coord) %>%
    as.matrix

  rownames(X) <- tree@phylo$tip.label

  #matrix with trait info for nodes
  A <- treeio::get.data(tree) %>%
    dplyr::arrange(node) %>%
    dplyr::filter(node %in% (n+1):(n+m)) %>%
    dplyr::select(X_coord, Y_coord) %>%
    as.matrix

  rownames(A) <- as.character((n+1):(n+m))

  #make continuous color scale

  color.gradient <- function(x, colors=c("red","yellow","springgreen","royalblue"), colsteps=100) {
    return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }

  #fit color scale to Euclidian distance from center
  col_vec <- c(color.gradient(sqrt(X[,1]^2 + X[,2]^2)), color.gradient(sqrt(A[,1]^2 + A[,2]^2)))
  names(col_vec) <- 1:(n+m)

  #make spatial only figure
  plot.new()
  png(filename = spatial_file)
  phytools::phylomorphospace(
    ape::as.phylo(tree),
    X,
    A=A,
    control=list(col.node=col_vec),
    xlab = "X_coord",
    ylab = "Y_coord",
    label='off'
  )
  dev.off()
  #make combined figure
  plot.new()

  reset.par()

  #dev.copy(png,"test.png")
  dev.new()
  par(xpd = NA, # switch off clipping, necessary to always see axis labels
      bg = "transparent", # switch off background to avoid obscuring adjacent plots
      oma = c(2, 2, 0, 0)) # move plot to the right and up

  phytools::phylomorphospace(
    ape::as.phylo(tree),
    X,
    A=A,
    control=list(col.node=col_vec),
    xlab = "X_coord",
    ylab = "Y_coord",
    label='off'
  )

  p1 <- recordPlot()

  dev.off()

}


#' Calculate LBI
#'
#' Calculate LBI statistic on tree nodes.
#'
#' @param tree treedata
#' @param muts data.frame
#' @param tau numeric
#' @importFrom ape as.phylo
#' @importFrom dplyr filter mutate percent_rank
#' @importFrom magrittr %>%
#' @importFrom treeImbalance lbi
#'
#' @return treedata
#' @export
calc_lbi <- function(tree, muts, tau = 0){
  if (tau == 0) {
    tau <- 1/nrow(muts)
  }
  lbi_values <- treeImbalance::lbi(ape::as.phylo(tree), tau = tau)
  tree@data$lbi <- lbi_values
  leaf_rank <- tree@data %>%
    dplyr::filter(type == "leaf") %>%
    dplyr::mutate("lbi_rank" = dplyr::percent_rank(lbi_values))
  node_rank <- tree@data %>%
    dplyr::filter(type == "node") %>%
    dplyr::mutate("lbi_rank"= dplyr::percent_rank(lbi_values))
  tree@data$lbi_rank <- c(leaf_rank$lbi_rank, node_rank$lbi_rank)
  return(tree)
}

#' Plot MCBD Tree
#'
#' Plot an MCC tree with median lambda / mu values on a coloured scale
#' @param treefile file containing the MCC tree
#' @param plotfile output file. If NULL, plots are sent to the R console.
#' @param scalel min and max values of lambda to use for the colour scale. If NULL, calculated from the tree values.
#' @param scalem min and max values of mu to use for the colour scale. If NULL, calculated from the tree values.
#' @param colour_gradient colours used on the scale, in order from low to high values.
#' @importFrom ape plot.phylo
#' @importFrom treeio read.beast
#'
#' @export
MCC_colour_plot = function(treefile, plotfile = NULL, scalel = NULL, scalem = NULL,
                           colour_gradient = c("red","yellow","green")) {f

  treedata = treeio::read.beast(treefile)
  lambdas = treedata@data$lambda_median[order(as.integer(treedata@data$node))]
  mus = treedata@data$mu_median[order(as.integer(treedata@data$node))]

  if(is.null(scalel)) {
    maxl = signif(1.05*max(lambdas),2)
    minl = signif(0.95*min(lambdas),2)
  }
  else {
    maxl = scalel[2]
    minl = scalel[1]
  }
  intervalsl = seq(minl, maxl, 0.1*(maxl - minl))

  if(is.null(scalem)) {
    maxm = signif(1.05*max(mus),2)
    minm = signif(0.95*min(mus),2)
  }
  else {
    maxm = scalem[2]
    minm = scalem[1]
  }
  intervalsm = seq(minm, maxm, 0.1*(maxm - minm))

  if(!is.null(plotfile)) pdf(plotfile, width = 17/grDevices::cm(1), height = 10/grDevices::cm(1))
  layout(mat = matrix(1:2, byrow = T, nrow = 1), widths = c(lcm(15),lcm(2)), heights = lcm(10))

  colsl = colour.gradient(lambdas[treedata@phylo$edge[,2]], intervals = intervalsl, colours = colour_gradient)
  ape::plot.phylo(treedata@phylo,type = "fan",edge.color = colsl,show.tip.label = F, edge.width = 2, no.margin = T)

  legendplot(intervalsl, colours = colour_gradient)

  colsm = colour.gradient(mus[treedata@phylo$edge[,2]], intervals = intervalsm, colours = colour_gradient)
  ape::plot.phylo(treedata@phylo,type = "fan",edge.color = colsm,show.tip.label = F, edge.width = 2, no.margin = T)

  legendplot(intervalsm, colours = colour_gradient)

  if(!is.null(plotfile)) dev.off()
}

legendplot = function(intervals = seq(0,11,0.5), axis.int = NULL, colours = c("red","yellow","green")) {
  lut = colorRampPalette(colours) (length(intervals))
  if(is.null(axis.int)) axis.int = intervals

  scale = (length(lut)-1)/(max(intervals)-min(intervals))
  plot(c(0,10), c(min(intervals),max(intervals)), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")
  axis(2, axis.int, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min(intervals)
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

colour.gradient <- function(x, intervals = seq(0,11,0.1), colours=c("red","yellow","green")) {
  colfun = colorRampPalette(colours)
  return(  colfun(length(intervals)) [ findInterval(x, intervals) ] )
}

#' Calc distance from center
#'
#' Calculate distance of cell from tumor center.
#' @param tree treedata
#' @importFrom purrr map2_dbl
#' @importFrom stats dist
#' @importFrom dplyr mutate
#'
#' @return treedata
#'
#' @export
calc_dist_from_center <- function(tree) {




  tree@data <- tree@data %>%
    dplyr::mutate("distance_from_center" = purrr::map2_dbl(X_coord, Y_coord, function(x,y) dist(matrix(c(0,0,x,
                                                                                                  y),
                                                                                                nrow = 2,
                                                                                                byrow = TRUE)))
    )
  return(tree)
}

#' Calc distance from center
#'
#' Calculate distance of cell from tumor center.
#' @param tree treedata
#' @param sp1 integer
#' @param sp2 integer
#'
#' @return numeric
#'
#' @export
calc_spatial_distance_given_leaves <- function(tree,sp1, sp2) {
  row1 <- which(tree@data$index == sp1) #match index to nodes
  row2 <- which(tree@data$index == sp2)
  x1 <- tree@data$X_coord[row1] #get spatial coordinates
  x2 <- tree@data$X_coord[row2]

  y1 <- tree@data$Y_coord[row1]
  y2 <- tree@data$Y_coord[row2]

  return(sqrt((x1 - x2)^2 + (y1 - y2)^2)) #return euclidian distance
}

#' Make tree from simulated cells
#'
#' Make tree from simulated cells including spatial constraints.
#' @param sampled_cells data.frame
#' @param sim_cells data.frame
#' @param add_cell_to_label logical. Default TRUE
#' @param branch_unit string. Default "time". Other options are "molecular", "none".
#' @param include_edge_states string
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join
#'
#' @return numeric
#'
#' @export
make_tree_from_simulated_cells <- function(sampled_cells,
                                           sim_cells,
                                           add_cell_to_label = TRUE,
                                           branch_unit = "time",
                                           include_edge_states = TRUE) {
  sim_cells <- sim_cells %>%
    normalize_locs()


  sampled_cells <- dplyr::left_join(sampled_cells, sim_cells, by = c("alpha", "birthdate", "cells_free",
                                                                     "deathdate", "index", "locx", "locy",
                                                                     "mutation_rate", "mutations", "parent_index",
                                                                     "proliferation_rate", "sequence"))
  print("Making newick tree")
  sim_nwk_list <- convert_nodes_to_string(sampled_cells$index,
                                          sim_cells,
                                          branch_unit = branch_unit)
  print("Adding metadata")
  tree <- convert_nwk_to_treedata(nwk_list = sim_nwk_list,
                                  all_cells = sim_cells,
                                  sampled_cells = sampled_cells,
                                  include_edge_states = include_edge_states)

  #tree <- as.treedata(read.tree(text = sim_nwk_list$tree.text))

  if (add_cell_to_label) {
    tree@phylo$tip.label <- paste0("cell_", tree@phylo$tip.label)
    tree@phylo$node.label <- paste0("cell_", tree@phylo$node.label)
  }

  return(tree)
}

#' Cluster cells
#'
#' Cluster cells from genetic data
#' @param sampled_cells data.frame
#' @param k data.frame
#'
#' @importFrom magrittr %>%
#' @importFrom stats kmeans
#' @importFrom tibble add_column
#'
#' @return numeric
#'
#' @export
cluster_cells <- function(sampled_cells, k = 10) {

  pa <- as.data.frame(define_mut_presence_absence(sampled_cells))

  clusters <- stats::kmeans(pa, centers = k)$cluster
  clustered_cells <- sampled_cells %>%
    tibble::add_column("cluster" = clusters)

  return(clustered_cells)
}


convert_nodes_to_string_binary_multitype<- function(orphan_cells,
                                                        all_cells) {


  #should be already sorted by index, but just to make sure
  all_cells <- all_cells %>% dplyr::arrange(index)


  #start out with each leaf as its own subtree
  subtrees <- as.character(orphan_cells)

  #keep track of all nodes in tree
  node_vec <- c()
  n <- length(orphan_cells)

  #also keep track of time spent in state 1

  times_in_state_1_vec <- c()
  total_times_vec <- c()

  node_states <- c()

  while(length(subtrees) > 1) { #will connect subtrees by MCRA until reach single tree

    #make all pariwise combinations of cells without mrca
    all_combs <- as.data.frame(t(combn(orphan_cells, 2)))

    #find ancestors of candidate cells
    common_ancestors <- purrr::map2_dbl(all_combs$V1, all_combs$V2, function(index_1, index_2) find_MRCA(index_1 = index_1,
                                                                                                         index_2 = index_2,
                                                                                                         all_cells = all_cells))
    common_ancestors_data <- all_cells[common_ancestors,]

    #find youngest ancestor and two children
    youngest_ancestor <- common_ancestors_data$index[common_ancestors_data$birthdate == max(common_ancestors_data$birthdate)]

    #in case more than 2 cells have same ancestor,
    #choose one to start, we will come back to the other
    youngest_ancestor <- youngest_ancestor[1]


    #keep track of children descended from that ancestor

    children <- unique(unlist(all_combs[which(common_ancestors == youngest_ancestor),]))


    children_substrings <- c()
    for (child in children) {

      #c <- tail(which(node_vec == child), n = 1)


      t_a <- all_cells[youngest_ancestor, "deathdate"]
      t_c <- all_cells[child, "deathdate"] #use deathdate


      #correct branch lengths for state-dependent clock rate

      ## find all ancestors and filter for those within branch (between identified mrca and child node)
      child_all_ancestors <- c(collect_all_ancestors(child, all_cells), child)
      all_cells_on_branch_df <- all_cells[all_cells$index %in% child_all_ancestors & all_cells$birthdate >= t_a, ]

      #calculate start and end points
      branch_start_time <- min(all_cells_on_branch_df$birthdate)
      branch_end_time <- max(all_cells_on_branch_df$deathdate)
      total_branch_time <- branch_end_time - branch_start_time


      #extract record of states for cells on branch
      state_rec <- as.integer(unlist(strsplit(purrr::map_chr(all_cells_on_branch_df$state_rec, function(sr) gsub("[^0-9]", "", sr)),split = "")))


      #get starting state
      start_state <- state_rec[1]




      state_transitions <- unlist(purrr::map(all_cells_on_branch_df$state_transitions, function(st) as.numeric(gsub("\\[|\\]", "", strsplit(as.character(st),",")[[1]]))))
      state_transitions <- state_transitions[!is.na(state_transitions)]

      time_period_states <- (1- start_state) * rep_len(c(0,1), length.out = (length(state_transitions) + 1)) +
        start_state * rep_len(c(1,0), length.out = (length(state_transitions) + 1))


      node_states <- c(node_states, time_period_states)

      time_breaks <- c(branch_start_time, state_transitions, branch_end_time)
      time_periods <- diff(time_breaks)


      if (clock_scaled) {
        genetic_distance <- state_dependent_clock_rates[2] * time_periods * time_period_states +
          state_dependent_clock_rates[1] * time_periods * (1 - time_period_states)
      } else {

        genetic_distance <- time_periods * time_period_states + time_periods * (1 - time_period_states)
      }


      time_spent_in_state1 = sum(time_periods * time_period_states)


      #fraction_time_spent_state1 <- time_spent_in_state1 / total_branch_time

      times_in_state_1_vec <- c(times_in_state_1_vec, time_spent_in_state1)
      total_times_vec <- c(total_times_vec, total_branch_time )

      #branch length closest to child is last branch
      child_bl <- round(genetic_distance[length(genetic_distance)],4)


      #branch length for children in tree
      #children_branch_lengths <- c(children_branch_lengths, as.character(round(genetic_distance[length(genetic_distance)],4)))

      #if branches have more than one type, need to re-name node into subnodes and keep branch lengths separately

      if (length(time_period_states) > 1) {

        multi_type_nodes <- c(paste0(child, letters[length(time_period_states):2]), child)

        node_vec <- c(node_vec, multi_type_nodes)


        child_substring <- paste(subtrees[which(orphan_cells == child)], as.character(child_bl), sep = ":")

        for (j in 2:length(multi_type_nodes)) {
          child_substring <- paste(paste("(", child_substring, ")", rev(multi_type_nodes)[j], sep = "") , round(rev(genetic_distance)[j],4), sep = ":")
        }

      } else {

        node_vec <- c(node_vec, child)

        child_substring <- paste(subtrees[which(orphan_cells == child)], as.character(child_bl), sep = ":")
      }

      children_substrings <- c(children_substrings, child_substring)



    }

    #subtrees <- c(subtrees, paste0("(",paste(paste(subtrees[which(orphan_cells %in% children)],branch_lengths,sep=":"), collapse = ","),")",as.character(youngest_ancestor)))
    subtrees <- c(subtrees, paste("(", paste(children_substrings, collapse = ","),")",as.character(youngest_ancestor), sep = "", collapse = ""))





    #make new subtree from combining children and MCRA

    subtrees <- subtrees[-which(orphan_cells %in% children)] #remove child trees after used

    orphan_cells <- c(orphan_cells, youngest_ancestor) #add ancestor as orphan

    orphan_cells <- orphan_cells[-which(orphan_cells %in% children)] #remove children
  }

  #to include root

  #if (length(orphan_cells) == 1) {

  child <- youngest_ancestor #get last root cell

  t_root <- 0
  t_c <- all_cells[all_cells$index == child, "deathdate"] #use deathdate


  #correct branch lengths for state-dependent clock rate

  ## find all ancestors and filter for those within branch (between identified mrca and child node)
  child_all_ancestors <- c(collect_all_ancestors(child, all_cells), child)
  all_cells_on_branch_df <- all_cells[all_cells$index %in% child_all_ancestors & all_cells$birthdate >= t_root, ]

  #calculate start and end points
  branch_start_time <- min(all_cells_on_branch_df$birthdate)
  branch_end_time <- max(all_cells_on_branch_df$deathdate)
  total_branch_time <- branch_end_time - branch_start_time


  #extract record of states for cells on branch
  state_rec <- as.integer(unlist(strsplit(purrr::map_chr(all_cells_on_branch_df$state_rec, function(sr) gsub("[^0-9]", "", sr)),split = "")))


  #get starting state
  start_state <- state_rec[1]

  root_state <- start_state




  state_transitions <- unlist(purrr::map(all_cells_on_branch_df$state_transitions, function(st) as.numeric(gsub("\\[|\\]", "", strsplit(as.character(st),",")[[1]]))))
  state_transitions <- state_transitions[!is.na(state_transitions)]

  time_period_states <-  time_period_states <- (1- start_state) * rep_len(c(0,1), length.out = (length(state_transitions) + 1)) +
    start_state * rep_len(c(1,0), length.out = (length(state_transitions) + 1))




  node_states <- c(node_states, time_period_states)

  time_breaks <- c(branch_start_time, state_transitions, branch_end_time)
  time_periods <- diff(time_breaks)


  if (clock_scaled) {
    genetic_distance <- state_dependent_clock_rates[2] * time_periods * time_period_states +
      state_dependent_clock_rates[1] * time_periods * (1 - time_period_states)
  } else {
    genetic_distance <- time_periods * time_period_states + time_periods * (1 - time_period_states)
  }


  time_spent_in_state1 = sum(time_periods * time_period_states)


  #fraction_time_spent_state1 <- time_spent_in_state1 / total_branch_time

  times_in_state_1_vec <- c(times_in_state_1_vec, time_spent_in_state1)
  total_times_vec <- c(total_times_vec, total_branch_time )

  #branch length closest to child is last branch
  child_bl <- round(genetic_distance[length(genetic_distance)],4)


  #branch length for children in tree
  #children_branch_lengths <- c(children_branch_lengths, as.character(round(genetic_distance[length(genetic_distance)],4)))

  #if branches have more than one type, need to re-name node into subnodes and keep branch lengths separately

  if (length(time_period_states) > 1) {

    multi_type_nodes <- c(paste0(child, letters[length(time_period_states):2]), child)

    node_vec <- c(node_vec, multi_type_nodes)


    child_substring <- paste(subtrees[1], as.character(child_bl), sep = ":")

    for (j in 2:length(multi_type_nodes)) {
      child_substring <- paste(paste("(", child_substring, ")", rev(multi_type_nodes)[j], sep = "") , round(rev(genetic_distance)[j],4), sep = ":")
    }

  } else {

    node_vec <- c(node_vec, child)

    child_substring <- paste(subtrees[1], as.character(child_bl), sep = ":")
  }



  #add roots
  subtrees <-  paste("(", paste(child_substring, collapse = ","),")","root", sep = "", collapse = "")

  node_vec <- c(node_vec, "root")
  node_states <- c(node_states, root_state)









  return(list("tree.text"= paste0(tail(subtrees, n = 1),";"),
              "node_vec" = node_vec,
              "node_states" = node_states))


}

#modified function to inlude node states
convert_nwk_to_treedata_bd <- function(nwk_list) {


  tree <- treeio::as.treedata(ggtree::read.tree(text = nwk_list$tree.text))

  #just in case
  # all_cells <- all_cells %>%
  #     dplyr::arrange(index)

  tree_data <- data.frame("index" = as.character(nwk_list$node_vec),
                          "state" = nwk_list$node_states)

  tree_data <- tree_data[match(c(tree@phylo$tip.label,tree@phylo$node.label), tree_data$index),]
  tree_data$node <- 1:nrow(tree_data)

  tree@data <- as_tibble(tree_data)


  return(tree)
}

#' Make full tree from simulated cells fast
#'
#' Make tree from simulated cells including spatial constraints but faster by converted to nwk string first.
#' Use prune_simulated tree to convert to subsampled tree. 
#' 
#' @param all_cells data.frame of all simulated cells in tumor
#' @param add_all_states logical. Default FALSE. If included edge/center states, which takes more time. 
#' @param sampled_cells_indices Vector of indices of sampled cells. If NULL, all cells will be included in tree. 
#' @param branch_unit string. Default "time". Other options are "genetic" and "none".
#' @param cell_locations_df_file File to find locations through time for pushing simulation. 
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join mutate
#' @importFrom tibble tibble
#' @importFrom treeio as.treedata
#' @importFrom ape read.tree
#' @importFrom phangorn Ancestors
#'
#' @return numeric
#'
convert_all_cells_to_tree_fast <- function (all_cells, add_all_states = FALSE, sampled_cells_indices = NULL, 
          branch_unit = "time", cell_locations_df_file = NULL){
  
  get_genetic_branch_length <- function(cell_index, all_cells) {
    mut_curr_cell <- all_cells$collapsed_sequence[all_cells$index == 
                                                    cell_index]
    parent_cell_index <- all_cells$parent_index[all_cells$index == 
                                                  cell_index]
    mut_parent_cell <- all_cells$collapsed_sequence[all_cells$index == 
                                                      parent_cell_index]
    return(adist(mut_curr_cell, mut_parent_cell)[1])
  }
  if (branch_unit == "time") {
    all_cells <- all_cells %>% dplyr::mutate(cell_name = paste0("cell_", 
                                                                index, sep = ""), branch_length = round(deathdate - 
                                                                                                          birthdate, 2)) %>% dplyr::mutate(nwk_node = paste0(cell_name, 
                                                                                                                                                             ":", branch_length))
  }
  else if (branch_unit == "genetic") {
    
    #print("Computing genetic distance")
    all_cells$collapsed_sequence <- purrr::map_chr(all_cells$sequence,
                                                   extract_sequences)
    

    all_cells$branch_length <- purrr::map_dbl(all_cells$index, 
                                              function(i) get_genetic_branch_length(i, all_cells = all_cells))
    all_cells <- all_cells %>% dplyr::mutate(cell_name = paste0("cell_", 
                                                                index, sep = "")) %>% dplyr::mutate(nwk_node = paste0(cell_name, 
                                                                                                                      ":", branch_length))
  }
  else {
    all_cells <- all_cells %>% dplyr::mutate(cell_name = paste0("cell_", 
                                                                index, sep = ""), branch_length = NA) %>% dplyr::mutate(nwk_node = paste0(cell_name, 
                                                                                                                                          ":", branch_length))
  }
  curr_leaves <- c(1)
  curr_nwk_string <- paste0("(", all_cells$nwk_node[which(all_cells$index == 
                                                            1)], ");")
  if (!is.null(sampled_cells_indices)) {
    included_cell_indices <- c(unique(unlist(purrr::map(sampled_cells_indices, 
                                                        function(child) collect_all_ancestors(index = child, 
                                                                                              all_cells = all_cells)))), sampled_cells_indices)
  }
  else {
    included_cell_indices <- all_cells$index
  }
  while (length(curr_leaves) > 0) {
    new_leaves <- c()
    for (leaf in curr_leaves) {
      children_indices <- all_cells$index[which((all_cells$parent_index == 
                                                   leaf) & (all_cells$index %in% included_cell_indices))]
      if (length(children_indices) > 0) {
        parent_node_string <- all_cells$nwk_node[which(all_cells$index == 
                                                         leaf)]
        children_node_string <- paste0("(", paste(all_cells$nwk_node[which(all_cells$parent_index == 
                                                                             leaf & (all_cells$index %in% included_cell_indices))], 
                                                  collapse = ","), ")", sep = "")
        new_node <- paste0(children_node_string, parent_node_string)
        curr_nwk_string <- gsub(parent_node_string, new_node, 
                                curr_nwk_string, perl = TRUE)
        new_leaves <- c(new_leaves, children_indices)
      }
    }
    curr_leaves <- new_leaves
  }
  tree <- treeio::as.treedata(ape::read.tree(text = curr_nwk_string))
  endpoint <- max(all_cells$deathdate)
  tree@data <- tibble(node = 1:(length(tree@phylo$node.label) + 
                                  length(tree@phylo$tip.label)), cell_name = c(tree@phylo$tip.label, 
                                                                               tree@phylo$node.label)) %>% dplyr::left_join(all_cells, 
                                                                                                                            by = "cell_name") %>% dplyr::mutate(alive = deathdate == 
                                                                                                                                                                  endpoint)
  if (add_all_states) {
    alive_nodes <- tree@data$node[tree@data$deathdate == 
                                    max(tree@data$deathdate, na.rm = TRUE)]
    alive_nodes <- alive_nodes[!is.na(alive_nodes)]
    non_exinct_nodes <- c()
    for (n in alive_nodes) {
      non_exinct_nodes <- c(non_exinct_nodes,
                            phangorn::Ancestors(tree@phylo, 
                                                node = n, type = "all"))
    }
    states_df <- find_all_cell_state_at_death(all_cells,
                                              cell_locations_df_file = cell_locations_df_file) %>% 
      dplyr::select(index, state)
    tree@data <- tree@data %>% dplyr::left_join(., states_df, by = "index")
    tree@phylo$edge.length[is.na(tree@phylo$edge.length)] <- 0
  }
  return(tree)
}

#' Prune simulated tree
#'
#' Prune full simulated tree to subsampled tips
#' 
#' @param tree treeio S4 tree of simulated tumor 
#' @param sampled_cells_indices Vector of indices of sampled cells. If NULL, all cells will be included in tree. 
#' @param add_all_states logical. Default FALSE. If included edge/center states, which takes more time. 
#' @param all_cells data.frame of all simulated cells. Default NULL but necessary to get edge/center states. 
#' @param cell_locations_df_file File to find locations through time for pushing simulation. 
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join mutate
#' @importFrom tibble tibble
#' @importFrom treeio as.treedata
#' @importFrom ape read.tree
#' @importFrom phangorn Ancestors
#' @importFrom stringr str_count
#'
#' @return numeric

#' @importFrom ape read.tree drop.tip
prune_simulated_tree <- function (tree, sampled_cells_indices, add_all_states = FALSE, 
          all_cells = NULL, cell_locations_df_file = NULL, branch_unit = "time"){
  sampled_cells <- paste0("cell_", sampled_cells_indices, sep = "")
  pruned.phylo <- ape::drop.tip(tree@phylo, tree@phylo$tip.label[-match(sampled_cells, 
                                                                        tree@phylo$tip.label)],
                                collapse.singles = TRUE)
  pruned.tree <- treeio::as.treedata(pruned.phylo)
  pruned.tree@data <- tibble(cell_name = c(pruned.phylo$tip.label, 
                                           pruned.phylo$node.label)) %>% dplyr::left_join(., tree@data, 
                                                                                          by = "cell_name") %>% dplyr::mutate(node = 1:(length(pruned.phylo$node.label) + 
                                                                                                                                          length(pruned.phylo$tip.label)))
  
  #manually add root
  if (branch_unit == "time") {
    
    #get deathdate of MRCA of all sampled cells
    pruned.tree@phylo$root.edge <- min(pruned.tree@data$deathdate)
    
  } else {
    
    #get number of mutations accumulated by MRCA of all sampled cells
    pruned.tree@phylo$root.edge <- str_count(pruned.tree@data$mutations[which.min(pruned.tree@data$deathdate)], pattern = ",")
  }
  
  
  if (add_all_states) {
    if (is.null(all_cells)) {
      error("To add all states need data.frame of all cells in simulation")
    }
    states <- purrr::map_dbl(pruned.tree@data$index, function(index) tumortree::find_cell_state_at_death(index = index, 
                                                                                                         all_cells = all_cells,
                                                                                                         cell_locations_df_file = cell_locations_df_file))
    pruned.tree@data$state <- states
  }
  return(pruned.tree)
}

