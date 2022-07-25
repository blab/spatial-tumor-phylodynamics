#' Get free cells at time t
#'
#' Given index, get the number of free surrounding lattice spots of cells at time t.
#'
#' @param index integer
#' @param time numeric
#' @param all_cells data.frame
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#'
#' @return integer
#' @export
get_free_cells_at_t <- function(index,time, all_cells){
    alive_cells <- all_cells %>%
        dplyr::filter(birthdate <= time & deathdate >= time) #inclusive on both sides
    x_loc <- alive_cells[alive_cells$index == index,]$locx
    y_loc <- alive_cells[alive_cells$index == index,]$locy
    curr_free_cells <- get_n_adjacent_cells(x_loc, y_loc, alive_cells = alive_cells)
    return(curr_free_cells)
}

#' Get n adjacent cells
#'
#' Given x-y location and current alive cells, find how many free adjacent lattice spots there are.
#'
#' @param x_loc numeric
#' @param y_loc numeric
#' @param all_cells data.frame
#'
#' @return data.frame
#' @export
get_n_adjacent_cells <- function(x_loc, y_loc, alive_cells) {
    check_x_vec <- c(x_loc + 1, x_loc, x_loc - 1, x_loc - 1, x_loc + 1, x_loc - 1, x_loc, x_loc + 1)
    check_y_vec <- c(y_loc + 1, y_loc + 1, y_loc + 1, y_loc, y_loc, y_loc - 1, y_loc - 1, y_loc - 1)
    n_adjacent_cells <- 8
    for (i in 1:length(check_x_vec)) {
        check_x <- check_x_vec[i]
        check_y <- check_y_vec[i]
        if (sum(alive_cells$locx == check_x & alive_cells$locy == check_y) == 1) {
            n_adjacent_cells <-  n_adjacent_cells - 1
        }
    }
    return(n_adjacent_cells)
}

#' Count spatial populations at t
#'
#' Count population of spatial comparments in tumor at time t.
#'
#' @param all_cells data.frame
#' @param time numeric
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom purrr map2_dbl
#'
#' @return data.frame
#' @export
count_spatial_pops_at_t <- function(all_cells, time){
    alive_cells <- all_cells %>%
        dplyr::filter(birthdate <= time & deathdate > time)
    curr_free_cells <- purrr::map2_dbl(alive_cells$locx,
                                       alive_cells$locy,
                                       function(x_loc, y_loc) get_n_adjacent_cells(x_loc, y_loc, alive_cells = alive_cells))
    n_cells_per_pop <- purrr::map_dbl(0:7, function(fc) sum(curr_free_cells == fc))
    df <- data.frame("cells_free" = 0:7,
                     "time" = time,
                     "n_cells" = n_cells_per_pop)
    return(df)
}

#' Get spatial populations at over time
#'
#' Create time points and get populations over time
#'
#' @param all_cells data.frame
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows
#' @importFrom purrr map
#'
#' @return data.frame
#' @export
get_cell_pops_over_time <- function(all_cells) {
    #set up intervals for each time to count alive cells
    endpoint <- max(all_cells$deathdate)
    timestep <- 1/24
    time_seq <- seq(0, (endpoint - timestep), by = timestep)

    pop_counts_list <- purrr::map(time_seq, function(time) count_spatial_pops_at_t(all_cells, time))
    #count up total alive cells

    pop_counts_df <- pop_counts_list %>%
        dplyr::bind_rows

    return(pop_counts_df)
}

#' Get all cell states
#'
#' Get spatial state of cell over lifespan.
#'
#' @param index integer
#' @param all_cells data.frame
#' @importFrom purrr map_dbl
#'
#' @return data.frame
#' @export
get_all_states_cell <- function(index,all_cells){
    birthdate <- all_cells[all_cells$index == index, ]$birthdate
    deathdate <- all_cells[all_cells$index == index, ]$deathdate
    time_alive <- seq(birthdate, deathdate, by = 1/24)
    time_alive <- time_alive[-length(time_alive)] # don't include death date
    free_cells_over_time <- purrr::map_dbl(time_alive, function(t) get_free_cells_at_t(index= index, time = t, all_cells = all_cells))
    return(data.frame("state" = free_cells_over_time, "time" = time_alive))
}

#' Get all cell states of lineage
#'
#' Get spatial states of full lineage over history.
#'
#' @param leaf integer
#' @param ancestor integer
#' @param all_cells data.frame
#' @importFrom dplyr bind_rows
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom tibble add_column
#'
#' @return data.frame
#' @export
get_lineage_states <- function(leaf, ancestor = 1, all_cells) {

    #first get all ancestors for leaf (that define end of lineage from ancestor (default is root))
    #because a dividing cells is replace by two new indexed cells, just need to find all ancestors in between
    #desired ancestor and leaf and for each cell's lifespan calculate how much time it spend on edge and the center
    lineage_cells <- collect_all_ancestors(leaf, all_cells = all_cells)
    lineage_cells <- c(leaf, lineage_cells[1:which(lineage_cells == ancestor)])

    lineage_states <- purrr::map(lineage_cells,  function(index) get_all_states_cell(index = index, all_cells = all_cells)) %>%
        dplyr::bind_rows() %>%  #total time include for qc
        tibble::add_column("leaf_index" = leaf)

    return(lineage_states)

}

#' Check if cell is most extreme.
#'
#' Check if cell is at an extreme x or y location in the tumor to determine edge.
#'
#' @param x_loc numeric
#' @param y_loc numeric
#' @param alive_cells data.frame
#'
#' @return integer
#' @export
check_if_cell_most_extreme <- function(x_loc, y_loc, alive_cells) {

    #first fix x location and ask if it is the most extreme on y axis

    #are there cells above?
    cells_above <- sum((alive_cells$locx == x_loc) & (alive_cells$locy > y_loc)) > 0

    #are there cells below?
    cells_below <- sum((alive_cells$locx == x_loc) & (alive_cells$locy < y_loc)) > 0

    #it is most extreme on this axis as long as only one of these statements are true
    most_extreme_on_y_axis <- ! (cells_above & cells_below)

    #second, fix y location and ask if it is the most extreme on x axis

    #are there cells above?
    cells_left <- sum((alive_cells$locx < x_loc) & (alive_cells$locy == y_loc)) > 0

    #are there cells below?
    cells_right <- sum((alive_cells$locx > x_loc) & (alive_cells$locy == y_loc)) >0

    #it is most extreme on this axis as long as only one of these statements are true
    most_extreme_on_x_axis <- ! (cells_left & cells_right)

    #on edge if is most extreme on either x or y axes
    most_extreme_on_any_axis <- most_extreme_on_x_axis |  most_extreme_on_y_axis


    return(as.integer(most_extreme_on_any_axis))
}

#' Check if cell is adjacent boundary
#'
#' Check if cell is next to an extreme x or y location in the tumor to determine edge.
#'
#' @param x_loc numeric
#' @param y_loc numeric
#' @param alive_cells data.frame
#'
#' @return integer
#' @export
check_boundary_adjacent <- function(x_loc, y_loc, alive_cells) {

    check_x_vec <- c(x_loc, x_loc - 1, x_loc + 1, x_loc)
    check_y_vec <- c(y_loc + 1, y_loc, y_loc, y_loc - 1)
    neighbor_cell_states <- c()
    for (i in 1:length(check_x_vec)) {
        check_x <- check_x_vec[i]
        check_y <- check_y_vec[i]
        neighbor_cell_states <- c(neighbor_cell_states,
                                  alive_cells$most_extreme[alive_cells$locx == check_x & alive_cells$locy == check_y])

    }

    return(as.integer(sum(neighbor_cell_states) > 0))
}

#' Find cell state at death
#'
#' Find binary edge/center spatial state at time of cell division or death
#'
#' @param index integer of cell index
#' @param all_cells data.frame of all simulated cells
#' @param cell_locations_df_file CSV file for locations of cells through time. If NULL, which is default, assume current locations have not changed. 
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom purrr map2_dbl
#'
#' @return integer
#' @export
find_cell_state_at_death <- function (index,
                                      all_cells,
                                      cell_locations_df_file = NULL) {
    print(index)
    time <- all_cells$deathdate[all_cells$index == index]
    alive_cells <- all_cells %>% dplyr::filter(birthdate < time & 
                                                   deathdate >= time)
    if (!is.null(cell_locations_df_file)) {
        
        print("Using locations file")
        cell_locations_df <- read_csv(cell_locations_df_file)
        cell_locations_df_t <- cell_locations_df %>% dplyr::filter(round(t,2) == round(time, 2))
        alive_cells <- alive_cells %>% dplyr::select(-locx, -locy) %>% 
            dplyr::left_join(., cell_locations_df_t, by = "index")
    }
    alive_cells <- mark_boundary(alive_cells, alive_cells)
    state_at_death <- alive_cells$est_edge[alive_cells$index == 
                                               index, "est_edge"]
    return(as.numeric(state_at_death))
}


#' Mark boundary
#'
#' Mark if set of cells are on boundary given alive cells.
#'
#' @param cells_to_mark data.frame
#' @param alive_cells data.frame
#' @importFrom purrr map2_dbl
#'
#' @return data.frame
#' @export
mark_boundary <- function(cells_to_mark, alive_cells){

    cells_to_mark$most_extreme <- purrr::map2_dbl(cells_to_mark$locx, cells_to_mark$locy,
                                              function(x_loc, y_loc) check_if_cell_most_extreme(x_loc, y_loc, alive_cells = alive_cells))

    alive_cells$most_extreme <- purrr::map2_dbl(alive_cells$locx, alive_cells$locy,
                                            function(x_loc, y_loc) check_if_cell_most_extreme(x_loc, y_loc, alive_cells = alive_cells))

    cells_to_mark$edge_adjacent <- purrr::map2_dbl(cells_to_mark$locx, cells_to_mark$locy,
                                               function(x_loc, y_loc) check_boundary_adjacent(x_loc, y_loc, alive_cells = alive_cells))
    cells_to_mark$est_edge <- as.integer(cells_to_mark$most_extreme | cells_to_mark$edge_adjacent)

    return(cells_to_mark)
}


#' Get lineage segment
#'
#' Collect all cells of branch
#'
#' @param node integer
#' @param tree treedata
#' @importFrom phangorn Ancestors
#'
#' @return vector
#' @export
get_lineage_segment <- function(node, tree, sim_cells) {
    parent <- phangorn::Ancestors(tree@phylo, node, type = "parent")
    daughter_index <- tree@data$index[tree@data$node == node]

    parent_index <- tree@data$index[tree@data$node == parent]

    lineage_cells <- c(collect_all_ancestors(daughter_index, all_cells = sim_cells))
    segment_cells <- lineage_cells[1:which(lineage_cells == parent_index)]
    return(segment_cells)
}

#' Get cell time on est edge
#'
#' Get time cells spends on esimated edge
#'
#' @param cell_index integer
#' @param sim_cells data.frame
#' @param timestep numeric
#' @importFrom dplyr filter
#' @importFrom purrr map_dbl
#' @importFrom magrittr %>%
#'
#' @return data.frame
#' @export
get_cell_time_on_est_edge <- function(cell_index, sim_cells, timestep =1/24) {
    cell_birthdate <- sim_cells$birthdate[sim_cells$index == cell_index]
    cell_deathdate <- sim_cells$deathdate[sim_cells$index == cell_index]
    cell_xloc <- sim_cells$locx[sim_cells$index == cell_index]
    cell_yloc <- sim_cells$locy[sim_cells$index == cell_index]
    times <- seq(cell_birthdate, cell_deathdate - timestep, by = timestep)
    edge_stat_vec <- c()
    for (t in times) {
        alive_cells <- sim_cells %>%
            dplyr::filter((birthdate < t) & (deathdate >= t))

        alive_cells$most_extreme <- purrr::map2_dbl(alive_cells$locx,
                                                    alive_cells$locy,
                                                    function(x,y) check_if_cell_most_extreme(x_loc = x,
                                                                                             y_loc =y,
                                                                                             alive_cells = alive_cells))

        cell_on_edge <- (check_boundary_adjacent(cell_xloc, cell_yloc, alive_cells) == 1) | (alive_cells$most_extreme[alive_cells$index == cell_index] == 1)
        edge_stat_vec <- c(edge_stat_vec, cell_on_edge)
    }

    return(data.frame("edge_time" = sum(edge_stat_vec) * timestep,
                      "total_time" = length(times) * timestep))
}

#' Calc segment edge time
#'
#' Calculate time spend on edge of branch
#'
#' @param node integer
#' @param tree treedata
#' @param sim_cells data.frame
#' @importFrom dplyr filter
#' @importFrom purrr map
#' @importFrom magrittr %>%
#'
#' @return numeric
#' @export
calc_segment_edge_time <- function(node, tree, sim_cells) {
    if(node == tree@phylo$Nnode + 2) {

        return(NA)

    } else {

        curr_index <- tree@data$index[tree@data$node == node]

        lineage_cells <- get_lineage_segment(node = node, tree = tree, sim_cells = sim_cells)
        segment_edge_hx <- purrr::map(c(curr_index, lineage_cells),
                                      function(cell_index) get_cell_time_on_est_edge(cell_index,
                                                                                     sim_cells = sim_cells)) %>%
            bind_rows()
        frac_of_seg_on_edge = sum(segment_edge_hx$edge_time) / sum(segment_edge_hx$total_time)
        return(frac_of_seg_on_edge)
    }
}

#' Get spatial compartment growth rates
#'
#' Extract instantaneous edge and center growth rates over time
#'
#' @param all_cells data.frame
#' @param time numeric
#' @param time_interval numeric
#' @importFrom dplyr filter
#' @importFrom purrr map2_dbl
#' @importFrom magrittr %>%
#'
#' @return numeric
#' @export
get_compartment_growth_rates_at_t <- function(all_cells, time, time_interval = 1/24){

    successful_parents <- sort(unique(all_cells$parent_index))

    #get all alive cells in time interval
    alive_cells <- all_cells %>%
        dplyr::filter(birthdate <= time & deathdate > time)

    #get "most extreme" cells that found base of boundary estimation
    alive_cells$most_extreme <- purrr::map2_dbl(alive_cells$locx, alive_cells$locy,
                                                function(x_loc, y_loc) check_if_cell_most_extreme(x_loc, y_loc, alive_cells = alive_cells))

    #extend boundary cells to those directly adjacent to extreme cells
    alive_cells$edge_adjacent <- purrr::map2_dbl(alive_cells$locx, alive_cells$locy,
                                                 function(x_loc, y_loc) check_boundary_adjacent(x_loc, y_loc, alive_cells = alive_cells))

    #define edge and boundary cells
    alive_cells$est_edge <- as.integer(alive_cells$most_extreme | alive_cells$edge_adjacent)


    #get compartment totals
    N_center_cells <- sum(!alive_cells$est_edge)
    N_edge_cells <- sum(alive_cells$est_edge)

    #get divisions in time interval
    #in this simulation a cell 'dies' which it divides,
    #so a division is marked by a "deathdate" in the the interval
    #with daughter cells remaining in the population

    #cells born in time interval
    alive_cells$successful_parent <- alive_cells$index %in% successful_parents
    #alive_cells$divided <- (alive_cells$birthdate <= time) & (alive_cells$birthdate > time - time_interval)

    #alive_cells$divided <- (alive_cells$deathdate > time) & (alive_cells$deathdate <= (time + time_interval)) & alive_cells$successful_parent
    alive_cells$divided <- (alive_cells$deathdate <= (time + time_interval)) & alive_cells$successful_parent
    divisions_of_center_cells <- sum(!alive_cells$est_edge & alive_cells$divided)
    divisions_of_edge_cells <- sum(alive_cells$est_edge & alive_cells$divided)

    #get number of deaths in time interval
    #number of deaths without offspring

    alive_cells$not_successful_parent <- !(alive_cells$index %in% successful_parents)
    not_successful_parents <-alive_cells$index[alive_cells$not_successful_parent]
    #alive_cells$died <- (alive_cells$deathdate >= time) & (alive_cells$deathdate < (time + time_interval)) & alive_cells$not_successful_parent
    alive_cells$died <- (alive_cells$deathdate <= (time + time_interval)) & alive_cells$not_successful_parent
    deaths_of_center_cells <- sum(!alive_cells$est_edge & alive_cells$died)
    deaths_of_edge_cells <- sum(alive_cells$est_edge & alive_cells$died)


    return(data.frame("time" = time,
                      "birth_rate" = c(divisions_of_edge_cells/N_edge_cells/time_interval,
                                       divisions_of_center_cells/N_center_cells/time_interval),
                      "state" = c("edge","center"),
                      "divisions" = c(divisions_of_edge_cells,divisions_of_center_cells),
                      "deaths" = c(deaths_of_edge_cells, deaths_of_center_cells),
                      "death_rate" = c(deaths_of_edge_cells / N_edge_cells/time_interval,
                                       deaths_of_center_cells/ N_center_cells/time_interval),
                      "growth_rate" = c((divisions_of_edge_cells - deaths_of_edge_cells)/N_edge_cells/(time_interval),
                                        (divisions_of_center_cells - deaths_of_center_cells)/N_center_cells/time_interval),
                      "N" = c(N_edge_cells, N_center_cells)
    )
    )


}

#' Edge distances versus cell divisions over time
#'
#' Extract instantaneous fates of cells over simulation time and record distance from current tumor edge
#'
#' @param all_cells data.frame
#' @param time numeric
#' @param time_interval numeric
#' @importFrom dplyr filter
#' @importFrom purrr map2_dbl
#' @importFrom magrittr %>%
#' @importFrom spatstat.geom nndist
#' @return numeric
#' @export
get_edge_dist_versus_divisions_at_t <- function (all_cells, time, time_interval = 2/24, cell_locations_df_file = NULL) {
    
    successful_parents <- sort(unique(all_cells$parent_index))
    alive_cells <- all_cells %>% dplyr::filter(round(birthdate, 
                                                     2) <= round(time, 2) & round(deathdate, 2) > round(time, 
                                                                                                        2))
    if (!is.null(cell_locations_df_file)) {
        cell_locations_df <- read_csv(cell_locations_df_file)
        cell_locations_df_t <- cell_locations_df %>% filter(round(t, 
                                                                  2) == round(time, 2))
        alive_cells <- alive_cells %>% dplyr::select(-locx, -locy) %>% 
            dplyr::left_join(., cell_locations_df_t, by = "index")
    }
    alive_cells$most_extreme <- purrr::map2_dbl(alive_cells$locx, 
                                                alive_cells$locy, function(x_loc, y_loc) check_if_cell_most_extreme(x_loc, 
                                                                                                                    y_loc, alive_cells = alive_cells))
    alive_cells$edge_adjacent <- purrr::map2_dbl(alive_cells$locx, 
                                                 alive_cells$locy, function(x_loc, y_loc) check_boundary_adjacent(x_loc, 
                                                                                                                  y_loc, alive_cells = alive_cells))
    alive_cells$est_edge <- alive_cells$most_extreme | alive_cells$edge_adjacent
    nearest_points <- spatstat.geom::nndist(X = alive_cells$locx, 
                                            Y = alive_cells$locy, by = alive_cells$est_edge)
    alive_cells$dist_from_edge <- (1 - alive_cells$est_edge) * 
        nearest_points[, which(colnames(nearest_points) == "TRUE")]
    alive_cells$successful_parent <- alive_cells$index %in% successful_parents
    alive_cells$divided <- (alive_cells$deathdate <= (time + 
                                                          time_interval)) & alive_cells$successful_parent
    alive_cells$not_successful_parent <- !(alive_cells$index %in% 
                                               successful_parents)
    not_successful_parents <- alive_cells$index[alive_cells$not_successful_parent]
    alive_cells$died <- (alive_cells$deathdate <= (time + time_interval)) & 
        alive_cells$not_successful_parent
    cell_fate_df <- alive_cells %>% add_column(time = time)
    return(cell_fate_df)
}

#' Get n divisions
#'
#' Get number of divisions along lineage of extant cells
#'
#' @param index numerical index of extant cell
#' @param all_cells data.frame
#' @return number of cell divisions
#' @export
get_n_divisions <- function(index, all_cells) {
    
    ancestors <- collect_all_ancestors(index, all_cells)
    
    return(length(ancestors))
}

#' Get all cell states at death
#'
#' Get states of cells as death by time slice, this method is faster than going cell by cell. 
#'
#' @param all_cells data.frame
#' @param cell_locations_df_file CSV file of cell locations. Default: null
#' @return number of cell divisions
#' @export
find_all_cell_state_at_death <- function (all_cells, cell_locations_df_file = NULL) {
    deathdates <- sort(unique(all_cells$deathdate))
    all_cells$state <- NA
    if (!is.null(cell_locations_df_file)) {
        
        cell_locations_df <- read_csv(cell_locations_df_file)
    }
    for (time in deathdates) {
        alive_cells <- all_cells %>% dplyr::filter(birthdate < 
                                                       time & deathdate >= time)
        if (!is.null(cell_locations_df_file)) {
            
            
            cell_locations_df_t <- cell_locations_df %>% dplyr::filter(round(t, 
                                                                             2) == round(time, 2))
            alive_cells <- alive_cells %>% dplyr::select(-locx, 
                                                         -locy) %>% dplyr::left_join(., cell_locations_df_t, 
                                                                                     by = "index")
        }
        dead_cells <- alive_cells %>% dplyr::filter(round(deathdate, 
                                                          2) == round(time, 2))
        dead_cells <- mark_boundary(dead_cells, alive_cells)
        all_cells$state[match(dead_cells$index, all_cells$index)] <- dead_cells$est_edge
    }
    return(all_cells)
}


