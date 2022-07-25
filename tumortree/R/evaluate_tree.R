#' Match simlated and estimated node
#'
#' Match node in simulated tree to node in estimated tree
#' based on set of child leaves descended from node.
#' 
#' @param node integer
#' @param sim_tree treedata
#' @param est_tree treedata
#' @importFrom phangorn mrca.phylo Descendants
#' 
#' @return integer
#' @export
match_sim_est_node <- function(node, sim_tree, est_tree) {
    
    rev_bayes_node <- phangorn::mrca.phylo(est_tree@phylo,match(sim_tree@phylo$tip.label[phangorn::Descendants(sim_tree@phylo, node = node , type = "tips")[[1]]], est_tree@phylo$tip.label))
    return(rev_bayes_node)
}
#' Create node conversion chart
#'
#' Make data frame that matches simulated nodes to 
#' nodes in estimated tree
#' 
#' @param sim_tree treedata
#' @param est_tree treedata
#' @importFrom purrr map_dbl
#' 
#' @return data.frame
#' @export
create_node_conversion_chart <- function(sim_tree, est_tree) {
    sim_nodes <- sim_tree@data$node[sim_tree@data$node > length(sim_tree@phylo$tip.label)]
    est_nodes <- purrr::map_dbl(sim_nodes, function(node) match_sim_est_node(node = node, sim_tree = sim_tree, est_tree = est_tree))
    sim_leaves<- sim_tree@data$node[sim_tree@data$node <= length(sim_tree@phylo$tip.label)]
    est_leaves <- match(sim_tree@phylo$tip.label, est_tree@phylo$tip.label)
    
    conversion_df <- data.frame("sim" = c(sim_leaves, sim_nodes), "est" = c(est_leaves, est_nodes))
    return(conversion_df)
}
#' Make calibration dataframe
#'
#' Compare the probability of edge state with 
#' the true simulated states for internal nodes
#' of the esimated tree.  
#' 
#' @param sim_tree treedata
#' @param est_tree treedata
#' @param type string
#' @importFrom purrr map_dbl
#' 
#' @return data.frame
#' @export
make_calbriation_df <- function(sim_tree, est_tree, type = c("revbayes", "bdmm")) {
    
    sim_tree@data <- sim_tree@data %>% 
        dplyr::arrange(as.numeric(node))
    
    est_tree@data <- est_tree@data %>% 
        dplyr::arrange(as.numeric(node))
    
    calibration_df <- create_node_conversion_chart(sim_tree, est_tree)
    calibration_df$death_edge_state <- sim_tree@data$death_edge_state
    
    if (type == "revbayes") {
        
        calibration_df$predicted_states <- as.numeric(est_tree@data[calibration_df$est,]$anc_state_1)
        calibration_df$predicted_states_post <- as.numeric(est_tree@data[calibration_df$est,]$anc_state_1_pp)
        
    } else if (type == "bdmm") {
        
        calibration_df$predicted_states <- as.numeric(est_tree@data[calibration_df$est,]$type)
        calibration_df$predicted_states_post <- as.numeric(est_tree@data[calibration_df$est,]$type.prob)
        
    } else {
        
        stop("Valid tree type must be provided")
    }
    calibration_df$correct <- as.numeric(calibration_df$death_edge_state  == calibration_df$predicted_states)
    
    # filter to only estimated nodes, as leaves are given states
    calibration_df <- calibration_df %>% 
        dplyr::filter(sim > length(sim_tree@phylo$tip.label))
    
    #unwrap probabilities
    calibration_df$est_edge_prob <- calibration_df$predicted_states * calibration_df$predicted_states_post + (1 -calibration_df$predicted_states ) * (1 - calibration_df$predicted_states_post)
    return(calibration_df)
} 
