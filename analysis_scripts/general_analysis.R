# Functions gor general analysis of the scripts, not animal specific

get_graph_metrics <- function(g, results_list){
  #' Get Graph metrics to be applied generally across all networks. 
  #' 
  #' Not specific to an animal network 
  #' 
  #' @param g igraph Graph Object
  #' @param results_list named list, with the columns outline in the analysis template file. 
  #'
  #'@return results_list: named list with the additional network metrics
  results_list[["density"]] <- igraph::edge_density(g)
  results_list[["eignvector_centrality"]] <- mean(igraph::eigen_centrality(g, directed = T)$vector)
  results_list[["betweenness_centrality"]] <- mean(igraph::betweenness(g, directed = T))
  results_list[["n_nodes"]] <- length(V(g))
  
  return(results_list)
}