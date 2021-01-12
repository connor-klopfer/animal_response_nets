# Functions gor general analysis of the scripts, not animal specific

# Import the analysis scripts here, where you import, run the analysis, and return
# the final dataset
source("analysis_scripts/mice_analysis.R", echo = FALSE)
source("analysis_scripts/guppy_analysis.R", echo = FALSE)

get_graph_metrics <- function(g, results_list){
  #' Get Graph metrics to be applied generally across all networks. 
  #' 
  #' Not specific to an animal network 
  #' 
  #' @param g igraph Graph Object
  #' @param results_list named list, with the columns outline in the analysis template file. 
  #' 
  #' Create a dataframe with the network metrics. This data set will have the
  #' following columns:
  #' 
  #'  treatement: 
  #'  timepoint: 
  #'  network number 
  #'  population: 
  #'  average degree: 
  #'  clustering coefficient:
  #'  eigenvalue centrality: Pernode 
  #'       may need to be
  #' @return results_list: named list with the additional network metrics
  is_g_directed <- is_directed(g)
  results_list[["density"]] <- igraph::edge_density(g)
  results_list[["eignvector_centrality"]] <- mean(igraph::eigen_centrality(g, directed = is_g_directed)$vector)
  results_list[["betweenness_centrality"]] <- mean(igraph::betweenness(g, directed = is_g_directed))
  results_list[["n_nodes"]] <- length(V(g))
  
  return(results_list)
}

run_full_analysis <- function(){
  #' Run fulll analysis on all the animal networks
  #' 
  all_analysis <- list()
  
  #' Append the function to iterate through your data structure and get the 
  #' network metrics here:
  all_analysis[['guppy']] <- guppy_net_analysis()
  all_analysis[['mice']] <- mice_net_analysis()
  
  
  final_df <- do.call(rbind, all_analysis)
  return(final_df)
}

final_df <- run_full_analysis()