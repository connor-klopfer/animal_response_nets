# Functions gor general analysis of the scripts, not animal specific

# Import the analysis scripts here, where you import, run the analysis, and return
# the final dataset
library(dplyr)
library(igraph)
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
  results_list[["eigenvector_centrality"]] <- mean(igraph::eigen_centrality(g, directed = is_g_directed)$vector)
  results_list[["betweenness_centrality"]] <- mean(igraph::betweenness(g, directed = is_g_directed))
  results_list[["n_nodes"]] <- length(V(g))
  results_list[['effective_information']] <- einet::effective_information(g, effectiveness = F)
  # May need to normalise, if so, 
  results_list[['effective_information']] <- as.numeric(results_list[['effective_information']])/
    log2(as.numeric(length(V(g))))
  
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
  
  
  
  
  final_df <- do.call(rbind, all_analysis) %>% 
    rbind(get_processed_results())
  return(final_df)
}

get_processed_results <- function(){
  #' Add in processed results uploaded to gihub seperately and append to a final 
  #' dataset 
  
  all_result_dfs <- list()
  parent_dir <- "data_prep_scripts"
  results_filename <- list.files(path = parent_dir, pattern = ".csv")
  
  # FOr every CSV file in analysis scripts/
  for(r in 1:length(results_filename)){
    # Filepath
    results <- file.path(parent_dir, results_filename[r])
    results_df <- readr::read_csv(results) %>% select(-X1)
    # Convert the names so they're all the same. 
    results_df <- normalise_names(results_df)
    all_result_dfs[[r]] <- results_df
  }
  return(do.call(rbind, all_result_dfs))
}

normalise_names <- function(d_f){
  #' Normalise names for a given dataset. 
  names_dict <- list(
    'network' = c("Study", "study"), 
    'animal' = c("Animal", "species"), 
    'condition' = c('Treatment', "treatment"), 
    'replicate_num' = c("GroupID", "replicate"), 
    'n_nodes' = c("N_Nodes", "n_nodes"),
    'density' = c("Density"), 
    'eigenvector_centrality' = c('Mean_Eigenvector_Centrality'), 
    'betweenness_centrality' = c("Mean_Betweenness_Centrality"),
    "effective_information" = c("Effective_Information", "normal_ei")
  )
  
  for(n in names(names_dict)){
    for(nam in names(d_f)){
      if(nam %in% names_dict[[n]]){
        names(d_f)[which(names(d_f) == nam)] <- n
      }
    }
  }
  return(d_f)
}

final_df <- run_full_analysis()
