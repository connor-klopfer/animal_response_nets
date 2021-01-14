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
  results_list[["betweenness_centrality"]] <- get_normalised_betweenness(g, is_g_directed)
  results_list[["n_nodes"]] <- length(V(g))
  results_list[['effective_information']] <- einet::effective_information(g, effectiveness = F)
  # May need to normalise, if so, 
  results_list[['effective_information']] <- as.numeric(results_list[['effective_information']])/
    log2(as.numeric(length(V(g))))
  
  return(results_list)
}

get_graph_ind_metrics <- function(g){
  is_g_directed <- is_directed(g)
  
  ind_level_df <- convert_metric_to_df(
    igraph::eigen_centrality(g, directed = is_g_directed)$vector, "eigenvector_centrality") %>%
    left_join(
      get_normalised_betweenness(g, is_g_directed) %>%
        convert_metric_to_df("betweeness"), by = "node_ID") %>% 
    left_join(igraph::degree(g) %>% convert_metric_to_df("degree"), 
              by = "node_ID")
  
  return(ind_level_df)
}


convert_metric_to_df <- function(metric, metric_name){
  final <- data.frame(names(metric), metric, row.names = NULL)
  names(final)<- c("node_ID", metric_name)
  return(final)
}


get_graph_pop_metrics <- function(g, results_list){
  is_g_directed <- is_directed(g)
  results_list[["density"]] <- igraph::edge_density(g)
  results_list[["n_nodes"]] <- length(V(g))
  results_list[['effective_information']] <- einet::effective_information(g, effectiveness = F)
  # May need to normalise, if so, 
  results_list[['effective_information']] <- as.numeric(results_list[['effective_information']])/
    log2(as.numeric(length(V(g))))
  return(results_list)
}


get_normalised_betweenness <- function(g, is_directed){
  all_betweenness <- igraph::betweenness(g, directed = is_directed)
  normalised_betweenness <- (all_betweenness - min(all_betweenness)) / 
                               (max(all_betweenness) - min(all_betweenness))
  return(normalised_betweenness)
}

run_full_analysis <- function(){
  #' Run fulll analysis on all the animal networks
  #' 
  all_individual_analysis <- list()
  all_population_analysis <- list()
  #' Append the function to iterate through your data structure and get the 
  #' network metrics here:
  # all_analysis[['guppy']] <- guppy_net_analysis()
  # all_analysis[['mice']] <- mice_net_analysis()
  
  all_individual_analysis[['guppy']] <- guppy_ind_analysis()
  all_population_analysis[['guppy']] <- guppy_pop_analysis()
  
  all_individual_analysis[['mice']] <- mice_ind_analysis()
  all_population_analysis[["mice"]] <- mice_pop_analysis()
  
  do.call(rbind, all_individual_analysis) %>% 
    write.csv(file.path('analysis_results', "all_individual_results.csv"), 
              row.names = F)

  do.call(rbind, all_population_analysis) %>% 
    write.csv(file.path('analysis_results', "all_population_results.csv"), 
              row.names = F)
  
  # final_df <- do.call(rbind, all_analysis) %>% 
  #   rbind(get_processed_results()) %>% bind_rows(append_beetle_analysis())
  # return(final_df)
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

append_beetle_analysis <- function(){
  beetles <- readr::read_csv(file.path("analysis_results", "mean_measures_beetles.csv"))
  names(beetles) <- c("replicate_num", "condition", "n_nodes", "eigenvector_centrality", 
                      "betweenness_centrality")
  beetles$network <- "guppy"
  beetles$animal <- "beetle"
  return(beetles)
}

append_beetle_analysis()

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
