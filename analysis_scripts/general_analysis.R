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
  if(max(all_betweenness) == 0){
    normalised_betweenness <- (all_betweenness - min(all_betweenness)) / 
      (1)
  }else{
    normalised_betweenness <- (all_betweenness - min(all_betweenness)) / 
                                 (max(all_betweenness) - min(all_betweenness))
    
  }
  # normalised_betweenness <- all_betweenness
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
  
  all_individual_analysis[['guppy']] <- guppy_ind_analysis()# %>% mutate_all(as.character)
  all_population_analysis[['guppy']] <- guppy_pop_analysis()# %>% mutate_all(as.character)
  
  all_individual_analysis[['mice']] <- mice_ind_analysis()# %>% mutate_all(as.character)
  all_population_analysis[["mice"]] <- mice_pop_analysis()# %>% mutate_all(as.character)
  
  all_individual_analysis[["bats"]] <- import_analysis_results("bats_node_metrics.csv") %>% 
    select(-`Network ID`, -exp_treatment) %>% group_by(condition, replicate_num) %>%
    mutate(betweeness = as.numeric(betweeness)) %>% 
    mutate(betweeness = (betweeness - min(betweeness))/(max(betweeness) - min(betweeness))) %>% 
    ungroup() %>% mutate(betweeness = as.character(betweeness))
  
  
  all_population_analysis[["bats"]] <- import_analysis_results("bats_network_metrics.csv")
  
  all_population_analysis[['tree_lizards']] <- import_analysis_results("network_output_treelizards_0weightexcl.csv") %>%
    select(-X1)
  all_individual_analysis[['tree_lizards']] <- import_analysis_results("node_output_treelizards_0weightexcl.csv") %>% 
    select(-X1)
  
  all_population_analysis[['wood_ants']] <- import_analysis_results("network_output_woodants.csv") %>%
    select(-X1)
  all_individual_analysis[['wood_ants']] <- import_analysis_results("node_output_woodants.csv") %>% 
    select(-X1)
  
  all_population_analysis[['stickleback']] <- import_analysis_results("stickleback_networkdata.csv") %>%
    select(-X1)
  all_individual_analysis[['stickleback']] <- import_analysis_results("stickleback_individualdata.csv") %>% 
    select(-X1, -NetworkID)
  
  all_population_analysis[['great_tits']] <- import_analysis_results("tits_2015_network_metrics.csv") %>%
    select(-X1)
  all_individual_analysis[['great_tits']] <- import_analysis_results("tits_2015_individual_metrics.csv") %>%
    select(-X1, -network_id)
  
  all_population_analysis[['fairy_wrens']] <- import_analysis_results("fairywren_2013_2014_pre_post_population.csv") 
  all_individual_analysis[['fairy_wrens']] <- import_analysis_results("fairywren_2013_2014_pre_post_individual.csv")
  
  
  View(all_individual_analysis)
  View(all_population_analysis)

  # print(names(append_beetle_ind_analysis()))
  do.call(rbind, all_individual_analysis) %>% 
    bind_rows(append_beetle_ind_analysis() %>% mutate_all(as.character)) %>%
    write.csv(file.path('analysis_results', "all_individual_results.csv"),
              row.names = F)

  # print(names(append_beetle_population_analysis()))
  do.call(rbind, all_population_analysis) %>% 
    rbind(append_beetle_population_analysis()) %>%
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
  
  results_filename <- c(
    
  )

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
    'network' = c("Study", "study", "Network/Study"), 
    'animal' = c("Animal", "species"), 
    'condition' = c('Treatment', "treatment", "Condition"), 
    'replicate_num' = c("GroupID", "replicate", "Replicate #", "replicate_number", "Year"), 
    'n_nodes' = c("N_Nodes", "n_nodes", "N Nodes"),
    'density' = c("Density"), 
    'eigenvector_centrality' = c('Mean_Eigenvector_Centrality',  "Eignvector Centrality", 
                                 "EigenvectorCentrality", "Eigenvector_Centrality"), 
    'betweeness' = c("Mean_Betweenness_Centrality", "Betweenness Centrality", "betweenness", "Betweenness", 
                                 "betweeness", "betweenness_centrality", "betweenness_centrality", 
                     "BetweennessCentrality"),
    "effective_information" = c("Effective_Information", "normal_ei", "Effective Information"),
    "degree" = c("Degree"),
    "node_ID" = c("Node ID", "nodeID", "ID", "NodeID", "node_id")
    
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

import_analysis_results <- function(filename){
  
  imported_data <- readr::read_csv(file.path("analysis_results", filename)) %>% 
    normalise_names()
  
  if(filename == "bats_node_metrics.csv"){
    names(imported_data)[length(names(imported_data))] <- "exp_treatment"
    # imported_data <- mutate(betweenness = betweenness_centrality) %>% select(-betweenness_centrality)
  }
  return(imported_data %>% mutate_all(as.character))
}

run_full_analysis()
