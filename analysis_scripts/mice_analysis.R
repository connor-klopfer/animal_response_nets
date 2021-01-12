#' Analysis Scripts for comparing the change in network metric for analysing the
#' mouse nets to compare across different species. 
#' 
#' The metrics that are analysed will be used: 
#' 

# Import the mice networks. 
source("data_prep_scripts/import_mice_networks.R", echo = FALSE)
source("analysis_scripts/general_analysis.R", echo = FALSE)

mice_net_analysis <- function(m_df){
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
  #'      
  master_list <- list()
  master_idx <- 1
  for(tr in names(m_df)){
    for(tp in names(m_df[[tr]])){
      for(n in 1:length(m_df[[tr]][[tp]])){
        temp_net <- m_df[[tr]][[tp]][[n]]
        
        # NOTE: Temp fix for analysis, need to see whay these are being imported 
        # as null. 
        if(is.null(temp_net)){
          next()
        }
        # Metadata
        master_list[[master_idx]] <- c(
          "network" = paste0("mouse_", tr),
          "animal" = "mouse",
          "condition" = tp,
          "replicate_num" = n
          )
        master_idx <- master_idx + 1
      }
    }
  }
  # Append to dataframe. 
  final <- data.frame(do.call(rbind, master_list))
  return(final)
}

# Dataframe of the mice analysis results. 
results_df <- get_mice_metrics(all_mice_nets)
