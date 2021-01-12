#' Analysis Scripts for comparing the change in network metric for analysing the
#' mouse nets to compare across different species. 
#' 
#' The metrics that are analysed will be used: 
#' 

# Import the mice networks. 
source("data_prep_scripts/import_mice_networks.R", echo = FALSE)

mice_net_analysis <- function(){
  #' Iterate through the list of lists in to form analysis on the mice networks. 
  m_df <- import_mice_networks()
  
  master_list <- list()
  master_idx <- 1
  for(tr in names(m_df)){
    for(tp in names(m_df[[tr]])){
      for(n in 1:length(m_df[[tr]][[tp]])){
        temp_net <- m_df[[tr]][[tp]][[n]]

        # Metadata
        temp_list <- c(
          "network" = paste0("mouse_", tr),
          "animal" = "mouse",
          "condition" = tp,
          "replicate_num" = n
          )
        temp_list <- get_graph_metrics(temp_net, temp_list)
        master_list[[master_idx]] <- temp_list
        master_idx <- master_idx + 1
      }
    }
  }
  # Append to dataframe. 
  final <- data.frame(do.call(rbind, master_list))
  return(final)
}