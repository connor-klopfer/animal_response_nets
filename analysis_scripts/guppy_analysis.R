#' Iterate through the guppy datasets, and perform the ananlysis on all the 
#' guppy networks and append into a large dataset. 

source("data_prep_scripts/import_guppy_networks.R", echo = FALSE)

guppy_net_analysis <- function(){
  #' Iterate through the list of lists in to form analysis on the mice networks.
  
  g_df <- import_guppy_networks()
  
  master_list <- list()
  master_idx <- 1
  for(r in names(g_df)){
    for(n in 1:length(g_df[[r]])){
      temp_net <- g_df[[r]][[n]]
      
      # Metadata
       temp_list <- c(
        "network" = "guppy",
        "animal" = "guppy",
        "condition" = r,
        "replicate_num" = n
      )
       
      # NOTE: Warnings are emitted because there is a weighted directed graphs 
      # for the eignevector centrality caluculation. Might need to adress to see
      # if: acyclic, component size < 2
      # warnings for all nets
       
      temp_list <- get_graph_metrics(temp_net, temp_list)
      master_list[[master_idx]] <- temp_list
      master_idx <- master_idx + 1
    }
  }
  # Append to dataframe. 
  final <- data.frame(do.call(rbind, master_list))
  return(final)
}

# guppy_results <- guppy_net_analysis()
