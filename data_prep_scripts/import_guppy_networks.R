#' Import guppy networks in R 
#' 
#' Imports from the clean association networks. Attaches some of the fish 
#'  attributes to the nodes. The format will be a list, with one segment thAT
#'   goes to a list of high risk networks and low risk groups, where the matching 
#'   experiement groups are matched by index in the list. 

import_guppy_networks <- function(){
  
  filepath <- file.path("data", "guppies", "reformatted_guppy_networks.csv")
  all_nets <- list()
  guppy_edges <- readr::read_csv(filepath)
  
  for(r in unique(guppy_edges$risk)){
    all_nets[[r]] <- vector(mode = "list", length = 16)
    for(g in unique(guppy_edges$group)){
      all_nets[[r]][[g]] <- guppy_edges %>% 
        filter(risk == r, group == g, weight > 0) %>% 
        graph_from_data_frame()
    }
  }
  return(all_nets)

}