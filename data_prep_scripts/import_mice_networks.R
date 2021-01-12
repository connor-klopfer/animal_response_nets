#' Import Mice Netorks from the partitioned dataset. 
#' 
#' By Connor Klopfer 10Jan21
#' 
#' Imports network into a list of igraph networks. The structure is similar to 
#' the structurr outlined in 'import_mice_networks.py' but with a list in place 
#' if a python dictionary.  

import_mice_networks <- function(){
  raw_mice_data <- readr::read_csv(
    file.path("data", 'mice', "partitioned_dryad_encounters.csv"))
  
  # Flip the columns to be igraph compatible
  all_mice_data <- raw_mice_data %>% 
    relocate(`IDA[inside_first]`, .before = box) %>% 
    relocate(`IDB[inside_after_IDA]`, .before = box)
  
  # Create a master list to hold all the networks 
  all_nets = list()
  
  # For all treatment levels (LPS, Control, on_injection)
  for(tr in unique(all_mice_data$treatment)){
    all_nets[[tr]] = list()
    
    for(tp in unique(all_mice_data$timepoint)){
      # Dataset subset 
      condition_subset <- all_mice_data %>% 
        filter(treatment == tr, timepoint == tp)
      
      # Inititalise list
      all_nets[[tr]][[tp]] = vector(
        mode = "list",
        length = length(unique(condition_subset$component_num)))
      
      # For all the components within that condition
      for(c in unique(condition_subset$component_num)){
        
        # Make the graph of the component from the dataframe subset
        all_nets[[tr]][[tp]][[c + 1]] <- condition_subset %>% 
          filter(component_num == c) %>% 
          select(-c("treatment", "timepoint", "component_num")) %>%
          graph_from_data_frame()
      }
    }
  }
  
  return(all_nets)
}