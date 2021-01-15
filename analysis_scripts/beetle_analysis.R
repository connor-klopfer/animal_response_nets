#' Append the beetle analysis for both the population and individual level metrics. 
#' Run the effective information criterion. 


append_beetle_population_analysis <- function(){
  beetles <- readr::read_csv(file.path("analysis_results", "Population_metrics_beetle.csv"))
  names(beetles) <- c("replicate_num", "condition", "n_nodes", "density")
  beetles$network <- "beetle"
  beetles$animal <- "beetle"
  
  beetle_edges <- readr::read_csv(file.path("data", 
                                            "Beetles_FormicaEtAl2017","SRI_edgelist_no_control.csv")) %>% 
    select(-X1)

  beetle_nets <- list()
  for(treatment in unique(beetle_edges$Group.ID.Populations)){
    beetle_nets[[treatment]] <- list()
    for(tp in unique(beetle_edges$Treatment.Period)){
      data_subset <- beetle_edges %>% 
        filter(Group.ID.Populations == treatment, Treatment.Period == tp)
      
      beetle_nets[[treatment]][[tp]] <- igraph::graph_from_data_frame(data_subset)
    }
  }
  
  beetles$effective_information <- NA
  for(c in unique(beetles$replicate_num)){
    if(c %in% c("C1", "C2", "C3", "C4")){
      next()
    }
    for(tp in unique(beetles$condition)){
      beetles$effective_information[which((beetles$replicate_num == c) & 
                                            beetles$condition == tp)] <- as.character(
                                              einet::effective_information(beetle_nets[[c]][[tp]]))
    }
  }
  return(beetles %>% 
           mutate(condition = as.character(condition), 
                  n_nodes = as.character(n_nodes), 
                  animal = "beetle", 
                  network = "beetle") %>% 
           filter(!(replicate_num %in% c("C1", "C2", "C3", "C4"))))
}

append_beetle_population_analysis()


append_beetle_ind_analysis <- function(){
  beetle <- readr::read_csv(file.path("analysis_results", "Individual_metrics_beetles.csv")) %>% 
    select(-CLOSENESS_CENTRALITY, -DEGREE_NORMALIZED)
  names(beetle) <- c("replicate_num", "condition", "node_ID", 'degree', "eigenvector_centrality", "betweeness")
  return(beetle %>% mutate(network = "beetle", animal = "beetle") %>% 
           filter(!(replicate_num %in% c("C1", "C2", "C3", "C4"))))
  }



View(append_beetle_ind_analysis())


get_beetle_nets <- function(){
  beetle_edges <- readr::read_csv(file.path("data", 
                            "Beetles_FormicaEtAl2017","SRI_edgelist_no_control.csv")) %>% select(-X1)
  # View(beetle_edges)
  
  all_beetles <- list()
  for(treatment in unique(beetle_edges$Group.ID.Populations)){
    all_beetles[[treatment]] <- list()
    for(tp in unique(beetle_edges$Treatment.Period)){
      data_subset <- beetle_edges %>% 
        filter(Group.ID.Populations == treatment, Treatment.Period == tp)
      
      all_beetles[[treatment]][[tp]] <- 2
      #igraph::graph_from_data_frame(data_subset)
    }
  }
  
  return(all_beetles)
}

# beetle_sample <- read_in_beetle_graphs()


# append_beetle_analysis()
