#' Cleaning the grate_tits data into a clean edgelist from the simplical 
#' #' edgelist version. 

reformat_fairywrens_dataset <- function(original){
  #' Should give all the permutations seen in the dataset into a long, edgelist format. 
  bird_data <- original %>% select(contains("Bird"))
  
  id_data <- original %>% select(-contains("Bird"))
  bird_combos <- combn(bird_data, 2, simplify = FALSE)
  
  all_dfs <- list()
  for(x in 1:length(bird_combos)){
    temp_dataset <- cbind(id_data, bird_combos[[x]])
    n_cols <- ncol(temp_dataset)
    names(temp_dataset)[(n_cols - 1):n_cols] <- c("Bird1", "Bird2")
    all_dfs[[x]] <- temp_dataset
  }
  
  final <- do.call(rbind, all_dfs) %>% filter(!is.na(Bird1), !is.na(Bird2))
  
  write.csv(final, "data/Fairywrens_LantzAndKarubian2017/fairwren_edgelist.csv", row.names = F)
  return(final)
}

reformatted_fairywrens <- reformat_fairywrens_dataset(S1Table)
