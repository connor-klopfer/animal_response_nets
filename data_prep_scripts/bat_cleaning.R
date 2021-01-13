#' Data Converting hour-list of bat networks (67 hours)into ta dataframe


write_bat_edgelist_tofile <- function(hour_list){
  # Save all df's to a file
  all_dfs <- list()
  
  # For every hour
  for(h in 1:length(hour_list)){
    # Write to edgelist
    temp_edgelist <- igraph::as_data_frame(hour_list[[h]])
    temp_edgelist$hour <- h
    all_dfs[[h]] <- temp_edgelist
  }
  
  final <- do.call(rbind, all_dfs)
  
  write.csv(final, file = "data/bats/bats_byhour_edgelist.R", row.names = F)

}

write_bat_edgelist_tofile(net.list)
