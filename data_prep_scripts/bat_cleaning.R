#' Data Converting hour-list of bat networks (67 hours)into ta dataframe


write_bat_edgelist_tofile <- function(hour_list){
  # Save all df's to a file
  all_dfs <- list()
  
  # For every hour
  for(h in 1:length(hour_list)){
    # Write to edgelist
    temp_edgelist <- igraph::as_data_frame(hour_list[[h]], what = "edges")
    temp_edgelist$hour <- h
    all_dfs[[h]] <- temp_edgelist
  }
  
  final <- do.call(rbind, all_dfs)
  
  write.csv(final, file = "data/bats/bats_byhour_edgelist.R", row.names = F)

}

write_bat_edgelist_tofile(net.list)


# loop for generating multiple networks
# subset by hour, create an igraph for each hour, calculate metrics for each network hour

g <- read.table("bats_byhour_edgelist.csv", header = TRUE, sep = ",")

keys<-unique(g$hour)

new_df<-do.call(bind_rows, lapply(1:length(keys), function(i){
  print(paste0(i," of ",length(keys)))
  #subset out each hour
  temp<- g %>% 
    filter(hour==keys[i])
  #make igraph object
  g2<-graph.data.frame(temp) #this is the igraph object
  #compute metrics (e.g., density)
  g3<- length(V(g2))
  g4<-igraph::edge_density(g2)
  g5<-mean(igraph::eigen_centrality(g2, directed = FALSE)$vector)
  g6<-mean(igraph::betweenness(g2, directed = FALSE))
  g7<-igraph::mean_distance(g2)
  g8 <- effective_information(g2, effectiveness = FALSE)
  g9 <- g8/log2(g3) # normalized value to control for network size 
  #write igraph object to file so that we can plot network later for each hour
  save(g2, file=paste0("/Users/laurenstanton/Desktop/CNWW/Bats/igraph_objects_4_plotting/","hour_",keys[i],".RData"))
  #create summary dataframe
  dat<-data.frame(hour=keys[i], n_nodes =g3, density=g4, mean_eignvector_centrality=g5, mean_betweenness_centrality=g6, average.path.length = g7, effective_information=g8,ie_adjusted=g9)
return(dat)
}))

write.table(new_df,file="/Users/laurenstanton/Desktop/CNWW/Bats/bat_metrics.csv",sep=",", row.names=FALSE)

#to plot all 67 networks
for(i in 1:length(keys)){
  load(paste0("/Users/laurenstanton/Desktop/CNWW/Bats/igraph_objects_4_plotting/","hour_",keys[i],".RData"))
  plot(g2) #need to figure out what parameters to use for plotting networks - LAS
}

#to plot one at a time by changing value of i
i=15
load(paste0("/Users/laurenstanton/Desktop/CNWW/Bats/igraph_objects_4_plotting/","hour_",keys[i],".RData"))
plot(g2)  # need to figure out what parameters to use for plotting networks - LAS
