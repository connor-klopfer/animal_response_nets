rm(list=ls())

# Great tits #############

#setwd("C:/Users/alexis/Google Drive (1)/Between Computers/CNWW2021/Great_tits_Firth et al 2017")
load("C:/Users/alexis/Google Drive (1)/Between Computers/CNWW2021/Great_tits_Firth et al 2017/Raw Data For Firth et al.RData")

library(asnipe)
library(dplyr)
library(igraph)
library(einet)

# time must be numeric for asnipe package
records$TimeStamp<-as.numeric(as.POSIXct(records$time)) 

for(i in unique(records$week)){
  print(i)
  df1<-subset(records,week==i)
  
  df_assoc_pts<-(df1 %>%
                   select(ring, TimeStamp,location))
  
  try(Location2 <- df_assoc_pts %>% count(location) %>% filter(n > 2))
  
  try(df_assoc_pts_filtered<- df_assoc_pts %>% filter(location %in% Location2$location)) #all had more than 2 points anyway so wasn't really necessary
  
  try(df_assoc_pts_filtered$ring<-as.factor(df_assoc_pts_filtered$ring))
  
  try(global_ids<-levels(df_assoc_pts_filtered$ring))
  
  try(gmm_data<- asnipe::gmmevents(time = df_assoc_pts_filtered$TimeStamp,
                                   identity = df_assoc_pts_filtered$ring,
                                   location = df_assoc_pts_filtered$location,
                                   global_ids=global_ids))
  
  
  
  try(gbi<-(gmm_data$gbi))
  try(events<- gmm_data$metadata %>% mutate(time_window = (End - Start)))
  try(observations_per_event <- gmm_data$B)
  
  try(adjm_gmm <- asnipe::get_network(gbi, data_format = "GBI", association_index = "SRI"))
  try(net_gmm <- igraph::graph.adjacency(adjm_gmm, mode = "undirected", diag = F, weighted = T))
  
  try(save(adjm_gmm, file = paste("adjm_gmm_",i,".RData")))
  try( write.csv(adjm_gmm, paste("adjm_gmm_",i,".csv")))
  try(write.csv(gbi, paste("gbi_gmm_",i,".csv")))
  
  
  try( png(paste("plot", i, ".png", sep = ""), width = 1000, height = 1000))
  
  try(plot(net_gmm))
  
  dev.off()
  

  try(between<-igraph::betweenness(net_gmm, directed = TRUE,weights=E(net_gmm)$weight))
  
  eigen <- eigen_centrality(net_gmm, directed = TRUE, weights = E(net_gmm)$weight)
  ec <- eigen$vector
  ec.value <- eigen$value
  
  
 # pool individual
  pool_id <- cbind.data.frame(between, 
                              eigen, ec, ec.value)
  
  ### group-level ##
  
  density<-graph.density(net_gmm) #should be same as:  density=mean(igraph::degree(net_gmm))/(vcount(net_gmm)-1)
  
  apl <- igraph::mean_distance(net_gmm) #average.path.length(graph.ref.crowd) 
  
  ei <- effective_information(net_gmm, effectiveness = FALSE)
  
  n<-length(unique(df1$ring))
  
  eff <- ei/log2(n) # normalized value to control for network size
  
  
  #Find proportion unknown relationships, a measure of sparseness
  prunk <- EloRating::prunk(adjm_gmm)
  
  prunk.pu <- as.numeric(prunk[1])
  
  prunk.dyads <- as.numeric(prunk[2])
  
  # pool group
  pool_group <- cbind.data.frame(ec, ec.value,
                                 apl,
                                 ei, eff, 
                                 prunk.pu, prunk.dyads, density
  )
  
  #try(pool_id$ring <- rownames(pool_id))
  
  try(pool_id <-cbind(rownames(pool_id), data.frame(pool_id, row.names=NULL))) 
  
  names(pool_id)[names(pool_id) == "rownames(pool_id)"] <- "ring"
  
  try(pool_id$week<-i)
  
  
  #try(pool_group$ring <- rownames(pool_group))
  
  try(pool_group <-cbind(rownames(pool_group), data.frame(pool_group, row.names=NULL))) 
  
  names(pool_group)[names(pool_group) == "rownames(pool_group)"] <- "ring"
  
  try(pool_group$week<-i)
  
  print(i)
  
  try(EdgeList<-as_edgelist(net_gmm, names = TRUE))
  try(EdgeList<-as.data.frame(EdgeList))
  try(EdgeList$week<-i)
  
  
  
  assign( paste("pool_id", i, sep = "_") , pool_id, envir = globalenv() )      
  assign( paste("pool_group", i, sep = "_") , pool_group, envir = globalenv() )      
  assign( paste("EdgeList", i, sep = "_") , EdgeList, envir = globalenv() )      
  
  
}

ID_sheets<-grep("pool_id",names(.GlobalEnv),value=TRUE)
ID_sheets_list<-do.call("list",mget(ID_sheets))
poolID_allWeeks<-do.call(rbind, ID_sheets_list)
#write.csv(poolID_allWeeks,"great_tits_Firth.et.al.2017_pool_id.csv")


group_sheets<-grep("pool_group",names(.GlobalEnv),value=TRUE)
group_sheets_list<-do.call("list",mget(group_sheets))
pool_group_allWeeks<-do.call(rbind, group_sheets_list)
#write.csv(pool_group_allWeeks,"great_tits_Firth.et.al.2017_pool_group.csv")


Edge_sheets<-grep("EdgeList",names(.GlobalEnv),value=TRUE)
Edge_sheets_list<-do.call("list",mget(Edge_sheets))
Edge_sheets_allWeeks<-do.call(rbind, Edge_sheets_list)
write.csv(Edge_sheets_allWeeks,"EdgeLists_perWeek_Great_tits_Firth et al 2017.csv")

