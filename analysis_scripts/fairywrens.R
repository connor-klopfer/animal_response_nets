fairywren<-read.csv("C:/Users/alexi/Google Drive (1)/Between Computers/CNWW2021/animal_response_nets/data/Fairywrens_LantzAndKarubian2017/fairwren_edgelist.csv")

head(fairywren)

#SP = pre or post fire, pre/post = 2013 data, PRE/POST = 2014 data
#location used to determine whether affected by fire or not

fairywren_pre<-subset(fairywren, SP=="Pre"|SP=="PRE")
fairywren_post<-subset(fairywren, SP=="Post"|SP=="POST")

fairywren_pre_2013_affected<-subset(fairywren_pre,Fire=="Y")
fairywren_post_2013_affected<-subset(fairywren_post,Fire=="Y")

fairywren_PRE_2014<-subset(fairywren_pre,SP=="PRE")
fairywren_POST_2014<-subset(fairywren_post,SP=="POST")

library(dplyr)
PRE_2014_edge<-fairywren_PRE_2014 %>% select(Bird1,Bird2)
network2014_pre <- igraph::graph_from_data_frame(PRE_2014_edge, directed=FALSE, vertices=NULL)

POST_2014_edge<-fairywren_POST_2014 %>% select(Bird1,Bird2)
network2014_post <- igraph::graph_from_data_frame(POST_2014_edge, directed=FALSE, vertices=NULL)

#Individual level
between_PRE_2014<-igraph::betweenness(network2014_pre, directed = FALSE,weights=NULL)
between_POST_2014<-igraph::betweenness(network2014_post, directed = FALSE,weights=NULL)

EC_PRE_2014<-igraph::eigen_centrality(network2014_pre, directed = FALSE,weights=NULL)$vector
EC_POST_2014<-igraph::eigen_centrality(network2014_post, directed = FALSE,weights=NULL)$vector

degree_PRE_2014<-igraph::degree(network2014_pre, mode="all")
degree_POST_2014<-igraph::degree(network2014_post, mode="all")

#########combine ##########
library(data.table)


between_PRE_2014_df<-as.data.frame(between_PRE_2014)
between_PRE_2014_df$Condition<-"Pre"
between_PRE_2014_df$Year<-"2014"
names(between_PRE_2014_df)[names(between_PRE_2014_df) == "between_PRE_2014"] <- "Betweenness"
between_PRE_2014_df<-setDT(between_PRE_2014_df, keep.rownames = "Network_ID")[]
between_PRE_2014_df$ID<-paste(between_PRE_2014_df$Network_ID,between_PRE_2014_df$Condition, sep="_")


between_POST_2014_df<-as.data.frame(between_POST_2014)
between_POST_2014_df$Condition<-"Post"
between_POST_2014_df$Year<-"2014"
names(between_POST_2014_df)[names(between_POST_2014_df) == "between_POST_2014"] <- "Betweenness"
between_POST_2014_df<-setDT(between_POST_2014_df, keep.rownames = "Network_ID")[]
between_POST_2014_df$ID<-paste(between_POST_2014_df$Network_ID,between_POST_2014_df$Condition, sep="_")


EC_PRE_2014_df<-as.data.frame(EC_PRE_2014)
EC_PRE_2014_df$Condition<-"Pre"
EC_PRE_2014_df$Year<-"2014"
names(EC_PRE_2014_df)[names(EC_PRE_2014_df) == "EC_PRE_2014"] <- "Eigenvector_Centrality"
EC_PRE_2014_df<-setDT(EC_PRE_2014_df, keep.rownames = "Network_ID")[]
EC_PRE_2014_df$ID<-paste(EC_PRE_2014_df$Network_ID,EC_PRE_2014_df$Condition, sep="_")

EC_POST_2014_df<-as.data.frame(EC_POST_2014)
EC_POST_2014_df$Condition<-"Post"
EC_POST_2014_df$Year<-"2014"
names(EC_POST_2014_df)[names(EC_POST_2014_df) == "EC_POST_2014"] <- "Eigenvector_Centrality"
EC_POST_2014_df<-setDT(EC_POST_2014_df, keep.rownames = "Network_ID")[]
EC_POST_2014_df$ID<-paste(EC_POST_2014_df$Network_ID,EC_POST_2014_df$Condition, sep="_")


degree_PRE_2014_df<-as.data.frame(degree_PRE_2014)
degree_PRE_2014_df$Condition<-"Pre"
degree_PRE_2014_df$Year<-"2014"
names(degree_PRE_2014_df)[names(degree_PRE_2014_df) == "degree_PRE_2014"] <- "Degree"
degree_PRE_2014_df<-setDT(degree_PRE_2014_df, keep.rownames = "Network_ID")[]
degree_PRE_2014_df$ID<-paste(degree_PRE_2014_df$Network_ID,degree_PRE_2014_df$Condition, sep="_")


degree_POST_2014_df<-as.data.frame(degree_POST_2014)
degree_POST_2014_df$Condition<-"Post"
degree_POST_2014_df$Year<-"2014"
names(degree_POST_2014_df)[names(degree_POST_2014_df) == "degree_POST_2014"] <- "Degree"
degree_POST_2014_df<-setDT(degree_POST_2014_df, keep.rownames = "Network_ID")[]
degree_POST_2014_df$ID<-paste(degree_POST_2014_df$Network_ID,degree_POST_2014_df$Condition, sep="_")


degree_2014_PrePost<-rbind(degree_PRE_2014_df,degree_POST_2014_df)
between_2014_PrePost<-rbind(between_PRE_2014_df,between_POST_2014_df)
EC_2014_PrePost<-rbind(EC_PRE_2014_df,EC_POST_2014_df)


deg_betwn_2014<-merge(degree_2014_PrePost,between_2014_PrePost, by="ID")
deg_betwn_EC_2014<-merge(EC_2014_PrePost,deg_betwn_2014, by="ID")
write.csv(deg_betwn_EC_2014, "C:/Users/alexi/Google Drive (1)/Between Computers/CNWW2021/animal_response_nets/data/Fairywrens_LantzAndKarubian2017/deg_betwn_EC_2014_pre_post_fairywren.csv")

#individuals<-rbind()

#Group level


n_nodes_PRE_2014<-length(igraph::V(network2014_pre))
n_nodes_POST_2014<-length(igraph::V(network2014_post))

density_PRE_2014<-igraph::edge_density(network2014_pre)
density_POST_2014<-igraph::edge_density(network2014_post)

effective_info_PRE_2014<-einet::effective_information(network2014_pre, effectiveness = F)
EI_normal_pre<-effective_info_PRE_2014/n_nodes_PRE_2014
effective_info_POST_2014<-einet::effective_information(network2014_post, effectiveness = F)
EI_normal_post<-effective_info_POST_2014/n_nodes_POST_2014
####### 2013 #############################


fairywren_pre<-subset(fairywren, SP=="Pre"|SP=="PRE")
fairywren_post<-subset(fairywren, SP=="Post"|SP=="POST")

fairywren_pre_2013_affected<-subset(fairywren_pre,Fire=="Y")
fairywren_post_2013_affected<-subset(fairywren_post,Fire=="Y")

fairywren_PRE_2013<-fairywren_pre_2013_affected
fairywren_POST_2013<-fairywren_post_2013_affected

library(dplyr)
PRE_2013_edge<-fairywren_PRE_2013 %>% select(Bird1,Bird2)
network2013_pre <- igraph::graph_from_data_frame(PRE_2013_edge, directed=FALSE, vertices=NULL)

POST_2013_edge<-fairywren_POST_2013 %>% select(Bird1,Bird2)
network2013_post <- igraph::graph_from_data_frame(POST_2013_edge, directed=FALSE, vertices=NULL)

#Individual level
between_PRE_2013<-igraph::betweenness(network2013_pre, directed = FALSE,weights=NULL)
between_POST_2013<-igraph::betweenness(network2013_post, directed = FALSE,weights=NULL)

EC_PRE_2013<-igraph::eigen_centrality(network2013_pre, directed = FALSE,weights=NULL)$vector
EC_POST_2013<-igraph::eigen_centrality(network2013_post, directed = FALSE,weights=NULL)$vector

degree_PRE_2013<-igraph::degree(network2013_pre, mode="all")
degree_POST_2013<-igraph::degree(network2013_post, mode="all")

#########combine ##########
library(data.table)


between_PRE_2013_df<-as.data.frame(between_PRE_2013)
between_PRE_2013_df$Condition<-"Pre"
between_PRE_2013_df$Year<-"2013"
names(between_PRE_2013_df)[names(between_PRE_2013_df) == "between_PRE_2013"] <- "Betweenness"
between_PRE_2013_df<-setDT(between_PRE_2013_df, keep.rownames = "Network_ID")[]
between_PRE_2013_df$ID<-paste(between_PRE_2013_df$Network_ID,between_PRE_2013_df$Condition, sep="_")


between_POST_2013_df<-as.data.frame(between_POST_2013)
between_POST_2013_df$Condition<-"Post"
between_POST_2013_df$Year<-"2013"
names(between_POST_2013_df)[names(between_POST_2013_df) == "between_POST_2013"] <- "Betweenness"
between_POST_2013_df<-setDT(between_POST_2013_df, keep.rownames = "Network_ID")[]
between_POST_2013_df$ID<-paste(between_POST_2013_df$Network_ID,between_POST_2013_df$Condition, sep="_")


EC_PRE_2013_df<-as.data.frame(EC_PRE_2013)
EC_PRE_2013_df$Condition<-"Pre"
EC_PRE_2013_df$Year<-"2013"
names(EC_PRE_2013_df)[names(EC_PRE_2013_df) == "EC_PRE_2013"] <- "Eigenvector_Centrality"
EC_PRE_2013_df<-setDT(EC_PRE_2013_df, keep.rownames = "Network_ID")[]
EC_PRE_2013_df$ID<-paste(EC_PRE_2013_df$Network_ID,EC_PRE_2013_df$Condition, sep="_")

EC_POST_2013_df<-as.data.frame(EC_POST_2013)
EC_POST_2013_df$Condition<-"Post"
EC_POST_2013_df$Year<-"2013"
names(EC_POST_2013_df)[names(EC_POST_2013_df) == "EC_POST_2013"] <- "Eigenvector_Centrality"
EC_POST_2013_df<-setDT(EC_POST_2013_df, keep.rownames = "Network_ID")[]
EC_POST_2013_df$ID<-paste(EC_POST_2013_df$Network_ID,EC_POST_2013_df$Condition, sep="_")


degree_PRE_2013_df<-as.data.frame(degree_PRE_2013)
degree_PRE_2013_df$Condition<-"Pre"
degree_PRE_2013_df$Year<-"2013"
names(degree_PRE_2013_df)[names(degree_PRE_2013_df) == "degree_PRE_2013"] <- "Degree"
degree_PRE_2013_df<-setDT(degree_PRE_2013_df, keep.rownames = "Network_ID")[]
degree_PRE_2013_df$ID<-paste(degree_PRE_2013_df$Network_ID,degree_PRE_2013_df$Condition, sep="_")


degree_POST_2013_df<-as.data.frame(degree_POST_2013)
degree_POST_2013_df$Condition<-"Post"
degree_POST_2013_df$Year<-"2013"
names(degree_POST_2013_df)[names(degree_POST_2013_df) == "degree_POST_2013"] <- "Degree"
degree_POST_2013_df<-setDT(degree_POST_2013_df, keep.rownames = "Network_ID")[]
degree_POST_2013_df$ID<-paste(degree_POST_2013_df$Network_ID,degree_POST_2013_df$Condition, sep="_")


degree_2013_PrePost<-rbind(degree_PRE_2013_df,degree_POST_2013_df)
between_2013_PrePost<-rbind(between_PRE_2013_df,between_POST_2013_df)
EC_2013_PrePost<-rbind(EC_PRE_2013_df,EC_POST_2013_df)


deg_betwn_2013<-merge(degree_2013_PrePost,between_2013_PrePost, by="ID")
deg_betwn_EC_2013<-merge(EC_2013_PrePost,deg_betwn_2013, by="ID")
write.csv(deg_betwn_EC_2013, "C:/Users/alexi/Google Drive (1)/Between Computers/CNWW2021/animal_response_nets/data/Fairywrens_LantzAndKarubian2017/deg_betwn_EC_2013_pre_post_fairywren.csv")

#individuals<-rbind()

#Group level


n_nodes_PRE_2013<-length(igraph::V(network2013_pre))
n_nodes_POST_2013<-length(igraph::V(network2013_post))

density_PRE_2013<-igraph::edge_density(network2013_pre)
density_POST_2013<-igraph::edge_density(network2013_post)

effective_info_PRE_2013<-einet::effective_information(network2013_pre, effectiveness = F)
EI_normal_pre<-effective_info_PRE_2013/n_nodes_PRE_2013
effective_info_POST_2013<-einet::effective_information(network2013_post, effectiveness = F)
EI_normal_post<-effective_info_POST_2013/n_nodes_POST_2013


