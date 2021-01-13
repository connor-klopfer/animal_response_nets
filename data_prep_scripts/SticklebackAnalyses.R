# Lisa P. Barrett 1/12/21 CNWW Team Before and After - I'm sorry this is SO inefficient!

# # FYI "In all but one of the association networks obtained for
# our groups, all fish were seen to associate with all others
# in their group at least once (fig. 1D). For this reason, many
# of the standard network metrics, such as path length, betweenness,
# and clustering coefficient, which are commonly
# used to describe networks with incomplete connectedness,
# were unnecessary or inapplicable here (see Croft et al.
#                                        2008). Instead, we focused on just three metrics: mean
# subgroup size, network density, and network differentiation."

#Outside of R I put the first 5 pre association matrices and the 5 post association matrices
#into separate CSV files

#load libraries
install.packages("igraph")
library(igraph)
install.packages("einet")
library(einet)

## PRE - UNSTRUCTURED , just the first 5 since post has five
#group 1_unstructured_pre-foraging task imported, no heading  SHOULD HAVE IMPORTED W/ HEADING so now I have to delete first column and first row
Group1_Unstructured_Pre
Group1_Unstructured_Pre2<- Group1_Unstructured_Pre[-c(1),-c(1)]

#make matrix and graph
m1<-as.matrix(Group1_Unstructured_Pre2)
g1<-graph.adjacency(m,mode="undirected",weighted=TRUE,diag=FALSE)
plot.igraph(g1)

#do for the rest of the pre- groups
Group2_Unstructured_Pre
Group2_Unstructured_Pre2<- Group2_Unstructured_Pre[-c(1),-c(1)]
m2<- as.matrix(Group2_Unstructured_Pre2)
g2<-graph.adjacency(m2,mode="undirected",weighted=TRUE,diag=FALSE)
Group3_Unstructured_Pre
Group3_Unstructured_Pre2<- Group3_Unstructured_Pre[-c(1),-c(1)]
m3<-as.matrix(Group3_Unstructured_Pre2)
g3<-graph.adjacency(m3,mode="undirected",weighted=TRUE,diag=FALSE)
Group4_Unstructured_Pre
Group4_Unstructured_Pre2<-Group4_Unstructured_Pre[-c(1),-c(1)]
m4<-as.matrix(Group4_Unstructured_Pre2)
g4<-graph.adjacency(m4,mode="undirected",weighted=TRUE,diag=FALSE)
Group5_Unstructured_Pre
Group5_Unstructured_Pre2<-Group5_Unstructured_Pre[-c(1),-c(1)]
m5<-as.matrix(Group5_Unstructured_Pre2)
g5<-graph.adjacency(m5,mode="undirected",weighted=TRUE,diag=FALSE)


#calculate measures for pre- groups
#first make a list of their matrices and a list of their graphs
matrixListPre <- list(m1,m2,m3,m4,m5)
graphListPre<-list(g1,g2,g3,g4,g5)


sample_function <- function(g){
  all_betweenness <- igraph::betweenness(g,directed=FALSE)
  return(mean(all_betweenness))
}
pre_between<-lapply(graphListPre, sample_function)


between_function <- function(g){
  all_eigenvector <- igraph::eigen_centrality(g, directed=FALSE)$vector
  return(mean(all_eigenvector))
}
pre_eigen<-lapply(graphListPre, between_function)


pre_density<-lapply(graphListPre, igraph::edge_density)


length_function<-function(g){
  all_length<-length(V(g))
  return(all_length)
}
pre_length<-lapply(graphListPre,length_function)

ei_function<-function(g){
  all_ei<-einet::effective_information(g,effectiveness=FALSE)
  return(all_ei/log2(8)) #all my networks have 8 nodes
}
pre_ei<-lapply(graphListPre,ei_function)



## POST - UNSTRUCTURED 
Group1_Unstructured_Post
Group1_Unstructured_Post2<- Group1_Unstructured_Post[-c(1),-c(1)]
m1_post<-as.matrix(Group1_Unstructured_Post2)
g1_post<-graph.adjacency(m1_post,mode="undirected",weighted=TRUE,diag=FALSE)
Group2_Unstructured_Post
Group2_Unstructured_Post2<- Group2_Unstructured_Post[-c(1),-c(1)]
m2_post<-as.matrix(Group2_Unstructured_Post2)
g2_post<-graph.adjacency(m2_post,mode="undirected",weighted=TRUE,diag=FALSE)
Group3_Unstructured_Post
Group3_Unstructured_Post2<- Group3_Unstructured_Post[-c(1),-c(1)]
m3_post<-as.matrix(Group3_Unstructured_Post2)
g3_post<-graph.adjacency(m3_post,mode="undirected",weighted=TRUE,diag=FALSE)
Group4_Unstructured_Post
Group4_Unstructured_Post2<- Group4_Unstructured_Post[-c(1),-c(1)]
m4_post<-as.matrix(Group4_Unstructured_Post2)
g4_post<-graph.adjacency(m4_post,mode="undirected",weighted=TRUE,diag=FALSE)
Group5_Unstructured_Post
Group5_Unstructured_Post2<- Group5_Unstructured_Post[-c(1),-c(1)]
m5_post<-as.matrix(Group5_Unstructured_Post2)
g5_post<-graph.adjacency(m5_post,mode="undirected",weighted=TRUE,diag=FALSE)


matrixListPost <- list(m1_post,m2_post,m3_post,m4_post,m5_post)
graphListPost<-list(g1_post,g2_post,g3_post,g4_post,g5_post)

#calculate measures for post- groups

sample_function <- function(g){
  all_betweenness <- igraph::betweenness(g,directed=FALSE)
  return(mean(all_betweenness))
}
post_between<-lapply(graphListPost, sample_function)


between_function <- function(g){
  all_eigenvector <- igraph::eigen_centrality(g, directed=FALSE)$vector
  return(mean(all_eigenvector))
}
post_eigen<-lapply(graphListPost, between_function)


post_density<-lapply(graphListPost, igraph::edge_density)


length_function<-function(g){
  all_length<-length(V(g))
  return(all_length)
}
post_length<-lapply(graphListPost,length_function)

ei_function<-function(g){
  all_ei<-einet::effective_information(g,effectiveness=FALSE)
  return(all_ei/log2(8)) #all my networks have 8 nodes
}
post_ei<-lapply(graphListPost,ei_function)


##combine pre and post measures into data frame 
stickleback_data<- data.frame(matrix(NA, nrow=10, ncol=9))
stickleback_data$Study<-rep("Stickleback", 10)
stickleback_data$Animal<-rep("Stickleback",10)
stickleback_data$Treatment<-c(rep("Pre",5), rep("Post",5))
stickleback_data$GroupID<-c(1:5,1:5)
stickleback_data$N_Nodes<-c(pre_length,post_length)
stickleback_data$Density<-c(pre_density,post_density)
stickleback_data$Mean_Eigenvector_Centrality<-c(pre_eigen,post_eigen)
stickleback_data$Mean_Betweenness_Centrality<-c(pre_between,post_between)
stickleback_data$Effective_Information<-c(pre_ei, post_ei)
stickleback_data<-stickleback_data[,-c(1:9)]
as.data.frame(stickleback_data)
str(stickleback_data)
getwd()
df <- apply(stickleback_data,2,as.character)
write.csv(df, "stickleback_data_output.csv")
