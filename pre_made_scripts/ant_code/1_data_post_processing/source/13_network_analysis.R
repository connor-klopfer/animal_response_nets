####13_network_analysis.R#####

#### Takes an interaction list as an input, builds a network, and analyse its properties

###Created by Nathalie Stroeymeyt

####################################
to_keep_ori <- to_keep

#################################
options(digits=16) ; options(digits.secs=6) ; options("scipen" = 10)

###remove output file
if (file.exists(paste(data_path,"/processed_data/individual_behaviour/post_treatment/interactions_with_treated.txt",sep=""))){
  file.remove(paste(data_path,"/processed_data/individual_behaviour/post_treatment/interactions_with_treated.txt",sep=""))
}

#### get input file list
if (!grepl("survival",data_path)){
  input_path           <- paste(data_path,"/intermediary_analysis_steps/binned_interaction_lists",sep="")
  setwd(input_path)  
  input_folders        <- list.dirs(recursive=T,path="PreTreatment",full.names=F)
  input_folders        <- input_folders[which(input_folders!="")]
}else{
  input_path           <- paste(data_path,"/intermediary_analysis_steps/full_interaction_lists",sep="")
  setwd(input_path)  
  input_folders        <- list.dirs(recursive=T,path="PostTreatment",full.names=F)
  input_folders        <- input_folders[which(input_folders!="")]
}

queen_community_summary <- NULL
to_keep <- c(ls(),"to_keep","input_folder","network_files","options","option","summary_collective","summary_individual","outputfolder","network_file")
for (input_folder in input_folders){
  print(input_folder)
  setwd(input_path)
  network_files <- list.files(path=paste("PreTreatment/",input_folder,sep=""),full.names=T)
  if (input_folder=="observed"){
    network_files <- c(network_files,list.files(path=paste("PostTreatment/",input_folder,sep=""),full.names=T))
    if(grepl("main",data_path)){
      options <- c("all_workers","untreated_only")
    }else{
      options <- c("all_workers")
    }
  }else{
    options <- c("all_workers")
  }
  for (option in options){
    print(option)
    outputfolder <- paste(data_path,"/processed_data/network_properties",sep="")
    
    summary_collective <- NULL
    summary_individual <- NULL
    for (network_file in network_files){
      print(network_file)
      ####get file metadata
      root_name          <- gsub("_interactions.txt","",unlist(strsplit(network_file,split="/"))[grepl("colony",unlist(strsplit(network_file,split="/")))])
      components         <- unlist(strsplit(root_name,split="_"))
      colony             <- components[grepl("colony",components)]
      colony_number      <- as.numeric(gsub("colony","",colony))
      treatment          <- info[which(info$colony==colony_number),"treatment"]
      colony_size        <- info[which(info$colony==colony_number),"colony_size"]
      
      if (!all(!grepl("PreTreatment",components))){period <- "before"}else{period <- "after"}
      if(!grepl("survival",data_path)){
        time_hours         <- as.numeric(gsub("TH","",components[which(grepl("TH",components))]))
        time_of_day        <- as.numeric(gsub("TD","",components[which(grepl("TD",components))]))
      }else{
        if (period=="after"){
          time_hours   <- 0
          time_of_day <- 12
        }else{
          time_hours   <- -30
          time_of_day <- 6
        }
      }
      
      ####get appropriate task_group list, treated list and tag
      colony_treated     <- treated[which(treated$colony==colony_number),"tag"]
      colony_task_group  <- task_groups[which(task_groups$colony==colony),]
      tagfile            <- tag_list[which(grepl(colony,tag_list))]
      if (length(tagfile)>1){
        tagfile <- tagfile[grepl(unlist(strsplit(gsub("\\.txt","",root_name),"_"))[grepl("Treatment",unlist(strsplit(gsub("\\.txt","",root_name),"_")))],tagfile)]
      }
      
      ####read interactions
      interactions       <- read.table(network_file,header=T,stringsAsFactors = F)
      
      tag                <- read.tag(tagfile)$tag
      names(tag)[names(tag)=="#tag"] <- "tag"; tag <- tag[which(tag$tag!="#tag"),]
      tag[which(tag$age==0),"age"]   <- NA ###unknown ages are coded as 0 in the tag file
      if (!grepl("survival",data_path)){
        alive      <- tag[which(tag$final_status=="alive"),"tag"]
      }else{
        alive      <- tag[which(as.numeric(tag$death)==0|as.numeric(tag$death)>=max(interactions$Stopframe,na.rm=T)),"tag"]
      }
      tag <-tag[which(tag$tag%in%alive),] ###remove dead ants from tag file
      
      
      ####if untreated only, reduce tag , colony_task_group, and interactions
      if (option=="untreated_only"){
        colony_task_group <- colony_task_group[which(!colony_task_group$tag%in%colony_treated),]
        tag               <- tag[which(!tag$tag%in%colony_treated),]
        interactions      <- interactions[which(!((interactions$Tag1%in%colony_treated)|(interactions$Tag2%in%colony_treated))),]
      }
      
      actors             <- data.frame(name=as.character(tag$tag))
      
      #### add a column contaning interaction duration in min
      interactions["duration_min"] <- (interactions$Stoptime - interactions$Starttime + 0.5) /60 ###duration in minutes (one frame = 0.5 second)
      
      #### add a column containing the status of tag 1 and the status of tag2
      interactions[c("status_Tag1","status_Tag2")] <- "untreated"
      interactions[which(interactions$Tag1%in%colony_treated),"status_Tag1"] <- "treated"
      interactions[which(interactions$Tag2%in%colony_treated),"status_Tag2"] <- "treated"
      
      #### use this information to calculate, for each worker, the cumulated duration of interaction with treated workers
      if (input_folder=="observed"&option=="all_workers"){
        aggregated1                 <- aggregate(na.rm=T,na.action="na.pass",duration_min~Tag1+status_Tag2,FUN=sum,data=interactions[which(interactions$status_Tag2=="treated"),])
        names(aggregated1)          <- c("tag","partner_status","duration_min")
        aggregated2                 <- aggregate(na.rm=T,na.action="na.pass",duration_min~Tag2+status_Tag1,FUN=sum,data=interactions[which(interactions$status_Tag1=="treated"),])
        names(aggregated2)          <- c("tag","partner_status","duration_min")
        aggregated                  <- rbind(aggregated1,aggregated2)
        aggregated                  <- aggregate(na.rm=T,na.action="na.pass",duration_min~tag+partner_status,FUN=sum,data=aggregated)
        interactions_with_treated   <- merge(data.frame(tag=tag[which(tag$tag%in%alive),"tag"],stringsAsFactors = F),aggregated[c("tag","duration_min")],all.x=T)
        interactions_with_treated[is.na(interactions_with_treated$duration_min),"duration_min"] <- 0
        names(interactions_with_treated) <- c("tag","duration_of_contact_with_treated_min")
        interactions_with_treated["colony"] <- colony
        interactions_with_treated["time_hours"] <- time_hours
        ###write results
        ###individual behaviour: pre_vs_post treatment
        if (grepl("main",data_path)){
          behav <- read.table(paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",sep=""),header=T,stringsAsFactors = F)
          if (!"duration_of_contact_with_treated_min"%in%names(behav)){
            behav[c("duration_of_contact_with_treated_min")] <- NA
          }
          behav[match(as.character(interaction(interactions_with_treated$colony,interactions_with_treated$tag,interactions_with_treated$time_hours)),as.character(interaction(behav$colony,behav$tag,behav$time_hours))),c("duration_of_contact_with_treated_min")]  <- interactions_with_treated$duration_of_contact_with_treated_min
          options(digits=3)
          write.table(behav,file=paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",sep=""), row.names=F, col.names=T,append=F,quote=F)
          options(digits=16)
        }
        ###interactions with treated: post-treatment
        if (period=="after"){
          outputfoldy <- paste(data_path,"/processed_data/individual_behaviour/post_treatment",sep="")
          if(!file.exists(outputfoldy)){dir.create(outputfoldy,recursive=T)}
          int_with_treated <- data.frame(colony_size=colony_size,treatment=treatment,period=period,time_of_day=time_of_day,
                                         interactions_with_treated,stringsAsFactors = F)
          if (!file.exists(paste(data_path,"/processed_data/individual_behaviour/post_treatment/interactions_with_treated.txt",sep=""))){
            write.table(int_with_treated,file=paste(outputfoldy,"/interactions_with_treated.txt",sep=""),col.names=T,row.names=F,quote=F,append=F)
          }else{
            write.table(int_with_treated,file=paste(outputfoldy,"/interactions_with_treated.txt",sep=""),col.names=F,row.names=F,quote=F,append=T)
          }
        }
      }
      
      ###build NETWORK
      if (!grepl("survival",data_path)){
        net <- graph.data.frame(interactions[c("Tag1","Tag2")],directed=F,vertices=actors)
        ### add edge weights
        E(net)$weight <- interactions[,"duration_min"]
        ###simplify graph (merge all edges involving the same pair of ants into a single one whose weight = sum of these weights)
        net <- simplify(net,remove.multiple=TRUE,remove.loop=TRUE,edge.attr.comb="sum")
        ##################remove unconnected nodes
        unconnected <- actors[degree(net)==0,]
        net <- net - as.character(unconnected)
        ##################update actor list
        actors <- get.vertex.attribute(net,"name")
        
        ####Part 1: collective network properties ####
        ##Assortativity  - Age
        ####if age experiment, get colony ages
        if (grepl("age",data_path)){
          colony_ages <- ages [which(ages$colony==colony),]
          ####set queen age to NA as this would bias the result (the queen is the oldest individual and interacts mostly with the young nurses)
          colony_ages[which(colony_ages$tag==queenid),"age"] <- NA
          ####order the age acording to the order of the network's vertices
          ordered_ages        <- colony_ages[match(V(net)$name,as.character(colony_ages$tag)),"age"]
          #### calculate age assortativity
          age_assortativity <- assortativity(net-V(net)$name[is.na(ordered_ages)],types1=ordered_ages[!is.na(ordered_ages)],directed=F)
        }else{
          age_assortativity <- NA
        }
        ## Assortativity - Task
        ordered_task_groups <- colony_task_group[match(actors,as.character(colony_task_group$tag)),"task_group"]
        ordered_task_groups <- as.numeric(as.factor(ordered_task_groups))
        task_assortativity  <- assortativity_nominal(net,types=ordered_task_groups,directed=F)
        ##Clustering
        clustering <- mean(transitivity(net,type="barrat",weights=E(net)$weight,isolates = c("NaN")),na.rm=T)
        ##Degree mean and max
        degrees         <- degree(net,mode="all")
        degree_mean     <- mean(degrees,na.rm=T)
        degree_maximum  <- max(degrees,na.rm=T)
        ##Density
        density  <- igraph::edge_density(net)
        ##Diameter
        diameter <- igraph::diameter(net,directed=F,unconnected=TRUE,weights=(1/E(net)$weight)) ###here use the inverse of the weights, because the algorithm considers weights as distances rather than strengths of connexion
        ##Efficiency
        net_dist                    <- shortest.paths(net, weights=1/E(net)$weight, mode="all") ##again use the inverse of the weights, because the algorithm considers weights as distances rather than strengths of connexion
        net_dist[net_dist==0]       <- NA ##remove distances to self
        efficiency                  <- 1/net_dist ##transform each distance into an efficiency
        efficiency <- (1/((vcount(net)*(vcount(net)-1))))*(sum(efficiency,na.rm=TRUE))
        ## Modularity
        communities             <- cluster_louvain(net, weights = E(net)$weight)
        community_membership    <- communities$membership
        modularity              <- modularity(net,community_membership,weights=E(net)$weight)
        
        ###Add to data
        summary_collective <- rbind(summary_collective,data.frame(randy=input_folder,colony=colony,colony_size=colony_size,treatment=treatment,period=period,time_hours=time_hours,time_of_day=time_of_day,
                                                                  age_assortativity=age_assortativity,
                                                                  task_assortativity=task_assortativity,
                                                                  clustering=clustering,
                                                                  degree_mean=degree_mean,
                                                                  degree_maximum=degree_maximum,
                                                                  density=density,
                                                                  diameter=diameter,
                                                                  efficiency=efficiency,
                                                                  modularity=modularity,stringsAsFactors = F))
        
        ####Part 2: individual network properties ####
        ###prepare table
        tag["status"] <- "untreated"; tag[which(tag$tag%in%colony_treated),"status"] <- "treated"
        individual <- data.frame(randy=input_folder,colony=colony,colony_size=colony_size,treatment=treatment,tag=tag$tag,age=tag$age,status=tag$group,period=period,time_hours=time_hours,time_of_day=time_of_day,
                                 degree=NA,
                                 aggregated_distance_to_queen=NA,
                                 mean_aggregated_distance_to_treated=NA,
                                 same_community_as_queen=NA)
        ##degree
        individual[match(names(degrees),individual$tag),"degree"] <- degrees
        ##same community as queen
        queen_comm <- community_membership[which(V(net)$name==queenid)]
        community_membership <- community_membership==queen_comm
        individual[match(V(net)$name,individual$tag),"same_community_as_queen"] <- community_membership
        ##path length to queen
        if (queenid%in%actors){
          path_length_to_queen <- t(shortest.paths(net,v=actors,to="665",weights=1/E(net)$weight))
          individual[match(colnames(path_length_to_queen),individual$tag),"aggregated_distance_to_queen"] <- as.numeric(path_length_to_queen )
        }
        ########Mean path length to treated; aggregated_network
        if(option!="untreated_only"){
          path_length_to_treated                             <- as.data.frame(as.matrix(shortest.paths(net,v=actors,to=as.character(colony_treated)[as.character(colony_treated)%in%V(net)$name],weights=1/E(net)$weight)))
          path_length_to_treated["mean_distance_to_treated"] <- NA
          path_length_to_treated$mean_distance_to_treated    <- as.numeric(rowMeans(path_length_to_treated,na.rm=T))
          individual[match(rownames(path_length_to_treated),individual$tag),"mean_aggregated_distance_to_treated"] <- path_length_to_treated[,"mean_distance_to_treated"]
        }
        ###Add data to main data table
        summary_individual <- rbind(summary_individual,individual)
      }
      clean()
    }
    #####write #####
    if (!grepl("survival",data_path)){
      if (input_folder=="observed"){
        ####Main experiment: write pre_vs_post_treatment data
        if (!grepl("age",data_path)){
          outputfolder2 <- paste(outputfolder,"pre_vs_post_treatment",option,sep="/")
          if(!file.exists(outputfolder2)){dir.create(outputfolder2,recursive=T)}
          write.table(summary_collective[,names(summary_collective)!="randy"],file=paste(outputfolder2,"/colony_data.txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
          write.table(summary_individual[,names(summary_individual)!="randy"],file=paste(outputfolder2,"/individual_data.txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
        }
        ####All workers: write pre_treatment data into random_vs_observed folder
        if (option=="all_workers"){
          outputfolder3 <- paste(outputfolder,"random_vs_observed",sep="/")
          if(!file.exists(outputfolder3)){dir.create(outputfolder3,recursive=T)}
          write.table(summary_collective[which(summary_collective$period=="before"),],file=paste(outputfolder3,"/network_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
          write.table(summary_individual[which(summary_individual$period=="before"),],file=paste(outputfolder3,"/node_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
          
          outputfolder4 <- paste(outputfolder,"post_treatment",sep="/")
          if(!file.exists(outputfolder4)){dir.create(outputfolder4,recursive=T)}
          write.table(summary_collective[which(summary_collective$period=="after"),],file=paste(outputfolder4,"/network_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
          write.table(summary_individual[which(summary_individual$period=="after"),],file=paste(outputfolder4,"/node_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
        }
        ######Main experiment, All workers: add pre_treatment node_properties information to pre_treatment behaviour file
        if (!grepl("age",data_path)&option=="all_workers"){
          pre_treatment_behav_file <- paste(data_path,"/processed_data/individual_behaviour/pre_treatment/network_position_vs_time_outside.dat",sep="")
          pre_treatment_behav      <- read.table(pre_treatment_behav_file,header=T,stringsAsFactors = F)
          pre_treatment_behav      <- merge(pre_treatment_behav,summary_individual[which(summary_individual$period=="before"),c("colony","tag","time_hours","degree","aggregated_distance_to_queen")],all.x=T,all.y=T)
          pre_treatment_behav      <- pre_treatment_behav[order(pre_treatment_behav$colony,pre_treatment_behav$tag,pre_treatment_behav$time_hours),]
          write.table(pre_treatment_behav, file=pre_treatment_behav_file,col.names=T,row.names=F,quote=F,append=F)
        }
      }else{
        outputfolder3 <- paste(outputfolder,"random_vs_observed",sep="/")
        if(!file.exists(outputfolder3)){dir.create(outputfolder3,recursive=T)}
        write.table(summary_collective[which(summary_collective$period=="before"),],file=paste(outputfolder3,"/network_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
        write.table(summary_individual[which(summary_individual$period=="before"),],file=paste(outputfolder3,"/node_properties_",input_folder,".txt",sep=""),col.names = T,row.names=F,append=F,quote=F)
      }
    }
    ###Get characteristics of queen community vs. other communities (worker age, prop. of foragers)
    if (!grepl("survival",data_path)&option=="all_workers"){
      summary_individual_before <- read.table(paste(outputfolder,"/random_vs_observed/node_properties_",input_folder,".txt",sep=""),header = T,stringsAsFactors = F)
      summary_individual_before <- summary_individual_before[which(summary_individual_before$period=="before"),]
      ####if necessary: add age
      if (grepl("age",data_path)){
        summary_individual_before <- merge(summary_individual_before[,which(names(summary_individual_before)!="age")],ages,all.x=T,all.y=F)
      }else{
        summary_individual_before$age <- NA
      }
      ####add task_group
      summary_individual_before <- merge(summary_individual_before,task_groups,all.x=T,all.y=F)
      ###remove queen
      summary_individual_before <- summary_individual_before[which(summary_individual_before$tag!=queenid),]
      
      ###1. calculate mean proportion of foragers depending on with vs without queen
      summary_individual_before["forager"] <- 0
      summary_individual_before[which(summary_individual_before$task_group=="forager"),"forager"] <- 1
      prop_foragers <- aggregate(na.rm=T,na.action="na.pass",forager~colony+randy+colony_size+treatment+period+time_hours+same_community_as_queen,FUN=mean,data=summary_individual_before)
      prop_foragers <- aggregate(na.rm=T,na.action="na.pass",forager~colony+randy+colony_size+treatment+period+same_community_as_queen,FUN=mean,data=prop_foragers)
      names(prop_foragers)[names(prop_foragers)=="same_community_as_queen"] <- "in_queen_comm";names(prop_foragers)[names(prop_foragers)=="forager"] <- "proportion_of_foragers"
      ###2. calculate mean age of workers depending on with vs without queen
      if (grepl("age",data_path)){
        mean_age <- aggregate(na.rm=T,na.action="na.pass",age~colony+randy+colony_size+treatment+period+time_hours+same_community_as_queen,FUN=mean,data=summary_individual_before)
        mean_age <- aggregate(na.rm=T,na.action="na.pass",age~colony+randy+colony_size+treatment+period+same_community_as_queen,FUN=mean,data=mean_age)
        names(mean_age)[names(mean_age)=="same_community_as_queen"] <- "in_queen_comm"
        prop_foragers <- merge(prop_foragers,mean_age,all.x=T)
      }else{
        prop_foragers$age <- NA
      }
      prop_foragers[which(prop_foragers$in_queen_comm=="FALSE"),"in_queen_comm"] <- "not_with_queen"
      prop_foragers[which(prop_foragers$in_queen_comm=="TRUE"),"in_queen_comm"] <- "with_queen"
      if (grepl("random",input_folder)){
        prop_foragers["randy"] <- "random"
      }
      queen_community_summary <- rbind(queen_community_summary,prop_foragers)
    }
  }
}
if (!grepl("survival",data_path)){
  queen_community_summary <- aggregate(na.rm=T,na.action="na.pass",cbind(proportion_of_foragers,age)~.,FUN=mean,data=queen_community_summary)
  queen_community_summary <- queen_community_summary[order(queen_community_summary$randy,queen_community_summary$colony),]
  queen_community_summary$treatment <- queen_community_summary$randy
  if (!file.exists(paste(data_path,"/processed_data/network_properties/random_vs_observed",sep=""))){dir.create(paste(data_path,"/processed_data/network_properties/random_vs_observed",sep=""),recursive=T)}
  write.table(queen_community_summary,file=paste(data_path,"/processed_data/network_properties/random_vs_observed/queen_community.dat",sep=""),append=F,quote=F,row.names=F,col.names=T)
}
to_keep <- to_keep_ori
