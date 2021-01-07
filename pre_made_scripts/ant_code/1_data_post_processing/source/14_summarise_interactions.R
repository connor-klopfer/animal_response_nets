####14_summarise_interactions.R#####

####Takes an interaction list as an input, and calculates:
###                     - the total duration of interactions of each ant to each of four categories of individuals (queen, nurses, untreated foragers and treated foragers)
###                     - the overall distribution of interactions according to task groups and ages

###Created by Nathalie Stroeymeyt

####################################
to_keep_ori <- to_keep

#################################
options(digits=16) ; options(digits.secs=6) ; options("scipen" = 10)

#### get input file list
input_path           <- paste(data_path,"/intermediary_analysis_steps/binned_interaction_lists",sep="")
setwd(input_path)  
input_folders        <- list.dirs(recursive=T,path="PreTreatment",full.names=F)
input_folders        <- input_folders[which(input_folders!="")]

outputfolder1 <- paste(data_path,"/processed_data/individual_behaviour/random_vs_observed",sep="")
if (!file.exists(outputfolder1)){dir.create(outputfolder1,recursive = T)}

summary_dol <- NULL
to_keep <- c(ls(),"to_keep","input_folder","network_file","network_files","summary_interactions","summary_pairs","all_interactions")
for (input_folder in input_folders){
  print(input_folder)
  setwd(input_path)
  network_files <- list.files(path=paste("PreTreatment/",input_folder,sep=""),full.names=T)
  if (input_folder=="observed"&grepl("main",data_path)){
    network_files <- c(network_files,list.files(path=paste("PostTreatment/",input_folder,sep=""),full.names=T))
  }
  summary_interactions <- NULL
  summary_pairs        <- NULL
  all_interactions     <- NULL
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
    time_hours         <- as.numeric(gsub("TH","",components[which(grepl("TH",components))]))
    time_of_day        <- as.numeric(gsub("TD","",components[which(grepl("TD",components))]))
    if (grepl("age",data_path)){
      colony_ages <- ages[which(ages$colony==colony),]
    }
    
    ####get appropriate task_group list, treated list and tag
    colony_treated     <- treated[which(treated$colony==colony_number),"tag"]
    colony_task_group  <- task_groups[which(task_groups$colony==colony),]
    tagfile            <- tag_list[which(grepl(colony,tag_list))]
    if (length(tagfile)>1){
      tagfile <- tagfile[grepl(components[grepl("Treatment",components)],tagfile)]
    }
    tag                <- read.tag(tagfile)$tag
    names(tag)[names(tag)=="#tag"] <- "tag"; tag <- tag[which(tag$tag!="#tag"),]
    tag[which(tag$age==0),"age"]   <- NA ###unknown ages are coded as 0 in the tag file
    tag <-tag[which(tag$final_status=="alive"),] ###remove dead ants from tag file
    ####read interactions
    interactions       <- read.table(network_file,header=T,stringsAsFactors = F)  
    
    #### add a column contaning interaction duration in min
    interactions["duration_min"] <- (interactions$Stoptime - interactions$Starttime + 0.5) /60 ###duration in minutes (one frame = 0.5 second)
    
    #### add a column containing the status of tag 1 and the status of tag2
    foragers <- colony_task_group[which(colony_task_group$task_group=="forager"),"tag"]
    nurses   <- colony_task_group[which(colony_task_group$task_group=="nurse"),"tag"]
    
    #####1. calculate within/between caste interactions for before period
    interactions[c("status_Tag1","status_Tag2")] <- NA
    
    interactions[which(interactions$Tag1%in%foragers),"status_Tag1"] <- "forager"
    interactions[which(interactions$Tag2%in%foragers),"status_Tag2"] <- "forager"
    
    interactions[which(interactions$Tag1%in%nurses),"status_Tag1"] <- "nurse"
    interactions[which(interactions$Tag2%in%nurses),"status_Tag2"] <- "nurse"
    
    interactions[which(interactions$Tag1==queenid),"status_Tag1"] <- "queen"
    interactions[which(interactions$Tag2==queenid),"status_Tag2"] <- "queen"
    
    if (period=="before"){
      inter <- interactions
      inter <- inter[,sort(names(inter))]
      all_interactions <- rbind(all_interactions,data.frame(randy=input_folder,colony_size=colony_size,period=period,inter))
    }
    
    # #####2. continue calculations for pre vs post
    interactions[which(interactions$Tag1%in%colony_treated),"status_Tag1"] <- "treated"
    interactions[which(interactions$Tag2%in%colony_treated),"status_Tag2"] <- "treated"
    
    
    #### use this information to calculate, for each worker, the cumulated duration of interaction with treated workers
    aggregated1                 <- aggregate(na.rm=T,na.action="na.pass",duration_min~Tag1+status_Tag2,FUN=sum,data=interactions)
    names(aggregated1)          <- c("tag","partner_status","duration_min")
    aggregated2                 <- aggregate(na.rm=T,na.action="na.pass",duration_min~Tag2+status_Tag1,FUN=sum,data=interactions)
    names(aggregated2)          <- c("tag","partner_status","duration_min")
    aggregated                  <- rbind(aggregated1,aggregated2)
    aggregated                  <- aggregate(na.rm=T,na.action="na.pass",duration_min~tag+partner_status,FUN=sum,data=aggregated)
    
    full_table                  <- expand.grid(tag=tag[which(tag$final_status=="alive"),"tag"],partner_status=c("queen","nurse","forager","treated"),stringsAsFactors = F)
    full_table                  <- merge(full_table,aggregated[c("tag","partner_status","duration_min")],all.x=T,all.y=T)
    full_table[is.na(full_table$duration_min),"duration_min"] <- 0
    
    full_table                  <- merge(full_table,tag[c("tag","group")]); names(full_table)[names(full_table)=="group"] <- "status"
    full_table                  <- merge(full_table,colony_task_group[c("tag","task_group")]); full_table[which(full_table$status=="treated"),"task_group"] <- "treated"
    if (!grepl("age",data_path)){
      full_table$age <- NA
    }else{
      full_table                <- merge(full_table,colony_ages,all.x=T,all.y=F)
    }
    full_table           <- full_table[c("tag","age","task_group","status","partner_status","duration_min")]
    summary_interactions <- rbind(summary_interactions,data.frame(randy=input_folder,colony=colony,colony_size=colony_size,treatment=treatment,period=period,time_of_day=time_of_day,time_hours=time_hours,full_table,stringsAsFactors = F))
    
    
    #####2. continue calculations for pre vs post
    interactions[which(interactions$Tag1%in%colony_treated),"status_Tag1"] <- "treated"
    interactions[which(interactions$Tag2%in%colony_treated),"status_Tag2"] <- "treated"
    
    
    #### use this information to calculate, for each worker, the cumulated duration of interaction with treated workers
    aggregated1                 <- aggregate(na.rm=T,na.action="na.pass",duration_min~Tag1+status_Tag2,FUN=sum,data=interactions)
    names(aggregated1)          <- c("tag","partner_status","duration_min")
    aggregated2                 <- aggregate(na.rm=T,na.action="na.pass",duration_min~Tag2+status_Tag1,FUN=sum,data=interactions)
    names(aggregated2)          <- c("tag","partner_status","duration_min")
    aggregated                  <- rbind(aggregated1,aggregated2)
    aggregated                  <- aggregate(na.rm=T,na.action="na.pass",duration_min~tag+partner_status,FUN=sum,data=aggregated)
    
    full_table                  <- expand.grid(tag=tag[which(tag$final_status=="alive"),"tag"],partner_status=c("queen","nurse","forager","treated"),stringsAsFactors = F)
    full_table                  <- merge(full_table,aggregated[c("tag","partner_status","duration_min")],all.x=T,all.y=T)
    full_table[is.na(full_table$duration_min),"duration_min"] <- 0
    
    full_table                  <- merge(full_table,tag[c("tag","group")]); names(full_table)[names(full_table)=="group"] <- "status"
    full_table                  <- merge(full_table,colony_task_group[c("tag","task_group")]); full_table[which(full_table$status=="treated"),"task_group"] <- "treated"
    if (!grepl("age",data_path)){
      full_table$age <- NA
    }else{
      full_table                <- merge(full_table,colony_ages,all.x=T,all.y=F)
    }
    full_table           <- full_table[c("tag","age","task_group","status","partner_status","duration_min")]
    summary_interactions <- rbind(summary_interactions,data.frame(randy=input_folder,colony=colony,colony_size=colony_size,treatment=treatment,period=period,time_of_day=time_of_day,time_hours=time_hours,full_table,stringsAsFactors = F))
    
    
    clean()
  }
  #####if folder = observed, use the summary_interactions table to compute inter-caste contacts #####
  if (grepl("main",data_path)&input_folder=="observed"){
    summary_interactions <- summary_interactions [which(!summary_interactions$partner_status%in%c("treated","queen")),]
    summary_interactions <- summary_interactions [which(!summary_interactions$task_group%in%c("treated","queen")),]
    
    summary_interactions["within_vs_between"] <- summary_interactions$partner_status==summary_interactions$task_group
    summary_interactions <- aggregate(na.rm=T,na.action="na.pass",duration_min~.,FUN=sum,data=summary_interactions[c("colony","period","time_of_day","time_hours","tag","task_group","status","partner_status","duration_min","within_vs_between")])
    summary_interactions_between <- summary_interactions[which(!summary_interactions$within_vs_between),!names(summary_interactions)%in%c("partner_status","within_vs_between")];
    names(summary_interactions_between)[names(summary_interactions_between)=="duration_min"] <- "inter_caste_contact_duration"
    
    #####add information to individual behaviour file
    behav <- read.table(paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",sep=""),header=T,stringsAsFactors = F)
    behav <- merge(behav,summary_interactions_between[c("colony","tag","time_hours","inter_caste_contact_duration")],all.x=T,all.y=F,sort=F)
    behav <- behav[order(behav$colony,behav$tag,behav$time_hours),]
    write.table(behav,file=paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",sep=""), row.names=F, col.names=T,append=F,quote=F)
  }
  
  ###Use all_interactions to obtain information about age-based DoL before treatment
  ###first check that all_interactions only has Pre-treatment
  all_interactions <- all_interactions[which(all_interactions$period=="before"),]
  ###summ all interactions for each pair of ants
  all_interactions <- aggregate(na.rm=T,na.action="na.pass",duration_min~randy+colony+colony_size+period+Tag1+Tag2+status_Tag1+status_Tag2+treatment,FUN=sum,data=all_interactions)
  ###calculate intra_caste_over_inter_caste_WW_contact_duration
  all_interactions["same_caste"] <- all_interactions$status_Tag1==all_interactions$status_Tag2
  same_caste_interactions <- aggregate(na.rm=T,na.action="na.pass",duration_min~randy+colony+colony_size+treatment+period,FUN=sum,data=all_interactions[which((all_interactions$status_Tag1!="queen")&(all_interactions$status_Tag2!="queen")&(all_interactions$same_caste)),])
  names(same_caste_interactions)[names(same_caste_interactions)=="duration_min"] <- "duration_min_within"
  inter_caste_interactions <- aggregate(na.rm=T,na.action="na.pass",duration_min~randy+colony+colony_size+treatment+period,FUN=sum,data=all_interactions[which((all_interactions$status_Tag1!="queen")&(all_interactions$status_Tag2!="queen")&(!all_interactions$same_caste)),])
  names(inter_caste_interactions)[names(inter_caste_interactions)=="duration_min"] <- "duration_min_between"
  inter_intra_caste_interactions <- merge(same_caste_interactions,inter_caste_interactions,all.x=T,all.y=T)
  inter_intra_caste_interactions["intra_caste_over_inter_caste_WW_contact_duration"] <- inter_intra_caste_interactions$duration_min_within/inter_intra_caste_interactions$duration_min_between
  dol <-   inter_intra_caste_interactions[!names(inter_intra_caste_interactions)%in%c("duration_min_within","duration_min_between")]
  
  ###calculate queen contact with nurses vs. workers
  queen_interactions <- all_interactions[which(all_interactions$status_Tag1=="queen"|all_interactions$status_Tag2=="queen"),]
  queen_interactions[which(queen_interactions$status_Tag1=="queen"),"partner"] <- queen_interactions[which(queen_interactions$status_Tag1=="queen"),"Tag2"]
  queen_interactions[which(queen_interactions$status_Tag1=="queen"),"partner_status"] <- queen_interactions[which(queen_interactions$status_Tag1=="queen"),"status_Tag2"]
  queen_interactions[which(queen_interactions$status_Tag2=="queen"),"partner"] <- queen_interactions[which(queen_interactions$status_Tag2=="queen"),"Tag1"]
  queen_interactions[which(queen_interactions$status_Tag2=="queen"),"partner_status"] <- queen_interactions[which(queen_interactions$status_Tag2=="queen"),"status_Tag1"]
  
  interaction_with_nurses <-aggregate (na.rm=T,na.action="na.pass",duration_min~randy+colony+colony_size+period+treatment,FUN=sum,data=queen_interactions[which(queen_interactions$partner_status=="nurse"),])
  names(interaction_with_nurses)[names(interaction_with_nurses)=="duration_min"] <- "duration_min_with_nurses"
  interaction_with_forager <-aggregate (na.rm=T,na.action="na.pass",duration_min~randy+colony+colony_size+period+treatment,FUN=sum,data=queen_interactions[which(queen_interactions$partner_status=="forager"),])
  names(interaction_with_forager)[names(interaction_with_forager)=="duration_min"] <- "duration_min_with_foragers"
  queen_interac <- merge(interaction_with_nurses,interaction_with_forager,all.x=T,all.y=T)
  queen_interac["QNurse_over_QForager_contact_duration"] <- queen_interac$duration_min_with_nurses/queen_interac$duration_min_with_foragers
  dol <- merge(dol,queen_interac[c("randy","colony","period","QNurse_over_QForager_contact_duration")])
  ###if necessary: add age
  if (grepl("age",data_path)){
    partner_ages <- ages; names(partner_ages) <- c("colony","partner","partner_age")
    queen_interactions <- merge(queen_interactions,partner_ages,all.x=T,all.y=F)
    
    ages_Tag1                    <- ages; names(ages_Tag1) <- c("colony","Tag1","age_Tag1")
    ages_Tag2                    <- ages; names(ages_Tag2) <- c("colony","Tag2","age_Tag2")
    all_interactions             <- merge(all_interactions,ages_Tag1,all.x=T,all.y=F)
    all_interactions             <- merge(all_interactions,ages_Tag2,all.x=T,all.y=F)
    all_interactions["age_diff"] <- abs(all_interactions$age_Tag2-all_interactions$age_Tag1)
    
    ###write ordered pair of interacting ants
    all_interactions["ordered"]<- all_interactions[,"Tag1"]<all_interactions[,"Tag2"]
    all_interactions[which(all_interactions$ordered),"new_Tag1"] <-  all_interactions[which(all_interactions$ordered),"Tag1"]
    all_interactions[which(all_interactions$ordered),"new_Tag2"] <-  all_interactions[which(all_interactions$ordered),"Tag2"]
    all_interactions[which(!all_interactions$ordered),"new_Tag1"] <-  all_interactions[which(!all_interactions$ordered),"Tag2"]
    all_interactions[which(!all_interactions$ordered),"new_Tag2"] <-  all_interactions[which(!all_interactions$ordered),"Tag1"]
    all_interactions["pair"]   <- as.character(interaction(all_interactions$colony,all_interactions$new_Tag1,all_interactions$new_Tag2))
    ###to get slope: get each pair of possibly interacting ants
    for (colony in unique(all_interactions$colony)){
      colony_ages <- ages[which(ages$colony==colony),]
      colony_ages_Tag1                    <- colony_ages; names(colony_ages_Tag1) <- c("colony","Tag1","age_Tag1")
      colony_ages_Tag2                    <- colony_ages; names(colony_ages_Tag2) <- c("colony","Tag2","age_Tag2")
      tagfile            <- tag_list[which(grepl(colony,tag_list))][1]
      tag <- read.tag(tagfile)$tag
      names(tag)[names(tag)=="#tag"] <- "tag"; tag <- tag[which(tag$tag!="#tag"),]
      
      ####list ants of known ages
      known_ages                <- colony_ages[which(!is.na(colony_ages$age)),"tag"]
      ###remove queen fronm list
      known_ages                <- known_ages[known_ages!=queenid]
      ###remove dead ants from list
      known_ages                <- known_ages[which(known_ages%in%tag[which(tag$final_status=="alive"),"tag"])]
      
      ### prepare complete pair table
      summ_pairs           <- data.frame(combinations(n=length(known_ages),r=2,v=known_ages))
      summ_pairs["ordered"]<- summ_pairs[,1]<summ_pairs[,2]
      summ_pairs[which(summ_pairs$ordered),"Tag1"] <-  summ_pairs[which(summ_pairs$ordered),1]
      summ_pairs[which(summ_pairs$ordered),"Tag2"] <-  summ_pairs[which(summ_pairs$ordered),2]
      summ_pairs[which(!summ_pairs$ordered),"Tag1"] <-  summ_pairs[which(!summ_pairs$ordered),2]
      summ_pairs[which(!summ_pairs$ordered),"Tag2"] <-  summ_pairs[which(!summ_pairs$ordered),1]
      summ_pairs              <- merge(summ_pairs,colony_ages_Tag1,all.x=T,all.y=F)
      summ_pairs              <- merge(summ_pairs,colony_ages_Tag2,all.x=T,all.y=F)
      summ_pairs["age_diff"]  <- abs(summ_pairs$age_Tag2-summ_pairs$age_Tag1)
      summ_pairs["pair"]   <- as.character(interaction(colony,summ_pairs$Tag1,summ_pairs$Tag2))
      summ_pairs <- merge(summ_pairs,unique(all_interactions[which(all_interactions$colony==colony),c("colony","randy","colony_size","period","treatment")]))
      ###merge it with interactions whose age diff is known 
      summ_pairs <- merge(summ_pairs,all_interactions[which(all_interactions$colony==colony),c("colony","randy","colony_size","period","treatment","pair","duration_min")],all.x=T,all.y=F)
      ###and fill intreactions which did not happen with 0
      summ_pairs[which(is.na(summ_pairs$duration_min)),"duration_min"] <- 0
      model_WW <- lm(duration_min~age_diff,data=summ_pairs)
      dol[dol$colony==colony,"slope_WW_contact_duration_f_age_diff"] <- coef(model_WW)["age_diff"]
      
      ###do the same for queen interactions
      colony_ages <- colony_ages[which(!is.na(colony_ages$age)),]
      names(colony_ages) <-c("colony","partner","age")
      queen_W <- merge(colony_ages,unique(queen_interactions[which(queen_interactions$colony==colony),c("colony","randy","colony_size","period","treatment")]))
      ###merge it with interactions whose age diff is known 
      queen_W <- merge(queen_W,queen_interactions[which(queen_interactions$colony==colony),c("colony","randy","colony_size","period","treatment","partner","duration_min")],all.x=T,all.y=F)
      ###and fill intreactions which did not happen with 0
      queen_W[which(is.na(queen_W$duration_min)),"duration_min"] <- 0
      model_QW <- lm(duration_min~age,data=queen_W)
      dol[dol$colony==colony,"slope_QW_contact_duration_f_W_age"] <- coef(model_QW)["age"]
      
    }
  }
  summary_dol <- rbind(summary_dol,dol)  
  
}
if(!file.exists(paste(data_path,"/processed_data/collective_behaviour/random_vs_observed",sep=""))){dir.create(paste(data_path,"/processed_data/collective_behaviour/random_vs_observed",sep=""),recursive=T)}
if(!file.exists(paste(data_path,"/processed_data/collective_behaviour/random_vs_observed/interactions.dat",sep=""))){write.table(summary_dol,file=paste(data_path,"/processed_data/collective_behaviour/random_vs_observed/interactions.dat",sep=""),col.names=T,row.names=F,append=F,quote=F)}
to_keep <- to_keep_ori
