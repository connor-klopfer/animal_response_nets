####6_time_investment.R#####

### Calls C++ function time_investment (https://github.com/laurentkeller/anttrackingUNIL)
### Takes a datfile, a tagfile, a plumefile, and a time period as inputs, and returns the number of frames spent in each zone defined in the plume file

###Created by Nathalie Stroeymeyt

#################################3
####get dat list #####
input_dat <- paste(data_path,"intermediary_analysis_steps/dat_files",sep="/")
setwd(input_dat)
dat_list <- paste(input_dat,list.files(pattern="dat"),sep="/")
if (grepl("age",data_path)){
  dat_list <- dat_list[grepl("PreTreatment",dat_list)]
}

#### define outputfolders
outputfolder <- paste(data_path,"/processed_data/time_investment",sep="")
if (!file.exists(outputfolder)){dir.create(outputfolder,recursive=T)}
outputfile <- paste(outputfolder,"/temp.txt",sep="")
outputspace_file <- paste(outputfolder,"/time_investment.txt",sep="")

space_dat     <- NULL
behaviour_dat <- NULL
for (datfile in dat_list){
  root_name   <- gsub("\\.dat","",unlist(strsplit(datfile,split="/"))[grepl("colony",unlist(strsplit(datfile,split="/")))])
  colony      <- unlist(strsplit(root_name,split="_"))[grepl("colony",unlist(strsplit(root_name,split="_")))]
  tagfile     <- tag_list[grepl(colony,tag_list)]
  if (length(tagfile)>1){
    tagfile <- tagfile[grepl(unlist(strsplit(root_name,"_"))[grepl("Treatment",unlist(strsplit(root_name,"_")))],tagfile)]
  }
  
  split_file  <- split_list[grepl(root_name,split_list)]
  if (grepl("age",data_path)){
    colony_ages <- ages[which(ages$colony==colony),]
  }
  split_info  <- read.table(split_file,header=T,stringsAsFactors = F)
  
  tag   <- read.tag(tagfile)$tag
  names(tag)[names(tag)=="#tag"] <- "tag"
  tag <- tag[which(tag$tag!="#tag"),]
  tag[which(tag$age==0),"age"] <- NA

  
  plume_info  <- info_plume[which(info_plume$dat_file==paste(root_name,".dat",sep="")),]
  plume_file  <- paste(data_path,"/original_data/plume_files/",plume_info["plume_file"],sep="")
  box         <- info[which(info$colony==as.numeric(gsub("colony","",colony))),"box"]
  treatment   <- info[which(info$colony==as.numeric(gsub("colony","",colony))),"treatment"]
  colony_size <- info[which(info$colony==as.numeric(gsub("colony","",colony))),"colony_size"]
  final_frame <- info_datfile[which(info_datfile$dat_file==paste(root_name,".dat",sep="")),"end_frame"]
  if (grepl("PreTreatment",root_name)){
    period <- "before"
  }else {
    period <- "after"
  }
  
  
  
  for (i in 1:nrow(split_info)){
    start_frame <- split_info [i,"frames"]
    if (i<nrow(split_info)){
      duration  <- split_info [i+1,"frames"] - start_frame
    }else{
      duration  <- final_frame-start_frame+1
    }
    ###Run time investment
    command <- paste(executables_path,"/time_investment -i ",datfile," -t ",tagfile," -p ",plume_file," -z nest -z water -z sugar -b ",box," -s ",start_frame," -d ",duration," -o ",outputfile,";",sep="")
    print(command)
    system(command)
    Sys.sleep(1)
    
    ###Immediately re-read outputfile then delete it
    temp <- read.table(outputfile,header=T,stringsAsFactors = F)
    file.remove(outputfile)
    ###Modify temp 
    temp["outside"] <- temp$water+temp$sugar+temp$outzone
    temp <- temp[c("tag","nest","outside","undetected")]; temp["detected"] <- temp$nest+temp$outside
    temp["prop_time_outside"] <- temp$outside/temp$detected
    
    ### add information
    temp <- data.frame(colony=colony,colony_size=colony_size,treatment=treatment,period=period,time_hours=split_info[i,"time_hours"],time_of_day=split_info[i,"time_of_day"],temp,stringsAsFactors = F)
    ### remove ants that did not survive until the end
    temp <- merge(temp,tag[c("tag","final_status","group")]); names(temp)[names(temp)=="group"] <- "status"
    temp <- temp[which(temp$final_status=="alive"),]
    if (!grepl("age",data_path)){
      temp$age <- NA
    }else{
      temp <- merge(temp,colony_ages,all.x=T,all.y=F)
    }
    
    ###
    temp <- temp[c("colony","colony_size","treatment","tag","age","status","period","time_hours","time_of_day","nest","outside","detected","undetected","prop_time_outside")] 
    space_dat     <- rbind(space_dat,temp)
    behaviour_dat <- rbind(behaviour_dat,temp[c("colony","colony_size","treatment","tag","age","status","period","time_hours","time_of_day","prop_time_outside")])
    rm(list="temp")
  }
}
###Write space dat
write.table(space_dat,file=outputspace_file,col.names=T,row.names=F,quote=F,append = F)
if (!grepl("age",data_path)){
  output_behav <- paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
  if (!file.exists(output_behav)){dir.create(output_behav,recursive=T)}
  outputbehav_file <- paste(output_behav,"/individual_behavioural_data.txt",sep="")

  behaviour_dat <- behaviour_dat[c("colony","colony_size","treatment","tag","age","status","period","time_hours","time_of_day","prop_time_outside")]
  write.table(behaviour_dat,file=outputbehav_file,col.names=T,row.names=F,quote=F,append = F)
  
  output_pre_treatment <- paste(data_path,"/processed_data/individual_behaviour/pre_treatment",sep="")
  if (!file.exists(output_pre_treatment)){dir.create(output_pre_treatment,recursive=T)}
  write.table(space_dat[which(space_dat$period=="before"),],file=paste(output_pre_treatment,"/network_position_vs_time_outside.dat",sep=""),col.names=T,row.names=F,quote=F,append = F)
  
}