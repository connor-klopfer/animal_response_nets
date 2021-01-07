####10_process_interaction_files.R#####

#### Removes interactions involving ants that died before the end of the experiment from the interaction file
#### Adds metadata to the interaction file

###Created by Nathalie Stroeymeyt

#################################
options(digits=16) ; options(digits.secs=6) ; options("scipen" = 10)

####get interaction list #####
input_interac <- paste(data_path,"intermediary_analysis_steps/filtered_interactions",sep="/")
setwd(input_interac)
interac_list <- paste(input_interac,dir(),sep="/")

####define outputfolder
if (!grepl("survival",data_path)){
  outputfolder  <- paste(data_path,"/intermediary_analysis_steps/binned_interaction_lists",sep="")
}
outputfolder2 <- paste(data_path,"/intermediary_analysis_steps/full_interaction_lists",sep="")

##################################
for (interac in interac_list){
  root_name <- unlist(strsplit(interac,split="/"))[grepl("colony",unlist(strsplit(interac,split="/")))]
  print(paste(unlist(strsplit(root_name,split="_")),collapse="_"))
  colony    <- unlist(strsplit(root_name,split="_"))[grepl("colony",unlist(strsplit(root_name,split="_")))]
  tagfile   <- tag_list[grepl(colony,tag_list)]
  if (length(tagfile)>1){
    tagfile <- tagfile[grepl(unlist(strsplit(gsub("\\.txt","",root_name),"_"))[grepl("Treatment",unlist(strsplit(gsub("\\.txt","",root_name),"_")))],tagfile)]
  }
  
  outputfolder4 <- paste(outputfolder2,"/",gsub("\\.txt","",unlist(strsplit(root_name,split="_"))[grepl("Treatment",unlist(strsplit(root_name,split="_")))]),"/observed",sep="")
  if(!file.exists(outputfolder4)){dir.create(outputfolder4,recursive=T)}
  
  ### read in interaction file
  tab <- read.table(interac,header=T,stringsAsFactors = F,sep=",",comment.char="%")
  if (length(names(tab)[grepl("X\\.",names(tab))])>0){
    names(tab)[grepl("X\\.",names(tab))] <- gsub("X\\.","",names(tab)[grepl("X\\.",names(tab))])
    
    ####read-in tag list to define list of live ants
    tag        <- read.tag(tagfile)$tag; names(tag)[names(tag)=="#tag"] <- "tag"; tag <- tag[which(tag$tag!="#tag"),]
    if (!grepl("survival",data_path)){
      alive      <- tag[which(tag$final_status=="alive"),"tag"]
    }else{
      alive      <- tag[which(as.numeric(tag$death)==0|as.numeric(tag$death)>=max(tab$Stopframe,na.rm=T)),"tag"]
    }
    
    ###remove interactions involving ants that died before the end of the experiment
    tab <- tab[which( (tab$Tag1%in%alive)&(tab$Tag2%in%alive)),]
    
    ###Add metadata and for each of the pre-treatment and pre-treatment tabectories
    if(!grepl("survival",data_path)){
      splitfile   <- split_list[grepl(gsub("\\.txt","",root_name),split_list)]
      if (length(splitfile)>0){
        outputfolder3 <- paste(outputfolder,"/",gsub("\\.txt","",unlist(strsplit(root_name,split="_"))[grepl("Treatment",unlist(strsplit(root_name,split="_")))]),"/observed",sep="")
        if(!file.exists(outputfolder3)){dir.create(outputfolder3,recursive=T)}
        #### read in aggregation information
        splitinfo <- read.table(splitfile,header=T, stringsAsFactors = F)
        indices             <- unlist(lapply(tab$Starttime,function(x)max(which(x>=splitinfo$time))))
        tab["time_hours"]  <- splitinfo$time_hours[indices]
        tab["time_of_day"] <- splitinfo$time_of_day[indices]
        
      }else{
        tab["time_hours"]  <- 0
        tab["time_of_day"] <- 12
      }
      
    }else{
      tab["time_hours"]  <- 0
      tab["time_of_day"] <- 12
    }
    tab["colony"]      <- colony
    tab["treatment"]   <- info[which(info$colony==as.numeric(gsub("colony","",colony))),"treatment"]
    
    ### Write the interaction list corresponding to each time point into a folder
    if (length(unique(tab$time_hours))>1){
      for (time_hours in unique(tab$time_hours)){
        subtab      <- tab[which(tab$time_hours==time_hours),]
        time_of_day <- unique(subtab$time_of_day)
        outfile <- paste(outputfolder3,"/",gsub("\\.txt","",root_name),"_TH",time_hours,"_TD",time_of_day,"_interactions.txt",sep="")
        write.table(subtab,file=outfile,col.names=T,row.names=F,append=F,quote=F)
      }
    }
    ### Overwrite the full interaction list and copy it in another folder
    write.table(tab,file=interac,col.names=T,row.names=F,append=F,quote=F)
    outfile <- paste(outputfolder4,"/",gsub("\\.txt","",root_name),"_interactions.txt",sep="")
    command <- paste("cp ", interac, " ",outfile,sep="")
    system(command)
    
  }
}