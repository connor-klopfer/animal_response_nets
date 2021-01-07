####15_heatmaps_individual.R#####

### Calls C++ function heatmap3_tofile (https://github.com/laurentkeller/anttrackingUNIL)
### Takes a datfile, a tagfile, an ant, and a time period as inputs, and returns a file containing the number of visits made by the specified ant to each coordinate

###Created by Nathalie Stroeymeyt###

#################################
options(digits=16) ; options(digits.secs=6) ; options("scipen" = 10)

####get dat list #####
input_dat <- paste(data_path,"intermediary_analysis_steps/dat_files",sep="/")
setwd(input_dat)
dat_list <- paste(input_dat,list.files(pattern="dat"),sep="/")

#### define outputfolders
outputfolder <- paste(data_path,"/intermediary_analysis_steps/heatmaps/individual",sep="")
if (!file.exists(outputfolder)){dir.create(outputfolder,recursive=T)}

for (datfile in dat_list){
  root_name   <- gsub("\\.dat","",unlist(strsplit(datfile,split="/"))[grepl("colony",unlist(strsplit(datfile,split="/")))])
  outputfolder2 <- paste(outputfolder,root_name,sep="/");if (!file.exists(outputfolder2)){dir.create(outputfolder2,recursive=T)}

  colony      <- unlist(strsplit(root_name,split="_"))[grepl("colony",unlist(strsplit(root_name,split="_")))]
  tagfile     <- tag_list[grepl(colony,tag_list)]
  split_file  <- split_list[grepl(root_name,split_list)]
  
  split_info  <- read.table(split_file,header=T,stringsAsFactors = F)
  
  tag   <- read.tag(tagfile)$tag
  names(tag)[names(tag)=="#tag"] <- "tag"
  tag <- tag[which(tag$tag!="#tag"),]
  tag[which(tag$age==0),"age"] <- NA
  antlist <- tag[which(tag$final_status=="alive"),"tag"]
  
  box         <- info[which(info$colony==as.numeric(gsub("colony","",colony))),"box"]
  treatment   <- info[which(info$colony==as.numeric(gsub("colony","",colony))),"treatment"]
  colony_size <- info[which(info$colony==as.numeric(gsub("colony","",colony))),"colony_size"]
  final_time <- info_datfile[which(info_datfile$dat_file==paste(root_name,".dat",sep="")),"end_time"]
  if (grepl("PreTreatment",root_name)){
    period <- "before"
  }else {
    period <- "after"
  }
  
  for (i in 1:nrow(split_info)){
    print(paste(datfile,"; Time window =",i))
    start_time <- split_info [i,"time"]
    if (i<nrow(split_info)){
      duration  <- split_info [i+1,"time"] - start_time -0.25
    }else{
      duration  <- final_time-start_time-0.25
    }
    ###Run Heatmap
    for (ant in antlist){
      outputfile <- paste(outputfolder2,"/","ant",paste(rep(0,3-nchar(ant)),collapse=""),ant,"_TD",split_info[i,"time_of_day"],"_TH",split_info[i,"time_hours"],sep="")
      if (!file.exists(paste(outputfile,"_ant",ant,"1.txt",sep=""))){
        command <- paste(executables_path,"/heatmap3_tofile -i ",datfile," -t ",tagfile," -b ",box," -s ",start_time, " -d ", duration," -n 1 -a ",ant," -m ",outputfile,";",sep="")
        system(command)
      }
    }  
  }
}
