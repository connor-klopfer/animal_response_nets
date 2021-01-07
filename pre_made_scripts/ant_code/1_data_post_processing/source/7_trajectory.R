####7_trajectory.R#####

### Calls C++ functions trajectory (https://github.com/laurentkeller/anttrackingUNIL)

### Takes a datfile and a tagfile as inputs, and returns the trajectory of each ant in a separate file

#################################

####get dat list #####
input_dat <- paste(data_path,"intermediary_analysis_steps/dat_files_nest_zone",sep="/")
setwd(input_dat)
dat_list <- paste(input_dat,list.files(pattern="dat"),sep="/")

#### define outputfolders
outputfolder <- paste(data_path,"/intermediary_analysis_steps/trajectories",sep="")
if (!file.exists(outputfolder)){dir.create(outputfolder,recursive=T)}

for (datfile in dat_list){
  root_name   <- gsub("\\.dat","",unlist(strsplit(datfile,split="/"))[grepl("colony",unlist(strsplit(datfile,split="/")))])
  outputfolder2 <- paste(outputfolder,root_name,sep="/");if (!file.exists(outputfolder2)){dir.create(outputfolder2,recursive=T)}
  
  colony      <- unlist(strsplit(root_name,split="_"))[grepl("colony",unlist(strsplit(root_name,split="_")))]
  tagfile     <- tag_list[grepl(colony,tag_list)]
  
  command <- paste(executables_path,"/trajectory ",datfile," ",tagfile," ",outputfolder2,"/ant", ";",sep="")
  print(command)
  system(command)
}
