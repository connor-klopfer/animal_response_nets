####1_apply_rotation_to_datfiles.R#####

### Calls C++ functions datcorr created when compiling anttrackingUNIL/Antorient (https://github.com/laurentkeller/anttrackingUNIL)

### Takes tracking datfile and tagfile as an input and returns an oriented datfile

### Note that the tag files used as an input must have been oriented using Antorient (https://github.com/laurentkeller/anttrackingUNIL),
### had false detections removed, 
### and include the frame of death / tag loss information
### The final processed tag files are provided in the subfolder original_data/tag_files and will be used in all later steps of the analysis

###Created by Nathalie Stroeymeyt


###################################
### define directories
tag_path    <- paste(data_path,"original_data/tag_files",sep="/")
input_path  <- paste(data_path,"intermediary_analysis_steps/datfiles_uncorrected",sep="/")
output_path <- paste(data_path,"intermediary_analysis_steps/dat_files",sep="/")

####If necessary, create output folder
if(!exists(output_path)){
  dir.create(output_path,recursive = T,showWarnings = F)
}

####Navigate to input folder and list datfiles
setwd(input_path)
dat_list <- paste(input_path,list.files(pattern="\\.dat"),sep="/")

####Navigate to tag folder and list tagfiles
setwd(tag_path)
tag_list <- paste(tag_path,list.files(pattern="\\.tags"),sep="/")

for (dat_file in dat_list){
  name_root <- gsub("\\.dat","",unlist(strsplit(dat_file,split="\\/"))[grepl("\\.dat",unlist(strsplit(dat_file,split="\\/")))])
  tag_file <- paste(unlist(strsplit(name_root,"_"))[!grepl("Treatment",unlist(strsplit(name_root,"_")))],collapse="_")
  tag_file <- tag_list[grepl(tag_file,tag_list)]
  if (length(tag_file)>1){
    tag_file <- tag_file[grepl(unlist(strsplit(name_root,"_"))[grepl("Treatment",unlist(strsplit(name_root,"_")))],tag_file)]
  }
  new_dat_file <- paste(output_path,"/",name_root,".dat",sep="")
  log_file <- paste(input_path,"/",name_root,".log",sep="")
 
  command_line <- paste(paste(executables_path,"/datcorr",sep=""),dat_file,tag_file,new_dat_file,log_file,sep=" ")
  print(command_line)
  system(command_line)
}