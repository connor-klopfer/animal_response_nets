####1_trackconverter.R#####

### Calls C++ functions trackconverter and controldat (https://github.com/laurentkeller/anttrackingUNIL)

### Takes tracking csv files as an input and creates a binary data file of class 'datfile', a log file, a tags file containing the list of dead, and a double file containing information on double detections
### Checks the integrity of the file created

### Note that the tag files produced by trackconverter need to be modified to take into account the orientation of the ant compared to the tag using Antorient (https://github.com/laurentkeller/anttrackingUNIL),
### to remove false detections, 
### and to record the frame of death / tag loss if relevant (see 2_define_deaths.R)
### The final processed tag files are provided in the subfolder original_data/tag_files and will be used in all later steps of the analysis

###Created by Nathalie Stroeymeyt

###################################
### define directories
input_path  <- paste(data_path,"original_data",sep="/")
output_path <- paste(data_path,"intermediary_analysis_steps/datfiles_uncorrected",sep="/")

####If necessary, create output folder
if(!exists(output_path)){
  dir.create(output_path,recursive = T,showWarnings = F)
}

#####Navigate to tracking data folder and list data files
setwd(paste(input_path,"tracking",sep="/"))
tracking_files <- paste(input_path,"/tracking/",list.files(pattern="csv"),sep="")
tracking_files <- tracking_files[!grepl("corr",tracking_files)]

#####Navigate to output folder
setwd(output_path)

#####for each file
for (tracking_file in tracking_files){
  file_root_name <- gsub(paste(input_path,"/tracking/",sep=""),"",tracking_file)
  colony <- as.numeric(gsub("colony","",unlist(strsplit(file_root_name,split="_"))[1]))
  command_line <- paste(executables_path,"/trackconverter ",gsub("csv","dat",file_root_name), " ",info[info$colony==colony,"box"]," ",gsub("csv","log",file_root_name)," ", tracking_file,sep="")
  print(command_line)
  system(command_line)
  
  command_line <- paste(executables_path,"/controldat ",gsub("csv","dat",file_root_name),sep="")
  print(command_line)
  system(command_line)
}

#####Navigate to tracking data folder and remove temporary corrected files
setwd(paste(input_path,"tracking",sep="/"))
system("rm *corr*")
