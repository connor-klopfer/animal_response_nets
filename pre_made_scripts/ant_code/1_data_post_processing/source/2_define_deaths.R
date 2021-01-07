####2_define_deaths.R#####


### Calls C++ functions define_death (https://github.com/laurentkeller/anttrackingUNIL)
### Takes a datfile and a tagfile has input, infers the frame of death/tag loss of each ant, and fills that information in the relevant column of the tagfile

### Note that all death/tag loss times were further checked manually using the videos. All deaths times provided in the subfolder original_data/tag_files are therefore exact times rather than the times inferred using define_death

###Created by Nathalie Stroeymeyt

###################################
####Define directories
working_path <- paste(data_path,"intermediary_analysis_steps/datfiles_uncorrected",sep="/")
setwd(working_path)
datlist <-list.files(pattern="dat")

for (datfile in datlist){
  tagfile <- gsub("dat","tags",datfile)
  logfile <- gsub("dat","log",datfile)
  
  print(paste("finding deaths for",datfile))
  Command <- paste(executables_path,"/define_death ", 
                   datfile, " ",
                   tagfile, " ",
                   logfile,
                   sep=""
  )   
  system(Command)     
}

