####4_retagged_ant_modifications.R#####

### Calls C++ functions datcorr created when compiling Antorient (https://github.com/laurentkeller/anttrackingUNIL)
### Calls C++ functions change_tagid 

### Modifies dat files for treated ants that lost their tags and were given a new tag during treatment (modify tag id)
### Modifies dat files for treated ants that lost their tags and were retagged with the same tag during treatment (modify rotation)


###Created by Nathalie Stroeymeyt

####Navigate to tag folder and list tagfiles
tag_path    <- paste(data_path,"original_data/tag_files",sep="/")
setwd(tag_path)
tag_list <- paste(tag_path,list.files(pattern="\\.tags"),sep="/")

####Get information on retagged ants
setwd(paste(data_path,"original_data/retagged_ants",sep="/"))
retagged <- read.table("retagged_ants.txt",header = T,stringsAsFactors = F)

###1. Deal with ants that were given a new tag
newly_tagged <- retagged[which(retagged$comment%in%c("new_tag","old_tag")),]

###for each relevant replicate
for (colony_number in unique(newly_tagged$colony)){
  ###Get into the retagged ant folder, which contains the necessary input files to make the changes
  setwd(paste(data_path,"original_data/retagged_ants",sep="/"))
  colony <- paste("colony",paste(rep(0,3-nchar(colony_number)),collapse=""),colony_number,sep="")

  ###Get the tag file that contains information on the tag that still need to be rotated
  tag_files <- paste(data_path,"original_data/retagged_ants",list.files()[grepl(colony,list.files())],sep="/")
  sub_tag_file <- tag_files[which(grepl("fororient",tag_files))]
  tag_file <- tag_files[which(!grepl("fororient",tag_files))]
  ###Get further information on old and new tag file
  box <- info[info$colony==colony_number,"box"]
  post_treatment_tag <- newly_tagged[newly_tagged$comment=="new_tag","tag"]
  pre_treatment_tag <- newly_tagged[newly_tagged$comment=="old_tag","tag"]

  ###first, export taglist, which will be used as an argument for changeid
  ####if file exists delete it
    if (file.exists("taglist.txt")){
      file.remove("taglist.txt")
    }
  ###NB: the tag that appears in the final tag file is the post-treatment tag, as it is the one which corresponds to the qPCR data
  write.table(data.frame(from=pre_treatment_tag,to=post_treatment_tag),file="taglist.txt",append=F,row.names=F,col.names=F,sep=",")

  ###second, rewrite preTreatment dat file
  setwd(paste(data_path,"/intermediary_analysis_steps/dat_files",sep="/"))
  dat_list <- paste(data_path,"/intermediary_analysis_steps/dat_files",list.files(pattern="\\.dat"),sep="/")
  dat_file <- dat_list[which(grepl(colony,dat_list)&grepl("PreTreatment",dat_list))]
  name_root <- gsub("\\.dat","",unlist(strsplit(dat_file,split="\\/"))[grepl("\\.dat",unlist(strsplit(dat_file,split="\\/")))])
  log_file <- paste(paste(data_path,"intermediary_analysis_steps/datfiles_uncorrected",sep="/"),"/",name_root,".log",sep="")

  ###apply rotation for the tag contained in sub_tag_file.
  command_line <- paste(paste(executables_path,"/datcorr",sep=""),dat_file,sub_tag_file,gsub("\\.dat","_additional_rotation.dat",dat_file),log_file,sep=" ")
  print(command_line)
  system(command_line)

  ###change file names
  system(paste("mv ",dat_file," ",gsub("\\.dat","_original.dat",dat_file)))
  system(paste("mv ",gsub("\\.dat","_additional_rotation.dat",dat_file)," ",dat_file))

  ###and now modify tag id in the datfile
  command_line <- paste(executables_path,"/change_tagid -d ",dat_file," -t ",tag_file," -o " ,gsub("\\.dat","_old_tagid_changed_to_new_tagid.dat",dat_file)," -b ",box," -l ",paste(data_path,"original_data/retagged_ants","taglist.txt",sep="/"),sep="")
  print(command_line)
  system(command_line)

  system(paste("rm ",dat_file))
  system(paste("mv ",gsub("\\.dat","_old_tagid_changed_to_new_tagid.dat",dat_file)," ",dat_file))
  system(paste("rm ",gsub("\\.dat","_original.dat",dat_file)))
}

###2. Deal with ants that were given the same tag but with a different orientation
retagged <- retagged[which(retagged$comment%in%c("retagged")),]
###for each relevant replicate
for (colony_number in unique(retagged$colony)){
  ###Get into the retagged ant folder, which contains the necessary input files to make the changes
  setwd(paste(data_path,"original_data/retagged_ants",sep="/"))
  colony <- paste("colony",paste(rep(0,3-nchar(colony_number)),collapse=""),colony_number,sep="")

  ###Get the tag file that contains information on the tag that still need to be rotated
  sub_tag_file <- paste(data_path,"original_data/retagged_ants",list.files()[grepl(colony,list.files())],sep="/")

  box <- info[info$colony==colony_number,"box"]

  ###second, rewrite postTreatment dat file
  setwd(paste(data_path,"/intermediary_analysis_steps/dat_files",sep="/"))
  dat_list <- paste(data_path,"/intermediary_analysis_steps/dat_files",list.files(pattern="\\.dat"),sep="/")
  dat_file <- dat_list[which(grepl(colony,dat_list)&grepl("PostTreatment",dat_list))]
  name_root <- gsub("\\.dat","",unlist(strsplit(dat_file,split="\\/"))[grepl("\\.dat",unlist(strsplit(dat_file,split="\\/")))])
  log_file <- paste(paste(data_path,"intermediary_analysis_steps/datfiles_uncorrected",sep="/"),"/",name_root,".log",sep="")

  ###apply rotation for the tag contained in sub_tag_file.
  command_line <- paste(paste(executables_path,"/datcorr",sep=""),dat_file,sub_tag_file,gsub("\\.dat","_additional_rotation.dat",dat_file),log_file,sep=" ")
  print(command_line)
  system(command_line)

  ###change file names
  system(paste("mv ",dat_file," ",gsub("\\.dat","_original.dat",dat_file)))
  system(paste("mv ",gsub("\\.dat","_additional_rotation.dat",dat_file)," ",dat_file))
  system(paste("rm ",gsub("\\.dat","_original.dat",dat_file)))
}
