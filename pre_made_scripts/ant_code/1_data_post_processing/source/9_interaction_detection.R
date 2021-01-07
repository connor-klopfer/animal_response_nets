######9_interaction_detection.R

### Calls C++ functions interaction_close_front_contacts and filter_interactions_cut_immobile (https://github.com/laurentkeller/anttrackingUNIL)

### interaction_close_front_contacts takes datfile and tagfile as an input and returns the list of all pairs of interacting ants on each frame
### filter_interactions_cut_immobile takes the output from interaction_close_front_contacts, checks whether subsequent interactions between the same pair of ants belong to the same interaction event or 2 distinct events, and returns the list of all interaction events

###Created by Nathalie Stroeymeyt

###########################################################################################

####get dat list #####
input_dat <- paste(data_path,"intermediary_analysis_steps/dat_files",sep="/")
setwd(input_dat)
dat_list <- paste(input_dat,list.files(pattern="dat"),sep="/")

####output folders
output_unfiltered <- paste(data_path,"/intermediary_analysis_steps/unfiltered_interactions",sep=""); if (!file.exists(output_unfiltered)){dir.create(output_unfiltered,recursive=T)}
output_filtered <- paste(data_path,"/intermediary_analysis_steps/filtered_interactions",sep="");  if (!file.exists(output_filtered)){dir.create(output_filtered,recursive=T)}

######### parameter values for detect interactions
distance       <- 60   #####considers only ants who are distant by tl1/2+tl2/2+distance (speeds process up)
delta_angle    <- 0    ####tests points along an arc with angle = -delta_angle,+delta_angle to tag orientation
angle_interval <- 0    ####tests points on that arc every angle_interval degrees
angle_parallel <- 0    #####will use all contacts, even back to front (because all contact can lead to transmission event)

####filter interaction parameters
width_ratio  <- 0.7584637
width_factor <-0.7645475
filter_arguments <- " -t 20 -m 240 -d 25 -a 6 "
#### t: maximum time gap, in frames (10 sec); after this is considered a new interaction
#### m: maximum asleep duration, in frames . After that the interaction is either cut short, or if the ants move again, a new interaction is declared
#### d: maximum movement distance (25 pixels) for an interaction to be considered the same
#### a: asleep maximum movement. If movement less than that between successive frames then ant is considered asleep

for (datfile in dat_list){
  dat_root <- gsub("\\.dat","",unlist(strsplit(datfile,split="/"))[grepl("\\.dat",unlist(strsplit(datfile,split="/")))])
  colony   <- unlist(strsplit(dat_root,split="_"))[grepl("colony",unlist(strsplit(dat_root,split="_")))]
  tagfile  <- tag_list[which(grepl(colony,tag_list))]
  if (length(tagfile)>1){
    tagfile <- tagfile[grepl(unlist(strsplit(dat_root,"_"))[grepl("Treatment",unlist(strsplit(dat_root,"_")))],tagfile)]
  }
  
  ###detect interactions
  command <- paste(executables_path,"/interaction_close_front_contacts ",datfile," ",tagfile," ",output_unfiltered,"/",dat_root,".txt ",distance," ",angle_parallel," ",width_factor," ",width_ratio," ",delta_angle," ",angle_interval,";",sep="")
  print(command)
  system(command)
  
  ###filter interactions
  command <- paste(executables_path,"/filter_interactions_cut_immobile -i ",output_unfiltered,"/",dat_root,".txt -o ",output_filtered,"/",dat_root,".txt",filter_arguments, tagfile,sep="")
  print(command)
  system(command)
}