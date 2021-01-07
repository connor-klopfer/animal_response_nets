
clean <- function(){
  rm(list=ls(envir = .GlobalEnv)[!ls(envir = .GlobalEnv)%in%to_keep], envir = .GlobalEnv)
  no_print <- gc(verbose=F)
  Sys.sleep(1)
}

closest_match <- function(x,y){
  return(min(which(abs(x-y)==min(abs(x-y))),na.rm=T))
}
read.tag <- function(tagfile){
  tag <- read.table(tagfile,sep=",",comment.char="%",as.is=TRUE,fill=TRUE,stringsAsFactors=F)
  ##find out in which line the names are contained
  index_names <- match("#tag",tag[,1])
  ##remove the stuff above
  if (index_names > 1){
    header_part <- tag[1:(index_names-1),]
    tag <- data.frame(tag[index_names:nrow(tag),])
  }
  ##in case there was a problem with the reading, fix it
  if (ncol(tag)==1){
    ncols <- min(which(!is.na(as.numeric(as.character(tag[,1]))))) -1
    new_tag <- {}
    for (line in 1:(nrow(tag)/ncols)){
      new_tag <- rbind(new_tag,as.character(tag[(((line - 1) *ncols) +1 ):(((line - 1) *ncols) + ncols ),]))
    }
    tag <- data.frame(new_tag,stringsAsFactors=FALSE)
  }
  
  ##update the line in which the names are contained
  index_names <- match("#tag",tag[,1])
  ##get name list
  original_name_list <- as.character(tag[index_names,])
  names(tag) <- original_name_list
  return(list(tag=tag,header_part=header_part))
}

high_threshold <- 0.0411
info           <- read.table(paste(data_path,"/original_data/info.txt",sep=""),header=T,stringsAsFactors = F)
treated        <- read.table(paste(data_path,"/original_data/treated_worker_list.txt",sep=""),header=T,stringsAsFactors = F)
task_groups    <- read.table(paste(data_path,"original_data/task_groups.txt",sep="/"),header=T,stringsAsFactors = F)
if (grepl("age",data_path)){
  ages      <- read.table(paste(data_path,"original_data","ant_ages.txt",sep="/"),header=T,stringsAsFactors = F)
}

if(file.exists(paste(data_path,"/original_data/info_dat.txt",sep=""))){
  info_datfile   <- read.table(paste(data_path,"/original_data/info_dat.txt",sep=""),header=T,stringsAsFactors = F)
  info_plume     <- read.table(paste(data_path,"/original_data/info_plume.txt",sep=""),header=T,stringsAsFactors = F)
}

#### get time_aggregation_list
if (!grepl("survival",data_path)){
  input_aggregation_info        <- paste(data_path,"original_data/time_aggregation_info",sep="/")
  setwd(input_aggregation_info)
  split_list                    <- paste(input_aggregation_info,list.files(pattern="txt"),sep="/")
}

####get tag list #####
input_tag <- paste(data_path,"original_data/tag_files",sep="/")
setwd(input_tag)
tag_list <- paste(input_tag,list.files(pattern="tags"),sep="/")


queenid <- 665

