####8_process_trajectory_files.R#####

#### Defines activity bouts and calculates summary statistics for each trajectory 
#### Outputs: modified trajectory file and modifed individual behaviour file

#################################
options(digits=16) ; options(digits.secs=6) ; options("scipen" = 10)

#### Parameters
to_keep_ori <- to_keep

max_time <- 1.5 # max time in sec, for which it is acceptable to calculate a distance moved (e.g. 1.5 sec: allows for up to 2 missing points between two successive detections and still calculates distance moved)

###  Functions
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

####get trajectory list #####
input_traj <- paste(data_path,"intermediary_analysis_steps/trajectories",sep="/")
setwd(input_traj)
traj_list <- paste(input_traj,dir(),sep="/")
#### reduce traj list to only PreTreatment folders
#### for each ant we will successively consider the Pre-treatment, then the Post-treatment trajectory
traj_list <- traj_list[grepl("Pre",traj_list)]

TAB_ALL <- NULL

##################################
for (traj_folder in traj_list){
  root_name <- unlist(strsplit(traj_folder,split="/"))[grepl("colony",unlist(strsplit(traj_folder,split="/")))]
  print(paste(unlist(strsplit(root_name,split="_"))[!grepl("Treatment",unlist(strsplit(root_name,split="_")))],collapse="_"))
  colony    <- unlist(strsplit(root_name,split="_"))[grepl("colony",unlist(strsplit(root_name,split="_")))]
  tagfile   <- tag_list[grepl(colony,tag_list)]
  splitfiles   <- split_list[grepl(colony,split_list)]
  
  ####read-in tag list to define list of live ants
  tag        <- read.tag(tagfile)$tag; names(tag)[names(tag)=="#tag"] <- "tag"; tag <- tag[which(tag$tag!="#tag"),]
  alive      <- paste("ant_",tag[which(tag$final_status=="alive"),"tag"],".txt",sep="")
  
  splitinfo_Pre <- read.table(splitfiles[grepl("Pre",splitfiles)],header=T, stringsAsFactors = F)
  splitinfo_Post <- read.table(splitfiles[grepl("Post",splitfiles)],header=T, stringsAsFactors = F)
  
  ###Navigate in folder
  setwd(traj_folder)
  traj_files <- list.files()
  
  to_keep <- unique(c(to_keep,ls(),"traj_file","to_keep_ori"))
  ###Process each trajectory file
  for (traj_file in traj_files){
    print(gsub("_"," ",gsub("\\.txt","",traj_file)))
    traj_file_Pre  <- paste(traj_folder,traj_file,sep="/")
    traj_file_Post <- gsub("Pre","Post",traj_file_Pre)
    ####delete files corresponding to dead ants
    if (!traj_file %in% alive){
      print("Deleting file (ant died before the end)")
      file.remove(traj_file_Pre)
      file.remove(traj_file_Post)
    }else{# if (!traj_file %in% alive)
      print("Modifying trajectory...")
      ###check if analysis has already been done
      output <- system(paste( "head -1 ",traj_file_Post),intern=T)
      if (!grepl("bout_id",output)){
        bout_counter <- 0
        ####do the analysis for each of the Pre- and Post- files 
        for (suffix in c("_Pre","_Post")){
          
          output <- system(paste( "head -1 ",get(paste("traj_file",suffix,sep=""))),intern=T)
          if (!grepl("bout_id",output)){
            traj <- read.table(get(paste("traj_file",suffix,sep="")),header=T,sep=",",comment.char="%")
            names(traj)[names(traj)=="X.frame"] <- "frame"
            #####Remove occasional duplicates
            traj <- traj[!duplicated(traj$frame),]
            ##########################################################      
            ######1. calculate distance moved ########################
            ##########################################################
            ####Prepare two tables offset by one line to facilitate calculation for calculations
            traj_1 <- insertRow(traj,rep(NA,ncol(traj)),1)
            traj_2 <- insertRow(traj,rep(NA,ncol(traj)),(nrow(traj)+1))[1:(nrow(traj)+1),]
            names(traj_1) <- paste(names(traj_1),"_previousframe",sep="")
            ####merge them
            full_trajectory <- data.frame(traj_1,traj_2)
            ####calculate time difference between two successive coordinates
            full_trajectory["time_diff"] <- full_trajectory$time-full_trajectory$time_previousframe
            ####round time_diff to nearest 0.5
            full_trajectory$time_diff <- round(full_trajectory$time_diff/0.5)*0.5;full_trajectory$time_diff[full_trajectory$time_diff==0]<-0.5
            #### calculate distance moved between two successive frames
            full_trajectory[,"total_dist"] <- (sqrt((full_trajectory[,"xcoor"]-full_trajectory[,"xcoor_previousframe"])^2+(full_trajectory[,"ycoor"]-full_trajectory[,"ycoor_previousframe"])^2))
            #### calculate a distance per frame (i.e. speed) for intervals that are less than max_time 
            full_trajectory[which(full_trajectory$time_diff<=max_time),"dist_per_frame"] <- (sqrt((full_trajectory[which(full_trajectory$time_diff<=max_time),"xcoor"]-full_trajectory[which(full_trajectory$time_diff<=max_time),"xcoor_previousframe"])^2+(full_trajectory[which(full_trajectory$time_diff<=max_time),"ycoor"]-full_trajectory[which(full_trajectory$time_diff<=max_time),"ycoor_previousframe"])^2))/((full_trajectory[which(full_trajectory$time_diff<=max_time),"time_diff"]/0.5))
            #### round
            full_trajectory$dist_per_frame <- round(100*full_trajectory$dist_per_frame)/100
            full_trajectory$total_dist <- round(100*full_trajectory$total_dist)/100
            #### now use this new data to create a new traj file
            traj <- full_trajectory[(1:(nrow(full_trajectory)-1)),c("frame","time","box","xcoor","ycoor","dist_per_frame","total_dist")]
            
            #####calculate turn angles
            trajectory <- as.ltraj(traj[c("xcoor","ycoor")],date=as.POSIXct(traj$time,origin="1970-01-01"),id=gsub("\\.txt","",gsub("ant_","",traj_file)),typeII=T,slsp="missing")
            turn_angles <- abs(trajectory[[1]]$rel.angle)
            #####modify it slightly, because I defined distance as distance moved between current frame and the previous; whereas turn angles are defined between current frame and the next
            turn_angles <- c(NA,turn_angles[1:(length(turn_angles)-1)])
            traj["turn_angle"] <- turn_angles
            #####no movement = no turn. Therefore fill in NA turn angles when dist=0
            traj[which(is.na(traj$turn_angle)&traj$dist_per_frame==0),"turn_angle"] <- 0
            #####whenever the time separating 2 successive detections is greater than max_time, enter NA; because then the turn angle is not meaningful
            times_1 <- c(NA,traj$time);times_2 <- c(traj$time,NA);time_diff <- times_2-times_1;time_diff <- time_diff[1:(length(time_diff)-1)]
            time_diff <- round(time_diff/0.5)*0.5
            traj[which(time_diff>max_time),"turn_angle"] <- NA
          }else{
            traj <- read.table(get(paste("traj_file",suffix,sep="")),header=T,comment.char="%")
            traj <- traj[,which(!names(traj)%in%c("type","bout_id"))]
          }
          ####################################################################    
          #######2. cut trajectory into bouts of activity vs. inactivity #####
          ####################################################################      
          #####define paremeters
          min_segment_duration_sec <- 120; Lmin <- 2*min_segment_duration_sec #minimum duration of a bout, in frames
          max_gap <- 15 ##anything larger than that needs to be subdivided into two different bouts
          dist_threshold <- 15 #everything below that threshold (half a tag length) is considered inactive
          power <- 0.125 ###power used before clustering
          
          
          #####define an object that does not contain NAs in the dist_per_frame column
          traj_noNA <- traj[which(!is.na(traj$dist_per_frame)),]
          
          #####the analysis is only possible if there is enough data in traj_noNA
          if (nrow(traj_noNA)>2*Lmin){
            #####find segments based on distance moved
            #####apply distance threshold
            traj_noNA[traj_noNA$dist_per_frame<= dist_threshold,"dist_per_frame"] <- 0
            
            #####use function cpt.meanvar to find breakpoints in the trajectory data
            if (length(unique(traj_noNA$dist_per_frame))>2){
              segments_dist <- cpt.meanvar(traj_noNA$dist_per_frame, method = "PELT",minseglen=Lmin)#####best by far
              breakpoints <- segments_dist@cpts
              rm(list="segments_dist")
            }else{
              breakpoints <- nrow(traj_noNA)
            }
            
            #####use breakpoints to define start and end times for bouts, and add the information into thr traj file 
            breakpoints <- match(traj_noNA[breakpoints,"time"],traj$time)
            breakpoints[length(breakpoints)] <- max(nrow(traj),breakpoints[length(breakpoints)])
            breakpoints <- c(1,breakpoints)
            first_index <- breakpoints[1]
            bout_start_indices <- c(first_index,(breakpoints+1)[2:(length(breakpoints)-1)])
            bout_end_indices <- breakpoints[2:(length(breakpoints))]
            bout_indices <- data.frame(start=bout_start_indices,end=bout_end_indices)
            bout_count_index <- NULL
            for (i in 1:nrow(bout_indices)){
              bout_count_index <- c(  bout_count_index , rep(i,(bout_indices[i,"end"]-bout_indices[i,"start"]+1)))
            }
            bout_count_index <- bout_counter + bout_count_index
            traj["bout_index"] <- bout_count_index
            bout_counter <- max(bout_count_index,na.rm=T)
          }#if (nrow(traj_noNA)>20)
          assign(paste("traj",suffix,sep=""),traj)
          rm(list=c("traj"))
        }##for (suffix in c("_Pre","_Post"))
        ###Now combine the two traj objects into a single one for the meta-analysis of active vs. inactive
        if (ncol(traj_Pre)==ncol(traj_Post)){###do it only if bouts were defined for both periods
          trajectory_table <- rbind(data.frame(period="before",traj_Pre),data.frame(period="after",traj_Post))
          ###Perform cluster analysis on the whole table to determine active/inactive using the mean and standard deviation of speed (dist_per_frame) and turn_angle as an input
          for_multi <- aggregate(na.action="na.pass",cbind(turn_angle,(dist_per_frame)^power)~bout_index,function(x)cbind(mean(x,na.rm=T),sd(x,na.rm=T)),data=trajectory_table)
          for_multi <- data.frame(for_multi$bout_index,for_multi$turn_angle,for_multi$V2);names(for_multi) <- c("bout_index","turn_angle","turn_angle_SD","Dist","Dist_SD")
          cluster_BOTH <- kmeans( na.omit( for_multi   [c("turn_angle","Dist","turn_angle_SD","Dist_SD")]),  centers=2) ## clustering on ave & sd of relative turn angle & average distance
          ####Distinguish between active and inactive bout using the speed (active corresponds to higher speed)
          types_BOTH <- c("inactive","active")[order(data.frame(cluster_BOTH["centers"])$centers.Dist)]
          trajectory_table$type <- types_BOTH[cluster_BOTH$cluster[trajectory_table$bout_index]]
          
          ####Use the results from the cluster analysis to define bouts and interbouts
          if (length(unique(trajectory_table$type))==1){
            if (unique(trajectory_table$type)=="active"){to_fill <- "bout1"}else{to_fill <- "interbout1"}
            trajectory_table[1:nrow(trajectory_table),"bout_id"] <- to_fill
            
          }else{
            ##find indices of changes in activity type, in order to pool successive bouts of the same type together
            changes <- data.frame(type=c(NA,trajectory_table$type),type.2=c(trajectory_table$type,NA))
            changes$change <- changes$type==changes$type.2
            breakpoints <- which(!changes$change[1:(length(changes$change)-1)])
            
            start_indices <- c(1,breakpoints)
            end_indices <- c(breakpoints-1,nrow(trajectory_table))
            
            if(trajectory_table[1,"type"]!="active"){
              interbout_indices <- 1+2*c(0:(ceiling(length(start_indices)/2)-1))
              bout_indices <- 2*c(1:floor(length(start_indices)/2))
            }else{
              bout_indices <- 1+2*c(0:(ceiling(length(start_indices)/2)-1))
              interbout_indices <- 2*c(1:floor(length(start_indices)/2)) 
            }
            
            bout_indices <- data.frame(start=start_indices[bout_indices],end=end_indices[bout_indices])
            interbout_indices <- data.frame(start=start_indices[interbout_indices],end=end_indices[interbout_indices])
            
            ###copy the new bout information into the trajectory_table object
            bout_index_list <- NULL
            bout_count_index <- NULL
            for (i in 1:nrow(bout_indices )){
              bout_index_list <- c(bout_index_list, seq(bout_indices[i,"start"],bout_indices[i,"end"]))
              bout_count_index <- c(  bout_count_index , rep(i,(bout_indices[i,"end"]-bout_indices[i,"start"]+1)))
            }
            interbout_index_list <- NULL
            interbout_count_index <- NULL
            for (i in 1:nrow(interbout_indices )){
              interbout_index_list <- c(interbout_index_list, seq(interbout_indices[i,"start"],interbout_indices[i,"end"]))
              interbout_count_index <- c(  interbout_count_index , rep(i,(interbout_indices[i,"end"]-interbout_indices[i,"start"]+1)))
            }
            
            trajectory_table[bout_index_list,"bout_id"] <- paste("bout",bout_count_index,sep="")
            trajectory_table[interbout_index_list,"bout_id"] <- paste("interbout",interbout_count_index,sep="")
          }
          ####Finally, find large gaps in trajectory and subdivide each bout that was defined across the gap
          ####separate bout_id into root and index
          temp <- gregexpr("[0-9]+", trajectory_table$bout_id)
          trajectory_table["idx"] <- as.numeric(unlist(regmatches(trajectory_table$bout_id,temp)))
          trajectory_table[paste("bout",trajectory_table$idx,sep="")==trajectory_table$bout_id,"root"] <- "bout"
          trajectory_table[paste("interbout",trajectory_table$idx,sep="")==trajectory_table$bout_id,"root"] <- "interbout"
          
          ###get times where more than 15 min gap
          times_1 <- c(NA,trajectory_table$time);times_2 <- c(trajectory_table$time,NA);time_diff <- times_2-times_1;time_diff <- time_diff[1:(length(time_diff)-1)]
          indices <- which(time_diff>(60*max_gap))
          for (index in indices){
            trajectory_table[(index:nrow(trajectory_table)),][trajectory_table[(index:nrow(trajectory_table)),"root"]==trajectory_table[index,"root"],"idx"] <-       trajectory_table[(index:nrow(trajectory_table)),][trajectory_table[(index:nrow(trajectory_table)),"root"]==trajectory_table[index,"root"],"idx"] + 1
          }
          
          ###copy new bout id into the table
          trajectory_table <- within(trajectory_table, bout_id <- paste(root,idx,sep=""))
          
          ###Finally, divide the trajectory:table into pre- and post-treatment trajectories
          traj_Pre <-trajectory_table[which(trajectory_table$period=="before"),]
          traj_Post <-trajectory_table[which(trajectory_table$period=="after"),]
        }else{
          ###In case one of the periods did not have enough data for bouts to be defined, fill the columns with NA 
          traj_Pre  <- data.frame(period="before",traj_Pre); traj_Pre["bout_id"] <- NA; traj_Pre["type"] <- NA
          traj_Post <- data.frame(period="after",traj_Post); traj_Post["bout_id"] <- NA; traj_Post["type"] <- NA
        }
        
        for (suffix in c("_Pre","_Post")){
          ###Add metadata and for each of the pre-treatment and pre-treatment trajectories
          splitinfo           <- get(paste("splitinfo",suffix,sep=""))
          traj                <- get(paste("traj",suffix,sep=""))
          indices             <- unlist(lapply(traj$time,function(x)max(which(x>=splitinfo$time))))
          traj["time_hours"]  <- splitinfo$time_hours[indices]
          traj["time_of_day"] <- splitinfo$time_of_day[indices]
          traj["colony"]      <- colony
          traj["treatment"]   <- info[which(info$colony==as.numeric(gsub("colony","",colony))),"treatment"]
          traj                <- traj[c("colony","treatment","box","period","time_hours","time_of_day","time","frame","xcoor","ycoor","dist_per_frame","total_dist","turn_angle","type","bout_id")]
          
          ##################################
          ####### Now analyse trajectory ###
          ##################################
          ####prepare output
          tab <- expand.grid(colony=colony,
                             tag=gsub("\\.txt","",gsub("ant_","",traj_file)),
                             time_hours=splitinfo$time_hours,
                             proportion_time_active=NA,
                             average_bout_speed_pixpersec=NA,
                             total_distance_travelled_pix=NA)
          ###Analyse separately for each time point
          for (time_point in unique(tab$time_hours)){
            subtraj <- traj[which(traj$time_hours==time_point),]
            if (!all(is.na(subtraj$type))){
              bouts <- subtraj[which(subtraj$type=="active"),]
              interbouts <- subtraj[which(subtraj$type=="inactive"),]
              
              bout_number <- length(unique(bouts$bout_id))
              interbout_number <- length(unique(interbouts$bout_id))
              if (!(bout_number==0&interbout_number==0)){###if there is no bout data, leave everything as NA
                if (bout_number==0){###if ant completely inactive
                  tab[which(tab$time_hours==time_point),c("proportion_time_active","total_distance_travelled_pix")] <- 0
                }else {###if ant shows at least some activity
                  tab[which(tab$time_hours==time_point),"total_distance_travelled_pix"] <- sum(bouts$total_dist,na.rm=T)
                  bout_speed <- aggregate(na.rm=T,na.action="na.pass",dist_per_frame~bout_id,FUN=mean,data=bouts)
                  tab[which(tab$time_hours==time_point),"average_bout_speed_pixpersec"] <- mean(bout_speed$dist_per_frame,na.rm = T)*2 ####multiply by 2 because each frame lasts 0.5sec
                  if (interbout_number==0){###if ant completely active
                    tab[which(tab$time_hours==time_point),"proportion_time_active"] <- 1
                  }else{###if ant shows both activity and inactivity
                    ###calculate cumulated bout duration
                    bout_starts <- aggregate(na.rm=T,na.action="na.pass",time~bout_id+time_hours,FUN=min,data=bouts); names(bout_starts)[names(bout_starts)=="time"] <- "Starttime"
                    bout_ends <- aggregate(na.rm=T,na.action="na.pass",time~bout_id+time_hours,FUN=max,data=bouts); names(bout_ends)[names(bout_ends)=="time"] <- "Stoptime"
                    bout_durations <- merge(bout_starts,bout_ends)
                    bout_duration <-  sum(bout_durations$Stoptime-bout_durations$Starttime+0.5,na.rm=T)
                    ###calculate cumulated interbout duration                    
                    interbout_starts <- aggregate(na.rm=T,na.action="na.pass",time~bout_id+time_hours,FUN=min,data=interbouts); names(interbout_starts)[names(interbout_starts)=="time"] <- "Starttime"
                    interbout_ends <- aggregate(na.rm=T,na.action="na.pass",time~bout_id+time_hours,FUN=max,data=interbouts); names(interbout_ends)[names(interbout_ends)=="time"] <- "Stoptime"
                    interbout_durations <- merge(interbout_starts,interbout_ends)
                    interbout_duration <-  sum(interbout_durations$Stoptime-interbout_durations$Starttime+0.5,na.rm=T)
                    
                    tab[which(tab$time_hours==time_point),"proportion_time_active"] <- bout_duration/(bout_duration+interbout_duration)
                  }
                }
              }
            }
          }
          
          TAB_ALL <- rbind(TAB_ALL,tab)
          ###Finally, write results
          ###individual behaviour
          behav <- read.table(paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",sep=""),header=T,stringsAsFactors = F)
          if (!"proportion_time_active"%in%names(behav)){
            behav[c("proportion_time_active","average_bout_speed_pixpersec","total_distance_travelled_pix")] <- NA
          }
          behav[match(as.character(interaction(tab$colony,tab$tag,tab$time_hours)),as.character(interaction(behav$colony,behav$tag,behav$time_hours))),c("proportion_time_active","average_bout_speed_pixpersec","total_distance_travelled_pix")]  <- tab[c("proportion_time_active","average_bout_speed_pixpersec","total_distance_travelled_pix")]
          options(digits=3)
          write.table(behav,file=paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",sep=""), row.names=F, col.names=T,append=F,quote=F)
          options(digits=16)
          write.table(traj,file=get(paste("traj_file",suffix,sep="")), row.names=F, col.names=T,append=F,quote=F)
        }
      }else{#if (!grepl("bout_id",output))
        print("Ant already processed.")
      }
    }#else (from if (!traj_file %in% alive))
    clean()
  }#for (traj_file in traj_files)
}# for (traj_folder in traj_list)
to_keep <- to_keep_ori