####18_process_heatmaps.R#####

#### Takes heatmaps as an input, and calculates the part of the ant's home range that is included within in nest, 
####                                          the distance between the ant's center of gravity and the center of gravity of the untreated ants, 
####                                          and the degree of spatial overlap between the ant's space use and the brood pile

###Created by Nathalie Stroeymeyt###

#################################
to_keep_ori <- to_keep

options(digits=16) ; options(digits.secs=6) ; options("scipen" = 10)

#######get heatmap list
input_heatmaps <- paste(data_path,"/intermediary_analysis_steps/heatmaps/individual",sep="")
setwd(input_heatmaps)
folder_list  <- dir()
full_heatmap_list <- paste(input_heatmaps,list.files(recursive=T,pattern="ant"),sep="/")

#######get group heatmap list
input_group_heatmaps <- paste(data_path,"/intermediary_analysis_steps/heatmaps/group",sep="")
setwd(input_group_heatmaps)
full_group_heatmap_list <- paste(input_group_heatmaps,list.files(recursive=T,pattern="untreated"),sep="/")

#######get brood location list
input_brood <- paste(data_path,"/original_data/brood_coordinates",sep="")
setwd(input_brood)
full_brood_list <- paste(input_brood,list.files(recursive=T,pattern="brood"),sep="/")

#####Read behaviour file 
behav <- read.table(paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",sep=""),header=T,stringsAsFactors = F)
#####Create new columns
if (!"within_nest_home_range"%in%names(behav)){
  behav[c("within_nest_home_range","distance_antCoG_to_colonyCoG","BA_between_ant_and_brood")] <- NA
}

for (input_folder in folder_list){
  print(input_folder)
  heatmap_list <- data.frame(file=full_heatmap_list[grepl(input_folder,full_heatmap_list)],stringsAsFactors = F)
  for (i in 1:nrow(heatmap_list)){
    time_info <- paste(unlist(strsplit(heatmap_list[i,"file"],split="_"))[grepl("TD",unlist(strsplit(heatmap_list[i,"file"],split="_")))|grepl("TH",unlist(strsplit(heatmap_list[i,"file"],split="_")))],collapse="_")
    heatmap_list[i,"time_info"] <- time_info
  }
  whens <- unique(heatmap_list$time_info)###get the list of time bins
  for (when in whens){
    print(when)
    ####CALCULATE UNTREATED COG POSITION for that time
    reference_untreated <- full_group_heatmap_list[which((grepl(input_folder,full_group_heatmap_list))&(grepl(when,full_group_heatmap_list)))]
    t_untreated         <- read.table(reference_untreated, sep=",")
    centroid_untreated  <- data.frame(x=sum(t_untreated[,1]*t_untreated[,3])/sum(t_untreated[,3]),y=sum(t_untreated[,2]*t_untreated[,3])/sum(t_untreated[,3]))
    
    ####READ BROOD HEATMAP for that time
    reference_brood     <- full_brood_list[which((grepl(input_folder,full_brood_list))&(grepl(when,full_brood_list)))]
    t_brood <- read.table(reference_brood,header=T,stringsAsFactors = F)
    t_brood <- data.frame(X=t_brood$X,Y=t_brood$Y,obs=1)
    #######modify t_brood; indeed heatmaps and nest entrance position were calculated using cells of 5*5 pixels whereas t_brood currently holds the real coordinates in pixels
    t_brood$X <- floor(t_brood$X/5)
    t_brood$Y <- floor(t_brood$Y/5)
    observations_brood <- aggregate(obs~X+Y,FUN=sum,data=t_brood)
    locations_brood <- NULL
    if (nrow(observations_brood)>0){
      for (value in unique(observations_brood[,3])){
        subset <- observations_brood[observations_brood[,3]==value,1:2]
        locations_brood <- rbind( locations_brood,data.frame(X=rep(subset[,1],value),Y=rep(subset[,2],value)))
      }
    } 
    heatmap_files <- heatmap_list[which(heatmap_list$time_info==when),"file"]
    
    to_keep <- c(ls(),"to_keep","heatmap","customised_grid_all")
    for (heatmap in heatmap_files){
      ####get metadata
      colony_metadata <- unlist(strsplit(heatmap,split="/"))[grepl("colony", unlist(strsplit(heatmap,split="/")))]
      colony      <- unlist(strsplit(colony_metadata,split="_"))[grepl("colony",unlist(strsplit(colony_metadata,split="_")))]
      colony_number <- as.numeric(gsub("colony","",colony))
      treatment     <- info[which(info$colony==colony_number),"treatment"]
      colony_size   <- info[which(info$colony==colony_number),"colony_size"]
      source(paste(code_path,"/heatmap_to_homerange_parameters.R",sep=""))
      time_point  <- unlist(strsplit(colony_metadata,split="_"))[grepl("Treatment",unlist(strsplit(colony_metadata,split="_")))]
      if (time_point=="PreTreatment"){period <- "before"}else{period <- "after"}
      
      ant_metadata    <- unlist(strsplit(heatmap,split="/"))[grepl("txt", unlist(strsplit(heatmap,split="/")))]
      ant_metadata    <- unlist(strsplit(ant_metadata,split="_"))
      tag             <- as.numeric(gsub("ant","",ant_metadata[grepl("ant",ant_metadata)&!grepl("txt",ant_metadata)]))
      time_of_day     <- as.numeric(gsub("TD","",ant_metadata[grepl("TD",ant_metadata)]))
      time_hours      <- as.numeric(gsub("TH","",ant_metadata[grepl("TH",ant_metadata)]))
      
      #####Check if whether ant has already been processed, and do the rest only if it is the case
      if (is.na(behav[which(behav$colony==colony&behav$tag==tag&behav$time_hours==time_hours),"distance_antCoG_to_colonyCoG"])){
        print(paste("Analysing ",heatmap))
        ####read heatmap
        t   <- read.table(heatmap, sep=",")####read heatmap output file
        
        ###define grid (to do only the first time)
        if (!exists("customised_grid_all")){
          for_grid <- t[,1:2]###list of all x-y coordinates of the plot (each existing cell (even unvisited, empty cells) defined on 1 line by its x,y coordinates)
          names(for_grid) <- c("x","y")
          customised_grid_all <- (unique(round(for_grid/grid_parameter)))*grid_parameter
          coordinates(customised_grid_all)=c("x","y")
          gridded(customised_grid_all) <- TRUE
        }
        
        ###use to t to make a location matrix
        observations <- t[t[,3]!=0,]
        locations <- NULL
        if (nrow(observations)>0){
          count <- 0
          for (value in unique(observations[,3])){
            count <- count + 1
            subset <- observations[observations[,3]==value,1:2]
            locations  <- rbind(locations,data.frame(X=rep(subset[,1],value),Y=rep(subset[,2],value)))
          }
        } 
        
        
        ####STEP 1 - calculate part of home range located within the nest
        if (nrow(locations)>5){
          locations_sp <- SpatialPoints(locations)
          bandwidth <- kernelUD(locations_sp, h="href",grid=customised_grid_all,extent=0)@h$h###get bandwidth (without bounder - because there is a bug when using bounder and href)
          if (exists("kud")){rm(list=c("kud"))}
          try(kud <- kernelUD(locations_sp, h=bandwidth,boundary=bounder_all,grid=customised_grid_all,extent=0),silent=T)###calculate utilization distribution with border definition
          if (exists("nest_home_range")){rm(list=c("nest_home_range"))}
          if (exists("kud")){
            if(exists("vud")){rm(list=c("vud"))}
            try(vud <- getvolumeUD(kud,standardize=T),silent=T)
            if(exists("vud")){
              hr <- as.data.frame(vud)[,1]
              hrtemp <- data.frame(as.numeric(hr <= 95))
              coordinates(hrtemp) <- coordinates(vud)
              hrnest <- hrtemp[( hrtemp@coords[,2]>=ynestmin)&( hrtemp@coords[,2]<=ynestmax),]
              nest_home_range <- grid_parameter*grid_parameter*sum(hrnest@data)
              rm(list=c("hrnest"))
            }
          }
          if (exists("nest_home_range")){nest_home_range <-       nest_home_range *5*5}else{nest_home_range <- NA}    #####this corresponds to an are in (pixels)^2
        }else{
          nest_home_range <- 0
        }
        
        ####STEP 2 - calculate CoG location and overlap with brood
        ###############Centroid location
        centroid <- data.frame(x=mean(locations$X),y=mean(locations$Y))
        distance_to_colony <- euc.dist(centroid,centroid_untreated)   #####this corresponds to a distance in 5*pixels
        distance_to_colony <-  distance_to_colony * 5              ##### this corresponds to a distance in real pixels
        ###############Overlap with brood
        ####make a multi-animal spatial points
        dataset <- rbind(data.frame(Name="focal_ant",locations,stringsAsFactors = F),data.frame(Name="brood",locations_brood,stringsAsFactors = F))
        dataset$Name <- factor(dataset$Name)
        SPDF    <- SpatialPointsDataFrame(rbind(locations,locations_brood), dataset)
        ####calculate the utilization distributions with border definition
        if(exists("uds")){rm(list=("uds"))}
        try(uds     <- kernelUD(SPDF[,1], h=fixed_bandwidth, grid=customised_grid_all, boundary=bounder_all, extent=0),silent=T)
        if (exists("uds")){
          ####calculate overlap using Bhattacharyya's Affinity method
          nrows <- uds[[1]]@grid@cells.dim[2]
          ncols <- uds[[1]]@grid@cells.dim[1]
          densit_m_list <- list()
          for (r in 1:length(uds)){
            ####extract data
            densit <- data.frame(X=uds[[r]]@coords[,1],Y=uds[[r]]@coords[,2],ud=uds[[r]]@data)						  ## pixels denote density
            ####reorder
            densit <- densit[order(densit$X,densit$Y),]
            densit_m           <-t( matrix( (densit$ud/sum(densit$ud)) ,  ncol=ncols , nrow=nrows, byrow = F) )
            densit_m_list[[r]] <- densit_m   		  		  ## UPSIDE DOWN IMAGE	
          }
          overlap_map <- matrix (ncol=nrows,nrow=ncols) 
          for (e in 1:ncols){	   					   ## for each row of the distance map image matrix
            for (f in 1:nrows){					   ## for each column of the distance map image matrix
              ## for each pixel, find the min 'height' for the pair of ants.
              overlap_map[e,f]  <- sqrt(densit_m_list[[1]][e,f]) * sqrt(densit_m_list[[2]][e,f])
            }#f
          }#e
          overlap <- sum(overlap_map)
        }else{
          overlap <- NA
        }
        ####Now add data to behaviour file
        behav[which(behav$colony==colony&behav$tag==tag&behav$time_hours==time_hours),c("within_nest_home_range","distance_antCoG_to_colonyCoG","BA_between_ant_and_brood")]  <- c(nest_home_range,distance_to_colony,overlap)
      }
      clean()
    }
    options(digits=3) 
    write.table(behav,file=paste(data_path,"/processed_data/individual_behaviour/pre_vs_post_treatment/individual_behavioural_data.txt",sep=""), row.names=F, col.names=T,append=F,quote=F)
    rm(list=c("reference_untreated","t_untreated","centroid_untreated","reference_brood","t_brood","observations_brood","locations_brood"))
  }
  rm(list=c("heatmap_list","when"))
}
to_keep <- to_keep_ori