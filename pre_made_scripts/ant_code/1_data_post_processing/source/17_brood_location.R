####17_brood_location.R#####

####Calculates the distance between the center of gravity of the brood pile and the nest entrance based on the brood coordinates

###Created by Nathalie Stroeymeyt

###################################
### get input folder list
input_brood <- paste(data_path,"/original_data/brood_coordinates",sep="")
setwd(input_brood)
input_list  <- paste(input_brood,dir(),sep="/")

###define outputfolder
outputfolder <- paste(data_path,"/processed_data/collective_behaviour/pre_vs_post_treatment",sep="")
if (!file.exists(outputfolder)){dir.create(outputfolder,recursive = T)}

brood_COG <- NULL
for (folder in input_list){
  root_name     <- unlist(strsplit(folder,split="/"))[grepl("colony",unlist(strsplit(folder,split="/")))]
  colony        <- unlist(strsplit(root_name,split="_"))[grepl("colony",unlist(strsplit(root_name,split="_")))]
  colony_number <- as.numeric(gsub("colony","",colony))
  treatment     <- info[which(info$colony==colony_number),"treatment"]
  colony_size   <- info[which(info$colony==colony_number),"colony_size"]
  source(paste(code_path,"/heatmap_to_homerange_parameters.R",sep=""))
  
  if (grepl("PreTreatment",folder)){
    period <- "before"
  }else{
    period <- "after"
  }
  
  setwd(folder)
  brood_files <- list.files(pattern="brood")
  for (brood_file in brood_files){
    metadata     <- unlist(strsplit(gsub("\\.txt","",brood_file),split="_"))
    time_of_day  <- as.numeric(gsub("TD","",metadata[which(grepl("TD",metadata))]))
    time_hours   <- as.numeric(gsub("TH","",metadata[which(grepl("TH",metadata))]))
    
    ###read coordinate file
    t <- read.table(brood_file,header=T,stringsAsFactors = F)
    t <- data.frame(X=t$X,Y=t$Y,obs=1)
    
    #######modify t; indeed heatmaps and nest entrance position were calculated using cells of 5*5 pixels whereas t currently holds the real coordinates in pixels
    t$X <- floor(t$X/5)
    t$Y <- floor(t$Y/5)
    observations <- aggregate(obs~X+Y,FUN=sum,data=t)
    
    #######calculate centroid location
    centroid <- data.frame(x=mean(observations$X),y=mean(observations$Y))
    if ((centroid$y>=yforagingmin)&(centroid$y<=yforagingmax)){####if CoG is in foraging arena, then distance to nest entrance will be positive, whereas if it is inside the nest, it will be negative
      fac <- 1
    }else {
      fac <- -1
    }
    distance_to_entrance <- fac*euc.dist(centroid,nest_entrance)   #####this corresponds to a distance in 5*pixels
    distance_to_entrance <- distance_to_entrance * 5               ##### this corresponds to a distance in real pixels
    
    ####record data
    brood_COG <- rbind(brood_COG,data.frame(colony=colony,colony_size=colony_size,treatment=treatment,period=period,time_hours=time_hours,time_of_day=time_of_day,distance_from_nest_entrance_pix=distance_to_entrance,stringsAsFactors = F))
  }
}
brood_COG <- brood_COG[order(brood_COG$colony,brood_COG$time_hours),]
write.table(brood_COG,file=paste(outputfolder,"/brood_location.txt",sep=""),col.names=T,row.names=F,append=F,quote=F)