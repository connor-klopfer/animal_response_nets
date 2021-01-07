####11_randomise_interactions.R#####

### Sources C++ function randomised_edges 

#### Takes an interaction list as an input and returns a the same interaction interaction list in which the interaction partners have been randomised
#### Follows the 'Randomized edges (RE)' algorithm laid out by Holme and Saramäki 2012 (Physics Reports)

###Created by Nathalie Stroeymeyt

####################################
to_keep_ori <- to_keep

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp(paste(code_path,"/randomise_edges.cpp",sep=""))

###get input list
interac_list <- NULL
input_path1  <- paste(data_path,"/intermediary_analysis_steps/binned_interaction_lists/PreTreatment/observed",sep="")
setwd(input_path1)
interac_list <- c(interac_list,paste(input_path1,list.files(),sep="/"))
input_path2  <- paste(data_path,"/intermediary_analysis_steps/full_interaction_lists/PreTreatment/observed",sep="")
setwd(input_path2)
interac_list <- c(interac_list,paste(input_path2,list.files(),sep="/"))

to_keep <- c(ls(),"to_keep","i","interac")
for (i in 1:100){#####perform 100 randomisations
  print(paste("Performing randomisations",i,"out of 100..."))
  for (interac in interac_list){
    print(paste("Interaction file",which(interac==interac_list),"out of",length(interac_list)))
    root_name          <- unlist(strsplit(interac,split="/"))[grepl("colony",unlist(strsplit(interac,split="/")))]
    folder_name        <- gsub(paste("/",root_name,sep=""),"",interac)
    outputfolder       <- gsub("observed",paste("random_",paste(rep(0,3-nchar(i)),collapse=""),i,sep=""),folder_name)
    if (!file.exists(outputfolder)){dir.create(outputfolder)}
    outfile <- paste(outputfolder,root_name,sep="/")
    if (!file.exists(outfile)){
      #####read-in file
      interactions                 <- read.table(interac,header=T,stringsAsFactors=F)
      randomised_partners          <- randomise_edges(interactions[c("Tag1","Tag2","Startframe","Stopframe")])
      randomised_interactions      <- interactions
      randomised_interactions$Tag1 <- randomised_partners$Tag1;randomised_interactions$Tag2 <- randomised_partners$Tag2;
      write.table(randomised_interactions,file=outfile,col.names=T,row.names=F,quote=F,append=F)
    }
    clean();
  }
}
to_keep <- to_keep_ori