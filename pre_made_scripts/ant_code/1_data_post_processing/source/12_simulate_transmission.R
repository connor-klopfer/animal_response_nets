####12_simulate_transmission.R#####

#### Sources C++ function simulate_transmission from source folder

#### Takes an interaction list and simulates the transmission of that agent from a list of originally contaminated workers to the rest of the colony

###Created by Nathalie Stroeymeyt

####################################
to_keep_ori <- to_keep

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp(paste(code_path,"/simulate_transmission.cpp",sep=""))

#################################
options(digits=16) ; options(digits.secs=6) ; options("scipen" = 10)

#### get input file list
input_path         <- paste(data_path,"/intermediary_analysis_steps/full_interaction_lists",sep="")
setwd(input_path)  
input_files        <- paste(input_path,"/",list.files(recursive=T,pattern="colony"),sep="")

#### arguments 
N_SIM  <- 500
if (!grepl("survival",data_path)){
  seed_files <- c("treated_workers.txt","random_workers.txt","frequent_foragers.txt","occasional_foragers.txt","nurses.txt","low_degree.txt","high_degree.txt")
}else{
  seed_files <- c("treated_workers.txt")
}
to_keep <- c(ls(),"to_keep","seed_file","outputfolder","interac_folders","reorder","interac_folder","interac_list","summary_collective","summary_individual","interac")
for (seed_file in seed_files ){
  print(paste("Seeds =",gsub("\\.txt","",seed_file)))
  if (seed_file=="treated_workers.txt"){
    if (!grepl("survival",data_path)){
      outputfolder    <- paste(data_path,"/transmission_simulations/pre_vs_post_treatment/experimentally_exposed_seeds",sep="")
      interac_folders <- "observed"
      reorder         <- T
    }else{
      outputfolder    <- paste(data_path,"/transmission_simulations/post_treatment/experimentally_exposed_seeds",sep="")
      interac_folders <- "observed"
      reorder         <- F
    }
  }else{
    outputfolder    <- paste(data_path,"/transmission_simulations/random_vs_observed/",gsub("\\.txt","",seed_file),sep="")
    interac_folders <- c("observed",1:100)
    reorder         <- F
  }
  if (!file.exists(outputfolder)){dir.create(outputfolder,recursive = T)}  
  
  for (interac_folder in interac_folders){
    if (interac_folder!="observed"){interac_folder <- paste("random_",paste(rep(0,3-nchar(interac_folder)),collapse=""),interac_folder,sep="")}
    if (!file.exists(paste(outputfolder,"/individual_simulation_results_",interac_folder,".txt",sep=""))){
      summary_collective <- NULL
      summary_individual <- NULL

      interac_list <- input_files[grepl(interac_folder,input_files)]
      if (seed_file!="treated_workers.txt"){interac_list <- interac_list[grepl("PreTreatment",interac_list)]}
      for (interac in interac_list){
        print(interac)
        ####get colony info
        root_name          <- gsub("_interactions.txt","",unlist(strsplit(interac,split="/"))[grepl("colony",unlist(strsplit(interac,split="/")))])
        colony             <- unlist(strsplit(root_name,split="_"))[grepl("colony",unlist(strsplit(root_name,split="_")))]
        colony_number      <- as.numeric(gsub("colony","",colony))

        ####get period info
        if (grepl("PreTreatment",root_name)){period="before"}else{period="after"}

        ####read interactions
        interaction_table <- read.table(interac,header=T,stringsAsFactors=F)
        interaction_table <- interaction_table[order(interaction_table$Starttime),]

        ####read tag list to define list of live ants
        tagfile             <- tag_list[grepl(colony,tag_list)]
        tag        <- read.tag(tagfile)$tag; names(tag)[names(tag)=="#tag"] <- "tag"; tag <- tag[which(tag$tag!="#tag"),]
        if (!grepl("survival",data_path)){
          alive      <- tag[which(tag$final_status=="alive"),"tag"]
        }else{
          alive      <- tag[which(as.numeric(tag$death)==0|as.numeric(tag$death)>=max(interaction_table$Stopframe,na.rm=T)),"tag"]
        }

        ###read seeds
        seeds              <- read.table(paste(data_path,"/original_data/seeds/",seed_file,sep=""),header=T,stringsAsFactors = F)
        seeds              <- seeds[which(seeds$colony==colony_number),"tag"]
        seeds              <- seeds[which(seeds%in%alive)]


        ####read time aggregation info for both current file and PostTreatment file if relevant (i.e. if reorder=T), and get start time from that
        if (!grepl("survival",data_path)){
          splitfile          <- split_list[grepl(root_name,split_list)]
          splitinfo          <- read.table(splitfile,header=T,stringsAsFactors = F)
          if (reorder){
            splitfile_Post   <- split_list[grepl(gsub("Pre","Post",root_name),split_list)]
            splitinfo_Post   <- read.table(splitfile_Post,header=T,stringsAsFactors = F)
          }
          time_start         <- min(splitinfo$time,na.rm=T)
        }else{
          time_start         <- min(interaction_table$Starttime,na.rm=T)
        }

        ####for simulations to be comparable between the Pre-Treatment and Post-treatment periods, they should start at the same time of day
        ####so in that case, Pre-Treatment interactions should be reordered
        if (reorder & period=="before") {
          ####start time must be the same time of day as start time of after.
          time_start <- min(splitinfo_Post$time,na.rm=T) - (24*3600)
          ####get the corresonding row number in the interaction table
          start_index <- min(which(interaction_table$Starttime>=time_start))
          ####split the interaction table into two halves
          interactions_unchanged <- interaction_table[start_index:nrow(interaction_table),]
          interactions_to_change <- interaction_table[1:(start_index-1),]
          ####modify frame number and time of interactions_to_change
          interactions_to_change[c("Starttime","Stoptime")] <- interactions_to_change[c("Starttime","Stoptime")] + 24 * 3600
          frame_offset <- (interactions_unchanged[nrow(interactions_unchanged),"Startframe"]+floor((interactions_to_change[1,"Starttime"]-interactions_unchanged[nrow(interactions_unchanged),"Starttime"])*2))-interactions_to_change[1,"Startframe"]
          interactions_to_change[c("Startframe","Stopframe")] <- interactions_to_change[c("Startframe","Stopframe")] + frame_offset
          interaction_table <- rbind(interactions_unchanged,interactions_to_change)
        }
        ####Create antlist object
        antlist                                                 <- data.frame(tag=alive,status=0,load=0,stringsAsFactors = F)
        antlist[which(antlist$tag%in%seeds),c("status","load")] <- 1

        ####Change formats to match formats expected by the C++ program
        antlist$tag                                             <- as.integer(antlist$tag )
        antlist$status                                          <- as.integer(antlist$status )
        antlist$load                                            <- as.numeric(antlist$load )

        interaction_table$Tag1                                  <- as.integer(interaction_table$Tag1 )
        interaction_table$Tag2                                  <- as.integer(interaction_table$Tag2 )
        interaction_table$Startframe                            <- as.integer(interaction_table$Startframe )
        interaction_table$Stopframe                             <- as.integer(interaction_table$Stopframe )
        interaction_table$Starttime                             <- as.numeric(interaction_table$Starttime )
        interaction_table$Stoptime                              <- as.numeric(interaction_table$Stoptime )

        ####Perform simulations
        print("Performing simulations...")
        simulations <- NULL
        for (i in 1:N_SIM){
          simulations <- rbind(simulations,data.frame(sim_number=i,simulate_transmission(interaction_table,antlist,time_start)))
        }

        ####Summarise simulations
        ##Step 1: modify columns
        simulations$relative_contamination_time                      <- simulations$relative_contamination_time/3600 ###express contamination time in hours instead of seconds
        simulations                                                  <- simulations[,!names(simulations)%in%c("absolute_contamination_time","infectious")] ##remove unncessary columns
        simulations["contaminated"]                                  <- 1          ###add a column containing information on whether the ant was contaminated or not during the simulation
        simulations["status"]                                        <- "untreated"
        simulations[which(simulations$tag%in%seeds),"status"]        <- "treated"
        ##Step 2: add individuals that did not get contaminated during the simulations
        all_individuals <- expand.grid(sim_number=c(1:N_SIM),tag=alive)
        simulations     <- merge(all_individuals,simulations,all.x=T,all.y=T)
        if (nrow( simulations[which(is.na(simulations$relative_contamination_time)),])>0){
          simulations[which(is.na(simulations$relative_contamination_time)),"status"]     <- "untreated"
          simulations[which(is.na(simulations$relative_contamination_time)),c("relative_contamination_time","contaminated_by","final_load","contaminated")]     <- rep(c(24,-1,0,0) ,each=nrow(simulations[which(is.na(simulations$relative_contamination_time)),]))
        }
        simulations     <- simulations[order(simulations$sim_number,simulations$relative_contamination_time),]

        ##Step3: add further individual-level information
        ############# Add infection rank
        ranks <- aggregate(relative_contamination_time~sim_number,function(x)rank(x,ties.method="min"),data=simulations[which(simulations$status!="treated"),])
        for (sim_number in unique(ranks$sim_number)){
          simulations[which(simulations$sim_number==sim_number&simulations$status=="treated"),"rank"] <- 0
          simulations[which(simulations$sim_number==sim_number&simulations$status!="treated"),"rank"] <- as.numeric(ranks[which(ranks$sim_number==sim_number),c(2:ncol(ranks))])
        }
        ############# Add high load/low load
        simulations["high_level"] <- as.numeric(simulations$final_load>high_threshold)
        simulations["low_level"] <- as.numeric(simulations$final_load>0&simulations$final_load<=high_threshold)

        ##Step4: summarise simulations - colony-level data
        #########Step 4.1: Prevalence, Mean load, Prop. high level, Prop. low level, Mean load
        colony_level         <- aggregate(na.rm=T,na.action="na.pass",cbind(contaminated,final_load,high_level,low_level)~sim_number,FUN=mean,data=simulations[which(simulations$status!="treated"),])
        names(colony_level)  <- c("sim_number","Prevalence","Mean_load","Prop_high_level","Prop_low_level")
        #########Step 4.2: Load skewness
        skewness             <- aggregate(na.rm=T,na.action="na.pass",final_load~sim_number,FUN=skewness,data=simulations[which(simulations$status!="treated"),])
        names(skewness)      <- c("sim_number","Load_skewness")
        colony_level         <- merge(colony_level,skewness)
        #########Step 4.3: Queen load
        queen                <- simulations[which(simulations$tag==queenid),c("sim_number","final_load")]
        names(queen)         <- c("sim_number","Queen_load")
        colony_level         <- merge(colony_level,queen)
        #########Step 4.4: Transmission rate
        colony_level["logistic_r"] <- NA
        for (sim in unique(simulations$sim_number)){
          #######define table that will be used in the non linear fit
          subset <- simulations[which(simulations$sim_number==sim),]
          #######remove from subset the ants that were not infected during the simulation and were added later
          subset <- subset[which(subset$contaminated!=0),]
          #######sort by time
          subset <- subset[order(subset$relative_contamination_time,subset$status),]
          #######get population size
          pop_size <- sum(subset[1,c("nb_susceptible","nb_contaminated")])
          #######get relevant data
          spread                 <- subset[c("relative_contamination_time","rank")]
          spread["nb_seeds"]     <- nrow(subset[which(subset$status=="treated"),])
          spread["nb_contaminated"] <- spread$rank+spread$nb_seeds
          spread <- spread[which(!duplicated(spread[c("relative_contamination_time","nb_contaminated")])),]
          spread["proportion_contaminated"] <- spread$nb_contaminated/pop_size

          #######try logistic fit
          y <- spread$proportion_contaminated
          x <- spread$relative_contamination_time
          P_zero <- spread[which(spread$rank==0),"proportion_contaminated"]
          #K <- 1
          if (exists("fit_logistic")){rm(list=c("fit_logistic"))}
          try(fit_logistic <- nls(
            y ~ (K*P_zero*exp(r*x))
            /
              (K+(P_zero*(-1+exp(r*x))))
            ,
            start=list(r=1,K=1)
            #start=list(r=1)
            ,
            control=list(maxiter = 1000)
          )
          ,silent=T)
          if (exists("fit_logistic")){
            r <- summary(fit_logistic)$coefficients["r","Estimate"]
            K <- summary(fit_logistic)$coefficients["K","Estimate"]
            colony_level[which(colony_level$sim_number==sim),"logistic_r"] <- r
          }
          rm(list=c("fit_logistic","r","K"))
        }
        #########Step 4.5: Take average over all simulations
        colony_level         <- colMeans(colony_level[names(colony_level)!="sim_number"],na.rm=T)
        #########Step 4.6: Store
        summary_collective   <- rbind(summary_collective,cbind(data.frame(colony=colony,colony_size=info[which(info$colony==colony_number),"colony_size"],treatment=info[which(info$colony==colony_number),"treatment"],period=period,time_hours=interaction_table[1,"time_hours"],time_of_day=interaction_table[1,"time_of_day"],stringsAsFactors = F),t(colony_level)))

        ##Step5: summarise simulations - individual-level data
        #########Step 5.1: Final_load, Prob. contamination, Prob. high level, Prob. low level
        individual_level         <- aggregate(na.rm=T,na.action="na.pass", cbind(final_load,contaminated,high_level,low_level) ~ tag,FUN=mean,data=simulations)
        names(individual_level)  <- c("tag","simulated_load","probability_of_transmission","probability_high_level","probability_low_level")
        #########Step 5.2: Transmission latency and transmission rank
        individual_level[c("transmission_latency","transmission_rank")] <- NA
        for (ant in alive){
          if (all(simulations[which(simulations$tag==ant),"contaminated"]==1)){
            individual_level[which(individual_level$tag==ant),"transmission_latency"] <- mean(simulations[which(simulations$tag==ant),"relative_contamination_time"])
            individual_level[which(individual_level$tag==ant),"transmission_rank"]    <- mean(simulations[which(simulations$tag==ant),"rank"])

          }else{
            model                                                                  <- coxph(Surv(relative_contamination_time,contaminated)~1,data=simulations[which(simulations$tag==ant),])
            mean_data                                                              <- summary(survfit(model),rmean="common")$table
            individual_level[which(individual_level$tag==ant),"transmission_latency"] <-  mean_data["*rmean"]

            model                                                                  <- coxph(Surv(rank,contaminated)~1,data=simulations[which(simulations$tag==ant),])
            mean_data                                                              <- summary(survfit(model),rmean="common")$table
            individual_level[which(individual_level$tag==ant),"transmission_rank"] <-  mean_data["*rmean"]

          }
        }
        #########Step 5.3: Store
        individual_level$antid          <- as.character(interaction(colony,individual_level$tag))
        individual_level$status         <- "untreated";individual_level[which(individual_level$tag%in%seeds),"status"] <- "treated"
        individual_level["colony"]      <- colony
        individual_level["colony_size"] <- info[which(info$colony==colony_number),"colony_size"]
        individual_level["treatment"]   <- info[which(info$colony==colony_number),"treatment"]
        individual_level["period"]      <- period
        individual_level["time_hours"]  <- interaction_table[1,"time_hours"]
        individual_level["time_of_day"] <- interaction_table[1,"time_of_day"]

        individual_level       <- individual_level[c("colony","colony_size","treatment","tag","antid","status","period","time_hours","time_of_day","simulated_load","probability_of_transmission","probability_high_level","probability_low_level","transmission_latency","transmission_rank")]
        summary_individual     <- rbind(summary_individual,individual_level)
        clean()
      }
      summary_collective <- summary_collective[order(summary_collective$colony,summary_collective$time_hours),]
      summary_individual <- summary_individual[order(summary_individual$colony,summary_individual$tag,summary_individual$time_hours),]

      if (grepl("random_vs_observed",outputfolder)){
        summary_collective["randy"] <- interac_folder
        summary_individual["randy"] <- interac_folder
      }
      write.table(summary_collective,paste(outputfolder,"/collective_simulation_results_",interac_folder,".txt",sep=""),col.names=T,row.names=F,quote=F,append=F)
      write.table(summary_individual,paste(outputfolder,"/individual_simulation_results_",interac_folder,".txt",sep=""),col.names=T,row.names=F,quote=F,append=F)

    }
  }
  
  if (grepl("random_vs_observed",outputfolder )){
    setwd(outputfolder)
    file_list <- list.files(pattern="individual")
    all_indiv_results <- NULL
    for (file in file_list){
      temp <- read.table(file, header=T,stringsAsFactors = F)
      if (grepl("random",unique(temp$randy))){
        temp$randy <- "random"
      }
      all_indiv_results <- rbind(all_indiv_results,temp)
    }
    all_indiv_results <- aggregate(na.rm=T,na.action="na.pass",cbind(simulated_load,probability_of_transmission,probability_high_level,probability_low_level,transmission_latency,transmission_rank)~.,FUN=mean,data=all_indiv_results)
    all_indiv_results$treatment <- all_indiv_results$randy
    write.table(all_indiv_results,file=paste(outputfolder,"/summarised_individual_results.dat",sep=""),append=F,quote=F,row.names=F,col.names=T)
    
  }
  
  
  
}
to_keep <- to_keep_ori

