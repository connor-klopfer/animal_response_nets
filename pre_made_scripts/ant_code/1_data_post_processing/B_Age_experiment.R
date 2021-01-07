rm(list=ls())

########################################################################
###  INPUT TO DEFINE BY USER############################################
####Please fill in the path to the data folder, the code folder, and the c++ executables folder, as in the exemple below
data_path  <- "~/Repositories/Data_Repository/age_experiment"
code_path  <- "~/Repositories/Code_Repository/1_data_post_processing/source"
executables_path <- "~/executables"
########################################################################
###   END INPUT ########################################################
########################################################################

source(paste(code_path,"/libraries.R",sep=""))
source(paste(code_path,"/functions_and_parameters.R",sep=""))

to_keep <- c(ls(),"to_keep")

# #####Run analysis programs ####
source(paste(code_path,"/1_trackconverter.R",sep=""))
clean()
source(paste(code_path,"/2_define_deaths.R",sep=""))
clean()
source(paste(code_path,"/3_apply_rotation_to_datfiles.R",sep=""))
clean()
source(paste(code_path,"/5_zoneconverter_nest.R",sep=""))
clean()
source(paste(code_path,"/6_time_investment.R",sep=""))
clean()
source(paste(code_path,"/9_interaction_detection.R",sep=""))
clean()
source(paste(code_path,"/10_process_interaction_files.R",sep=""))
clean()
source(paste(code_path,"/11_randomise_interactions.R",sep=""))
clean()
source(paste(code_path,"/13_network_analysis.R",sep=""))
clean()
source(paste(code_path,"/14_summarise_interactions.R",sep=""))
clean()