####5_zoneconverter_nest.R#####

### Calls C++ functions zoneconverter (https://github.com/laurentkeller/anttrackingUNIL)

### Takes as input a datfile and a plume file created using Plume (https://github.com/laurentkeller/anttrackingUNIL)
### Creates a dat files containing spatial information (nest zone vs. all other zones)

###Created by Nathalie Stroeymeyt
input_plume      <- paste(data_path,"/original_data/plume_files",sep="")
input_dat      <- paste(data_path,"/intermediary_analysis_steps/dat_files",sep="")
output_path     <- paste(data_path,"/intermediary_analysis_steps/dat_files_nest_zone",sep="")

if (!file.exists(output_path)){
  dir.create(output_path,recursive=T)
}

for (dat_file in unique(info_plume$dat_file)){
  subset <- info_plume[info_plume$dat_file==dat_file,]
  for (line in 1:nrow(subset)){
    box <- subset[line,"box"]
    plume_file <- subset[line,"plume_file"]
    
    if (line==1){
      command <- paste(executables_path,"/zone_converter -d ",input_dat,"/",dat_file," -p ",input_plume,"/", plume_file," -z nest -i 111 -o ",output_path,"/",dat_file," -b ",box, ";",sep="")
      system(command)
    }else{
      command <- paste("mv ",output_path,"/",dat_file," ",output_path,"/",dat_file,"_temp ;",sep="") 
      system(command)
      command <- paste(executables_path,"/zone_converter -d ",output_path,"/",dat_file,"_temp"," -p ",input_plume,"/", plume_file," -z nest -i 111 -o ",output_path,"/",dat_file," -b ",box, ";",sep="")
      system(command)
      command <- paste("rm ",output_path,"/",dat_file,"_temp ;",sep="")
      system(command)
    }
  }
}#dat_file
