########Other non-journal specific plotting parameters ##########
task_group_file <- "task_groups.txt"
refine          <- task_group_file
plot_type <- "bars"  ####bars or means or bars_points or boxplot or violinplot
if (plot_type=="boxplot"){
  relative_function <- "median"
}else{
  relative_function <- "mean"
}
par(mgp=c(0.8,0.1,0),mar=c(2,2,0.85,0),tcl=-0.2,lend=2,xpd=T,lwd=line_max)
pars <- par(mgp=c(0.8,0.1,0),mar=c(2,2,0.85,0),tcl=-0.2,lend=2,xpd=T,lwd=line_max)

#######Define general parameters that will be use throughout the analysis
queenid <- 665
night_start <- 18; light_start <- 6
stat_line <- -0.2
pix_to_mm_ratio <- max(c(0.0225877193,0.0229658793))


colour_palette <- c(viridis(10, alpha = 1, begin = 0, end = 1, direction = -1, option = "D")[c(1,3,5,7,9)],"#9ad0f3", "#0072B2")
colour_palette_workers <- viridis(10, alpha = 1, begin = 0, end = 1, direction = -1, option = "D")[1:5]
colour_palette_age <- rep(GetColorHex("lightskyblue"),2)

forager_colour <- colour_palette[1]
occasional_forager_colour <-  colour_palette[2]
nurse_colour <- colour_palette[3]
untreated_colour <- GetColorHex("grey60")
treated_colour <- GetColorHex("grey20")
queen_colour <- colour_palette[5]#queen_colour <- "mediumorchid4"
worker_colour <- colour_palette[2]
control_colour <- GetColorHex("skyblue1")
pathogen_colour <- GetColorHex("royalblue2")
random_colour <- GetColorHex("rosybrown1")
observed_colour <- GetColorHex("red4")
high_load_colour <- GetColorHex("springgreen4")
low_load_colour  <- GetColorHex("springgreen2")
#####define treatment and status names and labels
treatments <- c("control","pathogen","random","observed")
treatment_colours <- 
  c(control_colour,pathogen_colour,random_colour,observed_colour)
names(treatment_colours) <- treatments

statuses <- c("treated","untreated","queen"
              ,
              "forager","occasional_forager","nurse","queen"
              ,
              "worker","queen"
              ,
              "control","pathogen"
              ,
              "random","observed"
              ,
              "with_queen","not_with_queen"
              ,
              "high_predicted_load","low_predicted_load"
)
statuses_colours <- 
  c(treated_colour,untreated_colour,queen_colour
    ,
    forager_colour,occasional_forager_colour,nurse_colour,queen_colour
    ,
    worker_colour,queen_colour
    ,
    control_colour,pathogen_colour
    ,
    random_colour,observed_colour
    ,
    queen_colour,worker_colour
    ,
    high_load_colour,low_load_colour
  )
names(statuses_colours) <- statuses
statuses_colours <- statuses_colours[!duplicated(names(statuses_colours))]
full_statuses_names <- c("Treated\nforagers","Untreated\nworkers","Queen\n"
                         ,
                         "Untreated\nforagers","Occasional\nforagers","Nurses\n","Queen\n"
                         ,
                         "Workers","Queen\n"
                         ,
                         "Sham","Path."
                         ,
                         "Rand.","Obs."
                         ,
                         "Q\ncomm.","Other\ncomm."
                         ,
                         "High predicted load","Low predicted load"
);names(full_statuses_names) <- statuses
full_statuses_names <- gsub("Path.","",full_statuses_names)
status_order <- c(
  "queen"
  ,
  "nurse","occasional_forager","forager"
  ,
  "worker"
  ,
  "control","pathogen"
  ,
  "treated","untreated"
  ,
  "random","observed"
  ,
  "with_queen","not_with_queen"
  ,
  "high_predicted_load","low_predicted_load"
)

alphas <- c(1,1,0.5,1)
names(alphas) <- c("control","pathogen","random","observed")

high_threshold <- 0.0411
sq_mean_to_spore_nb_ratio <- 804427.048


light_start <- 6
dark_start <- 18

time_ori <- 12;unit <-3

Extended <- F

####keep all those until the end
to_keep <- c(ls(),"to_keep")
clean()