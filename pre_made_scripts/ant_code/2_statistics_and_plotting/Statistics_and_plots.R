rm(list=ls())

#####Overall parameters and functions ##########
####define folders
figurefolder <- "~/figures"
disk_path    <- "~/Repositories/Data_Repository"
source_path  <- "~/Repositories/Code_Repository/2_statistics_and_plotting/source"

####source programs
source(paste(source_path,"/libraries.R",sep=""))
source(paste(source_path,"/plotting_parameters.R",sep=""))
source(paste(source_path,"/functions.R",sep=""))
source(paste(source_path,"/analysis_parameters.R",sep=""))

############## Figure 1##############
#########Define parameters for Figure 1 
figure_height <- page_height/1.8
#####at the top of figure 1 we will have 2 images. I now give their height in inches 
image_height <- 1.75
##### parameters defining the height of each plotting row in figure 1
a <- (double_col/2)/figure_height
d <- image_height/figure_height
c <- 0.05
b <- (1 -a -c -d)
to_keep <- c(to_keep,"a","b","c","d")
######Open pdf Figure 1 #####
pdf(file=paste(figurefolder,"/Figure1.pdf",sep=""),family=text_font,font=text_font,bg="white",width=double_col,height=figure_height,pointsize=pointsize_more_than_2row2col)
##### Set-up layout and plot parameters #####
par(pars)
ncoli <- 26
layout(matrix(c(rep(8,ncoli),
                rep(1,ncoli/2),rep(2,ncoli/2),####networks at the top
                rep(8,ncoli),
                rep(7,ncoli/13),rep(3,3*ncoli/13),rep(4,3*ncoli/13), rep(6,3*ncoli/13),rep(5,3*ncoli/13)###second row, first half
                
), 4, ncoli, byrow = TRUE),heights=c(d,a,c,b))
####First plot networks #######
root_path <- paste(disk_path,"/main_experiment",sep="")######linux laptop
plot_network(case="topology_comparison",which_to_draw=c("PreTreatment_observed","PostTreatment_observed"))
### clean before next step
clean(); 
Sys.sleep(2)

####Second, plot changes in networks ########
root_path <- paste(disk_path,"/main_experiment/processed_data",sep="")
variable_list <- c("modularity","clustering","task_assortativity","efficiency")
names(variable_list) <- c("modularity","clustering","assortativity","efficiency")
transf_variable_list <- c("none","none","none","none")
analysis <- list(variable_list=variable_list,
                 transf_variable_list=transf_variable_list,
                 violin_plot_param = list(c(1.5,0,-0.02,0.25,0.01),c(1.5,-0.02,-0.02,0.25,0.01),c(1.5,-0.02,-0.02,0.25,0.01),c(1.5,-0.02,-0.02,0.25,0.001)))

plot_before_vs_after(root_path=root_path, data_path=paste(root_path,"/network_properties/pre_vs_post_treatment/all_workers",sep=""),analysis=analysis,status="collective",collective=T,pool_plot=T,pattern="colony_data.txt")
######## clean before next step###
clean(); 
Sys.sleep(2)

####Third, add letters  #########
par(xpd=NA)
##LETTERS
y_text <- grconvertY((1-1/80), from='ndc')
x_text1 <- grconvertX(0+1/80, from='ndc')
x_text2 <- grconvertX(0.6+1/80, from='ndc')
text(x_text1,y_text,labels=panel_casse("a"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
text(x_text2,y_text,labels=panel_casse("b"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
x_text1 <- grconvertX(0+1/80, from='ndc'); x_text2 <- grconvertX(1/2+1/80, from='ndc')
y_text <- grconvertY((1-d-1/80), from='ndc')
text(x_text1,y_text,labels=panel_casse("c"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
y_text1 <- grconvertY((1-d-a-c/2-1/80), from='ndc')
text(x_text1,y_text1,labels=panel_casse("d"),font=panel_font, cex=panel_cex,adj=c(0,0),xpd=NA)
par(xpd=F)
####Fourth, Close figure 1 #######
dev.off()

############## ##############
########Extended Data 15: pathogen-induced changes for untreated-only networks########################
pdf(file=paste(figurefolder,"/Extended_Data_15.pdf",sep=""),family=text_font,font=text_font,bg="white",width=three_col,height=page_height*0.3,pointsize=pointsize_more_than_2row2col)
par(pars)
ncoli <- 26
layout(matrix(c(rep(5,ncoli),
                rep(5,ncoli/13),rep(1,3*ncoli/13),rep(2,3*ncoli/13), rep(4,3*ncoli/13),rep(3,3*ncoli/13)###second row, first half
), 2, ncoli, byrow = TRUE)
,heights=c(0.05,0.95)
)
####Plot changes in networks  - untreated workers only
root_path <- paste(disk_path,"/main_experiment/processed_data",sep="")
variable_list <- c("modularity","clustering","task_assortativity","efficiency")
names(variable_list) <- c("modularity","clustering","assortativity","efficiency")
transf_variable_list <- c("none","none","none","none")
analysis <- list(variable_list=variable_list,
                 transf_variable_list=transf_variable_list,
                 violin_plot_param = list(c(1.5,0,-0.02,0.25,0.01),c(1.5,-0.02,-0.02,0.25,0.01),c(1.5,-0.02,-0.02,0.25,0.01),c(1.5,-0.02,-0.02,0.25,0.001)))

plot_before_vs_after(root_path=root_path, data_path=paste(root_path,"/network_properties/pre_vs_post_treatment/untreated_only",sep=""),analysis=analysis,status="collective",collective=T,pool_plot=T,pattern="colony_data.txt")
dev.off()
######## clean before next step###
Sys.sleep(2)
clean(); 
Sys.sleep(2)

######Open pdf Figure 2 #####
pdf(file=paste(figurefolder,"/Figure2.pdf",sep=""),family=text_font,font=text_font,bg="white",width=double_col,height=page_height*0.25,pointsize=pointsize_less_than_2row2col)
######Plot qPCR data #####
full_statuses_names_ori <- full_statuses_names
full_statuses_names[full_statuses_names%in%c("Foragers","Untreated\nforagers")] <- "Foragers\n"

statuses_colours_ori <- statuses_colours
statuses_colours[names(statuses_colours)%in%c("queen","forager","nurse")] <- "black"

par(pars)
par_mar_ori <- par()$mar
par(mar=par_mar_ori+c(1,0,1,1))
widz <- c(2,1.5)
layout(matrix(c(1,2),nrow=1),widths=widz)
translated_high_threshold <- plot_qpcr(experiments=c("age_experiment","main_experiment"))
to_keep <- c(to_keep,"translated_high_threshold")
####Add letters ####
par(xpd=NA)
x_text1 <- grconvertX(1/80, from='ndc');x_text2 <- grconvertX(widz[1]/sum(widz)+1/80, from='ndc')
y_text <- grconvertY(0.97, from='ndc')
text(x_text1,y_text,labels=panel_casse("a"),font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)
text(x_text2,y_text,labels=panel_casse("b"),font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)
par(xpd=F)
full_statuses_names <- full_statuses_names_ori
statuses_colours <- statuses_colours_ori
par(mar=par_mar_ori)
####Close figure 2#####
dev.off()
######## clean before next step###
clean();
Sys.sleep(2)

######Open pdf Figure 3 #####
pdf(file=paste(figurefolder,"/Figure3.pdf",sep=""),family=text_font,font=text_font,bg="white",width=double_col,height=page_height*0.45,pointsize=pointsize_more_than_2row2col)
##### Set-up layout and plot parameters #####
par(pars)
ncoli <- 4
heits <- c(5/9,0.05,4/9)
widz <- c(0.075,0.45,0,0.45)
layout(matrix(c(rep(1,ncoli/2),rep(4,ncoli/2),
                rep(5,ncoli),
                rep(5,ncoli/4),rep(2,ncoli/4),rep(5,ncoli/4),rep(3,ncoli/4)
), 3, ncoli, byrow = TRUE),heights=heits,widths = widz)
######First, plot distribution and threshold identification#######
plot_distribution(experiments="main_experiment",desired_treatments=c("pathogen"))
####Second, individual simulation results - comparison #######
root_path <- paste(disk_path,"/main_experiment",sep="")
queen <- T; treated <- F; nurses <- T; foragers <- T; 
unit_ori <- unit; unit <- 24
time_window <- 24

variable_list <- c("probability_high_level","probability_low_level")
names(variable_list) <- c("prob. receiving high load","prob. receiving low load")
transf_variable_list <- c("power3","sqrt")
predictor_list <- c("task_group","task_group")
analysis <- list(variable_list=variable_list,
                 transf_variable_list=transf_variable_list,
                 predictor_list=predictor_list,
                 violin_plot_param = list(c(1.5,0,-0.02,0.35,0.11),c(1.5,-0.02,-0.02,0.35,0.11)))
plot_before_vs_after(root_path=root_path, data_path=paste(root_path,"/transmission_simulations/pre_vs_post_treatment/experimentally_exposed_seeds",sep=""),analysis=analysis,status="all",collective=F,pool_plot=F,pattern="individual_simulation_results_observed",plot_untransformed=T,aligned=T)
unit <- unit_ori
####Third, add survival curve  #######
par(pars)
survival_analysis(experiment="survival_experiment",which_to_plot="second_only")
####Fourth, add letters  #########
par(xpd=NA)
##LETTERS
x_text1 <- grconvertX(0+1/80, from='ndc'); x_text2 <- grconvertX((sum(widz[1:2]))/(sum(widz))+1/80, from='ndc')
y_text1 <- grconvertY((1-1/80), from='ndc')
y_text2 <- grconvertY((1-heits[1]-1/80), from='ndc')
text(x_text1,y_text1,labels=panel_casse("a"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
text(x_text1,y_text2,labels=panel_casse("b"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
text(x_text2,y_text1,labels=panel_casse("c"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
par(xpd=F)
####Fifth, Close figure 3 #######
dev.off()

######Open pdf Figure 4 #####
pdf(file=paste(figurefolder,"/Figure4.pdf",sep=""),family=text_font,font=text_font,bg="white",width=double_col,height=page_height/2.5,pointsize=pointsize_more_than_2row2col)
Extended <- T
##### Set-up layout and plot parameters #####
par(pars)
ncoli <- 13
layout(matrix(c(rep(5,ncoli),
                rep(5,ncoli/13),rep(1,6*ncoli/13),rep(2,6*ncoli/13),
                rep(5,ncoli),
                rep(5,ncoli/13),rep(3,6*ncoli/13),rep(4,6*ncoli/13) 
), 4, ncoli, byrow = TRUE),heights=c(0.05,0.425,0.1,0.425))
####Time spent in nest  ##########
root_path <- paste(disk_path,"/main_experiment",sep="")
queen <- F; treated <- T; nurses <- T; foragers <- T
variable_list <- c("prop_time_outside")
names(variable_list) <- c("prop. of time outside")
transf_variable_list <- c("power0.01")
predictor_list <- "task_group"

analysis <- list(variable_list=variable_list,
                 transf_variable_list=transf_variable_list,
                 predictor_list=predictor_list,
                 violin_plot_param = list(c(1.5,0,-0.02,1,0.1)))
plot_before_vs_after(root_path=root_path, data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep=""),analysis=analysis,status="all",collective=F,pool_plot=F,pattern="individual_behavioural_data",plot_untransformed = T,boldy=T)
######## clean before next step###
clean(); 
Sys.sleep(2)
####Distance to colony,Overlap_with_brood####
root_path <- paste(disk_path,"/main_experiment",sep="")
queen <- F; treated <- T; nurses <- T; foragers <- T
level <- "all"
variable_list <- c("distance_antCoG_to_colonyCoG","BA_between_ant_and_brood")
names(variable_list) <- c(paste("distance to colony","_changetomm",sep=""),"overlap with brood")
predictor_list <- c("task_group","task_group")
transf_variable_list <- c("log","none")

analysis <- list(variable_list=variable_list,
                 transf_variable_list=transf_variable_list,
                 predictor_list = predictor_list,
                 violin_plot_param = list(c(1.5,0,-0.02,1,0.15)))
plot_before_vs_after(root_path=root_path, data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep=""),analysis=analysis,status="all",collective=F,pool_plot=F,pattern="individual_behavioural_data",boldy=T)

######## clean before next step###
clean(); 
Sys.sleep(2)
####Interactions within vs. between, untreated only #####
root_path <- paste(disk_path,"/main_experiment",sep="")
queen <- F; treated <- F; nurses <- T; foragers <- T

variable_list <- c("inter_caste_contact_duration")
names(variable_list) <- c("Inter-task contact (min)")
transf_variable_list <- c("sqrt")
predictor_list <- "task_group"
analysis <- list(variable_list=variable_list,
                 transf_variable_list=transf_variable_list,
                 predictor_list=predictor_list,
                 violin_plot_param = list(c(1.5,0,-0.02,1,0.25),c(1.5,0,-0.02,1,0.25),c(1.5,0,-0.02,1,0.25),c(1.5,0,-0.02,1,0.25)))
plot_before_vs_after(root_path=root_path, data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep=""),analysis=analysis,status="all",collective=F,pool_plot=F,pattern="individual_behavioural_data",boldy=T)
####Fifth, add letters  #########
par(xpd=NA)
##LETTERS
x_text1 <- grconvertX(1/13+1/80, from='ndc'); x_text2 <- grconvertX(7/13+1/80, from='ndc')
y_text1 <- grconvertY((1-1/80), from='ndc')
y_text2 <- grconvertY((1-0.5-1/80), from='ndc')
text(x_text1,y_text1,labels=panel_casse("a"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
text(x_text1,y_text2,labels=panel_casse("c"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
text(x_text2,y_text1,labels=panel_casse("b"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
text(x_text2,y_text2,labels=panel_casse("d"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
par(xpd=F)
####Fifth, Close figure 4 #######
dev.off()
######Open pdf Extended data 2: obs. vs random #####
figure_height <- page_height*0.62
three_col  <- conv_unit(165, "mm","inch")
pdf(file=paste(figurefolder,"/Extended_data_2.pdf",sep=""),family=text_font,font=text_font,bg="white",width=three_col,height=figure_height,pointsize=pointsize_more_than_2row2col)
############## ##############
##### Set-up layout and plot parameters #####
a <- (double_col/2)/figure_height
c <- 0
b <- (1 -a)/2
par(pars)
ncoli <- 6
layout(matrix(c(rep(1,ncoli/2),rep(2,ncoli/2),####networks at the top
                rep(3,ncoli/3),rep(4,ncoli/3),rep(5,ncoli/3),
                rep(6,ncoli/3),rep(7,ncoli/3),rep(8,ncoli/3)
), 3, ncoli, byrow = TRUE),heights=c(a,b,b))
### clean before next step
clean(); 
Sys.sleep(2)
####First plot networks #######
#######chosen case: experiment 2, replicate 59, time_hours=9
root_path <- paste(disk_path,"/main_experiment",sep="")######linux laptop
plot_network(case="topology_comparison",which_to_draw=c("PreTreatment_random","PreTreatment_observed"))
### clean before next step
clean(); 
Sys.sleep(2)
####Second plot random vs observed #######
par_mar_ori <- par()$mar
par_mgp_ori <- par()$mgp
par(mar=par_mar_ori+c(0,1,0,1))
par(mgp=par_mgp_ori+c(0.1,0,0))
variable_list <- c("modularity","density","degree_mean", 
                   "task_assortativity","diameter","degree_maximum")
names(variable_list) <- c("Modularity","Density", "Mean degree",
                          "Task assortativity","Diameter","Max. degree")
full_statuses_names_ori <- full_statuses_names
full_statuses_names[full_statuses_names%in%c("Rand.","Obs.")] <- c("Random","Observed")
plot_observed_vs_random(experiments="main_experiment",variables=variable_list,pattern="network_properties",data_path="processed_data/network_properties/random_vs_observed")
par(mar=par_mar_ori)
par(mgp=par_mgp_ori)
full_statuses_names <- full_statuses_names_ori
### clean before next step

clean(); 
Sys.sleep(2)
####Third, add letters  #########
par(xpd=NA)
##LETTERS
x_text1 <- grconvertX(0+1/80, from='ndc'); 
y_text1 <- grconvertY((1-1/80), from='ndc')
y_text2 <- grconvertY((1-a-1/80), from='ndc')
text(x_text1,y_text1,labels=panel_casse("a"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
text(x_text1,y_text2,labels=panel_casse("b"),font=panel_font, cex=panel_cex,adj=c(0,0),xpd=NA)
par(xpd=F)
######## clean before next step###

clean(); 
Sys.sleep(2)
####Fourth, Close Extended Data 2 #######
dev.off()
######Extended data 3: simulation results from different seeds, obs vs. random, COLLECTIVE #####
pdf(file=paste(figurefolder,"/Extended_data_3.pdf",sep=""),family=text_font,font=text_font,bg="white",width=  three_col,height=page_height*0.4,pointsize=pointsize_more_than_2row2col)
seeds <- c("random_workers","frequent_foragers","occasional_foragers","nurses","low_degree","high_degree")
names(seeds) <- c("Seeds = random workers","Seeds = frequent foragers","Seeds = occasional foragers","Seeds = nurses","Seeds = low degree workers (LD)","Seeds = high degree workers (HD)" )
seeds <- seeds[c(2,3,4,1)]
nrep <- length(seeds)/2
ncols <- 7
d <- 0.05
variable_list <- c("logistic_r","Mean_load","Load_skewness")
names(variable_list) <- c("Transmission rate","Mean simulated load (W)","Simulated load skewness (W)")

par(pars)
par_mar_ori <- par()$mar
par(mar=par_mar_ori+c(-1,0,0,0))
layout(
  matrix(c(
    rep(13,ncols),
    c(1:3),13,c(4:6),
    rep(13,ncols),
    c(7:9),13,c(10:12)
  )
  ,
  nrow=nrep+2,
  ncol=ncols,
  byrow = T
  )
  ,widths=c(rep((1-d)/6,3),d,rep((1-d)/6,3))
  ,heights=c(d,(1-2*d)/2,d,(1-2*d)/2)
)

for (seed in seeds){
  full_statuses_names_ori <- full_statuses_names
  full_statuses_names[full_statuses_names%in%c("Rand.","Obs.")] <- c("Rd","Obs")
  print(seed)
  plot_observed_vs_random(experiments="main_experiment",variables=variable_list,data_path = paste("transmission_simulations/random_vs_observed/",seed,sep=""),pattern="collective_simulation_results_")
  
  
  ####plot title
  par(xpd=NA)
  x_text <- grconvertX(1/4+as.numeric(((which(seed==seeds)/2)-floor((which(seed==seeds)/2)))==0)/2, from='ndc')
  y_text <- grconvertY((1 - 2*d/3 - (ceiling((which(seed==seeds)/2))-1)*(d+(1-2*d)/2)), from='ndc')
  
  text(x_text,y_text,labels=names(seeds[seeds==seed]),font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)
  # text(x_text,y_text,labels="The observed pre-exposure networks protect colonies from outside pathogens",font=2, cex=max_cex,adj=0.5,xpd=NA)
  
  ##if necessary, plot line
  if (is.even(which(seed==seeds))){
    
    x_line <- grconvertX(1/2, from='ndc')
    y_line1 <- grconvertY((1 - d - (ceiling((which(seed==seeds)/2))-1)*(d+(1-2*d)/2)), from='ndc') ####lower left = 0-0;top_right=1-1
    y_line2 <- grconvertY((1 - d - (1-2*d)/2 - (ceiling((which(seed==seeds)/2))-1)*(d+(1-2*d)/2)), from='ndc')
    segments(x_line,y_line1,x_line,y_line2)
  }
  par(xpd=F)
  full_statuses_names <- full_statuses_names_ori
}
par(mar=par_mar_ori)
dev.off()
######## clean before next step###
clean(); 
Sys.sleep(2)

######Open pdf Extended data 4: degree and distance to queen vs time spent outside in observed networks #####
figure_height <- three_col/2.2
pdf(file=paste(figurefolder,"/Extended_data_4.pdf",sep=""),family=text_font,font=text_font,bg="white",width=three_col,height=figure_height,pointsize=pointsize_less_than_2row2col)
##### Set-up layout and plot parameters #####
par(pars)
par_mar_ori <- par()$mar
par_mgp_ori <- par()$mgp
par(mar=par_mar_ori+c(0.2,0,0.8,1))
par(mgp=par_mgp_ori+c(0.3,0,0))

ncoli <- 2
layout(matrix(c(1,2
), 1, ncoli, byrow = TRUE))

###read data
data           <- read.table(paste(disk_path,"/main_experiment/processed_data/individual_behaviour/pre_treatment/network_position_vs_time_outside.dat",sep=""),header=T,stringsAsFactors = F)
data["antid"]  <- as.character(interaction(data$colony,data$tag))
####list desired variables and transformations
variable_list <- c("degree","aggregated_distance_to_queen")
names(variable_list) <- c("degree","path length to queen")
predictor_list <- c("prop_time_outside","prop_time_outside")
names(predictor_list) <- c("Prop. of time outside","Prop. of time outside")
transf_variable_list <- c("none","log")
transf_predictor_list <- c("power2","none")
analysis <- list(variable_list=variable_list,
                 predictor_list=predictor_list,
                 transf_variable_list=transf_variable_list,
                 transf_predictor_list=transf_predictor_list,              
                 violin_plot_param = list(c(1,0,-0.025,0.07,5),c(1,0,-0.025,0.07,0.2)))
####plot
plot_regression(data=data,time_point="before",analysis=analysis,n_cat_horiz=15,n_cat_vertic=30,pool=c(F,F),collective=T)
par(mar=par_mar_ori)
par(mgp=par_mgp_ori)
######## clean before next step###

clean(); 
Sys.sleep(2)
####Second, add letters  #########
par(xpd=NA)
##LETTERS
x_text1 <- grconvertX(0+1/80, from='ndc'); 
x_text2 <- grconvertX(0.5+1/80, from='ndc'); 
y_text1 <- grconvertY((1-1/80), from='ndc')
text(x_text1,y_text1,labels=panel_casse("a"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
text(x_text2,y_text1,labels=panel_casse("b"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
par(xpd=F)
######## clean before next step###
clean(); 
Sys.sleep(2)
####Third, Close Extended Data 4 #######
dev.off()


######Extended data 5: simulation results from different seeds, observed #####
seeds <- c("random_workers","low_degree","high_degree","frequent_foragers","occasional_foragers","nurses")
names(seeds) <- c("R","LD","HD","FF","OF","N")
seeds <- seeds[c(2,3,4,5,6,1)]
color_pal <- c(GetColorHex("grey50"),GetColorHex("grey70"),GetColorHex("grey20"),statuses_colours["forager"],statuses_colours["occasional_forager"],statuses_colours["nurse"])
color_pal <- color_pal[c(2,3,4,5,6,1)]
names(color_pal) <- seeds

variables <- c("logistic_r","Prevalence","Mean_load","Load_skewness","Queen_load")
names(variables) <- c("Transmission rate","Prevalence","Mean simulated load (W)","Simulated load skewness (W)","Simulated load (Q)")
transf <- c("none","none","none","none","none")
names(transf) <- variables

pdf(file=paste(figurefolder,"/Extended_data_5.pdf",sep=""),family=text_font,font=text_font,bg="white",width=three_col,height=page_height*0.49,pointsize=pointsize_more_than_2row2col)
par(pars)
par(mar=c(2.6,2,1,1))
layout(matrix(c(1:(2*ceiling(length(variables)/2))), ceiling(length(variables)/2), 2, byrow = T))
plot_seeds(experiments="main_experiment",seeds=seeds,variables=variables,transf=transf,color_pal=color_pal)
dev.off()

######## clean before next step###

clean(); 
Sys.sleep(2)
######Extended data 6: age-and-caste organization, main experiment + age-marked experiment #####
full_statuses_names_ori <- full_statuses_names

full_statuses_names[full_statuses_names%in%c("Rand.","Obs.")] <- c("Rd","Obs")
full_statuses_names[full_statuses_names%in%c("Q\ncomm." ,"Other\ncomm." )] <- c("Q","W")
pdf(file=paste(figurefolder,"/Extended_data_6.pdf",sep=""),family=text_font,font=text_font,bg="white",width=three_col,height=page_height*0.55,pointsize=pointsize_more_than_2row2col)
til <- 0.1
par(pars)
par_mar_ori <- par("mar")
# par(mar=par_mar_ori+c(1,0,0,0))
ncols <- 6
nrows <- 3
heits <- c(til,til,3/4*(1-til),til,1-til,0,3/4*(1-til))
layout(matrix (
  c (rep(14,ncols),
     rep(14,ncols),
     rep(1,ncols/3),rep(8,ncols/6),rep(9,ncols/6),rep(13,ncols/3),
     rep(14,ncols),
     rep(2,ncols/6),rep(3,ncols/6), rep(4,ncols/6),rep(5,ncols/6), rep(10,ncols/6),rep(11,ncols/6),
     rep(14,ncols),
     rep(6,ncols/3), rep(7,ncols/3), rep(12,ncols/3)
  )
  ,
  nrow=1+2*nrows,
  ncol=ncols,
  byrow = T
  
),
heights=heits
)

plot_age_dol(experiments=c("age_experiment","main_experiment"))
par(mar=par_mar_ori)
dev.off()
full_statuses_names <- full_statuses_names_ori 
######## clean before next step###
clean();
Sys.sleep(2)
######Extended data 7: simulation results from different seeds, obs vs. random, INDIVIDUAL #####
pdf(file=paste(figurefolder,"/Extended_data_7.pdf",sep=""),family=text_font,font=text_font,bg="white",width=  three_col,height=page_height*0.4,pointsize=pointsize_more_than_2row2col)
root_path <- paste(disk_path,"/main_experiment",sep="")######linux laptop
seeds <- c("random_workers","frequent_foragers","occasional_foragers","nurses","low_degree","high_degree")
names(seeds) <- c("Seeds = random workers","Seeds = frequent foragers","Seeds = occasional foragers","Seeds = nurses","Seeds = low degree workers (LD)","Seeds = high degree workers (HD)" )
seeds <- seeds[c(2,3,4,1)]
nrep <- length(seeds)/2
ncols <- 5
d <- 0.04
variable_list <-  c("transmission_latency","simulated_load")
names(variable_list) <- c("Contamination latency (h)","Simulated load")
predictor_list <- c("task_group","task_group")
transf_variable_list <- c("none","sqrt")

par(pars)
layout(
  matrix(c(
    rep(9,ncols),
    c(1:2),9,c(3:4),
    rep(9,ncols),
    c(5:6),9,c(7:8)
  )
  ,
  nrow=nrep+2,
  ncol=ncols,
  byrow = T
  )
  ,widths=c(rep((1-d)/4,2),d,rep((1-d)/4,2))
  ,heights=c(d,(1-2*d)/2,d,(1-2*d)/2)
)

unit_ori <- unit; unit <- 24
time_window <- 24
full_statuses_names_ori <- full_statuses_names
full_statuses_names[full_statuses_names%in%c("Rand.","Obs.")] <- c("Random","Observed")
full_statuses_names[full_statuses_names%in%c("Queen","Queen\n")] <- "Q"
full_statuses_names[full_statuses_names%in%c("Nurses\n","Untreated\nnurses")] <- "N"
full_statuses_names[full_statuses_names%in%c("Foragers","Untreated\nforagers")] <- "F"


for (seed in seeds){
  print(seed)
  ####then plot, for each of observed and random, the simulated time to infection and simulated mean load of queen, nurses and foragers
  
  queen <- T; treated <- F; nurses <- T; foragers <- T
  
  analysis <- list(variable_list=variable_list,
                   transf_variable_list=transf_variable_list,
                   predictor_list=predictor_list,
                   violin_plot_param = list(c(1.5,0,-0.02,0.25,0.15),c(1.5,0,-0.02,0.25,0.15)))
  plot_before_vs_after(root_path=root_path, data_path=paste(root_path,"/transmission_simulations/random_vs_observed/",seed,sep=""),analysis=analysis,status="all",collective=F,pool_plot=F,pattern="summarised_individual_results.dat")
  
  ####plot title
  par(xpd=NA)
  x_text <- grconvertX(1/4+as.numeric(((which(seed==seeds)/2)-floor((which(seed==seeds)/2)))==0)/2, from='ndc')
  y_text <- grconvertY((1 - 2*d/3 - (ceiling((which(seed==seeds)/2))-1)*(d+(1-2*d)/2)), from='ndc')
  
  text(x_text,y_text,labels=names(seeds[seeds==seed]),font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)
  # text(x_text,y_text,labels="The observed pre-exposure networks protect colonies from outside pathogens",font=2, cex=max_cex,adj=0.5,xpd=NA)
  
  ##if necessary, plot line
  if (is.even(which(seed==seeds))){
    
    x_line <- grconvertX(1/2, from='ndc')
    y_line1 <- grconvertY((1 - d - (ceiling((which(seed==seeds)/2))-1)*(d+(1-2*d)/2)), from='ndc') ####lower left = 0-0;top_right=1-1
    y_line2 <- grconvertY((1 - d - (1-2*d)/2 - (ceiling((which(seed==seeds)/2))-1)*(d+(1-2*d)/2)), from='ndc')
    segments(x_line,y_line1,x_line,y_line2)
  }
  par(xpd=F)
  
}

unit <- unit_ori
full_statuses_names <- full_statuses_names_ori 

dev.off()
######## clean before next step###

clean(); 
Sys.sleep(2)
######Extended data 8: measured load on networks #####
pdf(file=paste(figurefolder,"/Extended_data_8.pdf",sep=""),family=text_font,font=text_font,bg="white",width=three_col,height=1.3*(three_col/3),pointsize=pointsize_more_than_2row2col)
par(pars)
layout(matrix(c(1,2,4,6,3,5),nrow=2,byrow = T),heights=c(0.7,0.2))
root_path <- paste(disk_path,"/main_experiment",sep="")######linux laptop
plot_network(case="node_color_f_property",which_to_draw=c("PostTreatment_observed"),colours=c("task_group","measured_load","simulated_load"))
par(mar=c(0.9,0,0,0))
plot(1:3,rep(1,3),pch=NA,bty='n',xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,3.5),ylim=c(0,3))
points(1:3,rep(1.5,3),pch=21,bg=statuses_colours[c("queen","nurse","forager")],cex=1.3,lwd=line_min)
text(1:3,rep(0.35,3),labels=c("Queen","Nurse","Forager"),cex=inter_cex,adj=0.5)
par(xpd=NA)
##LETTERS
x_text1 <- grconvertX(0+1/80, from='ndc'); x_text2 <- grconvertX(1/3+1/80, from='ndc');  x_text3 <- grconvertX(2/3++1/80, from='ndc')
y_text <- grconvertY((1-1/80), from='ndc')
text(x_text1,y_text,labels=panel_casse("a"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
text(x_text2,y_text,labels=panel_casse("b"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
text(x_text3,y_text,labels=panel_casse("c"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
par(xpd=F)
dev.off()
######Extended data 9: qPCR vs network properties ############
pdf(file=paste(figurefolder,"/Extended_data_9.pdf",sep=""),family=text_font,font=text_font,bg="white",width=  three_col,height=page_height*0.6,pointsize=pointsize_2row2col)
layout(matrix(c(3,1,4,2),byrow=T,nrow=2),heights=c(3,2))
partial_least_square_regression(experiments=c("age_experiment","main_experiment"))
x_text1 <- grconvertX(1/80, from='ndc');x_text2 <- grconvertX(0.5+1/80, from='ndc')
y_text1 <- grconvertY(1-2/80, from='ndc')
y_text2<- grconvertY(1-1/2-1/80, from='ndc')
y_text3<- grconvertY(1-3/5-1/80, from='ndc')
text(x_text1,y_text1,labels=panel_casse("a"),font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)
text(x_text1,y_text3,labels=panel_casse("b"),font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)

dev.off()
clean()
######Extenced Data 10: Distance to treated ############
pdf(file=paste(figurefolder,"/Extended_Data_10.pdf",sep=""),family=text_font,font=text_font,bg="white",width=three_col,height=page_height*0.7,pointsize=pointsize_more_than_2row2col)
full_statuses_names_ori <- full_statuses_names
full_statuses_names[full_statuses_names%in%c("Foragers","Untreated\nforagers")] <- "Foragers\n"
statuses_colours_ori <- statuses_colours
statuses_colours[names(statuses_colours)%in%c("queen","outdoor_ant","not_outdoor_ant")] <- "black"

plot_type_ori <- plot_type
plot_type <- "means"
par(pars)
par_mar_ori <- par()$mar
par(mar=par_mar_ori+c(1,0,1,1))
heits <- c(1,10,1,10,1,10)
layout(matrix(c(7,1,7,3,7,5,7,2,7,4,7,6),nrow=6),heights=heits)
plot_qpcr_vs_distance_to_treated(experiments=c("age_experiment","main_experiment"))
####Add letters 
par(xpd=NA)
x_text1 <- grconvertX(1/80, from='ndc');x_text2 <- grconvertX(0.5+1/80, from='ndc')
y_text1 <- grconvertY(1-1/80, from='ndc')
y_text2<- grconvertY(1-sum(heits[1:2])/sum(heits)-1/80, from='ndc')
y_text3 <- grconvertY(1-sum(heits[1:4])/sum(heits)-1/80, from='ndc')
text(x_text1,y_text1,labels=panel_casse("a"),font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)
# text(x_text2,y_text1,labels=panel_casse("b"),font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)
text(x_text1,y_text2,labels=panel_casse("b"),font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)
# text(x_text2,y_text2,labels=panel_casse("d"),font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)
text(x_text1,y_text3,labels=panel_casse("c"),font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)
# text(x_text2,y_text3,labels=panel_casse("f"),font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)

x_text2 <- grconvertX(0.5, from='ndc')
text(x_text2,y_text1,labels="All workers",font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)
text(x_text2,y_text2,labels="Nurses",font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)
text(x_text2,y_text3,labels="Low contact with treated",font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)

par(xpd=F)
plot_type <- plot_type_ori
full_statuses_names <- full_statuses_names_ori
statuses_colours <- statuses_colours_ori
par(mar=par_mar_ori)
dev.off()
clean()


######Extended data 11: survival experiment #####
pdf(file=paste(figurefolder,"/Extended_data_11.pdf",sep=""),family=text_font,font=text_font,bg="white",width=three_col,height=page_height*0.54,pointsize=pointsize_more_than_2row2col)
  par(pars)
  layout(
    matrix(
      c(
        1,5,2,
        rep(5,3),
        3,5,4
      )
      ,nrow=3, byrow=T)
    ,heights=c(0.45,0.01,0.45),widths=c(0.45,0.03,0.45)
  )
  
  survival_analysis(experiment="survival_experiment",which_to_plot = "detailed")
  par(xpd=NA)
  ##LETTERS
  x_text1 <- grconvertX(0+1/80, from='ndc'); 
  x_text2 <- grconvertX((0.45+0.05)/(0.45+0.05+0.45)+1/80, from='ndc'); 
  y_text1 <- grconvertY((1-1/80), from='ndc')
  y_text2 <- grconvertY((1-((0.45+0.01)/(0.45+0.01+0.45))-1/80), from='ndc')
  
  text(x_text1,y_text1,labels=panel_casse("a"),font=panel_font, cex=panel_cex,xpd=NA)
  text(x_text2,y_text1,labels=panel_casse("b"),font=panel_font, cex=panel_cex,xpd=NA)
  text(x_text1,y_text2,labels=panel_casse("c"),font=panel_font, cex=panel_cex,xpd=NA)
  text(x_text2,y_text2,labels=panel_casse("d"),font=panel_font, cex=panel_cex,xpd=NA)
  par(xpd=F)
 dev.off()
######## clean before next step###
clean();
Sys.sleep(2)

######Extended data 12: pathogen-induced network changes; individual node properties #####
pdf(file=paste(figurefolder,"/Extended_data_12.pdf",sep=""),family=text_font,font=text_font,bg="white",width=three_col,height=page_height/5,pointsize=pointsize_more_than_2row2col)
par(pars)
ncoli <- 13
layout(matrix(c(rep(3,ncoli),
                rep(3,ncoli/13),rep(1,6*ncoli/13),rep(2,6*ncoli/13)
), 2, ncoli, byrow = TRUE),heights=c(0.05,0.45))

queen <- F; treated <- T; nurses <- T; foragers <- T
root_path <- paste(disk_path,"/main_experiment",sep="")
variable_list <- c("degree","aggregated_distance_to_queen")
names(variable_list) <- c("degree","path length to queen")
transf_variable_list <- c("none","log")
predictor_list <- c("task_group","task_group")
analysis <- list(variable_list=variable_list,
                 transf_variable_list=transf_variable_list,
                 predictor_list=predictor_list,
                 violin_plot_param = list(c(1.5,0,-0.02,0.25,0.15),c(1.5,0,-0.02,0.25,0.15)))
plot_before_vs_after(root_path=root_path, data_path=paste(root_path,"/processed_data/network_properties/pre_vs_post_treatment/all_workers",sep=""),analysis=analysis,status="all",collective=F,pool_plot=F,pattern="individual_data",boldy=T)
dev.off()
######## clean before next step###
clean();
Sys.sleep(2)
######Extended data 13: pathogen-induced simulation changes; colony-level plot#####
pdf(file=paste(figurefolder,"/Extended_data_13.pdf",sep=""),family=text_font,font=text_font,bg="white",width=three_col,height=page_height/5,pointsize=pointsize_more_than_2row2col)
unit_ori <- unit; unit <- 24
time_window <- 24
par(pars)
ncoli <- 13
layout(matrix(c(rep(5,ncoli),
                rep(5,ncoli/13),rep(1,3*ncoli/13),rep(2,3*ncoli/13),rep(3,3*ncoli/13),rep(4,3*ncoli/13)
), 2, ncoli, byrow = TRUE),heights=c(0.05,0.45))
root_path <- paste(disk_path,"/main_experiment",sep="")
variable_list <- c("probability_of_transmission")
names(variable_list) <- c("prob. transmission")
transf_variable_list <- c("power8")
analysis <- list(variable_list=variable_list,
                 transf_variable_list=transf_variable_list,
                 violin_plot_param = list(c(1.5,0,-0.02,0.25,0.15),c(1.5,0,-0.02,0.25,0.15)))
plot_before_vs_after(root_path=root_path, data_path=paste(root_path,"/transmission_simulations/pre_vs_post_treatment/experimentally_exposed_seeds",sep=""),analysis=analysis,status="untreated",collective=F,pool_plot=F,pattern="individual_simulation_results",boldy=T,plot_untransformed = T)


variable_list <- c("simulated_load","probability_high_level","probability_low_level")
names(variable_list) <- c("simulated load","prob. high load","prob. low load")
transf_variable_list <- c("sqrt","none","none")
analysis <- list(variable_list=variable_list,
                 transf_variable_list=transf_variable_list,
                 violin_plot_param = list(c(1.5,0,-0.02,0.25,0.15),c(1.5,0,-0.02,0.25,0.15)))
plot_before_vs_after(root_path=root_path, data_path=paste(root_path,"/transmission_simulations/pre_vs_post_treatment/experimentally_exposed_seeds",sep=""),analysis=analysis,status="untreated",collective=F,pool_plot=F,pattern="individual_simulation_results",boldy=T)
unit <- unit_ori
dev.off()
######## clean before next step###
clean();
Sys.sleep(2)

######Extended data 14: treatment-induced changes in worker behaviour #####
# pdf(file=paste(figurefolder,"/Extended_data_14.pdf",sep=""),family=text_font,font=text_font,bg="white",width=three_col,height=page_height*0.4,pointsize=pointsize_more_than_2row2col)
par(pars)
par(mar=par()$mar+c(0,0,0,0))
ncoli <- 10
heits <- c(0.08,0.45,0.08,0.45)
heits <- heits/sum(heits)
d <- 0.06
to_keep <- c(to_keep,"heits")
layout(matrix(c(
  rep(7,ncoli),
  rep(7,ncoli/10),rep(1,3*ncoli/10),rep(2,3*ncoli/10),rep(3,3*ncoli/10),
  rep(7,ncoli),
  rep(7,ncoli/10),rep(4,3*ncoli/10),rep(5,3*ncoli/10),rep(6,3*ncoli/10)
), 4, ncoli, byrow = TRUE),heights=heits,widths=c(d,rep((1-d)/9,9)))

experiment <- "main_experiment"
####1. individual behaviour ######
root_path <- paste(disk_path,experiment,sep="/")
queen <- F; treated <- T; nurses <- T; foragers <- T
variable_list <- c("proportion_time_active","average_bout_speed_pixpersec","total_distance_travelled_pix","within_nest_home_range","duration_of_contact_with_treated_min")
names(variable_list) <- c("Prop. time active","Speed while active (mm/sec)_changetomm","Total distance moved (mm)_changetomm","area visited within the nest_changetomm2","contact with treated (min)")
transf_variable_list <- c("none","log","sqrt","sqrt","log")
predictor_list <- rep("task_group",length(variable_list))
analysis <- list(variable_list=variable_list,
                 transf_variable_list=transf_variable_list,
                 predictor_list=predictor_list,
                 violin_plot_param = list(c(1.5,0,-0.02,0.25,0.15),c(1.5,0,-0.02,0.25,0.15)))
plot_before_vs_after(root_path=root_path, data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep=""),analysis=analysis,status="all",collective=F,pool_plot=F,pattern="individual_behavioural_data",boldy=T)

###2. Brood displacement #####
root_path <- paste(disk_path,experiment,sep="/")

variable_list <- c("distance_from_nest_entrance_pix")
names(variable_list) <- c("brood location (mm)_changetomm")
transf_variable_list <- c("none")

analysis <- list(variable_list=variable_list,
                 transf_variable_list=transf_variable_list,
                 violin_plot_param = list(c(1.5,0,-0.02,0.25,0.15),c(1.5,0,-0.02,0.25,0.15),c(1.5,0,-0.02,0.25,0.15),c(1.5,0,-0.02,0.25,0.15)))

plot_before_vs_after(root_path=root_path, data_path=paste(root_path,"/processed_data/collective_behaviour/pre_vs_post_treatment",sep=""),analysis=analysis,status="collective",collective=T,pool_plot=T,pattern="brood_location.txt",boldy=T)

clean();
####3. Add letters ####
par(xpd=NA)
x_text1 <- grconvertX(1/80, from='ndc');x_text2 <- grconvertX(d+3*(1-d)/9+1/80, from='ndc');x_text2bis <- grconvertX(d+2*(1-d)/9, from='ndc');x_text3 <- grconvertX(d+6*(1-d)/9+1/80, from='ndc')
y_text1 <- grconvertY(1-1/80, from='ndc');y_text2 <- grconvertY(1-1/80-sum(heits[1:2]), from='ndc');
text(x_text1,y_text1,labels=panel_casse("a"),font=2, cex=max_cex,adj=c(0.5,1),xpd=NA)
text(x_text2,y_text1,labels=panel_casse("b"),font=2, cex=max_cex,adj=c(0.5,1),xpd=NA)
text(x_text3,y_text1,labels=panel_casse("c"),font=2, cex=max_cex,adj=c(0.5,1),xpd=NA)
text(x_text1,y_text2,labels=panel_casse("d"),font=2, cex=max_cex,adj=c(0.5,1),xpd=NA)
text(x_text2,y_text2,labels=panel_casse("e"),font=2, cex=max_cex,adj=c(0.5,1),xpd=NA)
text(x_text3,y_text2,labels=panel_casse("f"),font=2, cex=max_cex,adj=c(0.5,1),xpd=NA)
par(xpd=F)
####7. Close pdf ####
# dev.off()