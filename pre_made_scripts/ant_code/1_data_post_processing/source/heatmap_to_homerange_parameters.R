euc.dist <- function(A, B) {
  return( sqrt(((A$x-B$x)^2)+((A$y-B$y)^2)))
}
###analysis parameters#######
grid_parameter <- 6 ####the lower the number, the finer the grid (e.g. if fine_grid_parameter=4 the program will regroup 16 pixels into 1; if fine_grid_parameter=2 the program will regroup 4 pixels into 1 , etc. ), the longer the processing time
###top left corner coordinates
tlc_x <- -1
tlc_y <- -1
###top right corner coordinates
trc_x <- 609
trc_y <- -1
###bottom right corner coordinates (remember: Y is positive towards the bottom)
brc_x <- 609
brc_y <- 912
###bottom left corner coordinates
blc_x <- -1
blc_y <- 912

#####read infos
########get information of where the nest is
ynestmin_px <- info[info$colony==colony_number,"Ynestmin"]
ynestmax_px <- info[info$colony==colony_number,"Ynestmax"]
if (ynestmin_px==0){
  ynestmin <- 0
  ynestmax <- ceiling(ynestmax_px/5)
  yforagingmin <- ynestmax + 1
  yforagingmax <- 911
  nest_entrance <- data.frame(x=(1+608/2),y=ynestmax)
}else{
  yforagingmin <- 0
  yforagingmax <- floor(ynestmin_px/5) - 1
  ynestmin <- floor(ynestmin_px/5)
  ynestmax <- 911
  nest_entrance <- data.frame(x=(1+608/2),y=ynestmin)
}
yallmin <- 0
yallmax <- 911
xallmin <- 0
xallmax <- 608
######### define boundaries to limit estimation of utilisation distributions at the nest walls- See Benhamou & Cornelis, 2010
######all
bound_all <- structure(list(X=c(tlc_x ,trc_x,brc_x ,blc_x ,tlc_x),Y=c(tlc_y ,trc_y,brc_y ,blc_y ,tlc_y)), .Names = c("X", "Y"))   ## Prep for boundary of HR area
bound_all <- do.call("cbind",bound_all)
Slo1    <- Line(bound_all)
Sli1    <- Lines(list(Slo1), ID="frontier1")
bounder_all <- SpatialLines(list(Sli1))
######nest
bound_nest <- structure(list(X=c(tlc_x ,trc_x,brc_x ,blc_x ,tlc_x),Y=c(ynestmin-1 ,ynestmin-1,ynestmax+1 ,ynestmax +1,ynestmin-1)), .Names = c("X", "Y"))   ## Prep for boundary of HR area
bound_nest <- do.call("cbind",bound_nest)
Slo1    <- Line(bound_nest)
Sli1    <- Lines(list(Slo1), ID="frontier1")
bounder_nest <- SpatialLines(list(Sli1))
######foraging
bound_foraging <- structure(list(X=c(tlc_x ,trc_x,brc_x ,blc_x ,tlc_x),Y=c(yforagingmin -1 ,yforagingmin-1 ,yforagingmax +1 ,yforagingmax +1 ,yforagingmin-1 )), .Names = c("X", "Y"))   ## Prep for boundary of HR area
bound_foraging <- do.call("cbind",bound_foraging)
Slo1    <- Line(bound_foraging)
Sli1    <- Lines(list(Slo1), ID="frontier1")
bounder_foraging <- SpatialLines(list(Sli1))

span <- expand.grid(X=0:608,Y=0:911,obs=0)
overlap_methods <- c("BA")
DENSITY_95 <- F
fixed_bandwidth <- 15
