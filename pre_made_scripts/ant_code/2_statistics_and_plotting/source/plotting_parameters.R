########Science-specific plotting parameters #####
####pdf width in inches
single_col <- conv_unit(55, "mm","inch")
double_col <- conv_unit(120, "mm","inch")
three_col  <- conv_unit(183, "mm","inch")
page_height <- conv_unit(297, "mm","inch")

####Text size
####"In a layout with exactly two rows and columns the base value of “cex” is reduced by a factor of 0.83: if there are three or more of either rows or columns, the reduction factor is 0.66.” 
ref_cex <- 11
pointsize_less_than_2row2col <- ref_cex
pointsize_2row2col <- ref_cex/0.83
pointsize_more_than_2row2col <- ref_cex/0.66
panel_cex <- ref_cex/ref_cex
max_cex <- (ref_cex-1)/ref_cex
inter_cex <- (ref_cex-2)/ref_cex
min_cex <- (ref_cex-3)/ref_cex

####Text font
text_font <- "Helvetica"
panel_font <- 2 ####panel letters must be bold upright
panel_casse <- function(x){
  x <- tolower(x)
  x <- capitalize(x)
  return(x)
  
}
casse <- function(x){
  ###remove final point
  if (substr(x,nchar(x),nchar(x))=="."){
    x <- substr(x,1,(nchar(x)-1))
  }
  x <- tolower(x)
  x <- capitalize(x)
  return(x)
}

####Line sizes
line_max <- 1
line_min <- 0.5
line_inter <- line_min +0.5*(line_max-line_min)

if(!file.exists(figurefolder)){dir.create(figurefolder,recursive = T)}