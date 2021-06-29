#ASNWG search terms using litsearchr
#June 28,29 2021

#install packages

install.packages("litsearchr")
install.packages("revtools")
library(litsearchr)


###########################
##THE NAIVE SEARCH STRING##
###########################

## Created a naive search string using terms from three different themes(from PECO framework)
## Population(Animal groups) + Exposure(Pertubation) + Outcome(change in the network shape)


((primate* OR animal* OR mammal*) OR (insect* OR bird* OR reptile OR bat* OR fish* OR amphibian*) ) 

AND (fire OR disturbance* OR sick OR pathog* OR add* OR remov* OR disrupt* OR manipu*) 

AND ((dynamic OR social OR behavio*r* OR spatial OR colon* ) NEAR/2 (Network OR connect*) AND (plasticit* OR change* OR reduce* OR increase* OR association* OR shape* OR cohesiv* OR adapt* OR segregat* OR diffusion OR structure) ) 


#############################
#Import naive search results#
#############################

library(litsearchr)
library(revtools)

setwd("C://Users//alexp//OneDrive - University of New Brunswick//Grad Studies//Complex Networks Winter Workshop//R_code")

#import articles from the naive search string
ASNWG_import <-
  import_results(directory = "./naive_results/", verbose = FALSE)


#remove duplicates
ASNWG_data1 <- 
  litsearchr::remove_duplicates(ASNWG_import, field = "title", method = "exact") 
ASNWG_data <-
  litsearchr::remove_duplicates(ASNWG_data1, field = "title", method = "string_osa")
ASNWG_data1$title
ASNWG_data$title

#view in a CSV format
###colnames(ASNWG_data1)
write.csv(ASNWG_data1$title, "./remove_dup_1")
write.csv(ASNWG_data$title,"./remove_dup_2")

#extract potential keywords from 