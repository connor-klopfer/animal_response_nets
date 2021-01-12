#### Import data Lattanzio & Miles paper ####

library(readxl)
LattanzioMiles_morph_network <- read_excel("~/Postdoc_UC/CNWW workshop/LattanzioMiles/LattanzioMilesJAEdata.xlsx", 
                                    sheet = "morph_networks")

LattanzioMiles_contests <- read_excel("~/Postdoc_UC/CNWW workshop/LattanzioMiles/LattanzioMilesJAEdata.xlsx", 
                                           sheet = "contests")

LattanzioMiles_strength <- read_excel("~/Postdoc_UC/CNWW workshop/LattanzioMiles/LattanzioMilesJAEdata.xlsx", 
                                      sheet = "strength_data")

LattanzioMiles_NBdiet <- read_excel("~/Postdoc_UC/CNWW workshop/LattanzioMiles/LattanzioMilesJAEdata.xlsx", 
                                     sheet = "NBdiet")

LattanzioMiles_LBdiet <- read_excel("~/Postdoc_UC/CNWW workshop/LattanzioMiles/LattanzioMilesJAEdata.xlsx", 
                                     sheet = "LBdiet")

LattanzioMiles_HBdiet <- read_excel("~/Postdoc_UC/CNWW workshop/LattanzioMiles/LattanzioMilesJAEdata.xlsx", 
                                     sheet = "HBdiet")


