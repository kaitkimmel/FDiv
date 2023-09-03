################################
#### FD METRIC CALCULATIONS ####
################################

# CREATED BY: KAITLIN KIMMEL

# load libraries
library(FD)
library(here)
source("http://www.sthda.com/upload/rquery_cormat.r")



## LOAD DATA

knz.trait <- read.csv(here("data/Cleaned/knz_traitsout.csv"))
knz.ann <- read.csv(here("data/Cleaned/knzannun_commout.csv"))
knz.twe <- read.csv(here("data/Cleaned/knztweun_commout.csv"))
knz.anngr <- read.csv(here("data/Cleaned/knzanngr_commout.csv"))
knz.twegr <- read.csv(here("data/Cleaned/knztwegr_commout.csv"))

### CHECK TRAIT COVERAGE 
# "Before analysis, we will remove species that have less than 100% trait coverage. 
# We will, however, make sure that the communities are still represented by at 
# least 80% of species abundance" 

# remove species with > 100% trait coverage
traits.comp <- knz.trait[complete.cases(knz.trait),] #2 species without all trait data
traits.comp$sp_name <- tolower(traits.comp$sp_name)
# match traits with species in communities
knz.annsub <- knz.ann[,which(names(knz.ann) %in% traits.comp$sp_name)]
knz.annsub$comp.tot <- rowSums(knz.annsub)
min(knz.annsub$comp.tot) # min is 44.89
max(knz.annsub$comp.tot) # max is 77.41

knz.anngrsub <- knz.anngr[,which(names(knz.anngr) %in% traits.comp$sp_name)]
knz.anngrsub$comp.tot <- rowSums(knz.anngrsub)
min(knz.anngrsub$comp.tot) # min is 21.33
max(knz.anngrsub$comp.tot) # max is 55.58

knz.twesub <- knz.twe[,which(names(knz.twe) %in% traits.comp$sp_name)]
knz.twesub$comp.tot <- rowSums(knz.twesub)
min(knz.twesub$comp.tot) # min is 8.53
max(knz.twesub$comp.tot) # max is 32.52

knz.twegrsub <- knz.twegr[,which(names(knz.twegr) %in% traits.comp$sp_name)]
knz.twegrsub$comp.tot <- rowSums(knz.twegrsub)
min(knz.twegrsub$comp.tot) # min is 13.94
max(knz.twegrsub$comp.tot) # max is 33.61

################# NO COMMUNITIES MEET THRESHOLD SET FOR ANALYSIS ############################



