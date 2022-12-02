########################################
#### COMMUNITY DATA CLEANING SCRIPT ####
#######################################

## CREATED BY: KAITLIN KIMMEL 

## LOAD LIBRARIES
library(here)
library(dplyr)
library(tidyr)

# BIOCON COMMUNITY DATA
cdr.biocon <- read.delim("https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/301/9/c41fc5de0beadc305668d5d2a4d40b7b", header = TRUE)

max(cdr.biocon$Year) # 2018 

# get data only from the 16 species plots from 2018 in the August sample (end of growing season)
cdr.biocon <- cdr.biocon[cdr.biocon$Year == 2018 & cdr.biocon$CountOfSpecies == 16 & 
                           cdr.biocon$Season == 'August',]
# get rid of unneccessary columns
cdr.biocon <- cdr.biocon[,c(4,5,6,7,15,16)]
# look at what categories are recorded & remove any categories that are not species

unique(cdr.biocon$Species)

cdr.biocon <- cdr.biocon[-which(cdr.biocon$Species %in% c("Bare ground", "Miscellaneous litter", "Real weeds", "Moss")),]


temp <- cdr.biocon %>% 
  pivot_wider(names_from = "Species", values_from = "Percent.cover")

temp$tot.cov <- rowSums(temp[,c(5:20)])
temp$multiplyer <- 100/temp$tot.cov

efor (i in 1:nrow(temp)){
  temp[i,c(5:20)] <- temp[i,c(5:20)]*temp$multiplyer[i]
}

temp$check <- rowSums(temp[,c(5:20)])
colSums(temp[,c(5:20)]) 
# note some species never recorded in 16 species plots in 2018
# Achillea, Anmone, & Koeleria

## For FD package: rows are 'sites' in this case our plots & columns are 'species'
cdr.biocon <- temp[,c(1:20)]

## pull out different treatment groups and save .csv files

cdrnaca <- cdr.biocon[cdr.biocon$Nitrogen.Treatment == "Namb" & cdr.biocon$CO2.Treatment == "Camb",]
cdrneca <- cdr.biocon[cdr.biocon$Nitrogen.Treatment == "Nenrich" & cdr.biocon$CO2.Treatment == "Camb",]
cdrnace <- cdr.biocon[cdr.biocon$Nitrogen.Treatment == "Namb" & cdr.biocon$CO2.Treatment == "Cenrich",]
cdrnece <- cdr.biocon[cdr.biocon$Nitrogen.Treatment == "Nenrich" & cdr.biocon$CO2.Treatment == "Cenrich",]

write.csv(cdrnaca, here("data/Cleaned/cdrnaca_commout.csv"), row.names = FALSE)
write.csv(cdrneca, here("data/Cleaned/cdrneca_commout.csv"), row.names = FALSE)
write.csv(cdrnace, here("data/Cleaned/cdrnace_commout.csv"), row.names = FALSE)
write.csv(cdrnece, here("data/Cleaned/cdrnece_commout.csv"), row.names = FALSE)

