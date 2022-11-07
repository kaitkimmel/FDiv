#####################################
#### TRAIT DATA CLEANING SCRIPT ####
####################################

## CREATED BY: KAITLIN KIMMEL 

## LIBRARIES
library(here)
library(dplyr)
library(stringr)

## CEDAR CREEK TRAIT DATA
cdr.SLA <- read.delim("https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/322/8/5400b565d3fc39248e91632603db36ef", header = TRUE)
cdr.SLA2 <- read.delim("https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/538/9/d6a22544411064d397041a9548f73074", header = TRUE)
cdr.istar <- read.delim("https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/342/9/efdc9389203f05db5ddc3978471e8dc4", header = TRUE)
cdr.rstar <- read.delim("https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/324/8/854eb92af9a89cc75ae67bb02720ace3", header = TRUE)
cdr.rootcn <- read.delim("https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/298/10/9fa0cda40a2ed18481d0f1b8ba6d1b36", header = TRUE)
cdr.rootbio <- read.delim("https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/325/10/cf6dadc4d72a6a37711c8e9637fe4890", header = TRUE)
cdr.shootcn <- read.delim("https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/305/10/c95883ce27d8d3ee916d6c574a5a7d8b", header=TRUE)
cdr.seed <- read.delim("https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/329/8/a375253093fc75c213f5a6ca2c4c20df", header = TRUE)
cdr.bigbiotraits <- read.delim("https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/414/8/f09cbe352ad7b89c4bd60f5a772dee10", header = TRUE)

### Averages by species over all observations in monoculture

## SLA
# unique(cdr.SLA$Species) # only 11 species have SLA values
# 
# cdr.SLA <- cdr.SLA %>% 
#   group_by(monospecies, CO2.Treatment, Nitrogen.Treatment) %>% 
#   summarize(SLA = mean(SLA..cm.2.g.))

cdr.SLA2$SLA..cm2.g. <- as.numeric(cdr.SLA2$SLA..cm2.g.)
cdr.SLA2 <- cdr.SLA2[cdr.SLA2$CountOfSpecies ==1,]
cdr.SLA2 <- cdr.SLA2[-which(is.na(cdr.SLA2$SLA..cm2.g.)),]
cdr.SLA2 <- cdr.SLA2 %>%
  mutate(measured.species = replace(measured.species, measured.species == "Koleria cristata","Koeleria cristata")) %>%
  group_by(measured.species, CO2.Treatment, Nitrogen.Treatment) %>%
  summarize(SLA = mean(SLA..cm2.g.))
names(cdr.SLA2)[1] <- 'monospecies'


## I*
cdr.istar <- cdr.istar[cdr.istar$CountOfSpecies ==1 & cdr.istar$Clip.Strip. == "N",]
cdr.istar <- cdr.istar[-which(is.na(cdr.istar$proportionTrans)),]
unique(cdr.istar$monospecies)
cdr.istar <- cdr.istar %>%
  group_by(monospecies, CO2.Treatment, Nitrogen.Treatment) %>% 
  summarize(I_star = mean(proportionTrans))

## R*

cdr.rstar <- cdr.rstar[cdr.rstar$CountOfSpecies == 1 & cdr.rstar$Depth == "0-20",]
cdr.rstar <- cdr.rstar[-which(is.na(cdr.rstar$NO2.NO3.mg.kg.)),]
cdr.rstar <- cdr.rstar %>%
  group_by(monospecies, CO2.Treatment, Nitrogen.Treatment) %>% 
  summarize(R_star = mean(NO2.NO3.mg.kg.))

### Root %C & %N

cdr.rootcn <- cdr.rootcn[cdr.rootcn$CountOfSpecies == 1 & cdr.rootcn$Root.Sampling.Depth == "0-20cm",]
cdr.rootcn <- cdr.rootcn[-which(is.na(cdr.rootcn$Percent.Carbon)),]
cdr.rootcn <- cdr.rootcn %>%
  group_by(monospecies, CO2.Treatment, Nitrogen.Treatment) %>% 
  summarize(Root_C = mean(Percent.Carbon), Root_N = mean(Percent.Nitrogen))

## Root Biomass
cdr.rootbio <- cdr.rootbio[cdr.rootbio$CountOfSpecies == 1 & cdr.rootbio$Depth == "0-20" & cdr.rootbio$Sub.Sample != "Crowns",]
unique(cdr.rootbio$monospecies)# Need to clean names 
cdr.rootbio <- cdr.rootbio[-which(is.na(cdr.rootbio$Mass..g.m2.)),]

cdr.rootbio <- cdr.rootbio %>%
  mutate( monospecies = replace(monospecies, monospecies == "AchilleaMillefolium","Achillea millefolium")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "AgropyronRepens","Agropyron repens")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "AmorphaCanescens","Amorpha canescens")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "AndropogonGerardi","Andropogon gerardi")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "AnemoneCylindrica","Anemone cylindrica")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "AsclepiasTuberosa","Asclepias tuberosa")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "BoutelouaGracilis","Bouteloua gracilis")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "BromusInermis","Bromus inermis")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "KoeleriaCristata","Koeleria cristata")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "LespedezaCapitata","Lespedeza capitata")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "LupinusPerennis","Lupinus perennis")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "PetalostemumVillosum","Petalostemum villosum")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "PoaPratensis","Poa pratensis")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "SchizachyriumScoparium","Schizachyrium scoparium")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "SolidagoRigida","Solidago rigida")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "SorghastrumNutans","Sorghastrum nutans")) %>%
  group_by(monospecies, CO2.Treatment, Nitrogen.Treatment) %>% 
  summarize(Root_biomass = mean(Mass..g.m2.))

## Shoot %C and %N

cdr.shootcn <- cdr.shootcn[cdr.shootcn$CountOfSpecies == 1,]
cdr.shootcn <- cdr.shootcn[-which(is.na(cdr.shootcn$Carbon...)),]
unique(cdr.shootcn$monospecies)
cdr.shootcn <- cdr.shootcn %>%
  mutate( monospecies = replace(monospecies, monospecies == "AchilleaMillefolium","Achillea millefolium")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "AgropyronRepens","Agropyron repens")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "AmorphaCanescens","Amorpha canescens")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "AndropogonGerardi","Andropogon gerardi")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "AnemoneCylindrica","Anemone cylindrica")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "AsclepiasTuberosa","Asclepias tuberosa")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "BoutelouaGracilis","Bouteloua gracilis")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "BromusInermis","Bromus inermis")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "KoeleriaCristata","Koeleria cristata")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "LespedezaCapitata","Lespedeza capitata")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "LupinusPerennis","Lupinus perennis")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "PetalostemumVillosum","Petalostemum villosum")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "PoaPratensis","Poa pratensis")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "SchizachyriumScoparium","Schizachyrium scoparium")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "SolidagoRigida","Solidago rigida")) %>%
  mutate(monospecies = replace(monospecies, monospecies == "SorghastrumNutans","Sorghastrum nutans")) %>%
  group_by(monospecies, CO2.Treatment, Nitrogen.Treatment) %>% 
  summarize(Shoot_C = mean(Carbon...), Shoot_N = mean(Nitrogen...))

# seed mass
cdr.seed <- cdr.seed[cdr.seed$CountOfSpecies == 1,]
unique(cdr.seed$monospecies) # 11 species total, missing Achillea, Amorpha, Anemone, Asclepias, petalostemum

# Grasses = avg seed weight, other species = seed weight need to make one column

cdr.seed$seed_mass <- cdr.seed$Seed.wt..g.

for (i in 1:nrow(cdr.seed)){
  if(cdr.seed$monospecies[i] %in%
     c('Koeleria cristata', 'Bouteloua gracilis',	'Agropyron repens', 
       'Sorghastrum nutans',  'Andropogon gerardi', 'Schizachyrium scoparium',
       'Solidago rigida')){
    cdr.seed$seed_mass[i] <- cdr.seed$Avg.wt.per.seed..g.[i]
  }
}


cdr.seed <- cdr.seed %>% 
  group_by(monospecies, CO2.Treatment, Nitrogen.Treatment) %>%
  summarize(Seed_mass = mean(seed_mass))
### NEED TO LOOK INTO DATASET METHODS>>> DIFFERENCE BETWEEN AVG WEIGHT AND WEIGHT

##### NOTES:
## 3 trait datasets that do not span all species in all treatment combinations
## SLA (56/64) - missing 2 species (Poa, Bouteloua)
## Root %C & %N (61/64) - Anemone cylindrica only in Cenrich, Namb
#### Going to pull this number for the other treatment combinations
## Seed mass (44/64) - missing 5 species (Achillea, Amorpha, Anemone, Asclepias, Petalostemum)
## Pulled trait data from Big Bio experiment (neighboring BioCON) to grab traits
### Only measured in ambient conditions, but will apply to all treatments for this study
### We are more concerned with trait combinations than precise values at this point
### Think it is better to have a full data set - will update methods since 
### Changed a bit from original plan in Registered Report Stage 1



# Create single dataset of all cdr trait data

#put all data frames into list
newdf <- merge(cdr.SLA2, cdr.istar, all = TRUE)
newdf <- merge(newdf, cdr.rstar, all = TRUE)
newdf <- merge(newdf, cdr.rootcn, all = TRUE)
newdf <- merge(newdf, cdr.rootbio, all = TRUE)
newdf <- merge(newdf, cdr.shootcn, all = TRUE)
newdf <- merge(newdf, cdr.seed, all = TRUE)

# Copying Anemone data from Root %C & %N to other treatments
newdf$Root_C[newdf$monospecies == "Anemone cylindrica"] <- newdf$Root_C[newdf$monospecies == "Anemone cylindrica" & newdf$CO2.Treatment == "Cenrich" & newdf$Nitrogen.Treatment == "Namb"]
newdf$Root_N[newdf$monospecies == "Anemone cylindrica"] <- newdf$Root_N[newdf$monospecies == "Anemone cylindrica" & newdf$CO2.Treatment == "Cenrich" & newdf$Nitrogen.Treatment == "Namb"]

# Adding in Big Bio traits to fill in missing ones
newdf$SLA[newdf$monospecies == "Bouteloua gracilis"] <- cdr.bigbiotraits$SLA_cm2.g[cdr.bigbiotraits$Species == "Bouteloua gracilis"]
newdf$SLA[newdf$monospecies == "Poa pratensis"] <- cdr.bigbiotraits$SLA_cm2.g[cdr.bigbiotraits$Species == "Poa pratensis"]
newdf$Seed_mass[newdf$monospecies == "Achillea millefolium"] <- cdr.bigbiotraits$Seed_weight_.g.[cdr.bigbiotraits$Species == "Achillea millefolium"]
newdf$Seed_mass[newdf$monospecies == "Amorpha canescens"] <- cdr.bigbiotraits$Seed_weight_.g.[cdr.bigbiotraits$Species == "Amorpha canescens"]
newdf$Seed_mass[newdf$monospecies == "Anemone cylindrica"] <- cdr.bigbiotraits$Seed_weight_.g.[cdr.bigbiotraits$Species == "Anemone cylindrica"]
newdf$Seed_mass[newdf$monospecies == "Asclepias tuberosa"] <- cdr.bigbiotraits$Seed_weight_.g.[cdr.bigbiotraits$Species == "Asclepias tuberosa"]
newdf$Seed_mass[newdf$monospecies == "Petalostemum villosum"] <- cdr.bigbiotraits$Seed_weight_.g.[cdr.bigbiotraits$Species == "Petalostemum villosa"] # nomenclature differences b/t experiments

newdf$monospecies <- newdf$monospecies %>% str_replace_all(" ", ".")

### Split traits into treatments for future analysis & save .csv files
cdr.naca <- newdf[which(newdf$CO2.Treatment == "Camb" & newdf$Nitrogen.Treatment == "Namb"),]
cdr.nace <- newdf[which(newdf$CO2.Treatment == "Cenrich" & newdf$Nitrogen.Treatment == "Namb"),]
cdr.neca <- newdf[which(newdf$CO2.Treatment == "Camb" & newdf$Nitrogen.Treatment == "Nenrich"),]
cdr.nece <- newdf[which(newdf$CO2.Treatment == "Cenrich" & newdf$Nitrogen.Treatment == "Nenrich"),]

write.csv(cdr.naca, here("data/Cleaned/cdrnaca_traitsout.csv"), row.names = FALSE)
write.csv(cdr.nace, here("data/Cleaned/cdrnace_traitsout.csv"), row.names = FALSE)
write.csv(cdr.neca, here("data/Cleaned/cdrneca_traitsout.csv"), row.names = FALSE)
write.csv(cdr.nece, here("data/Cleaned/cdrnece_traitsout.csv"), row.names = FALSE)