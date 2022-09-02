#####################################
#### TRAIT DATA CLEANING SCRIPT ####
####################################

## CREATED BY: KAITLIN KIMMEL 

## LIBRARIES
library(here)
library(dplyr)

## CEDAR CREEK TRAIT DATA
cdr.SLA <- read.delim("https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/322/8/5400b565d3fc39248e91632603db36ef", header = TRUE)
cdr.istar <- read.delim("https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/342/9/efdc9389203f05db5ddc3978471e8dc4", header = TRUE)
cdr.rstar <- read.delim("https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/324/8/854eb92af9a89cc75ae67bb02720ace3", header = TRUE)
cdr.rootcn <- read.delim("https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/298/10/9fa0cda40a2ed18481d0f1b8ba6d1b36", header = TRUE)
cdr.rootbio <- read.delim("https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/325/10/cf6dadc4d72a6a37711c8e9637fe4890", header = TRUE)
cdr.shootcn <- read.delim("https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/305/10/c95883ce27d8d3ee916d6c574a5a7d8b", header=TRUE)
cdr.seed <- read.delim("https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/329/8/a375253093fc75c213f5a6ca2c4c20df", header = TRUE)


### Averages by species over all observations in monoculture

## SLA
unique(cdr.SLA$Species) # only 11 species have SLA values

cdr.SLA <- cdr.SLA %>% 
  group_by(monospecies, CO2.Treatment, Nitrogen.Treatment) %>% 
  summarize(SLA = mean(SLA..cm.2.g.))

## I*
cdr.istar <- cdr.istar[cdr.istar$CountOfSpecies ==1 & cdr.istar$Clip.Strip. == "N",]
cdr.istar <- cdr.istar[-which(is.na(cdr.istar$proportionTrans)),]
unique(cdr.istar$monospecies)
cdr.istar <- cdr.istar %>%
  group_by(monospecies, CO2.Treatment, Nitrogen.Treatment) %>% 
  summarize(I_star = mean(proportionTrans))

## R*

cdr.rstar <- cdr.rstar[cdr.rstar$CountOfSpecies == 1 & cdr.rstar$Depth == "0-20",]
cdr.rstar <- cdr.rstar %>%
  group_by(monospecies, CO2.Treatment, Nitrogen.Treatment) %>% 
  summarize(R_star = mean(NO2.NO3.mg.kg.))

### Root %C & %N

cdr.rootcn <- cdr.rootcn[cdr.rootcn$CountOfSpecies == 1 & cdr.rootcn$Root.Sampling.Depth == "0-20cm",]
cdr.rootcn <- cdr.rootcn[-which(is.na(cdr.rootcn$Percent.Carbon)),]
cdr.rootcn <- cdr.rootcn %>%
  group_by(monospecies, CO2.Treatment, Nitrogen.Treatment) %>% 
  summarize(Root.C = mean(Percent.Carbon), Root.N = mean(Percent.Nitrogen))

