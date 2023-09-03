#########################################
#### SEV TRAIT DATA CLEANING SCRIPT ####
########################################

## CREATED BY: KAITLIN KIMMEL 

## LIBRARIES
library(here)
library(dplyr)
library(stringr)

## SEV TRAIT DATA

sev.trait <- read.csv("https://pasta.lternet.edu/package/data/eml/knb-lter-sev/333/1/d6b91bac4e7d8e92e71c70f5d34d2af4", header = TRUE)
sev.height <- read.csv(here("data/Sev/sev_npp_height_biomass.csv")) # ASK TIM FOR PASTA LINK

# # TRAITS 
# 1. maximum plant height
# 2. leaf dry matter content 
# 3. specific leaf area 
# 4. d15N 
# 5. d13C
# 6. leaf %N 
# 7. leaf %C
# 8. stem dry matter content 
# 9. root dry matter content 
# 10. photosynthetic pathway

sev.trait <- sev.trait[,c(1:5,19,18,22,23,24,25,37,41,9)] # get necessary columns

# find mean value for each trait for each species
sev.means <- sev.trait %>% 
  group_by(Kartez, genus_new, species_new) %>% 
  summarize(ldmc = mean(ldmc, na.rm = TRUE), sla = mean(sla, na.rm = TRUE), 
            d15N = mean(lf_iso_d15N, na.rm = TRUE), d13C = mean(lf_iso_d13C, na.rm = TRUE),
            perc_N = mean(lf_iso_pN, na.rm = TRUE), perc_C = mean(lf_iso_pC, na.rm = TRUE),
            stem_dmc = mean(sdmc, na.rm = TRUE), root_dmc = mean(rdmc, na.rm = TRUE))

# pull out photosynthetic pathway
sev.photo <- unique(sev.trait[,c(2,14)])
# merge photosynthetic pathway to means
sev.means <- merge(sev.means, sev.photo)

#find mean height of species
sev.height <- sev.height %>%
  group_by(kartez) %>%
  summarize(height = mean(height))

names(sev.height)[1] = 'Kartez'

# add to mean trait data
sev.means <- merge(sev.means, sev.height, all.x = TRUE)

# change NaN to NA
sev.means$ldmc[which(is.nan(sev.means$ldmc))] <- NA
sev.means$sla[which(is.nan(sev.means$sla))] <- NA
sev.means$d15N[which(is.nan(sev.means$d15N))] <- NA
sev.means$d13C[which(is.nan(sev.means$d13C))] <- NA
sev.means$perc_N[which(is.nan(sev.means$perc_N))] <- NA
sev.means$perc_C[which(is.nan(sev.means$perc_C))] <- NA
sev.means$stem_dmc[which(is.nan(sev.means$stem_dmc))] <- NA
sev.means$root_dmc[which(is.nan(sev.means$root_dmc))] <- NA

# save data
write.csv(sev.means, here("data/Cleaned/sev_traitsout.csv"), row.names = FALSE)


