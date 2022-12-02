########################################
#### TKNZ RAIT DATA CLEANING SCRIPT ####
#######################################

## CREATED BY: KAITLIN KIMMEL 

## LIBRARIES
library(here)
library(dplyr)
library(stringr)
library(readxl)


### KONZA TRAIT DATA

kt1 <- read.csv(here("data/Konza/kon_firegraze_trait.csv"))
sp_names <- read.csv(here("data/Konza/SANA_grass_spp.csv"))
knz_sp <- read_xlsx(here("data/Konza/sp_list.xlsx"))
# trait data from watersheds that are burned annually with no grazers
#kt2 <- kt1[which(kt1$fire == 1 & kt1$trt == 'ungrazed'),]
# calculating based on all the watersheds to see if we have more sp. coverage

#TRAITS: 
# 1. plant height
# 2. leaf area
# 3. specific leaf area
# 4. leaf dry matter content 
# 5. stomatal length 
# 6. stomatal density
# 7. stomatal pore area index 
# 8. leaf %N 
# 9. leaf %C 
# 10. d13C
# 11. photosynthetic pathway
# 12. growth form

kt2 <- kt1[,c(1,5:13,14:17)]

trait_means <- kt2 %>%
  group_by(sp_code) %>%
  summarize(height = mean(height, na.rm = TRUE), leaf_area = mean(leaf_area, na.rm = TRUE),
            SLA = mean(SLA, na.rm = TRUE), LDMC = mean(LDMC, na.rm = TRUE), st_len = mean(stom_length, na.rm = TRUE),
            st_den = mean(stom_den, na.rm = TRUE), spi = mean(SPI, na.rm = TRUE),
            leaf_N = mean(N, na.rm = TRUE), leaf_C = mean(C, na.rm = TRUE), d13C = mean(C13, na.rm = TRUE))

cat_traits<- unique(kt2[,c(1,13,14)])
cat_traits <- cat_traits[-which(cat_traits$C3_C4 == ""),]

trait_means <- merge(trait_means, cat_traits, all.x = TRUE)

# ## Some species with missing traits. see if can fill in from rest of dataset
# 
# kt3 <- kt1[-which(kt1$fire == 1 & kt1$trt == 'ungrazed'),] # pull other fire & grazing treatments
# 
# # find mean values
# trait_means2 <- kt3 %>% 
#   group_by(sp_code) %>%
#   summarize(height = mean(height, na.rm = TRUE), leaf_area = mean(leaf_area, na.rm = TRUE),
#             SLA = mean(SLA, na.rm = TRUE), LDMC = mean(LDMC, na.rm = TRUE), st_len = mean(stom_length, na.rm = TRUE),
#             st_den = mean(stom_den, na.rm = TRUE), spi = mean(SPI, na.rm = TRUE),
#             leaf_N = mean(N, na.rm = TRUE), leaf_C = mean(C, na.rm = TRUE), d13C = mean(C13, na.rm = TRUE))
# 
# 
# # BOUCUR is the only one we can fill in completely. 
# trait_means$height[which(trait_means$sp_code == "BOUCUR")] <- trait_means2$height[which(trait_means$sp_code == "BOUCUR")]
# 
# # BROINE can fill in leaf N, leaf C, and d13C
# trait_means$leaf_N[which(trait_means$sp_code == "BROINE")] <- trait_means2$leaf_N[which(trait_means$sp_code == "BROINE")]
# trait_means$leaf_C[which(trait_means$sp_code == "BROINE")] <- trait_means2$leaf_C[which(trait_means$sp_code == "BROINE")]
# trait_means$d13C[which(trait_means$sp_code == "BROINE")] <- trait_means2$d13C[which(trait_means$sp_code == "BROINE")]
# 
#replace NaN with NA
trait_means$LDMC[which(is.nan(trait_means$LDMC))] <- NA
trait_means$st_len[which(is.nan(trait_means$st_len))] <- NA
trait_means$st_den[which(is.nan(trait_means$st_den))] <- NA
trait_means$spi[which(is.nan(trait_means$spi))] <- NA

# replace sp codes with sp names and numbers

trait_means<- merge(trait_means, sp_names, all.x = TRUE)

trait_means[which(is.na(trait_means$sp_name)),]
# Some speices not in the list - match from other konza sp data
# MUHRAC = muhlenbergia racemosa
# MUHPAN = ???
# PANOLI = dichanthelium oligosanthe
# SETVIR = setaria viridis
# SPOCRY = sporobolus cryptandrus

trait_means$sp_name[which(trait_means$sp_code == "MUHRAC")] <- "muhlenbergia_racemosa"
trait_means$sp_name[which(trait_means$sp_code == "PANOLI")] <- "dichanthelium_oligosanthe"
trait_means$sp_name[which(trait_means$sp_code == "SETVIR")] <- "setaria_viridis"
trait_means$sp_name[which(trait_means$sp_code == "SPOCRY")] <- "sporobolus_cryptandrus"
trait_means$sp_name[which(trait_means$sp_code == "MUHPAN")] <- "MUHPAN"

trait_means <- trait_means[,c(14,2:13)]

# save trait data

write.csv(trait_means, here("data/Cleaned/knz_traitsout.csv"), row.names = FALSE)

## ONLY 11 species with full species data

