########################################
#### TKNZ RAIT DATA CLEANING SCRIPT ####
#######################################

## CREATED BY: KAITLIN KIMMEL 

## LIBRARIES
library(here)
library(dplyr)
library(stringr)


### KONZA TRAIT DATA

kt1 <- read.csv(here("data/Konza/kon_firegraze_trait.csv"))

# trait data from watersheds that are burned annually with no grazers

kt2 <- kt1[which(kt1$fire == 1 & kt1$trt == 'ungrazed'),]

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

kt2 <- kt2[,c(1,5:13,14:17)]

trait_means <- kt2 %>%
  group_by(sp_code) %>%
  summarize(height = mean(height, na.rm = TRUE), leaf_area = mean(leaf_area, na.rm = TRUE),
            SLA = mean(SLA, na.rm = TRUE), LDMC = mean(LDMC, na.rm = TRUE), st_len = mean(stom_length, na.rm = TRUE),
            st_den = mean(stom_den, na.rm = TRUE), spi = mean(SPI, na.rm = TRUE),
            leaf_N = mean(N, na.rm = TRUE), leaf_C = mean(C, na.rm = TRUE), d13C = mean(C13, na.rm = TRUE))

cat_traits<- unique(kt2[,c(1,13,14)])
cat_traits <- cat_traits[-which(cat_traits$C3_C4 == ""),]

trait_means <- merge(trait_means, cat_traits, all.x = TRUE)

## Some species with missing traits. see if can fill in from rest of dataset

kt3 <- kt1[-which(kt1$fire == 1 & kt1$trt == 'ungrazed'),]

trait_means2 <- kt3 %>% 
  group_by(sp_code) %>%
  summarize(height = mean(height, na.rm = TRUE), leaf_area = mean(leaf_area, na.rm = TRUE),
            SLA = mean(SLA, na.rm = TRUE), LDMC = mean(LDMC, na.rm = TRUE), st_len = mean(stom_length, na.rm = TRUE),
            st_den = mean(stom_den, na.rm = TRUE), spi = mean(SPI, na.rm = TRUE),
            leaf_N = mean(N, na.rm = TRUE), leaf_C = mean(C, na.rm = TRUE), d13C = mean(C13, na.rm = TRUE))
