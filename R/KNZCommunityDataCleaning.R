############################################
#### KNZ COMMUNITY DATA CLEANING SCRIPT ####
############################################

## CREATED BY: KAITLIN KIMMEL 

## LOAD LIBRARIES
library(here)
library(dplyr)
library(tidyr)
library(readxl)


# Load Data
# community composition data
dat <- read.csv("https://pasta.lternet.edu/package/data/eml/knb-lter-knz/69/19/63768b48f41e790a40e7fa4f9267c3a2")
# species code to species name data
sp_dat <- read_xlsx(here("data/Konza/sp_list.xlsx"))

# using 2010 - the year the data was collected
dat <- dat[which(dat$RecYear == 2010),]

# selecting 4 watersheds for analysis
# watersheds: 1d (annual burn ungrazed) 
# 20b (20 yr burn ungrazed) 
# N1b (annual burn bison) 
# N20b (20 yr burn bison)
watersheds <- c("001d", "020b", "n01b", "n20b")
dat <- dat[which(dat$WaterShed %in% watersheds),]

# subset out columns needed
dat <- dat[,c(4,6,8,9,10,13)]
# get unique plot id for each watershed, transect, plot combo

uid <- unique(dat[,c(2:4)])
uid$uid <- 1:nrow(uid)
dat <- merge(dat, uid)

dat <- dat %>% group_by(WaterShed, uid, SpeCode) %>%
  summarize(cover = max(Cover))

# convert cover classes to midpoints
# cover classes & midpoints
# 1 - 0-1% cover 0.5%
# 2 - 2-5% cover 3.5%
# 3 - 5-25% cover 15%
# 4 - 25-50% cover 37.5%
# 5 - 50-75% cover 62.5%
# 6 - 75-95% cover 85%
# 7 - 95-100% cover 97.5%

dat$cover[which(dat$cover == 1)] <- 0.5
dat$cover[which(dat$cover == 2)] <- 3.5
dat$cover[which(dat$cover == 3)] <- 15
dat$cover[which(dat$cover == 4)] <- 37.5
dat$cover[which(dat$cover == 5)] <- 62.5
dat$cover[which(dat$cover == 6)] <- 85
dat$cover[which(dat$cover == 7)] <- 97.5

# add species names to dataset
sp_dat$sp_name <- paste(sp_dat$genus, sp_dat$species, sep = "_")
names(sp_dat)[1] <- "SpeCode"

dat <- merge(dat, sp_dat[,c(1,11)])


temp <- dat[,c(2:5)] %>% 
  pivot_wider(names_from = "sp_name", values_from = "cover", values_fill = 0)

temp$tot.cov <- rowSums(temp[,c(3:177)])
temp$multiplyer <- 100/temp$tot.cov
for (i in 1:nrow(temp)){
  temp[i,c(3:177)] <- temp[i,c(3:177)]*temp$multiplyer[i]
}

temp$check <- rowSums(temp[,c(3:177)])

fin.dat <- temp[,c(1:177)]

ann_ungraze <- fin.dat[which(fin.dat$WaterShed == "001d"),]
ann_graze <- fin.dat[which(fin.dat$WaterShed == "020b"),]
twenty_ungraze <- fin.dat[which(fin.dat$WaterShed == "n01b"),]
twenty_graze <- fin.dat[which(fin.dat$WaterShed == "n20b"),]

write.csv(ann_ungraze, here("data/Cleaned/knzannun_commout.csv"), row.names = FALSE)
write.csv(ann_graze, here("data/Cleaned/knzanngr_commout.csv"), row.names = FALSE)
write.csv(twenty_ungraze, here("data/Cleaned/knztweun_commout.csv"), row.names = FALSE)
write.csv(twenty_graze, here("data/Cleaned/knztwegr_commout.csv"), row.names = FALSE)

