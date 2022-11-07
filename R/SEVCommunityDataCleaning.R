############################################
#### SEV COMMUNITY DATA CLEANING SCRIPT ####
############################################

## CREATED BY: KAITLIN KIMMEL 

## LOAD LIBRARIES
library(here)
library(dplyr)
library(tidyr)


# load data
sev.comm <- read.csv(here('data/sevcommunity_KOtrait.csv'), row.names = 1)

## For FD package: rows are 'sites' in this case our plots & columns are 'species'

temp <- sev.comm[,c(1,3,4,7)] %>% 
  pivot_wider(names_from = "kartez", values_from = "cover", values_fill = 0)

temp$tot.cov <- rowSums(temp[,c(3:50)])
temp$multiplyer <- 100/temp$tot.cov
for (i in 1:nrow(temp)){
  temp[i,c(3:50)] <- temp[i,c(3:50)]*temp$multiplyer[i]
}

temp$check <- rowSums(temp[,c(3:50)]) # all = 100

#split communities
sev.blue <- temp[which(temp$site== 'core_blue'),]
sev.black <- temp[which(temp$site== 'core_black'),]

#save community data
write.csv(sev.blue,(here('data/Cleaned/sevblue_commout.csv')), row.names = FALSE)
write.csv(sev.black,(here('data/Cleaned/sevblack_commout.csv')), row.names = FALSE)

