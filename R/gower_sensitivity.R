#####################################
#### SENSITIVITY ANALYSES GOWER ####
####################################

# CREATED BY: KAITLIN KIMMEL

# load libraries
library(here)
library(tidyr)


## LOAD DATA

sev.black <- read.csv(here('data/Cleaned/sevblack_sen.csv'))
sev.blue <- read.csv(here('data/Cleaned/sevblue_sen.csv'))
cdr.1 <- read.csv(here('data/Cleaned/cdr1_sen.csv'))


### long to wide format for analysis
sev.black <- sev.black %>% pivot_wider(names_from = m, values_from = FRich)
sev.blue <- sev.blue %>% pivot_wider(names_from = m, values_from = FRich)
cdr.1 <- sev.blue %>% pivot_wider(names_from = m, values_from = FRich)

# Correlations between FRich values with different m values
cor(sev.black$`2`, sev.black$`3`)

### values are exactly the same using different m values. Not going to go further. 