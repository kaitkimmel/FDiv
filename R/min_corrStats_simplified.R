################################################################
#### STATISTICAL ANALYSES  EUCLIDEAN for TRAIT CORRELATIONS ####
################################################################

# CREATED BY: KAITLIN KIMMEL

# load libraries
library(nlme)
library(here)
library(AICcmodavg)
library(ggplot2)
library(ggpubr)
library(dplyr)


# load data

cdr.1 <- read.csv(here("data/Cleaned/cdr1_sc.csv"))
cdr.1 <- cdr.1[-which(is.na(cdr.1$min_cor)),]
cdr.2 <- read.csv(here("data/Cleaned/cdr2_sc.csv"))
cdr.2 <- cdr.2[-which(is.na(cdr.2$min_cor)),]
cdr.3 <- read.csv(here("data/Cleaned/cdr3_sc.csv"))
cdr.3 <- cdr.3[-which(is.na(cdr.3$min_cor)),]
cdr.4 <- read.csv(here("data/Cleaned/cdr4_sc.csv"))
cdr.4 <- cdr.4[-which(is.na(cdr.4$min_cor)),]
sev.blue <- read.csv(here("data/Cleaned/sevblue_sc.csv"))
sev.blue <- sev.blue[-which(is.na(sev.blue$min_cor)),]
sev.black <- read.csv(here("data/Cleaned/sevblack_sc.csv"))
sev.black <- sev.black[-which(is.na(sev.black$min_cor)),]



#############################################################
################ SEV #######################################
###########################################################

###########################################################
################### BLUE ################################
##########################################################

###### Fit Metric ~ correlation models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ min_cor, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FRic ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots, sev.blue[-which(is.na(sev.blue$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(FRic~1, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) ### mod1
fricmod_blue <- lme(FRic ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE richness
mod <- lme(kde.alpha ~ min_cor, random = ~1|sev.blueplots,  data = sev.blue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.alpha ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots, sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(kde.alpha~1, random = ~1|sev.blueplots,  data = sev.blue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) ### mod1
kde.alphamod_blue <- lme(kde.alpha ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))


### FEve
mod <- lme(FEve ~ min_cor, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FEve ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(FEve~1, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) ## mod & 4 similar
fevemod_blue <- lme(FEve ~ 1, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE evenness
mod <- lme(kde.evenness ~ min_cor, random = ~1|sev.blueplots,  data = sev.blue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.evenness ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots, sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(kde.evenness~1, random = ~1|sev.blueplots,  data = sev.blue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) ### mod1
kde.evennessmod_blue <- lme(kde.evenness ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

## FDis
mod <- lme(FDis ~ min_cor, random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FDis ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(FDis~1, random = ~1|sev.blueplots,  data = sev.blue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) #mod1
fdismod_blue <- lme(FDis ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE dispersion
mod <- lme(kde.dispersion ~ min_cor, random = ~1|sev.blueplots,  data = sev.blue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.dispersion ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots, sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(kde.dispersion~1, random = ~1|sev.blueplots,  data = sev.blue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) ### mod1
kde.dispersionmod_blue <- lme(kde.dispersion ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

# FDiv
mod <- lme(FDiv ~ min_cor, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FDiv ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(FDiv~1, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) #mod4
fdivmod_blue <- lme(FDiv ~ 1, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

### Rao's Q
mod <- lme(RaoQ ~ min_cor, random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(RaoQ ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots),control = lmeControl(opt = "optim"))
mod4 <- lme(RaoQ~1, random = ~1|sev.blueplots,  data = sev.blue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) #mod1
raoqmod_blue <- lme(RaoQ ~ min_cor  + I(min_cor^2), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))


summary(fricmod_blue)
summary(fevemod_blue)
summary(fdismod_blue)
summary(fdivmod_blue)
summary(raoqmod_blue)
summary(kde.alphamod_blue)
summary(kde.evennessmod_blue)
summary(kde.dispersionmod_blue)


###########################################################
################### BLACK ################################
##########################################################


###### Fit Metric ~ min_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ min_cor, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FRic ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots, sev.black[-which(is.na(sev.black$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(FRic ~ 1, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
AIC(mod, mod1, mod4) ### mod1
fricmod_black <- lme(FRic ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE richness
mod <- lme(kde.alpha ~ min_cor, random = ~1|sev.blackplots,  data = sev.black, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.alpha ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots, sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(kde.alpha ~ 1, random = ~1|sev.blackplots,  data = sev.black, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1,mod4) ### mod1
anova(mod1, mod, mod4) #no dif mod4
kde.alphamod_black <- lme(kde.alpha ~ 1, random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))

### FEve
mod <- lme(FEve ~ min_cor, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FEve ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(FEve ~ 1, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) ## mod
fevemod_black <- lme(FEve ~ min_cor, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE evenness
mod <- lme(kde.evenness ~ min_cor, random = ~1|sev.blackplots,  data = sev.black, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.evenness ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots, sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(kde.evenness ~ 1, random = ~1|sev.blackplots,  data = sev.black, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) ### mod4
kde.evennessmod_black <- lme(kde.evenness ~ 1, random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))

## FDis
mod <- lme(FDis ~ min_cor, random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FDis ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(FDis ~ 1, random = ~1|sev.blackplots,  data = sev.black, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) #mod1
fdismod_black <- lme(FDis ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE dispersion
mod <- lme(kde.dispersion ~ min_cor, random = ~1|sev.blackplots,  data = sev.black, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.dispersion ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots, sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|sev.blackplots,  data = sev.black, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) ### mod1
anova(mod1, mod, mod4) #mod4
kde.dispersionmod_black <- lme(kde.dispersion ~ 1, random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))

# FDiv
mod <- lme(FDiv ~ min_cor, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FDiv ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(FDiv ~ 1, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) #mod4
fdivmod_black <- lme(FDiv ~ 1, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

### Rao's Q
mod <- lme(RaoQ ~ min_cor, random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(RaoQ ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots),control = lmeControl(opt = "optim"))
mod4 <- lme(RaoQ ~ 1, random = ~1|sev.blackplots,  data = sev.black, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) #mod1
raoqmod_black <- lme(RaoQ ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))


summary(fricmod_black)
summary(fevemod_black)
summary(fdismod_black)
summary(fdivmod_black)
summary(raoqmod_black)
summary(kde.alphamod_black)
summary(kde.evennessmod_black)
summary(kde.dispersionmod_black)

#########################################
######### CEDAR CREEK #####################
##########################################

###########################################################
################### CDR 1 ################################
##########################################################

###### Fit Metric ~ min_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ min_cor, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FRic ~ 1, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
fricmod_cdr1 <- lme(FRic ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ min_cor, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.alpha ~ 1, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
kde.alphamod_cdr1 <- lme(kde.alpha ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ min_cor, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FEve ~ 1, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ## mod4
fevemod_cdr1 <- lme(FEve ~ 1, random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ min_cor, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.evenness ~ 1, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod
kde.evennessmod_cdr1 <- lme(kde.evenness ~ min_cor, random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ min_cor, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDis ~ 1, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod
anova(mod1, mod) #no dif mod
fdismod_cdr1 <- lme(FDis ~ min_cor, random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

##kde dspersion
mod <- lme(kde.dispersion ~ min_cor, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
anova(mod1, mod) #no dif mod
kde.dispersionmod_cdr1 <- lme(kde.dispersion ~ min_cor, random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ min_cor, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDiv ~ 1, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod
fdivmod_cdr1 <- lme(FDiv ~ min_cor, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ min_cor, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(RaoQ ~ 1, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod1
raoqmod_cdr1 <- lme(RaoQ ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod_cdr1)
summary(fevemod_cdr1)
summary(fdismod_cdr1)
summary(fdivmod_cdr1)
summary(raoqmod_cdr1)
summary(kde.alphamod_cdr1)
summary(kde.evennessmod_cdr1)
summary(kde.dispersionmod_cdr1)

###########################################################
################### CDR 2 ################################
##########################################################
###### Fit Metric ~ min_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ min_cor, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FRic ~ 1, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
fricmod_cdr2 <- lme(FRic ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ min_cor, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.alpha ~ 1, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
kde.alphamod_cdr2 <- lme(kde.alpha ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ min_cor, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FEve ~ 1, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ## mod best fit
fevemod_cdr2 <- lme(FEve ~ min_cor, random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ min_cor, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.evenness ~ 1, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod
kde.evennessmod_cdr2 <- lme(kde.evenness ~ min_cor, random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ min_cor, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDis ~ 1, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod1
fdismod_cdr2 <- lme(FDis ~ min_cor+ I(min_cor^2), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

##KDE disperion
mod <- lme(kde.evenness ~ min_cor, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.evenness ~ 1, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod
kde.dispersionmod_cdr2 <- lme(kde.evenness ~ min_cor, random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ min_cor, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDiv ~ 1, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod1
fdivmod_cdr2 <- lme(FDiv ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ min_cor, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(RaoQ ~ 1, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod1
raoqmod_cdr2 <- lme(RaoQ ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod_cdr2)
summary(fevemod_cdr2)
summary(fdismod_cdr2)
summary(fdivmod_cdr2)
summary(raoqmod_cdr2)
summary(kde.alphamod_cdr2)
summary(kde.evennessmod_cdr2)
summary(kde.dispersionmod_cdr2)

###########################################################
################### CDR 3 ################################
##########################################################
###### Fit Metric ~ min_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ min_cor, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FRic ~ 1, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
fricmod_cdr3 <- lme(FRic ~ min_cor+ I(min_cor^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ min_cor, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.alpha ~ 1, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
kde.alphamod_cdr3 <- lme(kde.alpha ~ min_cor+ I(min_cor^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ min_cor, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FEve ~ 1, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ## mod4
fevemod_cdr3 <- lme(FEve ~ 1, random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ min_cor, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.evenness ~ 1, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
kde.evennessmod_cdr3 <- lme(kde.evenness ~ min_cor+ I(min_cor^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ min_cor, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDis ~ 1, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod1
fdismod_cdr3 <- lme(FDis ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

##KDE dispersion
mod <- lme(kde.dispersion ~ min_cor, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
kde.dispersionmod_cdr3 <- lme(kde.dispersion ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ min_cor, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDiv ~ 1, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod1
fdivmod_cdr3 <- lme(FDiv ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))


### Rao's Q
mod <- lme(RaoQ ~ min_cor, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(RaoQ ~ 1, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod1
raoqmod_cdr3 <- lme(RaoQ ~ min_cor+ I(min_cor^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod_cdr3)
summary(fevemod_cdr3)
summary(fdismod_cdr3)
summary(fdivmod_cdr3)
summary(raoqmod_cdr3)
summary(kde.alphamod_cdr3)
summary(kde.evennessmod_cdr3)
summary(kde.dispersionmod_cdr3)

###########################################################
################### CDR 4 ################################
##########################################################
###### Fit Metric ~ min_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ min_cor, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FRic ~ 1, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
fricmod_cdr4 <- lme(FRic ~ min_cor+ I(min_cor^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ min_cor, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.alpha ~ 1, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
kde.alphamod_cdr4 <- lme(kde.alpha ~ min_cor+ I(min_cor^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ min_cor, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FEve ~ 1, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ## mod
fevemod_cdr4 <- lme(FEve ~ min_cor, random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ min_cor, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.evenness ~ 1, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
kde.evennessmod_cdr4 <- lme(kde.evenness ~ min_cor+ I(min_cor^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ min_cor, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDis ~ 1, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod1
fdismod_cdr4 <- lme(FDis ~ min_cor+ I(min_cor^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE dispersion
mod <- lme(kde.dispersion ~ min_cor, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod
kde.dispersionmod_cdr4 <- lme(kde.dispersion ~ min_cor, random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))


# FDiv
mod <- lme(FDiv ~ min_cor, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDiv~ 1, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod4
fdivmod_cdr4 <- lme(FDiv ~ 1, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ min_cor, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(RaoQ ~ 1, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod
raoqmod_cdr4 <- lme(RaoQ ~ min_cor, random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod_cdr4)
summary(fevemod_cdr4)
summary(fdismod_cdr4)
summary(fdivmod_cdr4)
summary(raoqmod_cdr4)
summary(kde.alphamod_cdr4)
summary(kde.evennessmod_cdr4)
summary(kde.dispersionmod_cdr4)





summary(fricmod_blue)
summary(fricmod_black)
summary(fricmod_cdr1)
summary(fricmod_cdr2)
summary(fricmod_cdr3)
summary(fricmod_cdr4)

summary(kde.alphamod_blue)
summary(kde.alphamod_black)
summary(kde.alphamod_cdr1)
summary(kde.alphamod_cdr2)
summary(kde.alphamod_cdr3)
summary(kde.alphamod_cdr4)

summary(fevemod_blue)
summary(fevemod_black)
summary(fevemod_cdr1)
summary(fevemod_cdr2)
summary(fevemod_cdr3)
summary(fevemod_cdr4)

summary(kde.evennessmod_blue)
summary(kde.evennessmod_black)
summary(kde.evennessmod_cdr1)
summary(kde.evennessmod_cdr2)
summary(kde.evennessmod_cdr3)
summary(kde.evennessmod_cdr4)

summary(fdismod_blue)
summary(fdismod_black)
summary(fdismod_cdr1)
summary(fdismod_cdr2)
summary(fdismod_cdr3)
summary(fdismod_cdr4)

summary(kde.dispersionmod_blue)
summary(kde.dispersionmod_black)
summary(kde.dispersionmod_cdr1)
summary(kde.dispersionmod_cdr2)
summary(kde.dispersionmod_cdr3)
summary(kde.dispersionmod_cdr4)

summary(fdivmod_blue)
summary(fdivmod_black)
summary(fdivmod_cdr1)
summary(fdivmod_cdr2)
summary(fdivmod_cdr3)
summary(fdivmod_cdr4)

summary(raoqmod_blue)
summary(raoqmod_black)
summary(raoqmod_cdr1)
summary(raoqmod_cdr2)
summary(raoqmod_cdr3)
summary(raoqmod_cdr4)


################### GRAPHS ##############################


# create combined data set of raw data
cdr.1$community <- "CDR1"
cdr.2$community <- "CDR2"
cdr.3$community <- "CDR3"
cdr.4$community <- "CDR4"
cdr.full <- rbind(cdr.1, cdr.2, cdr.3, cdr.4)

cdr.full <- cdr.full %>% group_by(community, min_cor) %>%
  summarize (FRic = mean(FRic), FEve = mean(FEve), FDis = mean(FDis),
             FDiv = mean(FDiv), RaoQ = mean(RaoQ), RaoQ = mean(RaoQ), kde.alpha = mean(kde.alpha), kde.evenness  = mean(kde.evenness), kde.dispersion = mean(kde.dispersion))


##############################################
############## CDR 1 ########################
#############################################
range(cdr.1$min_cor)
new.df1 <- data.frame(min_cor = seq(min(cdr.1$min_cor), max(cdr.1$min_cor), by = 0.01))
fr.1 <- as.data.frame(predictSE.lme(fricmod_cdr1, new.df1))
fr.1$min_cor <- seq(min(cdr.1$min_cor), max(cdr.1$min_cor), by = 0.01)
fr.1$lwr <- fr.1$fit - fr.1$se.fit
fr.1$upr <- fr.1$fit + fr.1$se.fit
fe.1 <- as.data.frame(predictSE.lme(fevemod_cdr1, new.df1))
fe.1$min_cor <- seq(min(cdr.1$min_cor), max(cdr.1$min_cor), by = 0.01)
fe.1$lwr <- fe.1$fit - fe.1$se.fit
fe.1$upr <- fe.1$fit + fe.1$se.fit
fdis.1 <- as.data.frame(predictSE.lme(fdismod_cdr1, new.df1))
fdis.1$min_cor <- seq(min(cdr.1$min_cor), max(cdr.1$min_cor), by = 0.01)
fdis.1$lwr <- fdis.1$fit - fdis.1$se.fit
fdis.1$upr <- fdis.1$fit + fdis.1$se.fit
fdiv.1 <- as.data.frame(predictSE.lme(fdivmod_cdr1, new.df1))
fdiv.1$min_cor <- seq(min(cdr.1$min_cor), max(cdr.1$min_cor), by = 0.01)
fdiv.1$lwr <- fdiv.1$fit - fdiv.1$se.fit
fdiv.1$upr <- fdiv.1$fit + fdiv.1$se.fit
rq.1 <- as.data.frame(predictSE.lme(raoqmod_cdr1, new.df1))
rq.1$min_cor <- seq(min(cdr.1$min_cor), max(cdr.1$min_cor), by = 0.01)
rq.1$lwr <- rq.1$fit - rq.1$se.fit
rq.1$upr <- rq.1$fit + rq.1$se.fit
kde.alpha.1 <- as.data.frame(predictSE.lme(kde.alphamod_cdr1, new.df1))
kde.alpha.1$min_cor <- seq(min(cdr.1$min_cor), max(cdr.1$min_cor), by = .01)
kde.alpha.1$lwr <- kde.alpha.1$fit - kde.alpha.1$se.fit
kde.alpha.1$upr <- kde.alpha.1$fit + kde.alpha.1$se.fit
kde.evenness.1 <- as.data.frame(predictSE.lme(kde.evennessmod_cdr1, new.df1))
kde.evenness.1$min_cor <- seq(min(cdr.1$min_cor), max(cdr.1$min_cor), by = .01)
kde.evenness.1$lwr <- kde.evenness.1$fit - kde.evenness.1$se.fit
kde.evenness.1$upr <- kde.evenness.1$fit + kde.evenness.1$se.fit
kde.dispersion.1 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdr1, new.df1))
kde.dispersion.1$min_cor <- seq(min(cdr.1$min_cor), max(cdr.1$min_cor), by = .01)
kde.dispersion.1$lwr <- kde.dispersion.1$fit - kde.dispersion.1$se.fit
kde.dispersion.1$upr <- kde.dispersion.1$fit + kde.dispersion.1$se.fit

range(cdr.2$min_cor)
new.df2 <- data.frame(min_cor = seq(min(cdr.2$min_cor), max(cdr.2$min_cor), by = 0.01))
fr.2 <- as.data.frame(predictSE.lme(fricmod_cdr2, new.df2))
fr.2$min_cor <- seq(min(cdr.2$min_cor), max(cdr.2$min_cor), by = .01)
fr.2$lwr <- fr.2$fit - fr.2$se.fit
fr.2$upr <- fr.2$fit + fr.2$se.fit
fe.2 <- as.data.frame(predictSE.lme(fevemod_cdr2, new.df2))
fe.2$min_cor <- seq(min(cdr.2$min_cor), max(cdr.2$min_cor), by = .01)
fe.2$lwr <- fe.2$fit - fe.2$se.fit
fe.2$upr <- fe.2$fit + fe.2$se.fit
fdis.2 <- as.data.frame(predictSE.lme(fdismod_cdr2, new.df2))
fdis.2$min_cor <- seq(min(cdr.2$min_cor), max(cdr.2$min_cor), by = .01)
fdis.2$lwr <- fdis.2$fit - fdis.2$se.fit
fdis.2$upr <- fdis.2$fit + fdis.2$se.fit
fdiv.2 <- as.data.frame(predictSE.lme(fdivmod_cdr2, new.df2))
fdiv.2$min_cor <- seq(min(cdr.2$min_cor), max(cdr.2$min_cor), by = .01)
fdiv.2$lwr <- fdiv.2$fit - fdiv.2$se.fit
fdiv.2$upr <- fdiv.2$fit + fdiv.2$se.fit
rq.2 <- as.data.frame(predictSE.lme(raoqmod_cdr2, new.df2))
rq.2$min_cor <- seq(min(cdr.2$min_cor), max(cdr.2$min_cor), by = .01)
rq.2$lwr <- rq.2$fit - rq.2$se.fit
rq.2$upr <- rq.2$fit + rq.2$se.fit
kde.alpha.2 <- as.data.frame(predictSE.lme(kde.alphamod_cdr2, new.df2))
kde.alpha.2$min_cor <- seq(min(cdr.2$min_cor), max(cdr.2$min_cor), by = .01)
kde.alpha.2$lwr <- kde.alpha.2$fit - kde.alpha.2$se.fit
kde.alpha.2$upr <- kde.alpha.2$fit + kde.alpha.2$se.fit
kde.evenness.2 <- as.data.frame(predictSE.lme(kde.evennessmod_cdr2, new.df2))
kde.evenness.2$min_cor <- seq(min(cdr.2$min_cor), max(cdr.2$min_cor), by = .01)
kde.evenness.2$lwr <- kde.evenness.2$fit - kde.evenness.2$se.fit
kde.evenness.2$upr <- kde.evenness.2$fit + kde.evenness.2$se.fit
kde.dispersion.2 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdr2, new.df2))
kde.dispersion.2$min_cor <- seq(min(cdr.2$min_cor), max(cdr.2$min_cor), by = .01)
kde.dispersion.2$lwr <- kde.dispersion.2$fit - kde.dispersion.2$se.fit
kde.dispersion.2$upr <- kde.dispersion.2$fit + kde.dispersion.2$se.fit

range(cdr.3$min_cor)
new.df3 <- data.frame(min_cor = seq(min(cdr.3$min_cor), max(cdr.3$min_cor), by = 0.01))
fr.3 <- as.data.frame(predictSE.lme(fricmod_cdr3, new.df3))
fr.3$min_cor <- seq(min(cdr.3$min_cor), max(cdr.3$min_cor), by = .01)
fr.3$lwr <- fr.3$fit - fr.3$se.fit
fr.3$upr <- fr.3$fit + fr.3$se.fit
fe.3 <- as.data.frame(predictSE.lme(fevemod_cdr3, new.df3))
fe.3$min_cor <- seq(min(cdr.3$min_cor), max(cdr.3$min_cor), by = .01)
fe.3$lwr <- fe.3$fit - fe.3$se.fit
fe.3$upr <- fe.3$fit + fe.3$se.fit
fdis.3 <- as.data.frame(predictSE.lme(fdismod_cdr3, new.df3))
fdis.3$min_cor <- seq(min(cdr.3$min_cor), max(cdr.3$min_cor), by = .01)
fdis.3$lwr <- fdis.3$fit - fdis.3$se.fit
fdis.3$upr <- fdis.3$fit + fdis.3$se.fit
fdiv.3 <- as.data.frame(predictSE.lme(fdivmod_cdr3, new.df3))
fdiv.3$min_cor <- seq(min(cdr.3$min_cor), max(cdr.3$min_cor), by = .01)
fdiv.3$lwr <- fdiv.3$fit - fdiv.3$se.fit
fdiv.3$upr <- fdiv.3$fit + fdiv.3$se.fit
rq.3 <- as.data.frame(predictSE.lme(raoqmod_cdr3, new.df3))
rq.3$min_cor <- seq(min(cdr.3$min_cor), max(cdr.3$min_cor), by = .01)
rq.3$lwr <- rq.3$fit - rq.3$se.fit
rq.3$upr <- rq.3$fit + rq.3$se.fit
kde.alpha.3 <- as.data.frame(predictSE.lme(kde.alphamod_cdr3, new.df3))
kde.alpha.3$min_cor <- seq(min(cdr.3$min_cor), max(cdr.3$min_cor), by = .01)
kde.alpha.3$lwr <- kde.alpha.3$fit - kde.alpha.3$se.fit
kde.alpha.3$upr <- kde.alpha.3$fit + kde.alpha.3$se.fit
kde.evenness.3 <- as.data.frame(predictSE.lme(kde.evennessmod_cdr3, new.df3))
kde.evenness.3$min_cor <- seq(min(cdr.3$min_cor), max(cdr.3$min_cor), by = .01)
kde.evenness.3$lwr <- kde.evenness.3$fit - kde.evenness.3$se.fit
kde.evenness.3$upr <- kde.evenness.3$fit + kde.evenness.3$se.fit
kde.dispersion.3 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdr3, new.df3))
kde.dispersion.3$min_cor <- seq(min(cdr.3$min_cor), max(cdr.3$min_cor), by = .01)
kde.dispersion.3$lwr <- kde.dispersion.3$fit - kde.dispersion.3$se.fit
kde.dispersion.3$upr <- kde.dispersion.3$fit + kde.dispersion.3$se.fit

range(cdr.4$min_cor)
new.df4 <- data.frame(min_cor = seq(min(cdr.4$min_cor), max(cdr.4$min_cor), by = 0.01))
fr.4 <- as.data.frame(predictSE.lme(fricmod_cdr4, new.df4))
fr.4$min_cor <- seq(min(cdr.4$min_cor), max(cdr.4$min_cor), by = .01)
fr.4$lwr <- fr.4$fit - fr.4$se.fit
fr.4$upr <- fr.4$fit + fr.4$se.fit
fe.4 <- as.data.frame(predictSE.lme(fevemod_cdr4, new.df4))
fe.4$min_cor <- seq(min(cdr.4$min_cor), max(cdr.4$min_cor), by = .01)
fe.4$lwr <- fe.4$fit - fe.4$se.fit
fe.4$upr <- fe.4$fit + fe.4$se.fit
fdis.4 <- as.data.frame(predictSE.lme(fdismod_cdr4, new.df4))
fdis.4$min_cor <- seq(min(cdr.4$min_cor), max(cdr.4$min_cor), by = .01)
fdis.4$lwr <- fdis.4$fit - fdis.4$se.fit
fdis.4$upr <- fdis.4$fit + fdis.4$se.fit
fdiv.4 <- as.data.frame(predictSE.lme(fdivmod_cdr4, new.df4))
fdiv.4$min_cor <- seq(min(cdr.4$min_cor), max(cdr.4$min_cor), by = .01)
fdiv.4$lwr <- fdiv.4$fit - fdiv.4$se.fit
fdiv.4$upr <- fdiv.4$fit + fdiv.4$se.fit
rq.4 <- as.data.frame(predictSE.lme(raoqmod_cdr4, new.df4))
rq.4$min_cor <- seq(min(cdr.4$min_cor), max(cdr.4$min_cor), by = .01)
rq.4$lwr <- rq.4$fit - rq.4$se.fit
rq.4$upr <- rq.4$fit + rq.4$se.fit
kde.alpha.4 <- as.data.frame(predictSE.lme(kde.alphamod_cdr4, new.df4))
kde.alpha.4$min_cor <- seq(min(cdr.4$min_cor), max(cdr.4$min_cor), by = .01)
kde.alpha.4$lwr <- kde.alpha.4$fit - kde.alpha.4$se.fit
kde.alpha.4$upr <- kde.alpha.4$fit + kde.alpha.4$se.fit
kde.evenness.4 <- as.data.frame(predictSE.lme(kde.evennessmod_cdr4, new.df4))
kde.evenness.4$min_cor <- seq(min(cdr.4$min_cor), max(cdr.4$min_cor), by = .01)
kde.evenness.4$lwr <- kde.evenness.4$fit - kde.evenness.4$se.fit
kde.evenness.4$upr <- kde.evenness.4$fit + kde.evenness.4$se.fit
kde.dispersion.4 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdr4, new.df4))
kde.dispersion.4$min_cor <- seq(min(cdr.4$min_cor), max(cdr.4$min_cor), by = .01)
kde.dispersion.4$lwr <- kde.dispersion.4$fit - kde.dispersion.4$se.fit
kde.dispersion.4$upr <- kde.dispersion.4$fit + kde.dispersion.4$se.fit


A <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fr.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= min_cor, y = fit), data = fr.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fr.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= min_cor, y = fit), data = fr.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fr.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= min_cor, y = fit), data = fr.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fr.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= min_cor, y = fit), data = fr.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = min_cor, y = FRic, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "FRich") +
  theme_pubr()

C <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fe.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= min_cor, y = fit), data = fe.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fe.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= min_cor, y = fit), data = fe.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fe.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= min_cor, y = fit), data = fe.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fe.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= min_cor, y = fit), data = fe.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = min_cor, y = FEve, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "FEve") +
  theme_pubr()

E <-ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdis.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= min_cor, y = fit), data = fdis.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdis.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= min_cor, y = fit), data = fdis.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdis.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= min_cor, y = fit), data = fdis.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdis.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= min_cor, y = fit), data = fdis.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = min_cor, y = FDis, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "FDis") +
  theme_pubr()

F <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdiv.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= min_cor, y = fit), data = fdiv.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdiv.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= min_cor, y = fit), data = fdiv.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdiv.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= min_cor, y = fit), data = fdiv.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdiv.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= min_cor, y = fit), data = fdiv.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = min_cor, y = FDiv, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "FDiv") +
  theme_pubr()

G <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = rq.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= min_cor, y = fit), data = rq.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = rq.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= min_cor, y = fit), data = rq.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = rq.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= min_cor, y = fit), data = rq.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = rq.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= min_cor, y = fit), data = rq.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = min_cor, y = RaoQ, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "Rao Q") +
  theme_pubr()

B <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.alpha.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.alpha.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.alpha.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.alpha.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.alpha.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.alpha.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.alpha.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.alpha.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = min_cor, y = kde.alpha, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "KDE Richness") +
  theme_pubr()

D <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.evenness.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.evenness.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.evenness.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.evenness.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.evenness.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.evenness.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.evenness.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.evenness.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = min_cor, y = kde.evenness, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "KDE Evenness") +
  theme_pubr()

H <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.dispersion.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.dispersion.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.dispersion.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.dispersion.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.dispersion.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.dispersion.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.dispersion.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.dispersion.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = min_cor, y = kde.dispersion, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "KDE Dispersion") +
  theme_pubr()

png(here("Figures/cdr_mincor.png"), height = 5, width = 9, units = 'in', res = 300)
ggarrange(plotlist = list(A, B, C, D, E, F, G, H), common.legend = TRUE,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
dev.off()


###SEV####

sev.black$community <- "SEV1"
sev.blue$community <- "SEV2"
sev.full <- rbind(sev.blue[,-10], sev.black[,-10])

sev.full <- sev.full %>% group_by(community, min_cor) %>%
  summarize (FRic = mean(FRic, na.rm = TRUE), FEve = mean(FEve, na.rm = TRUE), 
             FDis = mean(FDis, na.rm = TRUE), FDiv = mean(FDiv, na.rm = TRUE), 
             RaoQ = mean(RaoQ, na.rm = TRUE), kde.alpha = mean(kde.alpha), kde.evenness  = mean(kde.evenness), kde.dispersion = mean(kde.dispersion))


range(sev.blue$min_cor)
new.df5 <- data.frame(min_cor = seq(min(sev.blue$min_cor), max(sev.blue$min_cor), by = 0.01))
fr.blue <- as.data.frame(predictSE.lme(fricmod_blue, new.df5))
fr.blue$min_cor <- seq(min(sev.blue$min_cor), max(sev.blue$min_cor), by = 0.01)
fr.blue$lwr <- fr.blue$fit - fr.blue$se.fit
fr.blue$upr <- fr.blue$fit + fr.blue$se.fit
fe.blue <- as.data.frame(predictSE.lme(fevemod_blue, new.df5))
fe.blue$min_cor <- seq(min(sev.blue$min_cor), max(sev.blue$min_cor), by = 0.01)
fe.blue$lwr <- fe.blue$fit - fe.blue$se.fit
fe.blue$upr <- fe.blue$fit + fe.blue$se.fit
fdis.blue <- as.data.frame(predictSE.lme(fdismod_blue, new.df5))
fdis.blue$min_cor <- seq(min(sev.blue$min_cor), max(sev.blue$min_cor), by = 0.01)
fdis.blue$lwr <- fdis.blue$fit - fdis.blue$se.fit
fdis.blue$upr <- fdis.blue$fit + fdis.blue$se.fit
fdiv.blue <- as.data.frame(predictSE.lme(fdivmod_blue, new.df5))
fdiv.blue$min_cor <- seq(min(sev.blue$min_cor), max(sev.blue$min_cor), by = 0.01)
fdiv.blue$lwr <- fdiv.blue$fit - fdiv.blue$se.fit
fdiv.blue$upr <- fdiv.blue$fit + fdiv.blue$se.fit
rq.blue <- as.data.frame(predictSE.lme(raoqmod_blue, new.df5))
rq.blue$min_cor <- seq(min(sev.blue$min_cor), max(sev.blue$min_cor), by = 0.01)
rq.blue$lwr <- rq.blue$fit - rq.blue$se.fit
rq.blue$upr <- rq.blue$fit + rq.blue$se.fit
kde.alpha.blue <- as.data.frame(predictSE.lme(kde.alphamod_blue, new.df5))
kde.alpha.blue$min_cor <- seq(min(sev.blue$min_cor), max(sev.blue$min_cor), by = .01)
kde.alpha.blue$lwr <- kde.alpha.blue$fit - kde.alpha.blue$se.fit
kde.alpha.blue$upr <- kde.alpha.blue$fit + kde.alpha.blue$se.fit
kde.evenness.blue <- as.data.frame(predictSE.lme(kde.evennessmod_blue, new.df5))
kde.evenness.blue$min_cor <- seq(min(sev.blue$min_cor), max(sev.blue$min_cor), by = .01)
kde.evenness.blue$lwr <- kde.evenness.blue$fit - kde.evenness.blue$se.fit
kde.evenness.blue$upr <- kde.evenness.blue$fit + kde.evenness.blue$se.fit
kde.dispersion.blue <- as.data.frame(predictSE.lme(kde.dispersionmod_blue, new.df5))
kde.dispersion.blue$min_cor <- seq(min(sev.blue$min_cor), max(sev.blue$min_cor), by = .01)
kde.dispersion.blue$lwr <- kde.dispersion.blue$fit - kde.dispersion.blue$se.fit
kde.dispersion.blue$upr <- kde.dispersion.blue$fit + kde.dispersion.blue$se.fit

range(sev.black$min_cor)
new.df6 <- data.frame(min_cor = seq(min(sev.black$min_cor), max(sev.black$min_cor), by = 0.01))
fr.black <- as.data.frame(predictSE.lme(fricmod_black, new.df6))
fr.black$min_cor <- seq(min(sev.black$min_cor), max(sev.black$min_cor), by = .01)
fr.black$lwr <- fr.black$fit - fr.black$se.fit
fr.black$upr <- fr.black$fit + fr.black$se.fit
fe.black <- as.data.frame(predictSE.lme(fevemod_black, new.df6))
fe.black$min_cor <- seq(min(sev.black$min_cor), max(sev.black$min_cor), by = .01)
fe.black$lwr <- fe.black$fit - fe.black$se.fit
fe.black$upr <- fe.black$fit + fe.black$se.fit
fdis.black <- as.data.frame(predictSE.lme(fdismod_black, new.df6))
fdis.black$min_cor <- seq(min(sev.black$min_cor), max(sev.black$min_cor), by = .01)
fdis.black$lwr <- fdis.black$fit - fdis.black$se.fit
fdis.black$upr <- fdis.black$fit + fdis.black$se.fit
fdiv.black <- as.data.frame(predictSE.lme(fdivmod_black, new.df6))
fdiv.black$min_cor <- seq(min(sev.black$min_cor), max(sev.black$min_cor), by = .01)
fdiv.black$lwr <- fdiv.black$fit - fdiv.black$se.fit
fdiv.black$upr <- fdiv.black$fit + fdiv.black$se.fit
rq.black <- as.data.frame(predictSE.lme(raoqmod_black, new.df6))
rq.black$min_cor <- seq(min(sev.black$min_cor), max(sev.black$min_cor), by = .01)
rq.black$lwr <- rq.black$fit - rq.black$se.fit
rq.black$upr <- rq.black$fit + rq.black$se.fit
kde.alpha.black <- as.data.frame(predictSE.lme(kde.alphamod_black, new.df6))
kde.alpha.black$min_cor <- seq(min(sev.black$min_cor), max(sev.black$min_cor), by = .01)
kde.alpha.black$lwr <- kde.alpha.black$fit - kde.alpha.black$se.fit
kde.alpha.black$upr <- kde.alpha.black$fit + kde.alpha.black$se.fit
kde.evenness.black <- as.data.frame(predictSE.lme(kde.evennessmod_black, new.df6))
kde.evenness.black$min_cor <- seq(min(sev.black$min_cor), max(sev.black$min_cor), by = .01)
kde.evenness.black$lwr <- kde.evenness.black$fit - kde.evenness.black$se.fit
kde.evenness.black$upr <- kde.evenness.black$fit + kde.evenness.black$se.fit
kde.dispersion.black <- as.data.frame(predictSE.lme(kde.dispersionmod_black, new.df6))
kde.dispersion.black$min_cor <- seq(min(sev.black$min_cor), max(sev.black$min_cor), by = .01)
kde.dispersion.black$lwr <- kde.dispersion.black$fit - kde.dispersion.black$se.fit
kde.dispersion.black$upr <- kde.dispersion.black$fit + kde.dispersion.black$se.fit

AA <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fr.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= min_cor, y = fit), data = fr.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fr.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= min_cor, y = fit), data = fr.black, lwd = 2, color = "black") + 
  geom_point(aes(x = min_cor, y = FRic, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#000000", "navyblue"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "FRich") +
  theme_pubr()

CC <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fe.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= min_cor, y = fit), data = fe.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fe.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= min_cor, y = fit), data = fe.black, lwd = 2, color = "black") + 
  geom_point(aes(x = min_cor, y = FEve, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#000000", "navyblue"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "FEve") +
  theme_pubr()

EE <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdis.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= min_cor, y = fit), data = fdis.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdis.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= min_cor, y = fit), data = fdis.black, lwd = 2, color = "black") + 
  geom_point(aes(x = min_cor, y = FDis, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "FDis") +
  theme_pubr()

FF <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdiv.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= min_cor, y = fit), data = fdiv.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdiv.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= min_cor, y = fit), data = fdiv.black, lwd = 2, color = "black") + 
  geom_point(aes(x = min_cor, y = FDiv, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "FDiv") +
  theme_pubr()

GG <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = rq.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= min_cor, y = fit), data = rq.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = rq.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= min_cor, y = fit), data = rq.black, lwd = 2, color = "black") + 
  geom_point(aes(x = min_cor, y = RaoQ, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "Rao Q") +
  theme_pubr()

BB <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.alpha.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.alpha.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.alpha.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.alpha.black, lwd = 2, color = "black") + 
  geom_point(aes(x = min_cor, y = kde.alpha, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "KDE Richness") +
  theme_pubr()

DD <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.evenness.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.evenness.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.evenness.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.evenness.black, lwd = 2, color = "black") + 
  geom_point(aes(x = min_cor, y = kde.evenness, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "KDE Evenness") +
  theme_pubr()

HH <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.dispersion.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.dispersion.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.dispersion.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.dispersion.black, lwd = 2, color = "black") + 
  geom_point(aes(x = min_cor, y = kde.dispersion, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "KDE Dispersion") +
  theme_pubr()

png(here("Figures/Sev_mincorr.png"), height = 5, width = 9, units = 'in', res = 300)
ggarrange(plotlist = list(AA, BB, CC, DD, EE, FF, GG, HH), common.legend = TRUE, 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
dev.off()

png(here("Figures/min_corr.png"), height = 8, width = 13, units = 'in', res = 300)
ggarrange(plotlist = list(A, AA, B, BB, C, CC, D, DD, E, EE, F, FF, G, GG, H, HH), 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", 
                     "O", "P"), ncol = 4, nrow = 4, common.legend = TRUE, align = "hv")
dev.off()

