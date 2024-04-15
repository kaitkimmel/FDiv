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

cdre.1 <- read.csv(here("data/Cleaned/cdr1_sc_euc.csv"))
cdre.1 <- cdre.1[-which(is.na(cdre.1$mean_cor)),]
cdre.2 <- read.csv(here("data/Cleaned/cdr2_sc_euc.csv"))
cdre.2 <- cdre.2[-which(is.na(cdre.2$mean_cor)),]
cdre.3 <- read.csv(here("data/Cleaned/cdr3_sc_euc.csv"))
cdre.3 <- cdre.3[-which(is.na(cdre.3$mean_cor)),]
cdre.4 <- read.csv(here("data/Cleaned/cdr4_sc_euc.csv"))
cdre.4 <- cdre.4[-which(is.na(cdre.4$mean_cor)),]
seve.eblue <- read.csv(here("data/Cleaned/sevblue_sc_euc.csv"))
seve.eblue <- seve.eblue[-which(is.na(seve.eblue$mean_cor)),]
seve.eblack <- read.csv(here("data/Cleaned/sevblack_sc_euc.csv"))
seve.eblack <- seve.eblack[-which(is.na(seve.eblack$mean_cor)),]



#############################################################
################ seve #######################################
###########################################################

###########################################################
################### eblue ################################
##########################################################

###### Fit Metric ~ correlation models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ mean_cor, random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots, seve.eblue[-which(is.na(seve.eblue$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(FRic ~ 1, random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
AIC(mod, mod1, mod4) ### mod1
fricmod_eblue <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE richness
mod <- lme(kde.alpha ~ mean_cor, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots, seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(kde.alpha ~ 1, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) ### mod4
kde.alphamod_eblue <- lme(kde.alpha ~ 1, random = ~1|sev.blueplots,  data = seve.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

### FEve
mod <- lme(FEve ~ mean_cor, random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FEve ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(FEve ~ 1, random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FEve)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) ## mod1
fevemod_eblue <- lme(FEve ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE evenness
mod <- lme(kde.evenness ~ mean_cor, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots, seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(kde.evenness ~ 1, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) ### mod1
kde.evennessmod_eblue <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots,  data = seve.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

## FDis
mod <- lme(FDis ~ mean_cor, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FDis ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(FDis ~ 1, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) #mod1
fdismod_eblue <- lme(FDis ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots,  data = seve.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE dispersion
mod <- lme(kde.dispersion ~ mean_cor, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots, seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) ### mod4
kde.dispersionmod_eblue <- lme(kde.dispersion ~ 1, random = ~1|sev.blueplots,  data = seve.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

# FDiv
mod <- lme(FDiv ~ mean_cor, random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FDiv ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(FDiv ~ 1, random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FDiv)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) #mod1
anova(mod1, mod, mod4) #use mod
fdivmod_eblue <- lme(FDiv ~ mean_cor, random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

### Rao's Q
mod <- lme(RaoQ ~ mean_cor, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(RaoQ ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots),control = lmeControl(opt = "optim"))
mod4 <- lme(RaoQ ~ 1, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) #mod1
anova(mod1, mod) #use mod
raoqmod_eblue <- lme(RaoQ ~ mean_cor, random = ~1|sev.blueplots,  data = seve.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))


summary(fricmod_eblue)
summary(fevemod_eblue)
summary(fdismod_eblue)
summary(fdivmod_eblue)
summary(raoqmod_eblue)
summary(kde.alphamod_eblue)
summary(kde.evennessmod_eblue)
summary(kde.dispersionmod_eblue)

###########################################################
################### eblack ################################
##########################################################


###### Fit Metric ~ mean_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ mean_cor, random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots, seve.eblack[-which(is.na(seve.eblack$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(FRic ~ 1, random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) ### mod1 lowest
fricmod_eblack <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE richness
mod <- lme(kde.alpha ~ mean_cor, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots, seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(kde.alpha ~ 1, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) ### mod1
kde.alphamod_eblack <- lme(kde.alpha ~ mean_cor+ I(mean_cor^2), random = ~1|sev.blackplots,  data = seve.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))

### FEve
mod <- lme(FEve ~ mean_cor, random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FEve ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(FEve ~ 1, random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FEve)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) ## mod, mod4 similar - use mod 4
fevemod_eblack <- lme(FEve ~ 1, random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE evenness
mod <- lme(kde.evenness ~ mean_cor, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots, seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(kde.evenness ~ 1, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) ### mod4
kde.evennessmod_eblack <- lme(kde.evenness ~ 1, random = ~1|sev.blackplots,  data = seve.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))

## FDis
mod <- lme(FDis ~ mean_cor, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FDis ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(FDis ~ 1, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) #mod1
fdismod_eblack <- lme(FDis ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots,  data = seve.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE dispersion
mod <- lme(kde.dispersion ~ mean_cor, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots, seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) ### mod1
kde.dispersionmod_eblack <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots,  data = seve.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))

# FDiv
mod <- lme(FDiv ~ mean_cor, random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FDiv ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(FDiv ~ 1, random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FDiv)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) #mod1
fdivmod_eblack <- lme(FDiv ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

### Rao's Q
mod <- lme(RaoQ ~ mean_cor, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(RaoQ ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots),control = lmeControl(opt = "optim"))
mod4 <- lme(RaoQ ~ 1, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) #mod1
raoqmod_eblack <- lme(RaoQ ~ mean_cor+ I(mean_cor^2), random = ~1|sev.blackplots,  data = seve.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))


summary(fricmod_eblack)
summary(fevemod_eblack)
summary(fdismod_eblack)
summary(fdivmod_eblack)
summary(raoqmod_eblack)
summary(kde.alphamod_eblack)
summary(kde.evennessmod_eblack)
summary(kde.dispersionmod_eblack)

#########################################
######### CEDAR CREEK #####################
##########################################

###########################################################
################### cdre 1 ################################
##########################################################

###### Fit Metric ~ mean_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ mean_cor, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FRic ~ 1,  random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
fricmod_cdre1 <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ mean_cor, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.alpha ~ 1,  random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ### mod1
kde.alphamod_cdre1 <- lme(kde.alpha ~  mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ mean_cor, random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FEve ~ 1,  random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ## mod1
anova(mod1, mod)
fevemod_cdre1 <- lme(FEve ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ mean_cor, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.evenness ~ 1,  random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ### mod1
kde.evennessmod_cdre1 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ mean_cor, random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDis ~ 1,  random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) #mod1
fdismod_cdre1 <- lme(FDis ~ mean_cor+ I(mean_cor^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

##KDE dispersion
mod <- lme(kde.dispersion ~ mean_cor, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.dispersion ~ 1,  random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ###  mod1
kde.dispersionmod_cdre1 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ mean_cor, random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDiv ~ 1,  random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) #mod1
fdivmod_cdre1 <- lme(FDiv ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ mean_cor, random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(RaoQ ~ 1,  random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) #mod1 best fit
raoqmod_cdre1 <- lme(RaoQ ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod_cdre1)
summary(fevemod_cdre1)
summary(fdismod_cdre1)
summary(fdivmod_cdre1)
summary(raoqmod_cdre1)
summary(kde.alphamod_cdre1)
summary(kde.evennessmod_cdre1)
summary(kde.dispersionmod_cdre1)


###########################################################
################### cdre 2 ################################
##########################################################
###### Fit Metric ~ mean_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ mean_cor, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FRic ~ 1,  random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ### mod1
fricmod_cdre2 <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ mean_cor, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.alpha ~ 1,  random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ### mod1
kde.alphamod_cdre2 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ mean_cor, random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FEve ~ 1,  random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ## mod1
fevemod_cdre2 <- lme(FEve ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ mean_cor, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.evenness ~ 1,  random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ### mod1
kde.evennessmod_cdre2 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ mean_cor, random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDis ~ 1,  random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) #mod1
fdismod_cdre2 <- lme(FDis ~ mean_cor+ I(mean_cor^2), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

##KDE dispersion
mod <- lme(kde.dispersion ~ mean_cor, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.dispersion ~ 1,  random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ### mod1
kde.dispersionmod_cdre2 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ mean_cor, random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDiv ~ 1,  random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) #mod1
fdivmod_cdre2 <- lme(FDiv ~ mean_cor+ I(mean_cor^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ mean_cor, random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(RaoQ ~ 1,  random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) #mod1
raoqmod_cdre2 <- lme(RaoQ ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod_cdre2)
summary(fevemod_cdre2)
summary(fdismod_cdre2)
summary(fdivmod_cdre2)
summary(raoqmod_cdre2)
summary(kde.alphamod_cdre2)
summary(kde.evennessmod_cdre2)
summary(kde.dispersionmod_cdre2)

###########################################################
################### cdre 3 ################################
##########################################################
###### Fit Metric ~ mean_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ mean_cor, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FRic ~ 1,  random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ### mod1
fricmod_cdre3 <- lme(FRic ~ mean_cor+ I(mean_cor^2), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ mean_cor, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.alpha ~ 1,  random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ### mod1
kde.alphamod_cdre3 <- lme(kde.alpha ~ mean_cor+ I(mean_cor^2), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ mean_cor, random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FEve ~ 1,  random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ## mod
fevemod_cdre3 <- lme(FEve ~ mean_cor, random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ mean_cor, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.evenness ~ 1,  random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ### mod4
kde.evennessmod_cdre3 <- lme(kde.evenness ~ 1, random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ mean_cor, random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDis ~ 1,  random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) #mod
anova(mod, mod4) 
fdismod_cdre3 <- lme(FDis ~ mean_cor, random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

### KDE dispersion
mod <- lme(kde.dispersion ~ mean_cor, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.dispersion ~ 1,  random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ### mod1
kde.dispersionmod_cdre3 <- lme(kde.dispersion ~ mean_cor+ I(mean_cor^2), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ mean_cor, random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDiv ~ 1,  random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) #mod1
fdivmod_cdre3 <- lme(FDiv ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ mean_cor, random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(RaoQ ~ 1,  random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) #mod
raoqmod_cdre3 <- lme(RaoQ ~ mean_cor, random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod_cdre3)
summary(fevemod_cdre3)
summary(fdismod_cdre3)
summary(fdivmod_cdre3)
summary(raoqmod_cdre3)
summary(kde.alphamod_cdre3)
summary(kde.evennessmod_cdre3)
summary(kde.dispersionmod_cdre3)

###########################################################
################### cdre 4 ################################
##########################################################
###### Fit Metric ~ mean_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ mean_cor, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FRic ~ 1,  random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ### mod1
fricmod_cdre4 <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ mean_cor, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.alpha ~ 1,  random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ### mod
kde.alphamod_cdre4 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ mean_cor, random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FEve ~ 1,  random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ## mod  
fevemod_cdre4 <- lme(FEve ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ mean_cor, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.evenness ~ 1,  random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ### mod1
anova(mod1, mod4)
kde.evennessmod_cdre4 <- lme(kde.evenness ~ 1, random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ mean_cor, random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDis ~ 1,  random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) #mod1
fdismod_cdre4 <- lme(FDis ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE dispersion
mod <- lme(kde.dispersion ~ mean_cor, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.dispersion ~ 1,  random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ### mod1
kde.dispersionmod_cdre4 <- lme(kde.dispersion ~  mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ mean_cor, random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDiv ~ 1,  random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) #mod
anova(mod, mod4) #use mod
fdivmod_cdre4 <- lme(FDiv ~ mean_cor, random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ mean_cor, random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(RaoQ ~ 1,  random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) #mod1
raoqmod_cdre4 <- lme(RaoQ ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod_cdre4)
summary(fevemod_cdre4)
summary(fdismod_cdre4)
summary(fdivmod_cdre4)
summary(raoqmod_cdre4)
summary(kde.alphamod_cdre4)
summary(kde.evennessmod_cdre4)
summary(kde.dispersionmod_cdre4)





summary(fricmod_eblue)
summary(fricmod_eblack)
summary(fricmod_cdre1)
summary(fricmod_cdre2)
summary(fricmod_cdre3)
summary(fricmod_cdre4)

summary(kde.alphamod_eblue)
summary(kde.alphamod_eblack)
summary(kde.alphamod_cdre1)
summary(kde.alphamod_cdre2)
summary(kde.alphamod_cdre3)
summary(kde.alphamod_cdre4)

summary(fevemod_eblue)
summary(fevemod_eblack)
summary(fevemod_cdre1)
summary(fevemod_cdre2)
summary(fevemod_cdre3)
summary(fevemod_cdre4)

summary(kde.evennessmod_eblue)
summary(kde.evennessmod_eblack)
summary(kde.evennessmod_cdre1)
summary(kde.evennessmod_cdre2)
summary(kde.evennessmod_cdre3)
summary(kde.evennessmod_cdre4)

summary(fdismod_eblue)
summary(fdismod_eblack)
summary(fdismod_cdre1)
summary(fdismod_cdre2)
summary(fdismod_cdre3)
summary(fdismod_cdre4)

summary(kde.dispersionmod_eblue)
summary(kde.dispersionmod_eblack)
summary(kde.dispersionmod_cdre1)
summary(kde.dispersionmod_cdre2)
summary(kde.dispersionmod_cdre3)
summary(kde.dispersionmod_cdre4)

summary(fdivmod_eblue)
summary(fdivmod_eblack)
summary(fdivmod_cdre1)
summary(fdivmod_cdre2)
summary(fdivmod_cdre3)
summary(fdivmod_cdre4)

summary(raoqmod_eblue)
summary(raoqmod_eblack)
summary(raoqmod_cdre1)
summary(raoqmod_cdre2)
summary(raoqmod_cdre3)
summary(raoqmod_cdre4)




################### GRAPHS ##############################


# create combined data set of raw data
cdre.1$community <- "CDR1"
cdre.2$community <- "CDR2"
cdre.3$community <- "CDR3"
cdre.4$community <- "CDR4"
cdre.full <- rbind(cdre.1, cdre.2, cdre.3, cdre.4)

cdre.full <- cdre.full %>% group_by(community, mean_cor) %>%
  summarize (FRic = mean(FRic), FEve = mean(FEve), FDis = mean(FDis),
             FDiv = mean(FDiv), RaoQ = mean(RaoQ), RaoQ = mean(RaoQ), kde.alpha = mean(kde.alpha), kde.evenness  = mean(kde.evenness), kde.dispersion = mean(kde.dispersion))

##############################################
############## cdre 1 ########################
#############################################
range(cdre.1$mean_cor)
new.df1 <- data.frame(mean_cor = seq(min(cdre.1$mean_cor), max(cdre.1$mean_cor), by = 0.01))
fr.e1 <- as.data.frame(predictSE.lme(fricmod_cdre1, new.df1))
fr.e1$mean_cor <- seq(min(cdre.1$mean_cor), max(cdre.1$mean_cor), by = 0.01)
fr.e1$lwr <- fr.e1$fit - fr.e1$se.fit
fr.e1$upr <- fr.e1$fit + fr.e1$se.fit
fe.e1 <- as.data.frame(predictSE.lme(fevemod_cdre1, new.df1))
fe.e1$mean_cor <- seq(min(cdre.1$mean_cor), max(cdre.1$mean_cor), by = 0.01)
fe.e1$lwr <- fe.e1$fit - fe.e1$se.fit
fe.e1$upr <- fe.e1$fit + fe.e1$se.fit
fdis.e1 <- as.data.frame(predictSE.lme(fdismod_cdre1, new.df1))
fdis.e1$mean_cor <- seq(min(cdre.1$mean_cor), max(cdre.1$mean_cor), by = 0.01)
fdis.e1$lwr <- fdis.e1$fit - fdis.e1$se.fit
fdis.e1$upr <- fdis.e1$fit + fdis.e1$se.fit
fdiv.e1 <- as.data.frame(predictSE.lme(fdivmod_cdre1, new.df1))
fdiv.e1$mean_cor <- seq(min(cdre.1$mean_cor), max(cdre.1$mean_cor), by = 0.01)
fdiv.e1$lwr <- fdiv.e1$fit - fdiv.e1$se.fit
fdiv.e1$upr <- fdiv.e1$fit + fdiv.e1$se.fit
rq.e1 <- as.data.frame(predictSE.lme(raoqmod_cdre1, new.df1))
rq.e1$mean_cor <- seq(min(cdre.1$mean_cor), max(cdre.1$mean_cor), by = 0.01)
rq.e1$lwr <- rq.e1$fit - rq.e1$se.fit
rq.e1$upr <- rq.e1$fit + rq.e1$se.fit
kde.alpha.e1 <- as.data.frame(predictSE.lme(kde.alphamod_cdre1, new.df1))
kde.alpha.e1$mean_cor <- seq(min(cdre.1$mean_cor), max(cdre.1$mean_cor), by = .01)
kde.alpha.e1$lwr <- kde.alpha.e1$fit - kde.alpha.e1$se.fit
kde.alpha.e1$upr <- kde.alpha.e1$fit + kde.alpha.e1$se.fit
kde.evenness.e1 <- as.data.frame(predictSE.lme(kde.evennessmod_cdre1, new.df1))
kde.evenness.e1$mean_cor <- seq(min(cdre.1$mean_cor), max(cdre.1$mean_cor), by = .01)
kde.evenness.e1$lwr <- kde.evenness.e1$fit - kde.evenness.e1$se.fit
kde.evenness.e1$upr <- kde.evenness.e1$fit + kde.evenness.e1$se.fit
kde.dispersion.e1 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdre1, new.df1))
kde.dispersion.e1$mean_cor <- seq(min(cdre.1$mean_cor), max(cdre.1$mean_cor), by = .01)
kde.dispersion.e1$lwr <- kde.dispersion.e1$fit - kde.dispersion.e1$se.fit
kde.dispersion.e1$upr <- kde.dispersion.e1$fit + kde.dispersion.e1$se.fit

range(cdre.2$mean_cor)
new.df2 <- data.frame(mean_cor = seq(min(cdre.2$mean_cor), max(cdre.2$mean_cor), by = 0.01))
fr.e2 <- as.data.frame(predictSE.lme(fricmod_cdre2, new.df2))
fr.e2$mean_cor <- seq(min(cdre.2$mean_cor), max(cdre.2$mean_cor), by = .01)
fr.e2$lwr <- fr.e2$fit - fr.e2$se.fit
fr.e2$upr <- fr.e2$fit + fr.e2$se.fit
fe.e2 <- as.data.frame(predictSE.lme(fevemod_cdre2, new.df2))
fe.e2$mean_cor <- seq(min(cdre.2$mean_cor), max(cdre.2$mean_cor), by = .01)
fe.e2$lwr <- fe.e2$fit - fe.e2$se.fit
fe.e2$upr <- fe.e2$fit + fe.e2$se.fit
fdis.e2 <- as.data.frame(predictSE.lme(fdismod_cdre2, new.df2))
fdis.e2$mean_cor <- seq(min(cdre.2$mean_cor), max(cdre.2$mean_cor), by = .01)
fdis.e2$lwr <- fdis.e2$fit - fdis.e2$se.fit
fdis.e2$upr <- fdis.e2$fit + fdis.e2$se.fit
fdiv.e2 <- as.data.frame(predictSE.lme(fdivmod_cdre2, new.df2))
fdiv.e2$mean_cor <- seq(min(cdre.2$mean_cor), max(cdre.2$mean_cor), by = .01)
fdiv.e2$lwr <- fdiv.e2$fit - fdiv.e2$se.fit
fdiv.e2$upr <- fdiv.e2$fit + fdiv.e2$se.fit
rq.e2 <- as.data.frame(predictSE.lme(raoqmod_cdre2, new.df2))
rq.e2$mean_cor <- seq(min(cdre.2$mean_cor), max(cdre.2$mean_cor), by = .01)
rq.e2$lwr <- rq.e2$fit - rq.e2$se.fit
rq.e2$upr <- rq.e2$fit + rq.e2$se.fit
kde.alpha.e2 <- as.data.frame(predictSE.lme(kde.alphamod_cdre2, new.df2))
kde.alpha.e2$mean_cor <- seq(min(cdre.2$mean_cor), max(cdre.2$mean_cor), by = .01)
kde.alpha.e2$lwr <- kde.alpha.e2$fit - kde.alpha.e2$se.fit
kde.alpha.e2$upr <- kde.alpha.e2$fit + kde.alpha.e2$se.fit
kde.evenness.e2 <- as.data.frame(predictSE.lme(kde.evennessmod_cdre2, new.df2))
kde.evenness.e2$mean_cor <- seq(min(cdre.2$mean_cor), max(cdre.2$mean_cor), by = .01)
kde.evenness.e2$lwr <- kde.evenness.e2$fit - kde.evenness.e2$se.fit
kde.evenness.e2$upr <- kde.evenness.e2$fit + kde.evenness.e2$se.fit
kde.dispersion.e2 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdre2, new.df2))
kde.dispersion.e2$mean_cor <- seq(min(cdre.2$mean_cor), max(cdre.2$mean_cor), by = .01)
kde.dispersion.e2$lwr <- kde.dispersion.e2$fit - kde.dispersion.e2$se.fit
kde.dispersion.e2$upr <- kde.dispersion.e2$fit + kde.dispersion.e2$se.fit

range(cdre.3$mean_cor)
new.df3 <- data.frame(mean_cor = seq(min(cdre.3$mean_cor), max(cdre.3$mean_cor), by = 0.01))
fr.e3 <- as.data.frame(predictSE.lme(fricmod_cdre3, new.df3))
fr.e3$mean_cor <- seq(min(cdre.3$mean_cor), max(cdre.3$mean_cor), by = .01)
fr.e3$lwr <- fr.e3$fit - fr.e3$se.fit
fr.e3$upr <- fr.e3$fit + fr.e3$se.fit
fe.e3 <- as.data.frame(predictSE.lme(fevemod_cdre3, new.df3))
fe.e3$mean_cor <- seq(min(cdre.3$mean_cor), max(cdre.3$mean_cor), by = .01)
fe.e3$lwr <- fe.e3$fit - fe.e3$se.fit
fe.e3$upr <- fe.e3$fit + fe.e3$se.fit
fdis.e3 <- as.data.frame(predictSE.lme(fdismod_cdre3, new.df3))
fdis.e3$mean_cor <- seq(min(cdre.3$mean_cor), max(cdre.3$mean_cor), by = .01)
fdis.e3$lwr <- fdis.e3$fit - fdis.e3$se.fit
fdis.e3$upr <- fdis.e3$fit + fdis.e3$se.fit
fdiv.e3 <- as.data.frame(predictSE.lme(fdivmod_cdre3, new.df3))
fdiv.e3$mean_cor <- seq(min(cdre.3$mean_cor), max(cdre.3$mean_cor), by = .01)
fdiv.e3$lwr <- fdiv.e3$fit - fdiv.e3$se.fit
fdiv.e3$upr <- fdiv.e3$fit + fdiv.e3$se.fit
rq.e3 <- as.data.frame(predictSE.lme(raoqmod_cdre3, new.df3))
rq.e3$mean_cor <- seq(min(cdre.3$mean_cor), max(cdre.3$mean_cor), by = .01)
rq.e3$lwr <- rq.e3$fit - rq.e3$se.fit
rq.e3$upr <- rq.e3$fit + rq.e3$se.fit
kde.alpha.e3 <- as.data.frame(predictSE.lme(kde.alphamod_cdre3, new.df3))
kde.alpha.e3$mean_cor <- seq(min(cdre.3$mean_cor), max(cdre.3$mean_cor), by = .01)
kde.alpha.e3$lwr <- kde.alpha.e3$fit - kde.alpha.e3$se.fit
kde.alpha.e3$upr <- kde.alpha.e3$fit + kde.alpha.e3$se.fit
kde.evenness.e3 <- as.data.frame(predictSE.lme(kde.evennessmod_cdre3, new.df3))
kde.evenness.e3$mean_cor <- seq(min(cdre.3$mean_cor), max(cdre.3$mean_cor), by = .01)
kde.evenness.e3$lwr <- kde.evenness.e3$fit - kde.evenness.e3$se.fit
kde.evenness.e3$upr <- kde.evenness.e3$fit + kde.evenness.e3$se.fit
kde.dispersion.e3 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdre3, new.df3))
kde.dispersion.e3$mean_cor <- seq(min(cdre.3$mean_cor), max(cdre.3$mean_cor), by = .01)
kde.dispersion.e3$lwr <- kde.dispersion.e3$fit - kde.dispersion.e3$se.fit
kde.dispersion.e3$upr <- kde.dispersion.e3$fit + kde.dispersion.e3$se.fit

range(cdre.4$mean_cor)
new.df4 <- data.frame(mean_cor = seq(min(cdre.4$mean_cor), max(cdre.4$mean_cor), by = 0.01))
fr.e4 <- as.data.frame(predictSE.lme(fricmod_cdre4, new.df4))
fr.e4$mean_cor <- seq(min(cdre.4$mean_cor), max(cdre.4$mean_cor), by = .01)
fr.e4$lwr <- fr.e4$fit - fr.e4$se.fit
fr.e4$upr <- fr.e4$fit + fr.e4$se.fit
fe.e4 <- as.data.frame(predictSE.lme(fevemod_cdre4, new.df4))
fe.e4$mean_cor <- seq(min(cdre.4$mean_cor), max(cdre.4$mean_cor), by = .01)
fe.e4$lwr <- fe.e4$fit - fe.e4$se.fit
fe.e4$upr <- fe.e4$fit + fe.e4$se.fit
fdis.e4 <- as.data.frame(predictSE.lme(fdismod_cdre4, new.df4))
fdis.e4$mean_cor <- seq(min(cdre.4$mean_cor), max(cdre.4$mean_cor), by = .01)
fdis.e4$lwr <- fdis.e4$fit - fdis.e4$se.fit
fdis.e4$upr <- fdis.e4$fit + fdis.e4$se.fit
fdiv.e4 <- as.data.frame(predictSE.lme(fdivmod_cdre4, new.df4))
fdiv.e4$mean_cor <- seq(min(cdre.4$mean_cor), max(cdre.4$mean_cor), by = .01)
fdiv.e4$lwr <- fdiv.e4$fit - fdiv.e4$se.fit
fdiv.e4$upr <- fdiv.e4$fit + fdiv.e4$se.fit
rq.e4 <- as.data.frame(predictSE.lme(raoqmod_cdre4, new.df4))
rq.e4$mean_cor <- seq(min(cdre.4$mean_cor), max(cdre.4$mean_cor), by = .01)
rq.e4$lwr <- rq.e4$fit - rq.e4$se.fit
rq.e4$upr <- rq.e4$fit + rq.e4$se.fit
kde.alpha.e4 <- as.data.frame(predictSE.lme(kde.alphamod_cdre4, new.df4))
kde.alpha.e4$mean_cor <- seq(min(cdre.4$mean_cor), max(cdre.4$mean_cor), by = .01)
kde.alpha.e4$lwr <- kde.alpha.e4$fit - kde.alpha.e4$se.fit
kde.alpha.e4$upr <- kde.alpha.e4$fit + kde.alpha.e4$se.fit
kde.evenness.e4 <- as.data.frame(predictSE.lme(kde.evennessmod_cdre4, new.df4))
kde.evenness.e4$mean_cor <- seq(min(cdre.4$mean_cor), max(cdre.4$mean_cor), by = .01)
kde.evenness.e4$lwr <- kde.evenness.e4$fit - kde.evenness.e4$se.fit
kde.evenness.e4$upr <- kde.evenness.e4$fit + kde.evenness.e4$se.fit
kde.dispersion.e4 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdre4, new.df4))
kde.dispersion.e4$mean_cor <- seq(min(cdre.4$mean_cor), max(cdre.4$mean_cor), by = .01)
kde.dispersion.e4$lwr <- kde.dispersion.e4$fit - kde.dispersion.e4$se.fit
kde.dispersion.e4$upr <- kde.dispersion.e4$fit + kde.dispersion.e4$se.fit


cdre.full <- cdre.full %>%
  filter(row_number() %% 5 == 1)

A <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fr.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= mean_cor, y = fit), data = fr.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fr.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= mean_cor, y = fit), data = fr.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fr.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= mean_cor, y = fit), data = fr.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fr.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= mean_cor, y = fit), data = fr.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = mean_cor, y = FRic, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Mean Trait-Trait Correlation", y = "FRich") +
  theme_pubr()

C <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fe.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= mean_cor, y = fit), data = fe.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fe.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= mean_cor, y = fit), data = fe.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fe.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= mean_cor, y = fit), data = fe.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fe.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= mean_cor, y = fit), data = fe.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = mean_cor, y = FEve, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Mean Trait-Trait Correlation", y = "FEve") +
  theme_pubr()

E <-ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdis.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdis.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdis.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdis.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdis.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdis.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdis.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdis.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = mean_cor, y = FDis, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Mean Trait-Trait Correlation", y = "FDis") +
  theme_pubr()

F <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdiv.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdiv.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdiv.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdiv.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdiv.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdiv.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdiv.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdiv.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = mean_cor, y = FDiv, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Mean Trait-Trait Correlation", y = "FDiv") +
  theme_pubr()

G <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = rq.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= mean_cor, y = fit), data = rq.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = rq.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= mean_cor, y = fit), data = rq.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = rq.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= mean_cor, y = fit), data = rq.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = rq.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= mean_cor, y = fit), data = rq.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = mean_cor, y = RaoQ, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Mean Trait-Trait Correlation", y = "Rao Q") +
  theme_pubr()

B <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.alpha.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.alpha.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.alpha.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.alpha.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.alpha.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.alpha.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.alpha.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.alpha.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = mean_cor, y = kde.alpha, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Mean Trait-Trait Correlation", y = "KDE Richness") +
  theme_pubr()

C <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.evenness.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.evenness.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.evenness.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.evenness.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.evenness.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.evenness.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.evenness.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.evenness.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = mean_cor, y = kde.evenness, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Mean Trait-Trait Correlation", y = "KDE Evenness") +
  theme_pubr()

H <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.dispersion.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.dispersion.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.dispersion.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.dispersion.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.dispersion.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.dispersion.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.dispersion.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.dispersion.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = mean_cor, y = kde.dispersion, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Mean Trait-Trait Correlation", y = "KDE Dispersion") +
  theme_pubr()

png(here("Figures/cdre_meancor_euc.png"), height = 5, width = 9, units = 'in', res = 300)
ggarrange(plotlist = list(A, B, C, D, E, F, G, H), common.legend = TRUE,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
dev.off()


###seve####


seve.eblack$community <- "SEV1"
seve.eblue$community <- "SEV2"
seve.full <- rbind(seve.eblue[,-10], seve.eblack[,-10])

seve.full <- seve.full %>% group_by(community, mean_cor) %>%
  summarize (FRic = mean(FRic, na.rm = TRUE), FEve = mean(FEve, na.rm = TRUE), 
             FDis = mean(FDis, na.rm = TRUE), FDiv = mean(FDiv, na.rm = TRUE), 
             RaoQ = mean(RaoQ, na.rm = TRUE), kde.alpha = mean(kde.alpha), kde.evenness  = mean(kde.evenness), kde.dispersion = mean(kde.dispersion))


range(seve.eblue$mean_cor)
new.df5 <- data.frame(mean_cor = seq(min(seve.eblue$mean_cor), max(seve.eblue$mean_cor), by = 0.01))
fr.eblue <- as.data.frame(predictSE.lme(fricmod_eblue, new.df5))
fr.eblue$mean_cor <- seq(min(seve.eblue$mean_cor), max(seve.eblue$mean_cor), by = 0.01)
fr.eblue$lwr <- fr.eblue$fit - fr.eblue$se.fit
fr.eblue$upr <- fr.eblue$fit + fr.eblue$se.fit
fe.eblue <- as.data.frame(predictSE.lme(fevemod_eblue, new.df5))
fe.eblue$mean_cor <- seq(min(seve.eblue$mean_cor), max(seve.eblue$mean_cor), by = 0.01)
fe.eblue$lwr <- fe.eblue$fit - fe.eblue$se.fit
fe.eblue$upr <- fe.eblue$fit + fe.eblue$se.fit
fdis.eblue <- as.data.frame(predictSE.lme(fdismod_eblue, new.df5))
fdis.eblue$mean_cor <- seq(min(seve.eblue$mean_cor), max(seve.eblue$mean_cor), by = 0.01)
fdis.eblue$lwr <- fdis.eblue$fit - fdis.eblue$se.fit
fdis.eblue$upr <- fdis.eblue$fit + fdis.eblue$se.fit
fdiv.eblue <- as.data.frame(predictSE.lme(fdivmod_eblue, new.df5))
fdiv.eblue$mean_cor <- seq(min(seve.eblue$mean_cor), max(seve.eblue$mean_cor), by = 0.01)
fdiv.eblue$lwr <- fdiv.eblue$fit - fdiv.eblue$se.fit
fdiv.eblue$upr <- fdiv.eblue$fit + fdiv.eblue$se.fit
rq.eblue <- as.data.frame(predictSE.lme(raoqmod_eblue, new.df5))
rq.eblue$mean_cor <- seq(min(seve.eblue$mean_cor), max(seve.eblue$mean_cor), by = 0.01)
rq.eblue$lwr <- rq.eblue$fit - rq.eblue$se.fit
rq.eblue$upr <- rq.eblue$fit + rq.eblue$se.fit
kde.alpha.eblue <- as.data.frame(predictSE.lme(kde.alphamod_eblue, new.df5))
kde.alpha.eblue$mean_cor <- seq(min(seve.eblue$mean_cor), max(seve.eblue$mean_cor), by = .01)
kde.alpha.eblue$lwr <- kde.alpha.eblue$fit - kde.alpha.eblue$se.fit
kde.alpha.eblue$upr <- kde.alpha.eblue$fit + kde.alpha.eblue$se.fit
kde.evenness.eblue <- as.data.frame(predictSE.lme(kde.evennessmod_eblue, new.df5))
kde.evenness.eblue$mean_cor <- seq(min(seve.eblue$mean_cor), max(seve.eblue$mean_cor), by = .01)
kde.evenness.eblue$lwr <- kde.evenness.eblue$fit - kde.evenness.eblue$se.fit
kde.evenness.eblue$upr <- kde.evenness.eblue$fit + kde.evenness.eblue$se.fit
kde.dispersion.eblue <- as.data.frame(predictSE.lme(kde.dispersionmod_eblue, new.df5))
kde.dispersion.eblue$mean_cor <- seq(min(seve.eblue$mean_cor), max(seve.eblue$mean_cor), by = .01)
kde.dispersion.eblue$lwr <- kde.dispersion.eblue$fit - kde.dispersion.eblue$se.fit
kde.dispersion.eblue$upr <- kde.dispersion.eblue$fit + kde.dispersion.eblue$se.fit

range(seve.eblack$mean_cor)
new.df6 <- data.frame(mean_cor = seq(min(seve.eblack$mean_cor), max(seve.eblack$mean_cor), by = 0.01))
fr.eblack <- as.data.frame(predictSE.lme(fricmod_eblack, new.df6))
fr.eblack$mean_cor <- seq(min(seve.eblack$mean_cor), max(seve.eblack$mean_cor), by = .01)
fr.eblack$lwr <- fr.eblack$fit - fr.eblack$se.fit
fr.eblack$upr <- fr.eblack$fit + fr.eblack$se.fit
fe.eblack <- as.data.frame(predictSE.lme(fevemod_eblack, new.df6))
fe.eblack$mean_cor <- seq(min(seve.eblack$mean_cor), max(seve.eblack$mean_cor), by = .01)
fe.eblack$lwr <- fe.eblack$fit - fe.eblack$se.fit
fe.eblack$upr <- fe.eblack$fit + fe.eblack$se.fit
fdis.eblack <- as.data.frame(predictSE.lme(fdismod_eblack, new.df6))
fdis.eblack$mean_cor <- seq(min(seve.eblack$mean_cor), max(seve.eblack$mean_cor), by = .01)
fdis.eblack$lwr <- fdis.eblack$fit - fdis.eblack$se.fit
fdis.eblack$upr <- fdis.eblack$fit + fdis.eblack$se.fit
fdiv.eblack <- as.data.frame(predictSE.lme(fdivmod_eblack, new.df6))
fdiv.eblack$mean_cor <- seq(min(seve.eblack$mean_cor), max(seve.eblack$mean_cor), by = .01)
fdiv.eblack$lwr <- fdiv.eblack$fit - fdiv.eblack$se.fit
fdiv.eblack$upr <- fdiv.eblack$fit + fdiv.eblack$se.fit
rq.eblack <- as.data.frame(predictSE.lme(raoqmod_eblack, new.df6))
rq.eblack$mean_cor <- seq(min(seve.eblack$mean_cor), max(seve.eblack$mean_cor), by = .01)
rq.eblack$lwr <- rq.eblack$fit - rq.eblack$se.fit
rq.eblack$upr <- rq.eblack$fit + rq.eblack$se.fit
kde.alpha.eblack <- as.data.frame(predictSE.lme(kde.alphamod_eblack, new.df6))
kde.alpha.eblack$mean_cor <- seq(min(seve.eblack$mean_cor), max(seve.eblack$mean_cor), by = .01)
kde.alpha.eblack$lwr <- kde.alpha.eblack$fit - kde.alpha.eblack$se.fit
kde.alpha.eblack$upr <- kde.alpha.eblack$fit + kde.alpha.eblack$se.fit
kde.evenness.eblack <- as.data.frame(predictSE.lme(kde.evennessmod_eblack, new.df6))
kde.evenness.eblack$mean_cor <- seq(min(seve.eblack$mean_cor), max(seve.eblack$mean_cor), by = .01)
kde.evenness.eblack$lwr <- kde.evenness.eblack$fit - kde.evenness.eblack$se.fit
kde.evenness.eblack$upr <- kde.evenness.eblack$fit + kde.evenness.eblack$se.fit
kde.dispersion.eblack <- as.data.frame(predictSE.lme(kde.dispersionmod_eblack, new.df6))
kde.dispersion.eblack$mean_cor <- seq(min(seve.eblack$mean_cor), max(seve.eblack$mean_cor), by = .01)
kde.dispersion.eblack$lwr <- kde.dispersion.eblack$fit - kde.dispersion.eblack$se.fit
kde.dispersion.eblack$upr <- kde.dispersion.eblack$fit + kde.dispersion.eblack$se.fit

seve.full <- seve.full %>%
  filter(row_number() %% 5 == 1)
AA <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fr.eblue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= mean_cor, y = fit), data = fr.eblue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fr.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= mean_cor, y = fit), data = fr.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = mean_cor, y = FRic, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("#000000", "navyblue"), name = "Community") +
  labs(x = "Mean Trait-Trait Correlation", y = "FRich") +
  theme_pubr()

CC <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fe.eblue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= mean_cor, y = fit), data = fe.eblue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fe.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= mean_cor, y = fit), data = fe.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = mean_cor, y = FEve, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("#000000", "navyblue"), name = "Community") +
  labs(x = "Mean Trait-Trait Correlation", y = "FEve") +
  theme_pubr()

EE <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdis.eblue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdis.eblue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdis.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdis.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = mean_cor, y = FDis, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Mean Trait-Trait Correlation", y = "FDis") +
  theme_pubr()

FF <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdiv.eblue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdiv.eblue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdiv.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdiv.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = mean_cor, y = FDiv, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Mean Trait-Trait Correlation", y = "FDiv") +
  theme_pubr()

GG <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = rq.eblue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= mean_cor, y = fit), data = rq.eblue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = rq.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= mean_cor, y = fit), data = rq.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = mean_cor, y = RaoQ, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Mean Trait-Trait Correlation", y = "Rao Q") +
  theme_pubr()

BB <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.alpha.eblue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.alpha.eblue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.alpha.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.alpha.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = mean_cor, y = kde.alpha, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Mean Trait-Trait Correlation", y = "KDE Richness") +
  theme_pubr()

DD <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.evenness.eblue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.evenness.eblue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.evenness.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.evenness.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = mean_cor, y = kde.evenness, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Mean Trait-Trait Correlation", y = "KDE Evenness") +
  theme_pubr()

HH <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.dispersion.eblue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.dispersion.eblue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.dispersion.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.dispersion.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = mean_cor, y = kde.dispersion, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Mean Trait-Trait Correlation", y = "KDE Dispersion") +
  theme_pubr()

png(here("Figures/seve_meancorr_euc.png"), height = 5, width = 9, units = 'in', res = 300)
ggarrange(plotlist = list(AA, BB, CC, DD, EE, FF, GG, HH), common.legend = TRUE, 
          labels = c("AA", "BB", "CC", "DD", "EE", "FF", "GG", "HH"))
dev.off()


png(here("Figures/mean_corr_euc.png"), height = 8, width = 13, units = 'in', res = 300)
ggarrange(plotlist = list(A, AA, B, BB, C, CC, D, DD, E, EE, F, FF, G, GG, H, HH), 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", 
                     "O", "P"), ncol = 4, nrow = 4, common.legend = TRUE, align = "hv")
dev.off()

