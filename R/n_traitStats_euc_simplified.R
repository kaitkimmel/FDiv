##########################################
#### STATISTICAL ANALYSES  EUCLIDEAN ####
#########################################

# CREATED BY: KAITLIN KIMMEL

# load libraries
library(nlme)
library(here)
library(AICcmodavg)
library(ggplot2)
library(ggpubr)
library(dplyr)


# load data

cdre.1 <- read.csv(here("data/Cleaned/cdr1_euc.csv"))
cdre.2 <- read.csv(here("data/Cleaned/cdr2_euc.csv"))
cdre.3 <- read.csv(here("data/Cleaned/cdr3_euc.csv"))
cdre.4 <- read.csv(here("data/Cleaned/cdr4_euc.csv"))
sev.eblue <- read.csv(here("data/Cleaned/sevblue_euc.csv"))
sev.eblack <- read.csv(here("data/Cleaned/sevblack_euc.csv"))

#############################################################
################ SEV #######################################
###########################################################

###########################################################
################### eblue ################################
##########################################################


###### Fit Metric ~ N_trait models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ n_trait, random = ~1|sev.blueplots,  data = sev.eblue[-which(is.na(sev.eblue$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots, sev.eblue[-which(is.na(sev.eblue$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(FRic ~ 1, random = ~1|sev.blueplots,  data = sev.eblue[-which(is.na(sev.eblue$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) ### mod
fricmod_eblue <- lme(FRic ~ n_trait, random = ~1|sev.blueplots,  data = sev.eblue[-which(is.na(sev.eblue$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE richness
mod <- lme(kde.alpha ~ n_trait, random = ~1|sev.blueplots,  data = sev.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots, sev.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(kde.alpha ~ 1, random = ~1|sev.blueplots,  data = sev.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) ### mod1
kde.alphamod_eblue <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|sev.blueplots,  data = sev.eblue[-which(is.na(sev.eblue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.eblue[-which(is.na(sev.eblue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(FEve ~ 1, random = ~1|sev.blueplots,  data = sev.eblue[-which(is.na(sev.eblue$FEve)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) ## mod1
anova(mod, mod1) #not sig diff - using mod
fevemod_eblue <- lme(FEve ~ n_trait, random = ~1|sev.blueplots,  data = sev.eblue[-which(is.na(sev.eblue$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE evenness
mod <- lme(kde.evenness ~ n_trait, random = ~1|sev.blueplots,  data = sev.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots, sev.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(kde.evenness ~ 1, random = ~1|sev.blueplots,  data = sev.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) ### mod1
kde.evennessmod_eblue <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

## FDis
mod <- lme(FDis ~ n_trait, random = ~1|sev.blueplots,  data = sev.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(FDis ~ 1, random = ~1|sev.blueplots,  data = sev.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) #mod1 best fit
fdismod_eblue <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

##KDE dispersion
mod <- lme(kde.dispersion ~ n_trait, random = ~1|sev.blueplots,  data = sev.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots, sev.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|sev.blueplots,  data = sev.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) ### mod1
anova(mod1, mod) # sig diff - use mod 1
kde.dispersionmod_eblue <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))


# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|sev.blueplots,  data = sev.eblue[-which(is.na(sev.eblue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.eblue[-which(is.na(sev.eblue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(FDiv ~ 1, random = ~1|sev.blueplots,  data = sev.eblue[-which(is.na(sev.eblue$FDiv)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) #mod1 best fit
anova(mod, mod1) #mod1 sig diff from mod
fdivmod_eblue <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.eblue[-which(is.na(sev.eblue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|sev.blueplots,  data = sev.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots),control = lmeControl(opt = "optim"))
mod4 <- lme(RaoQ ~ 1, random = ~1|sev.blueplots,  data = sev.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod4) #mod best fit
raoqmod_eblue <- lme(RaoQ ~ n_trait, random = ~1|sev.blueplots,  data = sev.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))


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


###### Fit Metric ~ N_trait models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ n_trait, random = ~1|sev.blackplots,  data = sev.eblack[-which(is.na(sev.eblack$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots, sev.eblack[-which(is.na(sev.eblack$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(FRic ~ 1, random = ~1|sev.blackplots,  data = sev.eblack[-which(is.na(sev.eblack$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
AIC(mod, mod1, mod4) ### mod, mod1 similar
anova(mod, mod1) #mod for simplicity
fricmod_eblack <- lme(FRic ~ n_trait, random = ~1|sev.blackplots,  data = sev.eblack[-which(is.na(sev.eblack$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE richness
mod <- lme(kde.alpha ~ n_trait, random = ~1|sev.blackplots,  data = sev.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots, sev.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(kde.alpha ~ 1, random = ~1|sev.blackplots,  data = sev.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) ### mod1
anova(mod1, mod) #sig diff use mod 1
kde.alphamod_eblack <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|sev.blackplots,  data = sev.eblack[-which(is.na(sev.eblack$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.eblack[-which(is.na(sev.eblack$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(FEve ~ 1, random = ~1|sev.blackplots,  data = sev.eblack[-which(is.na(sev.eblack$FEve)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) ## mod
fevemod_eblack <- lme(FEve ~ n_trait, random = ~1|sev.blackplots,  data = sev.eblack[-which(is.na(sev.eblack$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE evenness
mod <- lme(kde.evenness ~ n_trait, random = ~1|sev.blackplots,  data = sev.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots, sev.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(kde.evenness ~ 1, random = ~1|sev.blackplots,  data = sev.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) ### mod
kde.evennessmod_eblack <- lme(kde.evenness ~ n_trait, random = ~1|sev.blackplots,  data = sev.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))



## FDis
mod <- lme(FDis ~ n_trait, random = ~1|sev.blackplots,  data = sev.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(FDis ~ 1, random = ~1|sev.blackplots,  data = sev.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) #mod1 best fit
fdismod_eblack <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE dispersion
mod <- lme(kde.dispersion ~ n_trait, random = ~1|sev.blackplots,  data = sev.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots, sev.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|sev.blackplots,  data = sev.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1,mod4) ### mod1
kde.dispersionmod_eblack <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))

# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|sev.blackplots,  data = sev.eblack[-which(is.na(sev.eblack$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.eblack[-which(is.na(sev.eblack$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(FDiv ~ 1, random = ~1|sev.blackplots,  data = sev.eblack[-which(is.na(sev.eblack$FDiv)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) #mod1 best fit
anova(mod, mod1) #mod1 sig diff 
fdivmod_eblack <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.eblack[-which(is.na(sev.eblack$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|sev.blackplots,  data = sev.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots),control = lmeControl(opt = "optim"))
mod4 <- lme(RaoQ ~ 1, random = ~1|sev.blackplots,  data = sev.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod4) #mod
raoqmod_eblack <- lme(RaoQ ~ n_trait, random = ~1|sev.blackplots,  data = sev.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))


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

###### Fit Metric ~ N_trait models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ n_trait, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FRic ~ 1, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
fricmod_cdre1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ n_trait, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.alpha ~ 1, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
kde.alphamod_cdre1 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FEve ~ 1, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod4) ## mod1 best fit
fevemod_cdre1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ n_trait, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.evenness ~ 1, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
kde.evennessmod_cdre1 <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))


## FDis
mod <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDis ~ 1, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod1 best fit
fdismod_cdre1 <- lme(FDis ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

##KDE dispersion
mod <- lme(kde.dispersion ~ n_trait, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
kde.dispersionmod_cdre1 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDiv ~ 1, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod1 best fit
fdivmod_cdre1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(RaoQ ~ 1, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod best fit
raoqmod_cdre1 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))


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
###### Fit Metric ~ N_trait models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ n_trait, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FRic ~ 1, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod4) ### mod1
fricmod_cdre2 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ n_trait, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.alpha ~ 1, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
kde.alphamod_cdre2 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FEve ~ 1, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ## mod best fit
fevemod_cdre2 <- lme(FEve ~ n_trait, random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ n_trait, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.evenness ~ 1, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod
kde.evennessmod_cdre2 <- lme(kde.evenness ~ n_trait, random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDis ~ 1, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod4) #mod1 best fit
fdismod_cdre2 <- lme(FDis ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

##KDE dispersion
mod <- lme(kde.dispersion ~ n_trait, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
anova(mod, mod1) ## sig diff use mod1
kde.dispersionmod_cdre2 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDiv ~ 1, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod1 best fit
fdivmod_cdre2 <- lme(FDiv ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))



### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(RaoQ ~ 1, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod
raoqmod_cdre2 <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))


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
###### Fit Metric ~ N_trait models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ n_trait, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FRic ~ 1, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod4) ### mod1
fricmod_cdre3 <- lme(FRic ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ n_trait, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.alpha ~ 1, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
kde.alphamod_cdre3 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FEve~ 1, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ## mod1 best fit
fevemod_cdre3 <- lme(FEve ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ n_trait, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.evenness ~ 1, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
kde.evennessmod_cdre3 <- lme(kde.evenness ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDis ~ 1, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod1
fdismod_cdre3 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

##KDE dispersion
mod <- lme(kde.dispersion ~ n_trait, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
kde.dispersionmod_cdre3 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDiv ~ 1, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod1 best fit
fdivmod_cdre3 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(RaoQ ~ 1, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod 
raoqmod_cdre3 <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))


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
###### Fit Metric ~ N_trait models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ n_trait, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FRic ~ 1, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod
fricmod_cdre4 <- lme(FRic ~ n_trait, random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ n_trait, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.alpha ~ 1, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
anova(mod, mod1) # sig diff use mod1
kde.alphamod_cdre4 <- lme(kde.alpha ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FEve ~ 1, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ## mod1
fevemod_cdre4 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ n_trait, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.evenness ~ 1, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod4
kde.evennessmod_cdre4 <- lme(kde.evenness ~ 1, random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDis ~ 1, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod1 
fdismod_cdre4 <- lme(FDis ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

##KDE dispersion
mod <- lme(kde.dispersion ~ n_trait, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) ### mod1
anova(mod, mod1) # sig diff use mod1
kde.dispersionmod_cdre4 <- lme(kde.dispersion ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDiv ~ 1, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod4) #mod1 best fit
fdivmod_cdre4 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(RaoQ ~ 1, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod4) #mod 
raoqmod_cdre4 <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod_cdre4)
summary(fevemod_cdre4)
summary(fdismod_cdre4)
summary(fdivmod_cdre4)
summary(raoqmod_cdre4)
summary(kde.alphamod_cdre4)
summary(kde.evennessmod_cdre4)
summary(kde.dispersionmod_cdre4)

################### GRAPHS ##############################

# create combined data set of raw data
cdre.1$community <- "ambient"
cdre.2$community <- "elevated_N"
cdre.3$community <- "elevated_CO2"
cdre.4$community <- "elevated_N_CO2"
cdre.full <- rbind(cdre.1, cdre.2, cdre.3, cdre.4)

cdre.full <- cdre.full %>% group_by(community, n_trait) %>%
  summarize (FRic = mean(FRic), FEve = mean(FEve), FDis = mean(FDis),
             FDiv = mean(FDiv), RaoQ = mean(RaoQ), kde.alpha = mean(kde.alpha), kde.evenness = mean(kde.evenness), kde.dispersion = mean(kde.dispersion))

new.df <- data.frame(n_trait = seq(2, 9, by = 1))


##############################################
############## cdre 1 ########################
#############################################
fr.e1 <- as.data.frame(predictSE.lme(fricmod_cdre1, new.df))
fr.e1$n_trait <- seq(2, 9, by = 1)
fr.e1$lwr <- fr.e1$fit - fr.e1$se.fit
fr.e1$upr <- fr.e1$fit + fr.e1$se.fit
fe.e1 <- as.data.frame(predictSE.lme(fevemod_cdre1, new.df))
fe.e1$n_trait <- seq(2, 9, by = 1)
fe.e1$lwr <- fe.e1$fit - fe.e1$se.fit
fe.e1$upr <- fe.e1$fit + fe.e1$se.fit
fdis.e1 <- as.data.frame(predictSE.lme(fdismod_cdre1, new.df))
fdis.e1$n_trait <- seq(2, 9, by = 1)
fdis.e1$lwr <- fdis.e1$fit - fdis.e1$se.fit
fdis.e1$upr <- fdis.e1$fit + fdis.e1$se.fit
fdiv.e1 <- as.data.frame(predictSE.lme(fdivmod_cdre1, new.df))
fdiv.e1$n_trait <- seq(2, 9, by = 1)
fdiv.e1$lwr <- fdiv.e1$fit - fdiv.e1$se.fit
fdiv.e1$upr <- fdiv.e1$fit + fdiv.e1$se.fit
rq.e1 <- as.data.frame(predictSE.lme(raoqmod_cdre1, new.df))
rq.e1$n_trait <- seq(2, 9, by = 1)
rq.e1$lwr <- rq.e1$fit - rq.e1$se.fit
rq.e1$upr <- rq.e1$fit + rq.e1$se.fit
kde.alpha.e1 <- as.data.frame(predictSE.lme(kde.alphamod_cdre1, new.df))
kde.alpha.e1$n_trait <- seq(2, 9, by = 1)
kde.alpha.e1$lwr <- kde.alpha.e1$fit - kde.alpha.e1$se.fit
kde.alpha.e1$upr <- kde.alpha.e1$fit + kde.alpha.e1$se.fit
kde.evenness.e1 <- as.data.frame(predictSE.lme(kde.evennessmod_cdre1, new.df))
kde.evenness.e1$n_trait <- seq(2, 9, by = 1)
kde.evenness.e1$lwr <- kde.evenness.e1$fit - kde.evenness.e1$se.fit
kde.evenness.e1$upr <- kde.evenness.e1$fit + kde.evenness.e1$se.fit
kde.dispersion.e1 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdre1, new.df))
kde.dispersion.e1$n_trait <- seq(2, 9, by = 1)
kde.dispersion.e1$lwr <- kde.dispersion.e1$fit - kde.dispersion.e1$se.fit
kde.dispersion.e1$upr <- kde.dispersion.e1$fit + kde.dispersion.e1$se.fit

fr.e2 <- as.data.frame(predictSE.lme(fricmod_cdre2, new.df))
fr.e2$n_trait <- seq(2, 9, by = 1)
fr.e2$lwr <- fr.e2$fit - fr.e2$se.fit
fr.e2$upr <- fr.e2$fit + fr.e2$se.fit
fe.e2 <- as.data.frame(predictSE.lme(fevemod_cdre2, new.df))
fe.e2$n_trait <- seq(2, 9, by = 1)
fe.e2$lwr <- fe.e2$fit - fe.e2$se.fit
fe.e2$upr <- fe.e2$fit + fe.e2$se.fit
fdis.e2 <- as.data.frame(predictSE.lme(fdismod_cdre2, new.df))
fdis.e2$n_trait <- seq(2, 9, by = 1)
fdis.e2$lwr <- fdis.e2$fit - fdis.e2$se.fit
fdis.e2$upr <- fdis.e2$fit + fdis.e2$se.fit
fdiv.e2 <- as.data.frame(predictSE.lme(fdivmod_cdre2, new.df))
fdiv.e2$n_trait <- seq(2, 9, by = 1)
fdiv.e2$lwr <- fdiv.e2$fit - fdiv.e2$se.fit
fdiv.e2$upr <- fdiv.e2$fit + fdiv.e2$se.fit
rq.e2 <- as.data.frame(predictSE.lme(raoqmod_cdre2, new.df))
rq.e2$n_trait <- seq(2, 9, by = 1)
rq.e2$lwr <- rq.e2$fit - rq.e2$se.fit
rq.e2$upr <- rq.e2$fit + rq.e2$se.fit
kde.alpha.e2 <- as.data.frame(predictSE.lme(kde.alphamod_cdre2, new.df))
kde.alpha.e2$n_trait <- seq(2, 9, by = 1)
kde.alpha.e2$lwr <- kde.alpha.e2$fit - kde.alpha.e2$se.fit
kde.alpha.e2$upr <- kde.alpha.e2$fit + kde.alpha.e2$se.fit
kde.evenness.2 <- as.data.frame(predictSE.lme(kde.evennessmod_cdre2, new.df))
kde.evenness.2$n_trait <- seq(2, 9, by = 1)
kde.evenness.2$lwr <- kde.evenness.2$fit - kde.evenness.2$se.fit
kde.evenness.2$upr <- kde.evenness.2$fit + kde.evenness.2$se.fit
kde.dispersion.e2 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdre2, new.df))
kde.dispersion.e2$n_trait <- seq(2, 9, by = 1)
kde.dispersion.e2$lwr <- kde.dispersion.e2$fit - kde.dispersion.e2$se.fit
kde.dispersion.e2$upr <- kde.dispersion.e2$fit + kde.dispersion.e2$se.fit

fr.e3 <- as.data.frame(predictSE.lme(fricmod_cdre3, new.df))
fr.e3$n_trait <- seq(2, 9, by = 1)
fr.e3$lwr <- fr.e3$fit - fr.e3$se.fit
fr.e3$upr <- fr.e3$fit + fr.e3$se.fit
fe.e3 <- as.data.frame(predictSE.lme(fevemod_cdre3, new.df))
fe.e3$n_trait <- seq(2, 9, by = 1)
fe.e3$lwr <- fe.e3$fit - fe.e3$se.fit
fe.e3$upr <- fe.e3$fit + fe.e3$se.fit
fdis.e3 <- as.data.frame(predictSE.lme(fdismod_cdre3, new.df))
fdis.e3$n_trait <- seq(2, 9, by = 1)
fdis.e3$lwr <- fdis.e3$fit - fdis.e3$se.fit
fdis.e3$upr <- fdis.e3$fit + fdis.e3$se.fit
fdiv.e3 <- as.data.frame(predictSE.lme(fdivmod_cdre3, new.df))
fdiv.e3$n_trait <- seq(2, 9, by = 1)
fdiv.e3$lwr <- fdiv.e3$fit - fdiv.e3$se.fit
fdiv.e3$upr <- fdiv.e3$fit + fdiv.e3$se.fit
rq.e3 <- as.data.frame(predictSE.lme(raoqmod_cdre3, new.df))
rq.e3$n_trait <- seq(2, 9, by = 1)
rq.e3$lwr <- rq.e3$fit - rq.e3$se.fit
rq.e3$upr <- rq.e3$fit + rq.e3$se.fit
kde.alpha.e3 <- as.data.frame(predictSE.lme(kde.alphamod_cdre3, new.df))
kde.alpha.e3$n_trait <- seq(2, 9, by = 1)
kde.alpha.e3$lwr <- kde.alpha.e3$fit - kde.alpha.e3$se.fit
kde.alpha.e3$upr <- kde.alpha.e3$fit + kde.alpha.e3$se.fit
kde.evenness.3 <- as.data.frame(predictSE.lme(kde.evennessmod_cdre3, new.df))
kde.evenness.3$n_trait <- seq(2, 9, by = 1)
kde.evenness.3$lwr <- kde.evenness.3$fit - kde.evenness.3$se.fit
kde.evenness.3$upr <- kde.evenness.3$fit + kde.evenness.3$se.fit
kde.dispersion.e3 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdre3, new.df))
kde.dispersion.e3$n_trait <- seq(2, 9, by = 1)
kde.dispersion.e3$lwr <- kde.dispersion.e3$fit - kde.dispersion.e3$se.fit
kde.dispersion.e3$upr <- kde.dispersion.e3$fit + kde.dispersion.e3$se.fit

fr.e4 <- as.data.frame(predictSE.lme(fricmod_cdre4, new.df))
fr.e4$n_trait <- seq(2, 9, by = 1)
fr.e4$lwr <- fr.e4$fit - fr.e4$se.fit
fr.e4$upr <- fr.e4$fit + fr.e4$se.fit
fe.e4 <- as.data.frame(predictSE.lme(fevemod_cdre4, new.df))
fe.e4$n_trait <- seq(2, 9, by = 1)
fe.e4$lwr <- fe.e4$fit - fe.e4$se.fit
fe.e4$upr <- fe.e4$fit + fe.e4$se.fit
fdis.e4 <- as.data.frame(predictSE.lme(fdismod_cdre4, new.df))
fdis.e4$n_trait <- seq(2, 9, by = 1)
fdis.e4$lwr <- fdis.e4$fit - fdis.e4$se.fit
fdis.e4$upr <- fdis.e4$fit + fdis.e4$se.fit
fdiv.e4 <- as.data.frame(predictSE.lme(fdivmod_cdre4, new.df))
fdiv.e4$n_trait <- seq(2, 9, by = 1)
fdiv.e4$lwr <- fdiv.e4$fit - fdiv.e4$se.fit
fdiv.e4$upr <- fdiv.e4$fit + fdiv.e4$se.fit
rq.e4 <- as.data.frame(predictSE.lme(raoqmod_cdre4, new.df))
rq.e4$n_trait <- seq(2, 9, by = 1)
rq.e4$lwr <- rq.e4$fit - rq.e4$se.fit
rq.e4$upr <- rq.e4$fit + rq.e4$se.fit
kde.alpha.e4 <- as.data.frame(predictSE.lme(kde.alphamod_cdre4, new.df))
kde.alpha.e4$n_trait <- seq(2, 9, by = 1)
kde.alpha.e4$lwr <- kde.alpha.e4$fit - kde.alpha.e4$se.fit
kde.alpha.e4$upr <- kde.alpha.e4$fit + kde.alpha.e4$se.fit
kde.evenness.e4 <- as.data.frame(predictSE.lme(kde.evennessmod_cdre4, new.df))
kde.evenness.e4$n_trait <- seq(2, 9, by = 1)
kde.evenness.e4$lwr <- kde.evenness.e4$fit - kde.evenness.e4$se.fit
kde.evenness.e4$upr <- kde.evenness.e4$fit + kde.evenness.e4$se.fit
kde.dispersion.e4 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdre4, new.df))
kde.dispersion.e4$n_trait <- seq(2, 9, by = 1)
kde.dispersion.e4$lwr <- kde.dispersion.e4$fit - kde.dispersion.e4$se.fit
kde.dispersion.e4$upr <- kde.dispersion.e4$fit + kde.dispersion.e4$se.fit


A <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FRic, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "FRich") +
  theme_pubr()

C <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FEve, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "FEve") +
  theme_pubr()

E <-ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FDis, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "FDis") +
  theme_pubr()

F <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FDiv, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "FDiv") +
  theme_pubr()

G <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = RaoQ, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "Rao Q") +
  theme_pubr()

B <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.alpha, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Richness") +
  theme_pubr()

D <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.evenness, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Evenness") +
  theme_pubr()

H <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.dispersion, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Dispersion") +
  theme_pubr()

png(here("Figures/cdre_Ntraits_euc.png"))
ggarrange(plotlist = list(A, B, C, D, E, F, G, H), common.legend = TRUE,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
dev.off()

###SEV####


# create combined data set of raw data
sev.eblack$community <- "black"
sev.eblue$community <- "blue"
seve.full <- rbind(sev.eblue[,-10], sev.eblack[,-10])

seve.full <- seve.full %>% group_by(community, n_trait) %>%
  summarize (FRic = mean(FRic, na.rm = TRUE), FEve = mean(FEve, na.rm = TRUE), 
             FDis = mean(FDis, na.rm = TRUE), FDiv = mean(FDiv, na.rm = TRUE), 
             RaoQ = mean(RaoQ, na.rm = TRUE), kde.alpha = mean(kde.alpha), kde.evenness = mean(kde.evenness), kde.dispersion = mean(kde.dispersion))

new.df <- data.frame(n_trait = seq(2, 10, by = 1))

fr.eblue <- as.data.frame(predictSE.lme(fricmod_eblue, new.df))
fr.eblue$n_trait <- seq(2, 10, by = 1)
fr.eblue$lwr <- fr.eblue$fit - fr.eblue$se.fit
fr.eblue$upr <- fr.eblue$fit + fr.eblue$se.fit
fe.eblue <- as.data.frame(predictSE.lme(fevemod_eblue, new.df))
fe.eblue$n_trait <- seq(2, 10, by = 1)
fe.eblue$lwr <- fe.eblue$fit - fe.eblue$se.fit
fe.eblue$upr <- fe.eblue$fit + fe.eblue$se.fit
fdis.eblue <- as.data.frame(predictSE.lme(fdismod_eblue, new.df))
fdis.eblue$n_trait <- seq(2, 10, by = 1)
fdis.eblue$lwr <- fdis.eblue$fit - fdis.eblue$se.fit
fdis.eblue$upr <- fdis.eblue$fit + fdis.eblue$se.fit
fdiv.eblue <- as.data.frame(predictSE.lme(fdivmod_eblue, new.df))
fdiv.eblue$n_trait <- seq(2, 10, by = 1)
fdiv.eblue$lwr <- fdiv.eblue$fit - fdiv.eblue$se.fit
fdiv.eblue$upr <- fdiv.eblue$fit + fdiv.eblue$se.fit
rq.eblue <- as.data.frame(predictSE.lme(raoqmod_eblue, new.df))
rq.eblue$n_trait <- seq(2, 10, by = 1)
rq.eblue$lwr <- rq.eblue$fit - rq.eblue$se.fit
rq.eblue$upr <- rq.eblue$fit + rq.eblue$se.fit
kde.alpha.eblue <- as.data.frame(predictSE.lme(kde.alphamod_eblue, new.df))
kde.alpha.eblue$n_trait <- seq(2, 10, by = 1)
kde.alpha.eblue$lwr <- kde.alpha.eblue$fit - kde.alpha.eblue$se.fit
kde.alpha.eblue$upr <- kde.alpha.eblue$fit + kde.alpha.eblue$se.fit
kde.evenness.eblue <- as.data.frame(predictSE.lme(kde.evennessmod_eblue, new.df))
kde.evenness.eblue$n_trait <- seq(2, 10, by = 1)
kde.evenness.eblue$lwr <- kde.evenness.eblue$fit - kde.evenness.eblue$se.fit
kde.evenness.eblue$upr <- kde.evenness.eblue$fit + kde.evenness.eblue$se.fit
kde.dispersion.eblue <- as.data.frame(predictSE.lme(kde.dispersionmod_eblue, new.df))
kde.dispersion.eblue$n_trait <- seq(2, 10, by = 1)
kde.dispersion.eblue$lwr <- kde.dispersion.eblue$fit - kde.dispersion.eblue$se.fit
kde.dispersion.eblue$upr <- kde.dispersion.eblue$fit + kde.dispersion.eblue$se.fit

fr.eblack <- as.data.frame(predictSE.lme(fricmod_eblack, new.df))
fr.eblack$n_trait <- seq(2, 10, by = 1)
fr.eblack$lwr <- fr.eblack$fit - fr.eblack$se.fit
fr.eblack$upr <- fr.eblack$fit + fr.eblack$se.fit
fe.eblack <- as.data.frame(predictSE.lme(fevemod_eblack, new.df))
fe.eblack$n_trait <- seq(2, 10, by = 1)
fe.eblack$lwr <- fe.eblack$fit - fe.eblack$se.fit
fe.eblack$upr <- fe.eblack$fit + fe.eblack$se.fit
fdis.eblack <- as.data.frame(predictSE.lme(fdismod_eblack, new.df))
fdis.eblack$n_trait <- seq(2, 10, by = 1)
fdis.eblack$lwr <- fdis.eblack$fit - fdis.eblack$se.fit
fdis.eblack$upr <- fdis.eblack$fit + fdis.eblack$se.fit
fdiv.eblack <- as.data.frame(predictSE.lme(fdivmod_eblack, new.df))
fdiv.eblack$n_trait <- seq(2, 10, by = 1)
fdiv.eblack$lwr <- fdiv.eblack$fit - fdiv.eblack$se.fit
fdiv.eblack$upr <- fdiv.eblack$fit + fdiv.eblack$se.fit
rq.eblack <- as.data.frame(predictSE.lme(raoqmod_eblack, new.df))
rq.eblack$n_trait <- seq(2, 10, by = 1)
rq.eblack$lwr <- rq.eblack$fit - rq.eblack$se.fit
rq.eblack$upr <- rq.eblack$fit + rq.eblack$se.fit
kde.alpha.eblack <- as.data.frame(predictSE.lme(kde.alphamod_eblack, new.df))
kde.alpha.eblack$n_trait <- seq(2, 10, by = 1)
kde.alpha.eblack$lwr <- kde.alpha.eblack$fit - kde.alpha.eblack$se.fit
kde.alpha.eblack$upr <- kde.alpha.eblack$fit + kde.alpha.eblack$se.fit
kde.evenness.eblack <- as.data.frame(predictSE.lme(kde.evennessmod_eblack, new.df))
kde.evenness.eblack$n_trait <- seq(2, 10, by = 1)
kde.evenness.eblack$lwr <- kde.evenness.eblack$fit - kde.evenness.eblack$se.fit
kde.evenness.eblack$upr <- kde.evenness.eblack$fit + kde.evenness.eblack$se.fit
kde.dispersion.eblack <- as.data.frame(predictSE.lme(kde.dispersionmod_eblack, new.df))
kde.dispersion.eblack$n_trait <- seq(2, 10, by = 1)
kde.dispersion.eblack$lwr <- kde.dispersion.eblack$fit - kde.dispersion.eblack$se.fit
kde.dispersion.eblack$upr <- kde.dispersion.eblack$fit + kde.dispersion.eblack$se.fit
AA <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.eblue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.eblue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FRic, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FRich") +
  theme_pubr()

CC <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.eblue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.eblue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FEve, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FEve") +
  theme_pubr()

EE <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.eblue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.eblue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FDis, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FDis") +
  theme_pubr()

FF <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.eblue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.eblue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FDiv, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FDiv") +
  theme_pubr()

GG <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.eblue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.eblue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = RaoQ, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "Rao Q") +
  theme_pubr()

BB <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.eblue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.eblue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.alpha, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Richness") +
  theme_pubr()

DD <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.eblue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.eblue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.evenness, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Eveness") +
  theme_pubr()

HH <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.eblue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.eblue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.dispersion, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Dispersion") +
  theme_pubr()

png(here("Figures/Sev_Ntraits_euc.png"))
ggarrange(plotlist = list(AA, BB, CC, DD, EE, FF, GG, HH), common.legend = TRUE, 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
dev.off()



png(here("Figures/n_trait_euc.png"), height = 8, width = 13, units = 'in', res = 300)
ggarrange(plotlist = list(A, AA, B, BB, C, CC, D, DD, E, EE, F, FF, G, GG, H, HH), 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", 
                     "O", "P"), ncol = 4, nrow = 4, common.legend = TRUE, align = "hv")
dev.off()

