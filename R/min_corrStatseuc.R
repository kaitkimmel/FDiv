################################################################
#### STATISTICAL ANALYSES  EUCLIDEAN for  MIN TRAIT CORRELATIONS ####
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

cdre.1 <- read.csv(here("data/Cleaned/cdr1_euc.csv"))
cdre.1 <- cdre.1[-which(is.na(cdre.1$min_cor)),]
cdre.2 <- read.csv(here("data/Cleaned/cdr2_euc.csv"))
cdre.2 <- cdre.2[-which(is.na(cdre.2$min_cor)),]
cdre.3 <- read.csv(here("data/Cleaned/cdr3_euc.csv"))
cdre.3 <- cdre.3[-which(is.na(cdre.3$min_cor)),]
cdre.4 <- read.csv(here("data/Cleaned/cdr4_euc.csv"))
cdre.4 <- cdre.4[-which(is.na(cdre.4$min_cor)),]
seve.eblue <- read.csv(here("data/Cleaned/sevblue_euc.csv"))
seve.eblue <- seve.eblue[-which(is.na(seve.eblue$min_cor)),]
seve.eblack <- read.csv(here("data/Cleaned/sevblack_euc.csv"))
seve.eblack <- seve.eblack[-which(is.na(seve.eblack$min_cor)),]



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
mod <- lme(FRic ~ min_cor, random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FRic ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots, seve.eblue[-which(is.na(seve.eblue$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FRic ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FRic ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(FRic ~ 1, random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
AIC(mod, mod1, mod2, mod3, mod4) ### mod3
fricmod_eblue <- lme(FRic ~ min_cor + I(min_cor^2)+ I(min_cor^3) + I(min_cor^4), random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE richness
mod <- lme(kde.alpha ~ min_cor, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.alpha ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots, seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(kde.alpha ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(kde.alpha ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(kde.alpha ~ 1, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
AIC(mod, mod1, mod2, mod3, mod4) ### mod3
kde.alphamod_eblue <- lme(kde.alpha ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blueplots,  data = seve.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

### FEve
mod <- lme(FEve ~ min_cor, random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FEve ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FEve ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FEve ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(FEve ~ 1, random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FEve)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
AIC(mod, mod1,mod2, mod3, mod4) ## mod3
anova(mod3, mod4) # no sig diff - mod4
fevemod_eblue <- lme(FEve ~ 1, random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE evenness
mod <- lme(kde.evenness ~ min_cor, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.evenness ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots, seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(kde.evenness ~ 1, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
AIC(mod, mod1, mod2, mod3, mod4) ### mod3
kde.evennessmod_eblue <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blueplots,  data = seve.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

## FDis
mod <- lme(FDis ~ min_cor, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FDis ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FDis ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FDis ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(FDis ~ 1, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
AIC(mod, mod1, mod2, mod3, mod4) #mod3 best fit
fdismod_eblue <- lme(FDis ~ min_cor+ I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blueplots,  data = seve.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE dispersion
mod <- lme(kde.dispersion ~ min_cor, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.dispersion ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots, seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(kde.dispersion ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(kde.dispersion ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
AIC(mod, mod1, mod2, mod3, mod4) ### mod3
kde.dispersionmod_eblue <- lme(kde.dispersion ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blueplots,  data = seve.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

# FDiv
mod <- lme(FDiv ~ min_cor, random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FDiv ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FDiv ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FDiv ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(FDiv ~ 1, random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FDiv)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
AIC(mod, mod1, mod2, mod3, mod4) #mod2
fdivmod_eblue <- lme(FDiv ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blueplots,  data = seve.eblue[-which(is.na(seve.eblue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

### Rao's Q
mod <- lme(RaoQ ~ min_cor, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(RaoQ ~ min_cor + I(min_cor^2), random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots),control = lmeControl(opt = "optim"))
mod2 <- lme(RaoQ ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(RaoQ ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blueplots,  data = seve.eblue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod4 <- lme(RaoQ ~ 1, random = ~1|sev.blueplots,  data = seve.eblue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
AIC(mod, mod1, mod2, mod3, mod4) #mod3 best fit
raoqmod_eblue <- lme(RaoQ ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blueplots,  data = seve.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))


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


###### Fit Metric ~ min_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ min_cor, random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FRic ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots, seve.eblack[-which(is.na(seve.eblack$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FRic ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FRic ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(FRic ~ 1, random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
AIC(mod, mod1, mod2, mod3, mod4) ### mod4 lowest
fricmod_eblack <- lme(FRic ~ 1, random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE richness
mod <- lme(kde.alpha ~ min_cor, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.alpha ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots, seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(kde.alpha ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(kde.alpha ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(kde.alpha ~ 1, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod2, mod3, mod4) ### mod3 lowest
kde.alphamod_eblack <- lme(kde.alpha ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blackplots,  data = seve.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))

### FEve
mod <- lme(FEve ~ min_cor, random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FEve ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FEve ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FEve ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(FEve ~ 1, random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FEve)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod2, mod3, mod4) ## mod4 and mod1 similar - use mod 4
fevemod_eblack <- lme(FEve ~ 1, random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE evenness
mod <- lme(kde.evenness ~ min_cor, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.evenness ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots, seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(kde.evenness ~ 1, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod2, mod3, mod4) ### mod2 lowest
kde.evennessmod_eblack <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blackplots,  data = seve.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))

## FDis
mod <- lme(FDis ~ min_cor, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FDis ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FDis ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FDis ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(FDis ~ 1, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod2, mod3, mod4) #mod3 best fit
fdismod_eblack <- lme(FDis ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blackplots,  data = seve.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE dispersion
mod <- lme(kde.dispersion ~ min_cor, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.dispersion ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots, seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(kde.dispersion ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(kde.dispersion ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod2, mod3, mod4) ### mod3 lowest
kde.dispersionmod_eblack <- lme(kde.dispersion ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blackplots,  data = seve.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))

# FDiv
mod <- lme(FDiv ~ min_cor, random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FDiv ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FDiv ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FDiv ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(FDiv ~ 1, random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FDiv)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod2, mod3, mod4) #mod3 best fit
fdivmod_eblack <- lme(FDiv ~ min_cor+ I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blackplots,  data = seve.eblack[-which(is.na(seve.eblack$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

### Rao's Q
mod <- lme(RaoQ ~ min_cor, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(RaoQ ~ min_cor + I(min_cor^2), random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots),control = lmeControl(opt = "optim"))
mod2 <- lme(RaoQ ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(RaoQ ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blackplots,  data = seve.eblack, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod4 <- lme(RaoQ ~ 1, random = ~1|sev.blackplots,  data = seve.eblack, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod2, mod3, mod4) #mod3
raoqmod_eblack <- lme(RaoQ ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|sev.blackplots,  data = seve.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))


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

###### Fit Metric ~ min_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ min_cor, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FRic ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FRic ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FRic ~ 1, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
AIC(mod, mod1, mod2, mod3, mod4) ### mod23
fricmod_cdre1 <- lme(FRic ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ min_cor, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.alpha ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.alpha ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.alpha ~ 1, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ### mod3
kde.alphamod_cdre1 <- lme(kde.alpha ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ min_cor, random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FEve ~ 1, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ## mod2
fevemod_cdre1 <- lme(FEve ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ min_cor, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.evenness ~ 1, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ### mod3
kde.evennessmod_cdre1 <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ min_cor, random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDis ~ 1, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) #mod3
fdismod_cdre1 <- lme(FDis ~ min_cor+ I(min_cor^2)+ I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

##kde dspersion
mod <- lme(kde.dispersion ~ min_cor, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.dispersion ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.dispersion ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ### mod3
kde.dispersionmod_cdre1 <- lme(kde.dispersion ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ min_cor, random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDiv ~ 1, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) #mod3
fdivmod_cdre1 <- lme(FDiv ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ min_cor, random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(RaoQ ~ 1, random = ~1|Plot,  data = cdre.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) #mod3
raoqmod_cdre1 <- lme(RaoQ ~ min_cor + I(min_cor^2)+ I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))


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
###### Fit Metric ~ min_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ min_cor, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FRic ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FRic ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FRic ~ 1, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ### mod3
fricmod_cdre2 <- lme(FRic ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ min_cor, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.alpha ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.alpha ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.alpha ~ 1, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ### mod3
kde.alphamod_cdre2 <- lme(kde.alpha ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ min_cor, random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FEve ~ 1, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ## mod1
fevemod_cdre2 <- lme(FEve ~ min_cor+ I(min_cor^2), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ min_cor, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.evenness ~ 1, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ### mod3
kde.evennessmod_cdre2 <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ min_cor, random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDis ~ 1, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) #mod2
anova(mod2, mod1)
fdismod_cdre2 <- lme(FDis ~ min_cor+ I(min_cor^2)+ I(min_cor^3), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

##KDE disperion
mod <- lme(kde.dispersion ~ min_cor, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.dispersion ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.dispersion ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ### mod2
kde.dispersionmod_cdre2 <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ min_cor, random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDiv ~ 1, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) #mod3 best fit
fdivmod_cdre2 <- lme(FDiv ~ min_cor+ I(min_cor^2)+ I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ min_cor, random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(RaoQ ~ 1, random = ~1|Plot,  data = cdre.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) #mod3
raoqmod_cdre2 <- lme(RaoQ ~ min_cor+ I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))


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
###### Fit Metric ~ min_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ min_cor, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FRic ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FRic ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FRic ~ 1, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ### mod3
fricmod_cdre3 <- lme(FRic ~ min_cor+ I(min_cor^2)+ I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ min_cor, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.alpha ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.alpha ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.alpha ~ 1, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ### mod3
anova(mod3, mod2) #sig dif use mod3
kde.alphamod_cdre3 <- lme(kde.alpha ~ min_cor+ I(min_cor^2)+ I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ min_cor, random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FEve ~ 1, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ## mod2
fevemod_cdre3 <- lme(FEve ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ min_cor, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.evenness ~ 1, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ### mod3
kde.evennessmod_cdre3 <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ min_cor, random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDis ~ 1, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) #mod3 best fit
fdismod_cdre3 <- lme(FDis ~ min_cor + I(min_cor^2) + I(min_cor^3)+ I(min_cor^4), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

##KDE dispersion
mod <- lme(kde.dispersion ~ min_cor, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.dispersion ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.dispersion ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ### mod2
kde.dispersionmod_cdre3 <- lme(kde.dispersion ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ min_cor, random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDiv ~ 1, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) #mod3 best fit
fdivmod_cdre3 <- lme(FDiv ~ min_cor + I(min_cor^2) + I(min_cor^3)+ I(min_cor^4), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ min_cor, random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(RaoQ ~ 1, random = ~1|Plot,  data = cdre.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) #mod3
raoqmod_cdre3 <- lme(RaoQ ~ min_cor+ I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))


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
###### Fit Metric ~ min_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ min_cor, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FRic ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FRic ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FRic ~ 1, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ### mod3
fricmod_cdre4 <- lme(FRic ~ min_cor+ I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ min_cor, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.alpha ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.alpha ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.alpha ~ 1, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ### mod2
anova(mod2, mod1)
kde.alphamod_cdre4 <- lme(kde.alpha ~ min_cor+ I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ min_cor, random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FEve ~ 1, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ## mod2 similar to mod and mod3
anova(mod2, mod)# sig diff - use mod2
fevemod_cdre4 <- lme(FEve ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ min_cor, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.evenness ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.evenness ~ 1, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ### mod3
anova(mod3, mod2) ##no sig dif use mod2
kde.evennessmod_cdre4 <- lme(kde.evenness ~ min_cor+ I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ min_cor, random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDis ~ 1, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) #mod3
fdismod_cdre4 <- lme(FDis ~ min_cor+ I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

##KDE dispersion
mod <- lme(kde.dispersion ~ min_cor, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.dispersion ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.dispersion ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(kde.dispersion ~ 1, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) ### mod2
kde.dispersionmod_cdre4 <- lme(kde.dispersion ~ min_cor+ I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ min_cor, random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(FDiv ~ 1, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) #mod2
anova(mod2, mod1)
fdivmod_cdre4 <- lme(FDiv ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ min_cor, random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ min_cor + I(min_cor^2), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ min_cor + I(min_cor^2) + I(min_cor^3), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ min_cor + I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod4 <- lme(RaoQ ~ 1, random = ~1|Plot,  data = cdre.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3, mod4) #mod3
raoqmod_cdre4 <- lme(RaoQ ~ min_cor+ I(min_cor^2) + I(min_cor^3) + I(min_cor^4), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))


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

cdre.full <- cdre.full %>% group_by(community, min_cor) %>%
  summarize (FRic = mean(FRic), FEve = mean(FEve), FDis = mean(FDis),
             FDiv = mean(FDiv), RaoQ = mean(RaoQ), RaoQ = mean(RaoQ), kde.alpha = mean(kde.alpha), kde.evenness  = mean(kde.evenness), kde.dispersion = mean(kde.dispersion))


##############################################
############## cdre 1 ########################
#############################################
range(cdre.1$min_cor)
new.df7 <- data.frame(min_cor = seq(min(cdre.1$min_cor), max(cdre.1$min_cor), by = 0.01))
fr.e1 <- as.data.frame(predictSE.lme(fricmod_cdre1, new.df7))
fr.e1$min_cor <- seq(min(cdre.1$min_cor), max(cdre.1$min_cor), by = 0.01)
fr.e1$lwr <- fr.e1$fit - fr.e1$se.fit
fr.e1$upr <- fr.e1$fit + fr.e1$se.fit
fe.e1 <- as.data.frame(predictSE.lme(fevemod_cdre1, new.df7))
fe.e1$min_cor <- seq(min(cdre.1$min_cor), max(cdre.1$min_cor), by = 0.01)
fe.e1$lwr <- fe.e1$fit - fe.e1$se.fit
fe.e1$upr <- fe.e1$fit + fe.e1$se.fit
fdis.e1 <- as.data.frame(predictSE.lme(fdismod_cdre1, new.df7))
fdis.e1$min_cor <- seq(min(cdre.1$min_cor), max(cdre.1$min_cor), by = 0.01)
fdis.e1$lwr <- fdis.e1$fit - fdis.e1$se.fit
fdis.e1$upr <- fdis.e1$fit + fdis.e1$se.fit
fdiv.e1 <- as.data.frame(predictSE.lme(fdivmod_cdre1, new.df7))
fdiv.e1$min_cor <- seq(min(cdre.1$min_cor), max(cdre.1$min_cor), by = 0.01)
fdiv.e1$lwr <- fdiv.e1$fit - fdiv.e1$se.fit
fdiv.e1$upr <- fdiv.e1$fit + fdiv.e1$se.fit
rq.e1 <- as.data.frame(predictSE.lme(raoqmod_cdre1, new.df7))
rq.e1$min_cor <- seq(min(cdre.1$min_cor), max(cdre.1$min_cor), by = 0.01)
rq.e1$lwr <- rq.e1$fit - rq.e1$se.fit
rq.e1$upr <- rq.e1$fit + rq.e1$se.fit
kde.alpha.e1 <- as.data.frame(predictSE.lme(kde.alphamod_cdre1, new.df7))
kde.alpha.e1$min_cor <- seq(min(cdre.1$min_cor), max(cdre.1$min_cor), by = .01)
kde.alpha.e1$lwr <- kde.alpha.e1$fit - kde.alpha.e1$se.fit
kde.alpha.e1$upr <- kde.alpha.e1$fit + kde.alpha.e1$se.fit
kde.evenness.e1 <- as.data.frame(predictSE.lme(kde.evennessmod_cdre1, new.df7))
kde.evenness.e1$min_cor <- seq(min(cdre.1$min_cor), max(cdre.1$min_cor), by = .01)
kde.evenness.e1$lwr <- kde.evenness.e1$fit - kde.evenness.e1$se.fit
kde.evenness.e1$upr <- kde.evenness.e1$fit + kde.evenness.e1$se.fit
kde.dispersion.e1 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdre1, new.df7))
kde.dispersion.e1$min_cor <- seq(min(cdre.1$min_cor), max(cdre.1$min_cor), by = .01)
kde.dispersion.e1$lwr <- kde.dispersion.e1$fit - kde.dispersion.e1$se.fit
kde.dispersion.e1$upr <- kde.dispersion.e1$fit + kde.dispersion.e1$se.fit

range(cdre.2$min_cor)
new.df8 <- data.frame(min_cor = seq(min(cdre.2$min_cor), max(cdre.2$min_cor), by = 0.01))
fr.e2 <- as.data.frame(predictSE.lme(fricmod_cdre2, new.df8))
fr.e2$min_cor <- seq(min(cdre.2$min_cor), max(cdre.2$min_cor), by = .01)
fr.e2$lwr <- fr.e2$fit - fr.e2$se.fit
fr.e2$upr <- fr.e2$fit + fr.e2$se.fit
fe.e2 <- as.data.frame(predictSE.lme(fevemod_cdre2, new.df8))
fe.e2$min_cor <- seq(min(cdre.2$min_cor), max(cdre.2$min_cor), by = .01)
fe.e2$lwr <- fe.e2$fit - fe.e2$se.fit
fe.e2$upr <- fe.e2$fit + fe.e2$se.fit
fdis.e2 <- as.data.frame(predictSE.lme(fdismod_cdre2, new.df8))
fdis.e2$min_cor <- seq(min(cdre.2$min_cor), max(cdre.2$min_cor), by = .01)
fdis.e2$lwr <- fdis.e2$fit - fdis.e2$se.fit
fdis.e2$upr <- fdis.e2$fit + fdis.e2$se.fit
fdiv.e2 <- as.data.frame(predictSE.lme(fdivmod_cdre2, new.df8))
fdiv.e2$min_cor <- seq(min(cdre.2$min_cor), max(cdre.2$min_cor), by = .01)
fdiv.e2$lwr <- fdiv.e2$fit - fdiv.e2$se.fit
fdiv.e2$upr <- fdiv.e2$fit + fdiv.e2$se.fit
rq.e2 <- as.data.frame(predictSE.lme(raoqmod_cdre2, new.df8))
rq.e2$min_cor <- seq(min(cdre.2$min_cor), max(cdre.2$min_cor), by = .01)
rq.e2$lwr <- rq.e2$fit - rq.e2$se.fit
rq.e2$upr <- rq.e2$fit + rq.e2$se.fit
kde.alpha.e2 <- as.data.frame(predictSE.lme(kde.alphamod_cdre2, new.df8))
kde.alpha.e2$min_cor <- seq(min(cdre.2$min_cor), max(cdre.2$min_cor), by = .01)
kde.alpha.e2$lwr <- kde.alpha.e2$fit - kde.alpha.e2$se.fit
kde.alpha.e2$upr <- kde.alpha.e2$fit + kde.alpha.e2$se.fit
kde.evenness.e2 <- as.data.frame(predictSE.lme(kde.evennessmod_cdre2, new.df8))
kde.evenness.e2$min_cor <- seq(min(cdre.2$min_cor), max(cdre.2$min_cor), by = .01)
kde.evenness.e2$lwr <- kde.evenness.e2$fit - kde.evenness.e2$se.fit
kde.evenness.e2$upr <- kde.evenness.e2$fit + kde.evenness.e2$se.fit
kde.dispersion.e2 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdre2, new.df8))
kde.dispersion.e2$min_cor <- seq(min(cdre.2$min_cor), max(cdre.2$min_cor), by = .01)
kde.dispersion.e2$lwr <- kde.dispersion.e2$fit - kde.dispersion.e2$se.fit
kde.dispersion.e2$upr <- kde.dispersion.e2$fit + kde.dispersion.e2$se.fit

range(cdre.3$min_cor)
new.df9 <- data.frame(min_cor = seq(min(cdre.3$min_cor), max(cdre.3$min_cor), by = 0.01))
fr.e3 <- as.data.frame(predictSE.lme(fricmod_cdre3, new.df9))
fr.e3$min_cor <- seq(min(cdre.3$min_cor), max(cdre.3$min_cor), by = .01)
fr.e3$lwr <- fr.e3$fit - fr.e3$se.fit
fr.e3$upr <- fr.e3$fit + fr.e3$se.fit
fe.e3 <- as.data.frame(predictSE.lme(fevemod_cdre3, new.df9))
fe.e3$min_cor <- seq(min(cdre.3$min_cor), max(cdre.3$min_cor), by = .01)
fe.e3$lwr <- fe.e3$fit - fe.e3$se.fit
fe.e3$upr <- fe.e3$fit + fe.e3$se.fit
fdis.e3 <- as.data.frame(predictSE.lme(fdismod_cdre3, new.df9))
fdis.e3$min_cor <- seq(min(cdre.3$min_cor), max(cdre.3$min_cor), by = .01)
fdis.e3$lwr <- fdis.e3$fit - fdis.e3$se.fit
fdis.e3$upr <- fdis.e3$fit + fdis.e3$se.fit
fdiv.e3 <- as.data.frame(predictSE.lme(fdivmod_cdre3, new.df9))
fdiv.e3$min_cor <- seq(min(cdre.3$min_cor), max(cdre.3$min_cor), by = .01)
fdiv.e3$lwr <- fdiv.e3$fit - fdiv.e3$se.fit
fdiv.e3$upr <- fdiv.e3$fit + fdiv.e3$se.fit
rq.e3 <- as.data.frame(predictSE.lme(raoqmod_cdre3, new.df9))
rq.e3$min_cor <- seq(min(cdre.3$min_cor), max(cdre.3$min_cor), by = .01)
rq.e3$lwr <- rq.e3$fit - rq.e3$se.fit
rq.e3$upr <- rq.e3$fit + rq.e3$se.fit
kde.alpha.e3 <- as.data.frame(predictSE.lme(kde.alphamod_cdre3, new.df9))
kde.alpha.e3$min_cor <- seq(min(cdre.3$min_cor), max(cdre.3$min_cor), by = .01)
kde.alpha.e3$lwr <- kde.alpha.e3$fit - kde.alpha.e3$se.fit
kde.alpha.e3$upr <- kde.alpha.e3$fit + kde.alpha.e3$se.fit
kde.evenness.e3 <- as.data.frame(predictSE.lme(kde.evennessmod_cdre3, new.df9))
kde.evenness.e3$min_cor <- seq(min(cdre.3$min_cor), max(cdre.3$min_cor), by = .01)
kde.evenness.e3$lwr <- kde.evenness.e3$fit - kde.evenness.e3$se.fit
kde.evenness.e3$upr <- kde.evenness.e3$fit + kde.evenness.e3$se.fit
kde.dispersion.e3 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdre3, new.df9))
kde.dispersion.e3$min_cor <- seq(min(cdre.3$min_cor), max(cdre.3$min_cor), by = .01)
kde.dispersion.e3$lwr <- kde.dispersion.e3$fit - kde.dispersion.e3$se.fit
kde.dispersion.e3$upr <- kde.dispersion.e3$fit + kde.dispersion.e3$se.fit

range(cdre.4$min_cor)
new.df10 <- data.frame(min_cor = seq(min(cdre.4$min_cor), max(cdre.4$min_cor), by = 0.01))
fr.e4 <- as.data.frame(predictSE.lme(fricmod_cdre4, new.df10))
fr.e4$min_cor <- seq(min(cdre.4$min_cor), max(cdre.4$min_cor), by = .01)
fr.e4$lwr <- fr.e4$fit - fr.e4$se.fit
fr.e4$upr <- fr.e4$fit + fr.e4$se.fit
fe.e4 <- as.data.frame(predictSE.lme(fevemod_cdre4, new.df10))
fe.e4$min_cor <- seq(min(cdre.4$min_cor), max(cdre.4$min_cor), by = .01)
fe.e4$lwr <- fe.e4$fit - fe.e4$se.fit
fe.e4$upr <- fe.e4$fit + fe.e4$se.fit
fdis.e4 <- as.data.frame(predictSE.lme(fdismod_cdre4, new.df10))
fdis.e4$min_cor <- seq(min(cdre.4$min_cor), max(cdre.4$min_cor), by = .01)
fdis.e4$lwr <- fdis.e4$fit - fdis.e4$se.fit
fdis.e4$upr <- fdis.e4$fit + fdis.e4$se.fit
fdiv.e4 <- as.data.frame(predictSE.lme(fdivmod_cdre4, new.df10))
fdiv.e4$min_cor <- seq(min(cdre.4$min_cor), max(cdre.4$min_cor), by = .01)
fdiv.e4$lwr <- fdiv.e4$fit - fdiv.e4$se.fit
fdiv.e4$upr <- fdiv.e4$fit + fdiv.e4$se.fit
rq.e4 <- as.data.frame(predictSE.lme(raoqmod_cdre4, new.df10))
rq.e4$min_cor <- seq(min(cdre.4$min_cor), max(cdre.4$min_cor), by = .01)
rq.e4$lwr <- rq.e4$fit - rq.e4$se.fit
rq.e4$upr <- rq.e4$fit + rq.e4$se.fit
kde.alpha.e4 <- as.data.frame(predictSE.lme(kde.alphamod_cdre4, new.df10))
kde.alpha.e4$min_cor <- seq(min(cdre.4$min_cor), max(cdre.4$min_cor), by = .01)
kde.alpha.e4$lwr <- kde.alpha.e4$fit - kde.alpha.e4$se.fit
kde.alpha.e4$upr <- kde.alpha.e4$fit + kde.alpha.e4$se.fit
kde.evenness.e4 <- as.data.frame(predictSE.lme(kde.evennessmod_cdre4, new.df10))
kde.evenness.e4$min_cor <- seq(min(cdre.4$min_cor), max(cdre.4$min_cor), by = .01)
kde.evenness.e4$lwr <- kde.evenness.e4$fit - kde.evenness.e4$se.fit
kde.evenness.e4$upr <- kde.evenness.e4$fit + kde.evenness.e4$se.fit
kde.dispersion.e4 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdre4, new.df10))
kde.dispersion.e4$min_cor <- seq(min(cdre.4$min_cor), max(cdre.4$min_cor), by = .01)
kde.dispersion.e4$lwr <- kde.dispersion.e4$fit - kde.dispersion.e4$se.fit
kde.dispersion.e4$upr <- kde.dispersion.e4$fit + kde.dispersion.e4$se.fit


A <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fr.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= min_cor, y = fit), data = fr.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fr.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= min_cor, y = fit), data = fr.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fr.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= min_cor, y = fit), data = fr.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fr.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= min_cor, y = fit), data = fr.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = min_cor, y = FRic, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "FRich") +
  theme_pubr()

C <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fe.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= min_cor, y = fit), data = fe.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fe.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= min_cor, y = fit), data = fe.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fe.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= min_cor, y = fit), data = fe.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fe.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= min_cor, y = fit), data = fe.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = min_cor, y = FEve, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "FEve") +
  theme_pubr()

E <-ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdis.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= min_cor, y = fit), data = fdis.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdis.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= min_cor, y = fit), data = fdis.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdis.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= min_cor, y = fit), data = fdis.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdis.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= min_cor, y = fit), data = fdis.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = min_cor, y = FDis, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "FDis") +
  theme_pubr()

F <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdis.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= min_cor, y = fit), data = fdis.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdis.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= min_cor, y = fit), data = fdis.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdis.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= min_cor, y = fit), data = fdis.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdis.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= min_cor, y = fit), data = fdis.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = min_cor, y = FDiv, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "FDiv") +
  theme_pubr()

G <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = rq.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= min_cor, y = fit), data = rq.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = rq.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= min_cor, y = fit), data = rq.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = rq.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= min_cor, y = fit), data = rq.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = rq.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= min_cor, y = fit), data = rq.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = min_cor, y = RaoQ, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "Rao Q") +
  theme_pubr()

B <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.alpha.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.alpha.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.alpha.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.alpha.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.alpha.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.alpha.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.alpha.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.alpha.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = min_cor, y = kde.alpha, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "KDE Richness") +
  theme_pubr()

D <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.evenness.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.evenness.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.evenness.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.evenness.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.evenness.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.evenness.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.evenness.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.evenness.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = min_cor, y = kde.evenness, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "KDE Evenness") +
  theme_pubr()

H <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.dispersion.e1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.dispersion.e1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.dispersion.e2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.dispersion.e2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.dispersion.e3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.dispersion.e3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.dispersion.e4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.dispersion.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = min_cor, y = kde.dispersion, color = community), data = cdre.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "KDE Dispersion") +
  theme_pubr()

png(here("Figures/cdre_mincor_euc.png"), height = 5, width = 9, units = 'in', res = 300)
ggarrange(plotlist = list(A, B, C, D, E, F, G, H), common.legend = TRUE,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
dev.off()


###seve####

seve.eblack$community <- "SEV1"
seve.eblue$community <- "SEV2"
seve.full <- rbind(seve.eblue[,-10], seve.eblack[,-10])

seve.full <- seve.full %>% group_by(community, min_cor) %>%
  summarize (FRic = mean(FRic, na.rm = TRUE), FEve = mean(FEve, na.rm = TRUE), 
             FDis = mean(FDis, na.rm = TRUE), FDiv = mean(FDiv, na.rm = TRUE), 
             RaoQ = mean(RaoQ, na.rm = TRUE), kde.alpha = mean(kde.alpha), kde.evenness  = mean(kde.evenness), kde.dispersion = mean(kde.dispersion))


range(seve.eblue$min_cor)
new.df11 <- data.frame(min_cor = seq(min(seve.eblue$min_cor), max(seve.eblue$min_cor), by = 0.01))
fr.eblue <- as.data.frame(predictSE.lme(fricmod_eblue, new.df11))
fr.eblue$min_cor <- seq(min(seve.eblue$min_cor), max(seve.eblue$min_cor), by = 0.01)
fr.eblue$lwr <- fr.eblue$fit - fr.eblue$se.fit
fr.eblue$upr <- fr.eblue$fit + fr.eblue$se.fit
fe.eblue <- as.data.frame(predictSE.lme(fevemod_eblue, new.df11))
fe.eblue$min_cor <- seq(min(seve.eblue$min_cor), max(seve.eblue$min_cor), by = 0.01)
fe.eblue$lwr <- fe.eblue$fit - fe.eblue$se.fit
fe.eblue$upr <- fe.eblue$fit + fe.eblue$se.fit
fdis.eblue <- as.data.frame(predictSE.lme(fdismod_eblue, new.df11))
fdis.eblue$min_cor <- seq(min(seve.eblue$min_cor), max(seve.eblue$min_cor), by = 0.01)
fdis.eblue$lwr <- fdis.eblue$fit - fdis.eblue$se.fit
fdis.eblue$upr <- fdis.eblue$fit + fdis.eblue$se.fit
fdiv.eblue <- as.data.frame(predictSE.lme(fdivmod_eblue, new.df11))
fdiv.eblue$min_cor <- seq(min(seve.eblue$min_cor), max(seve.eblue$min_cor), by = 0.01)
fdiv.eblue$lwr <- fdiv.eblue$fit - fdiv.eblue$se.fit
fdiv.eblue$upr <- fdiv.eblue$fit + fdiv.eblue$se.fit
rq.eblue <- as.data.frame(predictSE.lme(raoqmod_eblue, new.df11))
rq.eblue$min_cor <- seq(min(seve.eblue$min_cor), max(seve.eblue$min_cor), by = 0.01)
rq.eblue$lwr <- rq.eblue$fit - rq.eblue$se.fit
rq.eblue$upr <- rq.eblue$fit + rq.eblue$se.fit
kde.alpha.eblue <- as.data.frame(predictSE.lme(kde.alphamod_eblue, new.df11))
kde.alpha.eblue$min_cor <- seq(min(seve.eblue$min_cor), max(seve.eblue$min_cor), by = .01)
kde.alpha.eblue$lwr <- kde.alpha.eblue$fit - kde.alpha.eblue$se.fit
kde.alpha.eblue$upr <- kde.alpha.eblue$fit + kde.alpha.eblue$se.fit
kde.evenness.eblue <- as.data.frame(predictSE.lme(kde.evennessmod_eblue, new.df11))
kde.evenness.eblue$min_cor <- seq(min(seve.eblue$min_cor), max(seve.eblue$min_cor), by = .01)
kde.evenness.eblue$lwr <- kde.evenness.eblue$fit - kde.evenness.eblue$se.fit
kde.evenness.eblue$upr <- kde.evenness.eblue$fit + kde.evenness.eblue$se.fit
kde.dispersion.eblue <- as.data.frame(predictSE.lme(kde.dispersionmod_eblue, new.df11))
kde.dispersion.eblue$min_cor <- seq(min(seve.eblue$min_cor), max(seve.eblue$min_cor), by = .01)
kde.dispersion.eblue$lwr <- kde.dispersion.eblue$fit - kde.dispersion.eblue$se.fit
kde.dispersion.eblue$upr <- kde.dispersion.eblue$fit + kde.dispersion.eblue$se.fit

range(seve.eblack$min_cor)
new.df12 <- data.frame(min_cor = seq(min(seve.eblack$min_cor), max(seve.eblack$min_cor), by = 0.01))
fr.eblack <- as.data.frame(predictSE.lme(fricmod_eblack, new.df12))
fr.eblack$min_cor <- seq(min(seve.eblack$min_cor), max(seve.eblack$min_cor), by = .01)
fr.eblack$lwr <- fr.eblack$fit - fr.eblack$se.fit
fr.eblack$upr <- fr.eblack$fit + fr.eblack$se.fit
fe.eblack <- as.data.frame(predictSE.lme(fevemod_eblack, new.df12))
fe.eblack$min_cor <- seq(min(seve.eblack$min_cor), max(seve.eblack$min_cor), by = .01)
fe.eblack$lwr <- fe.eblack$fit - fe.eblack$se.fit
fe.eblack$upr <- fe.eblack$fit + fe.eblack$se.fit
fdis.eblack <- as.data.frame(predictSE.lme(fdismod_eblack, new.df12))
fdis.eblack$min_cor <- seq(min(seve.eblack$min_cor), max(seve.eblack$min_cor), by = .01)
fdis.eblack$lwr <- fdis.eblack$fit - fdis.eblack$se.fit
fdis.eblack$upr <- fdis.eblack$fit + fdis.eblack$se.fit
fdiv.eblack <- as.data.frame(predictSE.lme(fdivmod_eblack, new.df12))
fdiv.eblack$min_cor <- seq(min(seve.eblack$min_cor), max(seve.eblack$min_cor), by = .01)
fdiv.eblack$lwr <- fdiv.eblack$fit - fdiv.eblack$se.fit
fdiv.eblack$upr <- fdiv.eblack$fit + fdiv.eblack$se.fit
rq.eblack <- as.data.frame(predictSE.lme(raoqmod_eblack, new.df12))
rq.eblack$min_cor <- seq(min(seve.eblack$min_cor), max(seve.eblack$min_cor), by = .01)
rq.eblack$lwr <- rq.eblack$fit - rq.eblack$se.fit
rq.eblack$upr <- rq.eblack$fit + rq.eblack$se.fit
kde.alpha.eblack <- as.data.frame(predictSE.lme(kde.alphamod_eblack, new.df12))
kde.alpha.eblack$min_cor <- seq(min(seve.eblack$min_cor), max(seve.eblack$min_cor), by = .01)
kde.alpha.eblack$lwr <- kde.alpha.eblack$fit - kde.alpha.eblack$se.fit
kde.alpha.eblack$upr <- kde.alpha.eblack$fit + kde.alpha.eblack$se.fit
kde.evenness.eblack <- as.data.frame(predictSE.lme(kde.evennessmod_eblack, new.df12))
kde.evenness.eblack$min_cor <- seq(min(seve.eblack$min_cor), max(seve.eblack$min_cor), by = .01)
kde.evenness.eblack$lwr <- kde.evenness.eblack$fit - kde.evenness.eblack$se.fit
kde.evenness.eblack$upr <- kde.evenness.eblack$fit + kde.evenness.eblack$se.fit
kde.dispersion.eblack <- as.data.frame(predictSE.lme(kde.dispersionmod_eblack, new.df12))
kde.dispersion.eblack$min_cor <- seq(min(seve.eblack$min_cor), max(seve.eblack$min_cor), by = .01)
kde.dispersion.eblack$lwr <- kde.dispersion.eblack$fit - kde.dispersion.eblack$se.fit
kde.dispersion.eblack$upr <- kde.dispersion.eblack$fit + kde.dispersion.eblack$se.fit

AA <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fr.eblue, alpha = 0.5, fill = "navyeblue") + 
  geom_line(aes(x= min_cor, y = fit), data = fr.eblue, lwd = 2, color = "navyeblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fr.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= min_cor, y = fit), data = fr.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = min_cor, y = FRic, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("#000000", "navyeblue"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "FRich") +
  theme_pubr()

CC <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fe.eblue, alpha = 0.5, fill = "navyeblue") + 
  geom_line(aes(x= min_cor, y = fit), data = fe.eblue, lwd = 2, color = "navyeblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fe.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= min_cor, y = fit), data = fe.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = min_cor, y = FEve, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("#000000", "navyeblue"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "FEve") +
  theme_pubr()

EE <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdis.eblue, alpha = 0.5, fill = "navyeblue") + 
  geom_line(aes(x= min_cor, y = fit), data = fdis.eblue, lwd = 2, color = "navyeblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdis.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= min_cor, y = fit), data = fdis.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = min_cor, y = FDis, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyeblue"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "FDis") +
  theme_pubr()

FF <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdiv.eblue, alpha = 0.5, fill = "navyeblue") + 
  geom_line(aes(x= min_cor, y = fit), data = fdiv.eblue, lwd = 2, color = "navyeblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = fdiv.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= min_cor, y = fit), data = fdiv.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = min_cor, y = FDiv, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyeblue"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "FDiv") +
  theme_pubr()

GG <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = rq.eblue, alpha = 0.5, fill = "navyeblue") + 
  geom_line(aes(x= min_cor, y = fit), data = rq.eblue, lwd = 2, color = "navyeblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = rq.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= min_cor, y = fit), data = rq.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = min_cor, y = RaoQ, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyeblue"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "Rao Q") +
  theme_pubr()

BB <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.alpha.eblue, alpha = 0.5, fill = "navyeblue") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.alpha.eblue, lwd = 2, color = "navyeblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.alpha.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.alpha.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = min_cor, y = kde.alpha, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyeblue"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "KDE Richness") +
  theme_pubr()

DD <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.evenness.eblue, alpha = 0.5, fill = "navyeblue") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.evenness.eblue, lwd = 2, color = "navyeblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.evenness.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.evenness.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = min_cor, y = kde.evenness, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyeblue"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "KDE Evenness") +
  theme_pubr()

HH <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.dispersion.eblue, alpha = 0.5, fill = "navyeblue") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.dispersion.eblue, lwd = 2, color = "navyeblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = min_cor), data = kde.dispersion.eblack, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= min_cor, y = fit), data = kde.dispersion.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = min_cor, y = kde.dispersion, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("black", "navyeblue"), name = "Community") +
  labs(x = "Minimum Trait-Trait Correlation", y = "KDE Dispersion") +
  theme_pubr()

png(here("Figures/seve_mincorr_euc.png"), height = 5, width = 9, units = 'in', res = 300)
ggarrange(plotlist = list(AA, BB, CC, DD, EE, FF, GG, HH), common.legend = TRUE, 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
dev.off()

png(here("Figures/min_corr_euc.png"), height = 8, width = 13, units = 'in', res = 300)
ggarrange(plotlist = list(A, AA, B, BB, C, CC, D, DD, E, EE, F, FF, G, GG, H, HH), 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", 
                     "O", "P"), ncol = 4, nrow = 4, common.legend = TRUE, align = "hv")
dev.off()
