############################################################
#### STATISTICAL ANALYSES  GOWER for TRAIT CORRELATIONS ####
############################################################

# CREATED BY: KAITLIN KIMMEL

# load libraries
library(nlme)
library(here)
library(AICcmodavg)
library(ggplot2)
library(ggpubr)
library(dplyr)


# load data

cdr.1 <- read.csv(here("data/Cleaned/cdr1.csv"))
cdr.1 <- cdr.1[-which(is.na(cdr.1$mean_cor)),]
cdr.2 <- read.csv(here("data/Cleaned/cdr2.csv"))
cdr.2 <- cdr.2[-which(is.na(cdr.2$mean_cor)),]
cdr.3 <- read.csv(here("data/Cleaned/cdr3.csv"))
cdr.3 <- cdr.3[-which(is.na(cdr.3$mean_cor)),]
cdr.4 <- read.csv(here("data/Cleaned/cdr4.csv"))
cdr.4 <- cdr.4[-which(is.na(cdr.4$mean_cor)),]
sev.blue <- read.csv(here("data/Cleaned/sevblue.csv"))
sev.blue <- sev.blue[-which(is.na(sev.blue$mean_cor)),]
sev.black <- read.csv(here("data/Cleaned/sevblack.csv"))
sev.black <- sev.black[-which(is.na(sev.black$mean_cor)),]

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
mod <- lme(FRic ~ mean_cor, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots, sev.blue[-which(is.na(sev.blue$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FRic ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FRic ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod, mod1)
fricmod_blue <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE richness
mod <- lme(kde.alpha ~ mean_cor, random = ~1|sev.blueplots,  data = sev.blue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots, sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod, mod1) ## sig dif use mod1
kde.alphamod_blue <- lme(kde.alpha ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

### FEve
mod <- lme(FEve ~ mean_cor, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FEve ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FEve ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FEve ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1,mod2, mod3) ## mod
fevemod_blue <- lme(FEve ~ mean_cor, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE evenness
mod <- lme(kde.evenness ~ mean_cor, random = ~1|sev.blueplots,  data = sev.blue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots, sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod2, mod3) ### mod3
anova(mod2, mod3) ## sig dif use mod3
kde.evennessmod_blue <- lme(kde.evenness ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

## FDis
mod <- lme(FDis ~ mean_cor, random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FDis ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FDis ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FDis ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1,mod2, mod3) #mod best fit
fdismod_blue <- lme(FDis ~ mean_cor, random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE dispersion
mod <- lme(kde.dispersion ~ mean_cor, random = ~1|sev.blueplots,  data = sev.blue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots, sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod, mod1) #sig diff use mod1
kde.dispersionmod_blue <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

# FDiv
mod <- lme(FDiv ~ mean_cor, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FDiv ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FDiv ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FDiv ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1,mod2, mod3) #mod
fdivmod_blue <- lme(FDiv ~ mean_cor, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

### Rao's Q
mod <- lme(RaoQ ~ mean_cor, random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(RaoQ ~ mean_cor + I(mean_cor^2), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots),control = lmeControl(opt = "optim"))
mod2 <- lme(RaoQ ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(RaoQ ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1,mod2, mod3) #mod3
anova(mod3, mod1) ##sig dif use mod3
raoqmod_blue <- lme(RaoQ ~ mean_cor + I(mean_cor^2)+ I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))


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


###### Fit Metric ~ mean_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ mean_cor, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots, sev.black[-which(is.na(sev.black$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FRic ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FRic ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod2, mod3) ### mod1 lowest
anova(mod, mod1) # sig diff use mod1
fricmod_black <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE richness
mod <- lme(kde.alpha ~ mean_cor, random = ~1|sev.blackplots,  data = sev.black, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots, sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod2, mod3) ### mod1 lowest
anova(mod, mod1) # sig diff use mod1
kde.alphamod_black <- lme(kde.alpha ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))

### FEve
mod <- lme(FEve ~ mean_cor, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FEve ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FEve ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FEve ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1,mod2, mod3) ## mod
fevemod_black <- lme(FEve ~ mean_cor, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE evenness
mod <- lme(kde.evenness ~ mean_cor, random = ~1|sev.blackplots,  data = sev.black, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots, sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod2, mod3) ### mod lowest
kde.evennessmod_black <- lme(kde.evenness ~ mean_cor, random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))

## FDis
mod <- lme(FDis ~ mean_cor, random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FDis ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FDis ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FDis ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1,mod2, mod3) #mod best fit
fdismod_black <- lme(FDis ~ mean_cor, random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE dispersion
mod <- lme(kde.dispersion ~ mean_cor, random = ~1|sev.blackplots,  data = sev.black, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots, sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod2, mod3) ### mod1 lowest
anova(mod, mod1) # sig diff use mod1
kde.dispersionmod_black <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))

# FDiv
mod <- lme(FDiv ~ mean_cor, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FDiv ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FDiv ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FDiv ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1,mod2, mod3) #mod best fit
fdivmod_black <- lme(FDiv ~ mean_cor, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

### Rao's Q
mod <- lme(RaoQ ~ mean_cor, random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(RaoQ ~ mean_cor + I(mean_cor^2), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots),control = lmeControl(opt = "optim"))
mod2 <- lme(RaoQ ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(RaoQ ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1,mod2, mod3) #mod1
anova(mod, mod1) # all similar - using mod
raoqmod_black <- lme(RaoQ ~ mean_cor, random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))


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

###### Fit Metric ~ mean_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ mean_cor, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FRic ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FRic ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod2
fricmod_cdr1 <- lme(FRic ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ mean_cor, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod2
anova(mod1, mod2) ##sig dif use mod2
kde.alphamod_cdr1 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ mean_cor, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) ## mod
fevemod_cdr1 <- lme(FEve ~ mean_cor, random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ mean_cor, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod
kde.evennessmod_cdr1 <- lme(kde.evenness ~ mean_cor, random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ mean_cor, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod3
fdismod_cdr1 <- lme(FDis ~ mean_cor+ I(mean_cor^2)+ I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

##KDE dispersion
mod <- lme(kde.dispersion ~ mean_cor, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod3
kde.dispersionmod_cdr1 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ mean_cor, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod3 best fit
fdivmod_cdr1 <- lme(FDiv ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ mean_cor, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod3 best fit
raoqmod_cdr1 <- lme(RaoQ ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))


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
###### Fit Metric ~ mean_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ mean_cor, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FRic ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FRic ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod3
anova(mod2, mod3) ##sig dif use mod3
fricmod_cdr2 <- lme(FRic ~ mean_cor + I(mean_cor^2)+ I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ mean_cor, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod3
kde.alphamod_cdr2 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2)+ I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ mean_cor, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) ## mod best fit
fevemod_cdr2 <- lme(FEve ~ mean_cor, random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ mean_cor, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod, mod1) ##sig dif use mod1
kde.evennessmod_cdr2 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ mean_cor, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod3 best fit
fdismod_cdr2 <- lme(FDis ~ mean_cor+ I(mean_cor^2)+ I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

##KDE disperion
mod <- lme(kde.evenness ~ mean_cor, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod, mod1) ##sig dif use mod1
kde.dispersionmod_cdr2 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ mean_cor, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod3 best fit
anova(mod2, mod3) ##sig diff use mod3
fdivmod_cdr2 <- lme(FDiv ~ mean_cor+ I(mean_cor^2)+ I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ mean_cor, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod2
anova(mod1, mod2) ##sig dif use mod2
raoqmod_cdr2 <- lme(RaoQ ~ mean_cor+ I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))


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
###### Fit Metric ~ mean_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ mean_cor, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FRic ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FRic ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod2
fricmod_cdr3 <- lme(FRic ~ mean_cor+ I(mean_cor^2)+ I(mean_cor^3), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ mean_cor, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod3
anova(mod2, mod3) #sig dif use mod3
kde.alphamod_cdr3 <- lme(kde.alpha ~ mean_cor+ I(mean_cor^2)+ I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ mean_cor, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) ## mod best fit
fevemod_cdr3 <- lme(FEve ~ mean_cor, random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ mean_cor, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod2
anova(mod1, mod2) ##no sif dif use mod1
kde.evennessmod_cdr3 <- lme(kde.evenness ~ mean_cor+ I(mean_cor^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ mean_cor, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod1, mod2, mod3 similar
anova(mod1, mod3, mod2) # mod1 for simiplicity bc mod3 and mod1 not sig diff
fdismod_cdr3 <- lme(FDis ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

##KDE dispersion
mod <- lme(kde.dispersion ~ mean_cor, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod3
anova(mod1, mod3) #no sig dif use mod1
kde.dispersionmod_cdr3 <- lme(kde.dispersion ~ mean_cor+ I(mean_cor^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ mean_cor, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod1 best fit
fdivmod_cdr3 <- lme(FDiv ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ mean_cor, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod1-mod 3 similar - using mod1 for simplicity
raoqmod_cdr3 <- lme(RaoQ ~ mean_cor+ I(mean_cor^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))


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
###### Fit Metric ~ mean_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ mean_cor, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FRic ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FRic ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod3
anova(mod1, mod3) #no sig dif use mod1
fricmod_cdr4 <- lme(FRic ~ mean_cor+ I(mean_cor^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ mean_cor, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.alpha ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod3
anova(mod2, mod3) #sig dif use mod3
kde.alphamod_cdr4 <- lme(kde.alpha ~ mean_cor+ I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ mean_cor, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) ## mod  
fevemod_cdr4 <- lme(FEve ~ mean_cor, random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ mean_cor, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.evenness ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod3
anova(mod2, mod3) ##sig dif use mod3
kde.evennessmod_cdr4 <- lme(kde.evenness ~ mean_cor+ I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ mean_cor, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod1
fdismod_cdr4 <- lme(FDis ~ mean_cor+ I(mean_cor^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE dispersion
mod <- lme(kde.dispersion ~ mean_cor, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.dispersion ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod3
kde.dispersionmod_cdr4 <- lme(kde.dispersion ~ mean_cor+ I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ mean_cor, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod1 best fit
fdivmod_cdr4 <- lme(FDiv ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ mean_cor, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ mean_cor + I(mean_cor^2) + I(mean_cor^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ mean_cor + I(mean_cor^2) + I(mean_cor^3) + I(mean_cor^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod1-mod3 similar
anova(mod1, mod3, mod2) # all not sig diff - using mod1
raoqmod_cdr4 <- lme(RaoQ ~ mean_cor + I(mean_cor^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod_cdr4)
summary(fevemod_cdr4)
summary(fdismod_cdr4)
summary(fdivmod_cdr4)
summary(raoqmod_cdr4)
summary(kde.alphamod_cdr4)
summary(kde.evennessmod_cdr4)
summary(kde.dispersionmod_cdr4)

################### GRAPHS ##############################

# create combined data set of raw data
cdr.1$community <- "ambient"
cdr.2$community <- "elevated_N"
cdr.3$community <- "elevated_CO2"
cdr.4$community <- "elevated_N_CO2"
cdr.full <- rbind(cdr.1, cdr.2, cdr.3, cdr.4)

cdr.full <- cdr.full %>% group_by(community, mean_cor) %>%
  summarize (FRic = mean(FRic), FEve = mean(FEve), FDis = mean(FDis),
             FDiv = mean(FDiv), RaoQ = mean(RaoQ), kde.alpha = mean(kde.alpha), kde.evenness  = mean(kde.evenness), kde.dispersion = mean(kde.dispersion))

new.df <- data.frame(mean_cor = seq(-0.3, .8, by = .1))


##############################################
############## CDR 1 ########################
#############################################
fr.1 <- as.data.frame(predictSE.lme(fricmod_cdr1, new.df))
fr.1$mean_cor <- seq(-0.3, .8, by = .1)
fr.1$lwr <- fr.1$fit - fr.1$se.fit
fr.1$upr <- fr.1$fit + fr.1$se.fit
fe.1 <- as.data.frame(predictSE.lme(fevemod_cdr1, new.df))
fe.1$mean_cor <- seq(-0.3, .8, by = .1)
fe.1$lwr <- fe.1$fit - fe.1$se.fit
fe.1$upr <- fe.1$fit + fe.1$se.fit
fdis.1 <- as.data.frame(predictSE.lme(fdismod_cdr1, new.df))
fdis.1$mean_cor <- seq(-0.3, .8, by = .1)
fdis.1$lwr <- fdis.1$fit - fdis.1$se.fit
fdis.1$upr <- fdis.1$fit + fdis.1$se.fit
fdiv.1 <- as.data.frame(predictSE.lme(fdivmod_cdr1, new.df))
fdiv.1$mean_cor <- seq(-0.3, .8, by = .1)
fdiv.1$lwr <- fdiv.1$fit - fdiv.1$se.fit
fdiv.1$upr <- fdiv.1$fit + fdiv.1$se.fit
rq.1 <- as.data.frame(predictSE.lme(raoqmod_cdr1, new.df))
rq.1$mean_cor <- seq(-0.3, .8, by = .1)
rq.1$lwr <- rq.1$fit - rq.1$se.fit
rq.1$upr <- rq.1$fit + rq.1$se.fit
kde.alpha.1 <- as.data.frame(predictSE.lme(kde.alphamod_cdr1, new.df))
kde.alpha.1$mean_cor <- seq(-0.3, .8, by = .1)
kde.alpha.1$lwr <- kde.alpha.1$fit - kde.alpha.1$se.fit
kde.alpha.1$upr <- kde.alpha.1$fit + kde.alpha.1$se.fit
kde.evenness.1 <- as.data.frame(predictSE.lme(kde.evennessmod_cdr1, new.df))
kde.evenness.1$mean_cor <- seq(-0.3, .8, by = .1)
kde.evenness.1$lwr <- kde.evenness.1$fit - kde.evenness.1$se.fit
kde.evenness.1$upr <- kde.evenness.1$fit + kde.evenness.1$se.fit
kde.dispersion.1 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdr1, new.df))
kde.dispersion.1$mean_cor <- seq(-0.3, .8, by = .1)
kde.dispersion.1$lwr <- kde.dispersion.1$fit - kde.dispersion.1$se.fit
kde.dispersion.1$upr <- kde.dispersion.1$fit + kde.dispersion.1$se.fit

fr.2 <- as.data.frame(predictSE.lme(fricmod_cdr2, new.df))
fr.2$mean_cor <- seq(-0.3, .8, by = .1)
fr.2$lwr <- fr.2$fit - fr.2$se.fit
fr.2$upr <- fr.2$fit + fr.2$se.fit
fe.2 <- as.data.frame(predictSE.lme(fevemod_cdr2, new.df))
fe.2$mean_cor <- seq(-0.3, .8, by = .1)
fe.2$lwr <- fe.2$fit - fe.2$se.fit
fe.2$upr <- fe.2$fit + fe.2$se.fit
fdis.2 <- as.data.frame(predictSE.lme(fdismod_cdr2, new.df))
fdis.2$mean_cor <- seq(-0.3, .8, by = .1)
fdis.2$lwr <- fdis.2$fit - fdis.2$se.fit
fdis.2$upr <- fdis.2$fit + fdis.2$se.fit
fdiv.2 <- as.data.frame(predictSE.lme(fdivmod_cdr2, new.df))
fdiv.2$mean_cor <- seq(-0.3, .8, by = .1)
fdiv.2$lwr <- fdiv.2$fit - fdiv.2$se.fit
fdiv.2$upr <- fdiv.2$fit + fdiv.2$se.fit
rq.2 <- as.data.frame(predictSE.lme(raoqmod_cdr2, new.df))
rq.2$mean_cor <- seq(-0.3, .8, by = .1)
rq.2$lwr <- rq.2$fit - rq.2$se.fit
rq.2$upr <- rq.2$fit + rq.2$se.fit
kde.alpha.2 <- as.data.frame(predictSE.lme(kde.alphamod_cdr2, new.df))
kde.alpha.2$mean_cor <- seq(-0.3, .8, by = .1)
kde.alpha.2$lwr <- kde.alpha.2$fit - kde.alpha.2$se.fit
kde.alpha.2$upr <- kde.alpha.2$fit + kde.alpha.2$se.fit
kde.evenness.2 <- as.data.frame(predictSE.lme(kde.evennessmod_cdr2, new.df))
kde.evenness.2$mean_cor <- seq(-0.3, .8, by = .1)
kde.evenness.2$lwr <- kde.evenness.2$fit - kde.evenness.2$se.fit
kde.evenness.2$upr <- kde.evenness.2$fit + kde.evenness.2$se.fit
kde.dispersion.2 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdr2, new.df))
kde.dispersion.2$mean_cor <- seq(-0.3, .8, by = .1)
kde.dispersion.2$lwr <- kde.dispersion.2$fit - kde.dispersion.2$se.fit
kde.dispersion.2$upr <- kde.dispersion.2$fit + kde.dispersion.2$se.fit

fr.3 <- as.data.frame(predictSE.lme(fricmod_cdr3, new.df))
fr.3$mean_cor <- seq(-0.3, .8, by = .1)
fr.3$lwr <- fr.3$fit - fr.3$se.fit
fr.3$upr <- fr.3$fit + fr.3$se.fit
fe.3 <- as.data.frame(predictSE.lme(fevemod_cdr3, new.df))
fe.3$mean_cor <- seq(-0.3, .8, by = .1)
fe.3$lwr <- fe.3$fit - fe.3$se.fit
fe.3$upr <- fe.3$fit + fe.3$se.fit
fdis.3 <- as.data.frame(predictSE.lme(fdismod_cdr3, new.df))
fdis.3$mean_cor <- seq(-0.3, .8, by = .1)
fdis.3$lwr <- fdis.3$fit - fdis.3$se.fit
fdis.3$upr <- fdis.3$fit + fdis.3$se.fit
fdiv.3 <- as.data.frame(predictSE.lme(fdivmod_cdr3, new.df))
fdiv.3$mean_cor <- seq(-0.3, .8, by = .1)
fdiv.3$lwr <- fdiv.3$fit - fdiv.3$se.fit
fdiv.3$upr <- fdiv.3$fit + fdiv.3$se.fit
rq.3 <- as.data.frame(predictSE.lme(raoqmod_cdr3, new.df))
rq.3$mean_cor <- seq(-0.3, .8, by = .1)
rq.3$lwr <- rq.3$fit - rq.3$se.fit
rq.3$upr <- rq.3$fit + rq.3$se.fit
kde.alpha.3 <- as.data.frame(predictSE.lme(kde.alphamod_cdr3, new.df))
kde.alpha.3$mean_cor <- seq(-0.3, .8, by = .1)
kde.alpha.3$lwr <- kde.alpha.3$fit - kde.alpha.3$se.fit
kde.alpha.3$upr <- kde.alpha.3$fit + kde.alpha.3$se.fit
kde.evenness.3 <- as.data.frame(predictSE.lme(kde.evennessmod_cdr3, new.df))
kde.evenness.3$mean_cor <- seq(-0.3, .8, by = .1)
kde.evenness.3$lwr <- kde.evenness.3$fit - kde.evenness.3$se.fit
kde.evenness.3$upr <- kde.evenness.3$fit + kde.evenness.3$se.fit
kde.dispersion.3 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdr3, new.df))
kde.dispersion.3$mean_cor <- seq(-0.3, .8, by = .1)
kde.dispersion.3$lwr <- kde.dispersion.3$fit - kde.dispersion.3$se.fit
kde.dispersion.3$upr <- kde.dispersion.3$fit + kde.dispersion.3$se.fit

fr.4 <- as.data.frame(predictSE.lme(fricmod_cdr4, new.df))
fr.4$mean_cor <- seq(-0.3, .8, by = .1)
fr.4$lwr <- fr.4$fit - fr.4$se.fit
fr.4$upr <- fr.4$fit + fr.4$se.fit
fe.4 <- as.data.frame(predictSE.lme(fevemod_cdr4, new.df))
fe.4$mean_cor <- seq(-0.3, .8, by = .1)
fe.4$lwr <- fe.4$fit - fe.4$se.fit
fe.4$upr <- fe.4$fit + fe.4$se.fit
fdis.4 <- as.data.frame(predictSE.lme(fdismod_cdr4, new.df))
fdis.4$mean_cor <- seq(-0.3, .8, by = .1)
fdis.4$lwr <- fdis.4$fit - fdis.4$se.fit
fdis.4$upr <- fdis.4$fit + fdis.4$se.fit
fdiv.4 <- as.data.frame(predictSE.lme(fdivmod_cdr4, new.df))
fdiv.4$mean_cor <- seq(-0.3, .8, by = .1)
fdiv.4$lwr <- fdiv.4$fit - fdiv.4$se.fit
fdiv.4$upr <- fdiv.4$fit + fdiv.4$se.fit
rq.4 <- as.data.frame(predictSE.lme(raoqmod_cdr4, new.df))
rq.4$mean_cor <- seq(-0.3, .8, by = .1)
rq.4$lwr <- rq.4$fit - rq.4$se.fit
rq.4$upr <- rq.4$fit + rq.4$se.fit
kde.alpha.4 <- as.data.frame(predictSE.lme(kde.alphamod_cdr4, new.df))
kde.alpha.4$mean_cor <- seq(-0.3, .8, by = .1)
kde.alpha.4$lwr <- kde.alpha.4$fit - kde.alpha.4$se.fit
kde.alpha.4$upr <- kde.alpha.4$fit + kde.alpha.4$se.fit
kde.evenness.4 <- as.data.frame(predictSE.lme(kde.evennessmod_cdr4, new.df))
kde.evenness.4$mean_cor <- seq(-0.3, .8, by = .1)
kde.evenness.4$lwr <- kde.evenness.4$fit - kde.evenness.4$se.fit
kde.evenness.4$upr <- kde.evenness.4$fit + kde.evenness.4$se.fit
kde.dispersion.4 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdr4, new.df))
kde.dispersion.4$mean_cor <- seq(-0.3, .8, by = .1)
kde.dispersion.4$lwr <- kde.dispersion.4$fit - kde.dispersion.4$se.fit
kde.dispersion.4$upr <- kde.dispersion.4$fit + kde.dispersion.4$se.fit



A <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fr.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= mean_cor, y = fit), data = fr.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fr.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= mean_cor, y = fit), data = fr.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fr.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= mean_cor, y = fit), data = fr.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fr.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= mean_cor, y = fit), data = fr.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = mean_cor, y = FRic, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Mean correlation", y = "FRich") +
  theme_pubr()

B <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fe.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= mean_cor, y = fit), data = fe.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fe.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= mean_cor, y = fit), data = fe.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fe.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= mean_cor, y = fit), data = fe.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fe.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= mean_cor, y = fit), data = fe.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = mean_cor, y = FEve, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Mean correlation", y = "FEve") +
  theme_pubr()

C <-ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdis.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdis.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdis.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdis.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdis.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdis.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdis.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdis.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = mean_cor, y = FDis, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Mean correlation", y = "FDis") +
  theme_pubr()

D <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdiv.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdiv.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdiv.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdiv.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdiv.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdiv.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdiv.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdiv.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = mean_cor, y = FDiv, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Mean correlation", y = "FDiv") +
  theme_pubr()

E <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = rq.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= mean_cor, y = fit), data = rq.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = rq.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= mean_cor, y = fit), data = rq.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = rq.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= mean_cor, y = fit), data = rq.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = rq.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= mean_cor, y = fit), data = rq.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = mean_cor, y = RaoQ, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Mean correlation", y = "Rao Q") +
  theme_pubr()

F <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.alpha.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.alpha.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.alpha.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.alpha.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.alpha.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.alpha.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.alpha.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.alpha.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = mean_cor, y = kde.alpha, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Mean correlation", y = "kde.alpha") +
  theme_pubr()

G <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.evenness.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.evenness.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.evenness.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.evenness.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.evenness.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.evenness.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.evenness.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.evenness.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = mean_cor, y = kde.evenness, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Mean correlation", y = "kde.evenness") +
  theme_pubr()

H <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.dispersion.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.dispersion.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.dispersion.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.dispersion.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.dispersion.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.dispersion.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.dispersion.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.dispersion.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = mean_cor, y = kde.dispersion, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Mean correlation", y = "kde.dispersion") +
  theme_pubr()

png(here("Figures/cdr_meancor.png"), height = 5, width = 9, units = 'in', res = 300)
ggarrange(plotlist = list(A, B, C, D, E, F, G, H), common.legend = TRUE,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
dev.off()


###SEV####


# create combined data set of raw data
sev.black$community <- "black"
sev.blue$community <- "blue"
sev.full <- rbind(sev.blue[,-10], sev.black[,-10])

sev.full <- sev.full %>% group_by(community, mean_cor) %>%
  summarize (FRic = mean(FRic, na.rm = TRUE), FEve = mean(FEve, na.rm = TRUE), 
             FDis = mean(FDis, na.rm = TRUE), FDiv = mean(FDiv, na.rm = TRUE), 
             RaoQ = mean(RaoQ, na.rm = TRUE), kde.alpha = mean(kde.alpha), kde.evenness = mean(kde.evenness), kde.dispersion = mean(kde.dispersion))

new.df <- data.frame(mean_cor = seq(-0.3, .4, by = .1))

fr.blue <- as.data.frame(predictSE.lme(fricmod_blue, new.df))
fr.blue$mean_cor <- seq(-0.3, .4, by = .1)
fr.blue$lwr <- fr.blue$fit - fr.blue$se.fit
fr.blue$upr <- fr.blue$fit + fr.blue$se.fit
fe.blue <- as.data.frame(predictSE.lme(fevemod_blue, new.df))
fe.blue$mean_cor <- seq(-0.3, .4, by = .1)
fe.blue$lwr <- fe.blue$fit - fe.blue$se.fit
fe.blue$upr <- fe.blue$fit + fe.blue$se.fit
fdis.blue <- as.data.frame(predictSE.lme(fdismod_blue, new.df))
fdis.blue$mean_cor <- seq(-0.3, .4, by = .1)
fdis.blue$lwr <- fdis.blue$fit - fdis.blue$se.fit
fdis.blue$upr <- fdis.blue$fit + fdis.blue$se.fit
fdiv.blue <- as.data.frame(predictSE.lme(fdivmod_blue, new.df))
fdiv.blue$mean_cor <- seq(-0.3, .4, by = .1)
fdiv.blue$lwr <- fdiv.blue$fit - fdiv.blue$se.fit
fdiv.blue$upr <- fdiv.blue$fit + fdiv.blue$se.fit
rq.blue <- as.data.frame(predictSE.lme(raoqmod_blue, new.df))
rq.blue$mean_cor <- seq(-0.3, .4, by = .1)
rq.blue$lwr <- rq.blue$fit - rq.blue$se.fit
rq.blue$upr <- rq.blue$fit + rq.blue$se.fit
kde.alpha.blue <- as.data.frame(predictSE.lme(kde.alphamod_blue, new.df))
kde.alpha.blue$mean_cor <- seq(-0.3, .4, by = .1)
kde.alpha.blue$lwr <- kde.alpha.blue$fit - kde.alpha.blue$se.fit
kde.alpha.blue$upr <- kde.alpha.blue$fit + kde.alpha.blue$se.fit
kde.evenness.blue <- as.data.frame(predictSE.lme(kde.evennessmod_blue, new.df))
kde.evenness.blue$mean_cor <- seq(-0.3, .4, by = .1)
kde.evenness.blue$lwr <- kde.evenness.blue$fit - kde.evenness.blue$se.fit
kde.evenness.blue$upr <- kde.evenness.blue$fit + kde.evenness.blue$se.fit
kde.dispersion.blue <- as.data.frame(predictSE.lme(kde.dispersionmod_blue, new.df))
kde.dispersion.blue$mean_cor <- seq(-0.3, .4, by = .1)
kde.dispersion.blue$lwr <- kde.dispersion.blue$fit - kde.dispersion.blue$se.fit
kde.dispersion.blue$upr <- kde.dispersion.blue$fit + kde.dispersion.blue$se.fit


fr.black <- as.data.frame(predictSE.lme(fricmod_black, new.df))
fr.black$mean_cor <- seq(-0.3, .4, by = .1)
fr.black$lwr <- fr.black$fit - fr.black$se.fit
fr.black$upr <- fr.black$fit + fr.black$se.fit
fe.black <- as.data.frame(predictSE.lme(fevemod_black, new.df))
fe.black$mean_cor <- seq(-0.3, .4, by = .1)
fe.black$lwr <- fe.black$fit - fe.black$se.fit
fe.black$upr <- fe.black$fit + fe.black$se.fit
fdis.black <- as.data.frame(predictSE.lme(fdismod_black, new.df))
fdis.black$mean_cor <- seq(-0.3, .4, by = .1)
fdis.black$lwr <- fdis.black$fit - fdis.black$se.fit
fdis.black$upr <- fdis.black$fit + fdis.black$se.fit
fdiv.black <- as.data.frame(predictSE.lme(fdivmod_black, new.df))
fdiv.black$mean_cor <- seq(-0.3, .4, by = .1)
fdiv.black$lwr <- fdiv.black$fit - fdiv.black$se.fit
fdiv.black$upr <- fdiv.black$fit + fdiv.black$se.fit
rq.black <- as.data.frame(predictSE.lme(raoqmod_black, new.df))
rq.black$mean_cor <- seq(-0.3, .4, by = .1)
rq.black$lwr <- rq.black$fit - rq.black$se.fit
rq.black$upr <- rq.black$fit + rq.black$se.fit
kde.alpha.black <- as.data.frame(predictSE.lme(kde.alphamod_black, new.df))
kde.alpha.black$mean_cor <- seq(-0.3, .4, by = .1)
kde.alpha.black$lwr <- kde.alpha.black$fit - kde.alpha.black$se.fit
kde.alpha.black$upr <- kde.alpha.black$fit + kde.alpha.black$se.fit
kde.evenness.black <- as.data.frame(predictSE.lme(kde.evennessmod_black, new.df))
kde.evenness.black$mean_cor <- seq(-0.3, .4, by = .1)
kde.evenness.black$lwr <- kde.evenness.black$fit - kde.evenness.black$se.fit
kde.evenness.black$upr <- kde.evenness.black$fit + kde.evenness.black$se.fit
kde.dispersion.black <- as.data.frame(predictSE.lme(kde.dispersionmod_black, new.df))
kde.dispersion.black$mean_cor <- seq(-0.3, .4, by = .1)
kde.dispersion.black$lwr <- kde.dispersion.black$fit - kde.dispersion.black$se.fit
kde.dispersion.black$upr <- kde.dispersion.black$fit + kde.dispersion.black$se.fit

A <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fr.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= mean_cor, y = fit), data = fr.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fr.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= mean_cor, y = fit), data = fr.black, lwd = 2, color = "black") + 
  geom_point(aes(x = mean_cor, y = FRic, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#000000", "navyblue"), name = "Community") +
  labs(x = "Mean correlation", y = "FRich") +
  theme_pubr()

B <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fe.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= mean_cor, y = fit), data = fe.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fe.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= mean_cor, y = fit), data = fe.black, lwd = 2, color = "black") + 
  geom_point(aes(x = mean_cor, y = FEve, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#000000", "navyblue"), name = "Community") +
  labs(x = "Mean correlation", y = "FEve") +
  theme_pubr()

C <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdis.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdis.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdis.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdis.black, lwd = 2, color = "black") + 
  geom_point(aes(x = mean_cor, y = FDis, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Mean correlation", y = "FDis") +
  theme_pubr()

D <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdiv.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdiv.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = fdiv.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= mean_cor, y = fit), data = fdiv.black, lwd = 2, color = "black") + 
  geom_point(aes(x = mean_cor, y = FDiv, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Mean correlation", y = "FDiv") +
  theme_pubr()

E <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = rq.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= mean_cor, y = fit), data = rq.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = rq.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= mean_cor, y = fit), data = rq.black, lwd = 2, color = "black") + 
  geom_point(aes(x = mean_cor, y = RaoQ, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Mean correlation", y = "Rao Q") +
  theme_pubr()

F <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.alpha.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.alpha.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.alpha.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.alpha.black, lwd = 2, color = "black") + 
  geom_point(aes(x = mean_cor, y = kde.alpha, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Mean correlation", y = "kde.alpha") +
  theme_pubr()

G <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.evenness.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.evenness.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.evenness.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.evenness.black, lwd = 2, color = "black") + 
  geom_point(aes(x = mean_cor, y = kde.evenness, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Mean correlation", y = "kde.evenness") +
  theme_pubr()

H <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.dispersion.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.dispersion.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = mean_cor), data = kde.dispersion.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= mean_cor, y = fit), data = kde.dispersion.black, lwd = 2, color = "black") + 
  geom_point(aes(x = mean_cor, y = kde.dispersion, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Mean correlation", y = "kde.dispersion") +
  theme_pubr()


png(here("Figures/Sev_meancorr.png"))
ggarrange(plotlist = list(A, B, C, D, E, F, G, H), common.legend = TRUE, 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
dev.off()



