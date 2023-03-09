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

cdr.1 <- read.csv(here("data/Cleaned/cdr1_euc.csv"))
cdr.2 <- read.csv(here("data/Cleaned/cdr2_euc.csv"))
cdr.3 <- read.csv(here("data/Cleaned/cdr3_euc.csv"))
cdr.4 <- read.csv(here("data/Cleaned/cdr4_euc.csv"))
sev.blue <- read.csv(here("data/Cleaned/sevblue_euc.csv"))
#sev.blue <- sev.blue[-which(sev.blue$SR == 1),]
sev.black <- read.csv(here("data/Cleaned/sevblack_euc.csv"))
#sev.black <- sev.black[-which(sev.black$SR == 1),]


#############################################################
################ SEV #######################################
###########################################################

###########################################################
################### BLUE ################################
##########################################################


###### Fit Metric ~ N_trait models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ n_trait, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots, sev.blue[-which(is.na(sev.blue$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod2, mod3) ### mod
fricmod_blue <- lme(FRic ~ n_trait, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE richness
mod <- lme(kde.alpha ~ n_trait, random = ~1|sev.blueplots,  data = sev.blue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots, sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(kde.alpha ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(kde.alpha ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod2, mod3) ### mod1
anova( mod, mod1) # sig diff - use mod 1
kde.alphamod_blue <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1,mod2, mod3) ## mod1
anova(mod, mod1) #not sig diff - using mod
fevemod_blue <- lme(FEve ~ n_trait, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE evenness
mod <- lme(kde.evenness ~ n_trait, random = ~1|sev.blueplots,  data = sev.blue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots, sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(kde.evenness ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(kde.evenness ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod2, mod3) ### mod3
anova( mod3, mod1) # no sig diff - use mod 1
kde.evennessmod_blue <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

## FDis
mod <- lme(FDis ~ n_trait, random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1,mod2, mod3) #mod2 best fit
anova(mod2, mod1) # mod1 no sig
fdismod_blue <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

##KDE dispersion
mod <- lme(kde.dispersion ~ n_trait, random = ~1|sev.blueplots,  data = sev.blue, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots, sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(kde.dispersion ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(kde.dispersion ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod1, mod) # sig diff - use mod 1
kde.dispersionmod_blue <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))


# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1,mod2, mod3) #mod1 best fit
anova(mod, mod1) #mod1 sig diff from mod
fdivmod_blue <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots),control = lmeControl(opt = "optim"))
mod2 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1,mod2, mod3) #mod best fit
raoqmod_blue <- lme(RaoQ ~ n_trait, random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))


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


###### Fit Metric ~ N_trait models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ n_trait, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots, sev.black[-which(is.na(sev.black$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod2, mod3) ### mod, mod1, and mod 2 similar
anova(mod, mod1) #mod for simplicity
fricmod_black <- lme(FRic ~ n_trait, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE richness
mod <- lme(kde.alpha ~ n_trait, random = ~1|sev.blackplots,  data = sev.black, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots, sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(kde.alpha ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(kde.alpha ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod1, mod) #sig diff use mod 1
kde.alphamod_black <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1,mod2, mod3) ## mod
fevemod_black <- lme(FEve ~ n_trait, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE evenness
mod <- lme(kde.evenness ~ n_trait, random = ~1|sev.blackplots,  data = sev.black, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots, sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(kde.evenness ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(kde.evenness ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod2, mod3) ### mod
kde.evennessmod_black <- lme(kde.evenness ~ n_trait, random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))



## FDis
mod <- lme(FDis ~ n_trait, random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1,mod2, mod3) #mod2 best fit
anova( mod2, mod1) # mod2 sign diff
fdismod_black <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE dispersion
mod <- lme(kde.dispersion ~ n_trait, random = ~1|sev.blackplots,  data = sev.black, method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots, sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(kde.dispersion ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(kde.dispersion ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod, mod1) ## sig diff mod1
kde.dispersionmod_black <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))

# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))


AIC(mod, mod1,mod2, mod3) #mod1 best fit
anova(mod, mod1) #mod1 sig diff 
fdivmod_black <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots),control = lmeControl(opt = "optim"))
mod2 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1,mod2, mod3) #mod
raoqmod_black <- lme(RaoQ ~ n_trait, random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))


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

###### Fit Metric ~ N_trait models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod1,mod2, similar AIC values
anova(mod1, mod) #mod1 sig diff
fricmod_cdr1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.alpha ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.alpha ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod, mod1) ## sig diff use mod1
kde.alphamod_cdr1 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) ## mod1 best fit
anova(mod, mod1) # mod1 sig different
fevemod_cdr1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.evenness ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.evenness ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod, mod1) ##sig diff use mod1
kde.evennessmod_cdr1 <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))


## FDis
mod <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod2 best fit
anova(mod2, mod1) # mod2 sig diff
fdismod_cdr1 <- lme(FDis ~ n_trait+ I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

##KDE dispersion
mod <- lme(kde.dispersion ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.dispersion ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.dispersion ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod1, mod) ##sig diff use mod1
kde.dispersionmod_cdr1 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod1 best fit
anova(mod, mod1) #mod1 sig diff
fdivmod_cdr1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod best fit
raoqmod_cdr1 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))


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
###### Fit Metric ~ N_trait models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ n_trait, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod, mod1)
fricmod_cdr2 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ n_trait, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.alpha ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.alpha ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod, mod1) ##sig diff use mod1
kde.alphamod_cdr2 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) ## mod2 best fit
anova(mod, mod2) # sign diff. using mod2
fevemod_cdr2 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ n_trait, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.evenness ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.evenness ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod
kde.evennessmod_cdr2 <- lme(kde.evenness ~ n_trait, random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod2 best fit
anova(mod2, mod1) #mod1 similar to mod 2 - using 1
fdismod_cdr2 <- lme(FDis ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

##KDE dispersion
mod <- lme(kde.dispersion ~ n_trait, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.dispersion ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.dispersion ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod, mod1) ## sig diff use mod1
kde.dispersionmod_cdr2 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod1 best fit
anova(mod, mod1) # sign diff. using mod 1
fdivmod_cdr2 <- lme(FDiv ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))



### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod
raoqmod_cdr2 <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))


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
###### Fit Metric ~ N_trait models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ n_trait, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod3
anova( mod3, mod2) # sig diff
fricmod_cdr3 <- lme(FRic ~ n_trait+ I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ n_trait, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.alpha ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.alpha ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod1
anova( mod, mod1) # sig. diff. using mod1
kde.alphamod_cdr3 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) ## mod2 best fit
anova( mod2, mod1) # sig. diff. using mod2
fevemod_cdr3 <- lme(FEve ~ n_trait+ I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ n_trait, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.evenness ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.evenness ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod2
anova( mod2, mod1) # sig. diff. using mod2
kde.evennessmod_cdr3 <- lme(kde.evenness ~ n_trait+ I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod2 best fit
anova( mod2, mod1) # sig. diff. using mod2
fdismod_cdr3 <- lme(FDis ~ n_trait + I(n_trait^2)+ I(n_trait^3), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

##KDE dispersion
mod <- lme(kde.dispersion ~ n_trait, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.dispersion ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.dispersion ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod, mod1) # sig. diff. using mod1
kde.dispersionmod_cdr3 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod3 best fit
anova(mod, mod3, mod2, mod1) # mod2 and mod3 not sig diff - using mod2
fdivmod_cdr3 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod 
raoqmod_cdr3 <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))


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
###### Fit Metric ~ N_trait models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ n_trait, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod
fricmod_cdr4 <- lme(FRic ~ n_trait, random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness
mod <- lme(kde.alpha ~ n_trait, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.alpha ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.alpha ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod, mod1) # sig diff use mod1
kde.alphamod_cdr4 <- lme(kde.alpha ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) ## mod2  
anova( mod2, mod1) # sig. diff 
fevemod_cdr4 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness
mod <- lme(kde.evenness ~ n_trait, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.evenness ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.evenness ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod3
anova(mod, mod3) # sig diff use mod3
kde.evennessmod_cdr4 <- lme(kde.evenness ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

## FDis
mod <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod1 & mod2 similar - using mod1 for simplicity
fdismod_cdr4 <- lme(FDis ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

##KDE dispersion
mod <- lme(kde.dispersion ~ n_trait, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(kde.dispersion ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(kde.dispersion ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod, mod1) # sig diff use mod1
kde.dispersionmod_cdr4 <- lme(kde.dispersion ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod2 best fit
anova(mod, mod2, mod1) # sig diff. using mod2
fdivmod_cdr4 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod 
raoqmod_cdr4 <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))


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

cdr.full <- cdr.full %>% group_by(community, n_trait) %>%
  summarize (FRic = mean(FRic), FEve = mean(FEve), FDis = mean(FDis),
             FDiv = mean(FDiv), RaoQ = mean(RaoQ), kde.alpha = mean(kde.alpha), kde.evenness = mean(kde.evenness), kde.dispersion = mean(kde.dispersion))

new.df <- data.frame(n_trait = seq(2, 9, by = 1))


##############################################
############## CDR 1 ########################
#############################################
fr.1 <- as.data.frame(predictSE.lme(fricmod_cdr1, new.df))
fr.1$n_trait <- seq(2, 9, by = 1)
fr.1$lwr <- fr.1$fit - fr.1$se.fit
fr.1$upr <- fr.1$fit + fr.1$se.fit
fe.1 <- as.data.frame(predictSE.lme(fevemod_cdr1, new.df))
fe.1$n_trait <- seq(2, 9, by = 1)
fe.1$lwr <- fe.1$fit - fe.1$se.fit
fe.1$upr <- fe.1$fit + fe.1$se.fit
fdis.1 <- as.data.frame(predictSE.lme(fdismod_cdr1, new.df))
fdis.1$n_trait <- seq(2, 9, by = 1)
fdis.1$lwr <- fdis.1$fit - fdis.1$se.fit
fdis.1$upr <- fdis.1$fit + fdis.1$se.fit
fdiv.1 <- as.data.frame(predictSE.lme(fdivmod_cdr1, new.df))
fdiv.1$n_trait <- seq(2, 9, by = 1)
fdiv.1$lwr <- fdiv.1$fit - fdiv.1$se.fit
fdiv.1$upr <- fdiv.1$fit + fdiv.1$se.fit
rq.1 <- as.data.frame(predictSE.lme(raoqmod_cdr1, new.df))
rq.1$n_trait <- seq(2, 9, by = 1)
rq.1$lwr <- rq.1$fit - rq.1$se.fit
rq.1$upr <- rq.1$fit + rq.1$se.fit
kde.alpha.1 <- as.data.frame(predictSE.lme(kde.alphamod_cdr1, new.df))
kde.alpha.1$n_trait <- seq(2, 9, by = 1)
kde.alpha.1$lwr <- kde.alpha.1$fit - kde.alpha.1$se.fit
kde.alpha.1$upr <- kde.alpha.1$fit + kde.alpha.1$se.fit
kde.evenness.1 <- as.data.frame(predictSE.lme(kde.evennessmod_cdr1, new.df))
kde.evenness.1$n_trait <- seq(2, 9, by = 1)
kde.evenness.1$lwr <- kde.evenness.1$fit - kde.evenness.1$se.fit
kde.evenness.1$upr <- kde.evenness.1$fit + kde.evenness.1$se.fit
kde.dispersion.1 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdr1, new.df))
kde.dispersion.1$n_trait <- seq(2, 9, by = 1)
kde.dispersion.1$lwr <- kde.dispersion.1$fit - kde.dispersion.1$se.fit
kde.dispersion.1$upr <- kde.dispersion.1$fit + kde.dispersion.1$se.fit

fr.2 <- as.data.frame(predictSE.lme(fricmod_cdr2, new.df))
fr.2$n_trait <- seq(2, 9, by = 1)
fr.2$lwr <- fr.2$fit - fr.2$se.fit
fr.2$upr <- fr.2$fit + fr.2$se.fit
fe.2 <- as.data.frame(predictSE.lme(fevemod_cdr2, new.df))
fe.2$n_trait <- seq(2, 9, by = 1)
fe.2$lwr <- fe.2$fit - fe.2$se.fit
fe.2$upr <- fe.2$fit + fe.2$se.fit
fdis.2 <- as.data.frame(predictSE.lme(fdismod_cdr2, new.df))
fdis.2$n_trait <- seq(2, 9, by = 1)
fdis.2$lwr <- fdis.2$fit - fdis.2$se.fit
fdis.2$upr <- fdis.2$fit + fdis.2$se.fit
fdiv.2 <- as.data.frame(predictSE.lme(fdivmod_cdr2, new.df))
fdiv.2$n_trait <- seq(2, 9, by = 1)
fdiv.2$lwr <- fdiv.2$fit - fdiv.2$se.fit
fdiv.2$upr <- fdiv.2$fit + fdiv.2$se.fit
rq.2 <- as.data.frame(predictSE.lme(raoqmod_cdr2, new.df))
rq.2$n_trait <- seq(2, 9, by = 1)
rq.2$lwr <- rq.2$fit - rq.2$se.fit
rq.2$upr <- rq.2$fit + rq.2$se.fit
kde.alpha.2 <- as.data.frame(predictSE.lme(kde.alphamod_cdr2, new.df))
kde.alpha.2$n_trait <- seq(2, 9, by = 1)
kde.alpha.2$lwr <- kde.alpha.2$fit - kde.alpha.2$se.fit
kde.alpha.2$upr <- kde.alpha.2$fit + kde.alpha.2$se.fit
kde.evenness.2 <- as.data.frame(predictSE.lme(kde.evennessmod_cdr2, new.df))
kde.evenness.2$n_trait <- seq(2, 9, by = 1)
kde.evenness.2$lwr <- kde.evenness.2$fit - kde.evenness.2$se.fit
kde.evenness.2$upr <- kde.evenness.2$fit + kde.evenness.2$se.fit
kde.dispersion.2 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdr2, new.df))
kde.dispersion.2$n_trait <- seq(2, 9, by = 1)
kde.dispersion.2$lwr <- kde.dispersion.2$fit - kde.dispersion.2$se.fit
kde.dispersion.2$upr <- kde.dispersion.2$fit + kde.dispersion.2$se.fit

fr.3 <- as.data.frame(predictSE.lme(fricmod_cdr3, new.df))
fr.3$n_trait <- seq(2, 9, by = 1)
fr.3$lwr <- fr.3$fit - fr.3$se.fit
fr.3$upr <- fr.3$fit + fr.3$se.fit
fe.3 <- as.data.frame(predictSE.lme(fevemod_cdr3, new.df))
fe.3$n_trait <- seq(2, 9, by = 1)
fe.3$lwr <- fe.3$fit - fe.3$se.fit
fe.3$upr <- fe.3$fit + fe.3$se.fit
fdis.3 <- as.data.frame(predictSE.lme(fdismod_cdr3, new.df))
fdis.3$n_trait <- seq(2, 9, by = 1)
fdis.3$lwr <- fdis.3$fit - fdis.3$se.fit
fdis.3$upr <- fdis.3$fit + fdis.3$se.fit
fdiv.3 <- as.data.frame(predictSE.lme(fdivmod_cdr3, new.df))
fdiv.3$n_trait <- seq(2, 9, by = 1)
fdiv.3$lwr <- fdiv.3$fit - fdiv.3$se.fit
fdiv.3$upr <- fdiv.3$fit + fdiv.3$se.fit
rq.3 <- as.data.frame(predictSE.lme(raoqmod_cdr3, new.df))
rq.3$n_trait <- seq(2, 9, by = 1)
rq.3$lwr <- rq.3$fit - rq.3$se.fit
rq.3$upr <- rq.3$fit + rq.3$se.fit
kde.alpha.3 <- as.data.frame(predictSE.lme(kde.alphamod_cdr3, new.df))
kde.alpha.3$n_trait <- seq(2, 9, by = 1)
kde.alpha.3$lwr <- kde.alpha.3$fit - kde.alpha.3$se.fit
kde.alpha.3$upr <- kde.alpha.3$fit + kde.alpha.3$se.fit
kde.evenness.3 <- as.data.frame(predictSE.lme(kde.evennessmod_cdr3, new.df))
kde.evenness.3$n_trait <- seq(2, 9, by = 1)
kde.evenness.3$lwr <- kde.evenness.3$fit - kde.evenness.3$se.fit
kde.evenness.3$upr <- kde.evenness.3$fit + kde.evenness.3$se.fit
kde.dispersion.3 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdr3, new.df))
kde.dispersion.3$n_trait <- seq(2, 9, by = 1)
kde.dispersion.3$lwr <- kde.dispersion.3$fit - kde.dispersion.3$se.fit
kde.dispersion.3$upr <- kde.dispersion.3$fit + kde.dispersion.3$se.fit

fr.4 <- as.data.frame(predictSE.lme(fricmod_cdr4, new.df))
fr.4$n_trait <- seq(2, 9, by = 1)
fr.4$lwr <- fr.4$fit - fr.4$se.fit
fr.4$upr <- fr.4$fit + fr.4$se.fit
fe.4 <- as.data.frame(predictSE.lme(fevemod_cdr4, new.df))
fe.4$n_trait <- seq(2, 9, by = 1)
fe.4$lwr <- fe.4$fit - fe.4$se.fit
fe.4$upr <- fe.4$fit + fe.4$se.fit
fdis.4 <- as.data.frame(predictSE.lme(fdismod_cdr4, new.df))
fdis.4$n_trait <- seq(2, 9, by = 1)
fdis.4$lwr <- fdis.4$fit - fdis.4$se.fit
fdis.4$upr <- fdis.4$fit + fdis.4$se.fit
fdiv.4 <- as.data.frame(predictSE.lme(fdivmod_cdr4, new.df))
fdiv.4$n_trait <- seq(2, 9, by = 1)
fdiv.4$lwr <- fdiv.4$fit - fdiv.4$se.fit
fdiv.4$upr <- fdiv.4$fit + fdiv.4$se.fit
rq.4 <- as.data.frame(predictSE.lme(raoqmod_cdr4, new.df))
rq.4$n_trait <- seq(2, 9, by = 1)
rq.4$lwr <- rq.4$fit - rq.4$se.fit
rq.4$upr <- rq.4$fit + rq.4$se.fit
kde.alpha.4 <- as.data.frame(predictSE.lme(kde.alphamod_cdr4, new.df))
kde.alpha.4$n_trait <- seq(2, 9, by = 1)
kde.alpha.4$lwr <- kde.alpha.4$fit - kde.alpha.4$se.fit
kde.alpha.4$upr <- kde.alpha.4$fit + kde.alpha.4$se.fit
kde.evenness.4 <- as.data.frame(predictSE.lme(kde.evennessmod_cdr4, new.df))
kde.evenness.4$n_trait <- seq(2, 9, by = 1)
kde.evenness.4$lwr <- kde.evenness.4$fit - kde.evenness.4$se.fit
kde.evenness.4$upr <- kde.evenness.4$fit + kde.evenness.4$se.fit
kde.dispersion.4 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdr4, new.df))
kde.dispersion.4$n_trait <- seq(2, 9, by = 1)
kde.dispersion.4$lwr <- kde.dispersion.4$fit - kde.dispersion.4$se.fit
kde.dispersion.4$upr <- kde.dispersion.4$fit + kde.dispersion.4$se.fit


A <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FRic, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "FRich") +
  theme_pubr()

C <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FEve, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "FEve") +
  theme_pubr()

E <-ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FDis, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "FDis") +
  theme_pubr()

F <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FDiv, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "FDiv") +
  theme_pubr()

G <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = RaoQ, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "Rao Q") +
  theme_pubr()

B <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.alpha, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Richness") +
  theme_pubr()

D <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.evenness, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Evenness") +
  theme_pubr()

H <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.dispersion, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Dispersion") +
  theme_pubr()

png(here("Figures/cdr_Ntraits_euc.png"))
ggarrange(plotlist = list(A, B, C, D, E, F, G, H), common.legend = TRUE,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
dev.off()

###SEV####


# create combined data set of raw data
sev.black$community <- "black"
sev.blue$community <- "blue"
sev.full <- rbind(sev.blue[,-10], sev.black[,-10])

sev.full <- sev.full %>% group_by(community, n_trait) %>%
  summarize (FRic = mean(FRic, na.rm = TRUE), FEve = mean(FEve, na.rm = TRUE), 
             FDis = mean(FDis, na.rm = TRUE), FDiv = mean(FDiv, na.rm = TRUE), 
             RaoQ = mean(RaoQ, na.rm = TRUE), kde.alpha = mean(kde.alpha), kde.evenness = mean(kde.evenness), kde.dispersion = mean(kde.dispersion))

new.df <- data.frame(n_trait = seq(2, 10, by = 1))

fr.blue <- as.data.frame(predictSE.lme(fricmod_blue, new.df))
fr.blue$n_trait <- seq(2, 10, by = 1)
fr.blue$lwr <- fr.blue$fit - fr.blue$se.fit
fr.blue$upr <- fr.blue$fit + fr.blue$se.fit
fe.blue <- as.data.frame(predictSE.lme(fevemod_blue, new.df))
fe.blue$n_trait <- seq(2, 10, by = 1)
fe.blue$lwr <- fe.blue$fit - fe.blue$se.fit
fe.blue$upr <- fe.blue$fit + fe.blue$se.fit
fdis.blue <- as.data.frame(predictSE.lme(fdismod_blue, new.df))
fdis.blue$n_trait <- seq(2, 10, by = 1)
fdis.blue$lwr <- fdis.blue$fit - fdis.blue$se.fit
fdis.blue$upr <- fdis.blue$fit + fdis.blue$se.fit
fdiv.blue <- as.data.frame(predictSE.lme(fdivmod_blue, new.df))
fdiv.blue$n_trait <- seq(2, 10, by = 1)
fdiv.blue$lwr <- fdiv.blue$fit - fdiv.blue$se.fit
fdiv.blue$upr <- fdiv.blue$fit + fdiv.blue$se.fit
rq.blue <- as.data.frame(predictSE.lme(raoqmod_blue, new.df))
rq.blue$n_trait <- seq(2, 10, by = 1)
rq.blue$lwr <- rq.blue$fit - rq.blue$se.fit
rq.blue$upr <- rq.blue$fit + rq.blue$se.fit
kde.alpha.blue <- as.data.frame(predictSE.lme(kde.alphamod_blue, new.df))
kde.alpha.blue$n_trait <- seq(2, 10, by = 1)
kde.alpha.blue$lwr <- kde.alpha.blue$fit - kde.alpha.blue$se.fit
kde.alpha.blue$upr <- kde.alpha.blue$fit + kde.alpha.blue$se.fit
kde.evenness.blue <- as.data.frame(predictSE.lme(kde.evennessmod_blue, new.df))
kde.evenness.blue$n_trait <- seq(2, 10, by = 1)
kde.evenness.blue$lwr <- kde.evenness.blue$fit - kde.evenness.blue$se.fit
kde.evenness.blue$upr <- kde.evenness.blue$fit + kde.evenness.blue$se.fit
kde.dispersion.blue <- as.data.frame(predictSE.lme(kde.dispersionmod_blue, new.df))
kde.dispersion.blue$n_trait <- seq(2, 10, by = 1)
kde.dispersion.blue$lwr <- kde.dispersion.blue$fit - kde.dispersion.blue$se.fit
kde.dispersion.blue$upr <- kde.dispersion.blue$fit + kde.dispersion.blue$se.fit

fr.black <- as.data.frame(predictSE.lme(fricmod_black, new.df))
fr.black$n_trait <- seq(2, 10, by = 1)
fr.black$lwr <- fr.black$fit - fr.black$se.fit
fr.black$upr <- fr.black$fit + fr.black$se.fit
fe.black <- as.data.frame(predictSE.lme(fevemod_black, new.df))
fe.black$n_trait <- seq(2, 10, by = 1)
fe.black$lwr <- fe.black$fit - fe.black$se.fit
fe.black$upr <- fe.black$fit + fe.black$se.fit
fdis.black <- as.data.frame(predictSE.lme(fdismod_black, new.df))
fdis.black$n_trait <- seq(2, 10, by = 1)
fdis.black$lwr <- fdis.black$fit - fdis.black$se.fit
fdis.black$upr <- fdis.black$fit + fdis.black$se.fit
fdiv.black <- as.data.frame(predictSE.lme(fdivmod_black, new.df))
fdiv.black$n_trait <- seq(2, 10, by = 1)
fdiv.black$lwr <- fdiv.black$fit - fdiv.black$se.fit
fdiv.black$upr <- fdiv.black$fit + fdiv.black$se.fit
rq.black <- as.data.frame(predictSE.lme(raoqmod_black, new.df))
rq.black$n_trait <- seq(2, 10, by = 1)
rq.black$lwr <- rq.black$fit - rq.black$se.fit
rq.black$upr <- rq.black$fit + rq.black$se.fit
kde.alpha.black <- as.data.frame(predictSE.lme(kde.alphamod_black, new.df))
kde.alpha.black$n_trait <- seq(2, 10, by = 1)
kde.alpha.black$lwr <- kde.alpha.black$fit - kde.alpha.black$se.fit
kde.alpha.black$upr <- kde.alpha.black$fit + kde.alpha.black$se.fit
kde.evenness.black <- as.data.frame(predictSE.lme(kde.evennessmod_black, new.df))
kde.evenness.black$n_trait <- seq(2, 10, by = 1)
kde.evenness.black$lwr <- kde.evenness.black$fit - kde.evenness.black$se.fit
kde.evenness.black$upr <- kde.evenness.black$fit + kde.evenness.black$se.fit
kde.dispersion.black <- as.data.frame(predictSE.lme(kde.dispersionmod_black, new.df))
kde.dispersion.black$n_trait <- seq(2, 10, by = 1)
kde.dispersion.black$lwr <- kde.dispersion.black$fit - kde.dispersion.black$se.fit
kde.dispersion.black$upr <- kde.dispersion.black$fit + kde.dispersion.black$se.fit

AA <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FRic, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FRich") +
  theme_pubr()

CC <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FEve, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FEve") +
  theme_pubr()

EE <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FDis, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FDis") +
  theme_pubr()

FF <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FDiv, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FDiv") +
  theme_pubr()

GG <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = RaoQ, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "Rao Q") +
  theme_pubr()

BB <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.alpha, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Richness") +
  theme_pubr()

DD <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.evenness, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Eveness") +
  theme_pubr()

HH <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.dispersion, color = community), data = sev.full, size = 3) + 
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
