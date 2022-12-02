################################
#### STATISTICAL ANALYSES   ####
################################

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
cdr.2 <- read.csv(here("data/Cleaned/cdr2.csv"))
cdr.3 <- read.csv(here("data/Cleaned/cdr3.csv"))
cdr.4 <- read.csv(here("data/Cleaned/cdr4.csv"))
sev.blue <- read.csv(here("data/Cleaned/sevblue.csv"))
sev.blue <- sev.blue[-which(sev.blue$SR == 1),]
sev.black <- read.csv(here("data/Cleaned/sevblack.csv"))
sev.black <- sev.black[-which(sev.black$SR ==1),]


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

AIC(mod, mod1, mod2, mod3) ### mod3
anova(mod, mod3, mod2) # sig diff - use mod 3
fricmod_blue <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE richness

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1,mod2, mod3) ## mod and mod 2 similar - choosing mod because simpler 
fevemod_blue <- lme(FEve ~ n_trait, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE evenness

## FDis
mod <- lme(FDis ~ n_trait, random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1,mod2, mod3) #mod best fit
fdismod_blue <- lme(FDis ~ n_trait, random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))


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

AIC(mod, mod1,mod2, mod3) #mod2 best fit
anova(mod, mod2, mod1) # sig diff. using mod 2
raoqmod_blue <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))


summary(fricmod_blue)
summary(fevemod_blue)
summary(fdismod_blue)
summary(fdivmod_blue)
summary(raoqmod_blue)

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

AIC(mod, mod1, mod2, mod3) ### mod3
anova(mod, mod3, mod2)
fricmod_black <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE richness

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1,mod2, mod3) ## models all similar
anova(mod, mod1) # not sig diff, use mod
fevemod_black <- lme(FEve ~ n_trait, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE evenness

## FDis
mod <- lme(FDis ~ n_trait, random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1,mod2, mod3) #mod best fit
fdismod_black <- lme(FDis ~ n_trait, random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE dispersion

# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1,mod2, mod3) #mod2 best fit
anova(mod, mod2, mod1) #mod2 sig diff 
fdivmod_black <- lme(FDiv ~ n_trait + I(n_trait^2)+ I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots),control = lmeControl(opt = "optim"))
mod2 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1,mod2, mod3) #mod2 best fit
anova(mod, mod2, mod1) # sig diff. using mod 2
raoqmod_black <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))


summary(fricmod_black)
summary(fevemod_black)
summary(fdismod_black)
summary(fdivmod_black)
summary(raoqmod_black)

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

AIC(mod, mod1, mod2, mod3) ### mod1,mod2, mod3 very similar AIC values - choosing mod1 because it is simpler
anova(mod1, mod)
fricmod_cdr1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) ## mod1 best fit
anova(mod, mod1) # sig diff - using mod1
fevemod_cdr1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness

## FDis
mod <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod best fit
fdismod_cdr1 <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))


# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod best fit
fdivmod_cdr1 <- lme(FDiv ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod2 best fit
anova(mod, mod2, mod1) # sig diff. using mod 2
raoqmod_cdr1 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod_cdr1)
summary(fevemod_cdr1)
summary(fdismod_cdr1)
summary(fdivmod_cdr1)
summary(raoqmod_cdr1)


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

AIC(mod, mod1, mod2, mod3) ### mod
fricmod_cdr2 <- lme(FRic ~ n_trait, random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) ## mod2 best fit
anova(mod, mod2, mod1) # sig diff. using mod2
fevemod_cdr2 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness

## FDis
mod <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod best fit
fdismod_cdr2 <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))


# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod1 best fit
anova(mod, mod1) # sig diff. using mod 1
fdivmod_cdr2 <- lme(FDiv ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod1 and mod 2 similar - choosing mod1 because simpler
anova(mod, mod1, mod2)#mod1 and mod 2 similar - choosing mod1 because simpler
raoqmod_cdr2 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod_cdr2)
summary(fevemod_cdr2)
summary(fdismod_cdr2)
summary(fdivmod_cdr2)
summary(raoqmod_cdr2)

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

AIC(mod, mod1, mod2, mod3) ### mod2
anova(mod, mod2, mod1) # not sig. diff. using mod
fricmod_cdr3 <- lme(FRic ~ n_trait, random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) ## mod1 best fit
anova(mod, mod1) # sig. diff - using mod 1
fevemod_cdr3 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness

## FDis
mod <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod best fit
fdismod_cdr3 <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))


# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod2 best fit
anova(mod, mod2, mod1) # sig diff. using mod2
fdivmod_cdr3 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod 2 
anova(mod, mod2, mod1) # mod1 not sig diff from mod2, using mod1
raoqmod_cdr3 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod_cdr3)
summary(fevemod_cdr3)
summary(fdismod_cdr3)
summary(fdivmod_cdr3)
summary(raoqmod_cdr3)


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

AIC(mod, mod1, mod2, mod3) ### mod1
anova(mod, mod2, mod1) # mod 1 and mod2 not sig diff - using mod 1
fricmod_cdr4 <- lme(FRic ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) ## mod1  & mod 2 similiar
anova(mod, mod1, mod2) # sig. diff from mod- using mod 1
fevemod_cdr4 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness

## FDis
mod <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod best fit
fdismod_cdr4 <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))


# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod1 best fit
anova(mod, mod1) # sig diff. using mod1
fdivmod_cdr4 <- lme(FDiv ~ n_trait + I(n_trait^2) , random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod1 and mod2 similar 
anova(mod, mod1, mod2) # mod1 not sig diff from mod2, using mod1
raoqmod_cdr4 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod_cdr4)
summary(fevemod_cdr4)
summary(fdismod_cdr4)
summary(fdivmod_cdr4)
summary(raoqmod_cdr4)


################### GRAPHS ##############################

# create combined data set of raw data
cdr.1$community <- "ambient"
cdr.2$community <- "elevated_N"
cdr.3$community <- "elevated_CO2"
cdr.4$community <- "elevated_N_CO2"
cdr.full <- rbind(cdr.1, cdr.2, cdr.3, cdr.4)

cdr.full <- cdr.full %>% group_by(community, n_trait) %>%
  summarize (FRic = mean(FRic), FEve = mean(FEve), FDis = mean(FDis),
             FDiv = mean(FDiv), RaoQ = mean(RaoQ))

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


A <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.1, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.1, lwd = 2, color = "black") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.2, alpha = 0.5, fill = "goldenrod") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.2, lwd = 2, color = "goldenrod") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.3, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.3, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.4, alpha = 0.5, fill = "forestgreen") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.4, lwd = 2, color = "forestgreen") + 
  geom_point(aes(x = n_trait, y = FRic, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue", "goldenrod", "forestgreen"), name = "Community") +
  labs(x = "Number of Traits", y = "FRich") +
  theme_pubr()

B <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.1, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.1, lwd = 2, color = "black") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.2, alpha = 0.5, fill = "goldenrod") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.2, lwd = 2, color = "goldenrod") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.3, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.3, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.4, alpha = 0.5, fill = "forestgreen") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.4, lwd = 2, color = "forestgreen") + 
  geom_point(aes(x = n_trait, y = FEve, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue", "goldenrod", "forestgreen"), name = "Community") +
  labs(x = "Number of Traits", y = "FEve") +
  theme_pubr()

C <-ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.1, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.1, lwd = 2, color = "black") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.2, alpha = 0.5, fill = "goldenrod") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.2, lwd = 2, color = "goldenrod") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.3, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.3, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.4, alpha = 0.5, fill = "forestgreen") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.4, lwd = 2, color = "forestgreen") + 
  geom_point(aes(x = n_trait, y = FDis, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue", "goldenrod", "forestgreen"), name = "Community") +
  labs(x = "Number of Traits", y = "FDis") +
  theme_pubr()

D <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.1, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.1, lwd = 2, color = "black") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.2, alpha = 0.5, fill = "goldenrod") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.2, lwd = 2, color = "goldenrod") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.3, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.3, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.4, alpha = 0.5, fill = "forestgreen") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.4, lwd = 2, color = "forestgreen") + 
  geom_point(aes(x = n_trait, y = FDiv, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue", "goldenrod", "forestgreen"), name = "Community") +
  labs(x = "Number of Traits", y = "FDiv") +
  theme_pubr()

E <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.1, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.1, lwd = 2, color = "black") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.2, alpha = 0.5, fill = "goldenrod") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.2, lwd = 2, color = "goldenrod") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.3, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.3, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.4, alpha = 0.5, fill = "forestgreen") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.4, lwd = 2, color = "forestgreen") + 
  geom_point(aes(x = n_trait, y = RaoQ, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue", "goldenrod", "forestgreen"), name = "Community") +
  labs(x = "Number of Traits", y = "Rao Q") +
  theme_pubr()


ggarrange(plotlist = list(A, B, C, D, E), common.legend = TRUE,
          labels = c("A", "B", "C", "D", "E"))


###SEV####


# create combined data set of raw data
sev.black$community <- "black"
sev.blue$community <- "blue"
sev.full <- rbind(sev.blue[,-7], sev.black[,-7])

sev.full <- sev.full %>% group_by(community, n_trait) %>%
  summarize (FRic = mean(FRic, na.rm = TRUE), FEve = mean(FEve, na.rm = TRUE), 
             FDis = mean(FDis, na.rm = TRUE), FDiv = mean(FDiv, na.rm = TRUE), 
             RaoQ = mean(RaoQ, na.rm = TRUE))

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

A <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FRic, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FRich") +
  theme_pubr()

B <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FEve, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FEve") +
  theme_pubr()

C <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FDis, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FDis") +
  theme_pubr()

D <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FDiv, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FDiv") +
  theme_pubr()

E <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = RaoQ, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "Rao Q") +
  theme_pubr()

png()
ggarrange(plotlist = list(A, B, C, D, E), common.legend = TRUE, 
       labels = c("A", "B", "C", "D", "E"))


