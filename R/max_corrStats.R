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

cdr.1 <- read.csv(here("data/Cleaned/cdr1_euc.csv"))
cdr.1 <- cdr.1[-which(is.na(cdr.1$max_cor)),]
cdr.2 <- read.csv(here("data/Cleaned/cdr2_euc.csv"))
cdr.2 <- cdr.2[-which(is.na(cdr.2$max_cor)),]
cdr.3 <- read.csv(here("data/Cleaned/cdr3_euc.csv"))
cdr.3 <- cdr.3[-which(is.na(cdr.3$max_cor)),]
cdr.4 <- read.csv(here("data/Cleaned/cdr4_euc.csv"))
cdr.4 <- cdr.4[-which(is.na(cdr.4$max_cor)),]
sev.blue <- read.csv(here("data/Cleaned/sevblue_euc.csv"))
sev.blue <- sev.blue[-which(is.na(sev.blue$max_cor)),]
sev.black <- read.csv(here("data/Cleaned/sevblack_euc.csv"))
sev.black <- sev.black[-which(is.na(sev.black$max_cor)),]

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
mod <- lme(FRic ~ max_cor, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FRic ~ max_cor + I(max_cor^2), random = ~1|sev.blueplots, sev.blue[-which(is.na(sev.blue$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FRic ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FRic ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1, mod2, mod3) ### mod
fricmod_blue <- lme(FRic ~ max_cor, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE richness

### FEve
mod <- lme(FEve ~ max_cor, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FEve ~ max_cor + I(max_cor^2), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FEve ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FEve ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1,mod2, mod3) ## mod3
fevemod_blue <- lme(FEve ~ max_cor  + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))

### KDE evenness

## FDis
mod <- lme(FDis ~ max_cor, random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FDis ~ max_cor + I(max_cor^2), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FDis ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FDis ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1,mod2, mod3) #mod2 best fit 
fdismod_blue <- lme(FDis ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))


# FDiv
mod <- lme(FDiv ~ max_cor, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(FDiv ~ max_cor + I(max_cor^2), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod2 <- lme(FDiv ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(FDiv ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1,mod2, mod3) #mod3
fdivmod_blue <- lme(FDiv ~ max_cor  + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

### Rao's Q
mod <- lme(RaoQ ~ max_cor, random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod1 <- lme(RaoQ ~ max_cor + I(max_cor^2), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots),control = lmeControl(opt = "optim"))
mod2 <- lme(RaoQ ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
mod3 <- lme(RaoQ ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|sev.blueplots,  data = sev.blue, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))

AIC(mod, mod1,mod2, mod3) #mod2
raoqmod_blue <- lme(RaoQ ~ max_cor  + I(max_cor^2) + I(max_cor^3), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))


summary(fricmod_blue)
summary(fevemod_blue)
summary(fdismod_blue)
summary(fdivmod_blue)
summary(raoqmod_blue)

###########################################################
################### BLACK ################################
##########################################################


###### Fit Metric ~ max_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ max_cor, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FRic)),], method = "ML",correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FRic ~ max_cor + I(max_cor^2), random = ~1|sev.blackplots, sev.black[-which(is.na(sev.black$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FRic ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FRic ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FRic)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1, mod2, mod3) ### mod3 lowest
fricmod_black <- lme(FRic ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE richness

### FEve
mod <- lme(FEve ~ max_cor, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FEve ~ max_cor + I(max_cor^2), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FEve ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FEve ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1,mod2, mod3) ## all within 2 AIC
anova(mod1, mod, mod2) # mod1 sig diff from mod
fevemod_black <- lme(FEve ~ max_cor + I(max_cor^2), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE evenness

## FDis
mod <- lme(FDis ~ max_cor, random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FDis ~ max_cor + I(max_cor^2), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FDis ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FDis ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1,mod2, mod3) #mod3 best fit
fdismod_black <- lme(FDis ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))

### KDE dispersion

# FDiv
mod <- lme(FDiv ~ max_cor, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(FDiv ~ max_cor + I(max_cor^2), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod2 <- lme(FDiv ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(FDiv ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1,mod2, mod3) #mod3 best fit
fdivmod_black <- lme(FDiv ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

### Rao's Q
mod <- lme(RaoQ ~ max_cor, random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod1 <- lme(RaoQ ~ max_cor + I(max_cor^2), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots),control = lmeControl(opt = "optim"))
mod2 <- lme(RaoQ ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
mod3 <- lme(RaoQ ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|sev.blackplots,  data = sev.black, method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))

AIC(mod, mod1,mod2, mod3) #mod3
raoqmod_black <- lme(RaoQ ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))


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

###### Fit Metric ~ max_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ max_cor, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FRic ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FRic ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod2 
fricmod_cdr1 <- lme(FRic ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness

### FEve
mod <- lme(FEve ~ max_cor, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) ## all within 2 AIC
anova(mod3, mod, mod1) #using mod - not sig diff
fevemod_cdr1 <- lme(FEve ~ max_cor, random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness

## FDis
mod <- lme(FDis ~ max_cor, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod3 best fit
fdismod_cdr1 <- lme(FDis ~ max_cor+ I(max_cor^2) + I(max_cor^3)+ I(max_cor^4), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))


# FDiv
mod <- lme(FDiv ~ max_cor, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod3 best fit
fdivmod_cdr1 <- lme(FDiv ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ max_cor, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod3 
raoqmod_cdr1 <- lme(RaoQ ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod_cdr1)
summary(fevemod_cdr1)
summary(fdismod_cdr1)
summary(fdivmod_cdr1)
summary(raoqmod_cdr1)


###########################################################
################### CDR 2 ################################
##########################################################
###### Fit Metric ~ max_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ max_cor, random = ~1|Plot,  data = cdr.2, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FRic ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FRic ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod3
fricmod_cdr2 <- lme(FRic ~ max_cor + I(max_cor^2)+ I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness

### FEve
mod <- lme(FEve ~ max_cor, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) ## mod3 best fit
fevemod_cdr2 <- lme(FEve ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness

## FDis
mod <- lme(FDis ~ max_cor, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod3 best fit
fdismod_cdr2 <- lme(FDis ~ max_cor+ I(max_cor^2)+ I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))


# FDiv
mod <- lme(FDiv ~ max_cor, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod
fdivmod_cdr2 <- lme(FDiv ~ max_cor, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ max_cor, random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod3
raoqmod_cdr2 <- lme(RaoQ ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod_cdr2)
summary(fevemod_cdr2)
summary(fdismod_cdr2)
summary(fdivmod_cdr2)
summary(raoqmod_cdr2)

###########################################################
################### CDR 3 ################################
##########################################################
###### Fit Metric ~ max_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ max_cor, random = ~1|Plot,  data = cdr.3, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FRic ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FRic ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod3
fricmod_cdr3 <- lme(FRic ~ max_cor+ I(max_cor^2)+ I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness

### FEve
mod <- lme(FEve ~ max_cor, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) ## mod2 best fit
fevemod_cdr3 <- lme(FEve ~ max_cor + I(max_cor^2)  + I(max_cor^3), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness

## FDis
mod <- lme(FDis ~ max_cor, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod3 best fit
fdismod_cdr3 <- lme(FDis ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))


# FDiv
mod <- lme(FDiv ~ max_cor, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod1 - mod3 within 2 AIC best fit
anova(mod1, mod2, mod3) #using mod1 - no sig differences between models
fdivmod_cdr3 <- lme(FDiv ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ max_cor, random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod3
raoqmod_cdr3 <- lme(RaoQ ~ max_cor+ I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod_cdr3)
summary(fevemod_cdr3)
summary(fdismod_cdr3)
summary(fdivmod_cdr3)
summary(raoqmod_cdr3)

###########################################################
################### CDR 4 ################################
##########################################################
###### Fit Metric ~ max_cor models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ max_cor, random = ~1|Plot,  data = cdr.4, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FRic ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FRic ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1, mod2, mod3) ### mod3
fricmod_cdr4 <- lme(FRic ~ max_cor+ I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness

### FEve
mod <- lme(FEve ~ max_cor, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) ## mod 
fevemod_cdr4 <- lme(FEve ~ max_cor, random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness

## FDis
mod <- lme(FDis ~ max_cor, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod3
fdismod_cdr4 <- lme(FDis ~ max_cor+ I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))


# FDiv
mod <- lme(FDiv ~ max_cor, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod3 similar
fdivmod_cdr4 <- lme(FDiv ~ max_cor + I(max_cor^2)  + I(max_cor^3)+ I(max_cor^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ max_cor, random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ max_cor + I(max_cor^2), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ max_cor + I(max_cor^2) + I(max_cor^3), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod3
raoqmod_cdr4 <- lme(RaoQ ~ max_cor + I(max_cor^2) + I(max_cor^3) + I(max_cor^4), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))


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

cdr.full <- cdr.full %>% group_by(community, max_cor) %>%
  summarize (FRic = mean(FRic), FEve = mean(FEve), FDis = mean(FDis),
             FDiv = mean(FDiv), RaoQ = mean(RaoQ))

new.df <- data.frame(max_cor = seq(-.02, .92, by = 0.01))


##############################################
############## CDR 1 ########################
#############################################
fr.1 <- as.data.frame(predictSE.lme(fricmod_cdr1, new.df))
fr.1$max_cor <- seq(-.02, .92, by = 0.01)
fr.1$lwr <- fr.1$fit - fr.1$se.fit
fr.1$upr <- fr.1$fit + fr.1$se.fit
fe.1 <- as.data.frame(predictSE.lme(fevemod_cdr1, new.df))
fe.1$max_cor <- seq(-.02, .92, by = 0.01)
fe.1$lwr <- fe.1$fit - fe.1$se.fit
fe.1$upr <- fe.1$fit + fe.1$se.fit
fdis.1 <- as.data.frame(predictSE.lme(fdismod_cdr1, new.df))
fdis.1$max_cor <- seq(-.02, .92, by = 0.01)
fdis.1$lwr <- fdis.1$fit - fdis.1$se.fit
fdis.1$upr <- fdis.1$fit + fdis.1$se.fit
fdiv.1 <- as.data.frame(predictSE.lme(fdivmod_cdr1, new.df))
fdiv.1$max_cor <- seq(-.02, .92, by = 0.01)
fdiv.1$lwr <- fdiv.1$fit - fdiv.1$se.fit
fdiv.1$upr <- fdiv.1$fit + fdiv.1$se.fit
rq.1 <- as.data.frame(predictSE.lme(raoqmod_cdr1, new.df))
rq.1$max_cor <- seq(-.02, .92, by = 0.01)
rq.1$lwr <- rq.1$fit - rq.1$se.fit
rq.1$upr <- rq.1$fit + rq.1$se.fit

fr.2 <- as.data.frame(predictSE.lme(fricmod_cdr2, new.df))
fr.2$max_cor <- seq(-.02, .92, by = 0.01)
fr.2$lwr <- fr.2$fit - fr.2$se.fit
fr.2$upr <- fr.2$fit + fr.2$se.fit
fe.2 <- as.data.frame(predictSE.lme(fevemod_cdr2, new.df))
fe.2$max_cor <- seq(-.02, .92, by = 0.01)
fe.2$lwr <- fe.2$fit - fe.2$se.fit
fe.2$upr <- fe.2$fit + fe.2$se.fit
fdis.2 <- as.data.frame(predictSE.lme(fdismod_cdr2, new.df))
fdis.2$max_cor <- seq(-.02, .92, by = 0.01)
fdis.2$lwr <- fdis.2$fit - fdis.2$se.fit
fdis.2$upr <- fdis.2$fit + fdis.2$se.fit
fdiv.2 <- as.data.frame(predictSE.lme(fdivmod_cdr2, new.df))
fdiv.2$max_cor <- seq(-.02, .92, by = 0.01)
fdiv.2$lwr <- fdiv.2$fit - fdiv.2$se.fit
fdiv.2$upr <- fdiv.2$fit + fdiv.2$se.fit
rq.2 <- as.data.frame(predictSE.lme(raoqmod_cdr2, new.df))
rq.2$max_cor <- seq(-.02, .92, by = 0.01)
rq.2$lwr <- rq.2$fit - rq.2$se.fit
rq.2$upr <- rq.2$fit + rq.2$se.fit

fr.3 <- as.data.frame(predictSE.lme(fricmod_cdr3, new.df))
fr.3$max_cor <- seq(-.02, .92, by = 0.01)
fr.3$lwr <- fr.3$fit - fr.3$se.fit
fr.3$upr <- fr.3$fit + fr.3$se.fit
fe.3 <- as.data.frame(predictSE.lme(fevemod_cdr3, new.df))
fe.3$max_cor <- seq(-.02, .92, by = 0.01)
fe.3$lwr <- fe.3$fit - fe.3$se.fit
fe.3$upr <- fe.3$fit + fe.3$se.fit
fdis.3 <- as.data.frame(predictSE.lme(fdismod_cdr3, new.df))
fdis.3$max_cor <- seq(-.02, .92, by = 0.01)
fdis.3$lwr <- fdis.3$fit - fdis.3$se.fit
fdis.3$upr <- fdis.3$fit + fdis.3$se.fit
fdiv.3 <- as.data.frame(predictSE.lme(fdivmod_cdr3, new.df))
fdiv.3$max_cor <- seq(-.02, .92, by = 0.01)
fdiv.3$lwr <- fdiv.3$fit - fdiv.3$se.fit
fdiv.3$upr <- fdiv.3$fit + fdiv.3$se.fit
rq.3 <- as.data.frame(predictSE.lme(raoqmod_cdr3, new.df))
rq.3$max_cor <- seq(-.02, .92, by = 0.01)
rq.3$lwr <- rq.3$fit - rq.3$se.fit
rq.3$upr <- rq.3$fit + rq.3$se.fit

fr.4 <- as.data.frame(predictSE.lme(fricmod_cdr4, new.df))
fr.4$max_cor <- seq(-.02, .92, by = 0.01)
fr.4$lwr <- fr.4$fit - fr.4$se.fit
fr.4$upr <- fr.4$fit + fr.4$se.fit
fe.4 <- as.data.frame(predictSE.lme(fevemod_cdr4, new.df))
fe.4$max_cor <- seq(-.02, .92, by = 0.01)
fe.4$lwr <- fe.4$fit - fe.4$se.fit
fe.4$upr <- fe.4$fit + fe.4$se.fit
fdis.4 <- as.data.frame(predictSE.lme(fdismod_cdr4, new.df))
fdis.4$max_cor <- seq(-.02, .92, by = 0.01)
fdis.4$lwr <- fdis.4$fit - fdis.4$se.fit
fdis.4$upr <- fdis.4$fit + fdis.4$se.fit
fdiv.4 <- as.data.frame(predictSE.lme(fdivmod_cdr4, new.df))
fdiv.4$max_cor <- seq(-.02, .92, by = 0.01)
fdiv.4$lwr <- fdiv.4$fit - fdiv.4$se.fit
fdiv.4$upr <- fdiv.4$fit + fdiv.4$se.fit
rq.4 <- as.data.frame(predictSE.lme(raoqmod_cdr4, new.df))
rq.4$max_cor <- seq(-.02, .92, by = 0.01)
rq.4$lwr <- rq.4$fit - rq.4$se.fit
rq.4$upr <- rq.4$fit + rq.4$se.fit


A <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fr.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= max_cor, y = fit), data = fr.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fr.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= max_cor, y = fit), data = fr.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fr.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= max_cor, y = fit), data = fr.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fr.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= max_cor, y = fit), data = fr.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = max_cor, y = FRic, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "FRich") +
  theme_pubr()

B <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fe.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= max_cor, y = fit), data = fe.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fe.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= max_cor, y = fit), data = fe.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fe.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= max_cor, y = fit), data = fe.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fe.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= max_cor, y = fit), data = fe.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = max_cor, y = FEve, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "FEve") +
  theme_pubr()

C <-ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fdis.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= max_cor, y = fit), data = fdis.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fdis.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= max_cor, y = fit), data = fdis.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fdis.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= max_cor, y = fit), data = fdis.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fdis.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= max_cor, y = fit), data = fdis.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = max_cor, y = FDis, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "FDis") +
  theme_pubr()

D <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fdiv.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= max_cor, y = fit), data = fdiv.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fdiv.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= max_cor, y = fit), data = fdiv.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fdiv.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= max_cor, y = fit), data = fdiv.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fdiv.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= max_cor, y = fit), data = fdiv.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = max_cor, y = FDiv, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "FDiv") +
  theme_pubr()

E <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = rq.1, alpha = 0.5, fill = "#9D8F0F") + 
  geom_line(aes(x= max_cor, y = fit), data = rq.1, lwd = 2, color = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = rq.2, alpha = 0.5, fill = "#731279") + 
  geom_line(aes(x= max_cor, y = fit), data = rq.2, lwd = 2, color = "#731279") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = rq.3, alpha = 0.5, fill = "#00B7FF") + 
  geom_line(aes(x= max_cor, y = fit), data = rq.3, lwd = 2, color = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = rq.4, alpha = 0.5, fill = "#075A13") + 
  geom_line(aes(x= max_cor, y = fit), data = rq.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = max_cor, y = RaoQ, color = community), data = cdr.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F", "#00B7FF", "#731279", "#075A13"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "Rao Q") +
  theme_pubr()

png(here("Figures/cdr_maxcor.png"), height = 5, width = 9, units = 'in', res = 300)
ggarrange(plotlist = list(A, B, C, D, E), common.legend = TRUE,
          labels = c("A", "B", "C", "D", "E"))
dev.off()


###SEV####


# create combined data set of raw data
sev.black$community <- "black"
sev.blue$community <- "blue"
sev.full <- rbind(sev.blue[,-10], sev.black[,-10])

sev.full <- sev.full %>% group_by(community, max_cor) %>%
  summarize (FRic = mean(FRic, na.rm = TRUE), FEve = mean(FEve, na.rm = TRUE), 
             FDis = mean(FDis, na.rm = TRUE), FDiv = mean(FDiv, na.rm = TRUE), 
             RaoQ = mean(RaoQ, na.rm = TRUE))

new.df <- data.frame(max_cor = seq(-0.03, .6, by = 0.01))

fr.blue <- as.data.frame(predictSE.lme(fricmod_blue, new.df))
fr.blue$max_cor <- seq(-0.03, .6, by = 0.01)
fr.blue$lwr <- fr.blue$fit - fr.blue$se.fit
fr.blue$upr <- fr.blue$fit + fr.blue$se.fit
fe.blue <- as.data.frame(predictSE.lme(fevemod_blue, new.df))
fe.blue$max_cor <- seq(-0.03, .6, by = 0.01)
fe.blue$lwr <- fe.blue$fit - fe.blue$se.fit
fe.blue$upr <- fe.blue$fit + fe.blue$se.fit
fdis.blue <- as.data.frame(predictSE.lme(fdismod_blue, new.df))
fdis.blue$max_cor <- seq(-0.03, .6, by = 0.01)
fdis.blue$lwr <- fdis.blue$fit - fdis.blue$se.fit
fdis.blue$upr <- fdis.blue$fit + fdis.blue$se.fit
fdiv.blue <- as.data.frame(predictSE.lme(fdivmod_blue, new.df))
fdiv.blue$max_cor <- seq(-0.03, .6, by = 0.01)
fdiv.blue$lwr <- fdiv.blue$fit - fdiv.blue$se.fit
fdiv.blue$upr <- fdiv.blue$fit + fdiv.blue$se.fit
rq.blue <- as.data.frame(predictSE.lme(raoqmod_blue, new.df))
rq.blue$max_cor <- seq(-0.03, .6, by = 0.01)
rq.blue$lwr <- rq.blue$fit - rq.blue$se.fit
rq.blue$upr <- rq.blue$fit + rq.blue$se.fit

fr.black <- as.data.frame(predictSE.lme(fricmod_black, new.df))
fr.black$max_cor <- seq(-0.03, .6, by = 0.01)
fr.black$lwr <- fr.black$fit - fr.black$se.fit
fr.black$upr <- fr.black$fit + fr.black$se.fit
fe.black <- as.data.frame(predictSE.lme(fevemod_black, new.df))
fe.black$max_cor <- seq(-0.03, .6, by = 0.01)
fe.black$lwr <- fe.black$fit - fe.black$se.fit
fe.black$upr <- fe.black$fit + fe.black$se.fit
fdis.black <- as.data.frame(predictSE.lme(fdismod_black, new.df))
fdis.black$max_cor <- seq(-0.03, .6, by = 0.01)
fdis.black$lwr <- fdis.black$fit - fdis.black$se.fit
fdis.black$upr <- fdis.black$fit + fdis.black$se.fit
fdiv.black <- as.data.frame(predictSE.lme(fdivmod_black, new.df))
fdiv.black$max_cor <- seq(-0.03, .6, by = 0.01)
fdiv.black$lwr <- fdiv.black$fit - fdiv.black$se.fit
fdiv.black$upr <- fdiv.black$fit + fdiv.black$se.fit
rq.black <- as.data.frame(predictSE.lme(raoqmod_black, new.df))
rq.black$max_cor <- seq(-0.03, .6, by = 0.01)
rq.black$lwr <- rq.black$fit - rq.black$se.fit
rq.black$upr <- rq.black$fit + rq.black$se.fit

A <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fr.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= max_cor, y = fit), data = fr.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fr.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= max_cor, y = fit), data = fr.black, lwd = 2, color = "black") + 
  geom_point(aes(x = max_cor, y = FRic, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#000000", "navyblue"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "FRich") +
  theme_pubr()

B <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fe.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= max_cor, y = fit), data = fe.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fe.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= max_cor, y = fit), data = fe.black, lwd = 2, color = "black") + 
  geom_point(aes(x = max_cor, y = FEve, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#000000", "navyblue"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "FEve") +
  theme_pubr()

C <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fdis.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= max_cor, y = fit), data = fdis.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fdis.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= max_cor, y = fit), data = fdis.black, lwd = 2, color = "black") + 
  geom_point(aes(x = max_cor, y = FDis, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "FDis") +
  theme_pubr()

D <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fdiv.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= max_cor, y = fit), data = fdiv.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = fdiv.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= max_cor, y = fit), data = fdiv.black, lwd = 2, color = "black") + 
  geom_point(aes(x = max_cor, y = FDiv, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "FDiv") +
  theme_pubr()

E <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = rq.blue, alpha = 0.5, fill = "navyblue") + 
  geom_line(aes(x= max_cor, y = fit), data = rq.blue, lwd = 2, color = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = max_cor), data = rq.black, alpha = 0.5, fill = "black") + 
  geom_line(aes(x= max_cor, y = fit), data = rq.black, lwd = 2, color = "black") + 
  geom_point(aes(x = max_cor, y = RaoQ, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("black", "navyblue"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "Rao Q") +
  theme_pubr()

png(here("Figures/Sev_maxcorr.png"), height = 5, width = 9, units = 'in', res = 300)
ggarrange(plotlist = list(A, B, C, D, E), common.legend = TRUE, 
          labels = c("A", "B", "C", "D", "E"))
dev.off()



