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


# load data

cdr.1 <- read.csv(here("data/Cleaned/cdr1.csv"))
cdr.2 <- read.csv(here("data/Cleaned/cdr2.csv"))
cdr.3 <- read.csv(here("data/Cleaned/cdr3.csv"))
cdr.4 <- read.csv(here("data/Cleaned/cdr4.csv"))


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
anova(mod, mod1) # not sig different - using mod
fevemod_cdr1 <- lme(FEve ~ n_trait, random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

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

save(fricmod_cdr1, file = here("Models_out/cdr1_fric.rda"))
save(fevemod_cdr1, file = here("Models_out/cdr1_feve.rda"))
save(fdismod_cdr1, file = here("Models_out/cdr1_fdis.rda"))
save(fdivmod_cdr1, file = here("Models_out/cdr1_fdiv.rda"))
save(raoqmod_cdr1, file = here("Models_out/cdr1_raoq.rda"))


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
anova(mod, mod2, mod1) # sign diff. using mod2
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
anova(mod, mod1) # sign diff. using mod 1
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

save(fricmod_cdr2, file = here("Models_out/cdr2_fric.rda"))
save(fevemod_cdr2, file = here("Models_out/cdr2_feve.rda"))
save(fdismod_cdr2, file = here("Models_out/cdr2_fdis.rda"))
save(fdivmod_cdr2, file = here("Models_out/cdr2_fdiv.rda"))
save(raoqmod_cdr2, file = here("Models_out/cdr2_raoq.rda"))

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

save(fricmod_cdr3, file = here("Models_out/cdr3_fric.rda"))
save(fevemod_cdr3, file = here("Models_out/cdr3_feve.rda"))
save(fdismod_cdr3, file = here("Models_out/cdr3_fdis.rda"))
save(fdivmod_cdr3, file = here("Models_out/cdr3_fdiv.rda"))
save(raoqmod_cdr3, file = here("Models_out/cdr3_raoq.rda"))


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

save(fricmod_cdr4, file = here("Models_out/cdr4_fric.rda"))
save(fevemod_cdr4, file = here("Models_out/cdr4_feve.rda"))
save(fdismod_cdr4, file = here("Models_out/cdr4_fdis.rda"))
save(fdivmod_cdr4, file = here("Models_out/cdr4_fdiv.rda"))
save(raoqmod_cdr4, file = here("Models_out/cdr4_raoq.rda"))


################### GRAPHS ##############################
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


A <- ggplot() + 
  geom_point(aes(x = n_trait, y = FRic, color = as.factor(SR)), data = cdr.1) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.1, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = fr.1, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

B <- ggplot() + 
  geom_point(aes(x = n_trait, y = FEve, color = as.factor(SR)), data = cdr.1) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.1, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = fe.1, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

C <- ggplot() + 
  geom_point(aes(x = n_trait, y = FDis, color = as.factor(SR)), data = cdr.1) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.1, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.1, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

D <- ggplot() + 
  geom_point(aes(x = n_trait, y = FDiv, color = as.factor(SR)), data = cdr.1) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.1, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.1, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

E <- ggplot() + 
  geom_point(aes(x = n_trait, y = RaoQ, color = as.factor(SR)), data = cdr.1) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.1, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = rq.1, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

######################################
########## CDR 2 ######################
######################################

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


F <- ggplot() + 
  geom_point(aes(x = n_trait, y = FRic, color = as.factor(SR)), data = cdr.2) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.2, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = fr.2, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

G <- ggplot() + 
  geom_point(aes(x = n_trait, y = FEve, color = as.factor(SR)), data = cdr.2) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.2, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = fe.2, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

H <- ggplot() + 
  geom_point(aes(x = n_trait, y = FDis, color = as.factor(SR)), data = cdr.2) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.2, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.2, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

I <- ggplot() + 
  geom_point(aes(x = n_trait, y = FDiv, color = as.factor(SR)), data = cdr.2) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.2, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.2, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

J <- ggplot() + 
  geom_point(aes(x = n_trait, y = RaoQ, color = as.factor(SR)), data = cdr.2) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.2, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = rq.2, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

######################################
########## CDR 3 ######################
######################################

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


K <- ggplot() + 
  geom_point(aes(x = n_trait, y = FRic, color = as.factor(SR)), data = cdr.3) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.3, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = fr.3, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68", "#015B69"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

L <- ggplot() + 
  geom_point(aes(x = n_trait, y = FEve, color = as.factor(SR)), data = cdr.3) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.3, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = fe.3, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68", "#015B69"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

M <- ggplot() + 
  geom_point(aes(x = n_trait, y = FDis, color = as.factor(SR)), data = cdr.3) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.3, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.3, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68", "#015B69"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

N <- ggplot() + 
  geom_point(aes(x = n_trait, y = FDiv, color = as.factor(SR)), data = cdr.3) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.3, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.3, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68", "#015B69"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

O <- ggplot() + 
  geom_point(aes(x = n_trait, y = RaoQ, color = as.factor(SR)), data = cdr.3) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.3, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = rq.3, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68", "#015B69"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()


######################################
########## CDR 4 ######################
######################################

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


P <- ggplot() + 
  geom_point(aes(x = n_trait, y = FRic, color = as.factor(SR)), data = cdr.4) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.4, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = fr.4, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68", "#015B69"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

Q <- ggplot() + 
  geom_point(aes(x = n_trait, y = FEve, color = as.factor(SR)), data = cdr.4) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.4, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = fe.4, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68", "#015B69"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

R <- ggplot() + 
  geom_point(aes(x = n_trait, y = FDis, color = as.factor(SR)), data = cdr.4) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.4, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.4, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68", "#015B69"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

S <- ggplot() + 
  geom_point(aes(x = n_trait, y = FDiv, color = as.factor(SR)), data = cdr.4) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.4, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.4, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68", "#015B69"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

T <- ggplot() + 
  geom_point(aes(x = n_trait, y = RaoQ, color = as.factor(SR)), data = cdr.4) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.4, alpha = 0.5) + 
  geom_line(aes(x= n_trait, y = fit), data = rq.4, lwd = 2) + 
  scale_color_manual(values = c("#9D8F0F", "#5D8B33", "#197F54", "#006F68", "#015B69"), name = "Species\n Richness") +
  labs(x = "Number of Traits") +
  theme_pubr()

ggarrange(plotlist = list(A, B, C, D, E, F, G, H, I , J, K, L, M, N, O, P, Q, R, S, T),
          ncol = 5, nrow = 4, common.legend = TRUE,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I" , "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T"))
ggarrange
