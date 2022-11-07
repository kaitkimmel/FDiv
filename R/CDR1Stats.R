################################
#### STATISTICAL ANALYSES   ####
################################

# CREATED BY: KAITLIN KIMMEL

# load libraries
library(nlme)
library(AICcmodavg)


# load data

cdr.1 <- read.csv(here("data/Cleaned/cdr1.csv"))



###### Fit Metric ~ N_trait models
# Plot will be fit as a random effect
# Testing for best functional form between linear, quadratic, cubic, and quartic fits
### FRic
mod <- lme(FRic ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML",correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

anova(mod, mod1)
AIC(mod, mod1, mod2, mod3) ### mod1 and mod2 very similar AIC values - choosing mod1 because it is simpler
fricmod <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE richness

### FEve
mod <- lme(FEve ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) ## mod1 best fit
fevemod <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

### KDE evenness

## FDis
mod <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod 2 best fit
fdismod <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))


# FDiv
mod <- lme(FDiv ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod 1 best fit
fdivmod <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

### Rao's Q
mod <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod2 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
mod3 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))

AIC(mod, mod1,mod2, mod3) #mod best fit
raoqmod <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))


summary(fricmod)
summary(fevemod)
summary(fdismod)
summary(fdivmod)
summary(raoqmod)

