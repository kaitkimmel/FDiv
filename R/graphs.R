######## GRAPHS ##########

# libraries
library(ggplot2)
library(here)
library(nlme)
library(AICcmodavg)
library(dplyr)
library(ggpubr)

# data
# gower
cdr.1 <- read.csv(here("data/Cleaned/cdr1.csv"))
cdr.2 <- read.csv(here("data/Cleaned/cdr2.csv"))
cdr.3 <- read.csv(here("data/Cleaned/cdr3.csv"))
cdr.4 <- read.csv(here("data/Cleaned/cdr4.csv"))
sev.blue <- read.csv(here("data/Cleaned/sevblue.csv"))
sev.black <- read.csv(here("data/Cleaned/sevblack.csv"))

# euclidean
cdre.1 <- read.csv(here("data/Cleaned/cdr1_euc.csv"))
cdre.2 <- read.csv(here("data/Cleaned/cdr2_euc.csv"))
cdre.3 <- read.csv(here("data/Cleaned/cdr3_euc.csv"))
cdre.4 <- read.csv(here("data/Cleaned/cdr4_euc.csv"))
sev.eblue <- read.csv(here("data/Cleaned/sevblue_euc.csv"))
sev.eblack <- read.csv(here("data/Cleaned/sevblack_euc.csv"))

######################
####     Fig. 1   ####
## Number of traits ##
######################

# Pull in models
####### GOWER #############
# SEV BLUE
fricmod_blue <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))
kde.alphamod_blue <- lme(kde.alpha ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))
fevemod_blue <- lme(FEve ~ n_trait, random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))
kde.evennessmod_blue <- lme(kde.evenness ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))
fdismod_blue <- lme(FDis ~ n_trait, random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))
kde.dispersionmod_blue <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))
fdivmod_blue <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.blue[-which(is.na(sev.blue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
raoqmod_blue <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blueplots,  data = sev.blue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

# SEV BLACK
fricmod_black <- lme(FRic ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))
kde.alphamod_black <- lme(kde.alpha ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))
fevemod_black <- lme(FEve ~ n_trait, random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))
kde.evennessmod_black <- lme(kde.evenness ~ n_trait, random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))
fdismod_black <- lme(FDis ~ n_trait, random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))
kde.dispersionmod_black <- lme(kde.dispersion ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))
fdivmod_black <- lme(FDiv ~ n_trait + I(n_trait^2)+ I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black[-which(is.na(sev.black$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
raoqmod_black <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.black, correlation = corCompSymm(form = ~ 1|sev.blackplots))

# CDR1
fricmod_cdr1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))
kde.alphamod_cdr1 <- lme(kde.alpha ~ n_trait, random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))
fevemod_cdr1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))
kde.evennessmod_cdr1 <- lme(kde.evenness ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))
fdismod_cdr1 <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))
kde.dispersionmod_cdr1 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))
fdivmod_cdr1 <- lme(FDiv ~ n_trait, random = ~1|Plot,  data = cdr.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
raoqmod_cdr1 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.1, correlation = corCompSymm(form = ~ 1|Plot))

# CDR2
fricmod_cdr2 <- lme(FRic ~ n_trait, random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))
kde.alphamod_cdr2 <- lme(kde.alpha ~ n_trait, random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))
fevemod_cdr2 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))
kde.evennessmod_cdr2 <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))
fdismod_cdr2 <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))
kde.dispersionmod_cdr2 <- lme(kde.dispersion ~ n_trait, random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))
fdivmod_cdr2 <- lme(FDiv ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdr.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
raoqmod_cdr2 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.2, correlation = corCompSymm(form = ~ 1|Plot))

# CDR3
fricmod_cdr3 <- lme(FRic ~ n_trait, random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))
kde.alphamod_cdr3 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))
fevemod_cdr3 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))
kde.evennessmod_cdr3 <- lme(kde.evenness ~ n_trait, random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))
fdismod_cdr3 <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))
kde.dispersionmod_cdr3 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))
fdivmod_cdr3 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdr.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
raoqmod_cdr3 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.3, correlation = corCompSymm(form = ~ 1|Plot))

# CDR4
fricmod_cdr4 <- lme(FRic ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))
kde.alphamod_cdr4 <- lme(kde.alpha ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))
fevemod_cdr4 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))
kde.evennessmod_cdr4 <- lme(kde.evenness ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))
fdismod_cdr4 <- lme(FDis ~ n_trait, random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))
kde.dispersionmod_cdr4 <- lme(kde.dispersion ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))
fdivmod_cdr4 <- lme(FDiv ~ n_trait + I(n_trait^2) , random = ~1|Plot,  data = cdr.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
raoqmod_cdr4 <- lme(RaoQ ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdr.4, correlation = corCompSymm(form = ~ 1|Plot))

####### EUCLIDEAN #############

# BLUE
fricmod_eblue <- lme(FRic ~ n_trait, random = ~1|sev.blueplots,  data = sev.eblue[-which(is.na(sev.eblue$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))
kde.alphamod_eblue <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))
fevemod_eblue <- lme(FEve ~ n_trait, random = ~1|sev.blueplots,  data = sev.eblue[-which(is.na(sev.eblue$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blueplots))
kde.evennessmod_eblue <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))
fdismod_eblue <- lme(FDis ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))
kde.dispersionmod_eblue <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))
fdivmod_eblue <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|sev.blueplots,  data = sev.eblue[-which(is.na(sev.eblue$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blueplots))
raoqmod_eblue <- lme(RaoQ ~ n_trait, random = ~1|sev.blueplots,  data = sev.eblue, correlation = corCompSymm(form = ~ 1|sev.blueplots))

# BLACK
fricmod_eblack <- lme(FRic ~ n_trait, random = ~1|sev.blackplots,  data = sev.eblack[-which(is.na(sev.eblack$FRic)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))
kde.alphamod_eblack <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))
fevemod_eblack <- lme(FEve ~ n_trait, random = ~1|sev.blackplots,  data = sev.eblack[-which(is.na(sev.eblack$FEve)),], correlation = corCompSymm(form = ~ 1|sev.blackplots))
kde.evennessmod_eblack <- lme(kde.evenness ~ n_trait, random = ~1|sev.blackplots,  data = sev.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))
fdismod_eblack <- lme(FDis ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|sev.blackplots,  data = sev.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))
kde.dispersionmod_eblack <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))
fdivmod_eblack <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|sev.blackplots,  data = sev.eblack[-which(is.na(sev.eblack$FDiv)),], method = "ML", correlation = corCompSymm(form = ~ 1|sev.blackplots))
raoqmod_eblack <- lme(RaoQ ~ n_trait, random = ~1|sev.blackplots,  data = sev.eblack, correlation = corCompSymm(form = ~ 1|sev.blackplots))

# CDR1
fricmod_cdre1 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))
kde.alphamod_cdre1 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))
fevemod_cdre1 <- lme(FEve ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))
kde.evennessmod_cdre1 <- lme(kde.evenness ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))
fdismod_cdre1 <- lme(FDis ~ n_trait+ I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))
kde.dispersionmod_cdre1 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))
fdivmod_cdre1 <- lme(FDiv ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.1, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
raoqmod_cdre1 <- lme(RaoQ ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdre.1, correlation = corCompSymm(form = ~ 1|Plot))

# CDR2
fricmod_cdre2 <- lme(FRic ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))
kde.alphamod_cdre2 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))
fevemod_cdre2 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))
kde.evennessmod_cdre2 <- lme(kde.evenness ~ n_trait, random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))
fdismod_cdre2 <- lme(FDis ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))
kde.dispersionmod_cdre2 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))
fdivmod_cdre2 <- lme(FDiv ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdre.2, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
raoqmod_cdre2 <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdre.2, correlation = corCompSymm(form = ~ 1|Plot))

#CDR3
fricmod_cdre3 <- lme(FRic ~ n_trait+ I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))
kde.alphamod_cdre3 <- lme(kde.alpha ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))
fevemod_cdre3 <- lme(FEve ~ n_trait+ I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))
kde.evennessmod_cdre3 <- lme(kde.evenness ~ n_trait+ I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))
fdismod_cdre3 <- lme(FDis ~ n_trait + I(n_trait^2)+ I(n_trait^3), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))
kde.dispersionmod_cdre3 <- lme(kde.dispersion ~ n_trait + I(n_trait^2), random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))
fdivmod_cdre3 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdre.3, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
raoqmod_cdre3 <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdre.3, correlation = corCompSymm(form = ~ 1|Plot))

# CDR4
fricmod_cdre4 <- lme(FRic ~ n_trait, random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))
kde.alphamod_cdre4 <- lme(kde.alpha ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))
fevemod_cdre4 <- lme(FEve ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))
kde.evennessmod_cdre4 <- lme(kde.evenness ~ n_trait + I(n_trait^2) + I(n_trait^3) + I(n_trait^4), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))
fdismod_cdre4 <- lme(FDis ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))
kde.dispersionmod_cdre4 <- lme(kde.dispersion ~ n_trait+ I(n_trait^2), random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))
fdivmod_cdre4 <- lme(FDiv ~ n_trait + I(n_trait^2) + I(n_trait^3), random = ~1|Plot,  data = cdre.4, method = "ML", correlation = corCompSymm(form = ~ 1|Plot))
raoqmod_cdre4 <- lme(RaoQ ~ n_trait, random = ~1|Plot,  data = cdre.4, correlation = corCompSymm(form = ~ 1|Plot))



####### Predicted values #######
######### GOWER ################
# create combined data set of raw data
cdr.1$community <- "CDR1"
cdr.2$community <- "CDR2"
cdr.3$community <- "CDR3"
cdr.4$community <- "CDR4"
cdr.full <- rbind(cdr.1, cdr.2, cdr.3, cdr.4)

cdr.full <- cdr.full %>% group_by(community, n_trait) %>%
  summarize (FRic = mean(FRic), FEve = mean(FEve), FDis = mean(FDis),
             FDiv = mean(FDiv), RaoQ = mean(RaoQ), kde.alpha = mean(kde.alpha), kde.evenness = mean(kde.evenness), kde.dispersion = mean(kde.dispersion))

new.df <- data.frame(n_trait = seq(2, 9, by = 0.1))

#### CDR ####
fr.1 <- as.data.frame(predictSE.lme(fricmod_cdr1, new.df))
fr.1$n_trait <- seq(2, 9, by = 0.1)
fr.1$lwr <- fr.1$fit - fr.1$se.fit
fr.1$upr <- fr.1$fit + fr.1$se.fit
fe.1 <- as.data.frame(predictSE.lme(fevemod_cdr1, new.df))
fe.1$n_trait <- seq(2, 9, by = 0.1)
fe.1$lwr <- fe.1$fit - fe.1$se.fit
fe.1$upr <- fe.1$fit + fe.1$se.fit
fdis.1 <- as.data.frame(predictSE.lme(fdismod_cdr1, new.df))
fdis.1$n_trait <- seq(2, 9, by = 0.1)
fdis.1$lwr <- fdis.1$fit - fdis.1$se.fit
fdis.1$upr <- fdis.1$fit + fdis.1$se.fit
fdiv.1 <- as.data.frame(predictSE.lme(fdivmod_cdr1, new.df))
fdiv.1$n_trait <- seq(2, 9, by = 0.1)
fdiv.1$lwr <- fdiv.1$fit - fdiv.1$se.fit
fdiv.1$upr <- fdiv.1$fit + fdiv.1$se.fit
rq.1 <- as.data.frame(predictSE.lme(raoqmod_cdr1, new.df))
rq.1$n_trait <- seq(2, 9, by = 0.1)
rq.1$lwr <- rq.1$fit - rq.1$se.fit
rq.1$upr <- rq.1$fit + rq.1$se.fit
kde.alpha.1 <- as.data.frame(predictSE.lme(kde.alphamod_cdr1, new.df))
kde.alpha.1$n_trait <- seq(2, 9, by = 0.1)
kde.alpha.1$lwr <- kde.alpha.1$fit - kde.alpha.1$se.fit
kde.alpha.1$upr <- kde.alpha.1$fit + kde.alpha.1$se.fit
kde.evenness.1 <- as.data.frame(predictSE.lme(kde.evennessmod_cdr1, new.df))
kde.evenness.1$n_trait <- seq(2, 9, by = 0.1)
kde.evenness.1$lwr <- kde.evenness.1$fit - kde.evenness.1$se.fit
kde.evenness.1$upr <- kde.evenness.1$fit + kde.evenness.1$se.fit
kde.dispersion.1 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdr1, new.df))
kde.dispersion.1$n_trait <- seq(2, 9, by = 0.1)
kde.dispersion.1$lwr <- kde.dispersion.1$fit - kde.dispersion.1$se.fit
kde.dispersion.1$upr <- kde.dispersion.1$fit + kde.dispersion.1$se.fit

fr.2 <- as.data.frame(predictSE.lme(fricmod_cdr2, new.df))
fr.2$n_trait <- seq(2, 9, by = 0.1)
fr.2$lwr <- fr.2$fit - fr.2$se.fit
fr.2$upr <- fr.2$fit + fr.2$se.fit
fe.2 <- as.data.frame(predictSE.lme(fevemod_cdr2, new.df))
fe.2$n_trait <- seq(2, 9, by = 0.1)
fe.2$lwr <- fe.2$fit - fe.2$se.fit
fe.2$upr <- fe.2$fit + fe.2$se.fit
fdis.2 <- as.data.frame(predictSE.lme(fdismod_cdr2, new.df))
fdis.2$n_trait <- seq(2, 9, by = 0.1)
fdis.2$lwr <- fdis.2$fit - fdis.2$se.fit
fdis.2$upr <- fdis.2$fit + fdis.2$se.fit
fdiv.2 <- as.data.frame(predictSE.lme(fdivmod_cdr2, new.df))
fdiv.2$n_trait <- seq(2, 9, by = 0.1)
fdiv.2$lwr <- fdiv.2$fit - fdiv.2$se.fit
fdiv.2$upr <- fdiv.2$fit + fdiv.2$se.fit
rq.2 <- as.data.frame(predictSE.lme(raoqmod_cdr2, new.df))
rq.2$n_trait <- seq(2, 9, by = 0.1)
rq.2$lwr <- rq.2$fit - rq.2$se.fit
rq.2$upr <- rq.2$fit + rq.2$se.fit
kde.alpha.2 <- as.data.frame(predictSE.lme(kde.alphamod_cdr2, new.df))
kde.alpha.2$n_trait <- seq(2, 9, by = 0.1)
kde.alpha.2$lwr <- kde.alpha.2$fit - kde.alpha.2$se.fit
kde.alpha.2$upr <- kde.alpha.2$fit + kde.alpha.2$se.fit
kde.evenness.2 <- as.data.frame(predictSE.lme(kde.evennessmod_cdr2, new.df))
kde.evenness.2$n_trait <- seq(2, 9, by = 0.1)
kde.evenness.2$lwr <- kde.evenness.2$fit - kde.evenness.2$se.fit
kde.evenness.2$upr <- kde.evenness.2$fit + kde.evenness.2$se.fit
kde.dispersion.2 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdr2, new.df))
kde.dispersion.2$n_trait <- seq(2, 9, by = 0.1)
kde.dispersion.2$lwr <- kde.dispersion.2$fit - kde.dispersion.2$se.fit
kde.dispersion.2$upr <- kde.dispersion.2$fit + kde.dispersion.2$se.fit

fr.3 <- as.data.frame(predictSE.lme(fricmod_cdr3, new.df))
fr.3$n_trait <- seq(2, 9, by = 0.1)
fr.3$lwr <- fr.3$fit - fr.3$se.fit
fr.3$upr <- fr.3$fit + fr.3$se.fit
fe.3 <- as.data.frame(predictSE.lme(fevemod_cdr3, new.df))
fe.3$n_trait <- seq(2, 9, by = 0.1)
fe.3$lwr <- fe.3$fit - fe.3$se.fit
fe.3$upr <- fe.3$fit + fe.3$se.fit
fdis.3 <- as.data.frame(predictSE.lme(fdismod_cdr3, new.df))
fdis.3$n_trait <- seq(2, 9, by = 0.1)
fdis.3$lwr <- fdis.3$fit - fdis.3$se.fit
fdis.3$upr <- fdis.3$fit + fdis.3$se.fit
fdiv.3 <- as.data.frame(predictSE.lme(fdivmod_cdr3, new.df))
fdiv.3$n_trait <- seq(2, 9, by = 0.1)
fdiv.3$lwr <- fdiv.3$fit - fdiv.3$se.fit
fdiv.3$upr <- fdiv.3$fit + fdiv.3$se.fit
rq.3 <- as.data.frame(predictSE.lme(raoqmod_cdr3, new.df))
rq.3$n_trait <- seq(2, 9, by = 0.1)
rq.3$lwr <- rq.3$fit - rq.3$se.fit
rq.3$upr <- rq.3$fit + rq.3$se.fit
kde.alpha.3 <- as.data.frame(predictSE.lme(kde.alphamod_cdr3, new.df))
kde.alpha.3$n_trait <- seq(2, 9, by = 0.1)
kde.alpha.3$lwr <- kde.alpha.3$fit - kde.alpha.3$se.fit
kde.alpha.3$upr <- kde.alpha.3$fit + kde.alpha.3$se.fit
kde.evenness.3 <- as.data.frame(predictSE.lme(kde.evennessmod_cdr3, new.df))
kde.evenness.3$n_trait <- seq(2, 9, by = 0.1)
kde.evenness.3$lwr <- kde.evenness.3$fit - kde.evenness.3$se.fit
kde.evenness.3$upr <- kde.evenness.3$fit + kde.evenness.3$se.fit
kde.dispersion.3 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdr3, new.df))
kde.dispersion.3$n_trait <- seq(2, 9, by = 0.1)
kde.dispersion.3$lwr <- kde.dispersion.3$fit - kde.dispersion.3$se.fit
kde.dispersion.3$upr <- kde.dispersion.3$fit + kde.dispersion.3$se.fit

fr.4 <- as.data.frame(predictSE.lme(fricmod_cdr4, new.df))
fr.4$n_trait <- seq(2, 9, by = 0.1)
fr.4$lwr <- fr.4$fit - fr.4$se.fit
fr.4$upr <- fr.4$fit + fr.4$se.fit
fe.4 <- as.data.frame(predictSE.lme(fevemod_cdr4, new.df))
fe.4$n_trait <- seq(2, 9, by = 0.1)
fe.4$lwr <- fe.4$fit - fe.4$se.fit
fe.4$upr <- fe.4$fit + fe.4$se.fit
fdis.4 <- as.data.frame(predictSE.lme(fdismod_cdr4, new.df))
fdis.4$n_trait <- seq(2, 9, by = 0.1)
fdis.4$lwr <- fdis.4$fit - fdis.4$se.fit
fdis.4$upr <- fdis.4$fit + fdis.4$se.fit
fdiv.4 <- as.data.frame(predictSE.lme(fdivmod_cdr4, new.df))
fdiv.4$n_trait <- seq(2, 9, by = 0.1)
fdiv.4$lwr <- fdiv.4$fit - fdiv.4$se.fit
fdiv.4$upr <- fdiv.4$fit + fdiv.4$se.fit
rq.4 <- as.data.frame(predictSE.lme(raoqmod_cdr4, new.df))
rq.4$n_trait <- seq(2, 9, by = 0.1)
rq.4$lwr <- rq.4$fit - rq.4$se.fit
rq.4$upr <- rq.4$fit + rq.4$se.fit
kde.alpha.4 <- as.data.frame(predictSE.lme(kde.alphamod_cdr4, new.df))
kde.alpha.4$n_trait <- seq(2, 9, by = 0.1)
kde.alpha.4$lwr <- kde.alpha.4$fit - kde.alpha.4$se.fit
kde.alpha.4$upr <- kde.alpha.4$fit + kde.alpha.4$se.fit
kde.evenness.4 <- as.data.frame(predictSE.lme(kde.evennessmod_cdr4, new.df))
kde.evenness.4$n_trait <- seq(2, 9, by = 0.1)
kde.evenness.4$lwr <- kde.evenness.4$fit - kde.evenness.4$se.fit
kde.evenness.4$upr <- kde.evenness.4$fit + kde.evenness.4$se.fit
kde.dispersion.4 <- as.data.frame(predictSE.lme(kde.dispersionmod_cdr4, new.df))
kde.dispersion.4$n_trait <- seq(2, 9, by = 0.1)
kde.dispersion.4$lwr <- kde.dispersion.4$fit - kde.dispersion.4$se.fit
kde.dispersion.4$upr <- kde.dispersion.4$fit + kde.dispersion.4$se.fit

#### SEV ####
sev.black$community <- "SEV1"
sev.blue$community <- "SEV2"
sev.full <- rbind(sev.blue[,-10], sev.black[,-10])

sev.full <- sev.full %>% group_by(community, n_trait) %>%
  summarize (FRic = mean(FRic, na.rm = TRUE), FEve = mean(FEve, na.rm = TRUE), 
             FDis = mean(FDis, na.rm = TRUE), FDiv = mean(FDiv, na.rm = TRUE), 
             RaoQ = mean(RaoQ, na.rm = TRUE), kde.alpha = mean(kde.alpha), 
             kde.evenness = mean(kde.evenness), 
             kde.dispersion = mean(kde.dispersion))

new.df1 <- data.frame(n_trait = seq(2, 10, by = 0.1))

fr.blue <- as.data.frame(predictSE.lme(fricmod_blue, new.df1))
fr.blue$n_trait <- seq(2, 10, by = 0.1)
fr.blue$lwr <- fr.blue$fit - fr.blue$se.fit
fr.blue$upr <- fr.blue$fit + fr.blue$se.fit
fe.blue <- as.data.frame(predictSE.lme(fevemod_blue, new.df1))
fe.blue$n_trait <- seq(2, 10, by = 0.1)
fe.blue$lwr <- fe.blue$fit - fe.blue$se.fit
fe.blue$upr <- fe.blue$fit + fe.blue$se.fit
fdis.blue <- as.data.frame(predictSE.lme(fdismod_blue, new.df1))
fdis.blue$n_trait <- seq(2, 10, by = 0.1)
fdis.blue$lwr <- fdis.blue$fit - fdis.blue$se.fit
fdis.blue$upr <- fdis.blue$fit + fdis.blue$se.fit
fdiv.blue <- as.data.frame(predictSE.lme(fdivmod_blue, new.df1))
fdiv.blue$n_trait <- seq(2, 10, by = 0.1)
fdiv.blue$lwr <- fdiv.blue$fit - fdiv.blue$se.fit
fdiv.blue$upr <- fdiv.blue$fit + fdiv.blue$se.fit
rq.blue <- as.data.frame(predictSE.lme(raoqmod_blue, new.df1))
rq.blue$n_trait <- seq(2, 10, by = 0.1)
rq.blue$lwr <- rq.blue$fit - rq.blue$se.fit
rq.blue$upr <- rq.blue$fit + rq.blue$se.fit
kde.alpha.blue <- as.data.frame(predictSE.lme(kde.alphamod_blue, new.df1))
kde.alpha.blue$n_trait <- seq(2, 10, by = 0.1)
kde.alpha.blue$lwr <- kde.alpha.blue$fit - kde.alpha.blue$se.fit
kde.alpha.blue$upr <- kde.alpha.blue$fit + kde.alpha.blue$se.fit
kde.evenness.blue <- as.data.frame(predictSE.lme(kde.evennessmod_blue, new.df1))
kde.evenness.blue$n_trait <- seq(2, 10, by = 0.1)
kde.evenness.blue$lwr <- kde.evenness.blue$fit - kde.evenness.blue$se.fit
kde.evenness.blue$upr <- kde.evenness.blue$fit + kde.evenness.blue$se.fit
kde.dispersion.blue <- as.data.frame(predictSE.lme(kde.dispersionmod_blue, new.df1))
kde.dispersion.blue$n_trait <- seq(2, 10, by = 0.1)
kde.dispersion.blue$lwr <- kde.dispersion.blue$fit - kde.dispersion.blue$se.fit
kde.dispersion.blue$upr <- kde.dispersion.blue$fit + kde.dispersion.blue$se.fit

fr.black <- as.data.frame(predictSE.lme(fricmod_black, new.df1))
fr.black$n_trait <- seq(2, 10, by = 0.1)
fr.black$lwr <- fr.black$fit - fr.black$se.fit
fr.black$upr <- fr.black$fit + fr.black$se.fit
fe.black <- as.data.frame(predictSE.lme(fevemod_black, new.df1))
fe.black$n_trait <- seq(2, 10, by = 0.1)
fe.black$lwr <- fe.black$fit - fe.black$se.fit
fe.black$upr <- fe.black$fit + fe.black$se.fit
fdis.black <- as.data.frame(predictSE.lme(fdismod_black, new.df1))
fdis.black$n_trait <- seq(2, 10, by = 0.1)
fdis.black$lwr <- fdis.black$fit - fdis.black$se.fit
fdis.black$upr <- fdis.black$fit + fdis.black$se.fit
fdiv.black <- as.data.frame(predictSE.lme(fdivmod_black, new.df1))
fdiv.black$n_trait <- seq(2, 10, by = 0.1)
fdiv.black$lwr <- fdiv.black$fit - fdiv.black$se.fit
fdiv.black$upr <- fdiv.black$fit + fdiv.black$se.fit
rq.black <- as.data.frame(predictSE.lme(raoqmod_black, new.df1))
rq.black$n_trait <- seq(2, 10, by = 0.1)
rq.black$lwr <- rq.black$fit - rq.black$se.fit
rq.black$upr <- rq.black$fit + rq.black$se.fit
kde.alpha.black <- as.data.frame(predictSE.lme(kde.alphamod_black, new.df1))
kde.alpha.black$n_trait <- seq(2, 10, by = 0.1)
kde.alpha.black$lwr <- kde.alpha.black$fit - kde.alpha.black$se.fit
kde.alpha.black$upr <- kde.alpha.black$fit + kde.alpha.black$se.fit
kde.evenness.black <- as.data.frame(predictSE.lme(kde.evennessmod_black, new.df1))
kde.evenness.black$n_trait <- seq(2, 10, by = 0.1)
kde.evenness.black$lwr <- kde.evenness.black$fit - kde.evenness.black$se.fit
kde.evenness.black$upr <- kde.evenness.black$fit + kde.evenness.black$se.fit
kde.dispersion.black <- as.data.frame(predictSE.lme(kde.dispersionmod_black, new.df1))
kde.dispersion.black$n_trait <- seq(2, 10, by = 0.1)
kde.dispersion.black$lwr <- kde.dispersion.black$fit - kde.dispersion.black$se.fit
kde.dispersion.black$upr <- kde.dispersion.black$fit + kde.dispersion.black$se.fit


Frich_g <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.1, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.2, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.3, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.4, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.blue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.black, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = fr.1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FRic, color = community), data = cdr.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = fr.blue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FRic, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#731279", "#00B7FF", "#075A13","#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FRich") +
  theme_pubr()

FEve_g <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.1, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.2, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.3, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.4, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.blue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.black, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = fe.1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FEve, color = community), data = cdr.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = fe.blue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FEve, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#731279", "#00B7FF", "#075A13","#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FEve") +
  theme_pubr()

FDis_g <-ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.1, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.2, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.3, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.4, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.blue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.black, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = fdis.1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FDis, color = community), data = cdr.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.blue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FDis, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#731279", "#00B7FF", "#075A13","#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FDis") +
  theme_pubr()

FDiv_g <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.1, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.2, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.3, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.4, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.blue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.black, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = fdiv.1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FDiv, color = community), data = cdr.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.blue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FDiv, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#731279", "#00B7FF", "#075A13","#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FDiv") +
  theme_pubr()

RaoQ_g <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.1, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.2, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.3, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.4, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.blue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.black, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = rq.1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = RaoQ, color = community), data = cdr.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = rq.blue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = RaoQ, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#731279", "#00B7FF", "#075A13","#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "Rao Q") +
  theme_pubr()

kderichness_g <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.1, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.2, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.3, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.4, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.blue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.black, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.alpha, color = community), data = cdr.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.blue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.alpha, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#731279", "#00B7FF", "#075A13","#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Richness") +
  theme_pubr()

kdeevenness_g <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.1, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.2, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.3, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.4, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.blue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.black, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.evenness, color = community), data = cdr.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.blue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.evenness, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#731279", "#00B7FF", "#075A13","#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Evenness") +
  theme_pubr()

kdedispersion_g <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.1, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.2, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.3, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.4, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.blue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.black, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.dispersion, color = community), data = cdr.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.blue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.dispersion, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#731279", "#00B7FF", "#075A13","#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Dispersion") +
  theme_pubr()



####### Predicted values #######
######### ECULIDEAN  ###########
# create combined data set of raw data
# create combined data set of raw data
cdre.1$community <- "CDR1"
cdre.2$community <- "CDR2"
cdre.3$community <- "CDR3"
cdre.4$community <- "CDR4"
cdre.full <- rbind(cdre.1, cdre.2, cdre.3, cdre.4)

cdre.full <- cdre.full %>% group_by(community, n_trait) %>%
  summarize (FRic = mean(FRic), FEve = mean(FEve), FDis = mean(FDis),
             FDiv = mean(FDiv), RaoQ = mean(RaoQ), kde.alpha = mean(kde.alpha), kde.evenness = mean(kde.evenness), kde.dispersion = mean(kde.dispersion))

new.df2 <- data.frame(n_trait = seq(2, 9, by = .1))


############## cdr  ########################

fr.1e <- as.data.frame(predictSE.lme(fricmod_cdre1, new.df2))
fr.1e$n_trait <- seq(2, 9, by = .1)
fr.1e$lwr <- fr.1e$fit - fr.1e$se.fit
fr.1e$upr <- fr.1e$fit + fr.1e$se.fit
fe.1e <- as.data.frame(predictSE.lme(fevemod_cdre1, new.df2))
fe.1e$n_trait <- seq(2, 9, by = .1)
fe.1e$lwr <- fe.1e$fit - fe.1e$se.fit
fe.1e$upr <- fe.1e$fit + fe.1e$se.fit
fdis.1e <- as.data.frame(predictSE.lme(fdismod_cdre1, new.df2))
fdis.1e$n_trait <- seq(2, 9, by = .1)
fdis.1e$lwr <- fdis.1e$fit - fdis.1e$se.fit
fdis.1e$upr <- fdis.1e$fit + fdis.1e$se.fit
fdiv.1e <- as.data.frame(predictSE.lme(fdivmod_cdre1, new.df2))
fdiv.1e$n_trait <- seq(2, 9, by = .1)
fdiv.1e$lwr <- fdiv.1e$fit - fdiv.1e$se.fit
fdiv.1e$upr <- fdiv.1e$fit + fdiv.1e$se.fit
rq.1e <- as.data.frame(predictSE.lme(raoqmod_cdre1, new.df2))
rq.1e$n_trait <- seq(2, 9, by = .1)
rq.1e$lwr <- rq.1e$fit - rq.1e$se.fit
rq.1e$upr <- rq.1e$fit + rq.1e$se.fit
kde.alpha.1e <- as.data.frame(predictSE.lme(kde.alphamod_cdre1, new.df2))
kde.alpha.1e$n_trait <- seq(2, 9, by = .1)
kde.alpha.1e$lwr <- kde.alpha.1e$fit - kde.alpha.1e$se.fit
kde.alpha.1e$upr <- kde.alpha.1e$fit + kde.alpha.1e$se.fit
kde.evenness.1e <- as.data.frame(predictSE.lme(kde.evennessmod_cdre1, new.df2))
kde.evenness.1e$n_trait <- seq(2, 9, by = .1)
kde.evenness.1e$lwr <- kde.evenness.1$fit - kde.evenness.1$se.fit
kde.evenness.1e$upr <- kde.evenness.1$fit + kde.evenness.1$se.fit
kde.dispersion.1e <- as.data.frame(predictSE.lme(kde.dispersionmod_cdre1, new.df2))
kde.dispersion.1e$n_trait <- seq(2, 9, by = .1)
kde.dispersion.1e$lwr <- kde.dispersion.1e$fit - kde.dispersion.1e$se.fit
kde.dispersion.1e$upr <- kde.dispersion.1e$fit + kde.dispersion.1e$se.fit

fr.2e <- as.data.frame(predictSE.lme(fricmod_cdre2, new.df2))
fr.2e$n_trait <- seq(2, 9, by = .1)
fr.2e$lwr <- fr.2e$fit - fr.2e$se.fit
fr.2e$upr <- fr.2e$fit + fr.2e$se.fit
fe.2e <- as.data.frame(predictSE.lme(fevemod_cdre2, new.df2))
fe.2e$n_trait <- seq(2, 9, by = .1)
fe.2e$lwr <- fe.2e$fit - fe.2e$se.fit
fe.2e$upr <- fe.2e$fit + fe.2e$se.fit
fdis.2e <- as.data.frame(predictSE.lme(fdismod_cdre2, new.df2))
fdis.2e$n_trait <- seq(2, 9, by = .1)
fdis.2e$lwr <- fdis.2e$fit - fdis.2e$se.fit
fdis.2e$upr <- fdis.2e$fit + fdis.2e$se.fit
fdiv.2e <- as.data.frame(predictSE.lme(fdivmod_cdre2, new.df2))
fdiv.2e$n_trait <- seq(2, 9, by = 0.1)
fdiv.2e$lwr <- fdiv.2e$fit - fdiv.2e$se.fit
fdiv.2e$upr <- fdiv.2e$fit + fdiv.2e$se.fit
rq.2e <- as.data.frame(predictSE.lme(raoqmod_cdre2, new.df2))
rq.2e$n_trait <- seq(2, 9, by = 0.1)
rq.2e$lwr <- rq.2e$fit - rq.2e$se.fit
rq.2e$upr <- rq.2e$fit + rq.2e$se.fit
kde.alpha.2e <- as.data.frame(predictSE.lme(kde.alphamod_cdre2, new.df2))
kde.alpha.2e$n_trait <- seq(2, 9, by = 0.1)
kde.alpha.2e$lwr <- kde.alpha.2e$fit - kde.alpha.2e$se.fit
kde.alpha.2e$upr <- kde.alpha.2e$fit + kde.alpha.2e$se.fit
kde.evenness.2e <- as.data.frame(predictSE.lme(kde.evennessmod_cdre2, new.df2))
kde.evenness.2e$n_trait <- seq(2, 9, by = 0.1)
kde.evenness.2e$lwr <- kde.evenness.2e$fit - kde.evenness.2e$se.fit
kde.evenness.2e$upr <- kde.evenness.2e$fit + kde.evenness.2e$se.fit
kde.dispersion.2e <- as.data.frame(predictSE.lme(kde.dispersionmod_cdre2, new.df2))
kde.dispersion.2e$n_trait <- seq(2, 9, by = 0.1)
kde.dispersion.2e$lwr <- kde.dispersion.2e$fit - kde.dispersion.2e$se.fit
kde.dispersion.2e$upr <- kde.dispersion.2e$fit + kde.dispersion.2e$se.fit

fr.3e <- as.data.frame(predictSE.lme(fricmod_cdre3, new.df2))
fr.3e$n_trait <- seq(2, 9, by = 0.1)
fr.3e$lwr <- fr.3e$fit - fr.3e$se.fit
fr.3e$upr <- fr.3e$fit + fr.3e$se.fit
fe.3e <- as.data.frame(predictSE.lme(fevemod_cdre3, new.df2))
fe.3e$n_trait <- seq(2, 9, by = 0.1)
fe.3e$lwr <- fe.3e$fit - fe.3e$se.fit
fe.3e$upr <- fe.3e$fit + fe.3e$se.fit
fdis.3e <- as.data.frame(predictSE.lme(fdismod_cdre3, new.df2))
fdis.3e$n_trait <- seq(2, 9, by = 0.1)
fdis.3e$lwr <- fdis.3e$fit - fdis.3e$se.fit
fdis.3e$upr <- fdis.3e$fit + fdis.3e$se.fit
fdiv.3e <- as.data.frame(predictSE.lme(fdivmod_cdre3, new.df2))
fdiv.3e$n_trait <- seq(2, 9, by = 0.1)
fdiv.3e$lwr <- fdiv.3e$fit - fdiv.3e$se.fit
fdiv.3e$upr <- fdiv.3e$fit + fdiv.3e$se.fit
rq.3e <- as.data.frame(predictSE.lme(raoqmod_cdre3, new.df2))
rq.3e$n_trait <- seq(2, 9, by = 0.1)
rq.3e$lwr <- rq.3e$fit - rq.3e$se.fit
rq.3e$upr <- rq.3e$fit + rq.3e$se.fit
kde.alpha.3e <- as.data.frame(predictSE.lme(kde.alphamod_cdre3, new.df2))
kde.alpha.3e$n_trait <- seq(2, 9, by = 0.1)
kde.alpha.3e$lwr <- kde.alpha.3e$fit - kde.alpha.3e$se.fit
kde.alpha.3e$upr <- kde.alpha.3e$fit + kde.alpha.3e$se.fit
kde.evenness.3e <- as.data.frame(predictSE.lme(kde.evennessmod_cdre3, new.df2))
kde.evenness.3e$n_trait <- seq(2, 9, by = 0.1)
kde.evenness.3e$lwr <- kde.evenness.3e$fit - kde.evenness.3e$se.fit
kde.evenness.3e$upr <- kde.evenness.3e$fit + kde.evenness.3e$se.fit
kde.dispersion.3e <- as.data.frame(predictSE.lme(kde.dispersionmod_cdre3, new.df2))
kde.dispersion.3e$n_trait <- seq(2, 9, by = 0.1)
kde.dispersion.3e$lwr <- kde.dispersion.3e$fit - kde.dispersion.3e$se.fit
kde.dispersion.3e$upr <- kde.dispersion.3e$fit + kde.dispersion.3e$se.fit

fr.4e <- as.data.frame(predictSE.lme(fricmod_cdre4, new.df2))
fr.4e$n_trait <- seq(2, 9, by = 0.1)
fr.4e$lwr <- fr.4e$fit - fr.4e$se.fit
fr.4e$upr <- fr.4e$fit + fr.4e$se.fit
fe.4e <- as.data.frame(predictSE.lme(fevemod_cdre4, new.df2))
fe.4e$n_trait <- seq(2, 9, by = 0.1)
fe.4e$lwr <- fe.4e$fit - fe.4e$se.fit
fe.4e$upr <- fe.4e$fit + fe.4e$se.fit
fdis.4e <- as.data.frame(predictSE.lme(fdismod_cdre4, new.df2))
fdis.4e$n_trait <- seq(2, 9, by = 0.1)
fdis.4e$lwr <- fdis.4e$fit - fdis.4e$se.fit
fdis.4e$upr <- fdis.4$fit + fdis.4$se.fit
fdiv.4e <- as.data.frame(predictSE.lme(fdivmod_cdre4, new.df2))
fdiv.4e$n_trait <- seq(2, 9, by = 0.1)
fdiv.4e$lwr <- fdiv.4e$fit - fdiv.4e$se.fit
fdiv.4e$upr <- fdiv.4e$fit + fdiv.4e$se.fit
rq.4e <- as.data.frame(predictSE.lme(raoqmod_cdre4, new.df2))
rq.4e$n_trait <- seq(2, 9, by = 0.1)
rq.4e$lwr <- rq.4e$fit - rq.4e$se.fit
rq.4e$upr <- rq.4e$fit + rq.4e$se.fit
kde.alpha.4e <- as.data.frame(predictSE.lme(kde.alphamod_cdre4, new.df2))
kde.alpha.4e$n_trait <- seq(2, 9, by = 0.1)
kde.alpha.4e$lwr <- kde.alpha.4e$fit - kde.alpha.4e$se.fit
kde.alpha.4e$upr <- kde.alpha.4e$fit + kde.alpha.4e$se.fit
kde.evenness.4e <- as.data.frame(predictSE.lme(kde.evennessmod_cdre4, new.df2))
kde.evenness.4e$n_trait <- seq(2, 9, by = 0.1)
kde.evenness.4e$lwr <- kde.evenness.4e$fit - kde.evenness.4e$se.fit
kde.evenness.4e$upr <- kde.evenness.4e$fit + kde.evenness.4e$se.fit
kde.dispersion.4e <- as.data.frame(predictSE.lme(kde.dispersionmod_cdre4, new.df2))
kde.dispersion.4e$n_trait <- seq(2, 9, by = 0.1)
kde.dispersion.4e$lwr <- kde.dispersion.4e$fit - kde.dispersion.4e$se.fit
kde.dispersion.4e$upr <- kde.dispersion.4e$fit + kde.dispersion.4e$se.fit


#### SEV #########
# create combined data set of raw data
sev.eblack$community <- "SEV1"
sev.eblue$community <- "SEV2"
sev.efull <- rbind(sev.eblue[,-10], sev.eblack[,-10])

sev.efull <- sev.efull %>% group_by(community, n_trait) %>%
  summarize (FRic = mean(FRic, na.rm = TRUE), FEve = mean(FEve, na.rm = TRUE), 
             FDis = mean(FDis, na.rm = TRUE), FDiv = mean(FDiv, na.rm = TRUE), 
             RaoQ = mean(RaoQ, na.rm = TRUE), kde.alpha = mean(kde.alpha), kde.evenness = mean(kde.evenness), kde.dispersion = mean(kde.dispersion))

new.df3 <- data.frame(n_trait = seq(2, 9, by = 0.1))

fr.eblue <- as.data.frame(predictSE.lme(fricmod_eblue, new.df3))
fr.eblue$n_trait <- seq(2, 9, by = 0.1)
fr.eblue$lwr <- fr.eblue$fit - fr.eblue$se.fit
fr.eblue$upr <- fr.eblue$fit + fr.eblue$se.fit
fe.eblue <- as.data.frame(predictSE.lme(fevemod_eblue, new.df3))
fe.eblue$n_trait <- seq(2, 9, by = 0.1)
fe.eblue$lwr <- fe.eblue$fit - fe.eblue$se.fit
fe.eblue$upr <- fe.eblue$fit + fe.eblue$se.fit
fdis.eblue <- as.data.frame(predictSE.lme(fdismod_eblue, new.df3))
fdis.eblue$n_trait <- seq(2, 9, by = 0.1)
fdis.eblue$lwr <- fdis.eblue$fit - fdis.eblue$se.fit
fdis.eblue$upr <- fdis.eblue$fit + fdis.eblue$se.fit
fdiv.eblue <- as.data.frame(predictSE.lme(fdivmod_eblue, new.df3))
fdiv.eblue$n_trait <- seq(2, 9, by = 0.1)
fdiv.eblue$lwr <- fdiv.eblue$fit - fdiv.eblue$se.fit
fdiv.eblue$upr <- fdiv.eblue$fit + fdiv.eblue$se.fit
rq.eblue <- as.data.frame(predictSE.lme(raoqmod_eblue, new.df3))
rq.eblue$n_trait <- seq(2, 9, by = 0.1)
rq.eblue$lwr <- rq.eblue$fit - rq.eblue$se.fit
rq.eblue$upr <- rq.eblue$fit + rq.eblue$se.fit
kde.alpha.eblue <- as.data.frame(predictSE.lme(kde.alphamod_eblue, new.df3))
kde.alpha.eblue$n_trait <- seq(2, 9, by = 0.1)
kde.alpha.eblue$lwr <- kde.alpha.eblue$fit - kde.alpha.eblue$se.fit
kde.alpha.eblue$upr <- kde.alpha.eblue$fit + kde.alpha.eblue$se.fit
kde.evenness.eblue <- as.data.frame(predictSE.lme(kde.evennessmod_eblue, new.df3))
kde.evenness.eblue$n_trait <- seq(2, 9, by = 0.1)
kde.evenness.eblue$lwr <- kde.evenness.eblue$fit - kde.evenness.eblue$se.fit
kde.evenness.eblue$upr <- kde.evenness.eblue$fit + kde.evenness.eblue$se.fit
kde.dispersion.eblue <- as.data.frame(predictSE.lme(kde.dispersionmod_eblue, new.df3))
kde.dispersion.eblue$n_trait <- seq(2, 9, by = 0.1)
kde.dispersion.eblue$lwr <- kde.dispersion.eblue$fit - kde.dispersion.eblue$se.fit
kde.dispersion.eblue$upr <- kde.dispersion.eblue$fit + kde.dispersion.eblue$se.fit

fr.eblack <- as.data.frame(predictSE.lme(fricmod_eblack, new.df3))
fr.eblack$n_trait <- seq(2, 9, by = 0.1)
fr.eblack$lwr <- fr.eblack$fit - fr.eblack$se.fit
fr.eblack$upr <- fr.eblack$fit + fr.eblack$se.fit
fe.eblack <- as.data.frame(predictSE.lme(fevemod_eblack, new.df3))
fe.eblack$n_trait <- seq(2, 9, by = 0.1)
fe.eblack$lwr <- fe.eblack$fit - fe.eblack$se.fit
fe.eblack$upr <- fe.eblack$fit + fe.eblack$se.fit
fdis.eblack <- as.data.frame(predictSE.lme(fdismod_eblack, new.df3))
fdis.eblack$n_trait <- seq(2, 9, by = 0.1)
fdis.eblack$lwr <- fdis.eblack$fit - fdis.eblack$se.fit
fdis.eblack$upr <- fdis.eblack$fit + fdis.eblack$se.fit
fdiv.eblack <- as.data.frame(predictSE.lme(fdivmod_eblack, new.df3))
fdiv.eblack$n_trait <- seq(2, 9, by = 0.1)
fdiv.eblack$lwr <- fdiv.eblack$fit - fdiv.eblack$se.fit
fdiv.eblack$upr <- fdiv.eblack$fit + fdiv.eblack$se.fit
rq.eblack <- as.data.frame(predictSE.lme(raoqmod_eblack, new.df3))
rq.eblack$n_trait <- seq(2, 9, by = 0.1)
rq.eblack$lwr <- rq.eblack$fit - rq.eblack$se.fit
rq.eblack$upr <- rq.eblack$fit + rq.eblack$se.fit
kde.alpha.eblack <- as.data.frame(predictSE.lme(kde.alphamod_eblack, new.df3))
kde.alpha.eblack$n_trait <- seq(2, 9, by = 0.1)
kde.alpha.eblack$lwr <- kde.alpha.eblack$fit - kde.alpha.eblack$se.fit
kde.alpha.eblack$upr <- kde.alpha.eblack$fit + kde.alpha.eblack$se.fit
kde.evenness.eblack <- as.data.frame(predictSE.lme(kde.evennessmod_eblack, new.df3))
kde.evenness.eblack$n_trait <- seq(2, 9, by = 0.1)
kde.evenness.eblack$lwr <- kde.evenness.eblack$fit - kde.evenness.eblack$se.fit
kde.evenness.eblack$upr <- kde.evenness.eblack$fit + kde.evenness.eblack$se.fit
kde.dispersion.eblack <- as.data.frame(predictSE.lme(kde.dispersionmod_eblack, new.df3))
kde.dispersion.eblack$n_trait <- seq(2, 9, by = 0.1)
kde.dispersion.eblack$lwr <- kde.dispersion.eblack$fit - kde.dispersion.eblack$se.fit
kde.dispersion.eblack$upr <- kde.dispersion.eblack$fit + kde.dispersion.eblack$se.fit


Frich_e <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.1e, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.2e, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.3e, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.4e, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.eblue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.eblack, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = fr.1e, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.2e, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.3e, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.4e, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FRic, color = community), data = cdre.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = fr.eblue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FRic, color = community), data = sev.efull, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#731279", "#00B7FF", "#075A13","#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FRich") +
  theme_pubr()

FEve_e <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.1e, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.2e, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.3e, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.4e, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.eblue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.eblack, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = fe.1e, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.2e, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.3e, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.4e, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FEve, color = community), data = cdre.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = fe.eblue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FEve, color = community), data = sev.efull, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#731279", "#00B7FF", "#075A13","#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FEve") +
  theme_pubr()

FDis_e <-ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.1e, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.2e, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.3e, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.4e, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.eblue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.eblack, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = fdis.1e, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.2e, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.3e, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.4e, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FDis, color = community), data = cdre.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.eblue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FDis, color = community), data = sev.efull, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#731279", "#00B7FF", "#075A13","#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FDis") +
  theme_pubr()

FDiv_e <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.1e, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.2e, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.3e, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.4e, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.eblue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.eblack, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = fdiv.1e, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.2e, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.3e, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.4e, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FDiv, color = community), data = cdre.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.eblue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FDiv, color = community), data = sev.efull, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#731279", "#00B7FF", "#075A13","#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "FDiv") +
  theme_pubr()

RaoQ_e <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.1e, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.2e, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.3e, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.4e, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.eblue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.eblack, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = rq.1e, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.2e, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.3e, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.4e, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = RaoQ, color = community), data = cdre.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = rq.eblue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = RaoQ, color = community), data = sev.efull, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#731279", "#00B7FF", "#075A13","#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "Rao Q") +
  theme_pubr()

kderichness_e <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.1e, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.2e, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.3e, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.4e, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.eblue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.eblack, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.1e, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.2e, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.3e, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.4e, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.alpha, color = community), data = cdre.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.eblue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.alpha, color = community), data = sev.efull, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#731279", "#00B7FF", "#075A13","#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Richness") +
  theme_pubr()

kdeevenness_e <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.1e, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.2e, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.3e, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.4e, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.eblue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.eblack, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.1e, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.2e, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.3e, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.4e, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.evenness, color = community), data = cdre.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.eblue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.evenness, color = community), data = sev.efull, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#731279", "#00B7FF", "#075A13","#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Evenness") +
  theme_pubr()

kdedispersion_e <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.1e, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.2e, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.3e, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.4e, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.eblue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.eblack, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.1e, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.2e, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.3e, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.4e, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.dispersion, color = community), data = cdre.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.eblue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.dispersion, color = community), data = sev.efull, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#731279", "#00B7FF", "#075A13","#000000", "navyblue"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Dispersion") +
  theme_pubr()


png(here('Figures/n_trait_full.png'), height = 9, width = 12, units = 'in', res = 150)
ggarrange(plotlist = list(Frich_g, Frich_e, kderichness_g, kderichness_e, FEve_g, FEve_e, 
                          kdeevenness_g, kdeevenness_e, FDis_g, FDis_e, FDiv_g, FDiv_e, 
                          RaoQ_g, RaoQ_e, kdedispersion_g, kdedispersion_e), ncol = 4, nrow = 4, 
          labels = c("A", "B", 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P'), common.legend = TRUE)
dev.off()
