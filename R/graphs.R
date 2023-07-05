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


############# PREDICTED VALUES & AVERAGES ##################
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
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.1, alpha = 0.2, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.2, alpha = 0.2, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.3, alpha = 0.2, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.4, alpha = 0.2, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.blue, alpha = 0.2, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.black, alpha = 0.2, fill = "black") +
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

png(here("Figures/cdr_Ntraits.png"))
ggarrange(plotlist = list(A, B, C, D, E, F, G, H), common.legend = TRUE,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
dev.off()

###SEV####


# create combined data set of raw data


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

png(here("Figures/Sev_Ntraits.png"))
ggarrange(plotlist = list(A, B, C, D, E, F, G, H), common.legend = TRUE, 
          labels = c("AA", "BB", "CC", "DD", "EE", "FF", "GG", "HH"))
dev.off()

png(here("Figures/n_trait.png"), height = 8, width = 13, units = 'in', res = 300)
ggarrange(plotlist = list(A, AA, B, BB, C, CC, D, DD, E, EE, F, FF, G, GG, H, HH), 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", 
                     "O", "P"), ncol = 4, nrow = 4, common.legend = TRUE, align = "hv")
dev.off()





