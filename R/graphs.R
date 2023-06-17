######## GRAPHS ##########

# libraries
library(ggplot2)
library(here)
library(nlme)

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







