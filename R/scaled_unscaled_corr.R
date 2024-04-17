#Correlation between metrics with scaled and unscaled 

#load libraries
library(nlme)
library(here)
library(AICcmodavg)
library(dplyr)

#read in scaled
cdr.1.sc <- read.csv(here("data/Cleaned/cdr1_sc.csv"))%>%
  subset(SR != "NA")
cdr.2.sc <- read.csv(here("data/Cleaned/cdr2_sc.csv"))%>%
  subset(SR != "NA")
cdr.3.sc <- read.csv(here("data/Cleaned/cdr3_sc.csv"))%>%
  subset(SR != "NA")
cdr.4.sc <- read.csv(here("data/Cleaned/cdr4_sc.csv"))%>%
  subset(SR != "NA")
sev.blue.sc <- read.csv(here("data/Cleaned/sevblue_sc.csv"))
#sev.blue <- sev.blue[-which(sev.blue$SR == 1),]
sev.black.sc <- read.csv(here("data/Cleaned/sevblack_sc.csv"))
colnames(sev.black.sc) <- paste(colnames(sev.black.sc),"sc",sep="_") 
#sev.black <- sev.black[-which(sev.black$SR ==1),]


#read in unscaled data
cdr.1 <- read.csv(here("data/Cleaned/cdr1.csv"))%>%
  subset(SR != "NA")
cdr.2 <- read.csv(here("data/Cleaned/cdr2.csv"))%>%
  subset(SR != "NA")
cdr.3 <- read.csv(here("data/Cleaned/cdr3.csv"))%>%
  subset(SR != "NA")
cdr.4 <- read.csv(here("data/Cleaned/cdr4.csv"))%>%
  subset(SR != "NA")
sev.blue <- read.csv(here("data/Cleaned/sevblue.csv"))
#sev.blue <- sev.blue[-which(sev.blue$SR == 1),]
sev.black <- read.csv(here("data/Cleaned/sevblack.csv"))
#sev.black <- sev.black[-which(sev.black$SR ==1),]


##merge the two
cdr1.comb <- left_join(cdr.1, cdr.1.sc, by = c("Plot", "Ring", "n_trait", "traits"))
cdr2.comb <- left_join(cdr.2, cdr.2.sc, by = c("Plot", "Ring", "n_trait", "traits"))
cdr3.comb <- left_join(cdr.3, cdr.3.sc, by = c("Plot", "Ring", "n_trait", "traits"))
cdr4.comb <- left_join(cdr.4, cdr.4.sc, by = c("Plot", "Ring", "n_trait", "traits"))
sev1.comb <- left_join(sev.blue, sev.blue.sc, by = c("sev.blueplots", "n_trait", "traits"))
sev2.comb <- cbind(sev.black, sev.black.sc)


#correlation between the two by community
#CDR1
cdr1.FRic <- lm(FRic.y~FRic.x, data = cdr1.comb)
summary(cdr1.FRic)
cor(cdr1.comb$FRic.y, cdr1.comb$FRic.x)

cdr1.kde.alpha <- lm(kde.alpha.y~kde.alpha.x, data = cdr1.comb)
summary(cdr1.kde.alpha)

cdr1.FEve <- lm(FEve.y~FEve.x, data = cdr1.comb)
summary(cdr1.FEve)

cdr1.kde.eve <- lm(kde.evenness.y~kde.evenness.x, data = cdr1.comb)
summary(cdr1.kde.eve)

cdr1.FDis <- lm(FDis.y~FDis.x, data = cdr1.comb)
summary(cdr1.FDis)

cdr1.kde.dis <- lm(kde.dispersion.y~kde.dispersion.x, data = cdr1.comb)
summary(cdr1.kde.dis)

cdr1.FDiv <- lm(FDiv.y~FDiv.x, data = cdr1.comb)
summary(cdr1.FDiv)

cdr1.rao <- lm(RaoQ.y~RaoQ.x, data = cdr1.comb)
summary(cdr1.rao)

#CDR2

cdr2.FRic <- lm(FRic.y~FRic.x, data = cdr2.comb)
summary(cdr2.FRic)

cdr2.kde.alpha <- lm(kde.alpha.y~kde.alpha.x, data = cdr2.comb)
summary(cdr2.kde.alpha)

cdr2.FEve <- lm(FEve.y~FEve.x, data = cdr2.comb)
summary(cdr2.FEve)

cdr2.kde.eve <- lm(kde.evenness.y~kde.evenness.x, data = cdr2.comb)
summary(cdr2.kde.eve)

cdr2.FDis <- lm(FDis.y~FDis.x, data = cdr2.comb)
summary(cdr2.FDis)

cdr2.kde.dis <- lm(kde.dispersion.y~kde.dispersion.x, data = cdr2.comb)
summary(cdr2.kde.dis)

cdr2.FDiv <- lm(FDiv.y~FDiv.x, data = cdr2.comb)
summary(cdr2.FDiv)

cdr2.rao <- lm(RaoQ.y~RaoQ.x, data = cdr2.comb)
summary(cdr2.rao)

#CDR3
cdr3.FRic <- lm(FRic.y~FRic.x, data = cdr3.comb)
summary(cdr3.FRic)

cdr3.kde.alpha <- lm(kde.alpha.y~kde.alpha.x, data = cdr3.comb)
summary(cdr3.kde.alpha)

cdr3.FEve <- lm(FEve.y~FEve.x, data = cdr3.comb)
summary(cdr3.FEve)

cdr3.kde.eve <- lm(kde.evenness.y~kde.evenness.x, data = cdr3.comb)
summary(cdr3.kde.eve)

cdr3.FDis <- lm(FDis.y~FDis.x, data = cdr3.comb)
summary(cdr3.FDis)

cdr3.kde.dis <- lm(kde.dispersion.y~kde.dispersion.x, data = cdr3.comb)
summary(cdr3.kde.dis)

cdr3.FDiv <- lm(FDiv.y~FDiv.x, data = cdr3.comb)
summary(cdr3.FDiv)

cdr3.rao <- lm(RaoQ.y~RaoQ.x, data = cdr3.comb)
summary(cdr3.rao)

#CDR4
cdr4.FRic <- lm(FRic.y~FRic.x, data = cdr4.comb)
summary(cdr4.FRic)

cdr4.kde.alpha <- lm(kde.alpha.y~kde.alpha.x, data = cdr4.comb)
summary(cdr4.kde.alpha)

cdr4.FEve <- lm(FEve.y~FEve.x, data = cdr4.comb)
summary(cdr4.FEve)

cdr4.kde.eve <- lm(kde.evenness.y~kde.evenness.x, data = cdr4.comb)
summary(cdr4.kde.eve)

cdr4.FDis <- lm(FDis.y~FDis.x, data = cdr4.comb)
summary(cdr4.FDis)

cdr4.kde.dis <- lm(kde.dispersion.y~kde.dispersion.x, data = cdr4.comb)
summary(cdr4.kde.dis)

cdr4.FDiv <- lm(FDiv.y~FDiv.x, data = cdr4.comb)
summary(cdr4.FDiv)

cdr4.rao <- lm(RaoQ.y~RaoQ.x, data = cdr4.comb)
summary(cdr4.rao)


#SEV1
sev1.FRic <- lm(FRic.y~FRic.x, data = sev1.comb)
summary(sev1.FRic)

sev1.kde.alpha <- lm(kde.alpha.y~kde.alpha.x, data = sev1.comb)
summary(sev1.kde.alpha)

sev1.FEve <- lm(FEve.y~FEve.x, data = sev1.comb)
summary(sev1.FEve)

sev1.kde.eve <- lm(kde.evenness.y~kde.evenness.x, data = sev1.comb)
summary(sev1.kde.eve)

sev1.FDis <- lm(FDis.y~FDis.x, data = sev1.comb)
summary(sev1.FDis)

sev1.kde.dis <- lm(kde.dispersion.y~kde.dispersion.x, data = sev1.comb)
summary(sev1.kde.dis)

sev1.FDiv <- lm(FDiv.y~FDiv.x, data = sev1.comb)
summary(sev1.FDiv)

sev1.rao <- lm(RaoQ.y~RaoQ.x, data = sev1.comb)
summary(sev1.rao)

#SEV2
sev2.FRic <- lm(FRic_sc~FRic, data = sev2.comb)
summary(sev2.FRic)

sev2.kde.alpha <- lm(kde.alpha_sc~kde.alpha, data = sev2.comb)
summary(sev2.kde.alpha)

sev2.FEve <- lm(FEve_sc~FEve, data = sev2.comb)
summary(sev2.FEve)

sev2.kde.eve <- lm(kde.evenness_sc~kde.evenness, data = sev2.comb)
summary(sev2.kde.eve)

sev2.FDis <- lm(FDis_sc~FDis, data = sev2.comb)
summary(sev2.FDis)

sev2.kde.dis <- lm(kde.dispersion_sc~kde.dispersion, data = sev2.comb)
summary(sev2.kde.dis)

sev2.FDiv <- lm(FDiv_sc~FDiv, data = sev2.comb)
summary(sev2.FDiv)

sev2.rao <- lm(RaoQ_sc~RaoQ, data = sev2.comb)
summary(sev2.rao)



###EUCLIDEAN
#read in scaled
cdr.1.sc <- read.csv(here("data/Cleaned/cdr1_sc_euc.csv"))%>%
  subset(SR != "NA")
cdr.2.sc <- read.csv(here("data/Cleaned/cdr2_sc_euc.csv"))%>%
  subset(SR != "NA")
cdr.3.sc <- read.csv(here("data/Cleaned/cdr3_sc_euc.csv"))%>%
  subset(SR != "NA")
cdr.4.sc <- read.csv(here("data/Cleaned/cdr4_sc_euc.csv"))%>%
  subset(SR != "NA")
sev.blue.sc <- read.csv(here("data/Cleaned/sevblue_sc_euc.csv"))
#sev.blue <- sev.blue[-which(sev.blue$SR == 1),]
sev.black.sc <- read.csv(here("data/Cleaned/sevblack_sc_euc.csv"))
colnames(sev.black.sc) <- paste(colnames(sev.black.sc),"sc",sep="_") 
#sev.black <- sev.black[-which(sev.black$SR ==1),]


#read in unscaled data
cdr.1 <- read.csv(here("data/Cleaned/cdr1_euc.csv"))%>%
  subset(SR != "NA")
cdr.2 <- read.csv(here("data/Cleaned/cdr2_euc.csv"))%>%
  subset(SR != "NA")
cdr.3 <- read.csv(here("data/Cleaned/cdr3_euc.csv"))%>%
  subset(SR != "NA")
cdr.4 <- read.csv(here("data/Cleaned/cdr4_euc.csv"))%>%
  subset(SR != "NA")
sev.blue <- read.csv(here("data/Cleaned/sevblue_euc.csv"))
#sev.blue <- sev.blue[-which(sev.blue$SR == 1),]
sev.black <- read.csv(here("data/Cleaned/sevblack_euc.csv"))
#sev.black <- sev.black[-which(sev.black$SR ==1),]


##merge the two
cdr1.comb <- left_join(cdr.1, cdr.1.sc, by = c("Plot", "Ring", "n_trait", "traits"))
cdr2.comb <- left_join(cdr.2, cdr.2.sc, by = c("Plot", "Ring", "n_trait", "traits"))
cdr3.comb <- left_join(cdr.3, cdr.3.sc, by = c("Plot", "Ring", "n_trait", "traits"))
cdr4.comb <- left_join(cdr.4, cdr.4.sc, by = c("Plot", "Ring", "n_trait", "traits"))
sev1.comb <- left_join(sev.blue, sev.blue.sc, by = c("sev.blueplots", "n_trait", "traits"))
sev2.comb <- cbind(sev.black, sev.black.sc)


#correlation between the two by community
#CDR1
cdr1.FRic <- lm(FRic.y~FRic.x, data = cdr1.comb)
summary(cdr1.FRic)

cdr1.kde.alpha <- lm(kde.alpha.y~kde.alpha.x, data = cdr1.comb)
summary(cdr1.kde.alpha)

cdr1.FEve <- lm(FEve.y~FEve.x, data = cdr1.comb)
summary(cdr1.FEve)

cdr1.kde.eve <- lm(kde.evenness.y~kde.evenness.x, data = cdr1.comb)
summary(cdr1.kde.eve)

cdr1.FDis <- lm(FDis.y~FDis.x, data = cdr1.comb)
summary(cdr1.FDis)

cdr1.kde.dis <- lm(kde.dispersion.y~kde.dispersion.x, data = cdr1.comb)
summary(cdr1.kde.dis)

cdr1.FDiv <- lm(FDiv.y~FDiv.x, data = cdr1.comb)
summary(cdr1.FDiv)

cdr1.rao <- lm(RaoQ.y~RaoQ.x, data = cdr1.comb)
summary(cdr1.rao)

#CDR2

cdr2.FRic <- lm(FRic.y~FRic.x, data = cdr2.comb)
summary(cdr2.FRic)

cdr2.kde.alpha <- lm(kde.alpha.y~kde.alpha.x, data = cdr2.comb)
summary(cdr2.kde.alpha)

cdr2.FEve <- lm(FEve.y~FEve.x, data = cdr2.comb)
summary(cdr2.FEve)

cdr2.kde.eve <- lm(kde.evenness.y~kde.evenness.x, data = cdr2.comb)
summary(cdr2.kde.eve)

cdr2.FDis <- lm(FDis.y~FDis.x, data = cdr2.comb)
summary(cdr2.FDis)

cdr2.kde.dis <- lm(kde.dispersion.y~kde.dispersion.x, data = cdr2.comb)
summary(cdr2.kde.dis)

cdr2.FDiv <- lm(FDiv.y~FDiv.x, data = cdr2.comb)
summary(cdr2.FDiv)

cdr2.rao <- lm(RaoQ.y~RaoQ.x, data = cdr2.comb)
summary(cdr2.rao)

#CDR3
cdr3.FRic <- lm(FRic.y~FRic.x, data = cdr3.comb)
summary(cdr3.FRic)

cdr3.kde.alpha <- lm(kde.alpha.y~kde.alpha.x, data = cdr3.comb)
summary(cdr3.kde.alpha)

cdr3.FEve <- lm(FEve.y~FEve.x, data = cdr3.comb)
summary(cdr3.FEve)

cdr3.kde.eve <- lm(kde.evenness.y~kde.evenness.x, data = cdr3.comb)
summary(cdr3.kde.eve)

cdr3.FDis <- lm(FDis.y~FDis.x, data = cdr3.comb)
summary(cdr3.FDis)

cdr3.kde.dis <- lm(kde.dispersion.y~kde.dispersion.x, data = cdr3.comb)
summary(cdr3.kde.dis)

cdr3.FDiv <- lm(FDiv.y~FDiv.x, data = cdr3.comb)
summary(cdr3.FDiv)

cdr3.rao <- lm(RaoQ.y~RaoQ.x, data = cdr3.comb)
summary(cdr3.rao)

#CDR4
cdr4.FRic <- lm(FRic.y~FRic.x, data = cdr4.comb)
summary(cdr4.FRic)

cdr4.kde.alpha <- lm(kde.alpha.y~kde.alpha.x, data = cdr4.comb)
summary(cdr4.kde.alpha)

cdr4.FEve <- lm(FEve.y~FEve.x, data = cdr4.comb)
summary(cdr4.FEve)

cdr4.kde.eve <- lm(kde.evenness.y~kde.evenness.x, data = cdr4.comb)
summary(cdr4.kde.eve)

cdr4.FDis <- lm(FDis.y~FDis.x, data = cdr4.comb)
summary(cdr4.FDis)

cdr4.kde.dis <- lm(kde.dispersion.y~kde.dispersion.x, data = cdr4.comb)
summary(cdr4.kde.dis)

cdr4.FDiv <- lm(FDiv.y~FDiv.x, data = cdr4.comb)
summary(cdr4.FDiv)

cdr4.rao <- lm(RaoQ.y~RaoQ.x, data = cdr4.comb)
summary(cdr4.rao)


#SEV1
sev1.FRic <- lm(FRic.y~FRic.x, data = sev1.comb)
summary(sev1.FRic)

sev1.kde.alpha <- lm(kde.alpha.y~kde.alpha.x, data = sev1.comb)
summary(sev1.kde.alpha)

sev1.FEve <- lm(FEve.y~FEve.x, data = sev1.comb)
summary(sev1.FEve)

sev1.kde.eve <- lm(kde.evenness.y~kde.evenness.x, data = sev1.comb)
summary(sev1.kde.eve)

sev1.FDis <- lm(FDis.y~FDis.x, data = sev1.comb)
summary(sev1.FDis)

sev1.kde.dis <- lm(kde.dispersion.y~kde.dispersion.x, data = sev1.comb)
summary(sev1.kde.dis)

sev1.FDiv <- lm(FDiv.y~FDiv.x, data = sev1.comb)
summary(sev1.FDiv)

sev1.rao <- lm(RaoQ.y~RaoQ.x, data = sev1.comb)
summary(sev1.rao)

#SEV2
sev2.FRic <- lm(FRic_sc~FRic, data = sev2.comb)
summary(sev2.FRic)

sev2.kde.alpha <- lm(kde.alpha_sc~kde.alpha, data = sev2.comb)
summary(sev2.kde.alpha)

sev2.FEve <- lm(FEve_sc~FEve, data = sev2.comb)
summary(sev2.FEve)

sev2.kde.eve <- lm(kde.evenness_sc~kde.evenness, data = sev2.comb)
summary(sev2.kde.eve)

sev2.FDis <- lm(FDis_sc~FDis, data = sev2.comb)
summary(sev2.FDis)

sev2.kde.dis <- lm(kde.dispersion_sc~kde.dispersion, data = sev2.comb)
summary(sev2.kde.dis)

sev2.FDiv <- lm(FDiv_sc~FDiv, data = sev2.comb)
summary(sev2.FDiv)

sev2.rao <- lm(RaoQ_sc~RaoQ, data = sev2.comb)
summary(sev2.rao)

