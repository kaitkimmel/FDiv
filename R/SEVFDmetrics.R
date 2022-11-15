################################
#### FD METRIC CALCULATIONS ####
################################

# CREATED BY: KAITLIN KIMMEL

# load libraries
library(FD)
library(here)
source("http://www.sthda.com/upload/rquery_cormat.r")



## LOAD DATA

sev.trait <- read.csv(here("data/Cleaned/sev_traitsout.csv"))
sev.blue <- read.csv(here("data/Cleaned/sevblue_commout.csv"))
sev.black <- read.csv(here("data/Cleaned/sevblack_commout.csv"))


### CHECK TRAIT COVERAGE 
# "Before analysis, we will remove species that have less than 100% trait coverage. 
# We will, however, make sure that the communities are still represented by at 
# least 80% of species abundance" 

# remove species with > 100% trait coverage
traits.comp <- sev.trait[complete.cases(sev.trait),] #26 species without all trait data

# match traits with species in communities
sev.bluesub <- sev.blue[,which(names(sev.blue) %in% traits.comp$Kartez)]
sev.bluesub$comp.tot <- rowSums(sev.bluesub)
min(sev.bluesub$comp.tot) # min is 79.96 close enough to 80% to inlcude. 

sev.blacksub <- sev.black[,which(names(sev.black) %in% traits.comp$Kartez)]
sev.blacksub$comp.tot <- rowSums(sev.blacksub)
min(sev.blacksub$comp.tot) # 77.72 - we will remove this row

# merge plot data back with reduced community
sev.blacksub$plotID <- sev.black$plotID
sev.blacksub <- sev.blacksub[-which(sev.blacksub$comp.tot < 80),]

# get the species with cover in the community (get rid of species with all 0 cover)
sev.blackplots <- sev.blacksub[,c(39)]
sev.black2 <- sev.blacksub[,-c(38,39)]
sev.blacksp <- names(sev.black2[which(colSums(sev.black2) != 0)])
sev.black2 <- sev.black2[,which(colSums(sev.black2) != 0)]
sev.black2 <- sev.black2[ , order(names(sev.black2))]
# pull the traits associated with each species
sev.blacktr <- sev.trait[which(sev.trait$Kartez %in% sev.blacksp),]
row.names(sev.blacktr) <- sev.blacktr$Kartez
sev.blacktr <- sev.blacktr[,-c(1:3)]

# get all combinations of traits
trait_comb_list <- list()
for (i in 2:10){
  trait_comb_list[[i-1]] <- combn(sev.blacktr[,c(1:10)], i, simplify = FALSE)
}

# Create dataframe to store metrics
df.out1 <- data.frame(SR = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA,
                      sev.blackplots = NA, n_trait = NA, traits = NA, mean_cor = NA, 
                      min_cor = NA, max_cor = NA)
# Loop to run through different trait combinations
for(j in 1:length(trait_comb_list)){
  focal_list <- trait_comb_list[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = sev.black2
    out <- dbFD(x, a, m = 2) # calculating metrics m = 2
    if(is.null(out$FDiv)){
      out$FDiv = rep(NA,29)
    }
    temp <- data.frame(SR = out$nbsp, FRic = out$FRic, FEve = out$FEve, FDiv = out$FDiv,
                       FDis = out$FDis, RaoQ = out$RaoQ)
    temp <- cbind(temp, sev.blackplots)
    temp$n_trait = ncol(focal_list[[i]])
    temp$traits = i
    if(ncol(focal_list[[i]]) == 4){
      if(is.numeric(focal_list[[i]][,c(1:4)])){
        temp.cor <- rquery.cormat(focal_list[[i]], type="flatten", graph=FALSE, method = 'spearman')
        temp$mean_cor <- mean(temp.cor$r$cor) # USE ABSOLUTE VALUES? 
        temp$min_cor <- min(temp.cor$r$cor)
        temp$max_cor <- max(temp.cor$r$cor)
      }
      else {
        temp$mean_cor <- NA
        temp$min_cor <- NA
        temp$max_cor <- NA
      }
    } else {
      temp$mean_cor <- NA
      temp$min_cor <- NA
      temp$max_cor <- NA
    }
    df.out1 <- rbind(df.out1, temp)
  }
}
