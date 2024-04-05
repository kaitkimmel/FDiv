################################
#### FD METRIC CALCULATIONS ####
################################

# CREATED BY: KAITLIN KIMMEL

# load libraries
library(BAT)
library(FD)
library(here)
source("http://www.sthda.com/upload/rquery_cormat.r")



## LOAD DATA

sev.trait <- read.csv(here("data/Cleaned/sev_traitsout.csv"))
sev.blue <- read.csv(here("data/Cleaned/sevblue_commout.csv"))
sev.blue <- subset(sev.blue, BOGR2 != 100 & plotID != "1 N 3")#remove rows with single species. necessary to calculate metrics
sev.black <- read.csv(here("data/Cleaned/sevblack_commout.csv"))
sev.black <- subset(sev.black, BOER4 != 100)#remove rows with single species. necessary to calculate metrics

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
### BLACK ####
sev.blackplots <- sev.blacksub[,c(39)]
sev.black2 <- sev.blacksub[,-c(38,39)]
sev.blacksp <- names(sev.black2[which(colSums(sev.black2) != 0)])
sev.black2 <- sev.black2[,which(colSums(sev.black2) != 0)]
sev.black2 <- sev.black2[ , order(names(sev.black2))]
# pull the traits associated with each species
sev.blacktr <- traits.comp[which(traits.comp$Kartez %in% sev.blacksp),]
row.names(sev.blacktr) <- sev.blacktr$Kartez
sev.blacktr <- sev.blacktr[,-c(1:3)]

### BLUE ####
sev.blue2 <- sev.blue[,-c(1,2)]
sev.bluesp <- names(sev.blue2[which(colSums(sev.blue2) != 0)]) # get species names to pull from trait data
sev.blue2 <- sev.blue2[,which(colSums(sev.blue2) != 0)] # remove species from community data
sev.blue2 <- sev.blue2[ ,order(names(sev.blue2))] # alphabetical order for community data
# pull the traits associated with each species
sev.bluetr <- traits.comp[which(traits.comp$Kartez %in% sev.bluesp),]
row.names(sev.bluetr) <- sev.bluetr$Kartez
sev.bluetr <- sev.bluetr[,-c(1:3)]
# some species in community that are not in trait data - need to remove them from community data
sev.bluesp1 <- row.names(sev.bluetr)
# remove those species from the community
sev.blue2 <- sev.blue2[,which(names(sev.blue2) %in% sev.bluesp1)]

#########################################################################
######### CREATING TRAIT DATAFRAMES FOR USE IN ANALYSIS ################
########################################################################
# get all combinations of traits
trait_comb_list <- list()
for (i in 2:10){
  trait_comb_list[[i-1]] <- combn(sev.blacktr[,c(1:10)], i, simplify = FALSE)
}

#### BLACK ####
# Create dataframe to store metrics
df.outblack <- data.frame(SR = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA, kde.alpha = NA, kde.evenness = NA, kde.dispersion = NA, sev.blackplots = NA, n_trait = NA, traits = NA, mean_cor = NA, min_cor = NA, max_cor = NA)


# Loop to run through different trait combinations
for(j in 1:length(trait_comb_list)){
  focal_list <- trait_comb_list[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = sev.black2
    out <- dbFD(x, a, m = 2) # calculating metrics m = 2
    if(is.null(out$FDiv)){
      out$FDiv = rep(NA,27)
    }
    #####This is where Tim is going to do KDE stuff
    temp.to <- kernel.build(comm = a, trait = focal_list[[i]],abund = TRUE, distance = "gower", axes = 2)
    kde.alpha <- kernel.alpha(temp.to)
    kde.alpha <- data.frame(kde.alpha)
    
    kde.evenness <- kernel.evenness(temp.to)
    kde.evenness <- data.frame(kde.evenness)
    
    kde.dispersion <- kernel.dispersion(temp.to)
    kde.dispersion <- data.frame(kde.dispersion)
    
    #####END TIM'S BIT
    temp <- data.frame(SR = out$nbsp, FRic = out$FRic, FEve = out$FEve, FDiv = out$FDiv,
                       FDis = out$FDis, RaoQ = out$RaoQ, kde.alpha = kde.alpha$kde.alpha,
                       kde.evenness = kde.evenness$kde.evenness, kde.dispersion = kde.dispersion$kde.dispersion)
    temp <- cbind(temp, sev.blackplots)
    temp$n_trait = ncol(focal_list[[i]])
    temp$traits = i
    if(ncol(focal_list[[i]]) == 4){
      if(is.numeric(focal_list[[i]][,c(1)])==TRUE&is.numeric(focal_list[[i]][,c(2)])==TRUE&is.numeric(focal_list[[i]][,c(3)])==TRUE&is.numeric(focal_list[[i]][,c(4)])==TRUE){
        temp.cor <- rquery.cormat(focal_list[[i]], type="flatten", graph=FALSE, method = 'spearman')
        temp$mean_cor <- mean(abs(temp.cor$r$cor)) # USE ABSOLUTE VALUES? 
        temp$min_cor <- min(abs(temp.cor$r$cor))
        temp$max_cor <- max(abs(temp.cor$r$cor))
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
    df.outblack <- rbind(df.outblack, temp)
  }
}

df.outblack <- df.outblack[-which(is.na(df.outblack)),]

#### BLUE ####
# get all combinations of traits
trait_comb_list1 <- list()
for (i in 2:10){
  trait_comb_list1[[i-1]] <- combn(sev.bluetr[,c(1:10)], i, simplify = FALSE)
}
# Create dataframe to store metrics
df.outblue <- data.frame(SR = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA, 
                         kde.alpha = NA, kde.evenness = NA, kde.dispersion = NA, sev.blueplots = NA, 
                         n_trait = NA, traits = NA, mean_cor = NA, min_cor = NA, max_cor = NA)
# get plots
sev.blueplots <- sev.blue[,2]
# Loop to run through different trait combinations
for(j in 1:length(trait_comb_list1)){
  focal_list <- trait_comb_list1[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = sev.blue2
    out <- dbFD(x, a, m = 2) # calculating metrics m = 2
    if(is.null(out$FDiv)){
      out$FDiv = rep(NA,28)
    }
    #####This is where Tim is going to start trying out KDE stuff
    temp.to <- kernel.build(comm = a, trait = focal_list[[i]],abund = TRUE, distance = "gower", axes = 2)
    kde.alpha <- kernel.alpha(temp.to)
    kde.alpha <- data.frame(kde.alpha)
    
    kde.evenness <- kernel.evenness(temp.to)
    kde.evenness <- data.frame(kde.evenness)
    
    kde.dispersion <- kernel.dispersion(temp.to)
    kde.dispersion <- data.frame(kde.dispersion)
    
    #####END TIM'S EXPERIMENT
    temp <- data.frame(SR = out$nbsp, FRic = out$FRic, FEve = out$FEve, FDiv = out$FDiv,
                       FDis = out$FDis, RaoQ = out$RaoQ, kde.alpha = kde.alpha$kde.alpha, 
                       kde.evenness = kde.evenness$kde.evenness, kde.dispersion = kde.dispersion$kde.dispersion)
    temp <- cbind(temp, sev.blueplots)
    temp$n_trait = ncol(focal_list[[i]])
    temp$traits = i
    if(ncol(focal_list[[i]]) == 4){
      if(is.numeric(focal_list[[i]][,c(1)])==TRUE&is.numeric(focal_list[[i]][,c(2)])==TRUE&is.numeric(focal_list[[i]][,c(3)])==TRUE&is.numeric(focal_list[[i]][,c(4)])==TRUE){
        temp.cor <- rquery.cormat(focal_list[[i]], type="flatten", graph=FALSE, method = 'spearman')
        temp$mean_cor <- mean(abs(temp.cor$r$cor)) # USE ABSOLUTE VALUES? 
        temp$min_cor <- min(abs(temp.cor$r$cor))
        temp$max_cor <- max(abs(temp.cor$r$cor))
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
    df.outblue <- rbind(df.outblue, temp)
  }
}

df.outblue <- df.outblue[-which(is.na(df.outblue)),]


#########################################
##### SENSITIVITY ANALYSES ##############
########################################
trait.list <- trait_comb_list[c(3:9)]


#### BLACK ####
# Create dataframe to store metrics
df.outblack.sen <- data.frame(FRich = NA, m = NA, sev.blackplots = NA, n_trait = NA, traits = NA, SR = NA)


# Loop to run through different trait combinations
for(j in 1:length(trait.list)){
  focal_list <- trait.list[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = sev.black2
    for(k in c(2:4)){
      out <- dbFD(x, a, m = k) 
      temp <- data.frame(SR = out$nbsp, FRich = out$FRic)
      temp <- cbind(temp, sev.blackplots)
      temp$n_trait = ncol(focal_list[[i]])
      temp$traits = i
      temp$m = k
      df.outblack.sen <- rbind(df.outblack.sen, temp)
    }
    
  }
}

df.outblack.sen <- df.outblack.sen[-which(is.na(df.outblack.sen)),]


#### BLUE ####
trait.list1 <- trait_comb_list1[c(3:9)]
# Create dataframe to store metrics
df.outblue.sen <- data.frame(FRich = NA, m = NA, sev.blueplots = NA, n_trait = NA, traits = NA, SR = NA)


# Loop to run through different trait combinations
for(j in 1:length(trait.list1)){
  focal_list <- trait.list1[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = sev.blue2
    for(k in c(2:4)){
      out <- dbFD(x, a, m = k) 
      temp <- data.frame(SR = out$nbsp, FRich = out$FRic)
      temp <- cbind(temp, sev.blueplots)
      temp$n_trait = ncol(focal_list[[i]])
      temp$traits = i
      temp$m = k
      df.outblue.sen <- rbind(df.outblue.sen, temp)
    }
    
  }
}

df.outblue.sen <- df.outblue.sen[-which(is.na(df.outblue.sen)),]

#### SAVE OUTPUT FOR ANALYSIS

write.csv(df.outblack, here("data/Cleaned/sevblack.csv"), row.names = FALSE)
write.csv(df.outblue, here("data/Cleaned/sevblue.csv"), row.names = FALSE)

write.csv(df.outblack.sen, here("data/Cleaned/sevblack_sen.csv"), row.names = FALSE)
write.csv(df.outblue.sen, here("data/Cleaned/sevblue_sen.csv"), row.names = FALSE)

###############################################################################
# SCALE AND CENTER ANALYSIS : re: R1 comments

sev.blacktr.sc <- as.data.frame(scale(sev.blacktr[c(1:8,10)], center = TRUE, scale = TRUE))
sev.blacktr.sc <- cbind(sev.blacktr$PhotoPath, sev.blacktr.sc)
names(sev.blacktr.sc)[1] <- 'PhotoPath'

sev.bluetr.sc <- as.data.frame(scale(sev.bluetr[,c(1:8,10)], center = TRUE, scale = TRUE))
sev.bluetr.sc <- cbind(sev.bluetr$PhotoPath, sev.bluetr.sc)
names(sev.bluetr.sc)[1] <- 'PhotoPath'

### Make trait df into one list 
trait_list_SC <- list(sev.blacktr.sc, sev.bluetr.sc)


#############################################################
######## FD METRIC CALCULATION #############################
###########################################################
### COMMUNITY 1
### Get all combinations of traits for community 1
trait_comb_list_sc <- list()
for (i in 2:10){
  trait_comb_list_sc[[i-1]] <- combn(trait_list_SC[[1]], i, simplify = FALSE)
}
# Create dataframe to store metrics
df.outblack.sc <- data.frame(SR = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA, 
                            kde.alpha = NA, kde.evenness = NA, kde.dispersion = NA, sev.blackplots = NA, 
                            n_trait = NA, traits = NA, mean_cor = NA, min_cor = NA, max_cor = NA)
# get plots
sev.blackplots <- sev.black2[,2]
# Loop to run through different trait combinations
for(j in 1:length(trait_comb_list_sc)){
  focal_list <- trait_comb_list_sc[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = sev.black2
    out <- dbFD(x, a, m = 2) # calculating metrics m = 2
    if(is.null(out$FDiv)){
      out$FDiv = rep(NA,27)
    }
    #####This is where Tim is going to start trying out KDE stuff
    temp.to <- kernel.build(comm = a, trait = focal_list[[i]], abund = TRUE, distance = "gower", axes = 2)
    kde.alpha <- kernel.alpha(temp.to)
    kde.alpha <- data.frame(kde.alpha)
    
    kde.evenness <- kernel.evenness(temp.to)
    kde.evenness <- data.frame(kde.evenness)
    
    kde.dispersion <- kernel.dispersion(temp.to)
    kde.dispersion <- data.frame(kde.dispersion)
    
    #####END TIM'S EXPERIMENT
    temp <- data.frame(SR = out$nbsp, FRic = out$FRic, FEve = out$FEve, FDiv = out$FDiv,
                       FDis = out$FDis, RaoQ = out$RaoQ, kde.alpha = kde.alpha$kde.alpha, 
                       kde.evenness = kde.evenness$kde.evenness, kde.dispersion = kde.dispersion$kde.dispersion)
    temp <- cbind(temp, sev.blackplots)
    temp$n_trait = ncol(focal_list[[i]])
    temp$traits = i
    if(ncol(focal_list[[i]]) == 4){
      if(is.numeric(focal_list[[i]][,c(1)])==TRUE&is.numeric(focal_list[[i]][,c(2)])==TRUE&is.numeric(focal_list[[i]][,c(3)])==TRUE&is.numeric(focal_list[[i]][,c(4)])==TRUE){
        temp.cor <- rquery.cormat(focal_list[[i]], type="flatten", graph=FALSE, method = 'spearman')
        temp$mean_cor <- mean(abs(temp.cor$r$cor)) # USE ABSOLUTE VALUES? 
        temp$min_cor <- min(abs(temp.cor$r$cor))
        temp$max_cor <- max(abs(temp.cor$r$cor))
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
    df.outblack.sc <- rbind(df.outblack.sc, temp)
  }
}

df.outblack.sc <- df.outblack.sc[-which(is.na(df.outblack.sc)),]

#df.out1.sc <- df.out1.sc[-which(is.na(df.out1.sc)),]
write.csv(df.outblack.sc, here("data/Cleaned/sevblack_sc.csv"), row.names = FALSE)
#### COMMUNITY 2

trait_comb_list1_sc <- list()
for (i in 2:10){
  trait_comb_list1_sc[[i-1]] <- combn(trait_list_SC[[2]], i, simplify = FALSE)
}
# Create dataframe to store metrics
df.outblue.sc <- data.frame(SR = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA, 
                         kde.alpha = NA, kde.evenness = NA, kde.dispersion = NA, sev.blueplots = NA, 
                         n_trait = NA, traits = NA, mean_cor = NA, min_cor = NA, max_cor = NA)
# get plots
sev.blueplots <- sev.blue[,2]
# Loop to run through different trait combinations
for(j in 1:length(trait_comb_list1_sc)){
  focal_list <- trait_comb_list1_sc[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = sev.blue2
    out <- dbFD(x, a, m = 2) # calculating metrics m = 2
    if(is.null(out$FDiv)){
      out$FDiv = rep(NA,28)
    }
    #####This is where Tim is going to start trying out KDE stuff
    temp.to <- kernel.build(comm = a, trait = focal_list[[i]], abund = TRUE, distance = "gower", axes = 2)
    kde.alpha <- kernel.alpha(temp.to)
    kde.alpha <- data.frame(kde.alpha)
    
    kde.evenness <- kernel.evenness(temp.to)
    kde.evenness <- data.frame(kde.evenness)
    
    kde.dispersion <- kernel.dispersion(temp.to)
    kde.dispersion <- data.frame(kde.dispersion)
    
    #####END TIM'S EXPERIMENT
    temp <- data.frame(SR = out$nbsp, FRic = out$FRic, FEve = out$FEve, FDiv = out$FDiv,
                       FDis = out$FDis, RaoQ = out$RaoQ, kde.alpha = kde.alpha$kde.alpha, 
                       kde.evenness = kde.evenness$kde.evenness, kde.dispersion = kde.dispersion$kde.dispersion)
    temp <- cbind(temp, sev.blueplots)
    temp$n_trait = ncol(focal_list[[i]])
    temp$traits = i
    if(ncol(focal_list[[i]]) == 4){
      if(is.numeric(focal_list[[i]][,c(1)])==TRUE&is.numeric(focal_list[[i]][,c(2)])==TRUE&is.numeric(focal_list[[i]][,c(3)])==TRUE&is.numeric(focal_list[[i]][,c(4)])==TRUE){
        temp.cor <- rquery.cormat(focal_list[[i]], type="flatten", graph=FALSE, method = 'spearman')
        temp$mean_cor <- mean(abs(temp.cor$r$cor)) # USE ABSOLUTE VALUES? 
        temp$min_cor <- min(abs(temp.cor$r$cor))
        temp$max_cor <- max(abs(temp.cor$r$cor))
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
    df.outblue.sc <- rbind(df.outblue.sc, temp)
  }
}

df.outblue.sc <- df.outblue.sc[-which(is.na(df.outblue.sc)),]
write.csv(df.outblue.sc, here("data/Cleaned/sevblue_sc.csv"), row.names = FALSE)
