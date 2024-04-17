################################
#### FD METRIC CALCULATIONS ####
################################

# CREATED BY: KAITLIN KIMMEL

# load libraries
library(BAT)
library(FD)
library(here)
source("http://www.sthda.com/upload/rquery_cormat.r")

# load data
# Cedar Creek Traits
cdr.traits1 <- read.csv(here("data/Cleaned/cdrnaca_traitsout.csv"))
cdr.traits2 <- read.csv(here("data/Cleaned/cdrneca_traitsout.csv"))
cdr.traits3 <- read.csv(here("data/Cleaned/cdrnace_traitsout.csv"))
cdr.traits4 <- read.csv(here("data/Cleaned/cdrnece_traitsout.csv"))

# Cedar Creek Community
cdr.comm1 <- read.csv(here("data/Cleaned/cdrnaca_commout.csv"))
cdr.comm2 <- read.csv(here("data/Cleaned/cdrneca_commout.csv"))
cdr.comm3 <- read.csv(here("data/Cleaned/cdrnace_commout.csv"))
cdr.comm4 <- read.csv(here("data/Cleaned/cdrnece_commout.csv"))

###################################################################
#### COMM & TRAIT 1 #### Note: make this so that all communities for 1 site can be passed in at once

### Make community df in to one list
comm_list <- list(cdr.comm1, cdr.comm2, cdr.comm3, cdr.comm4)

### Make trait df into one list 
trait_list <- list(cdr.traits1, cdr.traits2, cdr.traits3, cdr.traits4)

enviro_list <- list()

for (i in 1:length(comm_list)){
  ### 1. get plot and ring for each community df
  enviro_list[[i]] <- data.frame(Plot = comm_list[[i]]$Plot, Ring = comm_list[[i]]$Ring)
  ### 2. get rid of unnecessary columns
  comm_list[[i]] <- comm_list[[i]][,c(5:20)]
  ### 3. Remove species not found in plots from both community and trait data
  miss.sp <- names(comm_list[[i]])[which(colSums(comm_list[[i]]) == 0)] # get list of species not in community data
  # Achillea, Agropyron, Anemone, Bouteloua, Koeleria, Petalostemum, and Solidago
  comm_list[[i]] <- comm_list[[i]][,-which(colSums(comm_list[[i]])==0)] # remove species from community data
  trait_list[[i]] <- trait_list[[i]][-which(trait_list[[i]]$monospecies %in% miss.sp),]
  ### 4. make rownames of trait list the species names
  rownames(trait_list[[i]]) <- trait_list[[i]]$monospecies
  trait_list[[i]] <- trait_list[[i]][,-c(1:3)]
}

#############################################################
######## FD METRIC CALCULATION #############################
###########################################################
### COMMUNITY 1
### Get all combinations of traits for community 1
trait_comb_list <- list()
for (i in 2:9){
  trait_comb_list[[i-1]] <- combn(trait_list[[1]], i, simplify = FALSE)
}

# Create dataframe to store metrics
df.out1 <- data.frame(SR = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA, kde.alpha = NA, kde.evenness = NA, kde.dispersion = NA, Plot = NA, Ring = NA, n_trait = NA, traits = NA, mean_cor = NA, min_cor = NA, max_cor = NA)

# Loop to run through different trait combinations
for(j in 1:length(trait_comb_list)){
  focal_list <- trait_comb_list[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = comm_list[[1]]
    out <- dbFD(x, a, m = 2) # calculating metrics m = 2
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
                       FDis = out$FDis, RaoQ = out$RaoQ, kde.alpha = kde.alpha$kde.alpha, kde.evenness = kde.evenness$kde.evenness, kde.dispersion = kde.dispersion$kde.dispersion)    
    temp <- cbind(temp, enviro_list[[1]])
    temp$n_trait = ncol(focal_list[[i]])
    temp$traits = i
    if(ncol(focal_list[[i]]) == 4){
      temp.cor <- rquery.cormat(focal_list[[i]], type="flatten", graph=FALSE, method = 'spearman')
      temp$mean_cor <- mean(abs(temp.cor$r$cor)) # USE ABSOLUTE VALUES? 
      temp$min_cor <- min(abs(temp.cor$r$cor))
      temp$max_cor <- max(abs(temp.cor$r$cor))
    } else {
      temp$mean_cor <- NA
      temp$min_cor <- NA
      temp$max_cor <- NA
    }
    df.out1 <- rbind(df.out1, temp)
  }
}

df.out1 <- df.out1[-which(is.na(df.out1)),]


#### COMMUNITY 2
trait_comb_list2 <- list()
for (i in 2:9){
  trait_comb_list2[[i-1]] <- combn(trait_list[[2]], i, simplify = FALSE)
}


df.out2 <- data.frame(SR = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA, kde.alpha = NA, kde.evenness = NA, kde.dispersion = NA, Plot = NA, Ring = NA, n_trait = NA, traits = NA, mean_cor = NA, min_cor = NA, max_cor = NA)

for(j in 1:length(trait_comb_list2)){
  focal_list <- trait_comb_list2[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = comm_list[[2]]
    out <- dbFD(x, a, m = 2)
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
                       FDis = out$FDis, RaoQ = out$RaoQ, kde.alpha = kde.alpha$kde.alpha, kde.evenness = kde.evenness$kde.evenness, kde.dispersion = kde.dispersion$kde.dispersion) 
    temp <- cbind(temp, enviro_list[[2]])
    temp$n_trait = ncol(focal_list[[i]])
    temp$traits = i
    if(ncol(focal_list[[i]]) == 4){
      temp.cor <- rquery.cormat(focal_list[[i]], type="flatten", graph=FALSE, method = 'spearman')
      temp$mean_cor <- mean(abs(temp.cor$r$cor)) # USE ABSOLUTE VALUES? 
      temp$min_cor <- min(abs(temp.cor$r$cor))
      temp$max_cor <- max(abs(temp.cor$r$cor))
    } else {
      temp$mean_cor <- NA
      temp$min_cor <- NA
      temp$max_cor <- NA
    }
    df.out2 <- rbind(df.out2, temp)
  }
}

df.out2 <- df.out2[-which(is.na(df.out2)),]

### COMMUNITY 3
### Get all combinations of traits for community 3
trait_comb_list3 <- list()
for (i in 2:9){
  trait_comb_list3[[i-1]] <- combn(trait_list[[3]], i, simplify = FALSE)
}

# Create dataframe to store metrics
df.out3 <- data.frame(SR = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA, kde.alpha = NA, kde.evenness = NA, kde.dispersion = NA,  Plot = NA, Ring = NA, n_trait = NA, traits = NA, mean_cor = NA, min_cor = NA, max_cor = NA)
# Loop to run through different trait combinations
for(j in 1:length(trait_comb_list3)){
  focal_list <- trait_comb_list3[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = comm_list[[3]]
    out <- dbFD(x, a, m = 2) # calculating metrics m = 2
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
                       FDis = out$FDis, RaoQ = out$RaoQ, kde.alpha = kde.alpha$kde.alpha, kde.evenness = kde.evenness$kde.evenness, kde.dispersion = kde.dispersion$kde.dispersion) 
    temp <- cbind(temp, enviro_list[[3]])
    temp$n_trait = ncol(focal_list[[i]])
    temp$traits = i
    if(ncol(focal_list[[i]]) == 4){
      temp.cor <- rquery.cormat(focal_list[[i]], type="flatten", graph=FALSE, method = 'spearman')
      temp$mean_cor <- mean(abs(temp.cor$r$cor)) # USE ABSOLUTE VALUES? 
      temp$min_cor <- min(abs(temp.cor$r$cor))
      temp$max_cor <- max(abs(temp.cor$r$cor))
    } else {
      temp$mean_cor <- NA
      temp$min_cor <- NA
      temp$max_cor <- NA
    }
    df.out3 <- rbind(df.out3, temp)
  }
}

df.out3 <- df.out3[-which(is.na(df.out3)),]


### COMMUNITY 4
### Get all combinations of traits for community 4
trait_comb_list4 <- list()
for (i in 2:9){
  trait_comb_list4[[i-1]] <- combn(trait_list[[4]], i, simplify = FALSE)
}

# Create dataframe to store metrics
df.out4 <- data.frame(SR = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA, kde.alpha = NA, kde.evenness = NA, kde.dispersion = NA, Plot = NA, Ring = NA, n_trait = NA, traits = NA, mean_cor = NA, min_cor = NA, max_cor = NA)
# Loop to run through different trait combinations
for(j in 1:length(trait_comb_list4)){
  focal_list <- trait_comb_list4[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = comm_list[[4]]
    out <- dbFD(x, a, m = 2) # calculating metrics m = 2
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
                       FDis = out$FDis, RaoQ = out$RaoQ, kde.alpha = kde.alpha$kde.alpha, kde.evenness = kde.evenness$kde.evenness, kde.dispersion = kde.dispersion$kde.dispersion) 
    temp <- cbind(temp, enviro_list[[4]])
    temp$n_trait = ncol(focal_list[[i]])
    temp$traits = i
    if(ncol(focal_list[[i]]) == 4){
      temp.cor <- rquery.cormat(focal_list[[i]], type="flatten", graph=FALSE, method = 'spearman')
      temp$mean_cor <- mean(abs(temp.cor$r$cor)) # USE ABSOLUTE VALUES? 
      temp$min_cor <- min(abs(temp.cor$r$cor))
      temp$max_cor <- max(abs(temp.cor$r$cor))
    } else {
      temp$mean_cor <- NA
      temp$min_cor <- NA
      temp$max_cor <- NA
    }
    df.out4 <- rbind(df.out4, temp)
  }
}

df.out4 <- df.out4[-which(is.na(df.out4)),]
###

#### SENSTITIVY ANALYSES ####

trait.list1 <- trait_comb_list[c(3:8)]
# Create dataframe to store metrics
df.out1.sen <- data.frame(FRich = NA, m = NA, Ring = NA, Plot = NA, n_trait = NA, traits = NA, SR = NA)


# Loop to run through different trait combinations
for(j in 1:length(trait.list1)){
  focal_list <- trait.list1[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = comm_list[[1]]
    for(k in c(2:4)){
      out <- dbFD(x, a, m = k) 
      temp <- data.frame(SR = out$nbsp, FRich = out$FRic)
      temp <- cbind(temp, enviro_list[[1]])
      temp$n_trait = ncol(focal_list[[i]])
      temp$traits = i
      temp$m = k
      df.out1.sen <- rbind(df.out1.sen, temp)
    }
    
  }
}

df.out1.sen <- df.out1.sen[-which(is.na(df.out1.sen)),]

trait.list2 <- trait_comb_list2[c(3:8)]
# Create dataframe to store metrics
df.out2.sen <- data.frame(FRich = NA, m = NA, Ring = NA, Plot = NA, n_trait = NA, traits = NA, SR = NA)


# Loop to run through different trait combinations
for(j in 1:length(trait.list2)){
  focal_list <- trait.list2[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = comm_list[[2]]
    for(k in c(2:4)){
      out <- dbFD(x, a, m = k) 
      temp <- data.frame(SR = out$nbsp, FRich = out$FRic)
      temp <- cbind(temp, enviro_list[[2]])
      temp$n_trait = ncol(focal_list[[i]])
      temp$traits = i
      temp$m = k
      df.out2.sen <- rbind(df.out2.sen, temp)
    }
    
  }
}

df.out2.sen <- df.out2.sen[-which(is.na(df.out2.sen)),]

trait.list1 <- trait_comb_list1[c(3:8)]
# Create dataframe to store metrics
df.out3.sen <- data.frame(FRich = NA, m = NA, Ring = NA, Plot = NA, n_trait = NA, traits = NA, SR = NA)


# Loop to run through different trait combinations
for(j in 1:length(trait.list3)){
  focal_list <- trait.list3[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = comm_list[[3]]
    for(k in c(2:4)){
      out <- dbFD(x, a, m = k) 
      temp <- data.frame(SR = out$nbsp, FRich = out$FRic)
      temp <- cbind(temp, enviro_list[[3]])
      temp$n_trait = ncol(focal_list[[i]])
      temp$traits = i
      temp$m = k
      df.out3.sen <- rbind(df.out3.sen, temp)
    }
    
  }
}

df.out3.sen <- df.out3.sen[-which(is.na(df.out3.sen)),]

trait.list4 <- trait_comb_list4[c(3:8)]
# Create dataframe to store metrics
df.out4.sen <- data.frame(FRich = NA, m = NA, Ring = NA, Plot = NA, n_trait = NA, traits = NA, SR = NA)


# Loop to run through different trait combinations
for(j in 1:length(trait.list4)){
  focal_list <- trait.list4[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = comm_list[[4]]
    for(k in c(2:4)){
      out <- dbFD(x, a, m = k) 
      temp <- data.frame(SR = out$nbsp, FRich = out$FRic)
      temp <- cbind(temp, enviro_list[[4]])
      temp$n_trait = ncol(focal_list[[i]])
      temp$traits = i
      temp$m = k
      df.out4.sen <- rbind(df.out4.sen, temp)
    }
    
  }
}

df.out4.sen <- df.out4.sen[-which(is.na(df.out4.sen)),]
#### SAVE OUTPUT FOR ANALYSIS

write.csv(df.out1, here("data/Cleaned/cdr1.csv"), row.names = FALSE)
write.csv(df.out2, here("data/Cleaned/cdr2.csv"), row.names = FALSE)
write.csv(df.out3, here("data/Cleaned/cdr3.csv"), row.names = FALSE)
write.csv(df.out4, here("data/Cleaned/cdr4.csv"), row.names = FALSE)

write.csv(df.out1.sen, here("data/Cleaned/cdr1_sen.csv"), row.names = FALSE)
write.csv(df.out2.sen, here("data/Cleaned/cdr2_sen.csv"), row.names = FALSE)
write.csv(df.out3.sen, here("data/Cleaned/cdr3_sen.csv"), row.names = FALSE)
write.csv(df.out4.sen, here("data/Cleaned/cdr4_sen.csv"), row.names = FALSE)

###############################################################################
# SCALE AND CENTER ANALYSIS : re: R1 comments

c.traits1.sc <- as.data.frame(scale(cdr.traits1[,c(4:12)], center = TRUE, scale = TRUE))
c.traits1.sc <- cbind(cdr.traits1[,c(1:3)], c.traits1.sc)

c.traits2.sc <- as.data.frame(scale(cdr.traits2[,c(4:12)], center = TRUE, scale = TRUE))
c.traits2.sc <- cbind(cdr.traits2[,c(1:3)], c.traits2.sc)

c.traits3.sc <- as.data.frame(scale(cdr.traits3[,c(4:12)], center = TRUE, scale = TRUE))
c.traits3.sc <- cbind(cdr.traits3[,c(1:3)], c.traits3.sc)

c.traits4.sc <- as.data.frame(scale(cdr.traits4[,c(4:12)], center = TRUE, scale = TRUE))
c.traits4.sc <- cbind(cdr.traits4[,c(1:3)],c.traits4.sc)

### Make community df in to one list
comm_list <- list(cdr.comm1, cdr.comm2, cdr.comm3, cdr.comm4)

### Make trait df into one list 
trait_list_SC <- list(c.traits1.sc, c.traits2.sc, c.traits3.sc, c.traits4.sc)

enviro_list <- list()

for (i in 1:length(comm_list)){
  ### 1. get plot and ring for each community df
  enviro_list[[i]] <- data.frame(Plot = comm_list[[i]]$Plot, Ring = comm_list[[i]]$Ring)
  ### 2. get rid of unnecessary columns
  comm_list[[i]] <- comm_list[[i]][,c(5:20)]
  ### 3. Remove species not found in plots from both community and trait data
  miss.sp <- names(comm_list[[i]])[which(colSums(comm_list[[i]]) == 0)] # get list of species not in community data
  # Achillea, Agropyron, Anemone, Bouteloua, Koeleria, Petalostemum, and Solidago
  comm_list[[i]] <- comm_list[[i]][,-which(colSums(comm_list[[i]])==0)] # remove species from community data
  trait_list_SC[[i]] <- trait_list_SC[[i]][-which(trait_list_SC[[i]]$monospecies %in% miss.sp),]
  ### 4. make rownames of trait list the species names
  rownames(trait_list_SC[[i]]) <- trait_list_SC[[i]]$monospecies
  trait_list_SC[[i]] <- trait_list_SC[[i]][,-c(1:3)]
}
#############################################################
######## FD METRIC CALCULATION #############################
###########################################################
### COMMUNITY 1
### Get all combinations of traits for community 1
trait_sc_comb_list <- list()
for (i in 2:9){
  trait_sc_comb_list[[i-1]] <- combn(trait_list_SC[[1]], i, simplify = FALSE)
}

# Create dataframe to store metrics
df.out1.sc <- data.frame(SR = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA, kde.alpha = NA, kde.evenness = NA, kde.dispersion = NA, Plot = NA, Ring = NA, n_trait = NA, traits = NA, mean_cor = NA, min_cor = NA, max_cor = NA)

# Loop to run through different trait combinations
for(j in 1:length(trait_sc_comb_list)){
  focal_list <- trait_sc_comb_list[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = comm_list[[1]]
    out <- dbFD(x, a, m = 2) # calculating metrics m = 2
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
                       FDis = out$FDis, RaoQ = out$RaoQ, kde.alpha = kde.alpha$kde.alpha, kde.evenness = kde.evenness$kde.evenness, kde.dispersion = kde.dispersion$kde.dispersion)    
    temp <- cbind(temp, enviro_list[[1]])
    temp$n_trait = ncol(focal_list[[i]])
    temp$traits = i
    if(ncol(focal_list[[i]]) == 4){
      temp.cor <- rquery.cormat(focal_list[[i]], type="flatten", graph=FALSE, method = 'spearman')
      temp$mean_cor <- mean(abs(temp.cor$r$cor)) # USE ABSOLUTE VALUES? 
      temp$min_cor <- min(abs(temp.cor$r$cor))
      temp$max_cor <- max(abs(temp.cor$r$cor))
    } else {
      temp$mean_cor <- NA
      temp$min_cor <- NA
      temp$max_cor <- NA
    }
    df.out1.sc <- rbind(df.out1.sc, temp)
  }
}

#df.out1.sc <- df.out1.sc[-which(is.na(df.out1.sc)),]
write.csv(df.out1.sc, here("data/Cleaned/cdr1_sc.csv"), row.names = FALSE)
#### COMMUNITY 2
trait_sc_comb_list2 <- list()
for (i in 2:9){
  trait_sc_comb_list2[[i-1]] <- combn(trait_list_SC[[2]], i, simplify = FALSE)
}

# Create dataframe to store metrics
df.out2.sc <- data.frame(SR = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA, kde.alpha = NA, kde.evenness = NA, kde.dispersion = NA, Plot = NA, Ring = NA, n_trait = NA, traits = NA, mean_cor = NA, min_cor = NA, max_cor = NA)

# Loop to run through different trait combinations
for(j in 1:length(trait_sc_comb_list2)){
  focal_list <- trait_sc_comb_list2[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = comm_list[[2]]
    out <- dbFD(x, a, m = 2) # calculating metrics m = 2
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
                       FDis = out$FDis, RaoQ = out$RaoQ, kde.alpha = kde.alpha$kde.alpha, kde.evenness = kde.evenness$kde.evenness, kde.dispersion = kde.dispersion$kde.dispersion)    
    temp <- cbind(temp, enviro_list[[2]])
    temp$n_trait = ncol(focal_list[[i]])
    temp$traits = i
    if(ncol(focal_list[[i]]) == 4){
      temp.cor <- rquery.cormat(focal_list[[i]], type="flatten", graph=FALSE, method = 'spearman')
      temp$mean_cor <- mean(abs(temp.cor$r$cor)) # USE ABSOLUTE VALUES? 
      temp$min_cor <- min(abs(temp.cor$r$cor))
      temp$max_cor <- max(abs(temp.cor$r$cor))
    } else {
      temp$mean_cor <- NA
      temp$min_cor <- NA
      temp$max_cor <- NA
    }
    df.out2.sc <- rbind(df.out2.sc, temp)
  }
}

#df.out2.sc <- df.out2.sc[-which(is.na(df.out2.sc)),]
write.csv(df.out2.sc, here("data/Cleaned/cdr2_sc.csv"), row.names = FALSE)
### COMMUNITY 3
### Get all combinations of traits for community 3
trait_sc_comb_list3 <- list()
for (i in 2:9){
  trait_sc_comb_list3[[i-1]] <- combn(trait_list_SC[[3]], i, simplify = FALSE)
}

# Create dataframe to store metrics
df.out3.sc <- data.frame(SR = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA, kde.alpha = NA, kde.evenness = NA, kde.dispersion = NA,  Plot = NA, Ring = NA, n_trait = NA, traits = NA, mean_cor = NA, min_cor = NA, max_cor = NA)
# Loop to run through different trait combinations
for(j in 1:length(trait_sc_comb_list3)){
  focal_list <- trait_sc_comb_list3[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = comm_list[[3]]
    out <- dbFD(x, a, m = 2) # calculating metrics m = 2
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
                       FDis = out$FDis, RaoQ = out$RaoQ, kde.alpha = kde.alpha$kde.alpha, kde.evenness = kde.evenness$kde.evenness, kde.dispersion = kde.dispersion$kde.dispersion) 
    temp <- cbind(temp, enviro_list[[3]])
    temp$n_trait = ncol(focal_list[[i]])
    temp$traits = i
    if(ncol(focal_list[[i]]) == 4){
      temp.cor <- rquery.cormat(focal_list[[i]], type="flatten", graph=FALSE, method = 'spearman')
      temp$mean_cor <- mean(abs(temp.cor$r$cor)) # USE ABSOLUTE VALUES? 
      temp$min_cor <- min(abs(temp.cor$r$cor))
      temp$max_cor <- max(abs(temp.cor$r$cor))
    } else {
      temp$mean_cor <- NA
      temp$min_cor <- NA
      temp$max_cor <- NA
    }
    df.out3.sc <- rbind(df.out3.sc, temp)
  }
}

#df.out3.sc <- df.out3.sc[-which(is.na(df.out3.sc)),]
write.csv(df.out3.sc, here("data/Cleaned/cdr3_sc.csv"), row.names = FALSE)

### COMMUNITY 4
### Get all combinations of traits for community 4
trait_sc_comb_list4 <- list()
for (i in 2:9){
  trait_sc_comb_list4[[i-1]] <- combn(trait_list_SC[[4]], i, simplify = FALSE)
}

# Create dataframe to store metrics
df.out4.sc <- data.frame(SR = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA, kde.alpha = NA, kde.evenness = NA, kde.dispersion = NA, Plot = NA, Ring = NA, n_trait = NA, traits = NA, mean_cor = NA, min_cor = NA, max_cor = NA)
# Loop to run through different trait combinations
for(j in 1:length(trait_sc_comb_list4)){
  focal_list <- trait_sc_comb_list4[[j]]
  for (i in 1:length(focal_list)){
    x = gowdis(focal_list[[i]])
    a = comm_list[[4]]
    out <- dbFD(x, a, m = 2) # calculating metrics m = 2
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
                       FDis = out$FDis, RaoQ = out$RaoQ, kde.alpha = kde.alpha$kde.alpha, kde.evenness = kde.evenness$kde.evenness, kde.dispersion = kde.dispersion$kde.dispersion) 
    temp <- cbind(temp, enviro_list[[4]])
    temp$n_trait = ncol(focal_list[[i]])
    temp$traits = i
    if(ncol(focal_list[[i]]) == 4){
      temp.cor <- rquery.cormat(focal_list[[i]], type="flatten", graph=FALSE, method = 'spearman')
      temp$mean_cor <- mean(abs(temp.cor$r$cor)) # USE ABSOLUTE VALUES? 
      temp$min_cor <- min(abs(temp.cor$r$cor))
      temp$max_cor <- max(abs(temp.cor$r$cor))
    } else {
      temp$mean_cor <- NA
      temp$min_cor <- NA
      temp$max_cor <- NA
    }
    df.out4.sc <- rbind(df.out4.sc, temp)
  }
}

#df.out4.sc <- df.out4.sc[-which(is.na(df.out4.sc)),]
write.csv(df.out4.sc, here("data/Cleaned/cdr4_sc.csv"), row.names = FALSE)

##################
# Fix cdr2 columns
test <- read.csv(here('Data/Cleaned/cdr2_sc.csv'))
cdr1_plots <- enviro_list[[1]]
cdr2_plots <- enviro_list[[2]]
names(cdr2_plots) <-c('plot_2', 'ring_2')
plots <- cbind(cdr1_plots, cdr2_plots)

test2 <- dplyr::left_join(test, plots, by = c("Plot", "Ring")) %>%
  select(-c('Plot', "Ring")) %>%
  rename(Plot = plot_2, Ring = ring_2)
write.csv(test2, here('Data/Cleaned/cdr2_sc.csv'), row.names = FALSE)
