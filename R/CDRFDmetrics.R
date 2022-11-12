################################
#### FD METRIC CALCULATIONS ####
################################

# CREATED BY: KAITLIN KIMMEL

# load libraries
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
trait_list <- list(cdr.traits1, cdr.traits2, cdr.traits3, cdr.traits3)

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
df.out1 <- data.frame(SR = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA,
                     Plot = NA, Ring = NA, n_trait = NA, traits = NA, mean_cor = NA, 
                     min_cor = NA, max_cor = NA)
# Loop to run through different trait combinations
for(j in 1:length(trait_comb_list)){
  focal_list <- trait_comb_list[[j]]
  for (i in 1:length(focal_list)){
    x = focal_list[[i]]
    a = comm_list[[1]]
    out <- dbFD(x, a, m = 2) # calculating metrics m = 2
    temp <- data.frame(SR = out$nbsp, FRic = out$FRic, FEve = out$FEve, FDiv = out$FDiv,
                       FDis = out$FDis, RaoQ = out$RaoQ)
    temp <- cbind(temp, enviro_list[[1]])
    temp$n_trait = ncol(focal_list[[i]])
    temp$traits = i
    if(ncol(x) == 4){
      temp.cor <- rquery.cormat(x, type="flatten", graph=FALSE, method = 'spearman')
      temp$mean_cor <- mean(temp.cor$r$cor) # USE ABSOLUTE VALUES? 
      temp$min_cor <- min(temp.cor$r$cor)
      temp$max_cor <- max(temp.cor$r$cor)
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


df.out2 <- data.frame(SR = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA,
                     Plot = NA, Ring = NA, n_trait = NA, traits = NA, mean_cor = NA, 
                     min_cor = NA, max_cor = NA)
for(j in 1:length(trait_comb_list2)){
  focal_list <- trait_comb_list2[[j]]
  for (i in 1:length(focal_list)){
    x = focal_list[[i]]
    a = comm_list[[2]]
    out <- dbFD(x, a, m = 2)
    temp <- data.frame(SR = out$nbsp, FRic = out$FRic, FEve = out$FEve, FDiv =out$FDiv,
                       FDis = out$FDis, RaoQ = out$RaoQ)
    temp <- cbind(temp, enviro_list[[1]])
    temp$n_trait = ncol(focal_list[[i]])
    temp$traits = i
    if(ncol(x) == 4){
      temp.cor <- rquery.cormat(x, type="flatten", graph=FALSE, method = 'spearman')
      temp$mean_cor <- mean(temp.cor$r$cor) # USE ABSOLUTE VALUES? 
      temp$min_cor <- min(temp.cor$r$cor)
      temp$max_cor <- max(temp.cor$r$cor)
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
df.out3 <- data.frame(SR = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA,
                      Plot = NA, Ring = NA, n_trait = NA, traits = NA, mean_cor = NA, 
                      min_cor = NA, max_cor = NA)
# Loop to run through different trait combinations
for(j in 1:length(trait_comb_list3)){
  focal_list <- trait_comb_list3[[j]]
  for (i in 1:length(focal_list)){
    x = focal_list[[i]]
    a = comm_list[[3]]
    out <- dbFD(x, a, m = 2) # calculating metrics m = 2
    temp <- data.frame(SR = out$nbsp, FRic = out$FRic, FEve = out$FEve, FDiv =out$FDiv,
                       FDis = out$FDis, RaoQ = out$RaoQ)
    temp <- cbind(temp, enviro_list[[1]])
    temp$n_trait = ncol(focal_list[[i]])
    temp$traits = i
    if(ncol(x) == 4){
      temp.cor <- rquery.cormat(x, type="flatten", graph=FALSE, method = 'spearman')
      temp$mean_cor <- mean(temp.cor$r$cor) # USE ABSOLUTE VALUES? 
      temp$min_cor <- min(temp.cor$r$cor)
      temp$max_cor <- max(temp.cor$r$cor)
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
df.out4 <- data.frame(SR = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA,
                      Plot = NA, Ring = NA, n_trait = NA, traits = NA, mean_cor = NA, 
                      min_cor = NA, max_cor = NA)
# Loop to run through different trait combinations
for(j in 1:length(trait_comb_list4)){
  focal_list <- trait_comb_list4[[j]]
  for (i in 1:length(focal_list)){
    x = focal_list[[i]]
    a = comm_list[[4]]
    out <- dbFD(x, a, m = 2) # calculating metrics m = 2
    temp <- data.frame(SR = out$nbsp, FRic = out$FRic, FEve = out$FEve, FDiv =out$FDiv,
                       FDis = out$FDis, RaoQ = out$RaoQ)
    temp <- cbind(temp, enviro_list[[1]])
    temp$n_trait = ncol(focal_list[[i]])
    temp$traits = i
    if(ncol(x) == 4){
      temp.cor <- rquery.cormat(x, type="flatten", graph=FALSE, method = 'spearman')
      temp$mean_cor <- mean(temp.cor$r$cor) # USE ABSOLUTE VALUES? 
      temp$min_cor <- min(temp.cor$r$cor)
      temp$max_cor <- max(temp.cor$r$cor)
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
#### SAVE OUTPUT FOR ANALYSIS

write.csv(df.out1, here("data/Cleaned/cdr1.csv"), row.names = FALSE)
write.csv(df.out2, here("data/Cleaned/cdr2.csv"), row.names = FALSE)
write.csv(df.out3, here("data/Cleaned/cdr3.csv"), row.names = FALSE)
write.csv(df.out4, here("data/Cleaned/cdr4.csv"), row.names = FALSE)