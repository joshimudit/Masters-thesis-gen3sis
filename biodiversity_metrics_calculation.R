# /cluster/scratch/mjoshi


# load libraries ##############
library(ape)
library(raster)
library(picante)
library(apTreeshape)
library(moments)
library(adephylo)
library(dplyr)
library(phylobase)
# library(pez) # comparative.com
# library(phytools)
# library(FD)
# library(SYNCSA)


#########
model <- "maxprob"

# summary statistics table from gen3sis output
params <- readRDS(paste0("/cluster/scratch/mjoshi/output_maxprob/sims_summary_maxprob.rds"))
setwd("/cluster/scratch/mjoshi/output_maxprob/")

##############


# Columns for the summary statistics ##############

# # time step
# params$running_time_step <- NA
# params$n_total_diversity <- NA
# params$n_extinct_diversity <- NA
# params$n_extant_diversity <- NA

# Specie richness
params$SpecRich_total  <- NA
params$SpecRich_mean   <- NA
params$SpecRich_median <- NA
params$SpecRich_sd <- NA
params$SpecRich_max <- NA
params$SpecRich_skewness <- NA
params$SpecRich_lat_cor <- NA

# Clade age
params$clade_age <- NA

# Specie range
params$range_size_mean <- NA
params$range_size_max <-NA
params$range_size_median <- NA
params$range_size_sd <- NA
params$range_size_skewness <- NA

# Weighted endemism
params$endemism_weighted_mean <- NA
params$endemism_weighted_max <-NA
params$endemism_weighted_median <- NA
params$endemism_weighted_sd <- NA
params$endemism_weighted_skewness <- NA
params$endemism_weighted_lat <- NA


# Phylogenetic diversity (Faith's)
params$pd_mean <- NA
params$pd_max <-NA
params$pd_median <- NA
params$pd_sd <- NA
params$pd_skewness <- NA
params$pd_lat_cor <- NA


#Relised temp niche
params$temp_mean <-  NA
params$temp_max <-  NA
params$temp_median  <-  NA
params$temp_sd <- NA
params$temp_skewness <- NA

# Temperature trait
params$Trait_mean <- NA
params$Trait_max <-NA
params$Trait_median <- NA
params$Trait_sd <- NA
params$Trait_skewness <- NA





#####################

## Extra functions############################


abundanceCell <- function(species){
  abundance <- c()
  for(sp in 1:length(species)){
    x        <- species[[sp]]$abundance
    abundance <- c(abundance,x)
  }
  
  # sum abundance per cell
  abundance_cell <- tapply(abundance,names(abundance),sum)
  
  return(abundance_cell)
}

rescaleTree<-function(tree,scale){
  tree$edge.length<-
    tree$edge.length/max(nodeHeights(tree)[,2])*scale
  return(tree)
}

scale.to <- function(vec,vec.sum) {
  #mat is a vector
  #this function rescales each the vector values to sum to 'vec.sum'
  vec.tot <- sum(vec,na.rm=TRUE)
  if (vec.tot > 0) {
    vec.out <- vec.sum*vec/vec.tot
  } else {
    vec.out <- rep(0,times=length(vec))  #should columns for tips / branches with no occurrences be removed?
  }
  return(vec.out)
}

total_range <- function(grid){
  geographic_area <- as.data.frame(colnames(grid))
  colnames(geographic_area) <- c("Species")
  geographic_area$range_size <- NA
  geographic_area$range_size <- colSums(grid)
  return(geographic_area)
}

weighted_endemism <- function(x){
  proportions_matrix <- apply(x, 2, FUN=function(y){y <- y/sum(y)})
  WE <- rowSums(proportions_matrix, na.rm = TRUE)
  return(WE)
}

calc_PE <- function(tree, sites_x_tips,presence=c("presence","abundance","probability")) {
  
  # add code to check that the values are correct for the presence type:
  # 0 or 1 for presence - this calculates PE (Rosauer et al 2009)
  # from 0 to 1 for probability - this calculates model weighted PE (Rosauer, in prep)
  # any value for abundance - this calculation is equivalent to BED (Cadotte & Davies 2010)
  
  #default value for presence
  if (is.na(presence)) {presence="presence"}
  
  # change to a phylobase phylo4 object
  if (class(tree) == "phylo") {tree <- phylo4(tree)}
  
  sites_x_branches <- data.frame(cbind(rep(0,nrow(sites_x_tips))))
  
  for (i in 1:nTips(tree)) {
    sites_x_branches[,i] <- sites_x_tips[,which(labels(tree)[i]==colnames(sites_x_tips))]
    names( sites_x_branches)[i] <- labels(tree)[i]
  }
  
  rm(sites_x_tips); gc()
  branch.labels <- as.character(labels(tree))
  branch.count <- length(labels(tree))
  
  # add names and occupancy columns for internal branches
  for (i in (nTips(tree)+1):branch.count) {
    branch.labels[i] <- paste("b",i,sep="")
    desc <- as.integer(descendants(tree,i, type="tips"))
    if (presence=="abundance") {
      branch_col <- as.numeric(apply(sites_x_branches[,desc],MARGIN=1,FUN=sum))
    } else if (presence=="presence") {
      branch_col <- as.numeric(apply(sites_x_branches[,desc],MARGIN=1,FUN=max))
    } else if (presence=="probability") {
      branch_col <- as.numeric(apply(sites_x_branches[,desc],MARGIN=1,FUN=parent.prob))
    }
    sites_x_branches[,i] <- branch_col
    names(sites_x_branches[i]) <- branch.labels[i]
    #cat(i,branch.labels[i],length(desc),"\n")
    gc(verbose=F)
  }
  
  #scale columns (branches) to sum to 1
  sites_x_branches <- apply(sites_x_branches,MARGIN=2,FUN=scale.to,1)
  
  #now scale branches to sum to their length
  branch.lengths <- as.numeric(edgeLength(tree,1:branch.count))
  branch.lengths[is.na(branch.lengths)] <- 0
  for (i in 1:length(branch.lengths)) {
    sites_x_branches[,i] <- sites_x_branches[,i] * branch.lengths[i]
  }
  
  PE.vec <- apply(sites_x_branches,MARGIN=1,FUN=sum,na.rm=T)
  PE <- data.frame(cbind(1:nrow(sites_x_branches),PE.vec))
  names(PE) <- c("site","PE")
  return(PE)
}

tree_metrics <- function(tree){
  ED <- evol.distinct(tree, type = c("fair.proportion"),
                      scale = TRUE, use.branch.lengths = TRUE)
  names(ED)[2] <- "ED"
  ES <- evol.distinct(tree, type = c("equal.splits"),
                      scale = TRUE, use.branch.lengths = TRUE)
  
  names(ES)[2] <- "ES"
  tree_metrics <- left_join(ED,ES, by= "Species")#, match = "all")
  tree_metrics$DivRate <- 1/tree_metrics$ES
  return(tree_metrics)
}



# Calculation ################

landscape <- readRDS("/cluster/scratch/mjoshi/landscape_maxprob/landscapes.rds")
landscape_t <- landscape$temp
landscape_t$site <- rownames(landscape_t)

setwd(paste0("/cluster/scratch/mjoshi/output_maxprob")) #output for the respective landscape folders
file_list <- paste('thesis_m1_config_',1:100)


for(i in c(1:100)){
  #browser()  
  try(setwd(paste0("/cluster/scratch/mjoshi/output_maxprob/thesis_m1_config_",i)))
  
  try(setwd(paste0("/cluster/scratch/mjoshi/output_maxprob/thesis_m1_config_",i,"/richness")))
  print(file_list[i])
  try(files <- list.files())
  
  try(most_recent <- min(as.numeric(sapply(files, FUN=function(x){strsplit(strsplit(x, "_")[[1]][3], ".rds")[[1]]}))))
  if(most_recent != 0){ most_recent<- most_recent+1}
  try(params$running_time_step[i] <- most_recent)
  
  
  # RICHNESS SUMMARY STATS ###########
  richness <- try(readRDS(paste0("richness_t_", most_recent,".rds")))
  
  try(params$SpecRich_mean[i] <- mean(richness, na.rm=T))
  try(params$SpecRich_max[i] <- max(richness, na.rm=T))
  try(params$SpecRich_median[i] <- median(richness, na.rm=T))
  try(params$SpecRich_sd[i] <- sd(richness, na.rm=T))
  try(params$SpecRich_skewness[i] <- skewness(richness, na.rm=T))
  try(params$SpecRich_total[i] <- sum(richness, na.rm=T))
  
  try(richness <- data.frame(richness))
  try(richness$site <- rownames(richness))
  
  try(coordinates <- landscape_t %>% 
        dplyr::select("x", "y", "site"))
  
  try(spatial_metrics <- left_join(richness, coordinates, by='site'))#, match='all'))
  
  # Next simulation (in case of extinction)
  if(sum(richness$richness)==0){
    params$n_extant_diversity[i] <- 0
    setwd("../../")
    next
  }
  
  try(params$SpecRich_lat_cor[i] <-  cor(spatial_metrics$richness, abs(spatial_metrics$y), use='complete.obs'))
  
  
  
  
  # SPECIES-LEVEL SUMMARY STATS #############
  try(setwd("../species"))
  try(species <- readRDS(paste0("species_t_", most_recent,".rds")))
  try(params$n_extant_diversity[i] <- sum(sapply(species, FUN=function(x){if(length(x$abundance)>0){return(1)}else{0}})))
  
  try(setwd("../landscapes"))
  if(most_recent > 0){
    setwd("../../")
    next
  }
  
  
  # presence absence matrix
  try(landscape_t_0 <- readRDS("landscape_t_0.rds"))
  try(all_cells <- rownames(landscape_t_0$coordinates))
  try(all_species_presence <- do.call( cbind, lapply(species, FUN = function(x) {ifelse(all_cells %in% names(x$abundance), 1, 0)})))
  try(colnames(all_species_presence ) <- unlist(lapply(species, function(x){x$id}))) # colnames= species names
  try(presence_absence_matrix <- cbind(landscape$coordinates, all_species_presence))
  
  try(colnames(presence_absence_matrix) <- paste0("species", colnames(presence_absence_matrix)))
  
  ##################
  
  
  
  
  #PHYLOGENY SUMMARY STATS ################
  try(setwd("../phylogeny"))
  phy <- try(read.nexus(paste0("phylogeny_t_", most_recent,".nex")))
  if(class(phy)=="try-error"){setwd("../../") ; next}  
  
  # phy's will differ based on age of clade - so try to make comparable 
  try(phy_og <- phy)
  
  # CLADE AGE
  try(params$clade_age[i] <- max(node.age(phy)$ages))
  
  #phy <- try(drop.fossil(phy))
  # tree manually dropping extincts
  if(class(phy) == "try-error"){
    extinct_species <- paste0("species",which(sapply(species, FUN=function(x){length(x$abundance)==0})))
    phy <- phy_og
    phy <- try(drop.tip(phy, extinct_species))
  }
  # if that fails too move on
  if(class(phy) == "try-error"){setwd("../../") ; next}
  
  try(phy <- rescaleTree(phy, 1))
  
  
  # can't do uch if no species so move on to next sim
  if(length(phy$tip) < 3){setwd("../../") ; next} else { 
    
    if(any(phy$edge.length < 0 )){ phy$edge.length <- phy$edge.length + abs(min(phy$edge.length))}}
  ##################
  
  

  
  try(phy <- drop.tip(phy, setdiff(colnames(presence_absence_matrix), phy$tip.label)))
  try(presence_absence_matrix <- presence_absence_matrix[,colnames(presence_absence_matrix) %in% phy$tip.label])
  
  
  # Species range #######
  try(grid_range <- total_range(presence_absence_matrix))
  try(params$range_size_mean[i] <- mean(grid_range$range_size, na.rm=T))
  try(params$range_size_max[i] <- max(grid_range$range_size, na.rm=T))
  try(params$range_size_median[i] <- median(grid_range$range_size, na.rm=T))
  try(params$range_size_sd[i] <-  sd(grid_range$range_size, na.rm=T))
  try(params$range_size_skewness[i] <- skewness(grid_range$range_size, na.rm=T))
  
  # Weighted endemism #######
  try(spatial_metrics$endemism <- weighted_endemism(presence_absence_matrix))
  try(params$endemism_weighted_mean[i] <- mean(spatial_metrics$endemism, na.rm=T)) 
  try(params$endemism_weighted_max[i] <- max(spatial_metrics$endemism, na.rm=T))
  try(params$endemism_weighted_median[i] <- median(spatial_metrics$endemism, na.rm=T))
  try(params$endemism_weighted_sd[i] <- sd(spatial_metrics$endemism, na.rm=T))
  try(params$endemism_weighted_skewness[i] <- skewness(spatial_metrics$endemism, na.rm=T))
  try(params$endemism_weighted_lat[i] <- cor(spatial_metrics$endemism, abs(spatial_metrics$y), use = "complete.obs"))
  
  ##################
  
  
  

  # Faith's Phylogeny diversity (using picante)
  try(pd_estimate <- pd(presence_absence_matrix, phy, include.root=TRUE))
  
  try(params$pd_mean[i] <- mean(pd_estimate$PD, na.rm=T))
  try(params$pd_max[i] <- max(pd_estimate$PD, na.rm=T))
  try(params$pd_median[i] <- median(pd_estimate$PD, na.rm=T))
  try(params$pd_sd[i] <- sd(pd_estimate$PD, na.rm=T))
  try(params$pd_skewness[i] <- skewness(as.data.frame(pd_estimate$PD, na.rm=T)))
  
  try(spatial_metrics$PD <- pd_estimate)
 
  
  
  
  
  #
  # Trait-based metrics
  #
  
  try(setwd("../traits"))
  try(traits <- readRDS(paste0("traits_t_", most_recent,".rds")))
  
  # realized temperature niche
  realised_temp_niche <-  do.call(rbind, lapply(species, FUN=function(x){
    realised_temp <- landscape_t[which(landscape_t$site %in% names(x$abundance)),"0"]
    return(data.frame(Species=paste0("species", x$id), niche_trait=mean(realised_temp, na.rm=T)))
  }))
  
  try(params$temp_mean[i]   <-  mean(realised_temp_niche$niche_trait, na.rm=T)) # mean temp niche
  try(params$temp_max[i]   <-  max(realised_temp_niche$niche_trait, na.rm=T)) # max temp niche
  try(params$temp_median[i]   <-  median(realised_temp_niche$niche_trait, na.rm=T)) # median temp niche
  try(params$temp_sd[i]     <- sd(realised_temp_niche$niche_trait, na.rm=T))# sd temp niche
  try(params$temp_skewness[i] <- skewness(realised_temp_niche$niche_trait, na.rm=T))
  
  
  
  # Temperature trait
  try(temp_trait_means <- unlist(lapply(traits, FUN=function(x){y=mean(x[,"temp"], na.rm=T);return(y)})))
  # Temperature trait value
  try(params$Trait_mean[i] <- mean(temp_trait_means, na.rm=T))
  try(params$Trait_max[i] <- max(temp_trait_means, na.rm=T))
  try(params$Trait_median[i] <- median(temp_trait_means, na.rm=T))
  try(params$Trait_sd[i] <- sd(temp_trait_means, na.rm=T))
  try(params$Trait_skewness[i] <- skewness(temp_trait_means, na.rm=T))
  
 
  
  
  ####################### 
  
  # save everything
  try(setwd("../"))
  try(results_list <- list("spatial_metrics"=spatial_metrics,
                           "traits"=traits,
                           "phylogeny"=phy,
                           "species" = species))
  
  try(saveRDS(results_list, file="results_list_maxprob.rds"))
  try(setwd("../"))
  
  try(saveRDS(results_list, file="results_list_maxprob.rds"))
  try(saveRDS(params, file = paste0("summary_stats_maxprob.rds")))
  try(saveRDS(params, file = paste0("summary_stats_maxprob.txt")))
}


