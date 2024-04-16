######################################
###            METADATA            ###
######################################
# Mudit Joshi
#
# Date: 11.07.2023
#
# Landscape: Scotese 2018, 400mya, 800 time steps (500kyr)
#
# Project: Masters thesis
######################################

# setwd("E:/thesis/gen3sis_run")

library(gen3sis)
library(randtoolbox)

########################
### General settings ###
########################


random_seed = params$seed
start_time = params$start_time
end_time = 0
max_number_of_species = params$max_species
max_number_of_coexisting_species = 1e6
initial_abundance = params$initial_abundance

# a list of traits to include with each species
trait_names = c("temp","aridity")

# ranges to scale the input environemts with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# lsited with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]
# Scotese 220km range: -37.50240, 26.67235

# Range of environmental teperatures between -20 and +40 degree Celsius
environmental_ranges = list("temp" = c(-20, 40), "aridity" = c(0, 1))


#########################
### Observer Function ###
#########################


end_of_timestep_observer = function(data, vars, config){
  save_richness()
  save_species()
  save_traits()
  save_phylogeny()
  save_divergence()

  
  plot_richness(data$all_species, data$landscape)

  
  # make p/a matrices
  if(!file.exists(file.path(config$directories$output, "pa_matrices"))){dir.create(file.path(config$directories$output, "pa_matrices"))}
  pa_matrix <- data.frame(matrix(0,nrow=length(data$landscape$coordinates[,1]), ncol=(length(data$all_species)+2)))
  pa_matrix[,1:2]<-data$landscape$coordinates
  rownames(pa_matrix) <- rownames(data$landscape$coordinates)
  names(pa_matrix)[1:2] <- c("x", "y")
  names(pa_matrix)[3:length(pa_matrix)] <-unlist(lapply(data$all_species, FUN=function(x){x$id}))
  for(i in 3:(length(pa_matrix[1,]))){
    pa_matrix[names(data$all_species[[i-2]]$abundance),i] <- 1
  }
  write.table(pa_matrix, file.path(config$directories$output,"pa_matrices",  paste0("species_pa_",vars$ti, ".txt")), row.names = F, col.names = T)
}



######################
### Initialization ###
######################

create_ancestor_species <- function(landscape, config) {
  
  # raster of global temperatures 
  # browser()
  t_world <- raster::rasterFromXYZ(cbind(landscape$coordinates,landscape$environment[, 2:3, drop = F]))
  t_world <- extend(t_world, landscape$extent)
  
  # Will the simulation sample the ancestor from a partcular continent or randomly from anywhere?
  start_cells <- as.character(Which(t_world$temp > 0 & t_world$temp < 1 & t_world$aridity > 0.1 & t_world$aridity < 1, cells=T))
  
  # the list of new species
  all_species <- list()
  
  new_species <- create_species(as.character(start_cells), config)
  new_species$traits[ , "dispersal"] <- 1
  new_species$traits[ , "temp"] <- landscape$environment[start_cells, "temp"]
  new_species$traits[ , "aridity"] <- landscape$environment[start_cells, "aridity"]
  
  all_species <- append(all_species, list(new_species))
  
  
  return(all_species)
  
}


#################
### Dispersal ###
#################



get_dispersal_values <- function(n, species, landscape, config) {
  values <- rweibull(n, shape=2, scale= 500 )
  return(values)
}


##################
### Divergence ###
##################


# threshold for genetic distance after which a speciation event takes place
divergence_threshold = params$divergence_threshold #range from 2 to 10
lambda = params$lambda # range from 2 to 5

get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  
  
  if(any(landscape$environment[,"temp"]>1)){landscape$environment[which(landscape$environment[,"temp"]>1), "temp"] <- 1}
  if(any(landscape$environment[,"temp"]<0)){landscape$environment[which(landscape$environment[,"temp"]<0), "temp"] <- 0}
  
  divergences <- vector("numeric", length(unique(cluster_indices)))
  
  divergence_matrix <- matrix(0, ncol=length(unique(cluster_indices)),  nrow=length(unique(cluster_indices)))
  cluster_indices_order <- unique(cluster_indices)[order(unique(cluster_indices))]
  # order them, find pairwise temp dist
  for(i_x in 1:length(cluster_indices_order)){
    mean_temp_i <- mean(landscape$environment[names(species$abundance[which(cluster_indices %in% cluster_indices_order[i_x])]),"temp"], na.rm=T)
    for(j_x in 1:length(cluster_indices_order)){
      mean_temp_j <- mean(landscape$environment[names(species$abundance[which(cluster_indices %in% cluster_indices_order[j_x])]),"temp"], na.rm=T)
      divergence_matrix[i_x, j_x] <- (sum(mean_temp_i, mean_temp_j)/(2))^lambda
    }
  }
  rownames(divergence_matrix) <- cluster_indices_order
  colnames(divergence_matrix) <- cluster_indices_order
  
  if(length(divergence_matrix) == 1){ divergence_matrix <- as.numeric(divergence_matrix)}else{
    
    diag(divergence_matrix) <- 0
  }
  
  return(divergence_matrix)
}



#######################
### Trait Evolution ###
#######################


sigma_t <- params$sigma_t #range from 0.001 to 0.015
psi <- params$psi

apply_evolution <- function(species, cluster_indices, landscape, config) {
  
  # cell names
  traits <- species[["traits"]]
  cells <- rownames( traits )
  
  # evolve traits for each cluster
  for(cluster_index in unique(cluster_indices)){
    cells_cluster <- cells[which(cluster_indices == cluster_index)]
    t_theta_cluster <- mean(landscape$environment[cells_cluster,"temp"], na.rm=T)
    # evolve temperature
    traits[cells_cluster, "temp"] <- traits[cells_cluster, "temp"] + ( psi * (t_theta_cluster - traits[cells_cluster, "temp"]) ) + rnorm(1, mean = 0, sd = sigma_t)
  }
  
  # set bounds between 0 and 1 so the species can;t evolve a niche beyond that present in the data (all temp is scaled between 0 and 1)
  if(any(traits[, "temp"] > 1)){traits[which(traits[,"temp"]> 1), "temp"] <- 1}
  if(any(traits[, "temp"] < 0)){traits[which(traits[,"temp"]< 0), "temp"] <- 0}
  
  if(any(traits[, "aridity"] > 1)){traits[which(traits[,"aridity"]> 1), "aridity"] <- 1}
  if(any(traits[, "aridity"] < 0.1)){traits[which(traits[,"aridity"] <- 0.1), "aridity"] <- 0.1}
  
  
  return(traits)
  
}


#################################################
### Environmental and Ecological Interactions ###
#################################################


K_opt_max <-  params$K    # environmental filtering parameter (maximum abundance of all individuals in a mesic cell)
omega <- params$omega #range from=0.01 to 0.035     # environmental filtering parameter (higher values equal less restrictive niche)
aridity_cost <- params$aridity_cost # environmental filtering parameter (higher values equal less abundance in arid cells)

x0 <- params$inflexion # extinction parameter: inflexion point of extinction probability curve
decay <- params$decay # extinction parameter: inflexion point of extinction probability curve

apply_ecology <- function(abundance, traits, landscape, config) {
  
  # Realised carrying capacity (reduced in arid cells)
  K <- K_opt_max * exp(1-aridity_cost*(1 - landscape[, "aridity"]))
  
  # difference between landscape and species optimum temperature niche
  diff_t <- abs(traits[, "temp"] - landscape[, "temp"])
  
  # potential population size based on resource use efficiency 
  # omega is strength of environmental filtering
  # will equal K when tthe species is perfectly adapted
  # modified from McPeek 2008/2007
  Nij <- K * exp(-(diff_t/omega)^2)
  
  # carrying capcaity as a zero-sum game - species abundance is scaled based on their resourse use efficiency (Nij)
  Nij_sum <- sum(Nij)
  Nij_hat <- Nij * (sort(c(Nij_sum, K),partial=1)[1] / Nij_sum)
  Nij_hat[which(is.na(Nij_hat))] <- 0
  # now do an extinction filter based on population size
  prob_extinction <- (1/(1+exp(-decay*(x0 -  Nij_hat))))
  Nij_hat[which(sapply(prob_extinction, FUN=function(x){sample(c(0,1),1, prob= c(x, 1-(x)))})==0)] <- 0 

  
  return(Nij_hat)
}







