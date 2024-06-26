######################################
###            METADATA            ###
######################################
#
# Mudit Joshi
#
# Date: 04.10.2023
#
# Project: Masters thesis
#
# Description: Function to generate 100 configs from the template configuration file
# (each representing a unique plant clade)
######################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Dirs and libraries

# set home directory to save configs
HOME <- "E:/thesis/gen3sis_run"
setwd(HOME)

# load packages
library(randtoolbox) # for sobol sequences

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS

# define function to linearly map sobol sequences to desired parameter range
linMap <- function(x, from, to) {
  (x - min(x)) / max(x - min(x)) * (to - from) + from
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Generate params

# number of configs to generate per model
n <- 100
# number of variable parameters to generate sobol sequences for
p <- 5

###### Generate parameter tables and config files ####### 
for(m in c("m")){
  
  # set up and name data frame filled with n sobol sequences for p parameters
  set.seed(666)
  params_table <- data.frame(sobol(n, p, init = T))
  colnames(params_table) <- c("divergence_threshold", "lambda", "omega", "sigma_t", "psi")
  
  # 5 custom parameters to vary
  params_table[,"divergence_threshold"] <- round(linMap(params_table[,"divergence_threshold"], from=2, to=10)) 
  params_table[,"lambda"]               <- round(linMap(params_table[,"lambda"], from=2, to=5),3) 
  params_table[,"omega"]                <- round(linMap(params_table[,"omega"], from=0.01, to=0.035),3) # 0.035 is a niche breadth of appoximately 6-7 degrees (median MAT for vertebrates)
  params_table[,"sigma_t"]              <- round(linMap(params_table[,"sigma_t"], from=0.001, to=0.015),3) 
  params_table[,"psi"]                  <- round(linMap(params_table[,"psi"], from=0.001, to=0.03),3)
  
  # add fixed parameters 
  params_table$seed <- 666
  params_table$start_time <- 800 #400mya and 500kyr timesteps
  params_table$model <- sapply(seq(1, n, 1), FUN=function(x){paste(m,"", x, sep="_")})
  params_table$initial_abundance <- 10000
  params_table$grid_cell_distance <- 220 
  params_table$max_species <- 7500
  params_table$ancestor_number <- 1
  params_table$inflexion <- 1000
  params_table$decay <- 0.02
  params_table$aridity_cost <- 2
  params_table$K <- 30000

  # Write table
  setwd(file.path(HOME))
  write.table(params_table, file=file.path(HOME, paste0(m, "_parameters.txt")), row.names = F, col.names = T)
  
  #setwd(file.path(HOME, "scripts", "gen3sis"))
  if(!dir.exists("config")){dir.create("config"); setwd("config");dir.create("m")}
  
  for(i in 1:nrow(params_table)){ 
    
    #setwd(file.path(HOME, "scripts", "gen3sis"))
    
    params <- params_table[i,]
    
    config_i <- readLines("config.R")
    config_i <- gsub('*.params\\$divergence_threshold', params$divergence_threshold, config_i)
    config_i <- gsub('*.params\\$lambda', params$lambda, config_i)
    config_i <- gsub('*.params\\$sigma_t', params$sigma_t, config_i)
    config_i <- gsub('*.params\\$dispersal', params$dispersal, config_i)
    config_i <- gsub('*.params\\$inflexion', params$inflexion, config_i)
    config_i <- gsub('*.params\\$omega', params$omega, config_i)
    config_i <- gsub('*.params\\$seed', params$seed, config_i)
    config_i <- gsub('*.params\\$start_time', params$start_time, config_i)
    config_i <- gsub('*.params\\$model', params$model, config_i)
    config_i <- gsub('*.params\\$initial_abundance', params$initial_abundance, config_i)
    config_i <- gsub('*.params\\$grid_cell_distance', params$grid_cell_distance, config_i)
    config_i <- gsub('*.params\\$divergence', params$divergence, config_i)
    config_i <- gsub('*.params\\$max_species', params$max_species, config_i)
    config_i <- gsub('*.params\\$ancestor_number', params$ancestor_number, config_i)
    config_i <- gsub('*.params\\$decay', params$decay, config_i)
    config_i <- gsub('*.params\\$K', params$K, config_i)
    config_i <- gsub('*.params\\$aridity_cost', params$aridity_cost, config_i)
    config_i <- gsub('*.params\\$psi', params$psi, config_i)
    
    #setwd(paste0("config/", m))
    writeLines(config_i, paste0('thesis_', m,'_config_', i, '.R'))
    
  }
}



