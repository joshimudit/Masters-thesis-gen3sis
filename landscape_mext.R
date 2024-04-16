# Preparing directories & libraries
#setwd("E:/thesis/thesis_materials/result_paleographies")
#setwd(getwd())
library(raster)
library(gen3sis)
library(rmatio)
library(dplyr)
library(geoscale)
library(readxl)
library(ncdf4)

climate = "E:/thesis/thesis_materials/PLASIM_Phanerozoic_V2/"
topography= "E:/thesis/thesis_materials/Scotese_paleogeographic_reconstructions_2018/Scotese_Wright_2018_Maps_1_88_1degX1deg_PaleoDEMS_nc_v2/"
landscapes_list <- list()

#Defining timesteps and CO2 levels of climate data
climate_steps <- seq(-400, 0, 5)*1e+6
co2levels <- c(50, seq(500, 3000, 500))


######### Reference CO2 curve ----------------------------------------------------
# pCO2 peaks during 4 mass extinctions
co2data <- read_excel("co2_loess_fit.xlsx")
colnames(co2data)
co2data$time <- co2data$age * (-1e+6)
co2data$pCO2_ext = co2data$pCO2_maxprob
#View(co2data)

#late Devonian peak (360 mya) ~ 5000ppm (Cui et.al. 2021)
#co2data$pCO2_ext[504 : 506] = c(2500, 5000, 2500)
co2data$pCO2_ext[721] = c(5000)
#late Permian peak (252 mya) ~ 5000ppm (Cui et.al. 2021)
co2data$pCO2_ext[505] = c(5000)
#late Triassic peak (200 mya) ~ 4500ppm (Beerling et.al. 2002)
co2data$pCO2_ext[401] = c(4500)
#late Cretaceous peak (65 mya) ~ 2000ppm (Keller et.al. 2005)
co2data$pCO2_ext[131] = c(2000)



co2curve <- smooth.spline(co2data$time, co2data$pCO2_ext, spar = NULL)
plot(co2data$time, co2data$pCO2_ext, type="o", col="light blue", pch="o",ylim=c(0, 5000), ylab="pC02 (ppm)", xlab= "Time (mya)", lty=1, cex= 0.2, xlim = c(-400003900,0), main = "Scenario II: pCO2 peaks during 4 mass extinctions")

lines(co2curve, col="dark red", lty=3)




for (t in seq(0, -400e+6, -5e+5)){

  
  ### pCO2 peaks during 4 mass extinctions
  #late Devonian peak (360 mya)
  #late Permian peak (252 mya)
  #late Triassic peak (200 mya)
  #late Cretaceous peak (65 mya)
  # UPTO 3000 ppm only, because climate data available for maximum of 3000 ppm of pCO@
  
  if(t %in% seq(-360e+6, -359e+6, 5e+5)){
    co2_t <- 3000
  } else if (t %in% seq(-252e+6, -251e+6, 5e+5)){
    co2_t <- 3000
  } else if (t %in% seq(-200e+6, -199e+6, 5e+5)){
    co2_t <- 3000
  } else if (t %in% seq(-65e+6, -64e+6, 5e+5)){
    co2_t <- 2000
  } else {
    co2_t <- predict(co2curve, t)$y
  }
  
  
  
  # defining time 't' and the contribution of early and later time step on the value of t
  if(t %in% climate_steps){
    early_step <- t
    late_step <- t
    early_weight <- 1
    later_weight <- 0
  } else {
    late_step <- min(climate_steps[climate_steps > t])
    early_step <- max(climate_steps[climate_steps < t])
    early_weight <- 1 - abs((early_step - t))/(late_step - early_step)
    later_weight <- 1 - early_weight
  }
  
  
  # Extracting the interpolated CO2 values at any time 't'
  if(co2_t %in% co2levels){
    lower_co2 <- co2_t
    upper_co2 <- co2_t
    contribution_co2lower <- 1
    contribution_co2upper <- 0
  } else {
    lower_co2 <- co2levels[which(co2_t - co2levels == min(co2_t - co2levels[co2_t - co2levels > 0]))]
    upper_co2 <- co2levels[which(co2_t - co2levels == max(co2_t - co2levels[co2_t - co2levels <= 0]))]
    contribution_co2lower <- 1 - (co2_t - lower_co2)/(upper_co2 - lower_co2)
    contribution_co2upper <- 1 - (upper_co2 - co2_t)/(upper_co2 - lower_co2)
  }
  
  
  
  ######### Interpolating Topographies ----------------------------------------------------
  topos <- list.files(topography, pattern="\\.nc$")
  
  index <- grep(paste("_",early_step/(-1e+6), "Ma.nc", sep = ""), topos)
  topo_early <- raster(paste(topography,topos[index], sep=""))
  index <- grep(paste("_",late_step/(-1e+6), "Ma.nc", sep = ""), topos)
  topo_late <- raster(paste(topography,topos[index], sep=""))
  
  topo_t <- topo_early*early_weight + topo_late*later_weight
  
  # resampling into target resolution
  r <- raster(nrows=48, ncols=96,
              crs="+proj=longlat +datum=WGS84 +no_defs", 
              resolution = c(3.75, 3.75))
  topo_t <- resample(topo_t, r, method="bilinear")
  
  landscape_t <- topo_t
  landscape_t[topo_t <= 0] <- NA 
  
  
  
  
  ######### Interpolating Temperature ----------------------------------------------------
  
  # Calling temperature of time step previous than 't'
  T_early_lower_co2 <- raster(paste(climate,"run_",early_step/(-1e+6),"Ma/processed_",lower_co2,"_ppm/timmean_",lower_co2,"_ppm.nc", sep=""), var = "ts", level = 1)
  T_early_lower_co2 <- rotate(T_early_lower_co2)
  T_early_higher_co2 <- raster(paste(climate,"run_",early_step/(-1e+6),"Ma/processed_",upper_co2,"_ppm/timmean_",upper_co2,"_ppm.nc", sep=""), var = "ts", level = 1)
  T_early_higher_co2 <- rotate(T_early_higher_co2)
  T_early_co2 <- (contribution_co2lower * T_early_lower_co2 + contribution_co2upper * T_early_higher_co2) - 273.15
  # Calling temperature of time step after than 't'
  T_later_lower_co2 <- raster(paste(climate,"run_",late_step/(-1e+6),"Ma/processed_",lower_co2,"_ppm/timmean_",lower_co2,"_ppm.nc", sep=""), var = "ts", level = 1)
  T_later_lower_co2 <- rotate(T_later_lower_co2)
  T_later_higher_co2 <- raster(paste(climate,"run_",late_step/(-1e+6),"Ma/processed_",upper_co2,"_ppm/timmean_",upper_co2,"_ppm.nc", sep=""), var = "ts", level = 1)
  T_later_higher_co2 <- rotate(T_later_higher_co2)
  T_later_co2 <- (contribution_co2lower * T_later_lower_co2 + contribution_co2upper * T_later_higher_co2) - 273.15
  
  # temperature value at time 't'
  temp_t <- T_early_co2*early_weight + T_later_co2*later_weight
  temp_t[topo_t <= 0] <- NA
  
  
  
  
  ######### Interpolating Aridity ----------------------------------------------------
  # Need 3 rasters: Temperature, Net Shortwave radiation, Precipitation
  
  ### TEMPERATURE (Latent heat)
  temp_t_K <- temp_t + 273.15 # [to Kelvin]
  #lambda is latent heat of evaporation in mJ m-3 and is a function of temperature lambda = (3.146 - 0.002361 Td) * 10^3 in MJ m-3
  latent_heat <- (3.146 - 0.002361*temp_t_K) * 10^3
  
  
  ### Net Shortwave radiation (R): [W m-2] = [J s-1 m-2]
  # Calling Net Shortwave radiation of time step previous than 't'
  R_early_lower_co2 <- raster(paste(climate,"run_",early_step/(-1e+6),"Ma/processed_",lower_co2,"_ppm/timmean_",lower_co2,"_ppm.nc", sep=""), var = "rss", level = 4)
  R_early_lower_co2 <- rotate(R_early_lower_co2)
  R_early_higher_co2 <- raster(paste(climate,"run_",early_step/(-1e+6),"Ma/processed_",upper_co2,"_ppm/timmean_",upper_co2,"_ppm.nc", sep=""), var = "rss", level = 4)
  R_early_higher_co2 <- rotate(R_early_higher_co2)
  R_early_co2 <- (contribution_co2lower * R_early_lower_co2 + contribution_co2upper * R_early_higher_co2)
  # Calling Net Shortwave radiation of time step after than 't'
  R_later_lower_co2 <- raster(paste(climate,"run_",late_step/(-1e+6),"Ma/processed_",lower_co2,"_ppm/timmean_",lower_co2,"_ppm.nc", sep=""), var = "rss", level = 4)
  R_later_lower_co2 <- rotate(R_later_lower_co2)
  R_later_higher_co2 <- raster(paste(climate,"run_",late_step/(-1e+6),"Ma/processed_",upper_co2,"_ppm/timmean_",upper_co2,"_ppm.nc", sep=""), var = "rss", level = 4)
  R_later_higher_co2 <- rotate(R_later_higher_co2)
  R_later_co2 <- (contribution_co2lower * R_later_lower_co2 + contribution_co2upper * R_later_higher_co2)
  
  # Net Shortwave radiation value at time 't'
  rss_t <- R_early_co2*early_weight + R_later_co2*later_weight
  rss_t <- rss_t * (365*24*60*60*1e-6) # from [J s-1 m-2] to [MJ year-1 m-2]
  rss_t[topo_t <= 0] <- NA
  
  
  ### PRECIPITATION (P)
  # Calling precipitation of time step previous than 't'
  P_early_lower_co2 <- raster(paste(climate,"run_",early_step/(-1e+6),"Ma/processed_",lower_co2,"_ppm/timmean_",lower_co2,"_ppm.nc", sep=""), var = "pr", level = 1)
  P_early_lower_co2 <- rotate(P_early_lower_co2)
  P_early_higher_co2 <- raster(paste(climate,"run_",early_step/(-1e+6),"Ma/processed_",upper_co2,"_ppm/timmean_",upper_co2,"_ppm.nc", sep=""), var = "pr", level = 1)
  P_early_higher_co2 <- rotate(P_early_higher_co2)
  P_early_co2 <- (contribution_co2lower * P_early_lower_co2 + contribution_co2upper * P_early_higher_co2)
  # Calling precipitation of time step after than 't'
  P_later_lower_co2 <- raster(paste(climate,"run_",late_step/(-1e+6),"Ma/processed_",lower_co2,"_ppm/timmean_",lower_co2,"_ppm.nc", sep=""), var = "pr", level = 1)
  P_later_lower_co2 <- rotate(P_later_lower_co2)
  P_later_higher_co2 <- raster(paste(climate,"run_",late_step/(-1e+6),"Ma/processed_",upper_co2,"_ppm/timmean_",upper_co2,"_ppm.nc", sep=""), var = "pr", level = 1)
  P_later_higher_co2 <- rotate(P_later_higher_co2)
  P_later_co2 <- (contribution_co2lower * P_later_lower_co2 + contribution_co2upper * P_later_higher_co2)
  
  # precipitation value at time 't'
  prec_t <- P_early_co2*early_weight + P_later_co2*later_weight
  prec_t <- prec_t * 365*24*60*60 # conversion from [m s-1] to [m year-1]
  prec_t[topo_t <= 0] <- NA
  
  
  #
  ### Budyko aridity index: AI = Rn/(Lamda * P)
  #
  #Rn->average daily radiation MJ m-2 day-1
  #P is average daily rainfall in m day-1
  #lambda is latent heat of evaporation in mJ m-3 and is a function of temperature lambda = (3.146 - 0.002361 Td) * 10^3 in MJ m-3
  
  budyko_aridity <- rss_t/(latent_heat*prec_t)
  # declare everyting above 7 as super arid
  budyko_aridity[budyko_aridity >= 7] <- 7
  # Normalize and inverse
  aridity_t <- 1-(budyko_aridity/7) # 0 --> full aridity, 1 --> no aridity
  
  
  
  
  
  ######### Results ----------------------------------------------------
  
  #Compiling all the interpolations inside a landscape list
  landscapes_list$topo <- c(landscapes_list$topo, landscape_t)
  landscapes_list$temp <- c(landscapes_list$temp, temp_t)
  landscapes_list$aridity <- c(landscapes_list$aridity, aridity_t)
  
  par(mfrow = c(1,3))
  plot(landscape_t, zlim=c(0, 5000))
  plot(temp_t, zlim=c(-50, 40), main = paste(t))
  plot(aridity_t,zlim=c(0, 1), main = paste(t))
  
  print(t)
  print(co2_t)
  print(mean(values(temp_t), na.rm=T))
  
  
  # #
  # #Saving plots
  # # Check if the current timestep matches the required ones (360 myrs, 250 myrs, 200 myrs, and 65 myrs)
  # if (t == -360e+6 || t == -252e+6 || t == -200e+6 || t == -65e+6) {
  #   # Generate plots
  #   par(mfrow = c(1, 3))
  #   plot(landscape_t, zlim = c(0, 5000))
  #   plot(temp_t, zlim = c(-50, 40), main = paste(t))
  #   plot(aridity_t, zlim = c(0, 1), main = paste(t))
  #   
  #   # Save the plots as image files (adjust the file paths and names as needed)
  #   png(paste0("mext_landscape_", abs(t / 1e+6), "Ma.png"))
  #   plot(landscape_t, zlim = c(0, 5000))
  #   dev.off()
  #   
  #   png(paste0("mext_temperature_", abs(t / 1e+6), "Ma.png"))
  #   plot(temp_t, zlim = c(-50, 40), main = paste(t))
  #   dev.off()
  #   
  #   png(paste0("mext_aridity_", abs(t / 1e+6), "Ma.png"))
  #   plot(aridity_t, zlim = c(0, 1), main = paste(t))
  #   dev.off()
  # }
}





# Cost Function
cost_function <- function(source, habitable_src, dest, habitable_dest) {
  if (!all(habitable_src, habitable_dest)) {
    return(4/1000) # Costs 4 times more to cross grid cells
  } else {
    return(1/1000)
  }
}


#CREATING FINAL INPUT LANDSCAPE
create_input_landscape(landscapes = landscapes_list,
                       timesteps = as.character(seq(0,length(seq(-400e+6, 0, 5e+5))-1, 1)),
                       cost_function = cost_function,
                       directions = 8,
                       output_directory = "landscape_mext",
                       overwrite = T,
                       crs = "+proj=longlat +datum=WGS84 +no_defs",
                       calculate_full_distance_matrices = T,
                       verbose = T)





dev.off()
#PLOT LINE CUREVES
plot(co2data$time, co2data$pCO2_ext, type="o", col="light blue", pch="o",ylim=c(0, 5000), ylab="pC02 (ppm)", xlab= "Time (mya)", lty=1, cex= 0.2, xlim = c(-400003900,0), main = "Scenario II: pCO2 peaks during 4 mass extinctions")
lines(co2data$time, co2data$pCO2_maxprob, col="dark red", lty=3)
points(co2data$time, predict(co2curve, co2data$time)$y)
lines(co2curve)
#text(locator(), labels = c("Late Devonian (360 ma)", "End Permian (252 ma)", "End Triassic (200 ma)", "End Cretaceous (65 ma)"), cex= 0.8, col = "dark blue")


legend(-116503900,4500,legend=c("pCO2_mext", "pCO2_maxprob"), col=c("light blue","dark red"),
       pch=c("",""),lty=c(1,3), ncol=1, cex= 0.9)


