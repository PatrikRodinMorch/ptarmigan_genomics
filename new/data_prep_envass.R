library(readr)
library(stringr)
library(dplyr)
library(rgdal)
library(raster)
library(maps)
library(mapdata)
library(rstatix)

# Reading in the coordinate data and clean up the species name column
ptarm_coord <- read_csv("LatLong.csv") %>%
  mutate(Spp = str_replace(Spp, "Lagopus ", "L.")) %>%
  mutate(Spp = str_replace(Spp, "L.lagopus scotica", "L.scotica"))
head(ptarm_coord)

# Create species specific df's
muta <- ptarm_coord %>%
  filter(Spp == "L.muta")

lago <- ptarm_coord %>%
  filter(Spp == "L.lagopus" | Spp == "L.scotica")

# Downloading current climate and elevation data at 5-arcmin resolution (10km x 10km)
current_clim <- getData("worldclim", var="bio", res = 5, lon = ptarm_coord[["Long"]], lat = ptarm_coord[["Lat"]])
elevation <- getData("worldclim", var = "alt", res = 5, lon = ptarm_coord[["Long"]], lat = ptarm_coord[["Lat"]])
future_clim_RCP6 <- getData("CMIP5", var="bio", res = 5, lon = ptarm_coord[["Long"]], lat = ptarm_coord[["Lat"]], model = "HE", rcp = 60, year = 70)
future_clim_RCP85 <- getData("CMIP5", var="bio", res = 5, lon = ptarm_coord[["Long"]], lat = ptarm_coord[["Lat"]], model = "HE", rcp = 85, year = 70)

# Number of layers
nlayers(current_clim)
nlayers(elevation)
nlayers(future_clim_RCP6)
nlayers(future_clim_RCP85)

# Focus on the areas around the sample coordinates for muta
sample_range_muta <- extent(min(muta[["Long"]]) - 5, max(muta[["Long"]]) + 5, min(muta[["Lat"]]) - 5, max(muta[["Lat"]]) + 5)
sample_currentenv_muta <- crop(current_clim, sample_range_muta)
sample_futureenv6_muta <- crop(future_clim_RCP6, sample_range_muta)
sample_futureenv85_muta <- crop(future_clim_RCP6, sample_range_muta)
sample_currentelev_muta <- crop(elevation, sample_range_muta)

# Focus on the areas around the sample coordinates for lago
sample_range_lago <- extent(min(lago[["Long"]]) - 5, max(lago[["Long"]]) + 5, min(lago[["Lat"]]) - 5, max(lago[["Lat"]]) + 5)
sample_currentenv_lago <- crop(current_clim, sample_range_lago)
sample_futureenv6_lago <- crop(future_clim_RCP6, sample_range_lago)
sample_futureenv85_lago <- crop(future_clim_RCP6, sample_range_lago)
sample_currentelev_lago <- crop(elevation, sample_range_lago)


# Test plot using bioclim1 and elevation
plot(sample_currentenv_muta[["bio1"]]/10, main="Annual Mean Temperature")
map('worldHires',xlim=c(min(muta[["Long"]])-5,max(muta[["Long"]])+5), ylim=c(min(muta[["Lat"]])-5,max(muta[["Lat"]])+5), fill=FALSE, add=TRUE)
points(muta[["Long"]], muta[["Lat"]], pch="+", cex=1.4, col = "red")

plot(sample_currentelev_muta[["alt"]]/10, main="Elevation")
map('worldHires',xlim=c(min(muta[["Long"]])-5,max(muta[["Long"]])+5), ylim=c(min(muta[["Lat"]])-5,max(muta[["Lat"]])+5), fill=FALSE, add=TRUE)
points(muta[["Long"]], muta[["Lat"]], pch="+", cex=1.4, col = "red")

plot(sample_currentenv_lago[["bio1"]]/10, main="Annual Mean Temperature")
map('worldHires',xlim=c(min(lago[["Long"]])-5,max(lago[["Long"]])+5), ylim=c(min(lago[["Lat"]])-5,max(lago[["Lat"]])+5), fill=FALSE, add=TRUE)
points(lago[["Long"]], lago[["Lat"]], pch="+", cex=1.4, col = "red")

plot(sample_currentelev_lago[["alt"]]/10, main="Elevation")
map('worldHires',xlim=c(min(lago[["Long"]])-5,max(lago[["Long"]])+5), ylim=c(min(lago[["Lat"]])-5,max(lago[["Lat"]])+5), fill=FALSE, add=TRUE)
points(lago[["Long"]], lago[["Lat"]], pch="+", cex=1.4, col = "red")


# Extract the climate variables from the sample points

extract_sp_env <- function(sp, env, elev) {
  coords <- data.frame(Long = sp[["Long"]], Lat = sp[["Lat"]])
  points <- SpatialPoints(coords, proj4string = env@crs)
  
  environment <- extract(env, points)
  elevation <- extract(elev, points)
  
  ind_id <- sp[["Sample_name_vcf"]]
  ptarmigan_environment <- cbind.data.frame(ind_id, coords, elevation, environment)
  
  return(ptarmigan_environment)
}

muta_current_df <- extract_sp_env(muta, sample_currentenv_muta, sample_currentelev_muta)
muta_future6_df <- extract_sp_env(muta, sample_futureenv6_muta, sample_currentelev_muta)
muta_future85_df <- extract_sp_env(muta, sample_futureenv85_muta, sample_currentelev_muta)

lago_current_df <- extract_sp_env(lago, sample_currentenv_lago, sample_currentelev_lago)
lago_future6_df <- extract_sp_env(lago, sample_futureenv6_lago, sample_currentelev_lago)
lago_future85_df <- extract_sp_env(lago, sample_futureenv85_lago, sample_currentelev_lago)

# Here we estimate correlations between current climate predictors and remove them from the df in order to reduce multicollinearity in the modeling step
# Here i will remove one of two correlated predictors if r > |0.7|

correlation_test <- function(df, ...) {
  predictors <- df %>% 
    dplyr::select(!c(ind_id, Long, Lat, ...)) %>%
    cor_mat()
    
    #correlation_matrix <- round(cor(predictors), 2)
    
    return(predictors)
}

# Filter env variables for sp 1. Since these birds are adapted to climate the temperature variables for the warmest and coldest month is more interesting than the mean annual temperature
(muta_current_correlation <- correlation_test(muta_current_df))
(lago_current_correlation <- correlation_test(lago_current_df))

# This results in 6 variables (elevation, max_temp_warmest_month(bio5), min_temp_coldest_month(bio6), temperature_annual_range(bio7), mean_temp_warmest_quarter(bio8), precipitation_seasonality(bio15))
predictors_remove_muta <- c("bio1", "bio2", "bio3", "bio4", "bio9",  "bio10", "bio11", "bio12", "bio13", "bio14", "bio16", "bio17", "bio18", "bio19")
(muta_current_correlation <- correlation_test(muta_current_df, predictors_remove_muta))

# 4 uncorrelated predictors (mean_diurnal_range(bio2), min_temp_coldest_monsth(bio6), precip_wettest_month(bio13), precip_warmest_quartet(bio18)
predictors_remove_lago <- c("elevation", "bio1", "bio3", "bio8", "bio4", "bio7", "bio5", "bio9", "bio10", "bio11","bio12", "bio14", "bio15","bio16", "bio17", "bio19")
(lago_current_correlation <- correlation_test(lago_current_df, predictors_remove_lago))

# Create final csv
muta_final <- muta_current_df %>% 
  dplyr::select(!predictors_remove_muta)

write_csv(muta_final, "muta_env.csv")

lago_final <- lago_current_df %>% 
  dplyr::select(!predictors_remove_lago)

write_csv(lago_final, "lago_env.csv")



