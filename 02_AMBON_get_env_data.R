###AMBON environmental data
###AVC Jan 2025

### set up environment ---------------------------------------------------------
library(tidyverse)

#spatial modeling
library(sdmpredictors)
library(mregions2)
library(raster)
library(terra)
library(sp)
library(rgdal)
library(viridis)

load("data products/AMBON_total_detections_MM.csv")

# set environmental data directory
options(sdmpredictors_datadir = "~/env-data")

### Get environmental data -----------------------------------------------------

datasets <- list_datasets(terrestrial = FALSE, marine = TRUE)

MARSPEC <- list_layers(datasets[2,1])

env_data <- load_layers(c("MS_bathy_5m", "MS_biogeo05_dist_shore_5m", 
                          "MS_biogeo06_bathy_slope_5m", "MS_sst09_5m", 
                          "MS_sst10_5m", "MS_sst11_5m"))

data_extent <- raster::crop(env_data, extent(min(detect_data_dleu$longitude),
                                             max(detect_data_dleu$longitude),
                                             min(detect_data_dleu$latitude),
                                             max(detect_data_dleu$latitude)))

sst_rast <- mean(data_extent[[4:6]])

data_extent <- addLayer(data_extent[[1:3]], sst_rast)
names(data_extent) <- c("bathy", "distShore", "slope", "sst")
plot(data_extent)

data_correlations <- pearson_correlation_matrix(data_extent)
plot_correlation(data_correlations)

env_df <- as.data.frame(data_extent, xy=TRUE)

save(data_extent, env_df, env_data, file = "env-data/env_data_rasterdf.Rdata")
