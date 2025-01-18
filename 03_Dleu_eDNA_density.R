###Dleu eDNA distribution model
###AVC Jan 2025

### set up environment ---------------------------------------------------------
library(tidyverse)

#spatial modeling
library(MuMIn)
library(mgcv)
library(visreg)
library(sdmpredictors)
library(mregions2)
library(raster)
library(terra)
library(sp)
library(rgdal)
library(viridis)

load("data products/AMBON_total_detections_MM.csv")
load("env-data/env_data_rasterdf.Rdata")

### Merge detections and env data ----------------------------------------------

env_df_merge <- env_df %>% 
  mutate(lat_factor = as.factor(round(y, digits = 1))) %>% 
  mutate(lon_factor = as.factor(round(x, digits = 1))) %>% 
  group_by(lat_factor, lon_factor) %>% 
  summarise_all(mean)

detect_data_dleu <- detect_data_long %>% 
  group_by(location1, depth, Species) %>% 
  mutate(totReads = sum(nReads)) %>% 
  mutate(detect = case_when(totReads > 1 ~ 1,
                            TRUE ~ 0)) %>% 
  slice_head() %>% 
  relocate(totReads:detect, .after = nReads) %>% 
  filter(Species == "beluga whale") %>% 
  filter(!is.na(sample_type)) %>% 
  mutate(lat_factor = as.factor(round(latitude, digits = 1))) %>% 
  mutate(lon_factor = as.factor(round(longitude, digits = 1))) %>% 
  left_join(env_df_merge, by = c("lon_factor" = "lon_factor",
                                 "lat_factor" = "lat_factor")) 
  # mutate(sst_collectionday = case_when(collection_month == 9~MS_sst09_5m,
  #                                      collection_month == 10~MS_sst10_5m,
  #                                      collection_month == 11~MS_sst11_5m,
  #                                      TRUE~NA))


### Model beluga density -------------------------------------------------------

dleu_full_model <- gam(detect~#s(x, bs = "ts")+
                         #s(y, bs = "ts")+
                         s(bathy)+
                         s(distShore)+
                         s(slope)+
                         s(sst),
                       data = detect_data_dleu,
                       family = "binomial",
                       method = "REML")

visreg(dleu_full_model)

#check output
summary(dleu_full_model)
gam.check(dleu_full_model)
plot(dleu_full_model)
qqnorm(residuals(dleu_full_model))
qqline(residuals(dleu_full_model))

#density predictions
predictions =  terra::predict(data_extent, model = dleu_full_model, type = "response")
predictions[is.na(predictions[])] <- 0
spplot(predictions, colorkey = list(space = "left") ,scales = list(draw = TRUE))

### Plot model predictions ----------------------------------------------------
#convert to dataframe
pred_spdf <- as(predictions, "SpatialPixelsDataFrame")
pred_df <- as.data.frame(pred_spdf)

#create coastline shapefile
world <- terra::vect("env-data/ne_10m_land/ne_10m_land.shp")
bs_land <- crop(world, predictions)
bs_land_sf <- st_as_sf(bs_land)

#detections by station
dleu_pos <- detect_data_dleu %>% 
  group_by(location1, depth) %>% 
  slice_head() %>% 
  group_by(location1) %>% 
  mutate(numDetect = sum(detect))

#plot
dleu_density_map <- ggplot(bs_land)+
  geom_tile(data=pred_df, aes(x=x,y=y, fill=layer), alpha=0.7) +
  geom_sf(fill = "grey50", colour = NA) +
  #scale_fill_gradient(low = "blue", high = "red", na.value = "NULL") +
  scale_fill_viridis(option = "inferno", guide = "none") +
  scale_shape_manual(values = c(1,16), name = "Detection") +
  scale_color_manual(values = c("#0f85a0","#4a9152"), name = "Collection year")+
  theme_minimal() +
  #theme(legend.position = "none") +
  ggspatial::geom_spatial_point(data = detect_data_dleu, 
                                aes(x = longitude, y = latitude, 
                                    color = as.factor(collection_year),
                                    shape = as.factor(detect)),
                                size = 3,
                                alpha = 0.5,
                                stroke = 1) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
  
save(dleu_density_map, dleu_full_model, file = "data products/Dleu_eDNA_density.Rdata")                   