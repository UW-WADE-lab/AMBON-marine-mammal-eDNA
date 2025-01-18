####AMBON marine mammal biodiversity model
###AVC January 2025

### Set up environment ---------------------------------------------------------

library(tidyverse)
library(vegan)
library(ggvegan)


load("data products/AMBON_total_detections_MM.csv")
load("env-data/env_data_rasterdf.Rdata")

### Biodiversity index ---------------------------------------------------------

detect_data_wide <- detect_data_long %>% 
  dplyr::select(location1, Species, nReads) %>% 
  group_by(location1, Species) %>% 
  mutate(meanReads = mean(nReads)) %>% 
  slice_head() %>% 
  #unite(location1:depth, col = "Sample", sep = "_") %>% 
  pivot_wider(id_cols = location1, names_from = Species, values_from = meanReads) %>% 
  column_to_rownames("location1")
  

div <- diversity(detect_data_wide)

detect_data_count <- ifelse(detect_data_wide > 0, 1, 0)

sp_rich <- specnumber(detect_data_count)

### Merge with environmental data ----------------------------------------------

env_df_merge <- env_df %>% 
  mutate(lat_factor = as.factor(round(y, digits = 1))) %>% 
  mutate(lon_factor = as.factor(round(x, digits = 1))) %>% 
  group_by(lat_factor, lon_factor) %>% 
  summarise_all(mean)

sp_rich_merge <- as.data.frame(sp_rich) %>% 
  rownames_to_column("location1") %>% 
  left_join(detect_data_long %>% dplyr::select(location1, longitude, latitude)) %>% 
  distinct() %>% 
  mutate(lat_factor = as.factor(round(latitude, digits = 1))) %>% 
  mutate(lon_factor = as.factor(round(longitude, digits = 1))) %>% 
  left_join(env_df_merge, by = c("lon_factor" = "lon_factor",
                                 "lat_factor" = "lat_factor"))

### Model biodiversity distribution --------------------------------------------

sprich_full_model <- gam(sp_rich~#s(x, bs = "ts")+
                         #s(y, bs = "ts")+
                         s(bathy)+
                         s(distShore)+
                         s(slope)+
                         s(sst),
                       data = sp_rich_merge,
                       family = "poisson",
                       method = "REML")

visreg(sprich_full_model)

#check output
summary(sprich_full_model)
gam.check(sprich_full_model)
plot(sprich_full_model)
qqnorm(residuals(sprich_full_model))
qqline(residuals(sprich_full_model))

#density predictions
predictions =  terra::predict(data_extent, model = sprich_full_model, type = "response")
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

#numSpecs by station
numDetect_station <- as.data.frame(detect_data_count) %>% 
  rownames_to_column("location1") %>%
  pivot_longer(cols = -location1, names_to = "Species", values_to = "detect") %>% 
  filter(detect > 0) %>% 
  group_by(location1,Species) %>% 
  slice_head() %>% 
  group_by(location1) %>% 
  summarise(nSpec = sum(detect)) %>% 
  left_join(detect_data_long, by =c("location1" = "location1")) %>% 
  dplyr::select(location1,nSpec,latitude,longitude) %>% 
  distinct()

#plot
sprich_dist_map <- ggplot(bs_land)+
  geom_tile(data=pred_df, aes(x=x,y=y, fill=layer), alpha=0.7) +
  geom_sf(fill = "grey50", colour = NA) +
  #scale_fill_gradient(low = "blue", high = "red", na.value = "NULL") +
  scale_fill_viridis(option = "inferno", guide = "none") +
  scale_shape_manual(values = c(1,16), name = "Detection") +
  scale_color_manual(values = c("#0f85a0","#4a9152"), name = "Collection year")+
  theme_minimal() +
  #theme(legend.position = "none") +
  ggspatial::geom_spatial_point(data = numDetect_station, 
                                aes(x = longitude, y = latitude, 
                                    #color = as.factor(collection_year),
                                    size = nSpec),
                                alpha = 0.5) +
  scale_size(name="Number of species", range = c(0.5,6)) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())

save(sprich_dist_map, sprich_full_model, file = "data products/sprich_eDNA_density.Rdata")                   
