#### AMBON marine mammal detections
#### Depth ridgeplot, sampling map
#### AVC January 2025

#### set up environment --------------------------------------------------

#general
library(tidyverse)
library(ggridges)
library(PNWColors)
#maps
library(ggOceanMaps)
library(scatterpie)
library(mapdata)


detect_data <- read.csv("../03 detection data/AMBON_MiFish_cetacean_detections_2021_2023.csv")

#### format data ----------------------------------------------------------

detect_data_long <- detect_data %>% 
  rename("bowhead whale" = "Balaena.mysticetus",
         "minke whale" = "Balaenoptera.acutorostrata",
         "fin whale" = "Balaenoptera.physalus",
         "beluga whale" = "Delphinapterus.leucas",
         "bearded seal" = "Erignathus.barbatus",
         "grey whale" = "Eschrichtius.robustus",
         "Pacific white-sided dolphin" = "Lagenorhynchus.obliquidens",
         "humpback whale" = "Megaptera.novaeangliae",
         "walrus" = "Odobenus.rosmarus",
         "harbor/spotted seal" = "Phoca",
         "ribbon seal" = "Phoca.fasciata",
         "harbor porpoise" = "Phocoena.phocoena",
         "ringed seal" = "Pusa.hispida") %>% 
  pivot_longer("bowhead whale":"ringed seal", 
               names_to = "Species", 
               values_to = "nReads") %>% 
  relocate(Species:nReads, .after = "ABL_ID") %>% 
  filter(sample_type == "sample")

detections_by_station <- detect_data_long %>% 
  mutate(detect = case_when(nReads > 0 ~ 1,
                            TRUE ~ 0)) %>% 
  relocate(detect, .after = nReads) %>% 
  group_by(location1, depth, Species) %>% 
  mutate(nDetect = sum(detect), nObs = n()) %>% 
  relocate(nDetect:nObs, .after = detect)

### Density ridgeplot ----------------------------------------------------------

positive_detections <- detections_by_station %>% 
  filter(detect == 1) %>% 
  filter(grepl("SKQ2021",ABL_ID))

rare_detections <- positive_detections %>% 
  group_by(Species) %>% 
  filter(n() < 3)
  
depth_detection <- ggplot(positive_detections, aes(y = Species, x = depth, 
                                                   fill = Species,
                                                   color = Species)) +
  geom_density_ridges(scale = 0.5, bandwidth = 15,
                      jittered_points = TRUE,
                      point_alpha = 1,
                      point_shape = 21,
                      alpha = 0.6) +
  theme_minimal() + 
  coord_flip() +
  scale_x_reverse() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = pnw_palette("Cascades",12, type = "continuous")) +
  scale_color_manual(values = pnw_palette("Cascades",12, type = "continuous")) +
  geom_point(data = rare_detections) +
  theme(legend.position = "none") 
  
#### Map of sampling locations -------------------------------------------------

detect_stations <- detect_data %>% 
  drop_na(longitude) %>% 
  group_by(location1, collection_year) %>% 
  slice_head()

sampling_map <- basemap(limits=c(-175,-150,60,75), bathymetry = TRUE, rotate = TRUE) +
  ggspatial::geom_spatial_point(data = detect_stations, 
                                aes(x = longitude, y = latitude, 
                                    color = as.factor(collection_year),
                                    shape = as.factor(collection_year)),
                                alpha = 0.4,
                                stroke = 1,
                                size = 3) +
  scale_color_manual(values = c("#0f85a0","#4a9152"))+
  ggspatial::stat_spatial_identity(position = "dodge") +
    theme(legend.position = "none")

### Detections by species -----------------------------------------------

positive_detect_species <- detections_by_station %>% 
  filter(detect == 1) %>% 
  group_by(ABL_ID, Species) %>% 
  slice_head()

species_detections <- positive_detect_species %>% 
  ungroup() %>% 
  group_by(Species) %>% 
  summarize(n())

positive_detect_species_map <- positive_detect_species %>% 
  group_by(location1, Species) %>% 
  mutate(nDetect = n()) %>% 
  dplyr::select(location1,Species, nDetect, latitude,longitude) %>% 
  group_by(location1) %>% 
  # mutate(totalDetect = sum(nDetect)) %>% 
  # mutate(propDetect = nDetect/totalDetect) %>% 
  # select(-nDetect, -totalDetect) %>% 
  group_by(location1, Species) %>% 
  slice_head() %>% 
  pivot_wider(names_from = Species, values_from = nDetect, values_fill = 0) %>% 
  mutate(total = rowSums(across(3:14))) %>% 
  mutate(total = total/25) %>% 
  mutate(total = case_when(total < 0.1~0.1,
                           TRUE~total)) %>% 
  mutate(longitude = longitude + 360)

ak<-map_data('world2Hires','USA:Alaska')

ggplot()+geom_polygon(data=ak,aes(long,lat,group=group),fill=8,color="black") +
  coord_map(xlim=c(-175, -155), ylim=c(60,72)) +
  geom_scatterpie(data = positive_detect_species_map,
                  aes(x=longitude, y=latitude, group=location1),
                  cols=unique(positive_detections$Species),
                              legend_name = "species",
                              sorted_by_radius = TRUE) + 

  scale_fill_manual(values = c(pnw_palette("Cascades",10, type = "continuous"),
                               pnw_palette("Shuksan",10, type = "continuous"))) +
  theme(legend.position = "none") +
  theme_minimal()


save(detect_data_long, file = "AMBON_total_detections_MM.csv")

save(depth_detection, sampling_map, file = "AMBON_sampling_detections.Rdata")
