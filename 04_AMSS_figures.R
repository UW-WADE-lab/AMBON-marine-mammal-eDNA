###AMBON AMSS 2025 Make Figures
###AVC Jan 2025

### Set up environment ---------------------------------------------------------

library(tidyverse)
library(patchwork)

load("data products/Phis_eDNA_density.Rdata")
load("data products/Dleu_eDNA_density.Rdata")
load("data products/Erob_eDNA_density.Rdata")
load("data products/rarefaction_plots.Rdata")
load("data products/sprich_eDNA_density.Rdata")
load("data products/AMBON_sampling_detections.Rdata")


### Density Figure -------------------------------------------------------------

sp_dist <- dleu_density_map + erob_density_map + 
  phis_density_map + sprich_dist_map +
  plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = 'A')

pdf("figures/sp_dist.pdf")
plot(sp_dist)
dev.off()

### Rarefaction Figure ---------------------------------------------------------

rare_plots <- rare_techrep_plot + rare_biorep_plot

detection_pars <- rare_plots / depth_detection

png("figures/detection_pars.png", width = 2000, height = 1500)
plot(detection_pars)
dev.off()

pdf("figures/detection_pars.pdf")
plot(detection_pars)
dev.off()

### Sampling Map ---------------------------------------------------------------

sample_map <- sampling_map +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

png("figures/sampling_map.png", width = 1500, height = 1800)
plot(sample_map)
dev.off()

pdf("figures/sampling_map.pdf")
plot(sample_map)
dev.off()
