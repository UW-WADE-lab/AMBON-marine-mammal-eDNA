###Species rarefaction
###AVC Jan 2024

### Set up environment ---------------------------------------------------------

library(tidyverse)
library(vegan)

load("data products/AMBON_total_detections_MM.csv")

### Format data by tech rep ----------------------------------------------------

detect_by_rep <- detect_data_long %>% 
  filter(!(is.na(replicate))) %>% 
  mutate(detect = case_when(nReads > 0~1,
                            TRUE~0)) %>% 
  dplyr::select(ABL_ID, depth, Species, replicate, detect) %>% 
  pivot_wider(id_cols = c(ABL_ID, replicate, depth), names_from = Species, values_from = detect) %>%  
  #unite(ABL_ID:replicate, col = "SampleID") %>% 
  #column_to_rownames("SampleID") %>% 
  group_by(replicate) %>% 
  summarize_at(unique(detect_data_long$Species),sum)

sample.data <- detect_by_rep %>% 
  column_to_rownames("replicate") 

enviro.data <- detect_by_rep %>% 
  mutate(rowname = replicate) %>% 
  column_to_rownames("rowname")

### Estimate species accumulation by technical replicate -----------------------
#aggregate across all depths

Accum.1 <- specaccum(sample.data, y=enviro.data, 
                     method='exact', conditioned=FALSE, plotit=FALSE)


plot(Accum.1)
spec_accum_techrep <- cbind(Accum.1[["sites"]], Accum.1[["richness"]],
                            Accum.1[["sd"]]) %>% 
  as.data.frame() %>% 
  rename("sites" = "V1", "richness" = "V2", "sd" = "V3")

#fit glm to rarefaction data
test <- as.data.frame(cbind(richness = c(5,8,9), sites = c(1,2,3)))
pred_rare_tech <- lm(richness ~ log(sites), 
                      data=test)


summary(pred_rare_tech)
coef(pred_rare_tech)

plot(test$sites, test$richness)
curve(coef(pred_rare_tech)[1]+coef(pred_rare_tech)[2]*log(x))

predictions_analog <- data.frame(x = seq(1,10,1)) %>% 
  mutate(y = coef(pred_rare_tech)[1]+coef(pred_rare_tech)[2]*log(x))


### Plot species accumulation by technical replicate ---------------------------

rare_techrep_plot <- ggplot(data = predictions_analog, aes(x=x,y=y)) +
  geom_line(color = pnw_palette("Cascades",6)[4],
            linetype = "dotdash") +
  geom_line(data=spec_accum_techrep, 
            aes(x = sites, y = richness), inherit.aes = FALSE,
            size=1,color = pnw_palette("Cascades",6)[5]) +
  geom_point(data=spec_accum_techrep, 
             aes(x = sites, y = richness), inherit.aes = FALSE,
             size=2, color = pnw_palette("Cascades",6)[5]) +
  geom_ribbon(data=spec_accum_techrep, aes(x = sites, ymax = richness+sd, ymin = richness-sd),
              alpha=0.2, show.legend=FALSE,fill = pnw_palette("Cascades",6)[5],
              inherit.aes = FALSE) + 
  labs(x = "Technical Replicates", y = "Mean Species Richness") +
  theme_minimal() +
  scale_x_continuous(breaks=seq(1, 10, 1))

### Format data by bio rep -----------------------------------------------------

detect_by_biorep <- detect_data_long %>% 
  #filter(!(is.na(replicate))) %>% 
  mutate(detect = case_when(nReads > 0~1,
                            TRUE~0)) %>% 
  dplyr::select(ABL_ID, location1, depth, Species, replicate, detect) %>% 
  filter(depth < 200) %>% 
  mutate(depth_bin = cut(depth, breaks = seq(0,100,10), labels = FALSE)) %>% 
  dplyr::select(-depth) %>% 
  filter(!is.na(depth_bin)) %>% 
  pivot_wider(id_cols = c(ABL_ID, replicate, depth_bin), 
              names_from = Species, values_from = detect) %>%  
  #unite(ABL_ID:replicate, col = "SampleID") %>% 
  #column_to_rownames("SampleID") %>% 
  group_by(depth_bin) %>% 
  summarize_at(unique(detect_data_long$Species),sum)

sample.data <- detect_by_biorep %>% 
  column_to_rownames("depth_bin") 

enviro.data <- detect_by_biorep %>% 
  mutate(rowname = depth_bin) %>% 
  column_to_rownames("rowname")

### Estimate species accumulation by biological replicate ----------------------
#aggregate across all replicates

Accum.1 <- specaccum(sample.data, y=enviro.data, 
                     method='exact', conditioned=FALSE, plotit=FALSE)


plot(Accum.1)
spec_accum_biorep <- cbind(Accum.1[["sites"]], Accum.1[["richness"]],
                            Accum.1[["sd"]]) %>% 
  as.data.frame() %>% 
  rename("sites" = "V1", "richness" = "V2", "sd" = "V3")

### Plot species accumulation by technical replicate ---------------------------

rare_biorep_plot <- ggplot(data=spec_accum_biorep, 
                            aes(x = sites, y = richness, 
                                ymax = richness+sd, ymin = richness-sd)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(size=1,color = pnw_palette("Cascades",6)[5]) +
  geom_point(size=2, color = pnw_palette("Cascades",6)[5]) +
  geom_ribbon(alpha=0.2, show.legend=FALSE,fill = pnw_palette("Cascades",6)[5]) + 
  labs(x = "Biological Replicates", y = "Mean species richness") +
  theme_minimal() +
  scale_x_continuous(breaks=seq(1, 10, 1))

save(rare_techrep_plot, rare_biorep_plot, file = "data products/rarefaction_plots.Rdata")
