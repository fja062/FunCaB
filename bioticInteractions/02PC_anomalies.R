# source species composition compilation
source("~/Documents/FunCaB/bioticInteractions/composition/dataProcessing/00funcab_data_processing.R")

# source trait imputation
source("~/Documents/FunCaB/bioticInteractions/traits/trait_imputation.R")

# load libraries
library(vegan)
library(ggvegan)
library(TAM)

# filter for funcab-only turfs
# create grouping variable before merge
comp2 <- comp2 %>% 
  filter(!Treatment == "XC", Year == 2017) %>%
  mutate(speciesID = paste0(siteID,"_", species))

# set covers and heights to zero in removed plots
community_cover <- comp2 %>% 
  mutate(mossCov = if_else(grepl("B", Treatment) & Year > 2015, 0, mossCov),
         forbCov = if_else(grepl("F", Treatment) & Year > 2015, 0, forbCov),
         graminoidCov = if_else(grepl("G", Treatment) & Year > 2015, 0, graminoidCov),
         vegetationHeight = if_else(Treatment == "FGB" & Year > 2015, 0, vegetationHeight),
         mossHeight = if_else(Treatment == "FGB" & Year > 2015, 0, mossHeight))



# join imputed traits with species cover, and filter for treatment
communityFD <- left_join(community_cover, Species_traits, by = c("speciesID", "species", "siteID")) %>%
  select(siteID, Treatment, blockID, turfID, Year, species, functionalGroup, cover, C, N, CN, SLA, Lth, LDMC, sqrtLA, logHeight) %>%
  gather(C:logHeight, key = trait, value = value) %>% 
  mutate(trStrat = case_when(
    trait %in% c("N", "SLA", "LDMC") ~ "LE",
    trait %in% c("logHeight", "sqrtLA", "Lth") ~ "Gro",
    TRUE ~ "Other"
  )) %>% 
  filter(!is.na(cover), cover > 0) %>% 
  spread(key = trait, value = value)

# filter out spp covers in wrong removal plots
communityFD <- communityFD %>% 
  filter(!(Treatment %in% c("F", "FB", "FGB") & functionalGroup == "forb"),
         !(Treatment %in% c("G", "GB", "FGB") & functionalGroup == "graminoid"))

# join to climate data
communityFD <- communityFD %>% 
  left_join(weather)

# pca scores for traits
spp_rda <- Species_traits %>% 
  gather(C:logHeight, key = "trait", value= "value") %>% 
  group_by(speciesID, trait) %>% 
  summarise(value = mean(value)) %>% 
  group_by(trait) %>% 
  mutate(value = scale(value)) %>% 
  spread(key = trait, value = value) %>% 
  select(C:sqrtLA, -CN) %>% 
  rda()

spp_rdaPlot <- autoplot(spp_rda, legend.position = "none", colour = "grey50") +
  axis.dimLarge
ggsave(spp_rdaPlot, filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 3/figures/supfig1.jpg", dpi = 300, height = 6, width = 6)

pca1 <- spp_rda$CA$u %>% 
  as.data.frame() %>% 
  select(PC1:PC4)

spp <- Species_traits %>% 
  gather(C:logHeight, key = "trait", value= "value") %>% 
  group_by(speciesID, trait) %>% 
  summarise(value = mean(value)) %>% 
  group_by(trait) %>% 
  mutate(value = scale(value)) %>% 
  spread(key = trait, value = value) %>% 
  select(speciesID) %>% 
  bind_cols(pca1)

#save(spp, file = "~/OneDrive - University of Bergen/Research/FunCaB/paper 3/pcaScores.RData")


communityFDrd <- communityFD %>% 
  mutate(speciesID = paste0(siteID, "_", species)) %>% 
  left_join(spp, by = "speciesID")

PCvariancePlot <- communityFDrd %>% 
  filter(!precipLevel == 2000,
         !Treatment  =="GF") %>% 
  ggplot((aes(Treatment, PC1, fill = factor(tempLevel)))) +
  geom_boxplot() +
  #geom_point() +
  scale_fill_manual(values = c("white", "grey80", "grey50")) +
  facet_grid(.~functionalGroup, scales = "free_x") +
  theme_classic() +
  labs(y = "PC axis 1 (23.3 %)") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  axis.dimLarge

ggsave(PCvariancePlot, filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 3/figures/supfig2.jpg", dpi = 300, height = 4.5, width = 8)


communityFDwm <- communityFD %>% 
  group_by(siteID, blockID, turfID, Treatment, Year) %>%
  mutate(richness = sum(n_distinct(species))) %>%
  mutate(diversity = diversity(cover, index = "shannon")) %>%
  mutate(evenness = (diversity/log(richness))) %>% 
  group_by(siteID, blockID, turfID, Treatment, functionalGroup) %>%
  summarise(sumcover = sum(cover),
            richness = mean(richness),
            diversity = mean(diversity),
            evenness = mean(evenness),
            Wmean_PC1 = weighted.mean(PC1, cover, na.rm = TRUE)) %>% 
  ungroup()

community_pcAnom %>% 
  filter(!Treatment == "GF") %>% 
  left_join(weather) %>% 
  filter(!is.na(functionalGroup)) %>% 
  ggplot(aes(x = factor(tempLevel), y = factor(precipLevel), fill = PCAnom)) +
  geom_tile() +
  scale_fill_gradient2(low = pal1[3], mid = "snow1", high = pal1[4]) +
  facet_grid(functionalGroup ~ Treatment)

community_pcAnom <- communityFDwm %>% 
  group_by(turfID, siteID, functionalGroup) %>% 
  left_join(communityFDwm %>% filter(Treatment == "C") %>% ungroup() %>% select(CWmean_PC1 = Wmean_PC1, siteID, blockID, functionalGroup)) %>%
  mutate(PCAnom = Wmean_PC1 - CWmean_PC1,
         char = case_when(
           Treatment %in% c("B", "F", "G") ~ "effect",
           Treatment %in% c("FB", "GF", "GB") ~ "response"
         )) %>%
  filter(!Treatment %in% c("C", "FGB"))

# attach weather data
community_pcAnom <- community_pcAnom %>% 
  left_join(weather)
