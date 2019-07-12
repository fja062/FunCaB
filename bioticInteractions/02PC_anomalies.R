# source species composition compilation
source("~/OneDrive - University of Bergen/Research/FunCaB/SeedClim-Climate-Data/funcab/vegetation/00funcab_data_processing.R")

# source trait imputation
source("~/OneDrive - University of Bergen/Research/FunCaB/SeedClim-Climate-Data/funcab/vegetation/trait_imputation.R")

# load libraries
library(vegan)
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


# calculation of CWM and FDvar; Community weighted mean and community-weighted variance of trait values
# join imputed traits with species cover, and filter for treatment
community_pc <- left_join(community_cover, Species_traits, by = c("speciesID", "species", "siteID")) %>%
  select(siteID, Treatment, blockID, turfID, Year, species, functionalGroup, cover, C, N, CN, SLA, Lth, LDMC, sqrtLA, logHeight) %>%
  filter(!is.na(cover), cover > 0)


# pca scores for traits
spp_rda <- Species_traits %>% 
  gather(C:logHeight, key = "trait", value= "value") %>% 
  group_by(species, trait) %>% 
  summarise(value = mean(value)) %>% 
  group_by(trait) %>% 
  mutate(value = scale(value)) %>% 
  spread(key = trait, value = value) %>% 
  select(C:sqrtLA, -CN) %>% 
  rda()

spp_rdaPlot <- plot(spp_rda, display = "species")

pca1 <- spp_rda$CA$u %>% 
  as.data.frame() %>% 
  select(PC1:PC4)

spp <- Species_traits %>% 
  gather(C:logHeight, key = "trait", value= "value") %>% 
  group_by(species, trait) %>% 
  summarise(value = mean(value)) %>% 
  group_by(trait) %>% 
  mutate(value = scale(value)) %>% 
  spread(key = trait, value = value) %>% 
  select(species) %>% 
  bind_cols(pca1)

save(spp, file = "~/OneDrive - University of Bergen/Research/FunCaB/paper 3/pcaScores.RData")


community_pc <- community_pc %>% left_join(spp, by = "species")

PCvariancePlot <- community_pc %>% filter(Treatment %in% c("C", "GB", "FB")) %>% ggplot((aes(PC1, fill = interaction(functionalGroup, Treatment)))) +
  geom_density(alpha = 0.5) +
  facet_grid(.~tempLevel)

ggsave(PCvariancePlot, filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 3/figures/fgVarPlotPCtemp.jpg", dpi = 300, height = 6, width = 13)


community_pcWM <- community_pc %>% 
  group_by(siteID, blockID, turfID, Treatment, Year) %>%
  mutate(richness = sum(n_distinct(species))) %>%
  mutate(diversity = diversity(cover, index = "shannon")) %>%
  mutate(evenness = (diversity/log(richness))) %>% 
  group_by(siteID, blockID, turfID, Treatment, functionalGroup) %>%
  summarise(sumcover = sum(cover),
            richness = mean(richness),
            diversity = mean(diversity),
            evenness = mean(evenness),
            Wmean_PC1 = weighted_mean(PC1, cover)) %>% 
  ungroup()

community_pcAnom <- community_pcWM %>% 
  group_by(turfID, siteID, functionalGroup) %>% 
  left_join(community_FD %>% filter(Treatment == "C") %>% ungroup() %>% select(CWmean_PC1 = Wmean_PC1, siteID, blockID, trait)) %>%
  mutate(PCAnom = Wmean_PC1 - CWmean_PC1,
         char = case_when(
           Treatment %in% c("B", "F", "G") ~ "effect",
           Treatment %in% c("FB", "GF", "GB") ~ "response"
         )) %>%
  filter(!Treatment %in% c("C", "FGB"))

# attach weather data
community_pcAnom <- community_pcAnom %>% 
  left_join(weather)