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
  select(siteID, Treatment, blockID, turfID, Year, species, functionalGroup, cover, N, SLA, LDMC, Lth, logHeight)

# filter out spp covers in wrong removal plots
communityFD <- communityFD %>% 
  filter(!Treatment == "FGB" | is.na(functionalGroup)) %>%
  filter(!Treatment == "F"   | !functionalGroup == "forb") %>%
  filter(!Treatment == "FB"  | !functionalGroup == "forb") %>%
  filter(!Treatment == "GB"  | !functionalGroup == "graminoid") %>%
  filter(!Treatment == "G"   | !functionalGroup == "graminoid")

# join to climate data
communityFD <- communityFD %>% 
  left_join(weather)


#save(communityFD, file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/community_FD.RData")

load("~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/community_FD.RData")

communityFDwm <- communityFD %>% 
  group_by(siteID, blockID, turfID, Treatment, Year) %>%
  mutate(richness = sum(n_distinct(species))) %>%
  mutate(diversity = diversity(cover, index = "shannon")) %>%
  mutate(evenness = (diversity/log(richness))) %>% 
  filter(!precipLevel == 2000) %>% 
  group_by(siteID, blockID, turfID, Treatment, functionalGroup) %>%
  summarise(sumcover = sum(cover),
            richness = mean(richness),
            diversity = mean(diversity),
            evenness = mean(evenness),
            Wmean_SLA = weighted.mean(SLA, cover, na.rm = TRUE),
            Wmean_Lth = weighted.mean(Lth, cover, na.rm = TRUE),
            Wmean_N = weighted.mean(N, cover, na.rm = TRUE),
            Wmean_logHeight = weighted.mean(logHeight, cover, na.rm = TRUE),
            ) %>% 
  ungroup()


B.flux_CI <- B.flux_CI %>% ungroup() %>% 
  filter(Cflux == "GPP") %>% 
  left_join(dict_Site, by = c("Site" = "old")) %>% 
  select(siteID = new, blockID = Block, Treatment, turfID, CI, bryophyteCov, forbCov, graminoidCov, flux)

trComp <- communityFDwm %>% 
  left_join(B.flux_CI) %>% 
  filter(!Treatment %in% c("C", "FGB")) %>% 
  left_join(weather) %>% 
  mutate(sWmean_SLA = scale(Wmean_SLA/100, scale = FALSE, center = TRUE),
         precipdiv = precip7010/1000,
         sprecip7010 = scale(precipdiv, scale = FALSE, center = TRUE),
         stemp7010 = scale(temp7010, scale = FALSE, center = TRUE)) %>% 
  filter(!Treatment  == "F" | bryophyteCov  > 1,
         !Treatment  == "G" | bryophyteCov  > 1,
         !Treatment  == "GF" | bryophyteCov > 1
         ) %>% 
         mutate(bryophyteCov = bryophyteCov/100)

trComp %>% ggplot(aes(x = bryophyteCov, y = CI, colour = Treatment)) +
  geom_point()
trComp %>% filter(Treatment == "GF") %>% ggplot(aes(x = bryophyteCov, y = CI)) +
  geom_point()

# 1st analysis: CLIMATE
lm1 <- trComp %>% 
  filter(!Treatment == "GF") %>% 
  group_by(functionalGroup) %>% 
  nest(.key = "dat") %>% 
  mutate(mod = map(dat, ~lmer(CI ~ stemp7010 + sprecip7010 + Treatment + (1|siteID), data = .)),
         coef = map(mod, tidy)) %>% 
  unnest(coef)
 


#### Fig of Controls vs treated in 2015 for cover to check for similarity.





# 2nd analysis: FORB
trComp %>% filter(Treatment %in% c("GB", "G")) %>% 
  ggplot(aes(x = temp7010, y = CI, colour = Treatment)) +
  geom_point()

trComp %>% filter(Treatment %in% c("FB", "F")) %>% 
  ggplot(aes(x = Wmean_SLA, y = CI, colour = Treatment)) +
  geom_point()

trComp %>% filter(Treatment %in% c("GB", "G")) %>% 
  lmer(CI ~ sWmean_SLA + Wmean_Lth + Wmean_logHeight + Wmean_N + bryophyteCov + (1|siteID), data = .) %>% 
  tidy()

# 3rd analysis: GRAMINOID
trComp %>% filter(Treatment %in% c("FB", "F")) %>% 
  ggplot(aes(x = sWmean_SLA, y = CI, colour = Treatment)) +
  geom_point()

trComp %>% filter(Treatment %in% c("FB", "F")) %>% 
  lmer(CI ~ sWmean_SLA + Wmean_Lth + Wmean_logHeight + Wmean_N + bryophyteCov + (1|siteID), data = .) %>% 
  tidy()



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