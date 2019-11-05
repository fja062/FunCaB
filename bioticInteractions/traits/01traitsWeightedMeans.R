#----- TRAITS AND WEIGHTED MEANS ----#
# source Ragnhild's trait data
source("~/Documents/SeedclimComm/ragnhild_trait_data/load_traits.R") # warning here is fine, it just means those spp didn't have CN data collected

traits <- tbl(con, "taxon") %>% 
  collect() %>% 
  left_join(tbl(con, "moreTraits"), copy = TRUE, by = "species") %>% 
  select(species, seedMass)

# adding traits to my.GR.data
composition <- comp2 %>%
  left_join(traitdata, by = c("species", "siteID"))
  #left_join(traits, by = "species")


library(vegan)

composition <- composition %>%
  group_by(turfID, Year) %>%
  mutate(richness = sum(n_distinct(species))) %>% 
  mutate(diversity = diversity(cover, index = "shannon")) %>% 
  mutate(evenness = (diversity/log(richness))) %>% 
  #filter(!is.na(cover)) %>%
  group_by(turfID, siteID, Year, functionalGroup) %>% 
  mutate(wmH = weighted.mean(Height_mean, cover, na.rm = TRUE),
         wmSLA = weighted.mean(SLA_mean, cover, na.rm = TRUE),
         wmLA = weighted.mean(LA_mean, cover, na.rm = TRUE),
         wmLDMC = weighted.mean(LDMC_mean, cover, na.rm = TRUE),
         wmLTH = weighted.mean(Lth_mean, cover, na.rm = TRUE),
         wmCN = weighted.mean(CN_mean, cover, na.rm = TRUE)) %>% # to here for compensation analysis
  filter(Treatment %in% c("FB", "GB")) %>% 
  filter(!(functionalGroup == "forb" & Treatment == "FB"), !(functionalGroup == "graminoid" & Treatment == "GB"))
  select(-Height_mean, -LA_mean, -SLA_mean, -LTH_mean, -LDMC_mean, -CN_mean, -TTtreat, -species, -cover) %>% 
  distinct(turfID, Year, functionalGroup, .keep_all = TRUE) %>% 
  group_by(turfID, siteID, Year) %>% 
  spread(key =functionalGroup, value = wmH) %>% 
  mutate(mossHeight = if_else(grepl("B", Treatment), 0, mossHeight),
         forb = if_else(grepl("F", Treatment), 0, forb),
         graminoid = if_else(grepl("G", Treatment), 0, graminoid)) %>% 
  ungroup()

save(composition, file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/cleanedVegComp.RData")


composition %>% left_join(weather) %>% 
  filter(Year == 2015) %>% 
  ggplot(aes(x = wmLA, fill = functionalGroup)) + 
  geom_density(alpha = 0.3) + facet_grid(.~tempLevel) +
  geom_rug() +
  scale_fill_manual("", values = c("goldenrod", "grey")) +
  labs(x = "Community-weighted mean leaf area") +
  axis.dimLarge

ggsave(file = "~/OneDrive - University of Bergen/Research/FunCaB/paper 2/figures/supFig13.jpg", dpi = 300, width = 8, height = 4)
