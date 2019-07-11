# source species composition compilation
source("~/OneDrive - University of Bergen/Research/FunCaB/SeedClim-Climate-Data/funcab/vegetation/00funcab_data_processing.R")

# source trait imputation
source("~/OneDrive - University of Bergen/Research/FunCaB/SeedClim-Climate-Data/funcab/vegetation/trait_imputation.R")

# source plotting code
source("~/OneDrive - University of Bergen/Research/FunCaB/SeedclimComm/inst/graminoidRemovals/plotting_dim.R")

# load packages
library("vegan")
library("ggvegan")

# create grouping variable before merge
community_cover_1516 <- comp2 %>%
  mutate(speciesID = paste0(siteID,"_", species))

# calculation of CWM and FDvar; Community weighted mean and community-weighted variance of trait values
# join imputed traits with species cover
# wt.var() calculate weighted variance


# merge traits with community data, and filter for treatment
community_FD <- left_join(community_cover_1516, Species_traits, by = c("speciesID", "siteID")) %>%
  select(siteID, Treatment, turfID, Year, species, functionalGroup, cover, C, N, CN, SLA, Lth, LDMC, sqrtLA, logHeight) %>%
  filter(!is.na(cover), Treatment %in% c("C", "XC"), cover > 0) %>% 
  filter(!turfID == "Lav1XC") # must fix - this is a problem in 2016 where Eric often didn't record cover


# calculate community weighted means and variances
community_FD <- community_FD %>% 
  group_by(siteID, turfID, Year) %>%
  mutate(richness = sum(n_distinct(species))) %>%
  mutate(diversity = diversity(cover, index = "shannon")) %>%
  mutate(evenness = (diversity/log(richness))) %>% 
  summarize(richness = mean(richness),
            diversity = mean(diversity),
            evenness = mean(evenness),
            Wmean_LDMC= weighted.mean(LDMC, cover, na.rm=TRUE),
            Wmean_Lth= weighted.mean(Lth, cover, na.rm=TRUE),
            Wmean_LA= weighted.mean(sqrtLA, cover, na.rm=TRUE),
            Wmean_SLA= weighted.mean(SLA, cover, na.rm=TRUE),
            Wmean_Height= weighted.mean(logHeight, cover, na.rm=TRUE),
            Wmean_CN = weighted.mean(CN, cover, na.rm=TRUE),
            Wmean_C = weighted.mean(C, cover, na.rm=TRUE),
            Wmean_N = weighted.mean(N, cover, na.rm=TRUE),
            Wvar_LDMC= wt.var(LDMC, cover),
            Wvar_Lth= wt.var(Lth, cover),
            Wvar_LA= wt.var(sqrtLA, cover),
            Wvar_SLA= wt.var(SLA, cover),
            Wvar_Height= wt.var(logHeight, cover),
            Wvar_C = wt.var(C, cover),
            Wvar_N = wt.var(N, cover),
            Wvar_CN = wt.var(CN, cover)) %>%
  #gather(key= Trait, value= value, -c(turfID:Year))%>%
  ungroup()

# add climate info
community_FD <- community_FD %>% 
  left_join(weather)

# filter away 2017
community_FD <- community_FD %>% 
  mutate(turfID = factor(turfID)) %>% 
  filter(!Year == 2017)


#------------ PLOTTING --------------


community_FD_spp <- community_FD %>% select(-(siteID:Year), -(tempLevel:precipLevel))

set.seed(52)
NMDS <- metaMDS(community_FD_spp, noshare = TRUE, try = 50)#DNC

fNMDS <- fortify(NMDS) %>% 
  filter(Score == "sites") %>%
  bind_cols(community_FD %>% select(siteID:Year, tempLevel:precipLevel))

sNMDS <- fortify(NMDS) %>% 
  filter(Score == "species")

treat_colours <- c("grey", "grey40", "orange", "purple")

g <- ggplot(fNMDS, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 3, aes(colour = factor(precipLevel), shape = factor(tempLevel))) +
  geom_text(data = sNMDS, aes(x = NMDS1, y = NMDS2, label = Label)) +
  coord_equal() +
  scale_shape_manual(values = c(24, 22, 23, 25)) +
  scale_colour_manual(values = pal2[c(7,4,1,3)]) +
  scale_fill_manual(values = pal2[c(7,4,1,3)]) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(colour = "Treatment", fill = "Treatment", shape = "Site", size = "Year") 

ggsave(g, filename = "~/OneDrive - University of Bergen/Research/FunCaB/figures/NMDStraits.jpg", dpi = 300)

set.seed(32)
HA <- two_sites_nmds("H", "A")
AM <- two_sites_nmds("A", "M")
ML <- two_sites_nmds("M", "L")
LM <- two_sites_nmds("L", "M") 

HA <- HA %>% mutate(Dim1 = -Dim1)
ML <- ML %>% mutate(Dim1 = -Dim1)
LM <- LM %>% mutate(Dim1 = -Dim1)


all_ord <- bind_rows(
  `H - A` = HA, 
  `A - M` = AM, 
  `M - L` = ML, 
  `L - ` = LM, .id = "which") %>% 
  mutate(which = factor(which, levels = c("H - A", "A - M", "M - L", "L - "), labels = c("High Alpine - Alpine", "Alpine - Middle", "Middle - Lowland", " - Lowland")))

OrdinationPlot <- g %+% all_ord +
  labs(x = "NMDS1", y = "NMDS2") +
  facet_wrap(~ which)

ggsave(OrdinationPlot, filename = "community/FinalFigures/OrdinationPlot.png", height = 7, width = 8, dpi = 300)
