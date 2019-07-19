# source species composition compilation
source("~/Documents/FunCaB/bioticInteractions/composition/dataProcessing/00funcab_data_processing.R")

# source trait imputation
source("~/Documents/FunCaB/bioticInteractions/traits/trait_imputation.R")

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

# filter out spp covers in wrong removal plots
community_FD <- community_FD %>% 
  filter(!(Treatment %in% c("F", "FB", "FGB") & functionalGroup == "forb"),
         !(Treatment %in% c("G", "GB", "FGB") & functionalGroup == "graminoid"))

# join imputed traits with species cover, and filter for treatment
community_FD <- left_join(community_cover, Species_traits, by = c("speciesID", "species", "siteID")) %>%
  select(siteID, Treatment, blockID, turfID, Year, species, functionalGroup, cover, C, N, CN, SLA, Lth, LDMC, sqrtLA, logHeight) %>%
  filter(!is.na(cover), cover > 0)


# calculate community weighted means and variances
# https://rdrr.io/cran/TAM/man/weighted_Stats.html
community_FDWM <- community_FD %>% 
  group_by(siteID, blockID, turfID, Treatment, Year) %>%
  mutate(richness = sum(n_distinct(species))) %>%
  mutate(diversity = diversity(cover, index = "shannon")) %>%
  mutate(evenness = (diversity/log(richness))) %>% 
  group_by(siteID, blockID, turfID, Treatment, functionalGroup) %>%
  summarise(sumcover = sum(cover),
            richness = mean(richness),
            diversity = mean(diversity),
            evenness = mean(evenness),
            Wmean_LDMC = weighted_mean(LDMC, cover),
            Wmean_Lth= weighted_mean(Lth, cover),
            Wmean_LA= weighted_mean(sqrtLA, cover),
            Wmean_SLA= weighted_mean(SLA, cover),
            Wmean_Height= weighted_mean(logHeight),
            Wmean_CN = weighted_mean(CN, cover),
            Wmean_C = weighted_mean(C, cover),
            Wmean_N = weighted_mean(N, cover),
            #WmeanPCA = weighted_mean(PC1, cover),
            Wvar_LDMC= wt.var(LDMC, cover),
            Wvar_Lth = sqrt(wt.var(Lth, cover)),
            Wvar_LA  = sqrt(wt.var(sqrtLA, cover)),
            Wvar_SLA = wt.var(SLA, cover),
            Wvar_Height = wt.var(logHeight, cover),
            Wvar_C = sqrt(wt.var(C, cover)),
            Wvar_N = sqrt(wt.var(N, cover)),
            Wvar_CN = sqrt(wt.var(CN, cover))) %>% 
  ungroup()

# gather traits
community_FDWM <- community_FDWM %>% 
  gather(c(sumcover:Wvar_CN), key = "trait", value = "value") 


community_FDWM <- community_FDWM %>% 
  group_by(turfID, trait, siteID, functionalGroup) %>% 
  left_join(community_FDWM %>% filter(Treatment == "C") %>% ungroup() %>% select(Cvalue = value, siteID, blockID, trait)) %>%
  mutate(valueAnom = value - Cvalue,
         char = case_when(
           Treatment %in% c("B", "F", "G") ~ "effect",
           Treatment %in% c("FB", "GF", "GB") ~ "response"
         )) %>%
  filter(!Treatment %in% c("C", "FGB"))

#save(community_FDWM, file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/community_FD.RData")

# attach weather data
community_FDWM <- community_FDWM %>% 
  left_join(weather)


#------------ analyses ---------------#
#load packages
library(lme4)
library(MuMIn)
library(broom)
library(broom.mixed)

# Scaling explanatory variables except for richness (poisson distribution)
# relevel treatment so that TTC is the intercept
community_FD_analysis <- community_FDWM %>% 
  mutate(Sprecip0916 = as.numeric(scale(precip0916)),
         Stemp0916 = as.numeric(scale(temp0916)),
         sPrecip70 = as.numeric(scale(precip7010)),
         sTemp70 = as.numeric(scale(temp7010))) %>% 
  mutate(Treatment = factor(Treatment, levels = c("C", "FGB", "GF", "GB", "FB", "G", "F", "B"))) %>% 
  group_by(trait, functionalGroup) %>% 
  mutate(value = if_else(trait == "richness", value, scale(value))) %>%
  filter(is.finite(value))


#model of effect of graminoids and response of forbs
mod1temp <- community_FD_analysis %>% 
  filter(Treatment %in% c("C", "G", "B", "GB"), functionalGroup == "forb") %>% 
  mutate(traitII = trait) %>% 
  group_by(trait) %>%
  do({if(.$traitII[1] == "richness"){
    mod <- glmer(value ~ Treatment*sTemp70*sPrecip70 + (1|siteID/blockID), family = "poisson", data = .)
  } else {
    mod <- lmer(value ~ Treatment*sTemp70*sPrecip70 + (1|siteID/blockID), REML = FALSE, data = .)
    }
    tidy(mod)}) %>% 
  mutate(lower = (estimate - std.error*1.96),
         upper = (estimate + std.error*1.96)) %>% 
  arrange(desc(trait)) %>% 
  ungroup()

# model of effect of forbs and response of graminoids
mod2temp <- community_FD_analysis %>% 
  filter(Treatment %in% c("C", "F", "B", "FB"), functionalGroup == "graminoid") %>% 
  mutate(traitII = trait) %>%
  group_by(trait) %>% 
  do({if(.$traitII[1] == "richness"){
    mod <- glmer(value ~ Treatment*sTemp70*sPrecip70 + (1|siteID/blockID), family = "poisson", data = .)
  } else {
    mod <- lmer(value ~ Treatment*sTemp70*sPrecip70 + (1|siteID/blockID), REML = FALSE, data = .)
  }
    tidy(mod)}) %>% 
  mutate(lower = (estimate - std.error*1.96),
         upper = (estimate + std.error*1.96)) %>% 
  arrange(desc(trait)) %>% 
  ungroup()


mods <- bind_rows("forb" = mod1temp,"graminoid" = mod2temp, .id = "model") %>% 
  mutate(model) %>% 
  filter(!grepl("^sd_", term),
         !term == "(Intercept)") %>% 
  mutate(term = gsub("sPrecip70", "P", term),
         term = gsub("sTemp70", "t", term),
         term = gsub("Treatment", "", term),
         #term = gsub("(Intercept)", "bare ground", term),
         term = gsub(":", " x ", term)
  ) %>% 
  mutate(estimate = round(estimate, 3), std.error = round(std.error, 3), statistic = round(statistic, 3))
write_csv(mods, path = "~/OneDrive - University of Bergen/Research/FunCaB/paper 3/mods.csv")


modOutput <- mods %>% 
  filter(!grepl("Wvar", trait)) %>% 
  ggplot(aes(x = term, y = estimate, ymin = lower, ymax = upper)) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  #scale_shape_manual(legend.title.climate, values = c(21,24, 23)) + #, labels = c("0.6","2.7","6.5","10.5")
  scale_colour_manual(legend.title.climate, values = cbPalette) + #guide = guide_legend(reverse=TRUE)
  scale_x_discrete(limits = c("t", "P",   "C", "C x t", "C x P", "C x t x P", 
                              "F", "F x t", "F x P", "F x t x P",
                              "G", "G x t", "G x P", "G x t x P",
                              "B", "B x t", "B x P", "B x t x P",
                              "GF", "GF x t", "GF x P", "GF x t x P",
                              "GB", "GB x t", "GB x P", "GB x t x P",
                              "FB", "FB x t", "FB x P", "FB x t x P")) +
  geom_vline(xintercept =  c(3.5, 10.5, 17.5, 24.5, 31.5, 38.5, 45.5), colour = "grey90") +
  facet_wrap(trait ~ model, scales = "free_y", ncol = 2) +
  coord_flip() +
  axis.dimLarge +
  labs(y = "standardised coefficients", x = "Functional group")


ggsave(modOutput, filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 3/figures/supfig1.jpg", dpi = 300, height = 19, width = 13)


Glancemod1temp <- community_FD_analysis %>% 
  filter(Treatment %in% c("C", "G", "B", "GB"), functionalGroup == "forb") %>% 
  nest(-trait) %>%
  mutate(fit = map(data, ~ lmer(value ~ Treatment*Stemp0916*Sprecip0916 + (1|siteID/blockID), REML = FALSE, data = .)),
         results = map(fit, glance)) %>% 
  unnest(results)


#------------ PLOTTING --------------#
# source plotting colours etc
source("~/Documents/FunCaB/figures/plotting_dim.R")

g <- community_FDWM %>% 
  filter(!Treatment == "GF",
         !grepl("Wvar", trait)) %>%
 # mutate(Treatment = case_when(
 #  Treatment == "FB" ~ "G",
 #  Treatment == "GB" ~ "F"
 #)) %>% 
  group_by(trait) %>% 
  mutate(scaleAnom = scale(valueAnom)) %>% 
  filter(!is.na(functionalGroup)) %>% 
  ggplot(aes(x = factor(tempLevel), y = factor(precipLevel), fill = scaleAnom)) +
  geom_tile() +
  scale_fill_gradient2(low = pal1[3], mid = "snow1", high = pal1[4]) +
  facet_grid(trait~char*Treatment*functionalGroup, scales = "free")

ggsave(g, filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 3/figures/fig1PCA.jpg", dpi = 300, height = 4, width = 8)



#---------- abundance/dominance analyses and figures -----------#
abund <- comp2 %>% 
  left_join(weather) %>% 
  #filter(functionalGroup == "forb") %>% 
  group_by(siteID, Year, species, tempLevel, precipLevel, Treatment, functionalGroup) %>% 
  summarise(meanCov = mean(cover)) %>% 
  ungroup() %>% 
  group_by(precipLevel, tempLevel, species) %>%
  mutate(dominance = case_when(
    mean(meanCov) > 16 ~ "dominant",
    mean(meanCov)  < 16 ~ "subordinate"
  )) %>% 
  #mutate(dominance = if_else(mean(meanCov, na.rm = TRUE) > 15, "dominant", "subordinate")) %>% 
  mutate(label = if_else(dominance == "dominant", as.character(species), NA_character_))


# plot to illustrate which species are changing in abundance, both as a response but also as an effect - FORBS
abundplotForb <- abund %>% filter(Treatment %in% c("G", "B", "GB", "C"), !functionalGroup == "graminoid", !is.na(species)) %>% 
  left_join(abund %>% filter(Treatment == "C") %>% ungroup() %>% select(CmeanCov = meanCov, siteID, species)) %>%
  mutate(covAnom = meanCov - CmeanCov,
         char = case_when(
           Treatment %in% c("B", "F", "G") ~ "effect",
           Treatment %in% c("FB", "GF", "GB") ~ "response"
         )) %>% 
  filter(!Treatment == "C") %>% 
  ggplot(aes(x = Treatment, y = covAnom, colour = char)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  geom_label_repel(aes(label = label),
             na.rm = TRUE) +
  scale_colour_manual("", values = c("Black", "darkgreen")) +
  facet_grid(precipLevel ~ tempLevel, scales = "free_x")
  
ggsave(abundplotForb, filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 3/figures/fig2a.jpg", dpi = 300, height = 12, width = 8)


# plot to illustrate which species are changing in abundance, both as a response but also as an effect - GRAMINOIDS
abundplotGram <- abund %>% filter(Treatment %in% c("F", "B", "FB", "C"), !functionalGroup == "forb", !is.na(species)) %>% 
  left_join(abund %>% filter(Treatment == "C") %>% ungroup() %>% select(CmeanCov = meanCov, siteID, species)) %>%
  mutate(covAnom = meanCov - CmeanCov,
         char = case_when(
           Treatment %in% c("B", "F", "G") ~ "effect",
           Treatment %in% c("FB", "GF", "GB") ~ "response"
         )) %>% 
  filter(!Treatment == "C") %>% 
  ggplot(aes(x = Treatment, y = covAnom, colour = char)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  geom_label_repel(aes(label = label),
                   na.rm = TRUE) +
  scale_colour_manual("", values = c("Black", "darkgreen")) +
  facet_grid(precipLevel ~ tempLevel, scales = "free_x")

ggsave(abundplotGram, filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 3/figures/fig2b.jpg", dpi = 300, height = 12, width = 8)
