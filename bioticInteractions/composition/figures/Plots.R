source("~/OneDrive - University of Bergen/Research/FunCaB/seedclimComm/inst/graminoidRemovals/plotting_dim.R")

#### COVER ####
# if we dig deeper, and look at the cover of the functional groups, things start to look interesting.

supfig1 <- wholecom %>% 
  filter(Year == 2011) %>% 
  ggplot(aes(x = factor(tempLevel), y = sumcover, fill = functionalGroup)) +
  geom_boxplot() +
  scale_fill_manual("",values = c("grey60", "grey90")) +
  facet_wrap(~precipLevel) +
  axis.dimLarge +
  labs(x = "Mean summer temperature (°C)", y = "Total vegetation cover")

ggsave(supfig1, file = "~/OneDrive - University of Bergen/Research/FunCaB/figures/supfig1.jpg")


###################

cover_temp.orig <- timedelta %>% 
  mutate(temp = if_else(grepl("6.5", tempLevel), "alpine", if_else(grepl("8.5", tempLevel), "sub-alpine", "boreal"))) %>% 
  mutate(temp = factor(temp, levels = c("alpine", "sub-alpine", "boreal"))) %>% 
  ggplot(aes(x = Year, y = deltasumcover, colour = TTtreat, shape = TTtreat, group = TTtreat)) +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.6), size = 0.75) +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.6), geom = "line", size = 0.75) +
  scale_colour_manual(legend.title.climate, values = c("grey80", "black"), labels = c("Removal", "Untreated")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ temp) +
  scale_shape_manual(legend.title.climate, values = c(1, 16), labels = c( "Removal", "Untreated")) +
  theme_classic() +
    labs(y = paste("Δ forb cover (%)")) +
  theme(legend.direction = "horizontal",
        legend.justification=c(-0.1, 1.5), 
        legend.position=c(0, 1),
        strip.background = element_blank(),
        axis.text.x  = element_text(angle = 90)) +
  axis.dimLarge


sla_temp.orig <- timedelta %>% 
  mutate(temp = if_else(grepl("6.5", tempLevel), "alpine", if_else(grepl("8.5", tempLevel), "sub-alpine", "boreal"))) %>% 
  mutate(temp = factor(temp, levels = c("alpine", "sub-alpine", "boreal"))) %>% 
  ggplot(aes(x = Year, y = deltadiversity, shape = TTtreat, colour = TTtreat, group = TTtreat)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.6), size = 0.75) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.6), geom = "line", size = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(. ~ temp) +
  #scale_colour_manual(legend.title.climate, values = c("Black", "grey80"), labels = c("Untreated", "Removal")) +
  #scale_shape_manual(legend.title.climate, values = c(16, 1), labels = c("Untreated", "Removal")) +
  theme_classic() +
  labs(y = paste("Δ diversity (Shannon)")) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.text.x  = element_text(angle = 90))  +
  axis.dimLarge


delta <- plot_grid(cover_temp.orig, sla_temp.orig, labels = c('A', 'B'), ncol = 2, align = 'h')
ggsave(delta, width = 11, height = 4.5, dpi = 300, filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 1/figures/fig2ab.jpg")


abund <- my.GR.data %>% 
  filter(functionalGroup == "forb", temp7010 > 10, TTtreat == "RTC") %>%
  group_by(siteID, Year, species) %>% 
  mutate(meanCov = mean(cover)) %>% 
  ungroup() %>% 
  group_by(precipLevel, species) %>%
  mutate(dominance = case_when(
    precipLevel == "2700" & mean(meanCov) > 20 ~ "dominant",
    precipLevel == "2000" & mean(meanCov) > 20 ~ "dominant",
    precipLevel == "1200" & mean(meanCov) > 20 ~ "dominant",
    precipLevel == "600" & mean(meanCov) >  20 ~ "dominant",
    precipLevel == "2700" & mean(meanCov) < 20 ~ "subordinate",
    precipLevel == "2000" & mean(meanCov) < 20 ~ "subordinate",
    precipLevel == "1200" & mean(meanCov) < 20 ~ "subordinate",
    precipLevel == "600" & mean(meanCov)  < 20 ~ "subordinate"
  )) %>% 
  #mutate(dominance = if_else(mean(meanCov, na.rm = TRUE) > 15, "dominant", "subordinate")) %>% 
  mutate(label = if_else(Year == 2016 & dominance == "dominant", as.character(species), NA_character_))


lm1 <- lm(meanCov ~ species + precipLevel + 0, data = abund) %>% 
  tidy() %>% 
  mutate(term = gsub("species", "", term)) %>% 
  filter(!term == "precipLevel")

summary(lm1)

abund %>%
  left_join(lm1, by = c("species" = "term")) %>% 
  ggplot(aes(x = Year, y = meanCov, colour = dominance, group = species)) +
  #stat_summary(fun.data = "mean_cl_boot", geom = "line") +
  geom_line() +
  geom_label(aes(label = label),
             nudge_x = -0.6,
             nudge_y = -2.45,
             na.rm = TRUE) +
  scale_colour_manual(values = c("Black", "grey80")) +
  facet_grid(. ~ precipLevel) +
  geom_point(size = 3) +
  labs(y = "Mean cover (%)") +
  axis.dimLarge +
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave(filename = paste0("fig4.jpg"), width = 11, height = 4.5, dpi = 300, path = "~/OneDrive - University of Bergen/Research/FunCaB/figures")


my.GR.data %>%
  mutate(mossHeight = if_else(Year == 2016, mossHeight/10, mossHeight)) %>% 
  ggplot(aes(x = Year, y = mossHeight, colour = TTtreat, shape = TTtreat)) +
  stat_summary(fun.data = "mean_cl_boot") +
  stat_summary(fun.data = "mean_cl_boot", geom = "line") +
  scale_colour_manual(legend.title.weat, values = cbPalette, labels = c("Control", "Removal")) +
  scale_shape_manual(legend.title.weat, values = c(1, 16), labels = c("Control", "Removal")) +
  theme_classic() +
  facet_grid(.~precip) +
  axis.dimLarge +
  labs(x = "", y = "Moss depth (cm)") +
  theme(axis.text.x  = element_text(angle = 90))
  ggsave(filename = "moss_depth_precip.jpg", path = "/Users/fja062/Documents/seedclimComm/figures", height = 4, width = 7.5)


my.GR.data %>% 
  ggplot(aes(x = Year, y = totalBryophytes, colour = TTtreat, shape = TTtreat)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.5), geom = "line") +
  facet_grid(. ~ precip) +
  scale_colour_manual(legend.title.weat, values = cbPalette, labels = c("Control", "Removal")) +
  scale_shape_manual(legend.title.weat, values = c(1, 16), labels = c("Control", "Removal")) +
  axis.dimLarge +
  theme_classic() +
  theme(axis.text.x  = element_text(angle = 90))
  ggsave(filename = "moss_coverPRECIP.jpg", width = 7.5, height = 4, path = "/Users/fja062/Documents/seedclimComm/figures")


wholecom %>%
  filter(Year == 2011) %>% 
  ggplot(aes(x = wmeanSLA, colour = as.factor(tempLevel), fill = as.factor(tempLevel), shape = as.factor(tempLevel))) +
geom_density(alpha = 0.3)  +
  scale_colour_manual("Mean summer\ntemperature (°C)", values = cbPalette) +
  scale_fill_manual("Mean summer\ntemperature (°C)", values = cbPalette) +
  facet_grid(.~functionalGroup) +
  theme_classic() +
  axis.dimLarge +
  labs(x = expression("Community weighted mean SLA "(cm^2/g)), y = "Frequency density")

ggsave(filename = paste0("figsupp1.jpg"), width = 11, height = 4.5, dpi = 300, path = "~/OneDrive - University of Bergen/Research/FunCaB/figures")

# Supplementary Figure 1
# graminoid and non-graminoid plant cover
wholecom %>% 
  filter(Year == 2011) %>% 
  ggplot(aes(x = factor(tempLevel), y = sumcover, fill = functionalGroup)) + 
  geom_boxplot() + 
  facet_grid(.~precipLevel, switch = "x") +
  scale_fill_grey() + 
  theme_classic() +
  labs(y = "Average cover (%)", x = "Mean summer temperature (ªC) and mean annual precipitation (mm/yr)") +
  axis.dim +
  theme(legend.title = element_blank(),
        legend.position = c(0.9,0.9),
        strip.placement = "outside",
        strip.background = element_rect(colour = "white"))

ggsave(filename = "supFig1.jpg", width = 8, height = 4, dpi = 300, path = "~/OneDrive - University of Bergen/Research/FunCaB/paper 1/figures/")

wholecom %>% 
  select(functionalGroup, sumcover, Year, tempLevel, precipLevel, turfID) %>% 
  pivot_wider(values_from = "sumcover", names_from = "functionalGroup") %>%
  rowwise() %>% 
  mutate(percentGram = (graminoid/(forb + graminoid))*100) %>% 
  filter(Year == 2011) %>% 
  ggplot(aes(x = factor(tempLevel), y = percentGram)) +
  facet_grid(.~precipLevel, switch =  "x") +
  geom_boxplot() +
  scale_fill_grey() + 
  theme_classic() +
  labs(y = "Average cover (%)", x = "Mean summer temperature (ºC) and mean annual precipitation (mm/yr)") +
  axis.dim +
  theme(legend.title = element_blank(),
        legend.position = c(0.9,0.9),
        strip.placement = "outside",
        strip.background = element_rect(colour = "white"))


# Supplementary Figure 2
# average number of species per trait and per plot
avSpN <- my.GR.data %>%
  filter(functionalGroup == "forb") %>% 
  gather(c(SLA_mean:N_mean), key = "trait", value = "value") %>% 
  group_by(ID, trait) %>% 
  mutate(nSpecies = n_distinct(species), nSpeciesWD = n_distinct(value, na.rm = TRUE),
         nSpeciesWDperc = (nSpeciesWD/nSpecies)*100) %>% 
  separate(trait, c("trait", "mean"), sep = "_") %>% 
  select(-mean)

avSpN %>% group_by(trait) %>% 
  summarise(meanPercent = mean(nSpeciesWDperc),
            minPercent = min(nSpeciesWDperc),
            maxPercent = max(nSpeciesWDperc),
            meanSpp = mean(nSpecies)) %>% 
  ungroup()

avSpN <- avSpN %>% group_by(tempLevel, precipLevel, trait) %>% 
  summarise(meanPercent = mean(nSpeciesWDperc),
            minPercent = min(nSpeciesWDperc),
            maxPercent = max(nSpeciesWDperc),
            meanSpp = mean(nSpecies)) %>% 
  ungroup()



avSpN %>% ggplot(aes(x = trait, y = meanPercent, fill = factor(tempLevel))) + 
  geom_boxplot() +
  scale_fill_grey() +
  theme_classic() +
  labs(y = "mean percentage of species with available trait data (%)") +
  theme(axis.title.x = element_blank())

ggsave(filename = "supFig2.jpg", width = 8, height = 4, dpi = 300, path = "~/OneDrive - University of Bergen/Research/FunCaB/paper 1/figures/")
