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
  scale_colour_manual(legend.title.climate, values = c("Black", "grey80"), labels = c("Removal", "Untreated")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ temp) +
  scale_shape_manual(legend.title.climate, values = c(16, 1), labels = c( "Removal", "Untreated")) +
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
  scale_colour_manual(legend.title.climate, values = c("Black", "grey80"), labels = c("Untreated", "Removal")) +
  scale_shape_manual(legend.title.climate, values = c(16, 1), labels = c("Untreated", "Removal")) +
  theme_classic() +
  labs(y = paste("Δ diversity (Shannon)")) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.text.x  = element_text(angle = 90))  +
  axis.dimLarge


delta <- plot_grid(cover_temp.orig, sla_temp.orig, labels = c('A', 'B'), ncol = 2, align = 'h')
ggsave(delta, width = 11, height = 4.5, dpi = 300, filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 1/figures/fig2ab.jpg")


abund <- my.GR.data %>% 
  filter(functionalGroup == "forb", temp > 10, TTtreat == "RTC") %>% 
  group_by(siteID, Year, species) %>% 
  mutate(meanCov = mean(cover)) %>% 
  ungroup() %>% 
  group_by(precip, species) %>%
  mutate(dominance = case_when(
    precip == "2700" & mean(meanCov) > 16 ~ "dominant",
    precip == "2000" & mean(meanCov) > 16 ~ "dominant",
    precip == "1200" & mean(meanCov) > 16 ~ "dominant",
    precip == "600" & mean(meanCov) > 16 ~ "dominant",
    precip == "2700" & mean(meanCov) < 16 ~ "subordinate",
    precip == "2000" & mean(meanCov) < 16 ~ "subordinate",
    precip == "1200" & mean(meanCov) < 16 ~ "subordinate",
    precip == "600" & mean(meanCov)  < 16 ~ "subordinate"
  )) %>% 
  #mutate(dominance = if_else(mean(meanCov, na.rm = TRUE) > 15, "dominant", "subordinate")) %>% 
  mutate(label = if_else(Year == 2016 & dominance == "dominant", as.character(species), NA_character_))


lm1 <- lm(meanCov ~ species + precip + 0, data = abund) %>% 
  tidy() %>% 
  mutate(term = gsub("species", "", term)) %>% 
  filter(!term == "precip")

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
  facet_grid(. ~ precip) +
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
  theme(axis.text.x  = element_text(angle = 90)) +
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
  theme(axis.text.x  = element_text(angle = 90)) +
  ggsave(filename = "moss_coverPRECIP.jpg", width = 7.5, height = 4, path = "/Users/fja062/Documents/seedclimComm/figures")


wholecom %>%
  filter(Year == 2011) %>% 
  ggplot(aes(x = wmeanSLA_local, colour = as.factor(temp), fill = as.factor(temp), shape = as.factor(temp))) +
geom_density(alpha = 0.3)  +
  scale_colour_manual("Mean summer\ntemperature (°C)", values = cbPalette) +
  scale_fill_manual("Mean summer\ntemperature (°C)", values = cbPalette) +
  facet_grid(.~functionalGroup) +
  theme_classic() +
  axis.dimLarge +
  labs(x = expression("Community weighted mean SLA "(cm^2/g)), y = "Frequency density")

ggsave(filename = paste0("figsupp1.jpg"), width = 11, height = 4.5, dpi = 300, path = "~/OneDrive - University of Bergen/Research/FunCaB/figures")

rtcLDMC <- forbcom %>%
  mutate(nYear = as.numeric(as.character(Year))) %>% 
  filter(!is.na(wmeanLDMC_local))

modLDMC <- rtcLDMC %>% 
  lmer(wmeanLDMC_local ~ TTtreat*scale(summer_temp)*scale(annPrecip)*scale(nYear) - TTtreat*scale(summer_temp)*scale(annPrecip)*scale(nYear) + (1|siteID), REML = FALSE, data = .)

LDMCpred <- predict(modLDMC, newdata = LDMCnewDat)

rtcLDMC <- rtcLDMC %>% 
  mutate(LDMCpredL = (LDMCpred - sqrt(wmeanLDMC_local)*1.96),
         LDMCpredH = (LDMCpred + sqrt(wmeanLDMC_local)*1.96)) %>% 
  group_by(Year, siteID, TTtreat) %>% 
  mutate(mLDMCpredL = mean(LDMCpredL),
         mLDMCpredH = mean(LDMCpredH))
  
rtcLDMC %>% 
  #gather(deltawmeanCN_local,deltawmeanheight_local,deltawmeanLA_local,deltawmeanLDMC_local,deltawmeanLTH_local,deltawmeanSLA_local, key = "trait", value = "value") %>% 
ggplot(aes(nYear, LDMCpred, colour = factor(temp), fill = factor(temp), group = factor(temp))) +
  geom_point(aes(nYear, wmeanLDMC_local)) +
  geom_line() +
  #stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.5), geom = "line") +
  geom_ribbon(aes(ymin = mLDMCpredL, ymax = mLDMCpredH), alpha=0.3) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  #geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", size = 1) #+
  facet_wrap(. ~ TTtreat, scales = "free_y")
