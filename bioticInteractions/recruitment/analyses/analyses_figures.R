# analyses and figures
hist(seedTot$seed)

glmSeed <- glmmTMB(seedSum ~ Treatment*Round*stemp7010*sprecip7010 - Treatment:Round:stemp7010:sprecip7010 - Treatment:stemp7010:sprecip7010 + (1|siteID/blockID), 
                   family = "poisson", 
                   zi = ~ Treatment + Round + sprecip7010 + stemp7010, 
                   data = seedTot)

simulationOutput <- simulateResiduals(glmSeed, n = 250)
plot(simulationOutput)
testZeroInflation(simulationOutput)
summary(glmSeed)

glmSeed <- tidy(glmSeed) %>% 
  mutate(lower = (estimate - std.error*1.96),
         upper = (estimate + std.error*1.96))

glmSeed <- glmSeed %>% 
  mutate(term = gsub("Treatment", "", term),
         term = gsub(":", " x ", term),
         term = gsub("Round2", "S", term),
         term = gsub("sprecip7010", "P", term),
         term = gsub("stemp7010", "t", term),
         term = gsub("FGB", "Gap", term),
         term = gsub("GB", "F", term),
         term = gsub("B", "GF", term),
         term = gsub("G", "FB", term),
         term = gsub("FBap", "Gap", term),
         term = gsub("FBF", "GF", term)
  ) %>% 
  filter(component == "cond", effect == "fixed") %>% 
  select(-group, -component) %>% 
  mutate(estimate = round(estimate, 3), std.error = round(std.error, 3), statistic = round(statistic, 3), p.value = round(p.value, 3))

glmSeed %>% ggplot(aes(x = term, y = estimate, ymin = lower, ymax = upper)) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #geom_vline(xintercept =  c(2.5, 6.5, 10.5), colour = pal1[c(5,2,4,5,2,4)], size = 27, alpha = 0.2) +
  geom_vline(xintercept = c(6.5, 14.5, 21.5), colour = "olivedrab4", size = 36.5, alpha = 0.2) +
  geom_point(position = position_dodge(width = 0.1), size = 3, fill = "grey") +
  scale_colour_manual("", values = pal1[c(2,4,5,1)]) + 
  scale_shape_manual("", values = c(21, 22, 23, 24, 25,1), limits = c("Graminoids", "Forbs", "Bryophytes", "Intact vegetation", "Climate")) + 
  scale_x_discrete(limits = c("Gap x S x t", "F x S x t", "FB x S x t", "GF x S x t",
                              "Gap x S x P", "F x S x P", "FB x S x P", "GF x S x P",
                              "Gap x t", "F x t", "FB x t", "GF x t",
                              "Gap x P", "F x P", "FB x P", "GF x P",
                              "t", "P", "S", "Gap", "F", "FB", "GF")) +
  coord_flip() +
  axis.dimLarge +
  labs(y = "standardised coefficients") +
  theme_classic() +
  axis.dimLarge +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom")

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig1Coef.jpg", dpi = 300, width = 5.5, height = 7)

totMods <- bind_rows("P" = glmSeedP, "T" = glmSeedT, .id = "model")


write_csv(glmSeed, path = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/data/totMods.csv")


hist(survival$survival)

# see https://stats.stackexchange.com/questions/304132/glmer-not-converging for info on nAGQ argument
glmerSurv <- survival %>% 
  filter(tot > 0) %>%
  glmer(survival ~ Treatment * sprecip7010 * stemp7010  + (1 | siteID / blockID),
        family = "binomial",
        weights = tot,
        nAGQ = 0,
        data = .
  )

simulationOutput <- simulateResiduals(glmerSurv, n = 250)
plot(simulationOutput)
testZeroInflation(simulationOutput)

glmerSurv <- tidy(glmerSurv) %>% 
  mutate(lower = (estimate - std.error*1.96),
         upper = (estimate + std.error*1.96))

glmerSurv <- glmerSurv %>% 
  mutate(term = gsub("Treatment", "", term),
         term = gsub(":", " x ", term),
         term = gsub("sprecip7010", "P", term),
         term = gsub("stemp7010", "t", term),
         term = gsub("FGB", "Gap", term),
         term = gsub("GB", "F", term),
         term = gsub("B", "GF", term),
         term = gsub("G", "FB", term),
         term = gsub("FBap", "Gap", term),
         term = gsub("FBF", "GF", term)
  ) %>% 
  filter(is.na(group)) %>% 
  select(-group) %>% 
  mutate(estimate = round(estimate, 3), std.error = round(std.error, 3), statistic = round(statistic, 3), p.value = round(p.value, 3))


glmerSurv %>% ggplot(aes(x = term, y = estimate, ymin = lower, ymax = upper)) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #geom_vline(xintercept =  c(2.5, 6.5, 10.5), colour = pal1[c(5,2,4,5,2,4)], size = 27, alpha = 0.2) +
  geom_vline(xintercept = c(6.5, 12.5), colour = "olivedrab4", size = 40.5, alpha = 0.2) +
  geom_point(position = position_dodge(width = 0.1), size = 3, fill = "grey") +
  scale_colour_manual("", values = pal1[c(2,4,5,1)]) + 
  scale_shape_manual("", values = c(21, 22, 23, 24, 25,1), limits = c("Graminoids", "Forbs", "Bryophytes", "Intact vegetation", "Climate")) + 
  scale_x_discrete(limits = c("Gap x t", "F x t", "FB x t", "GF x t",
                              "Gap x P", "F x P", "FB x P", "GF x P",
                              "t", "P", "Gap", "F", "FB", "GF")) +
  coord_flip() +
  axis.dimLarge +
  labs(y = "standardised coefficients") +
  theme_classic() +
  axis.dim +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 10))

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig2Coef.jpg", dpi = 300, width = 4, height = 5)



write_csv(survMods, path = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/data/survMods.csv")

#### Figures ####

seedTotPlot %>% 
  filter(Treatment %in% c("Intact", "Gap")) %>% 
  ggplot(aes(x = tempLevel, y = perChange, linetype = Treatment, shape = Treatment, colour = Treatment)) +
  stat_summary(fun.data = "mean_cl_boot", na.rm = TRUE, position = position_dodge(width = 0.3), size = 1) +
  scale_shape_manual(values = c(21,24)) +
  scale_colour_manual(values = c("grey40", "darkgoldenrod3")) +
  geom_hline(yintercept = 0, colour = "lightgrey") +
  labs(title = "", y = "Change in seedling abundance (%)") +
  axis.dimLarge +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())

ggsave(filename = "../paper 4/figures/fig5a.jpg", dpi = 300, width = 7, height = 4)

# change in seedling abundance from early to late growing season
rc_rtcSum %>% 
  ggplot(aes(x = tempLevel, y = perChange, linetype = Treatment, colour = Treatment, shape = Treatment)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.3), size = 1) +
  geom_hline(yintercept = 0, colour = "lightgrey") +
  scale_shape_manual(values = c(21,24)) +
  scale_colour_manual(values = c("grey40", "darkgoldenrod3")) +
  labs(title = "", y = "Change in seedling abundance (%)") +
  axis.dimLarge +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())

ggsave(filename = "../paper 4/figures/fig5b.jpg", dpi = 300, width = 7, height = 4)


# change in seedling abundance from early to late growing season with precipitation
rc_rtcSum %>% 
  ggplot(aes(x = seasonDate, y = seed, linetype = Treatment, colour = Treatment, shape = Treatment, group = Treatment)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.2), size = 1.1) +
  stat_summary(fun.data = "mean_cl_boot", geom = "line", position = position_dodge(width = 0.2), size = 1.1) +
  facet_grid(.~precipLevel) +
  scale_shape_manual(values = c(21,24)) +
  scale_colour_manual(values = c("grey40", "darkgoldenrod3")) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  geom_vline(xintercept = c(as_date(c("25-05-2010","25-05-2012"), tz = "Europe/Oslo", format = "%d-%m-%Y")), alpha = 0.17, colour = "olivedrab4", size = 41) +
  labs(y = "Seedling sum") +
  theme_classic() +
  axis.dimLarge +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig4a.jpg", dpi = 300, width = 16, height = 4)


# change in seedling abundance from early to late growing season with temperature
rc_rtcSum %>% ggplot(aes(x = season, y = seed, linetype = Treatment, colour = Treatment, shape = Treatment, group = Treatment)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.2), size = 1.1) +
  stat_summary(fun.data = "mean_cl_boot", geom = "line", position = position_dodge(width = 0.2), size = 1.1) +
  facet_grid(.~ tempLevel) +
  scale_shape_manual(values = c(21,24)) +
  scale_colour_manual(values = c("grey40", "darkgoldenrod3")) +
  labs(y = "Seedling sum") +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  geom_vline(xintercept = c(as_date(c("25-05-2010","25-05-2012"), tz = "Europe/Oslo", format = "%d-%m-%Y")), alpha = 0.17, colour = "olivedrab4", size = 41) +
  theme_classic() +
  axis.dimLarge +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())
ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig4b.jpg", dpi = 300, width = 13, height = 4.3)

survival <- survival %>% 
  mutate(Treatment = factor(Treatment, labels = c(aC = "Intact", GB = "F", G = "FB", B = "GF", FGB = "Gap")))

# seedling survival following FG removal along temperature gradient
survival %>% filter(tot > 0) %>% 
  ggplot(aes(x = Treatment, y = survival)) + 
  geom_point(aes(size = tot), shape = 21, fill = "goldenrod", alpha = 0.4) + 
  facet_grid(.~tempLevel) +
  theme_classic() +
  axis.dimLarge +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank())

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig2a.jpg", width = 8, height = 3.6, dpi = 300)


# seedling survival following FG removal along precipitation gradient
survival %>% filter(tot > 0) %>% 
  ggplot(aes(x = Treatment, y = survival)) + 
  geom_point(aes(size = tot), shape = 21, fill = "goldenrod", alpha = 0.4) + 
  facet_grid(.~precipLevel) +
  theme_classic() +
  axis.dimLarge +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank())

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig2b.jpg", width = 9, height = 3.5, dpi = 300)

survival %>% 
  ggplot(aes(x = tempLevel, y = survival, lty = Treatment, shape = Treatment, colour = Treatment)) +
  geom_point(position = position_jitter(width = 0.2), aes(size = tot), shape = 21, fill = "grey80", alpha = 0.4) +
  stat_summary(fun.data = "mean_cl_boot", geom = "line", size = 1.3, position = position_dodge(width = 0.25)) +
  stat_summary(fun.data = "mean_cl_boot", size = 1, position = position_dodge(width = 0.25)) +
  scale_colour_manual(values = c("#592520", "#8C3A32", "#9F6E69", "#D95A4E", "#D9958F")) +
  scale_shape_manual(values = c(21,22,23,24,25)) +
  labs(x = "Temperature (ºC)", y = "Survival") +
  theme_classic() +
  axis.dimLarge +
  theme(legend.title = element_blank(),
        legend.key.width=unit(1.7,"cm"))

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig2c.jpg", width = 8, height = 5, dpi = 300)

survival %>% 
  ggplot(aes(x = precipLevel, y = survival, lty = Treatment, shape = Treatment, colour = Treatment)) +
  geom_point(position = position_jitter(width = 80), aes(size = tot), shape = 21, fill = "grey80", alpha = 0.4) +
  stat_summary(fun.data = "mean_cl_boot", geom = "line", size = 1.3, position = position_dodge(width = 0.25)) +
  stat_summary(fun.data = "mean_cl_boot", size = 1, position = position_dodge(width = 80)) +
  scale_colour_manual(values = c("#592520", "#8C3A32", "#9F6E69", "#D95A4E", "#D9958F")) +
  scale_shape_manual(values = c(21,22,23,24,25)) +
  labs(x = "Annual precipitation (mm/yr)", y = "Survival") +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2500)) +
  theme_classic() +
  axis.dimLarge +
  theme(legend.title = element_blank(),
        legend.key.width=unit(1.7,"cm"))

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig2d.jpg", width = 8.6, height = 5, dpi = 300)

lm1 <- rc_rtcSumAv %>% 
  filter(!season == "spr_12") %>% 
  mutate(precip7010 = scale((precip7010/1000), scale = FALSE, center = TRUE),
         temp7010 = scale(temp7010, scale = FALSE, center = TRUE)) %>%
  glmmTMB(seed ~ Treatment*tAnom*precip7010*temp7010 - Treatment:tAnom:precip7010:temp7010 + (1|siteID),
        family = "poisson", 
        zi = ~ monthN + temp7010,
        data = .)

simulationOutput <- simulateResiduals(lm1, n = 250)
plot(simulationOutput)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

summary(lm1)

lm2 <- rc_rtcSumAv %>% 
  filter(!season == "spr_12") %>% 
  mutate(precip7010 = scale((precip7010/1000), scale = FALSE, center = TRUE),
         temp7010 = scale(temp7010, scale = FALSE, center = TRUE)) %>%
  glmer.nb(seed ~ Treatment + (1|siteID), data = .)
  #glmmTMB(seed ~ Treatment*pAnom*precip7010*temp7010 - Treatment:pAnom:precip7010:temp7010 + (1|siteID),
  #        family = "poisson", 
  #        zi = ~ monthN + temp7010,
  #        data = .)

simulationOutput <- simulateResiduals(lm2, n = 250)
plot(simulationOutput)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

summary(lm2)

# filtering out an extreme temperature at -4 which doesn't affect trend and enables clearer illustration in the figure
rc_rtcSumAv %>% 
  mutate(monthN = recode(monthN, aut = "Late", spr = "Early")) %>% 
  filter(tAnom > -4) %>% 
  ggplot(aes(x = tAnom, y = seed, linetype = Treatment, shape = Treatment,colour = precipLevel)) +
  geom_vline(xintercept = 0, colour = "lightgrey") +
  geom_point(position = position_jitterdodge(dodge.width = 0.1), size = 3, stroke = 1, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_shape_manual(values = c(21,24)) +
  scale_colour_manual(values = c("#3DBFBF","#339F9F", "#038C8C", "#025959")) +
  facet_grid(tempLevel~monthN) +
  labs(x = "Temperature deviation from 2009-2018 mean (Δ ºC)",
       y = "Seed number") +
  ylim(c(0,80)) +
  theme_classic() +
  axis.dimLarge +
  theme(legend.title = element_blank())

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig6a.jpg", width = 10, height = 7, dpi = 300)

rc_rtcSumAv %>% 
  mutate(monthN = recode(monthN, aut = "Late", spr = "Early")) %>% 
  ggplot(aes(x = pAnom, y = seed, linetype = Treatment, shape = Treatment, colour = precipLevel)) +
  geom_vline(xintercept = 0, colour = "lightgrey") +
  geom_point(position = position_jitterdodge(dodge.width = 0.1), size = 3, stroke = 1, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_shape_manual(values = c(21,24)) +
  scale_colour_manual(values = c("#3DBFBF","#339F9F", "#038C8C", "#025959")) +
  facet_grid(. ~ monthN) +
  labs(x = "Soil moisture deviation from 2009-2018 mean",
       y = "Seed number") +
  ylim(c(0,80)) +
  theme_classic() +
  axis.dimLarge +
  theme(legend.title = element_blank())

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig6b.jpg", width = 9, height = 4, dpi = 300)


lm3 <- seedTot %>% 
  mutate(monthN = recode(monthN, spr = "Early", aut = "Late")) %>% 
  glmmTMB(seed ~ Treatment*tAnom*sprecip7010*stemp7010 - Treatment:tAnom:sprecip7010:stemp7010 + (1|siteID/monthN),
          family = "poisson", 
          zi = ~ 1,
          data = .)

simulationOutput <- simulateResiduals(lm3, n = 250)
plot(simulationOutput)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

summary(lm3)

# seedling abundance following FG removal along temperature gradient
seedTot %>% ggplot(aes(x = Treatment, y = seed, fill = season)) +
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(values = c("goldenrod", "grey")) +
  facet_grid(~tempLevel) +
  labs(y = "Seedling abundance") +
  theme_classic() +
  axis.dimLarge +
  theme(axis.title.x = element_blank())

seedTot %>% 
  mutate(monthN = recode(monthN, spr = "Early", aut = "Late")) %>% 
  ggplot(aes(x = tAnom, y = seed, linetype = Treatment, shape = Treatment, colour = tempLevelPlot)) +
  geom_vline(xintercept = 0, colour = "lightgrey") +
  geom_point(position = position_jitterdodge(dodge.width = 0.1), size = 3, stroke = 1, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_shape_manual(values = c(21,22,23,24,25)) +
  scale_colour_manual(values = c("#D9958F", "#8C3A32", "#592520", "#9F6E69", "#D95A4E")) +
  facet_grid(. ~ monthN) +
  labs(x = "Soil temperature deviation from 2009-2018 mean (ºC)",
       y = "Seed number") +
  ylim(c(0,65)) +
  labs(y = "Seedling abundance") +
  theme_classic() +
  axis.dimLarge +
  theme(legend.title = element_blank(),
        legend.key.width=unit(1.7,"cm")) + 
  guides(linetype = guide_legend(override.aes = list(color = "black")))

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig3a.jpg", width = 8, height = 3.6, dpi = 300)


# seedling abundance following FG removal along precipitation gradient
seedTotPlot %>% ggplot(aes(x = Treatment, y = seed, fill = season)) +
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(values = c("goldenrod", "grey")) +
  facet_grid(~precipLevel) +
  labs(y = "Seedling abundance") +
  theme_classic() +
  axis.dimLarge +
  theme(axis.title.x = element_blank())



seedTot %>% 
  mutate(monthN = recode(monthN, spr = "Early", aut = "Late")) %>% 
  ggplot(aes(x = tAnom, y = seed, linetype = Treatment, shape = Treatment, colour = precipLevelPlot)) +
  geom_vline(xintercept = 0, colour = "lightgrey") +
  geom_point(position = position_jitterdodge(dodge.width = 0.1), size = 3, stroke = 1, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_shape_manual(values = c(21,22,23,24,25)) +
  scale_colour_manual(values = c("#3DBFBF","#339F9F", "#038C8C", "#025959")) +
  facet_grid(. ~ monthN) +
  labs(x = "Soil temperature deviation from 2009-2018 mean (ºC)",
       y = "Seed number") +
  ylim(c(0,80)) +
  theme_classic() +
  axis.dimLarge +
  theme(legend.title = element_blank(),
        legend.key.width=unit(1.7,"cm")) + 
  guides(linetype = guide_legend(override.aes = list(color = "black")))


ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig3b.jpg", width = 8.1, height = 3.5, dpi = 300)


# climate
soilMplot <- soilM %>% 
  filter(dplyr::between(date, ymd("2018-06-01"), ymd("2018-10-01")),
         tempLevel == 8.5) %>%     
  ggplot(aes(x = date, y = value, linetype = factor(precipLevel))) +
  stat_summary(data = LTprecip, fun.data = "mean_cl_boot", aes(y = value), linetype = "solid", colour = "grey", geom = "line") +
  stat_summary(fun.data = "mean_cl_boot", geom = "line") +
  labs(y = "Soil moisture", linetype = "Mean annual\nprecipitation (mm)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


tempA <- mSubClim %>% 
  filter(year == 2018, logger %in% c("temp1", "temp2", "temp200cm", "temperature", "temperature2"))

tempAplot <- tempA %>% 
  filter(dplyr::between(date, ymd("2018-06-01"), ymd("2018-10-01")), tempLevel == 8.5) %>%     
  ggplot(aes(x = date, y = value, linetype = factor(precipLevel))) +
  stat_summary(data = LTtempPlot, fun.data = "mean_cl_boot", aes(y = value), linetype = "solid", colour = "grey", geom = "line") +
  stat_summary(fun.data = "mean_cl_boot", geom = "line") +
  labs(y = "Air temperature (ºC)", linetype = "Mean annual\nprecipitation (mm)") +
  theme(axis.title.x = element_blank())

precipLeg <- get_legend(soilMplot)

figs <- plot_grid(soilMplot + theme(legend.position = "none"), tempAplot + theme(legend.position = "none"), ncol = 1, align = "hv")
figs <- plot_grid(figs, precipLeg, ncol = 2, rel_widths = c(0.88,0.22))

ggsave(figs, filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig2.jpg", height = 6, width = 8, dpi = 300)


LTAnClimPlot <- mSubClimPlot %>% 
  ggplot(aes(x = date, y = tAnom, fill = sign)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = as_date("25-06-2018", tz = "Europe/Oslo", format = "%d-%m-%Y"), alpha = 0.17, colour = "olivedrab4", size = 5) +
  scale_fill_manual(values = c("darkcyan", "coral4")) +
  geom_col() +
  labs(y = "Air temperature\n(Δ 2009-2018 mean, ºC)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  axis.dim


LTClimPlot <- mSubClimPlot %>% 
  ggplot(aes(x = date, y = value)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = as_date("25-06-2018", tz = "Europe/Oslo", format = "%d-%m-%Y"), alpha = 0.17, colour = "olivedrab4", size = 5) +
  stat_summary(fun.data = "mean_cl_boot", geom = "line") +
  labs(y = "Air temperature (ºC)") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  axis.dim

gridPLotLTclim <- plot_grid(LTClimPlot, LTAnClimPlot, ncol = 1, align = "hv")

ggsave(gridPLotLTclim, filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig1.jpg", height = 6, width = 9, dpi = 300)

