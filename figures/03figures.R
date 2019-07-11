source("~/OneDrive - University of Bergen/Research/FunCaB/SeedclimComm/inst/graminoidRemovals/plotting_dim.R")

# load libraries
library(cowplot)
library(wesanderson)

soilTempPlot <- left_join(soilTempPlot, weather)

############################ 
#### Figure for article ####

## Figure 2 
intermediate <- soilTempPlot %>%
  mutate(weather = case_when(sunniness > 0.66 ~ "sunny",
                             sunniness > 0.33 ~ "variable",
                             sunniness < 0.33 ~ "cloudy")) %>% 
  filter( weather %in% c("sunny", "cloudy"), between(date, ymd("2015-08-01"), ymd("2015-09-30")), tempLevel == 8.5) %>% 
  ggplot(aes(x = hour, y = Value, colour = Treatment, linetype = weather)) +
  geom_smooth(se = FALSE, method = "loess") +
  scale_colour_manual("Functional groups", values = c("black", cbPalette[c(7,5,6,2,4,3,1)]), limits = c("FGB","GF", "GB", "FB", "G", "F", "B", "C"), labels = c("Bare ground","Bryophytes", "Forbs","Graminoids", "Bryophytes and forbs", "Bryophytes and graminoids", "Forbs and graminoids", "Bryophytes, forbs and\ngraminoids")) +
  labs(x = "Time (hr)", y = "Soil temperature (ºC)") +
  #scale_linetype_manual(values = c("73",1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 19)) +
  #ylim(7,19) +
  axis.dimLarge

  
x <- FD %>% 
  left_join(weather) %>% 
  filter(between(date, left = dmy("01-07-2015"), right = dmy("15-06-2016")), tempLevel == 8.5) %>% 
  ggplot(aes(x = date, y = sum, colour = Treatment)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "line", size = 0.8, na.rm = TRUE) +
  scale_colour_manual("Vegetation", values = c("black", cbPalette[c(7,5,6,2,4,3,1)]), limits = c("FGB", "GF", "GB", "FB", "G", "F", "B", "C")) +
  scale_fill_manual("Vegetation", values = c("black", cbPalette[c(7,5,6,2,4,3,1)]), limits = c("FGB", "GF", "GB", "FB", "G", "F", "B", "C")) +
  scale_x_date(date_labels = "%b %Y") +
  coord_cartesian(ylim = c(0, 31),
                  xlim = c(dmy("01-07-2015"), dmy("30-07-2016"))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  axis.dimLarge +
  ylab("Cumulative frost\ndays") +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.title = element_text(size=15))

y <- soilTempPlot %>%
  filter(between(date, ymd("2015-07-01"), ymd("2016-06-15")), tempLevel == 8.5) %>% 
  #mutate(Treatment = factor(Treatment, levels = c("C", "FB", "GB", "GF", "FGB"))) %>% 
  ggplot(aes(x = date, y = Value, colour = Treatment)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "line", size = 0.8, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual("Functional groups", values = c("black", cbPalette[c(7,5,6,2,4,3,1)]), limits = c("FGB","GF", "GB", "FB", "G", "F", "B", "C"), labels = c("Bare ground","Bryophytes", "Forbs","Graminoids", "Bryophytes & forbs", "Graminoids & bryophytes", "Forbs & graminoids", "Bryophytes, forbs &\ngraminoids")) +
  ylab("Soil temperature (°C)") +
  axis.dimLarge +
  theme(axis.title = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

legend <- get_legend(intermediate)
z <- plot_grid(y + theme(legend.position = "none"), x, ncol = 1, rel_heights = c(1.5, 1), align = "hv", labels = c("B", "C"), label_x = 0.055, label_y = 0.97)
#eps <- plot_grid(legend, lav, ncol = 1, rel_heights = c(1,1.65))

frostSumAnn <- plot_grid(intermediate + theme(legend.position = "none"), z, legend, nrow = 1, rel_widths = c(0.92, 0.9, 0.50), labels = "A", label_x = 0.02, label_y =0.97)
ggsave(frostSumAnn, filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 2/figures/fig2.jpg", dpi = 300, width = 13, height = 5)


# Figure 3A
covAnomPlot <- Cover %>% 
  filter(between(date, ymd("2015-08-01"), ymd("2015-09-30")),
         Treatment %in% c("GF", "FB", "GB"),
         !(blockID == 1 & siteID =="Ulvhaugen")) %>% 
  mutate(temp = if_else(grepl("6.5", tempLevel), "alpine", if_else(grepl("8.5", tempLevel), "sub-alpine", "boreal"))) %>% 
  mutate(temp = factor(temp, levels = c("alpine", "sub-alpine", "boreal"))) %>% 
  gather(key = response, value = value, graminoidCov, forbCov, mossCov) %>% 
  filter(value > 0,
         weather %in% c("sunny")) %>% 
  mutate(response = factor(response, levels = c("graminoidCov", "forbCov", "mossCov"))) %>% 
  filter(response %in% c("graminoidCov", "forbCov", "mossCov")) %>% 
  mutate(response = recode(response, graminoidCov = "Graminoid", forbCov = "Forb", mossCov = "Bryophyte")) %>% 
  ggplot(aes(x = value, y = maxAnom, fill = response)) +
  geom_point(alpha = 0.5, shape = 21) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", size = 0.5) +
  scale_fill_manual("", values = pal1[c(2,4,5)]) +
  facet_grid(. ~ temp) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Cover (%)",
       y = "temperature anomaly from bare ground") +
  theme_classic() +
  axis.dimLarge

ggsave(covAnomPlot, file = "~/OneDrive - University of Bergen/Research/FunCaB/paper 2/figures/Fig3A.jpg", dpi = 300, width = 9, height = 4)

# Figure 3B
heiAnomPlot <- maxminAnom %>% 
  filter(between(date, ymd("2015-08-01"), ymd("2015-09-30")),
         Treatment %in% c("GF", "FB", "GB")) %>% 
  mutate(temp = if_else(grepl("6.5", tempLevel), "alpine", if_else(grepl("8.5", tempLevel), "sub-alpine", "boreal"))) %>% 
  mutate(temp = factor(temp, levels = c("alpine", "sub-alpine", "boreal"))) %>% 
  gather(key = response, value = value, forbHeight, graminoidHeight, mossHeight) %>% 
  filter(between(value, 1, 100),
         weather %in% c("cloudy", "sunny")) %>% 
  mutate(response = factor(response, levels = c("forbHeight", "graminoidHeight", "mossHeight"))) %>% 
  mutate(response = recode(response, vegetationHeight = "Vascular plant", mossHeight = "Bryophyte")) %>% 
  ggplot(aes(x = value, y = maxAnom, fill = response)) +
  geom_point(alpha = 0.5, shape = 21) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", size = 0.5) +
  #geom_point() +
  scale_fill_manual(values = pal1[c(1,3,5)]) +
  #scale_fill_manual(values = c("black", "grey60")) +
  facet_grid(. ~ temp) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Biomass-weighted vegetation height",
       y = "temperature anomaly from bare ground") +
  theme_classic() +
  axis.dimLarge


ggsave(heiAnomPlot, file = "~/OneDrive - University of Bergen/Research/FunCaB/paper 2/figures/Fig3B.jpg", dpi = 300, width = 8, height = 4)


### supplementary figures ###
# Graminoid, forb and bryophyte cover distributions
g <- Cover %>% filter(Treatment == "FB") %>% 
  ggplot(aes(x = graminoidCov)) + geom_density(alpha = 0.3, fill = "goldenrod") + 
  geom_rug() + 
  geom_vline(xintercept = 10, linetype = "dashed") + 
  facet_grid(.~tempLevel) +
  xlim(0,100) +
  labs(x = "graminoid cover (%)")
f <- Cover %>% filter(Treatment == "GB") %>% 
  ggplot(aes(x = forbCov)) + geom_density(alpha = 0.3, fill = "goldenrod") + 
  geom_rug() + 
  geom_vline(xintercept = 10, linetype = "dashed") + 
  facet_grid(.~tempLevel) +
  xlim(0,100) +
  labs(x = "forb cover (%)")
m <- Cover %>% filter(Treatment == "GF") %>% 
  ggplot(aes(x = mossCov)) + geom_density(alpha = 0.3, fill = "goldenrod") + 
  geom_rug() + 
  geom_vline(xintercept = 10, linetype = "dashed") + 
  facet_grid(.~tempLevel) +
  xlim(0,100) +
  labs(x = "bryophyte cover (%)")

fgmDens <- plot_grid(g + theme(axis.title.y = element_blank()),f,m + theme(axis.title.y = element_blank()), ncol = 1,align = "hv")
ggsave(fgmDens, file = "~/OneDrive - University of Bergen/Research/FunCaB/paper 2/figures/supFig9.jpg", dpi = 300, height = 7, width = 7)

# temperature range anomalies
magAmpCovPlot <- magAmpAnom %>% 
  filter(between(date, ymd("2015-07-01"), ymd("2015-09-30")),
         Treatment %in% c("GF", "FB", "GB")) %>% 
  gather(key = response, value = value, graminoidCov, vegetationHeight, forbCov, mossHeight, mossCov) %>% 
  filter(value > 0,
         weather %in% c("cloudy", "sunny")) %>% 
  mutate(response = factor(response, levels = c("graminoidCov", "forbCov", "mossCov", "vegetationHeight", "mossHeight"))) %>% 
  filter(response %in% c("graminoidCov", "forbCov", "mossCov")) %>% 
  ggplot(aes(x = value, y = magAnom, fill = response)) +
  #geom_point(alpha = 0.9, shape = 21) +
  stat_summary(geom = "point", fun.y = "mean", alpha = 0.9, shape = 21, size = 2) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", size = 0.5) +
  #geom_point() +
  scale_fill_manual(values = pal1[c(2,4,5,1)]) +
  #scale_fill_manual(values = c("black", "grey60")) +
  #facet_grid(tempLevel ~ weather, scales = "free_x") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Cover (%)",
       y = "temperature range anomaly from bare ground")
ggsave(magAmpCovPlot, file = "~/OneDrive - University of Bergen/Research/FunCaB/paper 2/figures/supFig8.jpg", dpi = 300, width = 9, height = 4.5)

magAmpHeightPlot <- magAmpAnom %>% 
  filter(between(date, ymd("2015-07-01"), ymd("2015-09-30")),
         Treatment %in% c("GF", "FB", "GB")) %>% 
  gather(key = response, value = value, graminoidCov, vegetationHeight, forbCov, mossHeight, mossCov) %>% 
  filter(value > 0,
         weather %in% c("cloudy", "sunny")) %>% 
  mutate(response = factor(response, levels = c("graminoidCov", "forbCov", "mossCov", "vegetationHeight", "mossHeight"))) %>% 
  filter(response %in% c("vegetationHeight", "mossHeight")) %>% 
  ggplot(aes(x = value, y = magAnom, fill = response)) +
  #geom_point(alpha = 0.9, shape = 21) +
  stat_summary(geom = "point", fun.y = "mean", alpha = 0.9, shape = 21, size = 2) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", size = 0.5) +
  #geom_point() +
  scale_fill_manual(values = pal1[c(2,4,5)]) +
  #scale_fill_manual(values = c("black", "grey60")) +
  facet_grid(. ~ weather, scales = "free_x") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Height",
       y = "temperature range anomaly from bare ground")

ggsave(magAmpHeightPlot, file = "~/OneDrive - University of Bergen/Research/FunCaB/paper 2/figures/supFig7.jpg", dpi = 300, width = 9, height = 4.5)


