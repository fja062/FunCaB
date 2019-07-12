# see https://cms.met.no/site/2/klimaservicesenteret/klima-i-norge-2100/_attachment/11592?_ts=15c10419731 page 17 for predictions

#### maximum temperatures ####
treat <- maxminANALYSIS %>% 
  ungroup() %>% 
  mutate(Year = year(date),
         Treatment = recode(Treatment, "FGB" = "aFGB"),
         sTemp = scale(Temp),
         ssunniness = scale(sunniness)) %>% 
  distinct(siteID, turfID, Treatment, Temp, maxTemp, minTemp, sunniness, Precip, Block, Temperature_level) %>%
  filter(Treatment %in% c("aFGB", "C", "FB", "GB", "GF")) %>% 
  na.omit()

modTreat <- lmer(maxTemp ~ Treatment*scale(sunniness) + 
                   Treatment*scale(Temp) + 
                   Treatment*scale(Precip) +
                   (1|siteID/Block), REML = FALSE, data = treat)

P85 <- treat %>% mutate(Precip = Precip*1.17)
t85 <- treat %>% mutate(Temp = Temp + 3.9)
modSun <- treat %>% mutate(sunniness = sunniness*0.9)
tP85 <- treat %>% mutate(Temp = Temp + 3.9,
                         Precip = Precip*1.17)
tPSun85 <- treat %>% mutate(Temp = Temp + 3.9,
                            Precip = Precip*1.17,
                            sunniness = sunniness*0.9)

treat$modPreds <- predict(modTreat, re.form = ~1|siteID/Block)
treat$P85 <-  predict(modTreat, newdata = P85, re.form = ~1|siteID/Block)
treat$t85 <-  predict(modTreat, newdata = t85, re.form = ~1|siteID/Block)
treat$tP85 <-  predict(modTreat, newdata = tP85, re.form = ~1|siteID/Block)
treat$modSun <-  predict(modTreat, newdata = modSun, re.form = ~1|siteID/Block)
treat$tPSun85 <-  predict(modTreat, newdata = tPSun85, re.form = ~1|siteID/Block)

# single remainers, plus sunnines*temp
treat %>% gather(maxTemp, modPreds, modSun, P85, t85, tP85, tPSun85, key = Model, value = value) %>% 
  #filter(Temperature_level == 6.5) %>% 
  ggplot(aes(x = Treatment, y = value)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  #geom_vline(xintercept = c(3.5, 7.5), linetype = "dashed", colour = "grey60") +
  facet_grid(Temperature_level~Model)
#stat_summary(fun.y = mean, geom = "boxplot", size = 1)

pal1 <- wes_palette(7, name = "Darjeeling2", type = "continuous")

maxTempPredictions <- treat %>% gather(maxTemp, modPreds, modSun, P85, t85, tP85, tPSun85, key = Model, value = value) %>% 
  ggplot(aes(x = Temperature_level, y = value, colour = Treatment)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.6)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.6), geom = "line") +
  facet_grid(.~Model) +
  scale_color_manual(values = pal1)
ggsave(maxTempPredictions, filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 2/figures/maxTempPredictions.jpg", dpi = 300, width = 13, height = 6)



#### minimum temperature ####
treatmin <- treat
modTreatmin <- lmer(minTemp ~ Treatment*scale(sunniness) + 
                      Treatment*scale(Temp) + 
                      Treatment*scale(Precip) +
                      (1|siteID/Block), REML = FALSE, data = treatmin)

P85 <- treatmin %>% mutate(Precip = Precip*1.17)
t85 <- treatmin %>% mutate(Temp = Temp + 3.9)
modSun <- treatmin %>% mutate(sunniness = sunniness*0.9)
tP85 <- treatmin %>% mutate(Temp = Temp + 3.9,
                            Precip = Precip*1.17)
tPSun85 <- treatmin %>% mutate(Temp = Temp + 3.9,
                               Precip = Precip*1.17,
                               sunniness = sunniness*0.9)

treatmin$modPreds <- predict(modTreatmin, re.form = ~1|siteID/Block)
treatmin$P85 <-  predict(modTreatmin, newdata = P85, re.form = ~1|siteID/Block)
treatmin$t85 <-  predict(modTreatmin, newdata = t85, re.form = ~1|siteID/Block)
treatmin$tP85 <-  predict(modTreatmin, newdata = tP85, re.form = ~1|siteID/Block)
treatmin$modSun <-  predict(modTreatmin, newdata = modSun, re.form = ~1|siteID/Block)
treatmin$tPSun85 <-  predict(modTreatmin, newdata = tPSun85, re.form = ~1|siteID/Block)

pal1 <- wes_palette(7, name = "Darjeeling2", type = "continuous")

minTempPredictions <- treatmin %>% gather(minTemp, modPreds, modSun, P85, t85, tP85, tPSun85, key = Model, value = value) %>% 
  ggplot(aes(x = Temperature_level, y = value, colour = Treatment)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.6)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.6), geom = "line") +
  facet_grid(.~Model) +
  scale_color_manual(values = pal1)

ggsave(minTempPredictions, filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 2/figures/minTempPredictions.jpg", dpi = 300, width = 13, height = 6)


#### temperature amplitude ####
treat <- maxminANALYSIS %>% 
  ungroup() %>% 
  mutate(tempAmp = maxTemp - minTemp,
         Year = year(date),
         Treatment = recode(Treatment, "FGB" = "aFGB"),
         sTemp = scale(Temp),
         ssunniness = scale(sunniness)) %>% 
  distinct(siteID, turfID, Treatment, tempAmp, Temp, maxTemp, minTemp, sunniness, Precip, Block, Temperature_level) %>%
  filter(Treatment %in% c("aFGB", "C", "FB", "GB", "GF")) %>% 
  na.omit()

modTreat <- lmer(tempAmp ~ Treatment*scale(sunniness) + 
                   Treatment*scale(Temp) + 
                   Treatment*scale(Precip) +
                   (1|siteID/Block), REML = FALSE, data = treat)

P85 <- treat %>% mutate(Precip = Precip*1.17)
t85 <- treat %>% mutate(Temp = Temp + 3.9)
modSun <- treat %>% mutate(sunniness = sunniness*0.9)
tP85 <- treat %>% mutate(Temp = Temp + 3.9,
                         Precip = Precip*1.17)
tPSun85 <- treat %>% mutate(Temp = Temp + 3.9,
                            Precip = Precip*1.17,
                            sunniness = sunniness*0.9)

treat$modPreds <- predict(modTreat, re.form = ~1|siteID/Block)
treat$P85 <-  predict(modTreat, newdata = P85, re.form = ~1|siteID/Block)
treat$t85 <-  predict(modTreat, newdata = t85, re.form = ~1|siteID/Block)
treat$tP85 <-  predict(modTreat, newdata = tP85, re.form = ~1|siteID/Block)
treat$modSun <-  predict(modTreat, newdata = modSun, re.form = ~1|siteID/Block)
treat$tPSun85 <-  predict(modTreat, newdata = tPSun85, re.form = ~1|siteID/Block)

# single remainers, plus sunnines*temp
treat %>% gather(tempAmp, modPreds, modSun, P85, t85, tP85, tPSun85, key = Model, value = value) %>% 
  #filter(Temperature_level == 6.5) %>% 
  ggplot(aes(x = Treatment, y = value)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  #geom_vline(xintercept = c(3.5, 7.5), linetype = "dashed", colour = "grey60") +
  facet_grid(Temperature_level~Model)
#stat_summary(fun.y = mean, geom = "boxplot", size = 1)

pal1 <- wes_palette(7, name = "Darjeeling2", type = "continuous")

tempAmpPredictions <- treat %>% gather(tempAmp, modPreds, modSun, P85, t85, tP85, tPSun85, key = Model, value = value) %>% 
  ggplot(aes(x = Temperature_level, y = value, colour = Treatment)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.6)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.6), geom = "line") +
  facet_grid(.~Model) +
  scale_color_manual(values = pal1)
ggsave(tempAmpPredictions, filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 2/figures/maxTempPredictions.jpg", dpi = 300, width = 13, height = 6)
