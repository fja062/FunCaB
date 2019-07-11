# load libraries
library(lubridate)

#source biomass and composition
source("~/Documents/SeedClim-Climate-Data/funcab/vegetation/00funcab_data_processing.R")
#source("~/OneDrive - University of Bergen/Research/FunCaB/SeedClim-Climate-Data/funcab/vegetation/00funcab_data_processing.R")

source("~/Documents/SeedClim-Climate-Data/funcab/vegetation/biomass_cleaning.R")
#source("~/OneDrive - University of Bergen/Research/FunCaB/SeedClim-Climate-Data/funcab/vegetation/biomass_cleaning.R")

# load soil temperature data
load("~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/soilTemp.RData")

#extract air temperature
airTemp <- soilTemp %>% 
  filter(!is.na(airTemp)) %>% 
  filter(between(date, ymd("2015-07-01"), ymd("2015-10-30"))) %>% 
  group_by(siteID, date) %>% 
  summarise(maxTemp = max(airTemp),
            meanTemp = mean(airTemp))

# use air temperature to create cut-off for summer temperature analyses, ie. below daily average of 5deg.
ggplot(airTemp, aes(x = date, y = meanTemp, colour = siteID)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "line") +
  geom_hline(yintercept = 5)

airTemp %>% filter(between(date, ymd("2015-07-01"), ymd("2015-11-30"))) %>%
  group_by(siteID) %>% 
  filter(maxTemp <= 5) %>% 
  filter(date == min(date)) %>% 
  ggplot(aes(date, maxTemp, label = siteID)) +
  geom_text(position = position_jitter(width = 1, height = 1))
# cut-off date: 20 Oct 2015.

#summarise to mean daily temp
soilTempPlot <- soilTemp %>%
  ungroup() %>% 
  filter(between(date, ymd("2015-07-01"), ymd("2016-08-30")), !Treatment == "temp200cm", !TOD == "spinup") %>% 
  mutate(turfID = recode(turfID, "GudNAGB" = "Gud13GB"),
         Block = if_else(turfID == "Gud13GB", "13", Block))

soilTemp <- soilTempPlot %>% 
  group_by(siteID, turfID, blockID = Block, Treatment, date, sunniness, TOD, ID) %>%
  summarise(meanTemp = mean(Value),
            maxTemp = max(Value),
            minTemp = min(Value),
            magTemp = (maxTemp - minTemp)) %>% 
  ungroup()

biomassReg <- biomassReg %>% 
  ungroup()
  
vegComp <- soilTemp  %>%
  left_join(biomassReg) %>%
  group_by(turfID) %>% 
  distinct(siteID, turfID, blockID, Treatment, date, sunniness, TOD, meanTemp, .keep_all = TRUE) %>% 
  filter(TOD == "day")

# turn covers to zero where FG has been removed
vegComp <- vegComp %>% 
  mutate(forbBiomass = if_else(Treatment %in% c("F", "FB", "GF", "FGB"), 0, forbBiomass),
         graminoidBiomass = if_else(Treatment %in% c("G", "GF", "GB", "FGB"), 0, graminoidBiomass),
         mossBiomass = if_else(Treatment %in% c("B", "GB", "FB", "FGB"), 0, mossBiomass),
         forbCov = if_else(Treatment %in% c("F", "FB", "GF", "FGB"), 0, forbCov),
         graminoidCov = if_else(Treatment %in% c("G", "GF", "GB", "FGB"), 0, graminoidCov),
         mossCov = if_else(Treatment %in% c("B", "GB", "FB", "FGB"), 0, mossCov),
         mossHeight = if_else(Treatment %in% c("B", "GB", "FB", "FGB"), 0, mossHeight),
         vegetationHeight = if_else(Treatment %in% c("GF", "FGB"), 0, vegetationHeight))

# categorise weather and filter for summer months
vegComp <- vegComp %>% 
  mutate(weather = case_when(
    sunniness > 0.66 ~ "sunny",
    sunniness > 0.33 ~ "variable",
    sunniness < 0.33 ~ "cloudy")) %>% 
  ungroup()

save(vegComp, file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/cleanedVegComp.RData")


#--- soil freezing ---#
FD <- vegComp %>% 
  filter(!turfID == "Ves3G",
         between(date, left = dmy("01-08-2015"), right = dmy("30-06-2016")),
         TOD == "day") %>% 
  mutate(x = minTemp < 0) %>%
  arrange(date) %>% 
  group_by(turfID) %>% 
  mutate(sum = cumsum(x)) %>%
  ungroup()

save(FD, file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/cleanedFrostDays.RData")


# --- temperature anomalies ---#
maxmin <- vegComp %>% 
  mutate(month = month(date)) %>%
  filter(month %in% c(6, 7, 8, 9)) %>% 
  left_join(weather)

maxminAnom <- maxmin %>%
  filter(TOD == "day") %>% 
  left_join(maxmin %>% filter(Treatment == "FGB") %>% ungroup() %>% select(FGBmaxTemp = maxTemp, date, siteID, blockID)) %>%
  mutate(maxAnom = maxTemp - FGBmaxTemp) %>% 
  ungroup()

magAmpAnom <- maxmin %>%
  group_by(turfID, date) %>% 
  mutate(magTemp = maxTemp - minTemp) %>% 
  left_join(maxmin %>% filter(Treatment == "FGB") %>% ungroup() %>% select(FGBmagTemp = magTemp, date, siteID, blockID)) %>%
  mutate(magAnom = magTemp - FGBmagTemp) %>% 
  ungroup()


# filter for 1st July - 15th september for analyses, gather temperature calculationg, and recode FGB to be intercept
maxmin <- maxmin %>% 
  filter(between(date, left = dmy("01-07-2015"), right = dmy("20-09-2015"))) %>%
  mutate(Year = year(date),
         Treatment = recode(Treatment, "FGB" = "aFGB"))

# standardise variables
Cover <- maxmin %>% 
  ungroup() %>% 
  mutate(Year = year(date),
         sTemp70 = scale(temp7010, scale = FALSE, center = TRUE),
         sPrecip70 = scale((precip7010/1000), scale = FALSE, center = TRUE),
         ssunniness = scale(sunniness),
         sforbBiomass = scale(forbBiomass, scale = FALSE, center = TRUE),
         smossBiomass = scale(mossBiomass, scale = FALSE, center = TRUE),
         sgraminoidBiomass = scale(graminoidBiomass, scale = FALSE, center = TRUE),
         sforbCov = scale(forbCov, scale = FALSE, center = TRUE),
         smossCov = scale(mossCov, scale = FALSE, center = TRUE),
         sgraminoidCov = scale(graminoidCov, scale = FALSE, center = TRUE),
         svegetationHeight = scale(vegetationHeight, scale = FALSE, center = TRUE),
         smossHeight = scale(mossHeight, scale = FALSE, center = TRUE))


save(Cover, file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/cleanedSoilTemp.RData")
