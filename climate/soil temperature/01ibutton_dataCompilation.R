# load libraries
library(lubridate)

#source biomass and composition
source("~/Documents/FunCaB/bioticInteractions/composition/dataProcessing/00funcab_data_processing.R")

source("~/Documents/FunCaB/bioticInteractions/biomass/dataProcessing/biomass_cleaning.R")

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
# use if using biomass
#vegComp <- biomassReg %>% 
#  ungroup()

# use this if only using cover
vegCov <- comp2 %>% 
  filter(!Treatment == "XC", Year == 2015) %>% 
  distinct(siteID, blockID, Treatment, turfID, graminoidCov, forbCov, mossCov, functionalGroup, sumcover)
  
vegComp <- soilTemp  %>%
  left_join(vegCov) %>%
  group_by(turfID) %>% 
  distinct(siteID, turfID, blockID, Treatment, date, sunniness, TOD, meanTemp, .keep_all = TRUE) %>% 
  filter(TOD == "day")

# turn covers to NA where FG has been removed
#vegComp <- vegComp %>% 
#  mutate(forbCov = case_when(
#      grepl("F", Treatment) ~ NA_real_,
#      TRUE ~ forbCov),
#    graminoidCov = case_when(
#      grepl("G", Treatment) ~ NA_real_,
#      TRUE ~ graminoidCov),
#    mossCov = case_when(
#      grepl("B", Treatment) ~ NA_real_,
#      TRUE ~ mossCov))

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
         sforbCov = scale(forbCov, scale = FALSE, center = TRUE),
         smossCov = scale(mossCov, scale = FALSE, center = TRUE),
         sgraminoidCov = scale(graminoidCov, scale = FALSE, center = TRUE))


CoverSingles <- Cover %>% 
  group_by(siteID, blockID, date) %>% 
  filter(Treatment %in% c("aFGB", "FB", "GF", "GB")) %>% 
  spread(key = Treatment, value = maxTemp) %>%
  pivot_longer(c("FB", "GB", "GF"), names_to = "Treatment", values_to = "maxTemp") %>% 
  mutate(maxTempDelta = maxTemp - aFGB)

CoverDoubles <- Cover %>% 
  filter(Treatment %in% c("FB", "GF", "GB", "F", "G", "B")) %>% 
  spread(key = Treatment, value = maxTemp) %>%
  mutate(maxTempDelta = maxTemp - aFGB)

save(Cover, file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/cleanedSoilTemp.RData")
