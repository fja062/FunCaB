# load packages
require(lubridate)

# source dictionaries and data
source("~/Documents/FunCaB/dictionaries.R")

load("~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/climate.Rdata")
load("~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/soilMoisture2018.RData")

subClim <- climate %>% 
  ungroup() %>% 
  mutate(year = year(date), date = date(date), hour = hour(date), day = lubridate::yday(date)) %>%
  filter(year > 2008, 
         logger %in% c("jordf1", "jordf2", "jordfukt2", "nedbor", "sm300 1", "sm300 2", "soil moisture 1", "soil moisture 2", "soil temp", "temp1", "temp2", "temp200cm", "temp30cm", "veg temp", "veg. temp", "rain", "temperature", "temperature2"))

subClim <- subClim %>% 
  left_join(dict_Site, by = c("site" = "v3")) %>% 
  rename(siteID = new)

# long-term daily temperature average
LTtemp <- subClim %>% 
  filter(logger %in% c("temp200cm")) %>% 
  group_by(day, siteID) %>% 
  summarise(value = mean(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(day = day - 1, 
         date = as.Date(day, origin = "2018-01-01")) %>% 
  left_join(weather)

LTtempPlot <- LTtemp %>% 
  filter(tempLevel == 8.5, between(day, 152, 274))

# long-term daily precipitation average
LTprecip <- subClim %>% 
  filter(logger %in% c("jordf1", "jordf2", "jordfukt2", "soil moisture 1", "soil moisture 2", "sm300 2", "sm300 1")) %>% 
  group_by(day, siteID) %>% 
  summarise(value = mean(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(day = day - 1, 
         date = as.Date(day, origin = "2018-01-01")) %>% 
  left_join(weather)

LTprecipPlot <- LTprecip %>% 
  filter(tempLevel == 8.5, between(day, 152, 274))


# 2018 daily precipitation and temperature
mSubClim <- subClim %>% 
  group_by(date, day, logger, siteID, year) %>% 
  summarise(value = mean(value, na.rm = TRUE)) %>% 
  left_join(weather) %>% 
  ungroup()

soilM <- mSubClim %>% 
  filter(year == 2018, logger %in% c("jordf1", "jordf2", "jordfukt2", "soil moisture 1", "soil moisture 2", "sm300 2", "sm300 1"))



# long-term temperature and precipitation trend

mSubClimPlot <- mSubClim %>%
  mutate(day = day - 1) %>% 
  left_join(LTtemp %>% select(-date, LTval = value)) %>% 
  mutate(tAnom = value - LTval,
         sign = if_else(tAnom > 0, "red", "blue")) %>% 
  filter(logger %in% c("temp200cm"), 
         tempLevel == 8.5,
         dplyr::between(date, ymd("2009-03-01"), ymd("2018-12-31")))



AnnualLTtemp <- LTtemp %>% mutate(month = month(date)) %>% 
  filter(month %in% c(6,7,8,9)) %>% 
  group_by(siteID) %>% 
  summarise(LTval = mean(value)) %>% ungroup()

AnnualLTprecip <- LTprecip %>% mutate(month = month(date)) %>% 
  filter(month %in% c(6,7,8,9)) %>%  
  group_by(siteID) %>% 
  summarise(LTPval = mean(value)) %>% ungroup()

LTtemp <- LTtemp %>% mutate(month = month(date),
                  monthN = if_else(month %in% c(6,7), "spr",
                                   if_else(month %in% c(8,9), "aut", NULL))) %>% 
  filter(!is.na(monthN)) %>% 
  group_by(monthN, siteID) %>% 
  summarise(value = mean(value)) %>% 
  select(LTtval = value, monthN, siteID)

LTprecip <- LTprecip %>% mutate(month = month(date),
                  monthN = if_else(month %in% c(6,7), "spr",
                                   if_else(month %in% c(8,9), "aut", NULL))) %>% 
  filter(!is.na(monthN)) %>% 
  group_by(monthN, siteID) %>% 
  summarise(value = mean(value)) %>% 
  select(LTPval = value, monthN, siteID)


## monthly temp anomalies
soilM <- mSubClim %>% filter(logger %in% c("jordf1", "jordf2")) %>% 
  group_by(date, siteID) %>% 
  summarise(soilM = mean(value, na.rm = TRUE))

sm2018 <- SM2018 %>% 
  filter(siteID %in% c("Fauske", "Lavisdalen", "Skjellingahaugen", "Veskre"), 
         Treatment == "C") %>%
  group_by(siteID) %>% 
  summarise(soilM = mean(SM)/100) %>% 
  mutate(year = 2018,
         monthN = case_when(
           siteID == "Fauske" ~ "spr",
           siteID == "Veskre" ~ "spr",
           siteID == "Lavisdalen" ~ "spr",
           siteID == "Skjellingahaugen" ~ "aut"
         ))

# calculate monthly averages for drought and non-drought conditions
monthAv <- mSubClim %>%
  filter(logger == "temp200cm") %>%
  rename(soilT = value) %>% 
  left_join(soilM) %>% 
  gather(soilM, soilT, key = logger, value = value) %>% 
  mutate(month = month(date),
         monthN = if_else(month %in% c(6,7), "spr",
                          if_else(month %in% c(8,9), "aut", NULL))) %>% 
  filter(!is.na(monthN)) %>% 
  group_by(monthN, year, siteID, logger) %>% 
  summarise(value = mean(value, na.rm = TRUE)) %>%
  ungroup() %>% 
  spread(key = logger, value = value) %>% 
  left_join(LTtemp) %>% 
  left_join(LTprecip) %>% 
  left_join(sm2018, by = c("siteID", "year", "monthN"), suffix = c("", ".new")) %>% 
  mutate(soilM = if_else(is.na(soilM), soilM.new, soilM)) %>% 
  mutate(tAnom = soilT - LTtval,
         pAnom = soilM - LTPval) %>% 
  select(monthN, year, siteID, tAnom, pAnom, soilT, soilM)


# calculate summer averages
AnnualMonthAv <- mSubClim %>%
  filter(logger == "temp200cm") %>%
  rename(soilT = value) %>% 
  left_join(soilM) %>% 
  gather(soilM, soilT, key = logger, value = value) %>% 
  mutate(month = month(date)) %>% 
  filter(month %in% c(6,7,8,9)) %>% 
  group_by(year, siteID, logger) %>% 
  summarise(value = mean(value, na.rm = TRUE)) %>%
  ungroup() %>% 
  spread(key = logger, value = value) %>% 
  left_join(AnnualLTtemp) %>% 
  left_join(AnnualLTprecip) %>% 
  mutate(tAnom = soilT - LTval,
         pAnom = soilM - LTPval) %>% 
  select(year, siteID, tAnom, pAnom, soilT, soilM)


#save(AnnualMonthAv, file = "~/Documents/FunCaB/climate/data/AnnualMonthAv.RData")
