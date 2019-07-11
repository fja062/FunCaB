#### FUNCAB iButton Data ######
library(tidyverse)
library(lubridate)
library(readr)
library(readxl)
library(quantreg)
library(broom)

source("~/OneDrive - University of Bergen/Research/FunCaB/SeedClim-Climate-Data/funcab/dictionaries.R")

#### Read in iButtons Function
ReadIniButtons <- function(textfile){
  print(textfile)
  ending <- substr(textfile, nchar(textfile)-4+1, nchar(textfile))
  # Extract Date, Unit and Value
  if(ending == ".txt"){
    dat <- read_delim(textfile, delim = ",", col_names = FALSE)
    dat$Date <- dmy_hms(dat$X1) # convert to date
    dat$Unit <- dat$X2
    dat$Value <- as.numeric(paste(dat$X3, dat$X4, sep = ".")) # merge col 3 and 4 to one number
    dat <- dat %>% select(-X1, -X2, -X3, -X4)
  } else if(ending == ".csv"){
    # import body of data
    dat <- read_csv(textfile, skip = 19)
    dat$Date <- mdy_hms(dat$`Date/Time`) # convert to date
    dat <- dat %>% select(Date, Unit, Value, -`Date/Time`)
  } else {
    warning(paste(textfile, "format not recognised")) # warning if logger not recognized
    dat <- NULL
  }
  
  # Extract siteID, iButtonID and Year from file name
  dat$iButtonID <- basename(textfile)
  textfile2 <- tibble(dirname(textfile))
  colnames(textfile2) <- "ID"
  textfile2 <- textfile2 %>% 
    separate(col = ID, into = c("a", "b", "c", "d","e","f", "g", "Year", "siteID"), sep = "/")
  dat$siteID <- textfile2$siteID
  dat$Year <- textfile2$Year
  dat$ID <- paste0(dat$iButtonID, "_", dat$Year)
  return(dat)
}


# read in iButtonID dictionary 2016
dictionary16 <- read_excel(path = "~/OneDrive - University of Bergen/Research/FunCaB/Data/dictionaries&corrections/iButtonDictionary.xlsx", sheet = 1, col_names = TRUE)

dictionary16 <- dictionary16 %>% 
  mutate(Block = case_when(Block == "FCIII"  ~ "3",
                           Block == "FCII"   ~ "2",
                           Block ==  "FCI"   ~ "1",
                           Block == "FCV"    ~ "5",
                           Block == "FCIV"   ~ "4",
                           Block == "FCXV"   ~ "15",
                           Block == "FCXII"  ~ "12",
                           Block == "FCXIII" ~ "13",
                           Block == "FCVI"   ~ "6",
                           Block == "FCVIII" ~ "8",
                           Block == "IX"     ~ "9",
                           Block == "II"     ~ "2",
                           Block == "V"      ~ "5",
                           Block == "E2"     ~ "E2",
                           Block == "E5"     ~ "E5",
                           Block == "E9"     ~ "E9")) %>% 
  mutate(ID = paste0(iButtonID, "_", Year)) %>% 
  select(-TTC)

# read in iButtonID dictionary 2017
dictionary17 <- read_excel(path = "~/OneDrive - University of Bergen/Research/FunCaB/Data/dictionaries&corrections/iButtonDictionary.xlsx", sheet = 2, col_names = TRUE)
dictionary17 <- dictionary17 %>% 
  mutate(siteID = case_when(siteID == "LAV" ~ "Lavisdalen",
                            siteID == "GUD" ~ "Gudmedalen", 
                            siteID == "SKJ"~ "Skjellingahaugen", 
                            siteID == "ULV" ~ "Ulvhaugen", 
                            siteID == "ALR" ~ "Alrust", 
                            siteID == "FAU" ~ "Fauske", 
                            siteID == "HOG" ~ "Hogsete",
                            siteID == "VIK" ~  "Vikesland", 
                            siteID == "RAM" ~ "Rambera", 
                            siteID == "VES" ~ "Veskre", 
                            siteID == "OVS" ~ "Ovstedal", 
                            siteID == "ARH" ~ "Arhelleren")) %>% 
  mutate(ID = paste0(iButtonID, "_", Year)) %>% 
  select(-TTC)


dictionary <- dictionary16 %>% 
  bind_rows(dictionary17) %>% 
  mutate(ID = gsub(pattern = ".csv", "", ID),
         ID = gsub(pattern = "000000", "", ID),
         ID = gsub(pattern = ".txt", "", ID),
         ID = gsub(pattern = " DS1922L", "", ID)) %>% 
  filter(!Comments %in% c("FAIL", "Lost", "NO DATA", "no data on iButton", "X", "X Dead"))


# MAKE LIST OF ALL TXT FILES AND MERGE THEM TO ONE DATA FRAME
myfiles <- dir(path = paste0("/Volumes/BIO/Felles/007_Funcab_Seedclim/SeedClimClimateData/iButtondata"), pattern = "csv|txt", recursive = TRUE, full.names = TRUE)
myfiles <- myfiles[!grepl("log", myfiles, ignore.case = TRUE)] # remove log files
myfiles <- myfiles[!grepl("Empty", myfiles, ignore.case = TRUE)] # remove log files
myfiles <- myfiles[!grepl("unknown", myfiles, ignore.case = TRUE)] # remove log files


# read in ibutton logger data
mdat <- map_df(myfiles, ReadIniButtons)
head(mdat)

# clean logger ID to allow merge with the dictionary
trial <- mdat %>% 
  mutate(ID = gsub(pattern = ".csv", "", ID),
         ID = gsub(pattern = "000000", "", ID),
         ID = gsub(pattern = ".txt", "", ID),
         ID = gsub(pattern = " DS1922L", "", ID)) %>% 
  filter(!ID %in% c("093E41FB41_2017", "5E3E332741_2017", "843E3EA341_2017", "263E0D9441_2017", "FC3E48CE41second_2017", "E43E3BA541 GUD TTCP8_2016")) %>%
  mutate(ID = gsub(pattern = "123E354941_2017", "123E354941_2016", ID),
         ID = gsub(pattern = "0440282241_2017", "440282241_2017", ID),
         ID = gsub(pattern = "9C3E407D41_2017", "3E407D41_2016", ID),
         ID = gsub(pattern = "C03BBEDA41_2017", "3BBEDA41_2016", ID),
         ID = gsub(pattern = "E13E37AA41_2017", "3E37AA41_2016", ID),
         ID = gsub(pattern = "E13E37AA41_2017", "3E37AA41_2016", ID))
  


# these guys don't appear in 2017 dictionary: 7b3e410b41_2017, 783E3DA541_2017, 843E3BDE41_2017, 763E483541_2017, 5B3E34CB41_2017, 253E438A41_2017, 213BC3A141_2017, 063e3b6941_2017, 9A3E3C7041_2017 (lavisdalen), 1B40100C41_2017 (fauske), 963E328441_2017 (alrust)

# setdiff(trial$ID, dictionary$ID) <- what's in trial that's not in dictionary
# setdiff(dictionary$ID, trial$ID) <- what's in dictionary that's not in trial
# need to fix 2017 data = seems to be a lot missing from the dictionary

soilTemp <- trial %>% 
  mutate(Year = as.numeric(Year)) %>% 
  left_join(dictionary, by = c("ID", "siteID")) %>%
  select(-iButtonID.x, -iButtonID.y, -Year.x) %>% 
  mutate(Treatment = if_else(ID == "3E369341_2016", "GB", Treatment)) %>% 
  rename(Year = Year.y) %>% 
  filter(!Block %in% c("E5", "E2", "E9"))

# find out who the problematic ones are
soilTemp %>% filter(is.na(Year)) %>%
  filter(siteID == "Skjellingahaugen") %>% #distinct(siteID, ID, .keep_all = TRUE) %>% 
  ggplot(aes(x = Date, y = Value)) + geom_point() + facet_wrap(~siteID)

# Check start and end date
soilTemp %>% 
  group_by(Year, siteID) %>% 
  summarise(min(Date, na.rm = TRUE), max(Date, na.rm = TRUE)) %>% 
  arrange(siteID) %>% print(n = Inf)

# Remove ibuttons that are apparently logging until December 2017
wonkyiButtons <- soilTemp %>% 
  filter(Date > "2017-09-01 00:07:01") %>% 
  distinct(ID) %>% 
  pull(ID)

#soilTemp %>% filter(ID %in% wonkyiButtons) %>% ggplot(aes(x = Date, y = Value)) + geom_point() + facet_wrap(~ siteID)

`3E47D941_2016` <- soilTemp %>% filter(ID == "3E47D941_2016") %>% 
  mutate(Year = 2016, Block = "1", Treatment = "GF")

#Clean logger data to remove broken loggers and values before loggers were placed
iButtonData <- soilTemp %>% 
  filter(!ID %in% wonkyiButtons, !ID == "3E47D941_2016") %>% 
  bind_rows(`3E47D941_2016`) %>% 
  #filter(!ID %in% c("3E369341.csv_2016", "3E3DF841.csv_2016")) %>% # remove 2 iButtons from Gudmeldalen; these loggers need to be checked!!!
  filter(!(Date < ReplacementDate & grepl("2016", ID))) %>%  # remove values in 2015 before ibuttons have been placed
  filter(!(Date < ReplacementDate & grepl("2017", ID))) %>% 
  filter(!(Date > RemovalDate & grepl("2016", ID))) %>% 
  filter(!(Date > RemovalDate & grepl("2017", ID))) %>% 
  filter(!is.na(Treatment)) %>% 
  mutate(DOY = format(Date, "%Y-%m-%d")) %>% 
  mutate(turfID = paste0(substr(siteID,1,3), Block, Treatment)) %>% 
  select(-Unit) %>% 
  filter(!(ID == "3E3DF841_2016" & format(Date, "%Y-%m-%d") > "2016-04-25" & format(Date, "%Y-%m-%d") < "2016-05-26"),
         !(ID == "3E369341_2016" & format(Date, "%Y-%m-%d") > "2016-05-05"),
         !(ID == "D33E336C41_2017" & format(Date, "%Y-%m-%d") > "2016-12-19"),
         !(ID  %in% c("333E3B7B41_2017", "A93BD3C041_2017") & format(Date, "%Y-%m-%d") > "2017-05-22"),
         !(ID == "333E3B7B41_2017" & format(Date, "%Y-%m-%d") > "2016-10-05" & format(Date, "%Y-%m-%d") < "2016-10-26"),
         !(siteID == "Gudmedalen" & format(Date, "%Y-%m-%d") < "2016-05-26" & Value < -5),
         !(siteID == "Gudmedalen" & format(Date, "%Y-%m-%d") < "2016-05-26" & Value > 30),
         !(siteID == "Skjellingahaugen" & format(Date, "%Y-%m-%d") < "2015-10-01" & Value > 15),
         !(siteID == "Veskre" & format(Date, "%Y-%m-%d") < "2015-07-30"),
         Value - dplyr::lag(Value) < 10,
         Value > -20,
         Value < 40)

iButtonData <- iButtonData %>% 
  mutate(shortsiteID = substr(siteID, 1, 3)) %>% 
  unite(turfID, shortsiteID, Block, Treatment, sep = "", remove = FALSE) %>% 
  select(-shortsiteID) %>% 
  left_join(dict_TTC_turf)

##### save and load ibutton data #######
save(iButtonData, file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/iButton2015-2017.RData")

# yday ===> DOY

#### --- UVB LOGGER --- ####

# read in UVB data
# source the Functions_ReadInSeedClimClimate script
source("~/OneDrive - University of Bergen/Research/FunCaB/SeedClim-Climate-Data/funcab/climate/uv/Functions_ReadInUVBData.R")

UVB <- ImportData(site = c("Ovstedal", "Arhelleren", "Veskre", "Skjellingahaugen", "Rambera", "Fauske", "Alrust", "Ulvhaugen", "Vikesland", "Hogsete", "Lavisdalen", "Gudmedalen"))


uvb <- UVB %>%
  filter(!logger %in% c("power", "voltage2", "resistance1", "resistance2", "counter1", "counter2")) %>% 
  filter(value < 2.6) %>%
  filter(!(site == "skj" & value > 1)) %>% 
  filter(format(date, "%Y-%m-%d") > "2015-06-15") 

# summarise data to get one uvb value per hour.
uvb <- uvb %>%
  distinct() %>% 
  filter(!site %in% c("Skj", "Ulv")) %>%
  group_by(logger, site) %>% 
  mutate(value = value - min(value)) %>% 
  group_by(site, hour = hour(date), date = date(date)) %>%
  summarise(value = mean(value))

uvb %>% 
  group_by(site, hour) %>% 
  summarise(var = var(value)) %>% 
  arrange(var)

#save(uvb, file = "funcab_uvb.RData")

# 95% quantile regression to determine sunniness for each site at each hour of DOY
predictions <- uvb %>% 
  group_by(site, hour) %>% 
  do({
    fit <- rq(value ~ splines::bs(date, df = 10), tau = 0.95, data = .)
    augment(fit, data = .)
  })

#predictions %>% filter(hour == 12) %>% ggplot(aes(date, .fitted)) + geom_point(data = filter(uvb, hour == 12), aes(x = date, y = value, colour = factor(hour))) + geom_line() + facet_wrap(~ site)

sunniness <- predictions %>% 
  mutate(prop = value/.fitted) %>% 
  group_by(site, date) %>% 
  summarise(sunniness = weighted.mean(prop, w = .fitted)) %>% 
  mutate(sunniness = if_else(sunniness > 1, 1, sunniness)) %>%
  ungroup() %>% 
  mutate(siteID = plyr::mapvalues(site, from = dict_Site$v3, to = dict_Site$new)) %>% 
  select(-site)

sunniness <- sunniness %>% 
  complete(siteID, date) %>% 
  left_join(filter(sitePref, pref == 1)) %>% 
  left_join(sunniness, by = c("alternative" = "siteID", "date" = "date"), suffix = c(".old", ".alt")) %>% 
  mutate(sunniness = if_else(is.na(sunniness.old), sunniness.alt, sunniness.old)) %>%
  select(- sunniness.alt, - pref, - alternative)
  
gud <- sunniness %>% 
  filter(siteID == "Lavisdalen", format(date, "%Y-%m-%d") < "2015-09-09") %>%
  mutate(siteID = case_when(
    siteID == "Lavisdalen" ~ "Gudmedalen"))

UlvSkj <- sunniness %>% 
  filter(siteID %in% c("Alrust", "Rambera")) %>% 
  mutate(siteID = case_when(
    siteID == "Alrust" ~ "Ulvhaugen",
    siteID == "Rambera" ~ "Skjellingahaugen"
  ))

sunniness <- sunniness %>% 
  filter(!siteID %in% c("Ulvhaugen", "Skjellingahaugen"), !(siteID == "Gudmedalen" & format(date, "%Y-%m-%d") < "2015-09-09")) %>% 
  bind_rows(UlvSkj) %>% 
  bind_rows(gud)



# Air temperature
load("~/OneDrive - University of Bergen/Research/FunCaB/SeedClim-Climate-Data/data/Temperature.RData")

temperature <- temperature %>% 
  filter(!is.na(value)) %>% 
  #filter(date > "2015-06-15") %>% 
  #filter(date < "2017-06-13") %>% 
  mutate(site = plyr::mapvalues(site, from = dict_Site$old, to = dict_Site$new)) %>%
  filter(is.na(flag)) %>%
  filter(logger == "temp200cm") %>% 
  mutate(hour = hour(date), date = date(date)) %>% 
  select(siteID = site, date, hour, airTemp = value)


#### gridded climate data ####
load("~/OneDrive - University of Bergen/Research/FunCaB/SeedClim-Climate-Data/data/GriddedDailyClimateData2009-2016.RData")

climate <- climate %>%
  filter(between(Date, ymd("2015-07-01"), ymd("2016-08-30"))) %>% 
  mutate(Site = plyr::mapvalues(Site, from = dict_Site$old, to = dict_Site$new)) %>%
  select(-c(Temperature, Year, Month, Day), siteID = Site, gridPrecipitation = Precipitation, date = Date)
  

#----------- Data compilation -------------#
soilTemp <- iButtonData %>% 
  mutate(hour = hour(Date), date = date(Date)) %>% 
  select(-Date, -DOY, -Comments, -ReplacementDate2018, -RemovalDate, -ReplacementDate) %>% 
  bind_rows(temperature) %>% 
  left_join(sunniness) %>% 
  select(-sunniness.old) %>% 
  mutate(TOD = case_when(
    hour %in% c(10:18) ~ "day",
    hour %in% c(21,22,23,0,1,2,3,4,5,6) ~ "night",
    hour %in% c(19,20,7,8,9) ~ "spinup"))

soilTemp <- climate %>% 
  left_join(soilTemp) #%>% 
  #filter(!turfID == "Lav3GF") # removed due to irregular logging

save(soilTemp, file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/soilTemp.RData")
