# seedling data cleaning #

## -----> extract gridded precipitation ; where is this? ask aud


# D = dead, assign alive if found alive in following year
# S = seedling
# TP = toothpick missing; essentially missing data. assign alive if found alive in following year
# E = established - two pairs of opposite leaves, or two alternate leaves
# + = fertile


#load packages
library(lubridate)

#
# source vegetation data, plotting code and soil temperature data
source("~/Documents/FunCaB/bioticInteractions/composition/dataProcessing/00funcab_data_processing.R")

source("~/Documents/FunCaB/figures/plotting_dim.R")

source("~/Documents/FunCaB/climate/LT_climate.R")

#load seedling data
seed <- read_csv2("~/OneDrive - University of Bergen/Research/FunCaB/Data/primary/veg_recruitment/2018_funcab_seedlings.csv")


# load seed mass trait data
con <- src_sqlite(path = "~/OneDrive - University of Bergen/Research/FunCaB/seedclim.sqlite", create = FALSE)



######### data preparation ###########
seedMass <- tbl(con, "numeric_traits") %>% 
  filter(trait == "seedMass") %>% 
  collect()

turfDict <- comp2 %>% 
  distinct(siteID, blockID, Treatment, turfID) %>% 
  filter(Treatment %in% c("G", "B", "GB", "C", "FGB")) %>% 
  mutate(Round = 1, Round2 = 2) %>% 
  gather(Round, Round2, key = "v", value = "Round") %>% 
  select(-v)


# merge composition data with seed mass data
seedcomp <- comp2 %>% 
  filter(Year == 2018) %>% 
  left_join(seedMass) %>% 
  group_by(siteID, blockID, turfID, Treatment, mossCov, vegetationHeight, mossHeight, forbCov, graminoidCov, functionalGroup) %>% 
  summarise(seedMass = weighted.mean(value, cover)) %>% 
  ungroup()


seed <- seed %>% 
  filter(is.na(NS)) %>% 
  filter(!Comment %in% c("out of plot")) %>% 
  mutate(Round = factor(Round),
         blockID = as.character(blockID),
         Date1 = dmy(Date1),
         Date2 = dmy(Date2),
         Leuc_sp = as.numeric(Leuc_sp), 
         Tara_sp = as.numeric(Tara_sp),
         survival = if_else(is.na(survival), "alive", survival))

# join onto complete turf list to catch turfs with zero seedlings


######### ------ abundance ------ ###########
seedTot <- seed %>% 
  group_by(siteID, Round, blockID, Treatment, turfID) %>% 
  summarise(seed = sum(n())) %>% 
  ungroup() %>% 
  distinct(siteID, Round, blockID, Treatment, turfID, seed) %>% 
  right_join(turfDict %>% mutate(Round = as.factor(Round))) %>% 
  mutate(seed = case_when(is.na(seed) ~ 0L, 
                             TRUE ~ seed)) %>% 
  group_by(turfID, Round) %>% 
  mutate(seedSum_m2 = seed/0.0625) %>% 
  ungroup() 

seedTot <- seedTot %>% 
  #left_join(seedcomp, by = c("siteID", "turfID", "blockID", "Treatment")) %>% 
  mutate(Treatment = case_when(
    Treatment == "FGB" ~ "Gap",
    Treatment == "C" ~ "Intact",
    Treatment == "B" ~ "GF",
    Treatment == "GB" ~ "F",
    Treatment == "G" ~ "FB",
    TRUE ~ Treatment
  )) %>%
  #spread(Round, seedSum) %>% 
  #mutate(perChange = ((`2`/sum(`1`, `2`, na.rm = TRUE)) - (`1`/sum(`1`, `2`, na.rm = TRUE)))*100,
  #       ) %>% 
  #gather(`1`, `2`, key = season, value = seed) %>% 
  mutate(monthN = recode(Round, `1` = "spr", `2` = "aut"),
         year = 2018) %>% 
  select(-Round, -seedSum_m2)

# weather stuff
#seedTot <- seedTot %>%
#  left_join(weather) %>% 
#  mutate(
#    precip7010 = precip7010 / 1000,
#    sprecip7010 = scale(precip7010, center = TRUE, scale = FALSE),
#    stemp7010 = scale(temp7010, center = TRUE, scale = FALSE)
#  )



######### ----- survival ----- ###########
survival <- seed %>%
  filter(Round == 1) %>% 
  select(siteID, blockID, Treatment, turfID, survival) %>% 
  mutate(survival = recode(survival, "dead" = 0, "alive" = 1, "not found" = 0, "1" = 1, "dezd" = 0)) %>%
  group_by(Treatment, turfID, blockID, siteID) %>% 
  summarise(tot = n(),
            totS = sum(survival),
            survival = (totS/tot)*100) %>% 
  right_join(turfDict %>% mutate(Round = as.factor(Round))) %>% 
  mutate(tot = case_when(is.na(tot) ~ 0L, 
                             TRUE ~ tot)) %>% 
  mutate(totS = case_when(is.na(totS) ~ 0, 
                             TRUE ~ totS)) %>% 
  mutate(survival = case_when(is.na(survival) ~ 0, 
                             TRUE ~ survival)) %>%
  left_join(weather) %>% 
  ungroup()

survival <- survival %>% 
  mutate(Treatment = recode(Treatment, "C" = "aC"),
         tempLevelPlot = factor(tempLevel, labels = c("6.5" = "Alpine", "8.5" = "Sub-alpine", "10.5" = "Boreal")),
         precipLevelPlot = factor(precipLevel, labels = c("600" = "Dry", "1200" = "Semi-dry", "2000" = "Semi-wet", "2700" = "Wet")),
         sprecip7010 = scale(precip7010/1000, center = TRUE, scale = FALSE),
         stemp7010 = scale(temp7010, center = TRUE, scale = FALSE),
         survival = survival/100)


######### ---- 2009/2010 ----- ###########

rc_rtc <- read_delim("~/OneDrive - University of Bergen/Research/FunCaB/Data/primary/veg_recruitment/PM_rawdata_1112.csv", delim = ";", col_types = cols(.default = "c"))

rc_rtc <- rc_rtc %>% 
  mutate(Site = case_when(
    grepl("\\vs", Site) ~ "Ovs", 
    grepl("^\\L", Site) ~ "Lav", 
    grepl("\\lr", Site) ~ "Alr",
    grepl("\\g", Site) ~ "Hog",
    TRUE ~ Site
    ),
    Site = trimws(Site)) %>% 
  mutate(Species_11_2 = if_else(is.na(Species_11_2), Species_11_1, Species_11_2),
         Species_11_2 = if_else(is.na(Species_11_2), Species_10_2_1, Species_11_2),
         Species_11_2 = if_else(is.na(Species_11_2), Species_09_2, Species_11_2)) %>% 
  select(siteID = Site, Date, Treatment = Treat, blockID = Block, ID, species = Species_11_2, 
         aut_09 = C_09_2, spr_10 = C_10_1, aut_10 = C_10_2, spr_11 = C_11_1, aut_11 = C_11_2, spr_12 = C_12_1) %>% 
  mutate(blockID = as.numeric(case_when(
    blockID == "I" ~ 1,
    blockID == "II" ~ 2,
    blockID == "III" ~ 3,
    blockID == "IV" ~ 4,
    blockID == "V" ~ 5
  ))) %>% 
  filter(Treatment %in% c("RTG", "RTC")) %>% 
  left_join(dict_Site, by = c("siteID"  = "old")) %>% 
  select(Date:spr_12, "siteID" = new)

# remove graminoid seedlings from analyses
rc_rtcSum <- rc_rtc %>% 
  mutate(species = trimws(species, which = "both"),
         species = gsub(" ", ".", species)) %>% 
  filter(!species %in% c("G", "Ant.odo", "Ave.fle"), !grepl("Agr", species), !grepl("Car", species), !grepl("Des", species), !grepl("Fes", species), !grepl("Poa", species), !grepl("Luz", species)) %>%
  mutate(species = gsub(".cf", "", species),
         species = gsub(".CF", "", species)) %>%
  gather(spr_12, aut_11, spr_11, aut_10, spr_10, aut_09, key = season, value = seed, na.rm = TRUE)

# recode seedling survival
rc_rtcSum <- rc_rtcSum %>% 
  mutate(seed = substr(seed, 1,1),
         seed = case_when(
           seed == "s" ~ "S",
           seed == "H" ~ "S",
           seed == "?" ~ "T",
           seed == "C" ~ "S",
           seed == "E" ~ "S",
           TRUE ~ seed)) %>%
  mutate(year = as.numeric(paste0("20", substr(season, 5,6))),
         Date = paste0(substr(Date, 1,6), year),
         Date = dmy(Date)) %>% 
  mutate(seed = as.numeric(case_when(
  seed == "S" ~ 1,
  seed == "T" ~ 0,
  seed == "D" ~ 0,
  seed == "R" ~ 0))) %>% 
  filter(!is.na(ID)) %>% 
  distinct(siteID, blockID, Treatment, ID, season, species, seed) %>% 
  group_by(siteID, blockID, Treatment, ID) %>% 
  mutate(sum = sum(seed)) %>% 
  ungroup() %>% 
  spread(key = season, value = seed) %>% 
  mutate(spr_12 = if_else(rowSums(.[c("aut_10", "aut_09", "spr_10", "spr_11", "aut_11")], na.rm = TRUE) > 0,  0, spr_12),
         aut_11 = if_else(rowSums(.[c("aut_10", "aut_09", "spr_10", "spr_11")], na.rm = TRUE) > 0,  0, aut_11),
         spr_11 = if_else(rowSums(.[c("aut_10", "aut_09", "spr_10")], na.rm = TRUE) > 0,  0, spr_11),
         aut_10 = if_else(rowSums(.[c("aut_09", "spr_10")], na.rm = TRUE) > 0,  0, aut_10),
         spr_10 = if_else(rowSums(.["spr_10"], na.rm = TRUE) > 0,  0, spr_10))


rc_rtcSum <- rc_rtcSum %>% 
  gather(spr_12, aut_11, spr_11, aut_10, spr_10, aut_09, key = season, value = seed, na.rm = TRUE) %>% 
  mutate(turfID = paste0(substr(siteID,1,3), blockID, Treatment),
         year = as.numeric(paste0("20",substr(season, 5,6))), 
         monthN = substr(season, 1, 3)) %>%
  group_by(turfID, siteID, blockID, Treatment, year, season, monthN) %>% 
  summarise(seed = sum(seed, na.rm = TRUE)) %>% 
  #group_by(year, turfID) %>% 
  #spread(seas, seedSum) %>% 
  #mutate(perChange = (aut/sum(spr, aut, na.rm = TRUE) - spr/sum(spr, aut, na.rm = TRUE))*100) %>% 
  #gather(aut, spr, key = season, value = seed) %>% 
  ungroup() %>% 
  mutate(season = dmy(recode(season, aut_09 = "01-09-2009", spr_10 = "01-04-2010", aut_10 = "01-09-2010", spr_11 = "01-04-2011", aut_11 = "01-09-2011", spr_12 = "01-04-2012")),
         Treatment = factor(Treatment, labels = c(RTC = "Intact", RTG = "Gap")),
         blockID = as.character(blockID)) %>% 
  select(-season)


rc_rtcSumAv <- rc_rtcSum %>%
  bind_rows(seedTot) %>% 
  left_join(monthAv) %>% 
  left_join(weather) %>% 
  mutate(sprecip7010 = as.vector(scale((precip7010 / 1000), scale = FALSE, center = TRUE)),
         stemp7010 = as.vector(scale(temp7010, scale = FALSE, center = TRUE)),
         tempLevelPlot = factor(tempLevel, labels = c("6.5" = "Alpine", "8.5" = "Sub-alpine", "10.5" = "Boreal")),
         precipLevelPlot = factor(precipLevel, labels = c("600" = "Dry", "1200" = "Semi-dry", "2000" = "Semi-wet", "2700" = "Wet"))
         )





save(rc_rtcSumAv, file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/cleanedAbundData.RData")

save(survival, file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/cleanedSurvData.RData")
