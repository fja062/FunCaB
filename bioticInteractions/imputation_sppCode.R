library(tidyverse)
library(DBI)
library(dbplyr)
library(SDMTools)
library(readxl)
library(RSQLite)

con <- src_sqlite(path = "~/OneDrive - University of Bergen/Research/FunCaB/seedclim.sqlite", create = FALSE)
#con <- src_mysql(group = "seedclim", dbname = "seedclimComm", password = "password")

source("~/OneDrive - University of Bergen/Research/FunCaB/SeedClim-Climate-Data/funcab/dictionaries.R")
source("~/OneDrive - University of Bergen/Research/FunCaB/SeedclimComm/inst/graminoidRemovals/weather.R")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ---- database.controls.import ---- 
FG <- tbl(con, "character_traits") %>% 
  filter(trait == "functionalGroup") %>% 
  select(species, functionalGroup = value) %>% 
  collect()

my.GR.data <-tbl(con, "subTurfCommunity") %>%
  group_by(turfID, year, species) %>% 
  summarise(n_subturf = n()) %>% 
  collect() %>% 
  full_join(tbl(con, "turfCommunity") %>% collect()) %>%
  left_join(tbl(con, "taxon"), copy = TRUE) %>%
  left_join(tbl(con, "turfs"), copy = TRUE) %>%
  left_join(tbl(con, "plots"), by = c("destinationPlotID" = "plotID"), copy = TRUE) %>%
  left_join(tbl(con, "blocks"), by = "blockID", copy = TRUE) %>%
  left_join(tbl(con, "sites"), by = "siteID", copy = TRUE) %>%
  left_join(tbl(con, "turfEnvironment"), copy = TRUE) %>%
  select(siteID, blockID, plotID = destinationPlotID, turfID, TTtreat, GRtreat, Year = year, species, cover, temperature_level, precipitation_level, recorder, totalVascular, totalBryophytes, vegetationHeight, mossHeight, litter) %>%
  mutate(TTtreat = factor(TTtreat), GRtreat = factor(GRtreat)) %>%
  ungroup() %>% 
  filter(Year > 2014, TTtreat == "TTC"|GRtreat == "TTC")

# correct inconsistencies among cm/mm
my.GR.data <- my.GR.data %>%
  mutate(vegetationHeight = if_else(Year %in% c(2015, 2017), vegetationHeight*10, vegetationHeight),
         mossHeight = if_else(Year %in% c(2015, 2017), mossHeight*10, mossHeight)) %>% 
  mutate(mossHeight = if_else(Year == 2017 & turfID == "301 TTC", 7.5, mossHeight)) # careful!!! <- must check if this is correct

levels(my.GR.data$TTtreat) <- c(levels(my.GR.data$TTtreat),levels(my.GR.data$GRtreat))
my.GR.data$TTtreat[my.GR.data$TTtreat == ""| is.na(my.GR.data$TTtreat)] <- my.GR.data$GRtreat[my.GR.data$TTtreat == ""| is.na(my.GR.data$TTtreat)] # merge the GRtreat and TTtreat into one column
my.GR.data$GRtreat <- NULL

my.GR.data$recorder[is.na(my.GR.data$recorder)] <- "unknown botanist"

###--- fixes for botanist biases ---###
# Pascale fix
my.GR.data$cover[my.GR.data$recorder == "PM"] <- my.GR.data$cover[my.GR.data$recorder=="PM"]*1.20

# Siri fix
siri <- my.GR.data %>%
  filter(recorder == "Siri") %>%
  group_by(turfID, Year) %>%
  mutate(SumOfcover = sum(cover)) %>%
  filter(SumOfcover/totalVascular < 1.35)

siri.fix <- paste(as.character(my.GR.data$turfID), my.GR.data$Year) %in% paste(siri$turfID, siri$Year)
my.GR.data$cover[siri.fix] <- my.GR.data$cover[siri.fix]*1.3

# Owen fix
owen <- my.GR.data %>% 
  filter(recorder == "Owen") %>% 
  group_by(turfID, Year) %>% 
  mutate(sumOfCover = sum(cover)) %>% 
  filter(sumOfCover/totalVascular > 1.5)

owen.fix <- paste(as.character(my.GR.data$turfID), my.GR.data$Year) %in% paste(owen$turfID, owen$Year)
my.GR.data$cover[owen.fix] <- my.GR.data$cover[owen.fix]/1.5

my.GR.data <- my.GR.data %>%
  filter(turfID %in% dict_TTC_turf$TTtreat) %>% # or semi_join()
  mutate(Treatment = "C", TTtreat = turfID) %>%
  left_join(dict_TTC_turf, by = "TTtreat", suffix = c(".new", "")) %>% 
  mutate(blockID = substr(turfID, 4, 4)) %>% 
  select(-c(plotID, temperature_level, precipitation_level, totalVascular, litter, turfID.new)) %>% 
  filter(!is.na(cover),
         !(TTtreat == "37 TTC" & Year > 2015))


####-------- load funcab data ---------####

gudfun2015 <- read_excel("~/OneDrive - University of Bergen/Research/FunCaB/Data/primary/veg_composition/funcab_Gudmedalen.xlsx", col_types = "text")

funcab_2015 <- read_delim("~/OneDrive - University of Bergen/Research/FunCaB/Data/primary/veg_composition/funcab_composition_2015-utenGud.csv", delim = ";", col_types = cols(.default = "c")) # ; or \t

funcab_2016 <- read_delim("~/OneDrive - University of Bergen/Research/FunCaB/Data/primary/veg_composition/funcab_composition_2016.csv", delim = ";", col_types = cols(.default = "c"))

funcab_2017 <- read_delim("~/OneDrive - University of Bergen/Research/FunCaB/Data/primary/veg_composition/funcab_composition_2017.csv", delim = ";", col_types = cols(.default = "c"))

funcab_2018 <- read_excel("~/OneDrive - University of Bergen/Research/FunCaB/Data/primary/veg_composition/funcab_composition_2018.xlsx", col_types = "text")

scBryo <- read_excel("~/OneDrive - University of Bergen/Research/FunCaB/Data/primary/veg_composition/2017seedclimBryophyte.xlsx")

# calculate mean veg and moss heights for 2018
funcab_2018 <- funcab_2018 %>% 
  filter(!grepl("TT1", TTtreat), !grepl("OUT", TTtreat), !grepl("OUT", turfID)) %>% 
  mutate_at(vars(vegetationHeight, mossHeight), as.numeric) %>% 
  group_by(turfID) %>% 
  mutate(vegetationHeight = as.character(mean(vegetationHeight, na.rm = TRUE)),
         mossHeight =  as.character(mean(mossHeight, na.rm = TRUE))) %>% 
  ungroup()

# replace species names where mistakes have been found in database

problems <- read.csv("~/OneDrive - University of Bergen/Research/FunCaB/Data/dictionaries&corrections/speciesCorrections.csv", sep = ";", stringsAsFactors = FALSE) %>%
  filter(!old %in% c("Vio.can", "Com.ten", "Sel.sel")) %>%
  filter(cover != "WHAT HAPPENED") %>%
  mutate(cover = as.numeric(cover))

# load the dictionary merger
mergedictionary <- read.csv2(file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/dictionaries&corrections/mergedictionary.csv")

prob.sp <- problems %>%
  filter(!is.na(Year)) %>% 
  select(-functionalGroup)

prob.sp.name <- problems %>% 
  filter(is.na(Year), !old %in% c("Eri.bor")) %>% 
  select(old, new) %>% 
  bind_rows(mergedictionary)

problems.cover <- filter(problems, !is.na(cover)) %>%
  select(turfID, year = Year, species = old, cover)


# bind composition data, replace _ with . for compatibility in spp names
composition <- funcab_2016 %>% 
  bind_rows(funcab_2015) %>% 
  bind_rows(gudfun2015) %>% 
  bind_rows(funcab_2017) %>% 
  bind_rows(funcab_2018) %>% 
  filter(subPlot %in% c("%", "T")) %>% 
  select(c(siteID:subPlot), Year = year, recorder, c(totalGraminoids:mossHeight), litter, acro, pleuro, c(`Ach mil`:`Vis vul`)) %>% 
  select_if(colSums(!is.na(.)) > 0) %>% 
  gather(c("Ach mil":"Vio sp"), key = "species", value = "cover")

subTurfFreq <- composition %>% filter(subPlot == "T", !is.na(cover)) %>% 
  select(siteID, Treatment, turfID, Year, species, presence = cover) %>% 
  mutate(presence = 1)

composition <- composition %>% 
  filter(subPlot == "%") %>% 
  left_join(subTurfFreq) %>% 
  mutate(species = gsub("\\ |\\_", ".", species)) %>% 
  left_join(dict_TTC_turf, by = c("turfID" = "TTtreat"), suffix = c(".old", ".new")) %>%
  mutate(turfID = if_else(!is.na(turfID.new), turfID.new, turfID)) %>% 
  mutate_at(vars(cover, Year, totalGraminoids:pleuro), as.numeric) %>% 
  mutate(turfID = if_else((blockID == 16 & siteID == "Gudmedalen"), gsub("16", "5", turfID), turfID),
         turfID = if_else((siteID == "Alrust" & blockID == "3" & Year == 2015 & Treatment == "C"), "Alr3C", turfID),
         turfID = recode(turfID, "Alr4FGB" = "Alr5C"),
         turfID = recode(turfID, "Lav1G " = "Lav1G"),
         blockID = if_else(blockID == 16 & siteID == "Gudmedalen", gsub("16", "5", blockID), blockID),
         blockID = if_else(turfID == "Gud12C", "12", blockID)) %>% 
  filter(!(blockID == "4" & Year == 2015 & siteID == "Alrust"),
         !(turfID =="Gud12C" & Year == 2015),
         !is.na(turfID),
         !(turfID == "Alr3C" & recorder == "Siri"))



# overwrite problem spp with their correct names and covers
composition <- composition %>% 
  left_join(prob.sp, by = c("Year", "turfID", "siteID", "species" = "old"), suffix = c("", ".new")) %>%
  mutate(species = coalesce(new, species),
         cover = coalesce(cover.new, cover)) %>% 
  select(-new, -cover.new, -subPlot, - turfID.new) %>% 
  left_join(prob.sp.name, by = c("species" = "old")) %>% 
  mutate(species = if_else(!is.na(new), new, species))

FGBs <- composition %>% 
  filter(Treatment %in% c("FGB", "GF")) %>% 
  select(-species, -cover) %>% 
  distinct() %>% 
  filter(Year > 2015)

# filter out funcab controls that are also TTCs in 2015 & 2016
ttcs1516 <- composition %>% 
  filter(Treatment == "C", !Year %in% c(2017, 2018), !is.na(Year)) %>% 
  right_join(dict_TTC_turf) %>%
  select(-species, -cover, -pleuro, -acro, -litter, -presence, -new, -recorder) %>% 
  distinct()

ttcs17 <- composition %>% 
  filter(Treatment == "C", Year == 2017) %>% 
  right_join(dict_TTC_turf) %>%
  select(-new, -cover, -presence, -species) %>% 
  distinct() %>% 
  full_join(scBryo, by = "turfID", suffix = c(".old", "")) %>% 
  select(-totalBryophytes.old, -mossHeight.old, -vegetationHeight.old, -TTtreat.old, -litter)

####################
#### clean data ####

# join with TTC data
comp2 <- composition %>% 
  mutate(blockID = if_else(nchar(blockID) > 1, gsub("[^[:digit:]]", "", blockID), blockID)) %>% 
  full_join(my.GR.data, by = c("siteID", "blockID", "turfID", "Treatment", "Year", "species", "recorder"), suffix = c("", ".new")) %>% 
  mutate(cover = coalesce(cover.new, cover),
         totalBryophytes = coalesce(totalBryophytes.new, totalBryophytes),
         vegetationHeight = coalesce(vegetationHeight.new, vegetationHeight),
         mossHeight = coalesce(mossHeight.new, mossHeight)) %>% 
  select(-cover.new, -totalBryophytes.new, -vegetationHeight.new, -mossHeight.new)


#mean of previous and next year
sampling_year <- comp2 %>% 
  group_by(turfID) %>% 
  distinct(turfID, Year) %>% 
  arrange(turfID, Year) %>% 
  mutate(sampling = 1:n())

missingCov <- comp2 %>% group_by(turfID, species, Treatment) %>% filter(!is.na(presence) & is.na(cover)) %>% 
  select(siteID, blockID, turfID, Year, species, Treatment)

# covers interpolated from cover in year before/after
missingCov <- missingCov %>% 
  left_join(sampling_year) %>% 
  left_join(
    left_join(filter(comp2, !is.na(cover)), sampling_year),
    by = c("turfID", "species"), 
    suffix = c("", "_cover")) %>% #join to other years
  filter(abs(sampling - sampling_cover) == 1) %>% #next/previous year
  group_by(siteID, blockID, Treatment, turfID, species, Year) %>% 
  filter(n() == 2) %>% #need before and after year
  summarise(cover = mean(cover), flag = "Subturf w/o cover. Imputed as mean of adjacent years")

misCovSpp <- comp2 %>% filter(Treatment == "XC") %>% 
  filter(is.na(cover) & !is.na(presence)) %>% 
  distinct(siteID, blockID, Treatment, turfID, Year, species)

# covers interpolated from site means
misCovSpp2 <- comp2 %>% 
  right_join(misCovSpp %>% select(siteID, Year, species)) %>% 
  group_by(siteID, species, Year) %>% 
  summarise(cover = mean(cover, na.rm = TRUE)) %>% 
  filter(!is.na(cover)) %>% 
  right_join(misCovSpp)

# adding cover corrections
comp2 <- comp2 %>% 
  left_join(missingCov %>% select(-flag), by = c("siteID", "blockID", "Treatment", "turfID", "species", "Year"), suffix = c("", ".new")) %>% 
  mutate(cover = coalesce(cover.new, cover)) %>% 
  select(-cover.new) %>% 
  left_join(misCovSpp2, by = c("siteID", "blockID", "Treatment", "turfID", "species", "Year"), suffix = c("", ".new")) %>% 
  mutate(cover = coalesce(cover.new, cover)) %>% 
  select(-cover.new, -presence) %>% 
  filter(cover > 0)

# rejoin funcab attributes of the TTCs in 2016 and 2017
comp2 <- comp2 %>% 
  left_join(ttcs1516, by = c("siteID", "blockID", "Treatment", "turfID", "Year", "TTtreat"), suffix = c("", ".new")) %>% 
  mutate(mossHeight = coalesce(mossHeight.new, mossHeight),
         vegetationHeight = coalesce(vegetationHeight.new, vegetationHeight),
         totalBryophytes = coalesce(totalBryophytes.new, totalBryophytes),
         totalGraminoids = coalesce(totalGraminoids.new, totalGraminoids),
         totalForbs = coalesce(totalForbs.new, totalForbs)) %>% 
  select(-totalBryophytes.new, -vegetationHeight.new, -mossHeight.new, -totalForbs.new, -totalGraminoids.new) %>%
  left_join(ttcs17, by = c("siteID", "blockID", "Treatment", "turfID", "Year", "TTtreat", "recorder", "acro", "pleuro"), suffix = c("", ".new")) %>% 
  mutate(mossHeight = coalesce(mossHeight.new, mossHeight),
         vegetationHeight = coalesce(vegetationHeight.new, vegetationHeight),
         totalBryophytes = coalesce(totalBryophytes.new, totalBryophytes),
         totalGraminoids = coalesce(totalGraminoids.new, totalGraminoids),
         totalForbs = coalesce(totalForbs.new, totalForbs)) %>% 
  select(-totalBryophytes.new, -vegetationHeight.new, -mossHeight.new, -totalForbs.new, -totalGraminoids.new, -TTtreat, -new)

comp2 <- comp2 %>% 
  group_by(turfID, Year) %>% 
  mutate(totalBryophytes = if_else(is.na(totalBryophytes) & !is.na(pleuro) & !is.na(acro), pleuro + acro, totalBryophytes)) %>% 
  ungroup() %>% 
  mutate(turfID = if_else(grepl("TTC", turfID), turfID, substring(turfID, 4, n())),
         Treatment = gsub(" ", "", Treatment),
         turfID = paste0(str_sub(siteID, 1, 3), turfID),
         species = gsub(" ", ".", species),
         blockID = if_else(turfID == "Gud12C", "12", blockID))


# functional groups
comp2 <- comp2 %>% 
  group_by(turfID, Year, species) %>% 
  mutate(cover = case_when(
    turfID == "Alr2XC" & Year == 2016 & species == "Agr.cap" ~ sum(cover),
    TRUE ~ cover
  )) %>% 
  ungroup() %>% 
  left_join(FG)

comp2 %>% group_by(turfID, Year, species) %>% summarise(n = n_distinct(cover)) %>% filter(n > 1)

compiledResponseData <- comp2 %>% 
  filter(!Year == 2018) %>%
  left_join(weather) %>% 
  mutate(Site = substr(siteID, 1,3),
         turfID_Y = paste0(turfID, "_", Year)) %>% 
  distinct(turfID, Year, siteID, Site, turfID_Y, species, functionalGroup, cover, temp.C = temp7010, P.mm = precip7010)

save(compiledResponseData, file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/compiledResponseData.RData")
