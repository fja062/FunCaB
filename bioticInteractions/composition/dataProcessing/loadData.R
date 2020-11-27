#################################################################
# Script to process community data associated with Graminoid removal Ecography MS
#################################################################

library(tidyverse)
library(DBI)
library(dbplyr)
#library(SDMTools)
library(RSQLite)

con <- src_sqlite(path = "~/OneDrive - University of Bergen/Research/FunCaB/seedclim.sqlite", create = FALSE)
#con <- src_mysql(group = "seedclim", dbname = "seedclimComm", password = "password")

#source joining dictionaries
source("~/Documents/FunCaB/dictionaries.R")

#~~~~~~~~~~ Cover data ~~~~~~~~~~#
## ---- my.GR.data.import ---- 

problems <- read.csv("~/OneDrive - University of Bergen/Research/FunCaB/Data/dictionaries&corrections/speciesCorrections.csv", sep = ";", stringsAsFactors = FALSE) %>%
  filter(!old %in% c("Vio.can", "Com.ten", "Sel.sel")) %>%
  filter(cover != "WHAT HAPPENED") %>%
  mutate(cover = as.numeric(cover))

problems.cover <- filter(problems, !is.na(cover)) %>%
  select(turfID, year = Year, species = old, cover)


# extract functional groups from DB
FG <- tbl(con, "character_traits") %>% 
  filter(trait == "functionalGroup") %>% 
  select(species, functionalGroup = value) %>% 
  collect()

# extract community data for 2010-2016
my.GR.data <- tbl(con, "subTurfCommunity") %>%
  group_by(turfID, year, species) %>% 
  summarise(n_subturf = n()) %>% 
  collect() %>% 
  full_join(tbl(con, "turfCommunity") %>% collect()) %>%
  full_join(problems.cover, by = c("year", "turfID", "species"), suffix = c(".community", ".problems")) %>%
  mutate(cover = if_else(is.na(cover.community),
                         cover.problems,
                         cover.community)) %>% 
  left_join(tbl(con, "taxon"), copy = TRUE) %>%
  left_join(tbl(con, "turfs"), copy = TRUE) %>%
  left_join(tbl(con, "plots"), by = c("destinationPlotID" = "plotID"), copy = TRUE) %>%
  left_join(tbl(con, "blocks"), by = "blockID", copy = TRUE) %>%
  left_join(tbl(con, "sites"), by = "siteID", copy = TRUE) %>%
  left_join(tbl(con, "turfEnvironment"), copy = TRUE) %>%
  select(siteID, blockID, plotID = destinationPlotID, turfID, TTtreat, GRtreat, Year = year, species, cover, recorder, totalVascular, totalBryophytes, vegetationHeight, mossHeight, litter, pleuro, acro, liver, lichen, soil, rock, totalLichen, comment, date) %>%
  mutate(TTtreat = factor(TTtreat), GRtreat = factor(GRtreat)) %>%
  ungroup() %>%
  filter(Year %in% c(2010, 2011, 2012, 2013, 2015, 2016), TTtreat == "TTC"|GRtreat == "RTC"|GRtreat == "TTC")

my.GR.data <- my.GR.data %>%
  left_join(FG) %>% 
  mutate(functionalGroup = if_else(species %in% c("Gen.sp.", "Cre.pal", "Frag.vir", "Sch.gig", "Ste.bor", "Hie.ore", "Sel.sel."), "forb", functionalGroup),
         functionalGroup = if_else(species %in% c("Agr.can", "Phl.sp"), "graminoid", functionalGroup))

levels(my.GR.data$TTtreat) <- c(levels(my.GR.data$TTtreat),levels(my.GR.data$GRtreat))
my.GR.data$TTtreat[my.GR.data$TTtreat == ""| is.na(my.GR.data$TTtreat)] <- my.GR.data$GRtreat[my.GR.data$TTtreat == ""| is.na(my.GR.data$TTtreat)] # merge the GRtreat and TTtreat into one column
my.GR.data$GRtreat <- NULL
my.GR.data <- my.GR.data[!(my.GR.data$blockID == "Gud5" & my.GR.data$Year == 2010), ]
#my.GR.data <- my.GR.data[!(my.GR.data$turfID == "Fau1RTC" & my.GR.data$Year == 2010), ]
my.GR.data <- my.GR.data[!(my.GR.data$functionalGroup == "graminoid" & my.GR.data$Year > 2011 & my.GR.data$TTtreat == "RTC"), ]
my.GR.data$Year[my.GR.data$Year == 2010] <- 2011


my.GR.data <- my.GR.data %>% 
  mutate(turfID = recode(turfID, "Ram4RTCx" = "Ram4RTC"),
         turfID = recode(turfID, "Ram5RTCx" = "Ram5RTC")) %>% 
  mutate(ID = factor(paste(turfID, Year, sep = "_"))) %>% 
  filter(!(blockID %in% c("Skj11", "Skj12", "Gud11", "Gud12", "Gud13")))


my.GR.data$recorder[is.na(my.GR.data$recorder)] <- "unknown botanist"
my.GR.data$cover[my.GR.data$recorder == "PM"] <- my.GR.data$cover[my.GR.data$recorder=="PM"]*1.20

siri <- my.GR.data %>%
  filter(recorder == "Siri") %>%
  group_by(turfID, Year) %>%
  mutate(SumOfcover = sum(cover)) %>%
  filter(SumOfcover/totalVascular < 1.35)

siri.fix <- paste(as.character(my.GR.data$turfID), my.GR.data$Year) %in% paste(siri$turfID, siri$Year)
my.GR.data$cover[siri.fix] <- my.GR.data$cover[siri.fix]*1.3

owen <- my.GR.data %>% 
  filter(recorder == "Owen") %>% 
  group_by(turfID, Year) %>% 
  mutate(sumOfCover = sum(cover)) %>% 
  filter(sumOfCover/totalVascular > 1.5)

owen.fix <- paste(as.character(my.GR.data$turfID), my.GR.data$Year) %in% paste(owen$turfID, owen$Year)
my.GR.data$cover[owen.fix] <- my.GR.data$cover[owen.fix]/1.5


#gridded temperature etc
load(file = "~/Documents/FunCaB/climate/data/AnnualMonthAv.RData")
source("~/Documents/FunCaB/climate/weather.R")


my.GR.data <- my.GR.data %>%
  left_join(AnnualMonthAv %>% rename(Year = year), by = c("siteID", "Year")) %>% 
  left_join(weather, by = "siteID") %>% 
  ungroup()


#### code to get rid of turfs that have been attacked by ants or cows ####
lowcover <- my.GR.data %>% 
  group_by(turfID, Year, functionalGroup) %>% 
  mutate(sumcover = sum(cover)) %>% 
  filter(sumcover < 25) %>% 
  distinct(siteID, turfID,Year,sumcover)

my.GR.data %>% 
  filter(turfID %in% c("Ovs2RTC", "Ovs3RTC", "126 TTC"), functionalGroup == "forb") %>% 
  group_by(turfID, Year) %>% 
  mutate(sumcover = sum(cover)) %>% 
  distinct(siteID, turfID,Year, sumcover)

#my.GR.data <- filter(my.GR.data, !siteID %in% c("Ovs2RTC", "Ovs3RTC"))
## ---- my.GR.data.end ----

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############### TRAITS ###############

## ---- Traits.data.import ---- 

# source Ragnhild's trait data
source("~/Documents/SeedclimComm/gramRem_load-traits.R")

#load from data base
traits <- tbl(con, "numeric_traits") %>% 
  collect() %>% 
  filter(trait %in% c("Lower", "Nem", "BNem", "SBor", "MBor", "Nbor", "LAlp", "MAlp", "HAlp", "Upper"))

traits <- traits %>% 
  group_by(species) %>% 
  mutate(abundance = sum(value)) %>%
  spread(trait, value) %>% 
  mutate(specialism = factor(case_when(
    Nem == 0 & BNem == 0 & SBor == 0 & LAlp == 1 ~ "alpine", 
    HAlp == 0 & MAlp == 0 & LAlp == 0 & Nem == 1 ~ "lowland", 
    abundance >= 5.5 ~ "generalist",
    abundance < 5.5 ~ "other"))) %>% 
  select(-c(Nem, BNem, SBor, MBor, Nbor, LAlp, MAlp, HAlp)) %>% 
  ungroup()


# adding traits to my.GR.data
my.GR.data <- my.GR.data %>%
  left_join(traits, by = "species") %>%
  left_join(traitdata, by = c("species","siteID"))
  #left_join(Species_traits, by = c("species", "siteID")) 
  
# renaming pteridophytes and woody species to be forbs
my.GR.data <- my.GR.data %>% 
  mutate(TTtreat = factor(TTtreat),
         functionalGroup = recode(functionalGroup, "pteridophyte" = "forb", "woody" = "forb"),
         specialism = recode(specialism, "lowland" = "other"))


## ---- Traits.data.end ---- 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############### DIVERSITY MEASURES ###############

## ---- Diversity.data.import ---- 

#Species richness
library(vegan)

my.GR.data <- my.GR.data %>%
  group_by(turfID, Year, functionalGroup) %>%
  mutate(richness = sum(n_distinct(species))) %>% 
  mutate(diversity = diversity(cover, index = "shannon")) %>% 
  mutate(evenness = (diversity/log(richness)))


## ---- Diversity.data.end ---- 

#~~~~~~~~~~~~ WEIGHTED MEANS ###############

## ---- Community.mean.weighting ---- 

###### weighted means for whole community, grouped by functional group and year
# hadley wickham's bigvis package
weighted.var <- function(x, w = NULL, na.rm = FALSE) {
  if (na.rm) {
    na <- is.na(x) | is.na(w)
    x <- x[!na]
    w <- w[!na]
  }
  
  sum(w * (x - weighted.mean(x, w)) ^ 2) / (sum(w) - 1)
}

wholecom <- my.GR.data %>% 
  group_by(ID, functionalGroup) %>% 
  mutate(sumcover = sum(cover),
         wmeanLDMC = weighted.mean(LDMC_mean, cover, na.rm = TRUE),
         wmeanSLA = weighted.mean(SLA_mean, cover, na.rm = TRUE),
         wmeanLTH = weighted.mean(Lth_mean, cover, na.rm = TRUE),
         wmeanLA = weighted.mean(LA_mean, cover, na.rm = TRUE),
         wmeanheight = weighted.mean(Height_mean, cover, na.rm = TRUE),
         wmeanCN = weighted.mean(CN_mean, cover, na.rm = TRUE),
         wmeanN = weighted.mean(N_mean, cover, na.rm = TRUE),
         cwvLDMC = weighted.var(LDMC_mean, cover),
         cwvSLA = weighted.var(SLA_mean, cover),
         cwvLTH = weighted.var(Lth_mean, cover),
         cwvLA = weighted.var(LA_mean, cover),
         cwvheight = weighted.var(Height_mean, cover),
         cwvN = weighted.var(N_mean, cover)) %>% 
  ungroup() %>%
  mutate(funYear = as.factor(paste(functionalGroup, Year, sep = "_"))) %>%
  select(siteID, blockID, ID, c(turfID:Year), totalVascular:mossHeight, soil, functionalGroup:precip7010, richness:cwvN) %>% 
  distinct(ID, functionalGroup, TTtreat, .keep_all = TRUE)

# create forb-only data frame
forbcom <- wholecom %>%
  filter(!functionalGroup == "graminoid") 

## ---- Community.mean.weighting.end ---- 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################## WHOLECOM:: CALCULATING TRAITS DELTA ###################

## ---- Explanatory.variable.deltas ---- 

deltacalc <- sapply(1:nrow(forbcom[forbcom$TTtreat == "RTC",]), function(i){
  R <- forbcom[forbcom$TTtreat == "RTC",][i,]
  #browser()
  cols <- c("sumcover", "diversity", "richness", "evenness", "wmeanLDMC", "wmeanSLA", "wmeanLTH", "wmeanLA", "wmeanheight", "wmeanCN", "wmeanN", "cwvLDMC", "cwvSLA", "cwvLTH", "cwvLA", "cwvheight", "cwvN")
  friend <- forbcom$Year == R$Year & forbcom$blockID == R$blockID & forbcom$functionalGroup == R$functionalGroup & forbcom$TTtreat == "TTC"
  if(all (!friend)) {print(R$turfID)
    return(rep(NA, length(cols)))}
  stopifnot(sum(friend) == 1)
  
  f <- forbcom[friend,]
  x <- R[,cols] - f[,cols] # create difference between removal and control
  unlist(x)
})
deltacalc <- as.data.frame(t(deltacalc))
colnames(deltacalc) <- paste0("delta", colnames(deltacalc))

# bind to forb-only wmean data
rtcmeta <- cbind((forbcom[forbcom$TTtreat == "RTC",]), deltacalc)

######### CALCULATING TIME DELTA #######
timedeltacalc <- sapply(1:nrow(forbcom[forbcom$Year != 2011,]), function(i){
 R <- forbcom[forbcom$Year != 2011,][i,]
   #browser()
  cols <- c("sumcover", "diversity", "richness", "evenness", "wmeanLDMC", "wmeanSLA", "wmeanLTH", "wmeanLA", "wmeanheight", "wmeanCN", "cwvLDMC", "cwvSLA", "cwvLTH", "cwvLA", "cwvheight", "cwvCN")
   friend <- forbcom$turfID == R$turfID & forbcom$Year == 2011
  if(all (!friend)) {print(R$turfID)
     return(rep(NA, length(cols)))}
   stopifnot(sum(friend) == 1)
   
   f <- forbcom[friend,]
   x <- R[,cols] - f[,cols]
   unlist(x)
 })
timedeltacalc <- as.data.frame(t(timedeltacalc))
colnames(timedeltacalc) <- paste0("delta", colnames(timedeltacalc))
timedelta <- cbind((forbcom[forbcom$Year != 2011,]), timedeltacalc)



# save files to data folders #
##############################

#save(my.GR.data.FERT, file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/gramRemFert_dataDoc_FJ_SLO.RData")
#save(rtcmeta, file = "~/Documents/FunCaB/bioticInteractions/composition/data/rtcmeta.RData")
#save(wholecom, file = "~/Documents/FunCaB/bioticInteractions/composition/data/wholecom.RData")