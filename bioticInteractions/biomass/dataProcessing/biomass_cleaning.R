#### #### #### #### #### 
#### biomass #### 
# load packages
library(broom)

# read in biomass files
bio15 <- read_excel("~/OneDrive - University of Bergen/Research/FunCaB/Data/primary/veg_biomass/FunCaB_biomass_2015-2018.xlsx", sheet = 1) %>% select(-Date)
bio16 <- read_excel("~/OneDrive - University of Bergen/Research/FunCaB/Data/primary/veg_biomass/FunCaB_biomass_2015-2018.xlsx", sheet = 2) %>% select(-Date) %>% mutate(Biomass = as.numeric(Biomass))
bio17 <- read_excel("~/OneDrive - University of Bergen/Research/FunCaB/Data/primary/veg_biomass/FunCaB_biomass_2015-2018.xlsx", sheet = 3) %>% select(-Date) %>% mutate(Biomass = as.numeric(Biomass))
bio18 <- read_excel("~/OneDrive - University of Bergen/Research/FunCaB/Data/primary/veg_biomass/FunCaB_biomass_2015-2018.xlsx", sheet = 4) %>% select(-Date)

biomass <- bind_rows(bio15, bio16, bio17, bio18) %>% 
  filter(Round == 1,
         !Treatment == "RTC") %>% 
  left_join(dict_Site, by = c("Site" = "v2")) %>% 
  mutate(siteID = new) %>% 
  select(-Site, -old, -v3, -new) %>% 
  group_by(Year, siteID, Block, Treatment, Removed_FG) %>% 
  summarise(Biomass = mean(Biomass, na.rm = TRUE)) %>% 
  left_join(weather) %>% 
  ungroup() %>% 
  mutate(Year = as.numeric(substr(Year, 3,4)),
         Temperature = factor(tempLevel, levels = c(6.5, 8.5, 10.5), labels = c("Alpine", "Sub-alpine", "Boreal"))) %>%
  filter(!Treatment %in% c("RTC", "RTC2nd"))

# vegetation data
comp2 <- comp2 %>% 
  filter(!is.na(Treatment), !Treatment == "XC")

# FIX THIS!!!
# set covers and heights to NA in removed plots
comp201718 <- comp2 %>% 
  filter(Year > 2015) %>%
  mutate(mossCov = case_when(
    grepl("B", Treatment) ~ NA_real_,
    TRUE ~ mossCov),
    forbCov = case_when(
      grepl("F", Treatment) ~ NA_real_, 
      TRUE ~ forbCov),
    graminoidCov = case_when(
      grepl("G", Treatment) ~ NA_real_, 
      TRUE ~ graminoidCov),
    vegetationHeight = case_when(
      Treatment == "FGB" ~ NA_real_, 
      TRUE ~ vegetationHeight),
    mossHeight = case_when(
      Treatment == "FGB" ~ NA_real_, 
      TRUE ~ mossHeight))

composition2015 <- comp2 %>% 
  filter(Year == 2015) %>% 
  left_join(mossHeight, by = "turfID", suffix = c("", ".new")) %>% 
  mutate(mossHeight = if_else(is.na(mossHeight), mossHeight.new, mossHeight),
         blockID = as.character(blockID)) %>%
  distinct(siteID, Treatment, turfID, blockID, vegetationHeight, mossHeight, litter, mossCov, forbCov, graminoidCov) %>% 
  ungroup()

composition2015 <- composition2015 %>%
  group_by(siteID, blockID) %>% 
  mutate(mossHeight = mean(mossHeight, na.rm = TRUE))

# gather vegetation into columns for join with biomass
biomassComp <- composition2015 %>% 
  gather(graminoidCov, forbCov, mossCov, key = "functionalGroup", value = "covValue") %>% 
  mutate(functionalGroup = recode(functionalGroup, "mossCov" = "bryophyte", "graminoidCov" = "graminoid", "forbCov" = "forb"))
  
# align naming conventions among vegetation and biomass datasets
# join biomass data to composition data
biomassReg <- biomass %>% 
  mutate(blockID = as.character(if_else(siteID == "Gudmedalen", recode(Block, "1" = 5, "2" = 12, "3" = 13, "4" = 15), Block)),
         Treatment = if_else(Treatment == "FG", "GF", Treatment),
         turfID = paste0(substr(siteID, 1, 3), Block, Treatment),
         turfID = if_else(siteID == "Gudmedalen", paste0(substr(siteID, 1, 3), blockID, Treatment), turfID),
         turfID = recode(turfID, "Alr4FGB" = "Alr5C"),
         functionalGroup = recode(Removed_FG, "B" = "bryophyte", "G" = "graminoid", "F" = "forb")) %>% 
  left_join(biomassComp)

# turn biomass from g to g/m^2
biomassReg <- biomassReg %>% 
  mutate(Biomass_gm = Biomass/0.0625) 

# run and extract regressions for each Functional group with zero intercept
# forb
forbRegCoef <- biomassReg %>% 
  filter(functionalGroup == "forb") %>% 
  lm(Biomass_gm ~ 0 + covValue:vegetationHeight, data = .)
forbRegCoef <- tidy(forbRegCoef)

# graminoid
gramRegCoef <- biomassReg %>% 
  filter(functionalGroup == "graminoid") %>% 
  lm(Biomass_gm ~ 0 + covValue:vegetationHeight, data = .)
gramRegCoef <- tidy(gramRegCoef)

# moss
mossRegCoef <- biomassReg %>% 
  filter(functionalGroup == "bryophyte") %>% 
  lm(Biomass_gm ~ 0 + covValue:mossHeight, data = .)
mossRegCoef <- tidy(mossRegCoef)

# bind model estimates together
regCoef <- bind_rows("forb" = forbRegCoef, "graminoid" = gramRegCoef, "bryophyte" = mossRegCoef, .id = "functionalGroup") %>% 
  select(estimate, functionalGroup) %>% 
  spread(estimate, key = functionalGroup)


biomassReg <- biomassReg %>% 
  select(-functionalGroup, - Biomass, -Biomass_gm, -litter) %>% 
  full_join(composition2015 %>% select(-litter)) %>% 
  group_by(siteID, blockID, turfID, Treatment) %>% 
  mutate(forbBiomass = forbCov*vegetationHeight*regCoef$forb,
         graminoidBiomass = graminoidCov*vegetationHeight*regCoef$graminoid,
         mossBiomass = mossCov*mossHeight*regCoef$bryophyte) %>% 
  distinct() %>% 
  select(blockID, Treatment, turfID, siteID, forbBiomass, graminoidBiomass, mossBiomass, forbCov, graminoidCov, mossCov, vegetationHeight, mossHeight)
