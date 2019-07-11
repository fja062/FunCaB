#### #### #### #### #### 
#### biomass #### 
# load packages
library(broom)

# read in biomass files
bio15 <- read_excel("~/OneDrive - University of Bergen/Research/FunCaB/Data/primary/veg_biomass/biomass_removals_2015.xlsx")

biomass <- bio15 %>%
  filter(!Treatment %in% c("RTC", "RTC2nd")) %>% 
  left_join(dict_Site, by = c("Site" = "old")) %>% 
  select(-c(Site, Round, Year, v2, v3), "siteID" = new, "blockID" = Block, "turfID" = TurfID, "functionalGroup" = Func_group)

# vegetation data
comp2 <- comp2 %>% 
  filter(!is.na(Treatment), !Treatment == "XC")

# FIX THIS!!!
# set covers and heights to zero in removed plots
comp2 <- comp2 %>% 
  mutate(mossCov = if_else(grepl("B", Treatment) & Year > 2015, 0, mossCov),
         forbCov = if_else(grepl("F", Treatment) & Year > 2015, 0, forbCov),
         graminoidCov = if_else(grepl("G", Treatment) & Year > 2015, 0, graminoidCov),
         vegetationHeight = if_else(Treatment == "FGB" & Year > 2015, 0, vegetationHeight),
         mossHeight = if_else(Treatment == "FGB" & Year > 2015, 0, mossHeight))

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
  mutate(blockID = as.character(if_else(siteID == "Gudmedalen", recode(blockID, "1" = 5, "2" = 12, "3" = 13, "4" = 15), blockID)),
         turfID = if_else(siteID == "Gudmedalen", paste0(substr(siteID, 1, 3), blockID, Treatment), turfID),
         turfID = recode(turfID, "Alr4FGB" = "Alr5C"),
         functionalGroup = recode(functionalGroup, "B" = "bryophyte", "G" = "graminoid", "F" = "forb")) %>% 
  left_join(composition2015)

# turn biomass from g to g/m^2
biomassReg <- biomassReg %>% 
  mutate(Biomass_gm = Biomass_g/0.0625) 

# run and extract regressions for each Functional group with zero intercept
# forb
forbRegCoef <- biomassReg %>% 
  filter(functionalGroup == "forb") %>% 
  lm(Biomass_gm ~ 0 + forbCov:vegetationHeight, data = .)
forbRegCoef <- tidy(forbRegCoef)

# graminoid
gramRegCoef <- biomassReg %>% 
  filter(functionalGroup == "graminoid") %>% 
  lm(Biomass_gm ~ 0 + graminoidCov:vegetationHeight, data = .)
gramRegCoef <- tidy(gramRegCoef)

# moss
mossRegCoef <- biomassReg %>% 
  filter(functionalGroup == "bryophyte") %>% 
  lm(Biomass_gm ~ 0 + mossCov:mossHeight, data = .)
mossRegCoef <- tidy(mossRegCoef)

# bind model estimates together
regCoef <- bind_rows("forb" = forbRegCoef, "graminoid" = gramRegCoef, "bryophyte" = mossRegCoef, .id = "functionalGroup") %>% 
  select(estimate, functionalGroup) %>% 
  spread(estimate, key = functionalGroup)


biomassReg <- biomassReg %>% 
  select(-functionalGroup, - Biomass_g, -Biomass_gm, -litter) %>% 
  full_join(composition2015 %>% select(-litter)) %>% 
  group_by(siteID, blockID, turfID, Treatment) %>% 
  mutate(forbBiomass = forbCov*vegetationHeight*regCoef$forb,
         graminoidBiomass = graminoidCov*vegetationHeight*regCoef$graminoid,
         mossBiomass = mossCov*mossHeight*regCoef$bryophyte) %>% 
  distinct() %>% 
  select(blockID, Treatment, turfID, siteID, forbBiomass, graminoidBiomass, mossBiomass, forbCov, graminoidCov, mossCov, vegetationHeight, mossHeight)

