#### community composition data ####
load(file = "~/Volumes/fja062/PhD/Data/funcab/funcabCompDataWithTTCs.RData")

## ---- Traits.data.import ---- 

# source Ragnhild's trait data
source("~/Documents/seedclimComm/seedclimComm/ragnhild_trait_data/load_traits.R") # warning here is fine, it just means those spp didn't have CN data collected

traitdata <- traitdata %>% 
  select(siteID, species, Height_mean, LA_mean)

#load from data base
con <- DBI::dbConnect(RMySQL::MySQL(), group = "seedclim")
traits <- DBI::dbGetQuery(con, paste('SELECT taxon.* , Lower, Nem, BNem, SBor, MBor, NBor, LAlp, MAlp, HAlp, Upper, Min_height, Max_height
                                FROM taxon LEFT JOIN moreTraits ON taxon.species = moreTraits.species
                                ORDER BY taxon.species;'))

traits <- traits[traits$species %in% comp2$species,]
traits$HAlp <- as.numeric(traits$HAlp)
traits[,(12:19)][is.na(traits[,12:19])] <- 0 #change this to something more interpretable!
traits <- traits %>% 
  rowwise() %>% 
  mutate(abundance = sum(Nem, BNem, SBor, MBor, NBor, LAlp, MAlp, HAlp, na.rm = TRUE)) %>%
  select(species, functionalGroup)

head(traits)

identical(as.character(traits$species), comp2$species) #this should be identical, but if it's not it just means we are lacking the trait information for some species
identical(as.character(traitdata$species), comp2$species) #this should be identical


# adding traits to my.GR.data
composition <- comp2 %>%
  left_join(traitdata, by = c("species", "siteID")) %>%
  mutate(functionalGroup = if_else(is.na(functionalGroup), "forb", functionalGroup))

#comment these out depending on the analysis you want to run
composition$functionalGroup <- plyr::mapvalues(composition$functionalGroup, from = "pteridophyte", to = "forb")
composition$functionalGroup <- plyr::mapvalues(composition$functionalGroup, from = "woody", to = "forb")

#Species richness
library(vegan)

comp2 <- comp2 %>%
  group_by(turfID, Year, functionalGroup) %>%
  mutate(richness = sum(n_distinct(species))) %>% 
  mutate(diversity = diversity(cover, index = "shannon")) %>% 
  mutate(evenness = (diversity/log(richness)))

composition <- composition %>%
  select(-mossHeight) %>%
  left_join(mossDepth, by = c("turfID", "siteID", "Treatment")) %>% 
  group_by(turfID, siteID, functionalGroup) %>% 
  filter(!is.na(cover)) %>%
  mutate(wmH= weighted.mean(Height_mean, cover, na.rm=TRUE)) %>% #, wmLA= weighted.mean(LA_mean, cover, na.rm=TRUE)
  select(-Height_mean, -LA_mean, -species, -cover, -TotalGraminoids, -totalForbs, -totalBryophytes) %>% 
  distinct(turfID, functionalGroup, .keep_all = TRUE) %>% 
  spread(key =functionalGroup, value = wmH)


## ---- Traits.data.end ---- 
