# compile CWM an Fdiv trait values
# load imputation files for each traits for all species

trait_C <- read.csv("~/OneDrive - University of Bergen/Research/FunCaB/Data/imputation_files/speciesSitePredictions_C.csv", header =TRUE, sep = ";", dec = ",") %>%
  mutate(Species = paste0(Site, "_", Species)) %>%
  select( Site, Species, predictValue, predictionSE) #, predictionSE
trait_N <- read.csv("~/OneDrive - University of Bergen/Research/FunCaB/Data/imputation_files/speciesSitePredictions_N.csv", header =TRUE, sep = ";", dec = ",") %>%
  select(predictValue, predictionSE)
trait_CN <- read.csv("~/OneDrive - University of Bergen/Research/FunCaB/Data/imputation_files/speciesSitePredictions_CN.ratio.csv", header =TRUE, sep = ";", dec = ",") %>%
  select(predictValue, predictionSE)
trait_SLA <- read.csv("~/OneDrive - University of Bergen/Research/FunCaB/Data/imputation_files/speciesSitePredictions_SLA.csv", header =TRUE, sep = ";", dec = ",") %>%
  select(predictValue, predictionSE)
trait_Lth <- read.csv("~/OneDrive - University of Bergen/Research/FunCaB/Data/imputation_files/speciesSitePredictions_Lth_ave.csv", header =TRUE, sep = ";", dec = ",") %>%
  select(predictValue, predictionSE)
trait_LDMC <- read.csv("~/OneDrive - University of Bergen/Research/FunCaB/Data/imputation_files/speciesSitePredictions_LDMC.csv", header =TRUE, sep = ";", dec = ",") %>%
  select(predictValue, predictionSE)
trait_logLA <- read.csv("~/OneDrive - University of Bergen/Research/FunCaB/Data/imputation_files/speciesSitePredictions_logLA.csv", header =TRUE, sep = ";", dec = ",") %>%
  select(predictValue, predictionSE)
trait_logHeight <- read.csv("~/OneDrive - University of Bergen/Research/FunCaB/Data/imputation_files/speciesSitePredictions_logHeight.csv", header =TRUE, sep = ";", dec = ",") %>%
  select(predictValue, predictionSE)

Species_traits <-bind_cols(trait_C, trait_N, trait_CN, trait_SLA,
                           trait_Lth, trait_LDMC, trait_logLA, trait_logHeight)%>%
  rename(C = predictValue, N = predictValue1, CN = predictValue2, SLA =
           predictValue3, Lth = predictValue4, LDMC = predictValue5, LA =
           predictValue6, Height = predictValue7) %>% 
  mutate(sqrtLA = sqrt(exp(LA))) %>% 
  rename(logHeight = Height)

## check distribution of imputation values of traits
ggplot(Species_traits, aes(logHeight)) +
  geom_density() +
  facet_wrap(~Site)
## Distribution of imputed values within range of measured traits

## Checked SE values for each trait- Emp_nig, Emp_her, Ver_alp and Sil_vul very large SE remove from dataset
Species_traits <- Species_traits %>%
  filter(!grepl("Emp_",  Species)) %>%
  filter(!grepl("Ver_alp",  Species)) %>%
  filter(!grepl("Sil_vul",  Species)) %>% 
  mutate(Site = plyr::mapvalues(Site, from = dict_Site$old, to = dict_Site$new)) %>%
  select(siteID = Site, species = Species, C, N, CN, SLA, Lth, LDMC, sqrtLA, logHeight) %>% 
  mutate(species = substr(species, 5, n()),
         species = gsub("_", ".", species),
         speciesID = paste0(siteID, "_", species))
