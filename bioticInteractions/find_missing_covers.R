spp2016 <- my.GR.data %>% 
  filter(cover >=0 & Year == 2016) %>% 
  distinct(turfID, species)

my.GR.data %>% 
  anti_join(spp2016) %>% 
  filter(functionalGroup %in% c("graminoid", "forb")) %>% 
  filter(siteID %in% c("Arhelleren")) %>%
  group_by(species, turfID) %>% 
  ggplot(aes(x = Year, y = cover, colour = species)) +
  geom_point() +
  geom_line() +
  facet_grid(turfID ~ functionalGroup)
