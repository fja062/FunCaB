# biomass figures

biomass %>% ggplot(aes(x = Year, y = Biomass, colour = Treatment, linetype = Temperature)) +
  geom_line(size = 0.8) +
  facet_grid(precipLevel ~ Removed_FG) +
  theme_classic() +
  scale_color_viridis_d() +
  axis.dim +
  labs(y = "Mean biomass removed (g)")

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/supfig7.jpg", dpi = 300, height = 9, width = 8)


funCov <- comp2 %>% filter(Treatment == "C") %>% 
  left_join(weather) %>%
  group_by(tempLevel, precipLevel, Year) %>% 
  summarise(mossCov = mean(mossCov, na.rm = TRUE),
            graminoidCov = mean(graminoidCov, na.rm = TRUE),
            forbCov = mean(forbCov, na.rm = TRUE)
  )

funCov %>% 
  gather(mossCov, graminoidCov, forbCov, key = FG, value = cov) %>% ggplot(aes(x = Year, y = cov, colour = FG)) + 
  stat_summary(fun.data = "mean_cl_boot") + 
  stat_summary(fun.data = "mean_cl_boot", geom = "line") + 
  facet_grid(.~precipLevel)

funCov %>% 
  gather(mossCov, graminoidCov, forbCov, key = FG, value = cov) %>% ggplot(aes(x = Year, y = cov, colour = FG)) + 
  stat_summary(fun.data = "mean_cl_boot") + 
  stat_summary(fun.data = "mean_cl_boot", geom = "line") + 
  facet_grid(.~tempLevel)
