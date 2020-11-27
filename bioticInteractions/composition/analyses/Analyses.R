# load packages and source plotting code
library("tidyverse")
library("lme4")
library("broom.mixed")
library("lmerTest")
library(sjPlot)
library(yhat)
require(car)
require(GGally)
source("~/Documents/FunCaB/figures/plotting_dim.R")

# load data
load("~/Documents/FunCaB/bioticInteractions/composition/data/rtcmeta.RData")
load("~/Documents/FunCaB/bioticInteractions/composition/data/wholecom.RData")

# prepare data
forbcomAnalysis <- rtcmeta %>% 
  gather(c(deltasumcover:deltawmeanCN, deltacwvLDMC:deltacwvheight), key = "trait", value = "value") %>%
  select(-c(richness:cwvN), - deltawmeanN, -deltacwvN) %>% 
  tibble()

# Scaling explanatory variables
forbcomAnalysis <- forbcomAnalysis %>% 
  mutate(tempLevel = as.factor(tempLevel),
         precipLevel = as.factor(precipLevel),
         Sprecip7010 = as.factor(precip7010/100),
         Year = as.numeric(as.character(Year))-2000) %>% 
  group_by(trait) %>% 
  filter(is.finite(value)) %>% 
  ungroup()

# test for collinearity
forbcom %>% 
  filter(siteID %in% c("Rambera", "Ovstedal", "Gudmedalen")) %>% 
  select(c(wmeanheight, wmeanLTH, wmeanLDMC, wmeanSLA, tempLevel)) %>% 
  GGally::ggpairs(mapping = aes(colour = factor(tempLevel)))


modCov <- lmer(value ~ tAnom*pAnom + (1|siteID/turfID), REML = TRUE, data = subset(forbcomAnalysis, trait == "deltasumcover"))
  
# sum cover
modCov <- lmer(value ~ tempLevel*Year*precipLevel - tempLevel:precipLevel:Year - Year:precipLevel - tempLevel:precipLevel + (1|siteID/turfID), REML = TRUE, data = subset(forbcomAnalysis, trait == "deltasumcover"))

anova(modCov)
summary(modCov)

modCov <- forbcomAnalysis %>% 
  filter(trait == "deltasumcover") %>% 
  nest(data = everything()) %>%
  mutate(fit = map(data, ~ lmer(value ~ tempLevel*Year*precipLevel - tempLevel:precipLevel:Year - Year:precipLevel - tempLevel:precipLevel + (1|siteID/turfID), REML = TRUE, data = .)),
         results = map(fit, tidy, conf.int = TRUE)) %>%
  unnest(results) %>% 
  select(-data, -fit)

modCov
plot_residuals(modCov)
plot_model(modCov, type = "re")
plot_model(modCov)

forbcomAnalysis %>% 
  filter(trait == "deltasumcover") %>% 
  ggplot(aes(x = Year, y = value, colour = factor(tempLevel), fill = factor(tempLevel))) +
  scale_colour_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") +
  stat_summary(fun.data = "mean_se", geom = "point") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  geom_hline(yintercept = 0, colour = "grey50") +
  facet_wrap(~tempLevel) +
  theme_classic() +
  theme(legend.position = "none")

# richness

modRich <- lmerTest::lmer(value ~ tempLevel*Year*precipLevel - tempLevel:precipLevel:Year - tempLevel:precipLevel - precipLevel:Year - tempLevel:Year + (1|siteID/turfID), REML = TRUE, data = subset(forbcomAnalysis, trait == "deltarichness"))

anova(modRich)
summary(modRich)

modRich <- forbcomAnalysis %>% 
  filter(trait == "deltarichness") %>% 
  nest(data = everything()) %>%
    mutate(fit = map(data, ~ lmer(value ~ tempLevel + precipLevel + Year + (1|siteID/turfID) + (1|Year), REML = TRUE, data = .)),
           results = map(fit, tidy, conf.int = TRUE)) %>%
    unnest(results) %>% 
    select(-data, -fit)

plot_model(modRich, show.values = TRUE)
plot_model(modRich, type = "re")
plot_model(modRich, type = "pred", pred.type = "re")
plot_residuals(modRich)

# evenness
modEv <- forbcomAnalysis %>% 
  filter(trait == "deltaevenness") %>% 
    nest(data = everything()) %>%
    mutate(fit = map(data, ~ lmer(value ~ temp7010 + Sprecip7010 + (1|siteID/turfID) + (1|Year), REML = TRUE, data = .)),
           results = map(fit, tidy, conf.int = TRUE)) %>%
    unnest(results) %>% 
    select(-data, -fit)
  
modEv <- lmerTest::lmer(value ~ tempLevel*Year*precipLevel - tempLevel:precipLevel:Year - tempLevel:precipLevel + (1|siteID/turfID), REML = TRUE, data = subset(forbcomAnalysis, trait == "deltaevenness"))

anova(modEv)
summary(modEv)

plot_model(modEv, show.values = TRUE)
plot_model(modEv, type = "re")
plot_model(modEv, type = "pred")
plot_residuals(modEv)


# SLA
modSLA <- forbcomAnalysis %>% 
  filter(trait == "deltawmeanSLA") %>% 
  nest(data = everything()) %>%
    mutate(fit = map(data, ~ lmer(value ~ precipLevel*Year + (1|siteID/turfID), REML = TRUE, data = .)),
           results = map(fit, tidy, conf.int = TRUE)) %>%
    unnest(results) %>% 
    select(-data, -fit)


modSLA <- lmerTest::lmer(value ~ tempLevel*precipLevel*Year - tempLevel:precipLevel:Year - tempLevel:precipLevel - tempLevel:Year + (1|siteID/turfID), REML = TRUE, data = subset(forbcomAnalysis, trait == "deltawmeanSLA"))

anova(modSLA)
summary(modSLA)

plot_model(modSLA, show.values = TRUE)
plot_model(modSLA, type = "re")
plot_model(modSLA, type = "pred")
plot_residuals(modSLA)

forbcomAnalysis %>% 
  filter(trait == "deltawmeanSLA") %>% 
  ggplot(aes(x = Year, y = value, colour = factor(precipLevel), fill = factor(precipLevel))) +
  scale_colour_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") +
  stat_summary(fun.data = "mean_se", geom = "point") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  geom_hline(yintercept = 0, colour = "grey50") +
  facet_wrap(~precipLevel, scales = "free") +
  theme_classic() +
  theme(legend.position = "none")

# LTH
modLTH <- forbcomAnalysis %>% 
  filter(trait == "deltawmeanLTH") %>% 
  nest(data = everything()) %>%
  mutate(fit = map(data, ~ lmer(value ~ temp7010+Sprecip7010+Year + (1|siteID), REML = TRUE, data = .)),
         results = map(fit, tidy, conf.int = TRUE)) %>%
  unnest(results) %>%
  select(-data, -fit)
modLTH

forbcomAnalysis %>% 
  filter(trait == "deltawmeanLTH") %>% 
  ggplot(aes(x = Year, y = value)) +
  scale_colour_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") +
  stat_summary(fun.data = "mean_se", geom = "point") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  geom_hline(yintercept = 0, colour = "grey50") +
  theme_classic()


modLTH <- lmerTest::lmer(value ~ tempLevel*precipLevel*Year - tempLevel:precipLevel:Year - tempLevel:precipLevel - tempLevel:Year - precipLevel:Year + (1|siteID), REML = TRUE, data = subset(forbcomAnalysis, trait == "deltawmeanLTH"))

anova(modLTH)
summary(modLTH)

plot_model(modLTH, show.values = TRUE)
plot_model(modLTH, type = "re")
plot_model(modLTH, type = "pred")
plot_residuals(modLTH)


# height
modheight <- forbcomAnalysis %>% 
  filter(trait == "deltawmeanheight") %>% 
  nest(data = everything()) %>%
  mutate(fit = map(data, ~ lmer(value ~ temp7010 + Sprecip7010 + Year + (1|siteID), REML = TRUE, data = .)),
         results = map(fit, tidy, conf.int = TRUE)) %>%
  unnest(results) %>%
  select(-data, -fit)
modheight

forbcomAnalysis %>% 
  filter(trait == "deltawmeanheight") %>% 
  ggplot(aes(x = Year, y = value, colour = tempLevel)) +
  scale_colour_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") +
  stat_summary(fun.data = "mean_se", geom = "point") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  geom_hline(yintercept = 0, colour = "grey50") +
  theme_classic() +
  facet_wrap(~tempLevel)


modheight <- lmerTest::lmer(value ~ tempLevel*precipLevel*Year - tempLevel:precipLevel:Year - precipLevel:Year - tempLevel:Year - tempLevel:precipLevel + (1|siteID), REML = TRUE, data = subset(forbcomAnalysis, trait == "deltawmeanheight"))

anova(modheight)
summary(modheight)

plot_model(modheight, show.values = TRUE)
plot_model(modheight, type = "re")
plot_model(modheight, type = "pred")
plot_residuals(modheight)

# LDMC
modLDMC <- forbcomAnalysis %>% 
  filter(trait == "deltawmeanLDMC") %>% 
  nest(data = everything()) %>%
  mutate(fit = map(data, ~ lmer(value ~ temp7010 + Sprecip7010 + Year + (1|siteID), REML = TRUE, data = .)),
         results = map(fit, tidy, conf.int = TRUE)) %>%
  unnest(results) %>%
  select(-data, -fit)

forbcomAnalysis %>% 
  filter(trait == "deltawmeanLDMC") %>% 
  ggplot(aes(x = Year, y = value, colour = tempLevel)) +
  scale_colour_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") +
  stat_summary(fun.data = "mean_se", geom = "point") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  geom_hline(yintercept = 0, colour = "grey50") +
  theme_classic() +
  facet_wrap(~tempLevel)


modLDMC <- lmerTest::lmer(value ~ tempLevel*precipLevel*Year - tempLevel:precipLevel:Year - tempLevel:precipLevel - tempLevel:Year - precipLevel:Year + (1|siteID), REML = TRUE, data = subset(forbcomAnalysis, trait == "deltawmeanLDMC"))

anova(modLDMC)
summary(modLDMC)

plot_model(modLDMC, show.values = TRUE)
plot_model(modLDMC, type = "re")
plot_model(modLDMC, type = "pred")
plot_residuals(modLDMC)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# minimum interesting effect size => variation across years in controls
sim_lm <- function(n, delta, sd, ...){
  # simulate means
  mu <- rep(c(0, delta), each = n)
  # add noise
  y <- mu + rnorm(length(mu), sd = sd)
  ## predictor
  test <- crossing(year = c(1,2,3,4,5),
                   temp = c(1,2,3,4,5,6,7,8,9,10,11,12),
                   precip = c(11,12,13,14,15,16,17,18,19,20,21,22)
  )  
  # run test
  test <- lm(y ~ temp*precip*year, data = test)
  broom::glance(test)
}

sim_lm(n = 360, delta = 1, sd = 2)

nrep = 100

control <- crossing(rep_no = 1:nrep, n = seq(10, 100, 20)) 

runs <- control %>%
  pmap_df(~sim_lm(n = 360, delta = 1, sd = 2, nrep = nrep))  %>%
  mutate(sig = p.value <= 0.05)

p <- runs  %>%
  group_by(term) %>%
  summarise(power = mean(sig)) %>%
  ggplot(aes(x = n, y = power)) +
  geom_line() +
  geom_point()

runs %>%
  filter(sig) %>%
  ggplot(aes(x = n, y = estimate, group = n)) +
  geom_hline(yintercept = -1, colour = "red") +
  geom_violin(draw_quantiles = 0.5, fill = "grey50", alpha = 0.6)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mod1temp <- forbcomAnalysis %>% 
  mutate(traitII = trait) %>%
  group_by(trait) %>%
  do({if(.$trait[1] == "richness"){
    # richness
    mod <- glmer(value ~ TTtreat*Stemp0916*Sprecip0916*SYear - TTtreat:Stemp0916:Sprecip0916:SYear + (1|siteID/blockID), family = "poisson", data = .)
  } else if (.$trait[1] %in% c("wmeanLTH", "wmeanLA", "wmeanheight", "cwvSLA", "cwvLTH", "cwvLDMC", "cwvCN", "cwvheight", "cwvLA")){
    # leaf thickness, leaf area, height, cwv
    mod <- lmer(log(value) ~ TTtreat*Stemp0916*Sprecip0916*SYear - TTtreat:Stemp0916:Sprecip0916:SYear + (1|siteID/blockID), REML = FALSE, data = .)
  } else if (.$trait[1] %in% c("evenness", "diversity")){
    mod <- lmer(value^2 ~ TTtreat*Stemp0916*Sprecip0916*SYear - TTtreat:Stemp0916:Sprecip0916:SYear + (1|siteID/blockID), REML = FALSE, data = .)
  } else {
    mod <- lmer(value ~ TTtreat*Stemp0916*Sprecip0916*SYear - TTtreat:Stemp0916:Sprecip0916:SYear + (1|siteID/blockID), REML = FALSE, data = .) 
  }
    tidy(mod)}) %>% 
  arrange(desc(trait)) %>% 
  mutate(lower = (estimate - std.error*1.96),
         upper = (estimate + std.error*1.96)) %>% 
  ungroup()


mod1temp <- mod1temp %>% 
  mutate(test = case_when(
    grepl("wmean", trait) ~ "Mean",
    grepl("cwv", trait) ~ "Variance",
    grepl("^s|^r|^e|^d", trait) ~ "Mean"),
    term = case_when(
      term == "(Intercept)" ~ "Control",
      term == "Stemp0916" ~ "t",
      term == "Sprecip0916" ~ "P",
      term == "SYear" ~ "year",
      term =="TTtreatRTC:Stemp0916" ~ "t x removal",
      term =="TTtreatRTC:Sprecip0916" ~ "P x removal",
      term =="TTtreatRTC:Stemp0916:SYear" ~ "t x year x removal",
      term =="TTtreatRTC:Sprecip0916:SYear" ~ "P x year x removal",
      term =="TTtreatRTC:Stemp0916:Sprecip0916" ~ "P x t x removal",
      term =="Stemp0916:SYear" ~ "t x year",
      term =="Sprecip0916:SYear" ~ "P x year",
      term =="Stemp0916:Sprecip0916" ~ "P x t",
      term =="Stemp0916:Sprecip0916:SYear" ~ "P x t x year",
      term =="TTtreatRTC:SYear" ~ "Year x removal",
      term == "TTtreatRTC" ~ "removal")) %>% 
  mutate(trait = if_else(grepl("wmean", trait), substr(trait, 6, n()),
                         if_else(grepl("cwv", trait), substr(trait, 4, n()), trait))) %>% 
  mutate(sign = recode(trait, sumcover = 1, evenness = 1, richness = 1, diversity = 1, height = 0, LA = 0, LTH = 0, LDMC = 0, CN = 1, SLA = 1))

mod1temp %>% 
  select(-(p.value:upper), -sign) %>% 
  #filter(test == "Mean") %>% 
  #filter(term %in% c("Control", "removal", "t", "P", "year", "P x t x year", "t x year x removal", "P x year x removal"), test == "Mean") %>% 
  filter(!is.na(term)) %>% 
  mutate(estimate = round(estimate, 3), std.error = round(std.error, 3), statistic = round(statistic, 3)) %>% 
  write_csv(., path = "~/OneDrive - University of Bergen/Research/FunCaB/paper 1/results/mod1tempOUT1.csv")


coefEst <- mod1temp %>%
  filter(term %in% c("P x year x removal", "t x year x removal", "Year x removal", "removal")) %>% 
  ggplot(aes(x = trait, y = estimate, ymin = lower, ymax = upper, fill = factor(term, levels = c("P x year x removal", "t x year x removal", "Year x removal", "removal")), shape = factor(term, levels = c("P x year x removal", "t x year x removal", "Year x removal", "removal")), alpha = as.factor(sign))) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.5), aes(colour = factor(term, levels = c("P x year x removal", "t x year x removal", "Year x removal", "removal")))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2.9) +
  coord_flip() +
  geom_vline(xintercept =  c(1.5,2.5,3.5,5.5,6.5,8.5,9.5), colour = "grey90") +
  geom_vline(xintercept =  7.5, colour = "black") +
  geom_vline(xintercept =  4.5, colour = "grey50") +
  scale_alpha_manual(values = c(0.6, 1), guide = FALSE) +
  scale_fill_manual(legend.title.climate, values = c("#1C9099", "#E69F00", "grey90", "black")) +
  scale_colour_manual(legend.title.climate, values = c("black", "black", "black", "black")) +
  scale_shape_manual(legend.title.climate, values = c(25, 24, 23, 21)) +
  #scale_linetype_manual(legend.title.climate, values = c(1,1,3, 21,21,23, 25)) +
  scale_x_discrete(limits = c("SLA", "CN", "LDMC", "LTH", "LA", "height", "diversity", "richness", "evenness", "sumcover"), labels = c("SLA", "C:N ratio", "Leaf dry \n matter content", "Leaf thickness", "Leaf area", "Height", "diversity", "Richness", "Evenness", "Cover")) +
  facet_wrap(~test, strip.position = "top", scales = "free_x") +
  labs(y = "Standardised coefficients", x = "Leaf economic traits                 Structural traits               Community structure") +
  theme_cowplot(font_family = "Helvetica") +
  #ylim(c(-0.4, 0.6)) +
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom",
        legend.justification = "centre",
        legend.background = element_rect(fill = "white"),
        strip.text = element_text(size = 14, hjust = 0.42),
        axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_text(colour = c("black", "black", "grey40", "grey40", "grey40", "grey40", "black", "black", "black", "black")))
  
coefEst <- plot_grid(coefEst, labels = c("B                                                       C"), label_x = 0)



