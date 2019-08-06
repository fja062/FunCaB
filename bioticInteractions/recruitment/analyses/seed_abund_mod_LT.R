#### Load packages and data ####
library("rjags")
library("R2jags")
library("tidyverse")
library("tidybayes")
library("DHARMa")
# library("nimble")
#library("lme4")
#library("broom.mixed")

# data
load(file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/cleanedAbundData.RData")

# 1) 
#### data preparation ####
abdat <- rc_rtcSumAv  %>% 
  filter(Treatment %in% c("Intact", "Gap")) %>%
  mutate(stempLevel = as.vector(scale(tempLevel, scale = FALSE, center =TRUE)),
         stemp7010 = as.vector(scale(temp7010, scale = FALSE, center =TRUE)),
         precipDiv7010 = precip7010/1000,
         sprecip7010 = as.vector(scale(precipDiv7010, scale = FALSE, center =TRUE)),
         Treatment = factor(Treatment, levels = c("Intact", "Gap")),
         monthN = factor(monthN, levels = c("spr", "aut")),
         precipDiv = precipLevel/1000,
         sprecipLevel = as.vector(scale(precipDiv, scale = FALSE, center = TRUE)),
         soilTs = as.vector(scale(soilT, scale = FALSE, center = TRUE)),
         soilMs = as.vector(scale(soilM, scale = FALSE, center = TRUE)),
         ) %>%
  filter(!is.na(soilM),
         !is.na(soilT))
#MASS::fitdistr(abdat %>% pull(seed), densfun = "negative binomial")


#### Non-bayesian analysis ####
nbGlmerAb <- abdat %>% 
  glmer.nb(sqrt(seed) ~ monthN + Treatment + soilTs + I(soilTs^2) + soilMs + I(soilMs^2) + soilTs:Treatment + I(soilTs^2):Treatment + I(soilMs^2):Treatment + soilMs:Treatment + (1|siteID), data = ., family = "gamma")

simOut <- simulateResiduals(fittedModel = nbGlmerAb, n = 250)
plot(simOut)
testResiduals(simOut)
testZeroInflation(simOut)

summary(nbGlmerAb)

nbGlmerAb %>%
  tidy()  %>% 
  mutate(lower = (estimate - std.error*1.96),
         upper = (estimate + std.error*1.96))%>% 
  ggplot(aes(x = estimate, y = term, xmin = lower, xmax = upper)) +
  geom_errorbarh() +
  geom_point() +
  geom_vline(xintercept = 0)


#### Bayesian analysis ####
# i) set up a model matrix to feed directly into the model. This avoids potential coding errors. -1 removes the intercept, which I set separately so it can be drawn from a normal distribution.
matab.t <- model.matrix(~ monthN + Treatment + tAnom + pAnom + tAnom:Treatment + pAnom:Treatment + stemp7010 + sprecip7010, data = abdat)[,-1]

# fake data for predictions
abdatY <- crossing(Treatment = unique(abdat$Treatment), # rep is slowest on inside
                   monthN = unique(abdat$monthN), 
                   temp7010 = unique(abdat$tempLevel),
                   precip7010 = quantile(abdat$precipDiv7010, prob = c(0.4, 0.5)),
                   pAnom = seq(min(abdat$pAnom), max(abdat$pAnom), length = 50),
                   tAnom = quantile(abdat$tAnom, prob = c(0.25, 0.75))) %>% 
  mutate(stemp7010 = as.vector(scale(temp7010, scale = FALSE, center =TRUE))
  )
abdatY

# model matrix for fake data predictions
matab.tY <- model.matrix(~ monthN + Treatment + tAnom + pAnom + tAnom:Treatment + pAnom:Treatment + stemp7010 + precip7010, data = abdatY)
# remove intercept
matab.tY <- matab.tY[,-1]

# ii) model
cat("model {
  # Likelihood
  for (i in 1:n.dat) {
    y[i] ~ dnegbin(p[i], r)
    p[i] <- r / (r + mu[i])
    log(mu[i]) <- beta.intercept + inprod(beta, matX[i, ]) + beta.site[siteID[i]]
    
    # predictions for model validation, using original data
    yPred[i] ~ dnegbin(p[i], r)

  }

   # derived predictions
for (k in 1:n.datY){
    #Pred[k] ~ dnegbin(pPred[k], r)          # new data for each MCMC iteration
    pPred[k] <- r / (r + muPred[k])

    log(muPred[k]) <- beta.intercept + inprod(beta, matY[k, ])
    }

  # Priors
  r ~ dgamma(0.01, 0.01)            # prior for the precision of the survival probability
  beta.intercept ~ dnorm(0, 0.001)    # intercept prior

  for (b in 1:nEff) {
    beta[b] ~ dnorm(0, 0.001)         # priors for the betas
  }
  
  # priors random effects
  randTau ~ dgamma(0.01, 0.01)
  for (m in 1:n.site) {
    beta.site[m] ~ dnorm(0, randTau)
  }


}
", fill = TRUE, file = "~/Documents/FunCaB/analyses/seedAbund_tAnom.txt")


# specify the parameters to watch
paraNames.ab <- c("beta.intercept", "beta", "beta.site", "r", "yPred", "muPred", "mu")

# iii) Set up a list that contains all the necessary data
n.treat <- nlevels(factor(abdat$Treatment))
n.season <- nlevels(factor(abdat$monthN))

abDat <- list(y = abdat$seed, 
            n.dat = nrow(abdat),
            n.datY = nrow(abdatY),
            matX = matab.t,
            matY = matab.tY,
            nEff = ncol(matab.t),
            siteID = as.numeric(factor(abdat$siteID)),
            n.site = nlevels(factor(abdat$siteID)))



# iv) Compile the model and run the MCMC for an adaptation/burn-in phase and sample from the posteriors
AbundtAnom.mod <- jags(
  model.file = "~/Documents/FunCaB/analyses/seedAbund_tAnom.txt",
  data = abDat,
  n.iter = 20000,
  n.chains = 4,
  parameters.to.save = paraNames.ab,
  progress.bar = "text"
)


#  vi) diagnostics
# create variables for model checking
simulations.abT <- AbundtAnom.mod$BUGSoutput$sims.list$yPred
predictions.abT <- apply(AbundtAnom.mod$BUGSoutput$sims.list$mu, 2, median)
jagsModabT.meanlist <- AbundtAnom.mod$BUGSoutput$mean
jagsModabT.paramlist <- AbundtAnom.mod$BUGSoutput$sims.list
drawsAB <- AbundtAnom.mod$BUGSoutput$sims.list$pPred %>% as.data.frame()
dim(simulations.abT)

# extract names from matrix
rNames.t <- colnames(matab.t) %>% 
  enframe(name = "i", value = "term") %>% 
  mutate(i = paste0("beta[",i,"]"))



sim.abT <- createDHARMa(
    simulatedResponse = t(simulations.abT),
    observedResponse = abdat$seed,
    fittedPredictedResponse = predictions.abT,
    integerResponse = TRUE
  )

# check model fit
plot(sim.abT)                               # looks ok, slightly downturned
testResiduals(sim.abT)                      # looks good, no outliers, no dispersion problems
testZeroInflation(sim.abT)                  # no zero-inflation problems
plot(AbundtAnom.mod)                        # I think this looks alright...
testTemporalAutocorrelation(sim.abT)
testSpatialAutocorrelation(sim.abT)
traceplot(AbundtAnom.mod, match.head = TRUE, varname = "beta", mfrow = c(3,3))
traceplot(AbundtAnom.mod, match.head = TRUE, varname = "r")
#AbundtAnom.mod$BUGSoutput$sims.list$beta.intercept %>% as_data_frame() %>% as_tibble() %$% acf(V1)


source(file = "~/Documents/FunCaB/figures/plotting_dim.R")

# coefficients plot
modCoefPlot <- AbundtAnom.mod$BUGSoutput$summary %>% 
  as.data.frame() %>% 
  as_tibble(rownames = "term") %>% 
  filter(grepl("beta\\[", term)) %>% 
  full_join(rNames.t, by = c(term = "i")) %>% 
  mutate(term = if_else(!is.na(term.y), term.y, term)) %>% 
  select(-term.y) %>% 
  mutate(term = case_when(
    term == "stemp7010" ~ "t",
    term == "sprecip7010" ~ "P",
    term == "TreatmentGap" ~ "Gap",
    term == "monthNaut" ~ "Autumn",
    term == "tAnom" ~ "t∆",
    term == "pAnom" ~ "SM∆",
    term == "TreatmentGap:tAnom" ~ "Gap:t∆",
    term == "TreatmentGap:pAnom" ~ "Gap:SM∆"
    ),
  term = factor(term, levels = rev(c("Gap", "t∆", "Gap:t∆", "SM∆", "Gap:SM∆", "Autumn", "t", "P"))))



modCoefPlot %>% ggplot(aes(x = mean, y = term)) +
  geom_vline(xintercept = 0, colour = "grey50", size = 0.4) +
  geom_pointintervalh(aes(xmin = `2.5%`, xmax = `97.5%`), size = 0.4) +
  geom_pointintervalh(aes(xmin = `25%`, xmax = `75%`), size = 4, ) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.4) +
  geom_hline(yintercept = c(3.5, 5.5, 7.5), colour = "grey80", size = 0.4) +
  xlab("Effect size") +
  theme(axis.title.y = element_blank()) +
  axis.dimLarge

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig7.jpg", dpi = 300, height = 5, width = 5)



# predictions plot II
AbundtAnom.mod$BUGSoutput$summary %>% 
  as.data.frame() %>% 
  as_tibble(rownames = "term") %>% 
  filter(grepl("muPred", term)) %>% 
  bind_cols(abdatY) %>% 
  mutate(stemp7010 = as.character(round(stemp7010, digits = 1)),
         tempLevel = temp7010) %>% 
  filter(precip7010 > 1.34,
         tAnom > 0) %>% 
  ggplot(aes(x = pAnom, y = mean)) +
  geom_vline(xintercept = 0, colour = "grey80", linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "grey80", linetype = "dashed") +
  geom_point(data = abdat, aes(y = seed, x = pAnom, colour = Treatment), shape = 21, alpha = 0.4) +
  geom_ribbon(alpha = 0.2, aes(ymax = `97.5%`, ymin = `2.5%`, fill = Treatment)) +
  geom_line(aes(colour = Treatment)) +
  facet_grid(monthN~tempLevel) +
  scale_color_manual(values = c("grey60", "Black")) +
  scale_fill_manual(values = c("grey60", "Black")) +
  labs(y = "seedling abundance") +
  theme_classic() +
  axis.dimLarge

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig10.jpg", dpi = 300, height = 5.5, width = 9)


abdat %>% ggplot(aes(x = soilT, y = tAnom, size = seed, colour =factor(year))) +
  geom_hline(yintercept = 0, colour = "grey50") + 
  geom_point(alpha = 0.2) + 
  facet_grid(Treatment~tempLevel) +
  theme_classic() +
  labs(x = "Air temperature at 2 m (ºC)") +
  axis.dimLarge

abdat %>% ggplot(aes(x = soilM, y = pAnom, size = seed, colour =factor(year))) +
  geom_hline(yintercept = 0, colour = "grey50") + 
  geom_point(alpha = 0.2) + 
  facet_grid(Treatment~precipLevel) +
  theme_classic() +
  axis.dimLarge



# make plots
precipdat <- abdat %>% distinct(precip7010, sprecip7010, precipLevel, precipLevelPlot) %>% as_tibble()
drawsNew <- add_draws(data = abdatY, draws = drawsAB) %>% 
  ungroup() %>% 
  filter(.draw < 200) %>%
  left_join(precipdat)

drawsNew %>% 
  #median_qi(.width = c(0.5, 0.7, 0.975)) %>% 
  ggplot(aes(x = pAnom, y = .value, colour = ordered(Treatment))) +
  geom_line(aes(group = paste(monthN, .draw)), alpha = .1) +
  #geom_point(data = abdat) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  facet_grid(precipLevelPlot~monthN)

# soil moisture deviation
abdat %>% 
  mutate(monthN = case_when(
    monthN == "spr" ~ "early",
    monthN == "aut" ~ "late"
  )) %>%
  ggplot(aes(x = pAnom, y = seed, colour = Treatment)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_point(shape = 21, size = 3) +
  geom_smooth(method = "lm") +
  facet_grid(.~monthN) +
  scale_color_manual(values = c("Black", "grey60")) +
  labs(x = "soil moisture deviation from 2009-2018 mean",
       y = "seedling number") +
  theme_classic() +
  axis.dim +
  theme(legend.title = element_blank())

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig9.jpg", dpi = 300, width = 7.5, height = 4.5)

# temperature deviation
abdat %>% 
  mutate(monthN = case_when(
    monthN == "spr" ~ "early",
    monthN == "aut" ~ "late"
  )) %>%
  ggplot(aes(x = tAnom, y = seed, colour = Treatment)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_point(shape = 21, size = 3) +
  geom_smooth(method = "lm") +
  facet_grid(.~monthN) +
  scale_color_manual(values = c("Black", "grey60")) +
  labs(x = "soil moisture deviation from 2009-2018 mean",
       y = "seedling number") +
  theme_classic() +
  axis.dim +
  theme(legend.title = element_blank())

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig9b.jpg", dpi = 300, width = 7.5, height = 4.5)
