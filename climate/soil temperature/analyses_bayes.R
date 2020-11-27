# load packages
library("tidyverse")
library("lme4")
library("broom")
library("broom.mixed")
library("rjags")
library("R2jags")
library("tidyverse")
library("tidybayes")
library("DHARMa")

# source plotting code and soil temperature data
source("~/Documents/FunCaB/figures/plotting_dim.R")
load(file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/cleanedSoilTemp.RData")

# filter for covers above 10 % (although we could consider taking this out)
fgbs <- Cover %>% filter(Treatment == "aFGB") %>% 
  gather(forbCov, graminoidCov, mossCov, key = functionalGroup, value = cover) %>%
  mutate(cover = 0)

stDat <- Cover %>% 
  gather(forbCov, graminoidCov, mossCov, key = functionalGroup, value = cover) %>% 
  filter(!is.na(cover), weather == "sunny") %>% 
  mutate(stempLevel = as.vector(scale(tempLevel, scale = FALSE, center = TRUE)),
         stemp7010 = as.vector(scale(temp7010, scale = FALSE, center = TRUE)),
         precipDiv7010 = precip7010/1000,
         sprecip7010 = as.vector(scale(precipDiv7010, scale = FALSE, center = TRUE)),
         precipDiv = precipLevel/1000,
         sprecipLevel = as.vector(scale(precipDiv, scale = FALSE, center = TRUE)),
         sCover = as.vector(scale(cover, scale = FALSE, center = TRUE))) %>% 
  select(-sforbCov, -sgraminoidCov, -smossCov, weather)


stdat <- Cover %>% 
  filter(weather == "sunny", 
         date < "2015-09-30",
         date > "2015-08-09") %>% 
  gather(forbCov, graminoidCov, mossCov, key = functionalGroup, value = cover)
stdat %>% left_join(stdat %>% filter(Treatment == "aFGB") %>% ungroup() %>% select(covFGB = cover, siteID, blockID)) %>%
  mutate(forbCov = case_when(
    grepl("F", Treatment) ~ 0,
    TRUE ~ forbCov),
    graminoidCov = case_when(
      grepl("G", Treatment) ~ 0,
      TRUE ~ graminoidCov),
    mossCov = case_when(
      grepl("B", Treatment) ~ 0,
      TRUE ~ mossCov),
    stemp7010 = as.vector(scale(temp7010, scale = FALSE, center =TRUE)),
         precipDiv7010 = precip7010/1000,
         sprecip7010 = as.vector(scale(precipDiv7010, scale = FALSE, center =TRUE)),
         sgraminoidCov = graminoidCov/100,
         sforbCov = forbCov/100,
         smossCov = mossCov/100) %>% 
  distinct(siteID, turfID, blockID, Treatment, date, maxTemp, sgraminoidCov, sforbCov, smossCov, stemp7010, sprecip7010) %>% 


## model ##
modX <- lmer(maxTemp ~ sgraminoidCov*stemp7010 + sforbCov*stemp7010 + smossCov*stemp7010 + sgraminoidCov*sprecip7010 + sforbCov*sprecip7010 + smossCov*sprecip7010 + (1|siteID), data = stdat)


testDispersion(modX)
simulationOutput <- simulateResiduals(fittedModel = modX, n = 250)
simulationOutputRE <- simulateResiduals(fittedModel = modX, n = 250, use.u = F)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))

plot(simulationOutput)
plotResiduals(simulationOutput = simulationOutput)
noNA <- stdat %>% filter(!is.na(sgraminoidCov))
plotResiduals(simulationOutput, form = noNA$sgraminoidCov)
plotResiduals(simulationOutput, form = noNA$sforbCov)
plotResiduals(simulationOutput, form = noNA$smossCov)
plotResiduals(simulationOutput, form = noNA$stemp7010)
plotResiduals(simulationOutput, form = noNA$sprecip7010)
hist(simulationOutput)

modX <- tidy(modX) %>% 
  mutate(lower = (estimate - std.error*1.96),
         upper = (estimate + std.error*1.96))

modX %>%
  filter(!grepl("^sd_", term),
         !term == "(Intercept)") %>% 
  ggplot(aes(x = term, y = estimate, ymin = lower, ymax = upper)) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #geom_vline(xintercept =  c(8.5, 16.5), colour = pal1[2], size = 46, alpha = 0.2) +
  #geom_vline(xintercept = 1.45, colour = pal1[2], size = 24, alpha = 0.2) +
  geom_point(position = position_dodge(width = 0.1), size = 4) +
  #scale_colour_manual("", values = pal1[c(2,4,5,1,3,4)]) + 
  #scale_shape_manual("", values = c(21, 22, 23, 24, 25,1), limits = c("Graminoids", "Forbs", "Bryophytes", "Intact vegetation", "Climate")) + 
  coord_flip() +
  axis.dimLarge +
  labs(y = "Soil temperature deviance from intercept (ºC)") +
  theme_classic() +
  axis.dimLarge +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom")

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 2/figures/figmodOutput.jpg", dpi = 300, width = 6, height = 4)

modX %>% 
  mutate(estimate = round(estimate, 2),
         `std.error` = round(`std.error`, 2),
         statistic = round(statistic, 2)) %>% 
  write_delim(path = "~/OneDrive - University of Bergen/Research/FunCaB/paper 2/data/modOutput.csv", delim = ";")

#############
# 
stDat <- stDat %>% 
  group_by(turfID, functionalGroup) %>% 
  mutate(presence = case_when(
    grepl("F", Treatment) & functionalGroup == "forbCov" ~ 1,
    !grepl("F", Treatment) & functionalGroup == "forbCov" ~ -1,
    grepl("G", Treatment) & functionalGroup == "graminoidCov" ~ 1,
    !grepl("G", Treatment) & functionalGroup == "graminoidCov" ~ -1,
    grepl("B", Treatment) & functionalGroup == "mossCov" ~ 1,
    !grepl("B", Treatment) & functionalGroup == "mossCov" ~ -1,
    grepl("C", Treatment) ~ 1
  ),
  deltaCov = cover*presence)


#############

lmMod1 <- lmer(maxTemp ~ deltaCov*functionalGroup*Treatment + (1|siteID), data = stDat)


#### Bayesian analysis ####
# i) set up a model matrix to feed directly into the model. This avoids potential coding errors. -1 removes the intercept, which I set separately so it can be drawn from a normal distribution.

stDat <- stDat %>% filter(presence == 1)

soilTempModMat <- model.matrix(~ deltaCov*Treatment*functionalGroup + stemp7010 + sprecip7010, data = stDat)[,-1]


stDatY <- crossing(Treatment = unique(stDat$Treatment),
                   precip7010 = mean(stDat$precip7010),
                   temp7010 = mean(stDat$temp7010),
                   deltaCov = seq(min(stDat$deltaCov) + 1, max(stDat$deltaCov) - 1, length = 50),
                   functionalGroup = unique(stDat$functionalGroup)
                   ) %>%
  mutate(stemp7010 = as.vector(scale(temp7010, scale = FALSE, center = TRUE
  )),
  precipDiv7010 = precip7010/1000,
  sprecip7010 = as.vector(scale(precipDiv7010, scale = FALSE, center = TRUE
  ))
  )

stDatY

# model matrix for fake data predictions
soilTempModMatY <- model.matrix(~ deltaCov*Treatment*functionalGroup + stemp7010 + sprecip7010, data = stDatY)

# remove intercept
soilTempModMatY <- soilTempModMatY[,-1]

# ii) model
cat("model {
  # Likelihood
  for (i in 1:n.dat) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- beta.intercept + inprod(beta, matX[i, ]) + beta.site[siteID[i]]
    
    # predictions for model validation, using original data
    yPred[i] ~ dnorm(mu[i], tau)
  }

   # derived predictions
for (k in 1:n.datY){
    muPred[k] <- beta.intercept + inprod(beta, matY[k, ])
    }

  # Priors
  tau  <- 1/(sigma*sigma)             # prior for the precision
  sigma ~ dunif(0, 100)               # residual standard deviation
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
", fill = TRUE, file = "~/Documents/FunCaB/analyses/soilTempMaxTemp.txt")

# specify the parameters to watch
paraNames.ab <- c("beta.intercept", "beta", "beta.site", "yPred", "muPred", "mu")

# iii) Set up a list that contains all the necessary data
n.treat <- nlevels(factor(stDat$Treatment))

abDat <- list(y = stDat$maxTemp, 
              n.dat = nrow(stDat),
              n.datY = nrow(stDatY),
              matX = soilTempModMat,
              matY = soilTempModMatY,
              nEff = ncol(soilTempModMat),
              siteID = as.numeric(factor(stDat$siteID)),
              n.site = nlevels(factor(stDat$siteID)))



# iv) Compile the model and run the MCMC for an adaptation/burn-in phase and sample from the posteriors
AbundtAnom.mod <- jags(
  model.file = "~/Documents/FunCaB/analyses/soilTempMaxTemp.txt",
  data = abDat,
  n.iter = 5000,
  n.chains = 3,
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
rNames.t <- colnames(soilTempModMat) %>% 
  enframe(name = "i", value = "term") %>% 
  mutate(i = paste0("beta[",i,"]"))


sim.abT <- createDHARMa(
  simulatedResponse = t(simulations.abT),
  observedResponse = stDat$maxTemp,
  fittedPredictedResponse = predictions.abT,
  integerResponse = TRUE
)

# check model fit
plot(sim.abT)                               # looks good
testResiduals(sim.abT)                      # looks good, no outliers, no dispersion problems
testZeroInflation(sim.abT)                  # no zero-inflation problems
plot(AbundtAnom.mod)                        # I think this looks alright
testTemporalAutocorrelation(sim.abT)
testSpatialAutocorrelation(sim.abT)
traceplot(AbundtAnom.mod, match.head = TRUE, varname = "beta", mfrow = c(3,3))
traceplot(AbundtAnom.mod, match.head = TRUE, varname = "r")

modCoefPlot <- AbundtAnom.mod$BUGSoutput$summary %>% 
  as.data.frame() %>% 
  as_tibble(rownames = "term") %>% 
  filter(grepl("beta\\[", term)) %>% 
  full_join(rNames.t, by = c(term = "i")) %>% 
  mutate(term = if_else(!is.na(term.y), term.y, term)) %>% 
  select(-term.y)


modCoefPlot %>% ggplot(aes(x = mean, y = term)) +
  geom_vline(xintercept = 0, colour = "grey50", size = 0.4) +
  geom_pointintervalh(aes(xmin = `2.5%`, xmax = `97.5%`), size = 0.4) +
  geom_pointintervalh(aes(xmin = `25%`, xmax = `75%`), size = 4, ) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.4) +
  geom_hline(yintercept = c(3.5, 5.5, 7.5), colour = "grey80", size = 0.4) +
  xlab("Effect size") +
  theme(axis.title.y = element_blank())

AbundtAnom.mod$BUGSoutput$summary %>% 
  as.data.frame() %>% 
  as_tibble(rownames = "term") %>% 
  filter(grepl("muPred", term)) %>% 
  bind_cols(stDatY) %>%
  ggplot(aes(x = deltaCov, y = mean)) +
  geom_point(data = stDat, aes(y = maxTemp, x = deltaCov, colour = functionalGroup), shape = 21, alpha = 0.4) +
  geom_ribbon(alpha = 0.2, aes(ymax = `97.5%`, ymin = `2.5%`, fill = functionalGroup)) +
  geom_line(aes(colour = functionalGroup)) +
  scale_color_manual(values = c("darkgoldenrod3","grey50", "#BF5C6A")) +
  scale_fill_manual(values = c("darkgoldenrod3", "grey50", "#BF5C6A")) +
  theme_classic() +
  facet_grid(functionalGroup~Treatment, scales = "free_x")

#ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 2/figures/newAnalysesII.jpg", width = 12, height = 8, dpi = 300)