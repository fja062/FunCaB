# load packages
library("rjags")
library("R2jags")
library("tidyverse")
library("tidybayes")
library("DHARMa")
source("~/Documents/FunCaB/figures/plotting_dim.R")


# load data
load("~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/cleanedAbundData2018.RData")

dat18 <- seed2018 %>% 
  filter(!is.na(seed)) %>%     # remove NAs otherwise model fails
  mutate(blockID = as.character(blockID),
         Treatment = factor(Treatment, levels = c("Intact", "Gap", "F", "B", "G", "FB", 'GF', "GB"))) %>% 
  left_join(monthAv) %>% 
  left_join(weather) %>% 
  mutate(stempLevel = as.vector(scale(tempLevel, scale = FALSE, center =TRUE)),
         stemp7010 = as.vector(scale(temp7010, scale = FALSE, center =TRUE)),
         precipDiv7010 = precip7010/1000,
         sprecip7010 = as.vector(scale(precipDiv7010, scale = FALSE, center =TRUE)),
         monthN = factor(monthN, levels = c("spr", "aut")),
         precipDiv = precipLevel/1000,
         sprecipLevel = as.vector(scale(precipDiv, scale = FALSE, center = TRUE)),
         soilTs = as.vector(scale(soilT, scale = FALSE, center = TRUE)),
         soilMs = as.vector(scale(soilM, scale = FALSE, center = TRUE)),
  )

dat18 %>% ggplot(aes(x = temp7010, y = seed)) +
  geom_point() +
  facet_grid(monthN~precipLevel)

#### Bayesian analysis ####
# i) set up a model matrix to feed directly into the model. This avoids potential coding errors. -1 removes the intercept, which I set separately so it can be drawn from a normal distribution.
matab.t18 <- model.matrix(~ monthN + Treatment + stemp7010 + sprecip7010 + stemp7010:monthN + sprecip7010:monthN + monthN:Treatment, data = dat18)[,-1]

# fake data for predictions
datY18 <- crossing(Treatment = unique(dat18$Treatment),
  # rep is slowest on inside
  monthN = unique(dat18$monthN),
  temp7010 = unique(dat18$tempLevel),
  precip7010 = seq(min(dat18$precipDiv7010), max(dat18$precipDiv7010), length = 50)
) %>%
  mutate(stemp7010 = as.vector(scale(temp7010, scale = FALSE, center = TRUE
  )),
  sprecip7010 = as.vector(scale(precip7010, scale = FALSE, center = TRUE
  )))
datY18

# model matrix for fake data predictions
matab.tY18 <- model.matrix(~ monthN + Treatment + stemp7010 + sprecip7010 + stemp7010:monthN + sprecip7010:monthN + monthN:Treatment, data = datY18)
# remove intercept
matab.tY18 <- matab.tY18[,-1]


## ---->
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
", fill = TRUE, file = "~/Documents/FunCaB/analyses/seedAbund_tAnom18.txt")


# specify the parameters to watch
paraNames.ab <- c("beta.intercept", "beta", "beta.site", "r", "yPred", "muPred", "mu")

# iii) Set up a list that contains all the necessary data
n.treat <- nlevels(factor(dat18$Treatment))
n.season <- nlevels(factor(dat18$monthN))

modat18 <- list(y = dat18$seed, 
              n.dat = nrow(dat18),
              n.datY = nrow(datY18),
              matX = matab.t18,
              matY = matab.tY18,
              nEff = ncol(matab.t18),
              siteID = as.numeric(factor(dat18$siteID)),
              n.site = nlevels(factor(dat18$siteID)))



# iv) Compile the model and run the MCMC for an adaptation/burn-in phase and sample from the posteriors
AbundtAnom.mod <- jags(
  model.file = "~/Documents/FunCaB/analyses/seedAbund_tAnom18.txt",
  data = modat18,
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
rNames.t <- colnames(matab.t18) %>% 
  enframe(name = "i", value = "term") %>% 
  mutate(i = paste0("beta[",i,"]"))



sim.abT <- createDHARMa(
  simulatedResponse = t(simulations.abT),
  observedResponse = dat18$seed,
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


modCoefPlot <- AbundtAnom.mod$BUGSoutput$summary %>% 
  as.data.frame() %>% 
  as_tibble(rownames = "term") %>% 
  filter(grepl("beta\\[", term)) %>% 
  full_join(rNames.t, by = c(term = "i")) %>% 
  mutate(term = if_else(!is.na(term.y), term.y, term)) %>% 
  select(-term.y) %>% 
  mutate(
    term = gsub("Treatment", "", term),
    term = gsub("stemp7010", "t", term),
    term = gsub("sprecip7010", "P", term),
    term = gsub("monthNaut", "late", term),
    term = gsub(":", " x ", term)
  ) %>% 
  mutate(term = factor(term, levels = rev(c("Gap","late x Gap", "G", "late x G", "F", "late x F", "B", "late x B", "GF", "late x GF", "FB", "late x FB", "GB", "late x GB", "late", "t","late x t", "P", "late x P"))))


modCoefPlot %>% ggplot(aes(x = mean, y = term)) +
  geom_vline(xintercept = 0, colour = "grey50", size = 0.4) +
  geom_pointintervalh(aes(xmin = `2.5%`, xmax = `97.5%`), size = 0.4) +
  geom_pointintervalh(aes(xmin = `25%`, xmax = `75%`), size = 4, ) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.4) +
  geom_hline(yintercept = c(5.5, 11.5), colour = "grey80", size = 0.4) +
  xlab("Effect size") +
  theme(axis.title.y = element_blank()) +
  axis.dimLarge

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig12.jpg", dpi = 300, height = 6, width = 5)

AbundtAnom.mod$BUGSoutput$summary %>% 
  as.data.frame() %>% 
  as_tibble(rownames = "term") %>% 
  filter(grepl("muPred", term)) %>% 
  bind_cols(datY18) %>% 
  filter(Treatment %in% c("Intact", "Gap", "F", "G", "B")) %>% 
  mutate(stemp7010 = as.character(round(stemp7010, digits = 1)),
         tempLevel = temp7010,
         monthN = factor(monthN, levels = c("spr", "aut"), labels = c("early", "late"))) %>% 
  ggplot(aes(x = precip7010, y = mean)) +
  geom_point(data = dat18 %>% 
               filter(Treatment %in% c("Intact", "Gap", "F", "G", "B")) %>% 
               mutate(monthN = factor(monthN, levels = c("spr", "aut"), labels = c("early", "late"))), aes(y = seed, x = precipDiv, colour = factor(tempLevel)), shape = 21, alpha = 0.4) +
  geom_ribbon(alpha = 0.2, aes(ymax = `97.5%`, ymin = `2.5%`, fill = factor(tempLevel))) +
  geom_line(aes(colour = factor(tempLevel))) +
  facet_grid(monthN~Treatment) +
  scale_color_manual(values = c("grey80" ,"grey50", "grey20")) +
  scale_fill_manual(values =  c("grey80" ,"grey50", "grey20")) +
  labs(y = "seedling abundance", x = "precipitation (m/yr)") +
  theme_classic() +
  axis.dimLarge +
  theme(legend.title = element_blank())

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig11.jpg", dpi = 300, height = 5.5, width = 10)
