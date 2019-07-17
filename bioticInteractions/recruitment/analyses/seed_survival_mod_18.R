# load data
load(file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/cleanedSurvData.RData")
load(file = "~/OneDrive - University of Bergen/Research/FunCaB/Data/secondary/soilMoisture2018.RData")

# load libraries and scripts
library("rjags")
library("R2jags")
library("tidyverse")
library("tidybayes")
library("DHARMa")
source("~/Documents/FunCaB/figures/plotting_dim.R")

## add temperature anomaly data!
datSurv <- survival %>% 
  filter(Round == 1) %>%     # remove NAs otherwise model fails
  mutate(blockID = as.character(blockID)) %>% 
  left_join(SM2018 %>% mutate(blockID = as.character(blockID), Treatment = if_else(grepl("C", Treatment), "aC", Treatment)))
datSurv <- as.data.frame(datSurv)


# model matrix
matSurv.t <- model.matrix(~ Treatment*stemp7010*sprecip7010 - Treatment:stemp7010:sprecip7010, data = datSurv)[,-1]

datSurvY <- crossing(Treatment = unique(datSurv$Treatment),
    # rep is slowest on inside
    precip7010 = seq(600, 3029, length = 100) ,
    temp7010 = c(6.5, 8.5, 10.5)) %>% 
      mutate(sprecip7010 = (precip7010/ 1000) - mean(datSurv$precip7010/1000),
             stemp7010 =  temp7010 - mean(datSurv$temp7010))

matSurv.tY <- model.matrix( ~ Treatment*stemp7010*sprecip7010 - Treatment:stemp7010:sprecip7010, data = datSurvY)

matSurv.tY <- matSurv.tY[,-1]


##
cat("model {
  # Likelihood
  for(i in 1:n.dat) {

# Distribution of the number of surviving seedlings (using the number that started at the beginning of the season and the beta binomial parameters)
    numSurvived[i] ~ dbetabin(nbAlpha[i], nbBeta[i], N[i])

# mean survival probability & precision make the paramters for beta binomial model
    nbAlpha[i] <- meanSurvProb[i] * survPrec
    nbBeta[i] <- survPrec * (1 - meanSurvProb[i])
    
    logit(meanSurvProb[i]) <- beta.intercept +
      inprod(beta, matX[i,]) +
      beta.site[siteID[i]]
    
    # predictions for model validation, using original data
    yPred[i] ~ dbetabin(nbAlpha[i], nbBeta[i], N[i])
  }
  
# predictions
  for(j in 1:n.datY){
    #numSurvivedPred[j] ~ dbetabin(nbAlphaPred[j], nbBetaPred[j], 100)
    nbAlphaPred[j] <- meanSurvProbPred[j] * survPrec
    nbBetaPred[j] <- survPrec * (1 - meanSurvProbPred[j])

    logit(meanSurvProbPred[j]) <- beta.intercept +
      inprod(beta, matY[j,])
  }

  
  # Priors
  survPrec ~ dgamma(0.001, 0.001)     # Prior for the precision of the survival probability
  beta.intercept ~ dnorm(0, 0.001)    # intercept prior

  for(k in 1:nEff){  
    beta[k] ~ dnorm(0, 0.001)         # priors for the remaining betas
  }

  # priors random effects
  randTau ~ dgamma(0.001, 0.001)
  for(m in 1:n.site){
    beta.site[m] ~ dnorm(0, randTau)
  }  

}
", fill = TRUE, file = "~/Documents/FunCaB/analyses/funcabSurvival.txt")

# specify the parameters to watch
para.names.su <- c("beta", "beta.intercept", "beta.site", "survPrec", "meanSurvProb", "yPred", "meanSurvProbPred", "nbAlphaPred", "nbBetaPred")

# 2) Set up a list that contains all the necessary data
n.treat <- nlevels(factor(datSurv$Treatment))
n.site <- nlevels(factor(datSurv$siteID))

DataSurv <- list(numSurvived = datSurv$totS, 
              n.dat = nrow(datSurv),
              n.datY = nrow(matSurv.tY),
              matX = matSurv.t,
              matY = matSurv.tY,
              N = datSurv$tot,
              nEff = ncol(matSurv.t),
              siteID = as.numeric(factor(datSurv$siteID)), 
              n.site = n.site)

# check levels
levels(factor(datSurv$Treatment))


# Compile the model and run the MCMC for an adaptation (burn-in) phase
survModt <- jags(model.file = "~/Documents/FunCaB/analyses/funcabSurvival.txt",
  data = DataSurv,
  n.iter = 20000,
  n.chains = 5,
  parameters.to.save = para.names.su,
  progress.bar = "text",
  jags.module = "mix"
)

# create variables for model checking
simulationsSurv <- survModt$BUGSoutput$sims.list$yPred
predictionsSurv <- apply(survModt$BUGSoutput$sims.list$meanSurvProb, 2, median)
drawsSU <- simulationsSurv %>% as.data.frame()

dim(simulationsSurv)
simSurv <- createDHARMa(
  simulatedResponse = t(simulationsSurv),
  observedResponse = datSurv$totS,
  fittedPredictedResponse = predictionsSurv,
  integerResponse = TRUE
)

plot(simSurv)
testResiduals(simSurv)

survModt$BUGSoutput$summary %>% 
  as.data.frame() %>% 
  as_tibble(rownames = "term") %>% 
  filter(grepl("meanSurvProbPred", term)) %>% 
  bind_cols(datSurvY) %>% 
  ggplot(aes(x = precip7010, y = mean, ymax = `97.5%`, ymin = `2.5%`, fill = factor(temp7010))) +
  geom_ribbon(alpha = 0.2) +
  geom_line(aes(colour = factor(temp7010))) +
  facet_wrap(~Treatment)+
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")


# coefficients plot
modCoefSurv <- survModt$BUGSoutput$summary %>%
  as.data.frame() %>%
  as_tibble(rownames = "term") %>%
  filter(grepl("beta\\[", term)) %>%
  full_join(Names.t, by = c(term = "i")) %>%
  mutate(term = if_else(!is.na(term.y), term.y, term)) %>%
  select(-term.y) %>%
  mutate(
    term = gsub("Treatment", "", term),
    term = gsub("stemp7010", "t", term),
    term = gsub("sprecip7010", "P", term),
    term = gsub("FGB", "Gap", term),
    term = gsub("GB", "F", term),
    term = gsub("B", "GF", term),
    term = gsub("^G$", "FB", term),
    term = gsub("^GF:", "FB:", term),
    term = gsub("^GF:", "FB:", term),
    term = gsub(":", " x ", term),
    
  )


modCoefSurv %>% ggplot(aes(x = mean, y = term)) +
  geom_vline(xintercept = 0, colour = "grey50", size = 0.4) +
  geom_pointintervalh(aes(xmin = `2.5%`, xmax = `97.5%`), size = 1) +
  geom_pointintervalh(aes(xmin = `25%`, xmax = `75%`), size = 4, ) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.4) +
  xlab("Effect size") +
  theme(axis.title.y = element_blank()) +
  axis.dimLarge

ggsave(filename = "~/OneDrive - University of Bergen/Research/FunCaB/paper 4/figures/fig8.jpg", dpi = 300, height = 5, width = 5)




add_draws(data = datSurvY, draws = drawsSU) %>% 
  ungroup() %>%
  filter(.draw < 100) %>%
  ggplot(aes(x = precip7010, color = ordered(Treatment))) +
  geom_line(aes(y = .value, group = paste(Treatment, temp7010, .draw)), alpha = 0.1) +
  geom_point(data = datSurv%>% select(-temp7010) %>% rename(temp7010 = tempLevel), aes(x = precip7010, y = survival, size = tot), shape = 21) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(.~temp7010) 


Names.t <- colnames(matSurv.tY) %>% 
  enframe(name = "i", value = "term") %>% 
  mutate(i = paste0("beta[",i,"]"))







simulations.SM <- jagsModel.SM$BUGSoutput$sims.list$numSurvivedPred
predictions.SM <- apply(jagsModel.SM$BUGSoutput$sims.list$meanSurvProb, 2, median)
beta.SM <- jagsModel.SM$BUGSoutput$sims.list$beta

sim.SM <- createDHARMa(simulatedResponse = t(simulations.SM), observedResponse = datSurv.SM$totS, fittedPredictedResponse = predictions.SM, integerResponse = TRUE)

plot(sim.SM)
testResiduals(sim.SM)

rNames.SM2 <- colnames(matSurv.SM) %>% 
  enframe(name = "i", value = "term")

samples.SM2 %>% 
  spread_draws(c(beta[i])) %>% 
  median_qi(.width = c(0.5, 0.7, 0.95)) %>% 
  left_join(rNames.SM2) %>% 
  ggplot(aes(y = term, x = beta)) +
  #geom_eyeh()
  geom_pointintervalh() +
  geom_vline(xintercept = 0)




### model 2

matSurv.t <- model.matrix(~ Treatment*stemp7010*sprecip7010 - Treatment:stemp7010:sprecip7010, data = datSurv)[,-1]

datSurvY <- crossing(Treatment = unique(datSurv$Treatment),
                     # rep is slowest on inside
                     precip7010 = seq(600, 3029, length = 100) ,
                     temp7010 = c(6.5, 8.5, 10.5)) %>% 
  mutate(sprecip7010 = (precip7010/ 1000) - mean(datSurv$precip7010/1000),
         stemp7010 =  temp7010 - mean(datSurv$temp7010))

matSurv.tY <- model.matrix( ~ Treatment*stemp7010*sprecip7010 - Treatment:stemp7010:sprecip7010, data = datSurvY)

matSurv.tY <- matSurv.tY[,-1]


# specify the parameters to watch
para.names.su <- c("beta", "beta.intercept", "beta.site", "survPrec", "meanSurvProb", "yPred", "meanSurvProbPred", "nbAlphaPred", "nbBetaPred")

# 2) Set up a list that contains all the necessary data
n.treat <- nlevels(factor(datSurv$Treatment))
n.site <- nlevels(factor(datSurv$siteID))

DataSurv <- list(numSurvived = datSurv$totS, 
                 n.dat = nrow(datSurv),
                 n.datY = nrow(matSurv.tY),
                 matX = matSurv.t,
                 matY = matSurv.tY,
                 N = datSurv$tot,
                 nEff = ncol(matSurv.t),
                 siteID = as.numeric(factor(datSurv$siteID)), 
                 n.site = n.site)

# check levels
levels(factor(datSurv$Treatment))


# Compile the model and run the MCMC for an adaptation (burn-in) phase
survModt <- jags(model.file = "~/Documents/FunCaB/analyses/funcabSurvival.txt",
                 data = DataSurv,
                 n.iter = 20000,
                 n.chains = 5,
                 parameters.to.save = para.names.su,
                 progress.bar = "text",
                 jags.module = "mix"
)

# create variables for model checking
simulationsSurv <- survModt$BUGSoutput$sims.list$yPred
predictionsSurv <- apply(survModt$BUGSoutput$sims.list$meanSurvProb, 2, median)
drawsSU <- simulationsSurv %>% as.data.frame()

dim(simulationsSurv)
simSurv <- createDHARMa(
  simulatedResponse = t(simulationsSurv),
  observedResponse = datSurv$totS,
  fittedPredictedResponse = predictionsSurv,
  integerResponse = TRUE
)



datSurv %>% ggplot(aes(x = precip7010, y = survival, colour = Treatment)) +
  geom_point() +
  geom_line(aes(group = Treatment)) +
  facet_grid(.~tempLevel)
