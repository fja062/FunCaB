# --- load packages --- #
library(here)
library(rjags)
library(ggplot2)
library(reshape2)
library(dplyr)
library(lme4)

# --- data --- #
sumCovmod <- filter(community_FD_analysis, trait == "sumcover")


# --- traditional model --- #
mod.lmer <- lmer(value ~ Treatment*Stemp0916*Sprecip0916 + (1|siteID/blockID), REML = FALSE, data = sumCovmod)

summary(mod.lmer)

model{
  #likelihood
  for(datIter in 1:nData){
  obs[datIter] ~ dnorm(mu[datIter], tau)
  mu[datIter] <- intercept + treatment[dataIter] + temp[datIter] + precip[dataIter] + eps[datIter]
  }
  
  
  
  #random effects
  for(site in 1:siteID){
  for(block in 1:blockID){
  eps[site[block]] ~ dnorm(0, tau.block)
  eps[site] ~ dnorm(0, tau.site)
  }
  }
  
  
  #priors
  tau.block ~ dnorm(0, 0.001)
  tau.site ~ dnorm(0, 0.001)
  tau.obs ~ dnorm(0, 0.001)
  sigma.block ~ dnorm(0, 0.001)
  sigma.site ~ dunif(0,100)
  sigma.obs ~ dunif(0,100)
  treatment ~ dnorm(0, 0.001)
  temp ~ dnorm(0, 0.001)
  precip ~ dnorm(0, 0.001)
  
  
  
  }