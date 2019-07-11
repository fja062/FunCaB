#
lmerMod <- lmer(value ~ Treatment*Stemp0916*Sprecip0916 + (1|siteID), REML = FALSE, data = dat)
summary(lmerMod)


sink("ANOVA.txt")
cat(
  "model{
  # Likelihood
  for(i in 1:n.dat){
  y[i] ~ dnorm(mu[i], tau)
  mu[i] <- mu.alpha +
    beta.treat[TreatID[i]] +
    beta.temp*temp[i] + 
    beta.precip*precip[i] + 
    beta.treat.temp[TreatID[i]]*temp[i] + 
    beta.treat.precip[TreatID[i]]*precip[i] + 
    beta.temp.precip*temp[i]*precip[i] + 
    #beta.treat.temp.precip[TreatID[i]]*precip[i]*temp[i] +
    beta.site[siteID[i]]
  }
  
  # Priors fixed effects
  for(k in 2:n.treat){
    beta.treat[k] ~ dnorm(0, 0.001)
    beta.treat.temp[k] ~ dnorm(0, 0.001)
    beta.treat.precip[k] ~ dnorm(0, 0.001)
    #beta.treat.temp.precip[k] ~ dnorm(0, 0.001)
  }
  beta.treat[1] <- 0
  beta.treat.temp[1] <- 0
  beta.treat.precip[1] <- 0
  #beta.treat.temp.precip[1] <- 0
  mu.alpha ~ dnorm(0, 0.001)
  beta.precip ~ dnorm(0, 0.001)
  beta.temp ~ dnorm(0, 0.001)
  beta.temp.precip ~ dnorm(0, 0.001)
  tau <- 1/(sigma*sigma)
  sigma ~ dunif(0, 100)


  # priors random effects
    for(m in 2:n.site){
        beta.site[m] ~ dnorm(0, 0.001)
    }
  beta.site[1] <- 0


  }

  ",fill = TRUE)

sink()

# 2) Set up a list that contains all the necessary data
n.treat <- nlevels(dat$Treatment)
n.site <- nlevels(factor(dat$siteID))
Data = list(y = dat$value, 
            TreatID = as.numeric(dat$Treatment),
            siteID = as.numeric(factor(dat$siteID)),
            n.dat = nrow(dat), 
            n.treat = n.treat,
            n.site = n.site,
            temp = dat$Stemp7010,
            precip = dat$Sprecip7010)

# check levels
levels(dat$Treatment)

# 3) Specify a function to generate inital values for the parameters
inits.fn <- function() list(alpha.treat = rnorm(n.treat),
                            beta.temp = rnorm(1),
                            sigma = runif(1, 1, 100))

# Compile the model and run the MCMC for an adaptation (burn-in) phase
jagsModel <- jags.model(file = "ANOVA.txt", data=Data,
                        n.chains = 4,
                        n.adapt = 50000
                        )

para.names <- c("mu.alpha", "beta.treat", "beta.temp", "beta.precip", "beta.treat.temp", "beta.treat.precip",  "beta.treat.temp.precip", "beta.site")

samples <- coda.samples(jagsModel,
            variable.names = para.names,
            inits = inits.fn,
            n.iter = 100000,
            thin = 4)

summary(samples)


# plot diagnostics
plot(samples)
gelman.plot(samples)
gelman.diag(samples)
