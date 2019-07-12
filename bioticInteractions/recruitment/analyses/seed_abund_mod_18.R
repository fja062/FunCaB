#
library(rjags)

# load data
dat <- seedTot %>% 
  filter(!is.na(seed)) %>%     # remove NAs otherwise model fails
  mutate(blockID = as.character(blockID))


## ---->
sink("ANOVAfuncab.txt")
cat(
  "model{
  # Likelihood
  for(i in 1:n.dat){
  y[i] ~ dnegbin(p[i], r)
  p[i] <- r/(r + lambda[i])
  log(lambda[i]) <- mu[i]
  mu[i] <- beta.treat[TreatID[i]] +
  beta.temp*temp[i] + 
  #beta.precip*precip[i] + 
  #beta.treat.temp[TreatID[i]]*temp[i] + 
  #beta.treat.precip[TreatID[i]]*precip[i] + 
  #beta.temp.precip*temp[i]*precip[i] + 
  #beta.treat.temp.precip[TreatID[i]]*precip[i]*temp[i] +
  beta.site[siteID[i]]
  
  }
  
  # Priors fixed effects
  for(k in 2:n.treat){
  beta.treat[k] ~ dnorm(0, 0.001)
  #beta.treat.temp[k] ~ dnorm(0, 0.001)
  #beta.treat.precip[k] ~ dnorm(0, 0.001)
  #beta.treat.temp.precip[k] ~ dnorm(0, 0.001)
  }
  r ~ dunif(0,50)
  
  beta.treat[1] <- 0
  beta.treat.temp[1] <- 0
  #beta.treat.precip[1] <- 0
  #beta.treat.temp.precip[1] <- 0
  #beta.precip ~ dnorm(0, 0.001)
  beta.temp ~ dnorm(0, 0.001)
  #beta.temp.precip ~ dnorm(0, 0.001)
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
n.treat <- nlevels(factor(dat$Treatment))
n.site <- nlevels(factor(dat$siteID))

dat <- as.data.frame(dat)
Data <- list(y = dat$seed, 
             TreatID = as.numeric(factor(dat$Treatment)),
             siteID = as.numeric(factor(dat$siteID)),
             n.dat = nrow(dat), 
             n.treat = n.treat,
             n.site = n.site,
             temp = dat$stemp7010) %>% 
  map(as.vector)

# check levels
levels(factor(dat$Treatment))

# 3) Specify a function to generate inital values for the parameters
inits.fn <- function() list(alpha.treat = rnorm(n.treat),
                            beta.temp = rnorm(1),
                            sigma = runif(1, 1, 100),
                            seed = 5)

# Compile the model and run the MCMC for an adaptation (burn-in) phase
jagsModel <- jags.model(file = "ANOVAfuncab.txt", data = Data,
                        n.chains = 4,
                        n.adapt = 5000
)

para.names <- c("beta.treat", "beta.temp", "beta.treat.temp", "beta.site", "r")

#, "beta.temp", "beta.precip", "beta.treat.temp", "beta.treat.precip",
samples <- coda.samples(jagsModel,
                        variable.names = para.names,
                        inits = inits.fn,
                        n.iter = 2000,
                        thin = 4
)

dim(dat)
str(samples)
summary(samples)

samples %>% as.data.frame() %>% fortify() %>% ggplot(aes(x = bet.temp, y = `beta.treat[2]`)) +geom_point()


# plot diagnostics
plot(samples)
gelman.plot(samples[,c("r", "beta.temp")])
gelman.diag(samples[,c("r", "beta.temp")])
