
# Install/load the necessary packages and 
tryCatch(library(devtools), 
         error = function(e) {
           install.packages("devtools")
           library(devtools)
         })
tryCatch(library(tidyverse), 
         error = function(e) {
           install.packages("tidyverse")
           library(tidyverse)
         })
tryCatch(library(nimble), 
         error = function(e) {
           install.packages("nimble")
           library(nimble)
         })



# Functions to generate count data for use in an N-mixture model. 
# Variables are purposefully given confusing names so it's hard to cheat!
sim_dataset <- function(case) {
  set.seed(7887)
  blah <- 120
  blah2 <- rnorm(blah)
  blah3 <- rnorm(blah)
  h <- blah/40
  dd <- switch(case, rpois, rnbinom)
  view <- switch(case, {
    temp <- dd(blah, exp(0.5 * blah3 + 3.2))
    temp2 <- suppressMessages(bind_cols(
      lapply(1:h, function(x) rbinom(blah, temp, plogis(0.5 * blah2 - 0.2)))
    ))
    names(temp2) <- paste0("V", 1:h)
  }, {
    temp <- dd(blah, mu = 3 + exp(0.5 * blah3 + 2.2), size = 0.1)
    temp2 <- suppressMessages(bind_cols(
      lapply(1:h, function(x) 6 + rbinom(max(0, blah - 10), temp, plogis(0.5 * blah2 - 0.2)))
    ))
    names(temp2) <- paste0("V", 1:h)
  })
  
  return(list(
    y = temp2,
    x1 = blah2,
    x2 = blah3
  ))
}


library(nimble, warn.conflicts = FALSE)

code <- nimbleCode({
  ## noninformative priors
  mu ~ dflat()
  sigma ~ dhalfflat()
  ## likelihood
  for(i in 1:n) {
    y[i] ~ dnorm(mu, sd = sigma)
  }
})

data <- list(y = MASS::newcomb)
inits <- list(mu = 0, sigma = 5)
constants <- list(n = length(data$y))

model <- nimbleModel(code = code, data = data, constants = constants, inits = inits)


simNodes <- model$expandNodeNames("y")
parentNodes <- c("mu", "sigma")

# Store the true values of all the nodes we'll eventually simulate
true_simNodes_vals <- values(model, simNodes)


cmodel 	<- compileNimble(model)
mcmc    <- buildMCMC(model, monitors = parentNodes)
cmcmc   <- compileNimble(mcmc, project = model)

samples <- runMCMC(cmcmc, niter = 1000, nburnin = 500)


nSamp <- nrow(samples)
n <- length(data$y)
ppSamples <- matrix(0, nSamp, n)

set.seed(1)
for(i in 1:nSamp){
  cmodel[["mu"]] <- samples[i, "mu"]             ## or cmodel$mu <- samples[i, "mu"]
  cmodel[["sigma"]] <- samples[i, "sigma"]
  cmodel$simulate(simNodes, includeData = TRUE)
  ppSamples[i, ] <- cmodel[["y"]]
}


ppSamples <- matrix(0, nrow = nSamp, ncol =
                      length(model$expandNodeNames(dataNodes, returnScalarComponents = TRUE)))
postNames <- colnames(samples)

set.seed(1)
system.time({
  for(i in seq_len(nSamp)) {
    values(cmodel, postNames) <- samples[i, ]  # assign 'flattened' values
    cmodel$simulate(simNodes, includeData = TRUE)
    ppSamples[i, ] <- values(cmodel, dataNodes)
  }
})


obsMin <- min(data$y)
ppMin <- apply(ppSamples, 1, min)

# ## Check with plot in Gelman et al. (3rd edition), Figure 6.3
{
  hist(ppMin, xlim = c(-50, 20),
       main = "Discrepancy = min(y)", 
       xlab = "min(y_rep)")
  abline(v = obsMin, col = 'red')
}


deviance_dist <- numeric(nSamp)

# Simulated data
set.seed(1)
system.time({
  for(i in seq_len(nSamp)) {
    values(cmodel, postNames) <- samples[i, ]  # assign 'flattened' values
    cmodel$simulate(simNodes, includeData = TRUE)
    
    deviance_dist[i] <- -2 * cmodel$calculate()
  }
})

# Real data
values(cmodel, simNodes) <- true_simNodes_vals
deviance_observed <- numeric(nSamp)
system.time({
  for(i in seq_len(nSamp)) {
    values(cmodel, postNames) <- samples[i, ]  # assign 'flattened' values
    deviance_observed[i] <- -2 * cmodel$calculate()
  }
})

# Mean observed deviance
deviance_obs_mean <- mean(deviance_observed)

{
  hist(deviance_dist,
       main = "Discrepancy = deviance", 
       xlab = "deviance")
  abline(v = deviance_obs_mean, col = 'red')
}


d1 <- sim_dataset(1)
d2 <- sim_dataset(2)

d1$x1[1:10]
d2$x1[1:10]

head(d1$y)
head(d2$y)


nmix_code <- nimbleCode({
  # Model
  for (i in 1:nsite) {
    logit(p[i])    <- alpha0 + alpha1 * x1[i] + alpha2 * x2[i]
    log(lambda[i]) <- beta0 + beta1 * x1[i] + beta2 * x2[i]
    
    N[i] ~ dpois(lambda[i])
    
    for (j in 1:nobs) {
      y[i, j] ~ dbinom(size = N[i], prob = p[i])
    }
  }
  
  # Priors
  alpha0 ~ dnorm(0, sd = 5)
  alpha1 ~ dnorm(0, sd = 5)
  alpha2 ~ dnorm(0, sd = 5)
  beta0 ~ dnorm(0, sd = 10)
  beta1 ~ dnorm(0, sd = 5)
  beta2 ~ dnorm(0, sd = 5)
})

# Choose which dataset to use
dat <- d2

# Build the model
mod <- nimbleModel(
  nmix_code,
  constants = list(
    nobs = ncol(dat$y),
    nsite = nrow(dat$y),
    x1 = dat$x1,
    x2 = dat$x2
  ), 
  data = list(
    y = dat$y
  ),
  inits = list(
    alpha0 = 0, alpha1 = 0, alpha2 = 0,
    beta0  = 0, beta1  = 0, beta2  = 0,
    N = rep(max(dat$y) + 10, nrow(dat$y))
  )
)

cmod <- compileNimble(mod)
mcmc <- buildMCMC(mod, monitors = c("alpha0", "alpha1", "alpha2",
                                    "beta0", "beta1", "beta2", "N"))
cmcmc <- compileNimble(mcmc)

samples <- runMCMC(cmcmc, niter = 5000, nburnin = 2000, nchains = 2,
                   thin = 10)
summary <- MCMCvis::MCMCsummary(samples)

samples_long <- do.call(rbind, samples)
nSamp <- nrow(samples_long)
postNames <- colnames(samples_long)
