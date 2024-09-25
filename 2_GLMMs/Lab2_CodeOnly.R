source("helper/Prepare_EcoData.R")

library(tidyverse)



# Set a seed for replicability
set.seed(59670)

# Number of plants, number of fields
n_plants <- 100
n_fields <- 5

# create a vector of IDs indicating which field contained each plant.
field <- rep(1:n_fields, each = n_plants / n_fields)

# Simulate a moisture level for each field. The range of moisture is not even across each field. 
# Keep in mind that a linear model does not make any assumptions about the distribution of covariate data!
moisture <- c(runif(n_plants / n_fields, 2, 12),
              runif(n_plants / n_fields, 5, 15),
              runif(n_plants / n_fields, 6, 18),
              runif(n_plants / n_fields, 8, 18),
              runif(n_plants / n_fields, 9, 18))

# Center and scale these for convenience
moisture <- as.numeric(scale(moisture))

# So far we've just simulated moisture_i arbitrarily. Now let's get to the model. 
# We need to pick values for true parameters. (Remember these are log-scale)
beta0 <- 1.5
beta1 <- 0.6
sigma <- 1.1

# We need to simulate a random effect level for each field
alpha <- rnorm(n = n_fields, mean = 0, sd = sigma)
alpha

# Now we can calculate lambda deterministically, following the model eqn
log_lambda <- beta0 + beta1 * moisture + alpha[field] # <- nested indexing
lambda <- exp(log_lambda)

# Finally, draw random fruit counts from the Poisson
n_fruit <- rpois(n = n_plants, lambda)

# Put these in a data.frame
dat <- data.frame(
  plant_ID = 1:n_plants,
  field = field,
  n_fruit = n_fruit,
  moisture_scaled = moisture,
  true_log_lambda = log_lambda
)

ggplot(dat, aes(moisture_scaled, n_fruit)) +
  geom_point() +
  xlab("Moisture") + ylab("Number of fruit") +
  theme_minimal() +
  geom_smooth(method = "lm")


ggplot(dat,
       aes(moisture_scaled, n_fruit, 
           col = as.factor(field)), group = field) +
  geom_point() +
  xlab("Moisture") + ylab("Number of fruit") +
  theme_minimal() +
  geom_smooth(method = "lm")


glm_fit <- glm(n_fruit ~ moisture_scaled, family = "poisson", data = dat)
summary(glm_fit)

library(lme4)
glmm_fit <- glmer(n_fruit ~ moisture_scaled + (1 | field), 
                  family = "poisson", data = dat)
summary(glmm_fit)

beetle_data <- EcoData::volcanoisland
beetle_data$altitude_scaled <- as.numeric(scale(beetle_data$altitude))

head(beetle_data)



library(nimble)

GLM_code <- nimbleCode({
  
  # Data likelihood
  for (i in 1:nobs) {
    log(mu[i]) <- beta0 + beta1 * altitude[i] + beta2 * altitude[i]^2
    y[i] ~ dpois(mu[i])
  }
  
  # Priors for our parameters
  beta0 ~ dnorm(0, sd = 10)
  beta1 ~ dnorm(0, sd = 10)
  beta2 ~ dnorm(0, sd = 10)
})

GLMM_code <- nimbleCode({
  # Data likelihood
  for (i in 1:nobs) {
    log(mu[i]) <- beta0 + beta1 * altitude[i] + beta2 * altitude[i]^2 +
      ranef_plot[plot[i]]
    y[i] ~ dpois(mu[i])
  }
  
  # Priors for our parameters
  beta0 ~ dnorm(0, sd = 10)
  beta1 ~ dnorm(0, sd = 10)
  beta2 ~ dnorm(0, sd = 10)
  sigma ~ dgamma(shape = 1, rate = 1)
  
  # Random effect distribution
  for (i in 1:n_plot) {
    ranef_plot[i] ~ dnorm(0, sd = sigma)
  }
})

my_GLMM <- nimbleModel(
  code = GLMM_code,
  data = list(
    altitude = beetle_data$altitude_scaled,
    y = beetle_data$beetles
  ),
  constants = list(
    nobs = nrow(beetle_data),
    n_plot = length(unique(beetle_data$plot)),
    plot = as.numeric(beetle_data$plot) # plot
  ),
  inits = list(
    beta0 = 0,
    beta1 = 1,
    beta2 = 1,
    sigma = 1, # new
    ranef_plot = rnorm(length(unique(beetle_data$plot)), 0, 1) # new
  )
)

mcmcConf <- configureMCMC(my_GLMM)
mcmc <- buildMCMC(mcmcConf)

compiled_list <- compileNimble(my_GLMM, mcmc)

samples <- runMCMC(compiled_list$mcmc, niter = 10000, nburnin = 2000,
                   thin = 1, nchains = 3, samplesAsCodaMCMC = TRUE)

summary <- MCMCvis::MCMCsummary(samples)

summary

inat_data <- read_csv("data/inat_EGSq_dat.csv")

head(inat_data)

iNat_GLMM_code <- nimbleCode({
  
  # Code goes here
  
  
  
  # Priors
  beta0 ~ dnorm(0, sd = 5)
  beta1 ~ dnorm(0, sd = 5)
  beta2 ~ dnorm(0, sd = 5)
  beta3 ~ dnorm(0, sd = 5)
  beta4 ~ dnorm(0, sd = 5)
  sigma_year ~ dgamma(1, 1)
  sigma_cell ~ dgamma(1, 1)
})

iNat_GLMM_model <- nimbleModel(
  iNat_GLMM_code,
  constants = list(
    
  ),
  data = list(
    
  ), 
  inits = list(
    
  )
)

