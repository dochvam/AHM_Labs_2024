source("aux/Prepare_EcoData.R")


beetle_data <- EcoData::volcanoisland


hist(beetle_data$beetles)

# Altitude:
hist(beetle_data$altitude)
# Mean wind speed:
hist(beetle_data$windMean)
# Habitat quality measure:
hist(beetle_data$habitatQuality)
# Year:
hist(beetle_data$year)


plot(beetle_data$altitude, beetle_data$beetles)

plot(beetle_data$habitatQuality, beetle_data$beetles)

# Make year a factor so we get boxplots
plot(as.factor(beetle_data$year), beetle_data$beetles)

beetle_data$altitude_scaled <- as.numeric(scale(beetle_data$altitude))

# Store the mean and sd of the originals so we can back-transform later
altitude_mean <- mean(beetle_data$altitude)
altitude_sd <-     sd(beetle_data$altitude)

hist(beetle_data$altitude_scaled)



fit <- glm(formula = beetles ~ altitude_scaled, 
           data = beetle_data, 
           family = poisson(link = "log"))
summary(fit)


library(nimble)
my_code <- nimbleCode({
  
  # Data likelihood
  for (i in 1:nobs) {
    log(mu[i]) <- beta0 + beta1 * altitude[i]
    y[i] ~ dpois(mu[i])
  }
  
  # Priors for our two parameters
  beta0 ~ dnorm(0, sd = 10)
  beta1 ~ dnorm(0, sd = 10)
  
})




my_model <- nimbleModel(
  code = my_code,
  data = list(
    altitude = beetle_data$altitude_scaled,
    y = beetle_data$beetles
  ),
  constants = list(
    nobs = nrow(beetle_data)
  ),
  inits = list(
    beta0 = 0,
    beta1 = 1
  )
)


my_model$y[1:10]
my_model$mu[1:10]

# Calculate the log probability of y | beta0, beta1
my_model$calculate("y")


mcmcConf <- configureMCMC(my_model)
mcmc <- buildMCMC(mcmcConf)

compiled_list <- compileNimble(my_model, mcmc)

samples <- runMCMC(compiled_list$mcmc, niter = 10000, nburnin = 2000,
                   thin = 1, nchains = 3, samplesAsCodaMCMC = TRUE)

summary <- MCMCvis::MCMCsummary(samples)


plot(samples[, 1])
plot(samples[, 2])


# Comparing MLE and Bayesian results
# MLE estimates
summary(fit)

# Bayesian estimates
summary[, 1]


beetle_abund_pred <- data.frame(
  altitude = seq(-2, 2, length.out = 13),
  mean = NA, Q025 = NA, Q095 = NA
)

# Chains don't matter here, so we combine all samples into a matrix:
samples_mtx <- as.matrix(samples)
# Make a matrix where each row is a MCMC iteration and each column is an altitude
lambda_mtx <- matrix(NA, nrow = nrow(samples_mtx), ncol = nrow(beetle_abund_pred))

# Loop over posterior samples and calculate expected abundance
for (i in 1:nrow(lambda_mtx)) {
  lambda_mtx[i, ] <- exp(samples_mtx[i, "beta0"] + 
                           samples_mtx[i, "beta1"] * beetle_abund_pred$altitude) 
}

# Summarize posteriors
beetle_abund_pred$mean <- colMeans(lambda_mtx)
beetle_abund_pred$Q025 <- apply(lambda_mtx, 2, quantile, probs = 0.025)
beetle_abund_pred$Q975 <- apply(lambda_mtx, 2, quantile, probs = 0.975)

# Back-transform altitude to the correct scale
beetle_abund_pred$altitude_unscaled <-
  beetle_abund_pred$altitude * altitude_sd + altitude_mean

# Plot the output
library(ggplot2)
ggplot(beetle_abund_pred, aes(altitude, mean, ymin = Q025, ymax = Q975)) +
  geom_ribbon(alpha = 0.3) +
  geom_line() +
  theme_minimal() + xlab("Altitude (m)") + ylab("Expected count of beetles")