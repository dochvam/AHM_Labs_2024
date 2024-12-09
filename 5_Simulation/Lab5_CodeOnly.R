
# Install/load the necessary packages and 
tryCatch(library(devtools), 
         error = function(e) {
           install.packages("devtools")
           library(devtools)
         })
tryCatch(library(lubridate), 
         error = function(e) {
           install.packages("lubridate")
           library(lubridate)
         })
tryCatch(library(EcoData), 
         error = function(e) {
           devtools::install_github(repo = "TheoreticalEcology/EcoData", 
                                    dependencies = T, build_vignettes = T)
           library(EcoData)
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
tryCatch(library(ubms), 
         error = function(e) {
           install.packages("ubms")
           library(ubms)
         })



set.seed(1001)

### Run where the true effect size is 0.5
true_beta <- 0.5
true_mean <- 1
true_error_sd <- 1

ndat <- c(5, 10, 15, 20, 25, 30, 500)
nrep <- 1000

power1 <- numeric(length(ndat))

for (i in 1:length(ndat)) {
  significant <- logical(nrep)
  for (j in 1:nrep) {
    x <- runif(ndat[i], -1, 1)
    y <- rnorm(ndat[i], true_mean + true_beta * x, true_error_sd)
    
    this_fit <- lm(y ~ x)
    s <- summary(this_fit)
    
    # is the p-value for the slope significant?
    significant[j] <- s$coefficients[2, 4] < 0.05
  }
  power1[i] <- mean(significant) # Power is the rate at which we find an effect
}

### Repeat with true_beta = 0.2
true_beta <- 5
power2 <- numeric(length(ndat))

for (i in 1:length(ndat)) {
  significant <- logical(nrep)
  for (j in 1:nrep) {
    x <- runif(ndat[i], -1, 1)
    y <- rnorm(ndat[i], true_mean + true_beta * x, true_error_sd)
    
    this_fit <- lm(y ~ x)
    s <- summary(this_fit)
    
    # is the p-value for the slope significant?
    significant[j] <- s$coefficients[2, 4] < 0.05
  }
  power2[i] <- mean(significant) # Power is the rate at which we find an effect
}

bind_rows(
  data.frame(
    power = power1,
    ndat = ndat,
    beta = 0.5
  ),
  data.frame(
    power = power2,
    ndat = ndat,
    beta = 5
  )) %>% 
  ggplot() +
  geom_line(aes(ndat, power, group = beta, color = as.factor(beta))) +
  theme_minimal()



simulate_data <- function(nsites_per_forest, 
                          nreps,
                          true_forest_effect, 
                          true_forest_sd,
                          logit_occupancy_intercept,
                          logit_detection_intercept) {
  
  dist_edge <- seq(from = 0, to = 2000, length.out = nsites_per_forest)
  
  siteCovs <- as.data.frame(expand.grid(
    dist_edge = dist_edge, 
    forest_ID = 1:5
  )) %>% 
    mutate(site_ID = row_number())
  
  # Scale distance-to-edge
  siteCovs$dist_edge = as.numeric(scale(siteCovs$dist_edge))
  
  # Simulate the forest random effect
  forest_ranef <- rnorm(n = 5, mean = 0, sd = true_forest_sd)
  
  # Calculate the occupancy prob for each site
  psi <- expit(
    logit_occupancy_intercept + true_forest_effect * siteCovs$dist_edge +
      forest_ranef[siteCovs$forest_ID]
  )
  
  p <- expit(logit_detection_intercept)
  
  z <- numeric(nrow(siteCovs))
  y <- matrix(NA, nrow = nrow(siteCovs), ncol = nreps)
  
  for (i in 1:nrow(y)) {
    # Simulate the latent state
    z[i] <- rbinom(n = 1, size = 1, prob = psi[i])
    
    y[i,] <- rbinom(n = nreps, size = 1, prob = z[i] * p)
  }
  
  return(list(
    # Return the observed data
    y = y,
    siteCovs = siteCovs,
    umf = unmarkedFrameOccu(y = y, siteCovs = siteCovs),
    
    # Also return the inputs so we can easily track them later
    sim_specification = data.frame(
      nsites_per_forest = nsites_per_forest,
      nreps = nreps,
      true_forest_effect = true_forest_effect,
      true_forest_sd = true_forest_sd,
      logit_occupancy_intercept,
      logit_detection_intercept
    )
  ))
}


set.seed(6644)

datlist <- simulate_data(nsites_per_forest = 5, 
                         nreps = 6, 
                         true_forest_effect = 1, 
                         true_forest_sd = 0.2, 
                         logit_occupancy_intercept = 0, 
                         logit_detection_intercept = 0)

# View the data matrix (red = detected)
image(datlist$y)


# ubms

# Fit an occupancy model with no det covariates, with dist. to forest edge
# and forest ID random effect on occu
silence <- capture.output(
  test_fit <- stan_occu(~ 1 ~ dist_edge + (1|forest_ID), data = datlist$umf,
                        iter = 5000, warmup = 2500, verbose = FALSE, refresh=-1)
)
summary(test_fit, submodel = "state")
summary(test_fit, submodel = "det")




fit_and_evaluate_one <- function(this_datlist) {
  silence <- suppressWarnings(capture.output(
    this_fit <- stan_occu(~ 1 ~ dist_edge + (1|forest_ID), data = this_datlist$umf,
                          iter = 4000, warmup = 2000, verbose = FALSE, refresh=-1)
  ))
  
  s <- summary(this_fit, "state")
  
  # Is the credible interval different from 0
  return(sign(s$`2.5%`)[2] == sign(s$`97.5%`)[2])
}



set.seed(229983)

conditions <- data.frame(
  nsites_per_forest = c(5, 6, 10, 15),
  nreps = c(6, 5, 3, 2)
) %>% 
  mutate(nsurvey_tot = nsites_per_forest * nreps * 5,
         power = NA)

conditions

# Keep this low. In reality, we need way more (>100), but it's slow to estimate these
nsim <- 5

for (i in 1:nrow(conditions)) {
  significant <- logical(nsim)
  for (j in 1:nsim) {
    datlist <- simulate_data(nsites_per_forest = conditions$nsites_per_forest[i],
                             nreps = conditions$nreps[i], 
                             true_forest_effect = 1, 
                             true_forest_sd = 0.25, 
                             logit_occupancy_intercept = 0, 
                             logit_detection_intercept = 0)
    
    significant[j] <- fit_and_evaluate_one(datlist)
  }
  conditions$power[i] <- mean(significant)
}

# We need more simulations for this to be stable
conditions 