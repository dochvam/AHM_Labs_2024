
knitr::opts_chunk$set(echo = TRUE)

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

download.file("https://www.stacksjournal.org/wp-content/uploads/Articles/Parker-Shames-24003/Parker-Shames-24003-Dataset.zip", destfile = "Lab4_data.zip")

unzip("Lab4_data.zip")

cannabis_data <- read_csv("Parker-Shames-24003-Dataset/Parker-Shames-24003-Data.csv") %>% 
  mutate(Date = mdy(Date),
         Site = as.numeric(as.factor(Site)))

# Number of visits to each site
head(count(cannabis_data, Site))

# Counts of fox detections per site
cannabis_data %>% 
  group_by(Site) %>% 
  summarize(fox = sum(fox)) %>% 
  ggplot() + 
  geom_col(aes(Site, fox))

# How many sites had at least one detection of fox?
cannabis_data %>% 
  group_by(Site) %>% 
  summarize(fox = sum(fox) > 0) %>% 
  .$fox %>% 
  table()


cannabis_data <- cannabis_data %>% 
  arrange(Site)

site_data <- cannabis_data %>% 
  distinct(
    Site, Reg1_intercept, Reg2_intercept, Reg3_intercept,
    dist_cannabis_sqrt, elev_scaled, forest_prop_scaled, 
    dist_paved_scaled
  )




# Define the model:
code <- nimbleCode({
  
  # Site-level values
  for (i in 1:nsite) {
    logit(psi[i]) <- beta0 + beta1 * dist_cannabis_sqrt[i]
    z[i] ~ dbern(psi[i])
  }
  
  # Observation-level values
  for (o in 1:nobs) {
    logit(p[o]) <- alpha0
    y[o] ~ dbern(z[site[o]] * p[o])
  }
  
  # Priors
  beta0 ~ dnorm(0, sd = 5)
  beta1 ~ dnorm(0, sd = 5)
  alpha0 ~ dnorm(0, sd = 5)
})

# Build the model:
occ_model <- nimbleModel(
  code = code,
  constants = list(
    nobs = nrow(cannabis_data), 
    nsite = nrow(site_data),
    site = cannabis_data$Site
  ), 
  data = list(
    y = cannabis_data$fox,
    dist_cannabis_sqrt = site_data$dist_cannabis_sqrt
  ),
  inits = list(
    beta0 = 0, beta1 = 0, alpha0 = 0, z = rep(1, nrow(site_data))
  )
)



mcmcConf <- configureMCMC(occ_model)
mcmc <- buildMCMC(mcmcConf)

complist <- compileNimble(occ_model, mcmc)

samples <- runMCMC(complist$mcmc, niter = 10000, nburnin = 5000,
                   nchains = 3, samplesAsCodaMCMC = TRUE)

summary <- MCMCvis::MCMCsummary(samples)

View(summary)