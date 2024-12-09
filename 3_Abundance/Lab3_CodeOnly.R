
# Install/load the necessary packages and 
tryCatch(library(devtools), 
         error = function(e) {
           install.packages("devtools")
           library(devtools)
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

dat <- read_csv("https://raw.githubusercontent.com/dochvam/AHM_Labs_2024/main/3_Abundance/woodpecker_data.csv")

summary(dat)

nmix_code <- nimbleCode({
  # Model code goes here:
  
  
  # Priors go here:
  alpha0 ~ dnorm(0, 1)
  #   What are the parameters to be estimated?
})


mod <- nimbleModel(
  nmix_code,
  constants = list(
    
  ),
  data = list(
    
  ),
  inits = list( 
    
  )
)
  
cmod <- compileNimble(mod)
mcmc <- buildMCMC(mod)
cMCMC <- compileNimble(mcmc)

samples <- runMCMC(cMCMC, niter = 5000, nburnin = 1000, nchains = 2, thin = 1)

summary1 <- MCMCvis::MCMCsummary(samples)

summary1


