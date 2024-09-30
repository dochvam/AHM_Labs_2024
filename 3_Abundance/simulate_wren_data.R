library(tidyverse)

set.seed(63240)

nvisit <- 4
nsite <- 32

site_info <- data.frame(
  site_ID = 1:nsite,
  elevation = runif(nsite, 50, 100),
  management_unit = rep(1:4, each = nsite/4),
  habitat = sample(c("forest", "mixed_forest_developed"), size = nsite,
                   replace = TRUE)
)

visit_info <- as.data.frame(expand.grid(
  site_ID = 1:nsite, visit_ID = 1:nvisit
)) %>% 
  mutate(temp_C = runif(nsite*nvisit, 10, 30))


# true parameters
# det ~ habitat_mixed + temp_C + (1 | unit)
alpha0 <- -0.5
alpha1 <- 1.5
alpha2 <- 0.1
det_ranef_sd <- 0.1

# abund ~ habitat_mixed + elevation + (1 | unit)
beta0 <- 1.7
beta1 <- -0.8
beta2 <- 0.5
abund_ranef_sd <- 0.6

abund_ranef <- rnorm(4, 0, sd = abund_ranef_sd)
site_info$true_lambda <- exp(
  beta0 + beta1 * as.numeric(site_info$habitat == "mixed_forest_developed") +
    beta2 * scale(site_info$elevation) + abund_ranef[site_info$management_unit]
)
site_info$true_N <- rpois(nsite, site_info$true_lambda)

all_dat <- left_join(site_info, visit_info)

det_ranef <- rnorm(4, 0, sd = det_ranef_sd)
all_dat$true_p <- nimble::expit(
  alpha0 + alpha1 * as.numeric(all_dat$habitat == "mixed_forest_developed") +
    alpha2 * scale(all_dat$temp_C) + det_ranef[all_dat$management_unit]
)
all_dat$y <- rbinom(nrow(all_dat), all_dat$true_N, all_dat$true_p)

all_dat %>% 
  select(site_ID, visit_ID, elevation, temp_C, habitat, management_unit, 
         observed_count = y) %>% 
  write_csv("3_Abundance/woodpecker_data.csv")
