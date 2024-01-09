# Prepping data and Running simple STAN SAM model for kelp nutrients
# Heili Lowman
# January 9, 2024

#### Setup ####

# Load packages.
library(here)
library(tidyverse)
library(rstan)
library(shinystan)
library(bayesplot)

# Load data.
data <- readRDS(here('dat_working',
                     "CN_temp_monthly_20230720.rds"))

# Run each time you load in "rstan"
rstan_options(auto_write=TRUE)
# auto-caches model results in the same directory
options(mc.cores=parallel::detectCores())
# during runs, each chain needs a dedicated core

#### Model Fit ####

# For this first test, simply removing NAs.
data_noNA <- drop_na(data) # removes 76 records

# Prep data for STAN as a list
data_stan <- list(
  nrows = nrow(data_noNA),
  Temp = data_noNA$mean_Temp_C,
  CN = data_noNA$mean_CN
)

# Linear model run
# .stan file must end in a BLANK LINE, otherwise it will spit out an error
# usually default is 4 chains, 2000 iterations
stan_test_run <- stan(file = "models/SAM_STAN_template_lm.stan",
                 data = data_stan,
                 chains = 3,
                 iter = 5000,
                 control = list(max_treedepth = 12))

# Examine model convergence
shinystan::launch_shinystan(stan_test_run)
# All chains appear well-mixed and no signs of divergences.

# Examine summaries of the estimates.
stan_test_data <- summary(stan_test_run,
                  pars = c("b", "b0", "tau"),
                  probs = c(0.025, 0.5, 0.975))$summary # 2.5% and 97.5% percentiles
# Rhat values also look pretty good.

# Fun plotting ideas at: https://mc-stan.org/bayesplot/reference/MCMC-intervals.html
color_scheme_set("teal")
mcmc_areas(stan_test_run,
           pars = c("b", "b0"),
           point_est = "median",
           prob = 0.95) +
  labs(
    title = "Posterior distributions",
    subtitle = "with medians and 95% intervals"
  )

# End of script.
