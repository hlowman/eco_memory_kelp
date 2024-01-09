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
library(data.table)

# Load raw data.
data <- readRDS(here("dat_working",
                     "CN_temp_monthly_20230720.rds"))

# Run each time you load in "rstan"
rstan_options(auto_write=TRUE)
# auto-caches model results in the same directory
options(mc.cores=parallel::detectCores())
# during runs, each chain needs a dedicated core

#### LM Model Fit ####

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

#### SAM Model Fit ####

# From Ana's JAGS script, need to make the necessary data transformations.
# And going to replace NAs with mean values for the time being.
data_noNA_sam <- data %>%
  # replace NAs with mean
  mutate(CN_gapfilled = replace_na(mean_CN, mean(mean_CN, na.rm = TRUE))) %>%
  mutate(Temp_gapfilled = replace_na(mean_Temp_C, mean(mean_Temp_C, na.rm = TRUE)))

# Create new dataset of lags using gapfilled data.
Temp_lags_sam <- data_noNA_sam %>%
  dplyr::select(SITE, year, month, Temp_gapfilled) %>%
  group_by(SITE) %>%
  arrange(SITE, year, month) %>%
  #this creates a column for every lag 1:12 months ago
  do(data.frame(., setNames(shift(.$Temp_gapfilled, 1:12), c("mean_Temp_Cl1",
                                                          "mean_Temp_Cl2",
                                                          "mean_Temp_Cl3",
                                                          "mean_Temp_Cl4",
                                                          "mean_Temp_Cl5",
                                                          "mean_Temp_Cl6",
                                                          "mean_Temp_Cl7",
                                                          "mean_Temp_Cl8",
                                                          "mean_Temp_Cl9",
                                                          "mean_Temp_Cl10",
                                                          "mean_Temp_Cl11",
                                                          "mean_Temp_Cl12")))) %>%
  ungroup()

#join back together so each observation has its lags
data2_sam <- data_noNA_sam %>%
  left_join(Temp_lags_sam, by = c("SITE", "year", "month",
                              "Temp_gapfilled"))  

data2_sam_scale <-  data2_sam %>%
  # log-transform C:N data
  mutate(lnCN = log(mean_CN)) %>%
  # pivot to display temperatures long-format
  pivot_longer(c(Temp_gapfilled,mean_Temp_Cl1:mean_Temp_Cl12),
               names_to = "lag",
               values_to = "temp") %>%
  # scale all temps to each other
  mutate(temp_scaled = scale(temp)) %>%
  # pivot back to wide-format
  pivot_wider(names_from = lag,
              values_from = temp_scaled)

# Separate out temperature data for easier compiling for the model
data2_temp_sam <- data2_sam_scale %>%
  dplyr::select(Temp_gapfilled:mean_Temp_Cl12)

# Prep data for STAN as a list
data_stan_sam <- list(
  #indices for loops
  nrows = nrow(data2_sam_scale),
  nlag = 13, # 12 months into the past
  #reponse data
  CN = data2_sam_scale$lnCN,
  #temperature data
  Temp = data2_temp_sam)

saveRDS(data_stan_sam, here("model_summaries",
                        "STAN_input_list.RDS"))

# Fit SAM model.
stan_sam_run <- stan(file = "models/SAM_STAN_template.stan",
                      data = data_stan_sam,
                      chains = 3,
                      iter = 5000,
                      control = list(max_treedepth = 12))

# End of script.
