# Prepping data and Running simple STAN SAM model for kelp nutrients
# Heili Lowman
# January 9, 2024

# Note - for extra reference when working between JAGS and STAN, please
# find the link to both user manuals below:

# JAGS: https://people.stat.sc.edu/hansont/stat740/jags_user_manual.pdf
# STAN: https://mc-stan.org/docs/stan-users-guide/index.html

#### Setup ####

# Load packages.
library(here)
library(tidyverse)
library(rstan)
library(shinystan)
library(tidybayes)
library(brms)
library(bayesplot)
library(modelsummary)
library(data.table)
library(naniar)

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

#### SAM Model Fit 12mo ####

# From Ana's JAGS script, need to make the necessary data transformations.
# And going to replace NAs with mean values for the time being.
data_noNA_sam <- data %>%
  # replace NAs with mean
  mutate(CN_gapfilled = replace_na(mean_CN, mean(mean_CN, na.rm = TRUE))) %>%
  mutate(Temp_gapfilled = replace_na(mean_Temp_C, mean(mean_Temp_C, na.rm = TRUE)))

# OPTION 2: Replace w/ values from normal distribution outside the model script,
# as Ana did in her model script. rnorm(mu, sigma)
# Could also consider differences among sites/within years.
cn_dist <- rnorm(100, 
                 mean(data$mean_CN, na.rm = TRUE),
                 sd(data$mean_CN, na.rm = TRUE))

temp_dist <- rnorm(100,
                   mean(data$mean_Temp_C, na.rm = TRUE),
                   sd(data$mean_CN, na.rm = TRUE))

data_repNA_sam <- data %>%
  # replace NAs with randomly generated values from normal dist.
  mutate(CN_gapfilled = replace_na(mean_CN, sample(cn_dist, 1))) %>%
  mutate(Temp_gapfilled = replace_na(mean_Temp_C, sample(temp_dist,1)))

# Create new dataset of lags using gapfilled data.
Temp_lags_sam <- data_repNA_sam %>%
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
data2_sam <- data_repNA_sam %>%
  left_join(Temp_lags_sam, by = c("SITE", "year", "month",
                              "Temp_gapfilled"))  

saveRDS(data2_sam, here("dat_working", "SAM_structured_input_data_STAN.RDS"))

# Assemble data for model

# Response data

# base dataset creation
# will make it easier to omit first year of data from CN and Temp data
data2_sam_scaled <- data2_sam %>%
  pivot_longer(Temp_gapfilled:mean_Temp_Cl12,
               names_to = "lag",
               values_to = "temp") %>%
  # scale all temps to each other
  mutate(temp = scale(temp)) %>%
  pivot_wider(names_from = lag,
              values_from = temp) %>%
  # remove first year of data from all three sites for which
  # we do not have proper lagged temperature data available
  drop_na(mean_Temp_Cl1:mean_Temp_Cl12)

# nutritional content (C:N) data
y <- as.vector(log(data2_sam_scaled$CN_gapfilled))
  
# temperature lag matrix 
Temp <- data2_sam_scaled %>%
  dplyr::select(Temp_gapfilled:mean_Temp_Cl12) %>%
  #make matrix
  as.matrix()

# Loop indices
# number of observations
nrows <- nrow(data2_sam_scaled)
# number of months included in each antecedent summation
nlag <- 13

# Prep data for STAN as a list
data_stan_sam <- list(
  #indices for loops
  nrows = nrows,
  nlag = nlag,
  #reponse data
  CN = y,
  #temperature data
  Temp = Temp,
  #prior
  alpha = c(rep(1/13, 13)))

saveRDS(data_stan_sam, here("model_summaries", "STAN_input_list.RDS"))

# Fit SAM model.
stan_sam_run <- stan(file = "models/SAM_STAN_template.stan",
                      data = data_stan_sam,
                      chains = 3,
                      iter = 5000,
                      control = list(max_treedepth = 12))

# Ran on Pinyon, started at 4:07PM.
# No signs of divergences, yay!

saveRDS(stan_sam_run, here("model_summaries", "STAN_output_list_011224.RDS"))

# Examine summaries of the estimates.
stan_sam_summary <- summary(stan_sam_run,
                          pars = c("b0", "b","wA[1]", "wA[2]", "wA[3]", 
                                   "wA[4]", "wA[5]", "wA[6]", 
                                   "wA[7]", "wA[8]", "wA[9]", "wA[10]", 
                                   "wA[11]", "wA[12]", "wA[13]", "sigma"),
                          probs = c(0.025, 0.5, 0.975))$summary # 2.5% and 97.5% percentiles
# Rhat values also look pretty good.

saveRDS(stan_sam_summary, here("model_summaries", "STAN_output_summary_011224.RDS"))

# Quick plot of model values.
plot(stan_sam_run, pars = c("b0", "b",
                      "wA[1]", "wA[2]", "wA[3]", "wA[4]", "wA[5]", "wA[6]", 
                      "wA[7]", "wA[8]", "wA[9]", "wA[10]", "wA[11]", "wA[12]",
                      "wA[13]", "sigma"))

# More fun plot.
color_scheme_set("brightblue")
mcmc_intervals(stan_sam_run,
           pars = c("b", "wA[1]", "wA[2]", "wA[3]", 
                    "wA[4]", "wA[5]", "wA[6]", 
                    "wA[7]", "wA[8]", "wA[9]", "wA[10]", 
                    "wA[11]", "wA[12]", "wA[13]"),
           point_est = "median",
           prob = 0.5,
           prob_outer = 0.95) +
  labs(
    title = "Posterior distributions",
    subtitle = "with medians and 95% intervals"
  )

#### SAM Model Fit 24mo ####

# Re-running much of the same code from the chunk above so these can
# stand on their own.

# From Ana's JAGS script, need to make the necessary data transformations.
# And going to replace NAs with mean values for the time being.
data_noNA_sam <- data %>%
  # replace NAs with mean
  mutate(CN_gapfilled = replace_na(mean_CN, mean(mean_CN, na.rm = TRUE))) %>%
  mutate(Temp_gapfilled = replace_na(mean_Temp_C, mean(mean_Temp_C, na.rm = TRUE)))

# OPTION 2: Replace w/ values from normal distribution outside the model script,
# as Ana did in her model script. rnorm(mu, sigma)
# Could also consider differences among sites/within years.
cn_dist <- rnorm(100, 
                 mean(data$mean_CN, na.rm = TRUE),
                 sd(data$mean_CN, na.rm = TRUE))

temp_dist <- rnorm(100,
                   mean(data$mean_Temp_C, na.rm = TRUE),
                   sd(data$mean_CN, na.rm = TRUE))

data_repNA_sam <- data %>%
  # replace NAs with randomly generated values from normal dist.
  mutate(CN_gapfilled = replace_na(mean_CN, sample(cn_dist, 1))) %>%
  mutate(Temp_gapfilled = replace_na(mean_Temp_C, sample(temp_dist,1)))

# Create new dataset of lags using gapfilled data.
Temp_lags_sam <- data_repNA_sam %>%
  dplyr::select(SITE, year, month, Temp_gapfilled) %>%
  group_by(SITE) %>%
  arrange(SITE, year, month) %>%
  #this creates a column for every lag 1:12 months ago
  do(data.frame(., setNames(shift(.$Temp_gapfilled, 1:24), c("mean_Temp_Cl1",
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
                                                             "mean_Temp_Cl12",
                                                             "mean_Temp_Cl13",
                                                             "mean_Temp_Cl14",
                                                             "mean_Temp_Cl15",
                                                             "mean_Temp_Cl16",
                                                             "mean_Temp_Cl17",
                                                             "mean_Temp_Cl18",
                                                             "mean_Temp_Cl19",
                                                             "mean_Temp_Cl20",
                                                             "mean_Temp_Cl21",
                                                             "mean_Temp_Cl22",
                                                             "mean_Temp_Cl23",
                                                             "mean_Temp_Cl24")))) %>%
  ungroup()

#join back together so each observation has its lags
data3_sam <- data_repNA_sam %>%
  left_join(Temp_lags_sam, by = c("SITE", "year", "month",
                                  "Temp_gapfilled"))  

saveRDS(data3_sam, here("dat_working", "SAM_structured_input_data_2ylag_STAN.RDS"))

# Assemble data for model

# Response data

# base dataset creation
# will make it easier to omit first two years of data from CN and Temp data
data3_sam_scaled <- data3_sam %>%
  pivot_longer(Temp_gapfilled:mean_Temp_Cl24,
               names_to = "lag",
               values_to = "temp") %>%
  # scale all temps to each other
  mutate(temp = scale(temp)) %>%
  pivot_wider(names_from = lag,
              values_from = temp) %>%
  # remove first two years of data from all three sites for which
  # we do not have proper lagged temperature data available for now
  drop_na(mean_Temp_Cl1:mean_Temp_Cl24)

# nutritional content (C:N) data
y <- as.vector(log(data3_sam_scaled$CN_gapfilled))

# temperature lag matrix 
Temp <- data3_sam_scaled %>%
  dplyr::select(Temp_gapfilled:mean_Temp_Cl24) %>%
  #make matrix
  as.matrix()

# Loop indices
# number of observations
nrows <- nrow(data3_sam_scaled)
# number of months included in each antecedent summation
nlag <- 25

# Prep data for STAN as a list
data_stan_sam <- list(
  #indices for loops
  nrows = nrows,
  nlag = nlag,
  #reponse data
  CN = y,
  #temperature data
  Temp = Temp,
  #prior
  alpha = c(rep(1/25, 25))) # prior is 0.04

saveRDS(data_stan_sam, here("model_summaries", "STAN_input_list_2ylag.RDS"))

# Fit SAM model.
stan_sam_run_2y <- stan(file = "models/SAM_STAN_template.stan",
                     data = data_stan_sam,
                     chains = 3,
                     iter = 5000,
                     control = list(max_treedepth = 12))

# Ran on Pinyon, started at 10:33AMPM, finished 
# 1823 divergences, eek!

saveRDS(stan_sam_run_2y, here("model_summaries", "STAN_output_list_011524.RDS"))

# Quick plot of model values.
plot(stan_sam_run_2y, pars = c("wA[1]", "wA[2]", "wA[3]", "wA[4]", "wA[5]", "wA[6]", 
                            "wA[7]", "wA[8]", "wA[9]", "wA[10]", "wA[11]", "wA[12]",
                            "wA[13]", "wA[14]", "wA[15]", "wA[16]", 
                            "wA[17]", "wA[18]", "wA[19]", "wA[20]", "wA[21]", 
                            "wA[22]", "wA[23]", "wA[24]", "wA[25]"))

# End of script.
