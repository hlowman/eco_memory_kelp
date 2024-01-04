# Prepping data and Running simple JAGS SAM model for kelp nutrients
# Ana Miller-ter Kuile
# September 6, 2023


# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse",
                  'data.table',
                  'jagsUI',
                  'mcmcplots')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Read data ---------------------------------------------------------------


data <- readRDS(here('dat_working',
                     "CN_temp_monthly_20230720.rds"))


#eventually - get Heili to pull in lags for the data that "starts" the time series
# e.g. are there temp data for prior to the first observations
#for 12 months prior? or length we decide to go with?



# Prep temp data ----------------------------------------------------------



Temp_lags <- data %>%
  dplyr::select(SITE, year, month, mean_Temp_C) %>%
  group_by(SITE) %>%
  arrange(SITE, year, month) %>%
  #this creates a column for every lag 1:12 months ago
  do(data.frame(., setNames(shift(.$mean_Temp_C, 1:12), c("mean_Temp_Cl1",
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
data2 <- data %>%
  left_join(Temp_lags, by = c("SITE", "year", "month",
                              "mean_Temp_C"))

saveRDS(data2, here("dat_working", "SAM_structured_input_data.RDS"))

# Data objects for model --------------------------------------------------

#Loop indices
#i loop
nrows <- nrow(data)
#number of covariates for priors
n.covs <- 1
#number of lags for weight priors
n.lag <- 13

#Response data

#look at distribution to see if transformation is necessary
hist(data$mean_CN)
hist(log(data$mean_CN))
#will want to log transform 
y <- as.vector(log(data$mean_CN))

#temperature lag matrix 
Temp <- data2 %>%
  pivot_longer(mean_Temp_C:mean_Temp_Cl12,
               names_to = "lag",
               values_to = "temp") %>%
  #scale all temps to each other
  mutate(temp = scale(temp)) %>%
  pivot_wider(names_from = lag,
              values_from = temp) %>%
  dplyr::select(mean_Temp_C:mean_Temp_Cl12) %>%
  #make matrix
  as.matrix()

#LAGS
#eventually:
#kelp lags == avg. lifespan of a kelp plant
#for now:
#12 months of lags


# List data ---------------------------------------------------------------

data_list <- list(#indicies for loops
                  nrows = nrows,
                  n.covs = n.covs,
                  n.lag = n.lag,
                  #reponse data
                  y = y,
                  #temperature data
                  Temp = Temp)

saveRDS(data_list, here("model_summaries",
                        "JAGS_input_list.RDS"))


# Run JAGS model ----------------------------------------------------------

#Get parameters to track
parms <- c("b0", 
           "b",
           'wA')

#load model
model_file <- here("models",
              "SAM_template.R")

#Run the jags model
model <- jagsUI::jags(data = data_list,
                      inits = NULL,
                      parameters.to.save = parms,
                      model.file = model_file,
                      n.chains = 3,
                      n.iter = 4000,
                      parallel = TRUE)


# Check convergence -------------------------------------------------------

mcmcplot(model$samples)



# Get summaries and samples -----------------------------------------------

sum <- summary(model$samples)

saveRDS(sum, here("model_summaries", "model_summary_list.RDS"))

sims <- model$sims.list

saveRDS(sims, here("model_summaries", "model_samples.RDS"))
