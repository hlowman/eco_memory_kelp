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


# Get summary and plot some stuff -----------------------------------------

sum <- summary(model$samples)


#TEMPERUTARE overall effect:
#get mean value for temperature effect
mean <- as.data.frame(sum$statistics) %>%
  rownames_to_column(var = "parameter") %>%
  #select only temp variable
  filter(parameter == "b") %>%
  dplyr::select(Mean) %>%
  as_vector()

#get temperature effect samples
samps <- model$sims.list[[2]]

#make a ggplot of these values
temp_graph <- ggplot() +
  geom_density(aes(x = samps), fill = "#e34a33", size = 1) +
  theme_bw() +
  labs(x = "Temperature effect") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 15)) +
  geom_vline(xintercept = mean, linetype = 2, size = 1) 


#export this figure
ggsave(filename = here("figures",
                       "temp_effect.png"),
       plot = temp_graph,
       height = 6.5,
       width = 6.5,
       units = "in")

#get weights summarised
sum1 <- as.data.frame(sum$quantiles) %>%
  rownames_to_column(var = "parameter") %>%
  filter(str_detect(parameter, "wA")) %>%
  #get the lags out of the parameter names
  mutate(lag = str_sub(parameter, 
                       start = 4, 
                       end = nchar(parameter) - 1)) %>%
  #clean up parameter names
  mutate(parameter = str_sub(parameter, 
                             start = 1, 
                             end= 2)) %>%
  #re-level lag as factor
  mutate(lag = factor(lag, levels = c("1", '2',
                                      '3', '4',
                                      '5', '6',
                                      '7', '8',
                                      '9', '10',
                                      "11", "12",
                                      "13")))

#get proportion if all is equal
prop <- 1/13

#plot weights
weights_graph <- ggplot(sum1) +
  #if all equal - they should fall across this line
  geom_hline(yintercept = prop, linetype = 2, size = 1) +
  geom_pointrange(aes(x = lag,
                     y = `50%`,
                     ymin = `2.5%`,
                     ymax = `97.5%`), size = 1) +
  #get so that 1 = this month (0 months ago)
  scale_x_discrete(labels = c("0", "1", "2",
                              "3", "4", "5",
                              "6", "7", "8",
                              "9" , "10" ,"11",
                              "12")) +
  labs(x = "Months into the past",
       y = "Temperature importance weight \n (Median and 95% BCI)") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
      axis.title = element_text(size = 15)) 

#save this file
ggsave(filename = here("figures",
                       "temp_weights.png"),
       plot = weights_graph,
       height = 6.5,
       width = 6.5,
       units = "in")


# Export overall summary --------------------------------------------------

sum_tab <- as.data.frame(sum$quantiles) %>%
  rownames_to_column(var = "parameter") %>%
  filter(parameter != "deviance")

write.csv(sum_tab, here("model_summaries",
                        "SAM_model_summary_9_6_23.csv"))
