#SAM model pictures and summaries
# Ana Miller-ter Kuile
# September 6, 2023


# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

sum <- readRDS(here("model_summaries", "model_summary_list.RDS"))

sims <- readRDS(here("model_summaries", "model_samples.RDS"))

data_list <- readRDS(here("model_summaries",
                        "JAGS_input_list.RDS"))

data2 <- readRDS(here("dat_working", "SAM_structured_input_data.RDS"))

# Temperature overall effect ----------------------------------------------

#TEMPERUTARE overall effect:
#get mean value for temperature effect
mean <- as.data.frame(sum$statistics) %>%
  rownames_to_column(var = "parameter") %>%
  #select only temp variable
  filter(parameter == "b") %>%
  dplyr::select(Mean) %>%
  as_vector()

#get temperature effect samples
samps <- sims[[2]]

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


# Temp partial plot -------------------------------------------------------


#pull out medians of weights as vector
#make climate a matrix of temp lags
#multiply each by the weight and take the sum for each row (observation)
#it would be a "monthly weighted z-score"
#beta for temp
blT <- as.data.frame(sum$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm == "b") %>%
  dplyr::select(`50%`) %>%
  as_vector()

#b0
b0 <- as.data.frame(sum$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm == "b0") %>%
  dplyr::select(`50%`) %>%
  as_vector()

#get temparutres on scaled scale
#use Temp object from model data_list

#make scaled data long format to get mean and sd
tscale <- data2 %>%
  dplyr::select(SITE, year, month, mean_Temp_C:mean_Temp_Cl12) %>% #adjust if needed
  pivot_longer(mean_Temp_C:mean_Temp_Cl12,
               names_to = "lag",
               values_to = "temp") 

#get mean and SD of OG data to back-transform stuff
mean <- mean(tscale$temp, na.rm = T)
sd <- sd(tscale$temp, na.rm = T)

#get weights per month
t_wt <- as.data.frame(sum$quantiles) %>%
  rownames_to_column(var = "parameter") %>%
  filter(str_detect(parameter, "wA")) %>%
  dplyr::select(`50%`) %>%
  as_vector()

#get tmax dataset
regT <- data2 %>%
  dplyr::select(SITE, year, month, mean_CN, mean_Temp_C:mean_Temp_Cl12)

#multiply months by their weights
regT$TAnt <- apply(data_list$Temp, MARGIN = 1, FUN = function(x){sum(x*t_wt)})

#revert Tmax to OG data scale
regT <- regT %>%
  mutate(Temp = TAnt*sd + mean)

#regression prediction for Temperature
regT <- regT %>%
  mutate(reg = b0 + blT*TAnt,
         exp_reg = exp(reg))
# xscaled = (x – x̄) / s
# xscaled*sd + mean = x
xlab <-  expression("Weighted temperature " ( degree*C))

temp_pred <- ggplot(regT) +
  geom_point(aes(x = Temp, y = mean_CN),
             fill = "#e34a33", color = "black", shape = 21, alpha = 0.5) +
  geom_line(aes(x = Temp, y = exp_reg),
            lwd = 2, color = "black") +
  geom_line(aes(x = Temp, y = exp_reg), 
            color = "#e34a33", lwd = 1.5) +
  labs(x = xlab,
       y = "Kelp C:N") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15)) 


#export this figure
ggsave(filename = here("figures",
                       "temp_predicted.png"),
       plot = temp_pred,
       height = 6.5,
       width = 6.5,
       units = "in")

# Weights -----------------------------------------------------------------


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

