library(here)
library(tidyverse)


data <- readRDS(here('dat_working',
                     "CN_temp_monthly_20230720.rds"))


#LAGS
#eventually:
#kelp lags == avg. lifespan of a kelp plant
#for now:
#12 months of lags