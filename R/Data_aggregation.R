# Data Aggregation Script
# July 20, 2023
# Heili Lowman

# This script will import the monthly bottom temperature dataset from the SBC LTER
# and combine with C:N data for use in SAM model structures.

# Load packages.
library(tidyverse)
library(here)
library(lubridate)
library(calecopal)

# Load datasets.

# Bottom temperature aggregated monthly (2002-2022)
dat_temp <- readRDS("dat_working/Bottom_temp_monthly_20230720.rds")

# Newest version of monthly kelp tissue C:N data (2002-2022)
dat_kelp <- read_csv("dat_raw/CHN_Field_all_years.csv")

# First, going to aggregate C:N data so there aren't multiple replicates per site
dat_kelp_agg <- dat_kelp %>%
  group_by(SITE, YEAR, MONTH) %>%
  summarize(mean_C = mean(C, na.rm = TRUE),
            mean_N = mean(N, na.rm = TRUE),
            mean_CN = mean(CN_RATIO, na.rm = TRUE)) %>%
  ungroup()

# Join datasets together.
dat_all <- full_join(dat_temp, dat_kelp_agg, by = c("SITE", "year" = "YEAR", 
                                                "month" = "MONTH")) %>%
  select(SITE, year, month, day, mean_C, mean_N, mean_CN, mean_Temp_C)

# And export for use in modeling.
saveRDS(dat_all,"dat_working/CN_temp_monthly_20230720.rds")

# End of script.
