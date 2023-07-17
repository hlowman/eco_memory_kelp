# Temperature Trimming Script
# July 17, 2023
# Heili Lowman

# This script will import the raw bottom temperature dataset from the SBC LTER
# and trim to locations of interest so the dataset will be small enough to load
# onto github.

# Load packages.
library(tidyverse)
library(here)

# Load dataset.
dat <- read_csv("dat_raw/Bottom_temp_all_years_20220729.csv")

# Trim dataset down to ABUR, AQUE, and MOHK.
dat_trim <- dat %>%
  filter(SITE %in% c("ABUR", "AQUE", "MOHK"))

# Export dataset.
saveRDS(dat_trim,"dat_raw/Bottom_temp_ABUR_AQUE_MOHK_20230717.rds")

# End of script.
