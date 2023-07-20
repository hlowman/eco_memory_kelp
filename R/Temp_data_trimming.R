# Temperature Trimming Script
# July 17, 2023
# Heili Lowman

# This script will import the raw bottom temperature dataset from the SBC LTER
# and trim to locations of interest so the dataset will be small enough to load
# onto github.

# Load packages.
library(tidyverse)
library(here)
library(lubridate)
library(calecopal)

# Load dataset.
dat <- read_csv("dat_raw/Bottom_temp_all_years_20220729.csv")

# Trim dataset down to ABUR, AQUE, and MOHK.
dat_trim <- dat %>%
  filter(SITE %in% c("ABUR", "AQUE", "MOHK"))

# Export 15 minute frequency dataset.
saveRDS(dat_trim,"dat_raw/Bottom_temp_ABUR_AQUE_MOHK_20230717.rds")

# Add new columns designating year and month.
dat_trim <- dat_trim %>%
  mutate(year = year(DATE_LOCAL),
         month = month(DATE_LOCAL))

# Aggregate on a monthly basis.
dat_monthly <- dat_trim %>%
  group_by(SITE, year, month) %>%
  summarize(mean_Temp_C = mean(TEMP_C, na.rm = TRUE)) %>%
  ungroup() %>%
  # and create new columns for "dates"
  mutate(day = 1) %>%
  mutate(date = make_date(year, month, day))

# Quick plot to be sure they're all there.
ggplot(dat_monthly, aes(x = date, y = mean_Temp_C)) +
  geom_line(aes(color = SITE), alpha = 0.75, size = 2) +
  scale_color_manual(values = cal_palette("kelp2")) +
  labs(x = "Date", y = "Temperature (deg C)") +
  theme_bw() # ok looks alright

# And export for use in modeling.
saveRDS(dat_monthly,"dat_working/Bottom_temp_monthly_20230720.rds")

# End of script.
