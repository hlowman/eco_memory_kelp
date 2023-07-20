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
library(patchwork)

# Load datasets.

# Bottom temperature aggregated monthly (2002-2022)
dat_temp <- readRDS("dat_working/Bottom_temp_monthly_20230720.rds")

# Newest version of monthly kelp tissue C:N data (2002-2022)
dat_kelp <- read_csv("dat_raw/CHN_Field_all_years.csv")

# First, going to aggregate C:N data so there aren't multiple replicates per site
dat_kelp_agg <- dat_kelp %>%
  # and need to remove -99999 values
  filter(CN_RATIO > 0) %>%
  # and need to rename "MORK" sites as "MOHK"
  mutate(site = case_when(SITE == "MORK" ~ "MOHK",
                          TRUE ~ SITE)) %>%
  group_by(site, YEAR, MONTH) %>%
  summarize(mean_C = mean(C, na.rm = TRUE),
            mean_N = mean(N, na.rm = TRUE),
            mean_CN = mean(CN_RATIO, na.rm = TRUE)) %>%
  ungroup()

# Join datasets together.
dat_agg <- full_join(dat_temp, dat_kelp_agg, by = c("SITE" = "site", "year" = "YEAR", 
                                                "month" = "MONTH")) 

dat_all <- dat_agg %>%
  select(SITE, year, month, day, mean_C, mean_N, mean_CN, mean_Temp_C)

# And export for use in modeling.
saveRDS(dat_all,"dat_working/CN_temp_monthly_20230720.rds")

# And making some quick plots for the README of the github page.
(fig1 <- ggplot(dat_agg, aes(x = date, y = mean_CN)) +
    geom_line(aes(color = SITE), alpha = 0.75, size = 2) +
    scale_color_manual(values = cal_palette("figmtn")) +
    labs(x = "DATE", y = "C:N") +
    theme_bw())

(fig2 <- ggplot(dat_agg, aes(x = date, y = mean_Temp_C)) +
    geom_line(aes(color = SITE), alpha = 0.75, size = 2) +
    scale_color_manual(values = cal_palette("kelp2")) +
    labs(x = "DATE", y = "Temperature") +
    theme_bw())

(fig3 <- ggplot(dat_agg, aes(x = mean_Temp_C, y = mean_C)) +
    geom_point(aes(color = SITE), alpha = 0.75, size = 2) +
    scale_color_manual(values = cal_palette("figmtn")) +
    labs(x = "Temperature", y = "%C") +
    theme_bw() +
    theme(legend.position = "none"))

(fig4 <- ggplot(dat_agg, aes(x = mean_Temp_C, y = mean_N)) +
    geom_point(aes(color = SITE), alpha = 0.75, size = 2) +
    scale_color_manual(values = cal_palette("figmtn")) +
    labs(x = "Temperature", y = "%N") +
    theme_bw() +
    theme(legend.position = "none"))

(fig5 <- ggplot(dat_agg, aes(x = mean_Temp_C, y = mean_CN)) +
    geom_point(aes(color = SITE), alpha = 0.75, size = 2) +
    scale_color_manual(values = cal_palette("figmtn")) +
    labs(x = "Temperature", y = "C:N") +
    theme_bw())

(fig_comp <- fig1 /
  fig2 /
  (fig3 + fig4 + fig5))

# Export figure.
# ggsave(("figures/data_summary_figure.png"),
#        width = 25,
#        height = 20,
#        units = "cm"
# )

# End of script.
