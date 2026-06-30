#!/usr/bin/env Rscript

# ============================================================================
# Evaluate the completeness of USGS daily stream gage data
# ============================================================================
#
# KEY FEATURES:
#   - Reads all data for all GAGES-II watersheds in Colorado
#   - Summarizes record length and completeness, overall and by water year
#   - Prioritizes watersheds with data through 2023, sorted by completeness
# ============================================================================

library(tidyverse)
library(sf)
library(arrow)
library(colorspace)

# ============================================================================
# FILE PATHS AND CONFIGURATION
# ============================================================================
watershed_dataset <- "GAGES-II"
watershed_meta <- "/davenport-scratch/DATA/GAGES-II/basinchar_and_report_sept_2011/spreadsheets-in-csv-format/conterm_basinid.txt"
usgs_data_directory <- "/davenport-scratch/DATA/USGS/raw/Discharge/daily"

# Command line arguments for chunking (optional)
# Usage: Rscript script.R [chunk_index] [num_chunks]
args <- commandArgs(trailingOnly = TRUE)
chunk_index <- ifelse(length(args) >= 1, as.integer(args[1]), 1)
num_chunks <- ifelse(length(args) >= 2, as.integer(args[2]), 1)

# ============================================================================
# LOAD WATERSHEDS AND FILTER TO COLORADO
# ============================================================================
basin_meta <- read_csv(watershed_meta, show_col_types = FALSE) %>%
  filter(STATE == "CO") %>%
  select(GAGE_ID = STAID) %>%
  distinct()

# ============================================================================
# READ DAILY USGS DATA
# ============================================================================
all_usgs_files <- list.files(usgs_data_directory, pattern = ".rds")
gage_filter <- paste0(basin_meta$GAGE_ID, ".rds")
usgs_files <- file.path(
  usgs_data_directory,
  all_usgs_files[all_usgs_files %in% gage_filter]
)

# Merge data and filter to "A" quality flags
usgs_data <- lapply(usgs_files, readRDS) %>%
  bind_rows() %>%
  filter(str_sub(X_00060_00003_cd, 1, 1) == "A") %>%
  select(site_no, Date, X_00060_00003)

names(usgs_data) <- c("site_no", "date", "q_cfs")

site_nos <- unique(usgs_data$site_no)
dates <- unique(usgs_data$date)

full_grid <- expand.grid(
  site_no = site_nos, 
  date = dates
)

usgs_data <- merge(full_grid, usgs_data, by = c("site_no", "date"), all = T)

# ============================================================================
# CALCULATE SUMMARY STATISTICS
# ============================================================================
data_summary_overall <- usgs_data %>%
  mutate(year = year(date)) %>%
  group_by(site_no) %>%
  summarize(
    valid_fxn_overall = sum(!is.na(q_cfs)) / n()
  )

data_summary_by_year <- usgs_data %>%
  mutate(year = year(date)) %>%
  group_by(site_no, year) %>%
  summarize(
    valid_fxn = sum(!is.na(q_cfs)) / n(),
    .groups = "drop"
  ) %>%
  # Flag watersheds where 2023 is complete
  mutate(valid_fxn_2023 = ifelse(year == 2023, valid_fxn, NA)) %>%
  # Merge with overall summary
  merge(data_summary_overall, by = "site_no") %>%
  arrange(-valid_fxn_2023, -valid_fxn_overall, -valid_fxn)

site_levels <- data_summary_by_year %>%
  select(site_no) %>%
  distinct(site_no)

data_summary_by_year <- data_summary_by_year %>%
  mutate(site_no = factor(site_no, site_levels$site_no))

# ============================================================================
# CREATE PLOT
# ============================================================================
ggplot(data_summary_by_year) +
  geom_raster(aes(x = site_no, y = year, fill = valid_fxn)) +
  scale_x_discrete(
    name = "Site ID", 
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    name = "Year",
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_continuous_sequential(name = "Completeness", palette = "RdPu", labels=scales::percent) +
  ggtitle("USGS Streamgage Daily Record Completeness") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 3, angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave("_scripts/USGS_daily_record_completeness_heatmap.png", dpi = 300, width = 25, height = 4, units = "in")

# ============================================================================
# HELPER FUNCTIONS FOR COMPLETENESS THRESHOLDS
# ============================================================================

# Check if a watershed has data for a specific year with minimum completeness
has_year_complete <- function(df, target_year, min_completeness = 1.0) {
  df %>%
    filter(year == target_year, valid_fxn >= min_completeness) %>%
    pull(site_no) %>%
    unique()
}

# Check if a watershed has N+ consecutive years of minimum completeness
has_consecutive_years <- function(df, min_years = 20, min_completeness = 1.0) {
  df %>%
    filter(valid_fxn >= min_completeness) %>%
    group_by(site_no) %>%
    mutate(
      year_lag = year - lag(year, default = year[1] - 1),
      consecutive_group = cumsum(year_lag != 1)
    ) %>%
    group_by(site_no, consecutive_group) %>%
    mutate(consecutive_count = n()) %>%
    ungroup() %>%
    filter(consecutive_count >= min_years) %>%
    pull(site_no) %>%
    unique()
}

# Check if a watershed has N+ years of data with minimum completeness (gaps allowed)
has_years_with_threshold <- function(df, min_years = 20, min_completeness = 1.0) {
  df %>%
    filter(valid_fxn >= min_completeness) %>%
    group_by(site_no) %>%
    summarize(n_years = n(), .groups = "drop") %>%
    filter(n_years >= min_years) %>%
    pull(site_no) %>%
    unique()
}

# ============================================================================
# GENERATE COMPLETENESS SUMMARY STATISTICS
# ============================================================================
completeness_stats <- bind_rows(
  # 40+ years thresholds
  tibble(
    metric = "Complete data for 40+ years",
    description = "Watersheds with 40+ years of completely continuous data (100% daily coverage)",
    value = length(has_years_with_threshold(data_summary_by_year, min_years = 40, min_completeness = 1.0))
  ),

  # 40+ years thresholds (>90% complete)
  tibble(
    metric = "Data for 40+ years (>90%)",
    description = "Watersheds with 40+ years of data (90% daily coverage)",
    value = length(has_years_with_threshold(data_summary_by_year, min_years = 40, min_completeness = 0.9))
  ),
  
  # 30+ years thresholds
  tibble(
    metric = "Complete data for 30+ years",
    description = "Watersheds with 30+ years of completely continuous data (100% daily coverage)",
    value = length(has_years_with_threshold(data_summary_by_year, min_years = 30, min_completeness = 1.0))
  ),

  # 30+ years thresholds (>90% complete)
  tibble(
    metric = "Data for 30+ years (>90%)",
    description = "Watersheds with 30+ years of data (90% daily coverage)",
    value = length(has_years_with_threshold(data_summary_by_year, min_years = 30, min_completeness = 0.9))
  ),
  
  # 20+ years thresholds
  tibble(
    metric = "Complete data for 20+ years",
    description = "Watersheds with 20+ years of completely continuous data (100% daily coverage)",
    value = length(has_years_with_threshold(data_summary_by_year, min_years = 20, min_completeness = 1.0))
  ),

  # 30+ years thresholds (>90% complete)
  tibble(
    metric = "Data for 20+ years (>90%)",
    description = "Watersheds with 20+ years of data (90% daily coverage)",
    value = length(has_years_with_threshold(data_summary_by_year, min_years = 20, min_completeness = 0.9))
  ),
  
  # 2023 data
  tibble(
    metric = "Complete data for 2023",
    description = "Watersheds with 100% complete data in 2023",
    value = length(has_year_complete(data_summary_by_year, target_year = 2023, min_completeness = 1.0))
  ),
  
  # 2023 & 20+ consecutive years
  tibble(
    metric = "Complete data for 2023 & 20+ years consecutive",
    description = "Watersheds with complete 2023 data AND 20+ consecutive years of 100% complete data",
    value = length(
      intersect(
        has_year_complete(data_summary_by_year, target_year = 2023, min_completeness = 1.0),
        has_consecutive_years(data_summary_by_year, min_years = 20, min_completeness = 1.0)
      )
    )
  ),
  
  # 2023 & 20+ consecutive years (>90% complete)
  tibble(
    metric = "Complete data for 2023 & 20+ years consecutive (>90%)",
    description = "Watersheds with complete 2023 data AND 20+ consecutive years with >90% completeness",
    value = length(
      intersect(
        has_year_complete(data_summary_by_year, target_year = 2023, min_completeness = 1.0),
        has_consecutive_years(data_summary_by_year, min_years = 20, min_completeness = 0.90)
      )
    )
  ),
  
  # 2023 & 20+ years (gaps allowed)
  tibble(
    metric = "Complete data for 2023 & 20+ years (gaps allowed)",
    description = "Watersheds with complete 2023 data AND 20+ years of 100% complete data (gaps allowed)",
    value = length(
      intersect(
        has_year_complete(data_summary_by_year, target_year = 2023, min_completeness = 1.0),
        has_years_with_threshold(data_summary_by_year, min_years = 20, min_completeness = 1.0)
      )
    )
  ),
  
  # 2023 & 20+ years (>90% complete, gaps allowed)
  tibble(
    metric = "Complete data for 2023 & 20+ years (>90%, gaps allowed)",
    description = "Watersheds with complete 2023 data AND 20+ years with >90% completeness (gaps allowed)",
    value = length(
      intersect(
        has_year_complete(data_summary_by_year, target_year = 2023, min_completeness = 1.0),
        has_years_with_threshold(data_summary_by_year, min_years = 20, min_completeness = 0.90)
      )
    )
  ),
  
  # 2022 data
  tibble(
    metric = "Complete data for 2022",
    description = "Watersheds with 100% complete data in 2022",
    value = length(has_year_complete(data_summary_by_year, target_year = 2022, min_completeness = 1.0))
  ),
  
  # 2022 & 20+ consecutive years
  tibble(
    metric = "Complete data for 2022 & 20+ years consecutive",
    description = "Watersheds with complete 2022 data AND 20+ consecutive years of 100% complete data",
    value = length(
      intersect(
        has_year_complete(data_summary_by_year, target_year = 2022, min_completeness = 1.0),
        has_consecutive_years(data_summary_by_year, min_years = 20, min_completeness = 1.0)
      )
    )
  ),
  
  # 2022 & 20+ consecutive years (>90% complete)
  tibble(
    metric = "Complete data for 2022 & 20+ years consecutive (>90%)",
    description = "Watersheds with complete 2022 data AND 20+ consecutive years with >90% completeness",
    value = length(
      intersect(
        has_year_complete(data_summary_by_year, target_year = 2022, min_completeness = 1.0),
        has_consecutive_years(data_summary_by_year, min_years = 20, min_completeness = 0.90)
      )
    )
  ),
  
  # 2022 & 20+ years (gaps allowed)
  tibble(
    metric = "Complete data for 2022 & 20+ years (gaps allowed)",
    description = "Watersheds with complete 2022 data AND 20+ years of 100% complete data (gaps allowed)",
    value = length(
      intersect(
        has_year_complete(data_summary_by_year, target_year = 2022, min_completeness = 1.0),
        has_years_with_threshold(data_summary_by_year, min_years = 20, min_completeness = 1.0)
      )
    )
  ),
  
  # 2022 & 20+ years (>90% complete, gaps allowed)
  tibble(
    metric = "Complete data for 2022 & 20+ years (>90%, gaps allowed)",
    description = "Watersheds with complete 2022 data AND 20+ years with >90% completeness (gaps allowed)",
    value = length(
      intersect(
        has_year_complete(data_summary_by_year, target_year = 2022, min_completeness = 1.0),
        has_years_with_threshold(data_summary_by_year, min_years = 20, min_completeness = 0.90)
      )
    )
  )
)

completeness_stats <- completeness_stats %>%
  mutate(metric = factor(metric, levels = rev(completeness_stats$metric)))

ggplot(completeness_stats) +
  geom_col(aes(x = value, y = metric)) +
  scale_x_continuous(name = str_glue("Number of Watersheds (out of {nrow(data_summary_overall)})")) +
  scale_y_discrete(name = NULL) +
  labs(subtitle = "USGS Streamgage Daily Record Completeness") +
  theme_bw()

ggsave("_scripts/USGS_daily_record_completeness_threshold_summary.png", dpi = 300, width = 8, height = 4, units = "in")
