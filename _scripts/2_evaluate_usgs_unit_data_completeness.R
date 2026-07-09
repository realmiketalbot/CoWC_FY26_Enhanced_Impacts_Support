#!/usr/bin/env Rscript

# ============================================================================
# Evaluate the completeness of USGS instantaneous stream gage data
# ============================================================================
#
# KEY FEATURES:
#   - Reads instantaneous discharge data organized by year
#   - Processes year-by-year to manage memory for large datasets
#   - Summarizes record length and completeness
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
usgs_unit_base_dir <- "/davenport-scratch/DATA/USGS/raw/Discharge/Unit"

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
  distinct() %>%
  mutate(chunk = ntile(row_number(), num_chunks)) %>%
  filter(chunk == chunk_index) %>%
  pull(GAGE_ID)

message(sprintf("Processing %d Colorado watersheds (chunk %d of %d)", 
                length(basin_meta), chunk_index, num_chunks))

# ============================================================================
# FIND AVAILABLE YEARS AND SITE DIRECTORIES
# ============================================================================
# Get list of site directories (USGS site numbers)
site_dirs <- list.dirs(usgs_unit_base_dir, recursive = FALSE, full.names = FALSE)

# Filter to Colorado watersheds
site_dirs <- site_dirs[site_dirs %in% basin_meta]

message(sprintf("Found %d Colorado site directories with instantaneous data", 
                length(site_dirs)))

# Get available years
available_years <- NULL
for (site_dir in site_dirs) {
  site_path <- file.path(usgs_unit_base_dir, site_dir)
  year_files <- list.files(site_path, pattern = "^\\d{4}\\.rds$")
  years <- as.integer(gsub("\\.rds$", "", year_files))
  available_years <- c(available_years, years)
}
available_years <- sort(unique(available_years))

message(sprintf("Year range: %d to %d", min(available_years), max(available_years)))

# ============================================================================
# INITIALIZE SUMMARY DATA STRUCTURES
# ============================================================================
# Store completeness info: site_no, year, count of observations, count of "A" quality
summary_by_year <- tibble(
  site_no = character(),
  year = integer(),
  total_obs = integer(),
  quality_a_obs = integer(),
  valid_fxn = numeric()
)

summary_overall <- tibble(
  site_no = character(),
  total_obs = integer(),
  quality_a_obs = integer(),
  valid_fxn_overall = numeric(),
  max_year = integer()
)

# ============================================================================
# PROCESS DATA YEAR BY YEAR (MEMORY EFFICIENT)
# ============================================================================
for (year in available_years) {
  message(sprintf("Processing year %d...", year))
  
  year_data <- NULL
  
  # Load data for this year from all sites
  for (site_no in site_dirs) {
    year_file <- file.path(usgs_unit_base_dir, site_no, sprintf("%d.rds", year))
    
    if (!file.exists(year_file)) next
    
    tryCatch({
      site_year_data <- readRDS(year_file) %>%
        as_tibble() %>%
        mutate(
          site_no = as.character(site_no),
          year = year,
          date = as.Date(dateTime),
          quality_flag = str_sub(X_00060_00000_cd, 1, 1),
          q_cfs = X_00060_00000
        ) %>%
        select(site_no, year, date, quality_flag, q_cfs)
      
      year_data <- bind_rows(year_data, site_year_data)
      
    }, error = function(e) {
      warning(sprintf("Error reading %s: %s", year_file, e$message))
    })
  }
  
  if (is.null(year_data) || nrow(year_data) == 0) {
    message(sprintf("  No data found for year %d", year))
    next
  }
  
  # Summarize by site and year
  year_summary <- year_data %>%
    group_by(site_no, year) %>%
    summarize(
      total_obs = n(),
      quality_a_obs = sum(quality_flag == "A", na.rm = TRUE),
      valid_fxn = quality_a_obs / total_obs,
      .groups = "drop"
    )
  
  summary_by_year <- bind_rows(summary_by_year, year_summary)
  
  # Clean up memory
  rm(year_data, year_summary)
  gc()
}

message("Calculating overall statistics...")

# Calculate overall statistics
summary_overall <- summary_by_year %>%
  group_by(site_no) %>%
  summarize(
    total_obs = sum(total_obs),
    quality_a_obs = sum(quality_a_obs),
    valid_fxn_overall = quality_a_obs / total_obs,
    max_year = max(year),
    .groups = "drop"
  )

# ============================================================================
# PREPARE DATA FOR PLOTTING
# ============================================================================
# Sort by: 2023 completeness (descending), then overall completeness (descending)
data_summary_by_year <- summary_by_year %>%
  left_join(summary_overall %>% select(site_no, valid_fxn_overall, max_year), 
            by = "site_no") %>%
  mutate(
    valid_fxn_2023 = ifelse(year == 2023, valid_fxn, NA),
    has_2023 = max_year >= 2023
  ) %>%
  arrange(desc(has_2023), desc(valid_fxn_overall), desc(valid_fxn))

site_levels <- data_summary_by_year %>%
  distinct(site_no) %>%
  pull(site_no)

data_summary_by_year <- data_summary_by_year %>%
  mutate(site_no = factor(site_no, levels = site_levels))

# ============================================================================
# CREATE HEATMAP PLOT
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
  ggtitle("USGS Streamgage Unit Record Completeness") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 3, angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave("_scripts/USGS_unit_record_completeness_heatmap.png", dpi = 300, width = 25, height = 4, units = "in")

message("Heatmap saved")

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
    description = "Watersheds with 40+ years of completely continuous data (100% hourly coverage)",
    value = length(has_years_with_threshold(summary_by_year, min_years = 40, min_completeness = 1.0))
  ),

  # 40+ years thresholds (>90% complete)
  tibble(
    metric = "Data for 40+ years (>90%)",
    description = "Watersheds with 40+ years of data (90% hourly coverage)",
    value = length(has_years_with_threshold(summary_by_year, min_years = 40, min_completeness = 0.9))
  ),
  
  # 30+ years thresholds
  tibble(
    metric = "Complete data for 30+ years",
    description = "Watersheds with 30+ years of completely continuous data (100% hourly coverage)",
    value = length(has_years_with_threshold(summary_by_year, min_years = 30, min_completeness = 1.0))
  ),

  # 30+ years thresholds (>90% complete)
  tibble(
    metric = "Data for 30+ years (>90%)",
    description = "Watersheds with 30+ years of data (90% hourly coverage)",
    value = length(has_years_with_threshold(summary_by_year, min_years = 30, min_completeness = 0.9))
  ),
  
  # 20+ years thresholds
  tibble(
    metric = "Complete data for 20+ years",
    description = "Watersheds with 20+ years of completely continuous data (100% hourly coverage)",
    value = length(has_years_with_threshold(summary_by_year, min_years = 20, min_completeness = 1.0))
  ),

  # 20+ years thresholds (>90% complete)
  tibble(
    metric = "Data for 20+ years (>90%)",
    description = "Watersheds with 20+ years of data (90% hourly coverage)",
    value = length(has_years_with_threshold(summary_by_year, min_years = 20, min_completeness = 0.9))
  ),
  
  # 2023 data
  tibble(
    metric = "Complete data for 2023",
    description = "Watersheds with 100% complete hourly data in 2023",
    value = length(has_year_complete(summary_by_year, target_year = 2023, min_completeness = 1.0))
  ),
  
  # 2023 & 20+ consecutive years
  tibble(
    metric = "Complete data for 2023 & 20+ years consecutive",
    description = "Watersheds with complete 2023 data AND 20+ consecutive years of 100% complete hourly data",
    value = length(
      intersect(
        has_year_complete(summary_by_year, target_year = 2023, min_completeness = 1.0),
        has_consecutive_years(summary_by_year, min_years = 20, min_completeness = 1.0)
      )
    )
  ),
  
  # 2023 & 20+ consecutive years (>90% complete)
  tibble(
    metric = "Complete data for 2023 & 20+ years consecutive (>90%)",
    description = "Watersheds with complete 2023 data AND 20+ consecutive years with >90% hourly completeness",
    value = length(
      intersect(
        has_year_complete(summary_by_year, target_year = 2023, min_completeness = 1.0),
        has_consecutive_years(summary_by_year, min_years = 20, min_completeness = 0.90)
      )
    )
  ),
  
  # 2023 & 20+ years (gaps allowed)
  tibble(
    metric = "Complete data for 2023 & 20+ years (gaps allowed)",
    description = "Watersheds with complete 2023 data AND 20+ years of 100% complete hourly data (gaps allowed)",
    value = length(
      intersect(
        has_year_complete(summary_by_year, target_year = 2023, min_completeness = 1.0),
        has_years_with_threshold(summary_by_year, min_years = 20, min_completeness = 1.0)
      )
    )
  ),
  
  # 2023 & 20+ years (>90% complete, gaps allowed)
  tibble(
    metric = "Complete data for 2023 & 20+ years (>90%, gaps allowed)",
    description = "Watersheds with complete 2023 data AND 20+ years with >90% hourly completeness (gaps allowed)",
    value = length(
      intersect(
        has_year_complete(summary_by_year, target_year = 2023, min_completeness = 1.0),
        has_years_with_threshold(summary_by_year, min_years = 20, min_completeness = 0.90)
      )
    )
  ),
  
  # 2022 data
  tibble(
    metric = "Complete data for 2022",
    description = "Watersheds with 100% complete hourly data in 2022",
    value = length(has_year_complete(summary_by_year, target_year = 2022, min_completeness = 1.0))
  ),
  
  # 2022 & 20+ consecutive years
  tibble(
    metric = "Complete data for 2022 & 20+ years consecutive",
    description = "Watersheds with complete 2022 data AND 20+ consecutive years of 100% complete hourly data",
    value = length(
      intersect(
        has_year_complete(summary_by_year, target_year = 2022, min_completeness = 1.0),
        has_consecutive_years(summary_by_year, min_years = 20, min_completeness = 1.0)
      )
    )
  ),
  
  # 2022 & 20+ consecutive years (>90% complete)
  tibble(
    metric = "Complete data for 2022 & 20+ years consecutive (>90%)",
    description = "Watersheds with complete 2022 data AND 20+ consecutive years with >90% hourly completeness",
    value = length(
      intersect(
        has_year_complete(summary_by_year, target_year = 2022, min_completeness = 1.0),
        has_consecutive_years(summary_by_year, min_years = 20, min_completeness = 0.90)
      )
    )
  ),
  
  # 2022 & 20+ years (gaps allowed)
  tibble(
    metric = "Complete data for 2022 & 20+ years (gaps allowed)",
    description = "Watersheds with complete 2022 data AND 20+ years of 100% complete hourly data (gaps allowed)",
    value = length(
      intersect(
        has_year_complete(summary_by_year, target_year = 2022, min_completeness = 1.0),
        has_years_with_threshold(summary_by_year, min_years = 20, min_completeness = 1.0)
      )
    )
  ),
  
  # 2022 & 20+ years (>90% complete, gaps allowed)
  tibble(
    metric = "Complete data for 2022 & 20+ years (>90%, gaps allowed)",
    description = "Watersheds with complete 2022 data AND 20+ years with >90% hourly completeness (gaps allowed)",
    value = length(
      intersect(
        has_year_complete(summary_by_year, target_year = 2022, min_completeness = 1.0),
        has_years_with_threshold(summary_by_year, min_years = 20, min_completeness = 0.90)
      )
    )
  )
)

completeness_stats <- completeness_stats %>%
  mutate(metric = factor(metric, levels = rev(completeness_stats$metric)))

# ============================================================================
# CREATE SUMMARY STATISTICS PLOT
# ============================================================================
ggplot(completeness_stats) +
  geom_col(aes(x = value, y = metric)) +
  scale_x_continuous(name = str_glue("Number of Watersheds (out of {nrow(summary_overall)})")) +
  scale_y_discrete(name = NULL) +
  labs(subtitle = "USGS Streamgage Unit Record Completeness") +
  theme_bw()

ggsave("_scripts/USGS_unit_record_completeness_threshold_summary.png", dpi = 300, width = 8, height = 4, units = "in")

message("Summary statistics plot saved")

# ============================================================================
# DISPLAY RESULTS
# ============================================================================
print(completeness_stats)

message(sprintf("\nTotal watersheds analyzed: %d", nrow(summary_overall)))
message(sprintf("Total observations: %d", sum(summary_overall$total_obs)))
message(sprintf("Total quality 'A' observations: %d", sum(summary_overall$quality_a_obs)))

# ============================================================================
# FINAL SELECTION
# ============================================================================
# Filter to stations meeting "Complete data for 2022 & 20+ years consecutive"
selected_stations <- intersect(
  has_year_complete(summary_by_year, target_year = 2022, min_completeness = 1.0),
  has_consecutive_years(summary_by_year, min_years = 20, min_completeness = 1.0)
)

message(sprintf("Selected %d stations", length(selected_stations)))

# Write to CSV
write_csv(
  tibble(site_no = selected_stations),
  "_data/selected_unit_stations.csv"
)

message("Station IDs written to _data/selected_unit_stations.csv")