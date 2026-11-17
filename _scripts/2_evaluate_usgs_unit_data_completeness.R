#!/usr/bin/env Rscript

# ============================================================================
# Script 2: Measure record coverage of USGS instantaneous (unit) discharge
#           data for all Colorado GAGES-II sites
# ============================================================================
#
# Writes a site x water year x month coverage table that Script 2b uses to
# select stations, and that later defines the effective record length used
# for the Poisson exceedance rate and the dispersion index.
#
# WHAT CHANGED (2026-11-17, milestone M1 re-screen):
# The previous version measured "completeness" as
#     valid_fxn = (observations flagged "A") / (observations present)
# on calendar years. Missing data never entered that ratio, so a year with
# only a few months of approved data scored 100%. The station rule built
# on it ("100% in 2022 and 20+ consecutive 100% years") was an approval
# screen, not a coverage screen. This version measures coverage directly,
# and keeps the approval share as a separate column. Station selection moved
# to Script 2b so the criteria can be tuned without re-reading unit data.
#
# DEFINITION OF AN OBSERVED DAY
# A day is observed if it has at least one non-missing, non-negative
# discharge value. This is exactly the rule Script 3 uses to build the
# daily-max series, so coverage here matches the series peaks come from.
#
# INPUTS:
#   - GAGES-II conterm_basinid.txt (Colorado sites)
#   - /davenport-scratch/DATA/USGS/raw/Discharge/Unit/<site_no>/<year>.rds
#
# OUTPUTS:
#   - _data/usgs_unit_coverage_co.rds : one row per site x water year x month
#       over each site's span (first to last water year with any data),
#       including months with no data. Columns:
#         site_no, water_year, month, cal_year,
#         n_days_obs     : observed days in the month
#         days_in_month
#         frac_obs       : n_days_obs / days_in_month
#         n_obs          : valid unit values in the month
#         n_obs_approved : of those, values flagged "A" (approved)
#   - _scripts/USGS_unit_record_coverage_heatmap.png : site x water year
#
# Runtime: reads every year file once per site. Sites are independent, so
# this parallelises with furrr when N_WORKERS > 1.
# ============================================================================

library(tidyverse)
library(colorspace)
library(furrr)

# ============================================================================
# CONFIGURATION
# ============================================================================
watershed_meta     <- "/davenport-scratch/DATA/GAGES-II/basinchar_and_report_sept_2011/spreadsheets-in-csv-format/conterm_basinid.txt"
usgs_unit_base_dir <- "/davenport-scratch/DATA/USGS/raw/Discharge/Unit"
output_file        <- "_data/usgs_unit_coverage_co.rds"
heatmap_file       <- "_scripts/USGS_unit_record_coverage_heatmap.png"

N_WORKERS <- max(1L, min(parallelly::availableCores() - 1L, 16L))

# ============================================================================
# CANDIDATE SITES: Colorado GAGES-II with unit-data directories
# ============================================================================
co_sites <- read_csv(watershed_meta, show_col_types = FALSE) %>%
  filter(STATE == "CO") %>%
  distinct(STAID) %>%
  pull(STAID) %>%
  as.character()

site_dirs <- list.dirs(usgs_unit_base_dir, recursive = FALSE, full.names = FALSE)
candidates <- sort(intersect(site_dirs, co_sites))

message(sprintf("Colorado GAGES-II sites: %d | with unit-data directories: %d",
                length(co_sites), length(candidates)))

# ============================================================================
# PER-SITE COVERAGE
# ============================================================================
site_coverage <- function(site_no) {

  year_files <- list.files(file.path(usgs_unit_base_dir, site_no),
                           pattern = "^\\d{4}\\.rds$", full.names = TRUE)
  if (length(year_files) == 0) return(NULL)

  daily <- map_dfr(year_files, function(f) {
    tryCatch({
      readRDS(f) %>%
        as_tibble() %>%
        transmute(
          date     = as.Date(dateTime),
          q_cfs    = X_00060_00000,
          approved = str_sub(X_00060_00000_cd, 1, 1) == "A"
        ) %>%
        filter(!is.na(q_cfs), q_cfs >= 0) %>%          # same rule as Script 3
        group_by(date) %>%
        summarize(n_obs = n(), n_obs_approved = sum(approved, na.rm = TRUE),
                  .groups = "drop")
    }, error = function(e) {
      warning(sprintf("Error reading %s: %s", f, e$message))
      NULL
    })
  })

  if (is.null(daily) || nrow(daily) == 0) return(NULL)

  # A date can appear in two year files if files overlap; keep one row per date
  daily <- daily %>%
    group_by(date) %>%
    summarize(n_obs = sum(n_obs), n_obs_approved = sum(n_obs_approved), .groups = "drop")

  daily %>%
    mutate(
      month      = as.integer(lubridate::month(date)),
      water_year = as.integer(lubridate::year(date)) + as.integer(month >= 10L)
    ) %>%
    group_by(water_year, month) %>%
    summarize(n_days_obs     = n(),
              n_obs          = sum(n_obs),
              n_obs_approved = sum(n_obs_approved),
              .groups = "drop") %>%
    complete(water_year = full_seq(water_year, 1), month = 1:12,
             fill = list(n_days_obs = 0L, n_obs = 0L, n_obs_approved = 0L)) %>%
    mutate(
      site_no       = site_no,
      water_year    = as.integer(water_year),
      cal_year      = water_year - as.integer(month >= 10L),
      days_in_month = lubridate::days_in_month(
                        as.Date(sprintf("%d-%02d-01", cal_year, month))),
      frac_obs      = n_days_obs / days_in_month
    )
}

message(sprintf("Measuring coverage with %d worker(s)...", N_WORKERS))
if (N_WORKERS > 1) plan(multisession, workers = N_WORKERS) else plan(sequential)

coverage <- future_map(candidates, site_coverage, .progress = TRUE) %>%
  bind_rows() %>%
  select(site_no, water_year, month, cal_year, n_days_obs, days_in_month,
         frac_obs, n_obs, n_obs_approved) %>%
  arrange(site_no, water_year, cal_year, month)

plan(sequential)

saveRDS(coverage, output_file)
message(sprintf("Saved %d rows for %d sites to %s",
                nrow(coverage), n_distinct(coverage$site_no), output_file))

# ============================================================================
# DIAGNOSTICS
# ============================================================================
wy_cov <- coverage %>%
  group_by(site_no, water_year) %>%
  summarize(frac_obs = sum(n_days_obs) / sum(days_in_month), .groups = "drop")

# Mean coverage by calendar month across all site-years: where do gaps fall?
by_month <- coverage %>%
  group_by(month) %>%
  summarize(mean_frac_obs = mean(frac_obs), .groups = "drop")
message("\nMean coverage by calendar month (all sites, all years in span):")
message(paste(sprintf("  %s: %5.1f%%", month.abb[by_month$month],
                      100 * by_month$mean_frac_obs), collapse = "\n"))

site_order <- wy_cov %>%
  group_by(site_no) %>%
  summarize(mean_cov = mean(frac_obs), .groups = "drop") %>%
  arrange(desc(mean_cov)) %>%
  pull(site_no)

ggplot(wy_cov %>% mutate(site_no = factor(site_no, levels = site_order))) +
  geom_raster(aes(x = site_no, y = water_year, fill = frac_obs)) +
  scale_x_discrete(name = "Site ID", expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(name = "Water year", expand = expansion(mult = c(0, 0))) +
  scale_fill_continuous_sequential(name = "Days observed", palette = "RdPu",
                                   labels = scales::percent, limits = c(0, 1)) +
  ggtitle("USGS unit discharge record coverage (observed days / days in water year)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 3, angle = 90, hjust = 1, vjust = 0.5))

ggsave(heatmap_file, dpi = 300, width = 25, height = 4, units = "in")
message(sprintf("Heatmap saved to %s", heatmap_file))
message("\nScript 2 complete. Next: Script 2b (station selection).")
