#!/usr/bin/env Rscript

# ============================================================================
# Script 6a: Flood Type Classification — Depth Threshold (Davenport et al. 2020 fixed-depth rules)
# ============================================================================
# Classifies each POT peak by dominant flood-generating process using the
# fixed antecedent-window thresholds of Davenport et al. (2020):
#   - Rainfall:     total_rain_mm > 10 AND total_melt_mm < 5
#   - Snowmelt:     total_melt_mm > 10 AND total_rain_mm < 5
#   - Rain-on-Snow: total_rain_mm > 5  AND total_melt_mm > 5
#   - Unclassified: does not meet any of the above criteria
#
# Antecedent window: 8 days ending on peak date (lookback_day 1–8),
# matching Davenport et al. Fig. 2b. Window totals are summed from
# rain_mm_day and melt_mm_day in accum_curves_60d.rds (Script 5).
#
# Snowmelt is qsm_flux (direct WLDAS snowmelt flux, mm/day) summed over
# the window, not net SWE decrease. This avoids masking of melt by
# concurrent mid-window accumulation.
#
# Completeness: peaks with fewer than 8 days of WLDAS data in the window
# are marked unclassified (complete_window = FALSE).
#
# Peak selection follows Yu et al. (2022):
#   - Per-type selection and demotion applied first; types with fewer than
#     MIN_TYPE_SAMPLE events are demoted to "unclassified" before
#     single-sample selection.
#   - Single sample: largest m independent events where m = max(25, n_water_years),
#                    ranked directly by peak flow with no per-water-year cap.
#   - Per type:      largest n independent events where n = max(20, n_water_years),
#                    ranked directly by peak flow with no per-water-year cap.
#   Selection flags stored in in_single_sample and in_type_sample columns.
#
# INPUTS:
#   - _data/accum_curves_60d.rds  : daily rain/melt per peak (Script 5)
#   - _data/pot_peaks_all.rds     : POT peaks with pot_threshold_cfs
#   - _data/site_metadata.rds     : regulated flag, basin attributes
#
# OUTPUT:
#   - _data/pot_peaks_classified_depth_thresh.rds : peaks with flood_type,
#       flood_type_raw, pot_threshold_cfs, water_year,
#       in_single_sample, in_type_sample, complete_window,
#       total_rain_mm, total_melt_mm, n_days_wldas
#       (flood_type_raw retains the threshold-based label before any
#       sample-size demotion; flood_type is the final classification)
# ============================================================================

library(tidyverse)

source("functions.R")

# ============================================================================
# CONFIGURATION
# ============================================================================
accum_file   <- "_data/accum_curves_60d.rds"
peaks_file   <- "_data/pot_peaks_all.rds"
meta_file    <- "_data/site_metadata.rds"
output_file  <- "_data/pot_peaks_classified_depth_thresh.rds"

# Davenport et al. (2020) thresholds (mm over 8-day window)
ANTECEDENT_DAYS  <-  8   # lookback_day 1 (peak day) through 8
THRESH_RAIN_HIGH <- 10
THRESH_RAIN_LOW  <-  5
THRESH_MELT_HIGH <- 10
THRESH_MELT_LOW  <-  5

# Yu et al. (2022) top-n sample size thresholds
MIN_SINGLE_SAMPLE <- 25
MIN_TYPE_SAMPLE   <- 20

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading data...")
accum <- readRDS(accum_file)
peaks <- readRDS(peaks_file) %>%
  mutate(peak_date = as.Date(peak_date))
meta  <- readRDS(meta_file)

if (!"pot_threshold_cfs" %in% names(peaks)) {
  stop("pot_threshold_cfs column not found in pot_peaks_all.rds. ",
       "Re-run Script 3 to regenerate.")
}

# Deduplicate to one row per peak (pot_peaks_all is one row per peak-window day)
peaks <- peaks %>%
  distinct(site_no, peak_id, peak_date, peak_flow, pot_threshold_cfs)

message(sprintf("Peaks: %d across %d sites",
                nrow(peaks), n_distinct(peaks$site_no)))
message(sprintf("Accum curve rows: %d", nrow(accum)))

# ============================================================================
# COMPUTE 8-DAY WINDOW TOTALS
# ============================================================================
message("\nComputing Davenport 8-day antecedent window totals...")

window_totals <- accum %>%
  filter(lookback_day <= ANTECEDENT_DAYS) %>%
  group_by(site_no, peak_id) %>%
  summarise(
    total_rain_mm = sum(rain_mm_day, na.rm = TRUE),
    total_melt_mm = sum(melt_mm_day, na.rm = TRUE),
    n_days_wldas  = sum(!is.na(rain_mm_day)),
    .groups = "drop"
  ) %>%
  mutate(complete_window = n_days_wldas == ANTECEDENT_DAYS)

n_incomplete <- sum(!window_totals$complete_window)
if (n_incomplete > 0) {
  message(sprintf("  Warning: %d peaks have incomplete 8-day windows and will be unclassified",
                  n_incomplete))
}

# ============================================================================
# CLASSIFY
# ============================================================================
message("Applying Davenport et al. classification thresholds...")

classified <- peaks %>%
  left_join(window_totals, by = c("site_no", "peak_id")) %>%
  left_join(meta %>% select(site_no, regulated), by = "site_no") %>%
  mutate(
    flood_type_raw = case_when(
      !complete_window | is.na(total_rain_mm)                                        ~ "unclassified",
      total_rain_mm > THRESH_RAIN_HIGH & total_melt_mm < THRESH_MELT_LOW            ~ "rainfall",
      total_melt_mm > THRESH_MELT_HIGH & total_rain_mm < THRESH_RAIN_LOW            ~ "snowmelt",
      total_rain_mm > THRESH_RAIN_LOW  & total_melt_mm > THRESH_MELT_LOW            ~ "rain_on_snow",
      TRUE                                                                           ~ "unclassified"
    ),
    flood_type = factor(flood_type_raw,
                        levels = c("rainfall", "snowmelt", "rain_on_snow", "unclassified"))
  )

# ============================================================================
# PRE-SELECTION SUMMARY
# ============================================================================
message("\n=== CLASSIFICATION SUMMARY (pre-selection) ===")
type_counts <- classified %>%
  dplyr::count(flood_type) %>%
  mutate(pct = n / sum(n) * 100)

walk2(type_counts$flood_type, type_counts$n, function(t, n) {
  pct <- type_counts$pct[type_counts$flood_type == t]
  message(sprintf("  %-14s : %5d  (%5.1f%%)", t, n, pct))
})

message(sprintf("\nTotal peaks: %d | Sites: %d",
                nrow(classified), n_distinct(classified$site_no)))

# ============================================================================
# YU ET AL. (2022) TOP-N PEAK SELECTION
# ============================================================================
message("\nApplying Yu et al. (2022) top-n peak selection...")

# water_year() is defined in functions.R

classified <- classified %>%
  mutate(water_year = water_year(peak_date))

# Single-sample selection precedes demotion, following Yu et al. (2022):
# the single-sample pool is drawn from all peaks irrespective of flood type,
# before any per-type MIN_TYPE_SAMPLE filtering. Only peaks unclassified by
# data insufficiency are excluded.
single_sample_ids <- classified %>%
  filter(complete_window) %>%
  group_by(site_no) %>%
  mutate(m    = pmax(MIN_SINGLE_SAMPLE, n_distinct(water_year)),
         rank = row_number(desc(peak_flow))) %>%
  filter(rank <= m) %>%
  ungroup() %>%
  select(site_no, peak_id) %>%
  mutate(in_single_sample = TRUE)

# Per-type selection: top-n per (site, flood_type, water_year), retained only
# for site-types reaching MIN_TYPE_SAMPLE.
type_sample_ids <- classified %>%
  filter(flood_type != "unclassified") %>%
  group_by(site_no, flood_type) %>%
  mutate(n    = pmax(MIN_TYPE_SAMPLE, n_distinct(water_year)),
         rank = row_number(desc(peak_flow))) %>%
  filter(rank <= n) %>%
  filter(n() >= MIN_TYPE_SAMPLE) %>%
  ungroup() %>%
  select(site_no, peak_id, flood_type) %>%
  mutate(in_type_sample = TRUE)

retained_site_types <- type_sample_ids %>% distinct(site_no, flood_type)

# Demote types below MIN_TYPE_SAMPLE to "unclassified". in_single_sample has
# already been assigned and is unaffected.
classified <- classified %>%
  mutate(flood_type = if_else(
    flood_type != "unclassified" &
      !paste(site_no, flood_type) %in%
        paste(retained_site_types$site_no, retained_site_types$flood_type),
    "unclassified",
    as.character(flood_type)
  )) %>%
  mutate(flood_type = factor(flood_type,
                             levels = c("rainfall", "snowmelt", "rain_on_snow", "unclassified")))

classified <- classified %>%
  left_join(type_sample_ids %>% select(site_no, peak_id, in_type_sample),
            by = c("site_no", "peak_id")) %>%
  left_join(single_sample_ids, by = c("site_no", "peak_id")) %>%
  mutate(
    in_single_sample = replace_na(in_single_sample, FALSE),
    in_type_sample   = replace_na(in_type_sample,   FALSE)
  )

# ============================================================================
# POST-SELECTION SUMMARY
# ============================================================================
message("\n=== CLASSIFICATION SUMMARY (post-selection) ===")
type_counts_post <- classified %>%
  dplyr::count(flood_type) %>%
  mutate(pct = n / sum(n) * 100)

walk2(type_counts_post$flood_type, type_counts_post$n, function(t, n) {
  pct <- type_counts_post$pct[type_counts_post$flood_type == t]
  message(sprintf("  %-14s : %5d  (%5.1f%%)", t, n, pct))
})

selection_summary <- classified %>%
  group_by(site_no) %>%
  summarise(
    n_single_sample = sum(in_single_sample),
    n_rainfall      = sum(flood_type == "rainfall"     & in_type_sample),
    n_snowmelt      = sum(flood_type == "snowmelt"     & in_type_sample),
    n_rain_on_snow  = sum(flood_type == "rain_on_snow" & in_type_sample),
    .groups = "drop"
  )

message(sprintf("\nSites with >= %d single-sample peaks : %d",
                MIN_SINGLE_SAMPLE,
                sum(selection_summary$n_single_sample >= MIN_SINGLE_SAMPLE)))
message(sprintf("Sites with >= %d rainfall peaks      : %d",
                MIN_TYPE_SAMPLE, sum(selection_summary$n_rainfall    >= MIN_TYPE_SAMPLE)))
message(sprintf("Sites with >= %d snowmelt peaks      : %d",
                MIN_TYPE_SAMPLE, sum(selection_summary$n_snowmelt    >= MIN_TYPE_SAMPLE)))
message(sprintf("Sites with >= %d rain_on_snow peaks  : %d",
                MIN_TYPE_SAMPLE, sum(selection_summary$n_rain_on_snow >= MIN_TYPE_SAMPLE)))

# ============================================================================
# SAVE
# ============================================================================
saveRDS(classified, output_file)
message(sprintf("\nSaved %d classified peaks to %s", nrow(classified), output_file))
message(sprintf("  MIN_TYPE_SAMPLE  : %d", MIN_TYPE_SAMPLE))
message(sprintf("  Columns: %s", paste(names(classified), collapse = ", ")))