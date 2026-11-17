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
# This script only LABELS peaks. It does no sample selection: the Yu et al.
# (2022) top-n selection and the MIN_TYPE_SAMPLE demotion were removed for
# milestone M1 (see _planning/agu2026_interim_plan.md, item U2). Magnitude
# filtering happens in Script 7, through each flood type's fitting threshold.
#
# fit_eligible marks peaks that could have been classified (complete WLDAS
# window). It defines the single-sample pool in Script 7, so that pool
# excludes the same data-insufficient peaks the per-type pools exclude.
#
# Thresholds live in CLASS_CONFIG and are written to the output as the
# "class_config" attribute, so sensitivity runs can be told apart.
#
# INPUTS:
#   - _data/accum_curves_60d.rds  : daily rain/melt per peak (Script 5)
#   - _data/pot_peaks_all.rds     : candidate peaks (Script 3)
#   - _data/site_metadata.rds     : regulated flag, basin attributes
#
# OUTPUT:
#   - _data/pot_peaks_classified_depth_thresh.rds : one row per peak with
#       site_no, peak_id, peak_date, peak_flow, candidate_threshold_cfs,
#       total_rain_mm, total_melt_mm, n_days_wldas, complete_window,
#       regulated, flood_type_raw, flood_type, water_year, fit_eligible
#       flood_type is an exact factor copy of flood_type_raw, kept for one
#       release so existing consumers keep working. To be removed in M2.
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

# Davenport et al. (2020) thresholds (mm over 8-day window). Defaults
# reproduce the initial results; change them here for sensitivity runs.
CLASS_CONFIG <- list(
  antecedent_days  =  8,   # lookback_day 1 (peak day) through 8
  thresh_rain_high = 10,
  thresh_rain_low  =  5,
  thresh_melt_high = 10,
  thresh_melt_low  =  5
)

ANTECEDENT_DAYS  <- CLASS_CONFIG$antecedent_days
THRESH_RAIN_HIGH <- CLASS_CONFIG$thresh_rain_high
THRESH_RAIN_LOW  <- CLASS_CONFIG$thresh_rain_low
THRESH_MELT_HIGH <- CLASS_CONFIG$thresh_melt_high
THRESH_MELT_LOW  <- CLASS_CONFIG$thresh_melt_low

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading data...")
accum <- readRDS(accum_file)
peaks <- readRDS(peaks_file) %>%
  mutate(peak_date = as.Date(peak_date))
meta  <- readRDS(meta_file)

if (!"candidate_threshold_cfs" %in% names(peaks)) {
  stop("candidate_threshold_cfs column not found in pot_peaks_all.rds. ",
       "Re-run Script 3 to regenerate.")
}

# Deduplicate to one row per peak (pot_peaks_all is one row per peak-window day)
peaks <- peaks %>%
  distinct(site_no, peak_id, peak_date, peak_flow, candidate_threshold_cfs)

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

# water_year() and fit_eligible
# water_year() is defined in functions.R
classified <- classified %>%
  mutate(
    water_year   = water_year(peak_date),
    fit_eligible = coalesce(complete_window, FALSE)
  )

# ============================================================================
# SUMMARY
# ============================================================================
message("\n=== CLASSIFICATION SUMMARY ===")
type_counts <- classified %>%
  dplyr::count(flood_type) %>%
  mutate(pct = n / sum(n) * 100)

walk2(type_counts$flood_type, type_counts$n, function(t, n) {
  pct <- type_counts$pct[type_counts$flood_type == t]
  message(sprintf("  %-14s : %5d  (%5.1f%%)", t, n, pct))
})

message(sprintf("\nTotal peaks: %d | Sites: %d",
                nrow(classified), n_distinct(classified$site_no)))

# Per-site type counts, for reference only. The minimum-count rule is
# applied in Script 7 to exceedances of each type's fitting threshold.
site_type_counts <- classified %>%
  filter(flood_type_raw != "unclassified") %>%
  count(site_no, flood_type_raw)

for (t in c("rainfall", "snowmelt", "rain_on_snow")) {
  message(sprintf("Sites with >= 20 %-12s peaks (all magnitudes): %d", t,
                  sum(site_type_counts$flood_type_raw == t & site_type_counts$n >= 20)))
}

# ============================================================================
# SAVE
# ============================================================================
attr(classified, "class_config") <- CLASS_CONFIG

saveRDS(classified, output_file)
message(sprintf("\nSaved %d classified peaks to %s", nrow(classified), output_file))
message(sprintf("  Config  : %s",
                paste(names(CLASS_CONFIG), unlist(CLASS_CONFIG), sep = "=", collapse = ", ")))
message(sprintf("  Columns : %s", paste(names(classified), collapse = ", ")))
