#!/usr/bin/env Rscript

# ============================================================================
# Script 6b: Flood Type Classification — Fractional Threshold
# ============================================================================
# Classifies each data-sufficient POT peak into flood type using only
# rain_frac (the weighted rain fraction over the basin-adaptive antecedent
# window, search_radius) as the classification feature. SWE is not used.
#
# Classification rules (evaluated in order):
#   1. Unclassified  : n_days_data < completeness threshold, rain_frac NA, or
#                      total antecedent-window flux below NOISE_FLOOR_MM
#                      (flux_sufficient == FALSE)
#   2. Rainfall      : rain_frac >= RAIN_FRAC_THRESH (= 0.99)
#   3. Snowmelt      : rain_frac <  MIXED_SPLIT (population median of mixed events)
#   4. Rain-on-snow  : rain_frac >= MIXED_SPLIT
#
# "Snowmelt" and "rain_on_snow" are operationally defined: they describe
# position in the rain_frac distribution of mixed events (rain_frac < 0.99),
# not direct physical process attribution. The labels are retained for
# compatibility with the dashboard's type_labels / type_colors machinery.
#
# The RAIN_FRAC_THRESH tolerance (0.99 rather than 1.0) absorbs floating-point
# artifacts in the rain_frac ratio itself. See phase_0_signoff.md for
# rationale.
#
# MIXED_SPLIT is the population-level median of rain_frac over all
# data-sufficient peaks with rain_frac < RAIN_FRAC_THRESH, pooled across all
# sites. It is computed once from the full classified population and applied
# uniformly. Its value is saved as an attribute on the output for
# reproducibility.
#
# This script produces the 3-class output only. The 2-class output
# (frac_thresh_2class: rainfall vs. snow_influenced) is derived downstream
# with remap_2class() in functions.R.
#
# This script only LABELS peaks. It does no sample selection: the Yu et al.
# (2022) top-n selection and the MIN_TYPE_SAMPLE demotion were removed for
# milestone M1 (see _planning/agu2026_interim_plan.md, item U2). Magnitude
# filtering happens in Script 7, through each flood type's fitting threshold.
#
# fit_eligible = data_sufficient & !is.na(rain_frac) & flux_sufficient.
# It defines the single-sample pool in Script 7, so that pool excludes the
# same peaks the per-type pools exclude (the rule the old single-sample
# selection used).
#
# Thresholds live in CLASS_CONFIG and are written to the output as the
# "class_config" attribute. Defaults reproduce the initial results.
#
# INPUTS:
#   - _data/accum_curves_60d.rds  : daily rain/melt per peak (Script 5)
#   - _data/pot_peaks_all.rds     : candidate peaks with search_radius (Script 3)
#   - _data/site_metadata.rds     : regulated flag, basin attributes
#
# OUTPUT:
#   - _data/pot_peaks_classified_frac_thresh_3class.rds
#     flood_type is an exact factor copy of flood_type_raw, kept for one
#     release so existing consumers keep working. To be removed in M2.
#
# Author : Mike Talbot, PhD Student, Colorado State University
# Date   : 2026-05-21
# ============================================================================

library(tidyverse)

source("functions.R")

# ============================================================================
# CONFIGURATION
# ============================================================================
accum_file  <- "_data/accum_curves_60d.rds"
peaks_file  <- "_data/pot_peaks_all.rds"
meta_file   <- "_data/site_metadata.rds"
output_file <- "_data/pot_peaks_classified_frac_thresh_3class.rds"

CLASS_CONFIG <- list(
  # alpha_r. Peaks with rain_frac >= this value are classified as rainfall.
  # The 0.99 tolerance absorbs floating-point artifacts in the rain_frac
  # ratio itself.
  rain_frac_thresh = 0.99,

  # alpha_s. NULL = population median of rain_frac among mixed events
  # (the proposal's rule). A number fixes the split instead, for
  # sensitivity runs.
  mixed_split      = NULL,

  # Noise floor: a peak-level eligibility threshold on total antecedent-
  # window flux (total_flux = summed rain + melt over the window). Peaks
  # whose total flux falls below this value are excluded from analysis
  # entirely as suspected model/observation mismatches.
  noise_floor_mm   = 0.1,

  # Completeness: minimum fraction of search_radius days that must have
  # WLDAS data for a peak to be classified.
  min_days_frac    = 0.90
)

RAIN_FRAC_THRESH <- CLASS_CONFIG$rain_frac_thresh
NOISE_FLOOR_MM   <- CLASS_CONFIG$noise_floor_mm
MIN_DAYS_SHORT   <- CLASS_CONFIG$min_days_frac

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading data...")
accum <- readRDS(accum_file)
peaks <- readRDS(peaks_file) %>%
  mutate(peak_date = as.Date(peak_date))
meta  <- readRDS(meta_file)

# pot_peaks_all.rds is one row per (peak, window day) — deduplicate to one
# row per peak for use as the classification lookup table.
# ASSUMPTION: peak_id uniquely identifies each independent peak within a site.
peaks <- peaks %>%
  dplyr::distinct(site_no, peak_id, peak_date, peak_flow,
                  candidate_threshold_cfs, search_radius,
                  antecedent_runoff_mm, antecedent_days_used)

message(sprintf("Peaks: %d across %d sites",
                nrow(peaks), dplyr::n_distinct(peaks$site_no)))
message(sprintf("Accum curve rows: %d", nrow(accum)))

# ============================================================================
# COMPUTE RAIN_FRAC
# ============================================================================
message("\nComputing rain_frac (search_radius window)...")

# Weighted rain fraction over each peak's basin-adaptive antecedent window.
# total_rain_flux / total_melt_flux are plain sums of the raw daily values —
# no per-day activity floor. rain_frac = total rain / total flux; NA is
# returned for peaks where total flux == 0. Peaks whose total_flux falls
# below NOISE_FLOOR_MM are excluded downstream via the explicit
# flux_sufficient flag (see below), not filtered out here.
peaks_with_windows <- peaks %>%
  dplyr::mutate(short_window = search_radius)

feat <- accum %>%
  dplyr::inner_join(
    peaks_with_windows %>% dplyr::select(site_no, peak_id, short_window),
    by = c("site_no", "peak_id")
  ) %>%
  dplyr::filter(lookback_day <= short_window) %>%
  dplyr::group_by(site_no, peak_id) %>%
  dplyr::summarise(
    total_rain_flux = sum(rain_mm_day, na.rm = TRUE),
    total_melt_flux = sum(melt_mm_day, na.rm = TRUE),
    total_flux      = total_rain_flux + total_melt_flux,
    n_days_data     = dplyr::first(n_days_data),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    rain_frac = dplyr::if_else(
      total_flux > 0,
      total_rain_flux / total_flux,
      NA_real_
    )
  )

# ============================================================================
# JOIN AND APPLY DATA-SUFFICIENCY FLAG
# ============================================================================
message("Joining features and applying data-sufficiency flag...")

classified <- peaks %>%
  dplyr::left_join(
    peaks_with_windows %>% dplyr::select(site_no, peak_id, short_window),
    by = c("site_no", "peak_id")
  ) %>%
  dplyr::left_join(feat, by = c("site_no", "peak_id")) %>%
  dplyr::left_join(meta %>% dplyr::select(site_no, regulated), by = "site_no") %>%
  dplyr::mutate(
    data_sufficient = n_days_data >= ceiling(short_window * MIN_DAYS_SHORT),
    flux_sufficient = total_flux >= NOISE_FLOOR_MM
  )

# ============================================================================
# COMPUTE POPULATION MEDIAN SPLIT
# ============================================================================
message("\nComputing population median split for mixed events...")

# Mixed events: data-sufficient peaks with rain_frac < RAIN_FRAC_THRESH.
# MIXED_SPLIT is computed once from all sites pooled, then applied uniformly.
mixed_rain_frac <- classified %>%
  dplyr::filter(data_sufficient, !is.na(rain_frac), flux_sufficient,
                rain_frac < RAIN_FRAC_THRESH) %>%
  dplyr::pull(rain_frac)

MIXED_SPLIT <- if (is.null(CLASS_CONFIG$mixed_split)) {
  median(mixed_rain_frac)
} else {
  CLASS_CONFIG$mixed_split
}

message(sprintf("  Mixed events (n)  : %d", length(mixed_rain_frac)))
message(sprintf("  MIXED_SPLIT       : %.4f", MIXED_SPLIT))
message(sprintf("  rain_frac range   : [%.4f, %.4f]",
                min(mixed_rain_frac), max(mixed_rain_frac)))

# ============================================================================
# CLASSIFY
# ============================================================================
message("\nClassifying peaks...")

classified <- classified %>%
  dplyr::mutate(
    # Decision tree (evaluated in order):
    #   1. Unclassified : data insufficient, rain_frac NA, or total
    #                     antecedent-window flux below NOISE_FLOOR_MM
    #                     (flux_sufficient == FALSE)
    #   2. Rainfall     : rain_frac >= RAIN_FRAC_THRESH (= 0.99)
    #   3. Snowmelt     : rain_frac < MIXED_SPLIT (lower half of mixed events)
    #   4. Rain-on-snow : rain_frac >= MIXED_SPLIT (upper half of mixed events)
    #
    # "Snowmelt" and "rain_on_snow" are operationally defined — see header.
    flood_type_raw = dplyr::case_when(
      !data_sufficient | is.na(rain_frac) | !flux_sufficient ~ "unclassified",
      rain_frac >= RAIN_FRAC_THRESH       ~ "rainfall",
      rain_frac <  MIXED_SPLIT            ~ "snowmelt",
      TRUE                                ~ "rain_on_snow"
    ),
    flood_type = factor(flood_type_raw,
                        levels = c("rainfall", "snowmelt",
                                   "rain_on_snow", "unclassified")),
    # water_year() is defined in functions.R
    water_year   = water_year(peak_date),
    fit_eligible = dplyr::coalesce(
      data_sufficient & !is.na(rain_frac) & flux_sufficient, FALSE)
  )

# ============================================================================
# SUMMARY
# ============================================================================
message("\n=== CLASSIFICATION SUMMARY ===")
type_counts <- classified %>%
  dplyr::count(flood_type) %>%
  dplyr::mutate(pct = n / sum(n) * 100)

purrr::walk2(type_counts$flood_type, type_counts$n, function(t, n) {
  pct <- type_counts$pct[type_counts$flood_type == t]
  message(sprintf("  %-14s : %5d  (%5.1f%%)", t, n, pct))
})

message(sprintf("\nTotal peaks: %d | Sites: %d",
                nrow(classified), dplyr::n_distinct(classified$site_no)))

reg_summary <- classified %>%
  dplyr::count(regulated, flood_type) %>%
  dplyr::mutate(pct = n / sum(n) * 100)

message("\n=== BY REGULATION STATUS ===")
message(paste(capture.output(print(reg_summary)), collapse = "\n"))

# Per-site type counts, for reference only. The minimum-count rule is
# applied in Script 7 to exceedances of each type's fitting threshold.
site_type_counts <- classified %>%
  dplyr::filter(flood_type_raw != "unclassified") %>%
  dplyr::count(site_no, flood_type_raw)

for (t in c("rainfall", "snowmelt", "rain_on_snow")) {
  message(sprintf("Sites with >= 20 %-12s peaks (all magnitudes): %d", t,
                  sum(site_type_counts$flood_type_raw == t & site_type_counts$n >= 20)))
}

# ============================================================================
# SAVE
# ============================================================================

# Attach MIXED_SPLIT and the config as attributes so downstream scripts and
# diagnostics can retrieve the thresholds without re-reading the raw data.
attr(classified, "mixed_split")       <- MIXED_SPLIT
attr(classified, "rain_frac_thresh")  <- RAIN_FRAC_THRESH
attr(classified, "n_mixed_events")    <- length(mixed_rain_frac)
attr(classified, "class_config")      <- CLASS_CONFIG

saveRDS(classified, output_file)
message(sprintf("\nSaved %d classified peaks to %s", nrow(classified), output_file))
message(sprintf("  RAIN_FRAC_THRESH : %.2f", RAIN_FRAC_THRESH))
message(sprintf("  MIXED_SPLIT      : %.4f%s", MIXED_SPLIT,
                if (is.null(CLASS_CONFIG$mixed_split)) " (population median)" else " (fixed)"))
message(sprintf("  NOISE_FLOOR_MM   : %.2f", NOISE_FLOOR_MM))
message(sprintf("  Columns          : %s",   paste(names(classified), collapse = ", ")))
