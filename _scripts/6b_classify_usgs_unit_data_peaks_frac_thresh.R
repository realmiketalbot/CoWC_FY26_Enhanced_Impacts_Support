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
# (frac_thresh_2class: rainfall vs. snow_influenced) is derived in Script 7 as
# a downstream remap.
#
# INPUTS:
#   - _data/accum_curves_60d.rds  : daily rain/melt per peak (Script 5)
#   - _data/pot_peaks_all.rds     : POT peaks with search_radius per peak
#   - _data/site_metadata.rds     : regulated flag, basin attributes
#
# OUTPUT:
#   - _data/pot_peaks_classified_frac_thresh_3class.rds
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

# Rainfall threshold: peaks with rain_frac >= this value are classified as
# rainfall. The 0.99 tolerance absorbs floating-point artifacts in the
# rain_frac ratio itself.
RAIN_FRAC_THRESH <- 0.99

# Noise floor: a peak-level eligibility threshold on total antecedent-window
# flux (total_flux = summed rain + melt over the window). Peaks whose total
# flux falls below this value are excluded from analysis entirely (both
# type-conditioned classification and the single-sample pool) as suspected
# model/observation mismatches, not merely zeroed out day-by-day.
NOISE_FLOOR_MM <- 0.1    # mm, applied to total_flux

# Completeness: minimum fraction of search_radius days that must have WLDAS
# data for a peak to be classified (rather than marked unclassified).
MIN_DAYS_SHORT <- 0.90

# Yu et al. (2022) top-n sample size thresholds
MIN_SINGLE_SAMPLE <- 25   # minimum peaks for single-sample GPD/LP3 fit
MIN_TYPE_SAMPLE   <- 20   # minimum peaks per type for per-type GPD/LP3 fit

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
                  pot_threshold_cfs, search_radius,
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

MIXED_SPLIT <- median(mixed_rain_frac)

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
                                   "rain_on_snow", "unclassified"))
  )

# ============================================================================
# PRE-SELECTION SUMMARY
# ============================================================================
message("\n=== CLASSIFICATION SUMMARY (pre-selection) ===")
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

message("\n=== BY REGULATION STATUS (pre-selection) ===")
message(paste(capture.output(print(reg_summary)), collapse = "\n"))

# ============================================================================
# YU ET AL. (2022) TOP-N PEAK SELECTION
# ============================================================================
message("\nApplying Yu et al. (2022) top-n peak selection...")

# water_year() is defined in functions.R

classified <- classified %>%
  dplyr::mutate(water_year = water_year(peak_date))

# Single-sample: largest m independent peaks per site, ranked directly by
# peak flow with no per-water-year cap.
# m = max(MIN_SINGLE_SAMPLE, number of unique water years at site).
# Shares the flux_sufficient eligibility criterion with the type-conditioned
# classification, so both the single-sample and per-type pools exclude the
# same suspected model/observation-mismatch peaks.
single_sample_ids <- classified %>%
  dplyr::filter(data_sufficient & !is.na(rain_frac) & flux_sufficient) %>%
  dplyr::group_by(site_no) %>%
  dplyr::mutate(m    = pmax(MIN_SINGLE_SAMPLE, dplyr::n_distinct(water_year)),
                rank = dplyr::row_number(dplyr::desc(peak_flow))) %>%
  dplyr::filter(rank <= m) %>%
  dplyr::ungroup() %>%
  dplyr::select(site_no, peak_id) %>%
  dplyr::mutate(in_single_sample = TRUE)

# Per-type: largest n independent peaks per (site, type), ranked directly by
# peak flow with no per-water-year cap.
# Types with fewer than MIN_TYPE_SAMPLE events at a site are dropped entirely.
type_sample_ids <- classified %>%
  dplyr::filter(flood_type != "unclassified") %>%
  dplyr::group_by(site_no, flood_type) %>%
  dplyr::mutate(n    = pmax(MIN_TYPE_SAMPLE, dplyr::n_distinct(water_year)),
                rank = dplyr::row_number(dplyr::desc(peak_flow))) %>%
  dplyr::filter(rank <= n) %>%
  dplyr::filter(dplyr::n() >= MIN_TYPE_SAMPLE) %>%
  dplyr::ungroup() %>%
  dplyr::select(site_no, peak_id, flood_type) %>%
  dplyr::mutate(in_type_sample = TRUE)

classified <- classified %>%
  dplyr::left_join(single_sample_ids, by = c("site_no", "peak_id")) %>%
  dplyr::left_join(
    type_sample_ids %>% dplyr::select(site_no, peak_id, in_type_sample),
    by = c("site_no", "peak_id")
  ) %>%
  dplyr::mutate(
    in_single_sample = tidyr::replace_na(in_single_sample, FALSE),
    in_type_sample   = tidyr::replace_na(in_type_sample,   FALSE)
  )

# Demote types that did not meet MIN_TYPE_SAMPLE at a given site to
# "unclassified". in_single_sample is unaffected by this demotion.
retained_site_types <- type_sample_ids %>% dplyr::distinct(site_no, flood_type)

classified <- classified %>%
  dplyr::mutate(flood_type = dplyr::if_else(
    flood_type != "unclassified" &
      !paste(site_no, flood_type) %in%
        paste(retained_site_types$site_no, retained_site_types$flood_type),
    "unclassified",
    as.character(flood_type)
  )) %>%
  dplyr::mutate(flood_type = factor(flood_type,
                                    levels = c("rainfall", "snowmelt",
                                               "rain_on_snow", "unclassified")))

# ============================================================================
# POST-SELECTION SUMMARY
# ============================================================================
message("\n=== CLASSIFICATION SUMMARY (post-selection) ===")
type_counts_post <- classified %>%
  dplyr::count(flood_type) %>%
  dplyr::mutate(pct = n / sum(n) * 100)

purrr::walk2(type_counts_post$flood_type, type_counts_post$n, function(t, n) {
  pct <- type_counts_post$pct[type_counts_post$flood_type == t]
  message(sprintf("  %-14s : %5d  (%5.1f%%)", t, n, pct))
})

selection_summary <- classified %>%
  dplyr::group_by(site_no) %>%
  dplyr::summarise(
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
                MIN_TYPE_SAMPLE, sum(selection_summary$n_rainfall     >= MIN_TYPE_SAMPLE)))
message(sprintf("Sites with >= %d snowmelt peaks      : %d",
                MIN_TYPE_SAMPLE, sum(selection_summary$n_snowmelt     >= MIN_TYPE_SAMPLE)))
message(sprintf("Sites with >= %d rain-on-snow peaks  : %d",
                MIN_TYPE_SAMPLE, sum(selection_summary$n_rain_on_snow >= MIN_TYPE_SAMPLE)))

# ============================================================================
# SAVE
# ============================================================================

# Attach MIXED_SPLIT as an attribute so downstream scripts and diagnostics
# can retrieve the threshold value without re-reading the raw data.
attr(classified, "mixed_split")       <- MIXED_SPLIT
attr(classified, "rain_frac_thresh")  <- RAIN_FRAC_THRESH
attr(classified, "n_mixed_events")    <- length(mixed_rain_frac)

saveRDS(classified, output_file)
message(sprintf("\nSaved %d classified peaks to %s", nrow(classified), output_file))
message(sprintf("  RAIN_FRAC_THRESH : %.2f", RAIN_FRAC_THRESH))
message(sprintf("  MIXED_SPLIT      : %.4f", MIXED_SPLIT))
message(sprintf("  MIN_TYPE_SAMPLE  : %d",   MIN_TYPE_SAMPLE))
message(sprintf("  Columns          : %s",   paste(names(classified), collapse = ", ")))
