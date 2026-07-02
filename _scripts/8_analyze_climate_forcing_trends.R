#!/usr/bin/env Rscript

# ============================================================================
# Script 8: Annual Climate Forcing Trend Analysis
# ============================================================================
#
# PURPOSE:
#   Derives five POT-based annual site-level time series and one annual
#   maximum (AM) rain fraction series from WLDAS basin-mean outputs and
#   USGS POT peaks, then tests each for monotonic trends using the
#   Mann-Kendall / Theil-Sen framework. The AM rain fraction (am_rain_frac)
#   is the rain_frac of the single highest-flow classified peak per site-year;
#   am_flood_type records that peak's flood type label. Diagnostic and
#   publication-quality figures only — no Shiny pre-computation outputs.
#
# METHODOLOGICAL NOTES:
#
#   Why peak SWE (max), not mean SWE?
#     SWE is a storage state. The ecologically and hydrologically relevant
#     quantity is the peak snowpack that accumulated over the winter, which
#     drives spring melt potential. Mean SWE confounds accumulation and melt
#     timing within the same year. Peak SWE is the standard metric in
#     snowpack trend literature (Mote et al. 2018, Fyfe et al. 2017).
#
#   Why median rain fraction (not mean)?
#     rain_frac is bounded [0, 1] and within-year distributions are often
#     skewed. The median is more robust to a single anomalous peak (e.g.,
#     one large rainfall event during a predominantly snowmelt year) and
#     avoids leverage from boundary values.
#
#   Why all POT peaks for event rate (including "unclassified")?
#     The event rate measures whether the number of high-flow events per year
#     is changing, independent of process type. "Unclassified" peaks are those
#     above the PoT threshold where WLDAS data were insufficient to determine
#     the generating mechanism — they are real flood events. Excluding them
#     would artificially deflate the rate at sites with shorter WLDAS records.
#
# INPUTS:
#   _data/WLDAS/processed/daily_by_date/        (Script 4)
#   _data/pot_peaks_classified_frac_thresh_3class.rds (Script 6b)
#   _data/site_metadata.rds                     (Script 0)
#   _data/station_spatial.rds                   (Script 1, fallback for coords)
#   functions.R
#
# OUTPUTS:
#   _data/annual_climate_forcing_summary.rds    one row per (site_no, water_year);
#                                               includes am_rain_frac, am_flood_type
#   _data/climate_forcing_trends.rds            one row per (site_no, variable);
#                                               includes am_rain_frac trend rows
#   _diagnostics/climate_forcing_trends/        figure files
#
# Author : Mike Talbot, PhD Student, Colorado State University
# Date   : 2026-05-12
# ============================================================================

## ── Required-package check ────────────────────────────────────────────────────
# 'maps' is needed for ggplot2::map_data("state") / maps::map("state", ...).
# It is not in the original spec list but is required for the Colorado boundary.
required_pkgs <- c("tidyverse", "arrow", "Kendall", "trend", "sf",
                   "patchwork", "scales", "ggrepel", "maps")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      quietly = FALSE, FUN.VALUE = logical(1L))]
if (length(missing_pkgs) > 0L) {
  stop(
    "\nThe following packages are required but not installed:\n  ",
    paste(missing_pkgs, collapse = ", "),
    "\n\nInstall with:\n  install.packages(c(",
    paste0('"', missing_pkgs, '"', collapse = ", "), "))\n"
  )
}

## ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(patchwork)
library(scales)
library(ggrepel)

source("functions.R")

## ── Configuration ─────────────────────────────────────────────────────────────

WLDAS_DIR    <- "_data/WLDAS/processed/daily_by_date"
PEAKS_FILE   <- "_data/pot_peaks_classified_frac_thresh_3class.rds"
META_FILE    <- "_data/site_metadata.rds"
SPATIAL_FILE <- "_data/station_spatial.rds"
DIAG_DIR     <- "_diagnostics/climate_forcing_trends"
OUT_SUMMARY  <- "_data/annual_climate_forcing_summary.rds"
OUT_TRENDS   <- "_data/climate_forcing_trends.rds"

MIN_YEARS_ELIGIBLE <- 20L   # non-NA years required per variable for trend testing
MIN_SITES_REGIONAL <- 10L   # minimum sites needed to include a WY in regional agg
MK_ALPHA           <- 0.05
MIN_WY_DAYS_FRAC   <- 0.90  # minimum fraction of expected WLDAS days for a valid annual sum

VAR_KEYS <- c("total_precip_mm", "peak_swe_mm", "max_peak_discharge_cfs",
              "pot_event_rate", "median_rain_frac", "am_rain_frac")

VAR_LABELS <- c(
  total_precip_mm        = "Total Precip. (mm)",
  peak_swe_mm            = "Peak SWE (mm)",
  max_peak_discharge_cfs = "Max Peak Q (cfs)",
  pot_event_rate         = "POT Event Count",
  median_rain_frac       = "Median Rain Fraction",
  am_rain_frac           = "AM Rain Fraction"
)

VAR_SLOPE_UNITS <- c(
  total_precip_mm        = "mm / yr",
  peak_swe_mm            = "mm / yr",
  max_peak_discharge_cfs = "cfs / yr",
  pot_event_rate         = "events / yr",
  median_rain_frac       = "fraction / yr",
  am_rain_frac           = "fraction / yr"
)

# Neutral sequential palette — one color per variable
VAR_COLORS <- c(
  total_precip_mm        = "#4e79a7",
  peak_swe_mm            = "#76b7b2",
  max_peak_discharge_cfs = "#59a14f",
  pot_event_rate         = "#edc948",
  median_rain_frac       = "#b07aa1",
  am_rain_frac           = "#f28e2b"
)

# Trend direction colors — matches convention in 8_analyze_flood_type_temporal_trends.R
SLOPE_COLORS <- c(positive = "#b2182b", negative = "#2166ac", zero = "#888888")

## ── Output directory ──────────────────────────────────────────────────────────
if (!dir.exists(DIAG_DIR)) dir.create(DIAG_DIR, recursive = TRUE)

## ── Water year helper ─────────────────────────────────────────────────────────
# WY Y = Oct 1 of (Y-1) through Sep 30 of Y.
water_year <- function(d) {
  yr <- as.integer(format(d, "%Y"))
  mo <- as.integer(format(d, "%m"))
  ifelse(mo >= 10L, yr + 1L, yr)
}

## ── Safe max (avoids -Inf for all-NA groups) ──────────────────────────────────
safe_max <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

## ── Load POT peaks ────────────────────────────────────────────────────────────
message("Loading POT peaks...")

peaks_raw <- readRDS(PEAKS_FILE)

# One row per (peak, window day) — deduplicate to one row per independent peak.
# ASSUMPTION: 'peak_id' uniquely identifies each independent peak within a site.
peaks <- dplyr::distinct(peaks_raw, site_no, peak_id, .keep_all = TRUE)

n_dropped <- nrow(peaks_raw) - nrow(peaks)
if (n_dropped > 0L)
  message(sprintf("  Deduplication: dropped %d duplicate rows", n_dropped))

# ASSUMPTION: 'water_year', 'flood_type_raw', and 'rain_frac' are present
# in pot_peaks_classified_frac_thresh_3class.rds (added by Script 6b).
required_peak_cols <- c("site_no", "peak_id", "peak_date", "peak_flow",
                        "water_year", "flood_type_raw", "rain_frac")
missing_cols <- setdiff(required_peak_cols, names(peaks))
if (length(missing_cols) > 0L)
  stop("Missing required columns in peaks file: ",
       paste(missing_cols, collapse = ", "))

all_sites <- sort(unique(peaks$site_no))

message(sprintf("  Peaks: %d rows | %d sites | WY %d–%d",
                nrow(peaks),
                length(all_sites),
                min(peaks$water_year, na.rm = TRUE),
                max(peaks$water_year, na.rm = TRUE)))

## ── Site metadata and spatial coordinates ─────────────────────────────────────
message("Loading metadata and spatial coordinates...")

meta <- readRDS(META_FILE)

if (!all(c("site_no", "regulated") %in% names(meta)))
  stop("site_metadata.rds must contain 'site_no' and 'regulated' columns.")

reg_lookup <- dplyr::select(meta, site_no, regulated)

# Drop sf geometry if present — sfc_POINT is not vctrs-compatible in dplyr
spatial_raw <- readRDS(SPATIAL_FILE)
spatial <- if (inherits(spatial_raw, "sf")) {
  sf::st_drop_geometry(spatial_raw)
} else {
  as.data.frame(spatial_raw)
}

# ASSUMPTION: station_spatial.rds contains 'site_no', 'lat', 'lon' columns.
if (!all(c("site_no", "lat", "lon") %in% names(spatial)))
  stop("station_spatial.rds must contain 'site_no', 'lat', and 'lon' columns.")

spatial_pts <- dplyr::select(spatial, site_no, lon, lat)

## ── Open WLDAS dataset and verify schema ──────────────────────────────────────
message("Opening WLDAS parquet dataset...")

# ASSUMPTION: all parquet files in WLDAS_DIR share a uniform schema matching
# the output of Script 4: GAGE_ID (chr), date (Date), rainf_flux, snowf_flux,
# swe, qsm_flux (numeric, post-conversion to mm/day or mm).
wldas_ds <- arrow::open_dataset(WLDAS_DIR)

wldas_cols <- names(wldas_ds)
required_wldas_cols <- c("GAGE_ID", "date", "rainf_flux", "snowf_flux", "swe")
missing_wldas <- setdiff(required_wldas_cols, wldas_cols)
if (length(missing_wldas) > 0L)
  stop(
    "Required WLDAS columns absent from parquet schema: ",
    paste(missing_wldas, collapse = ", "),
    "\nAvailable columns: ", paste(wldas_cols, collapse = ", ")
  )

message(sprintf("  Schema OK (%d columns). Querying for %d sites...",
                length(wldas_cols), length(all_sites)))

## ── Collect WLDAS: filter to sites of interest, select needed columns ─────────
# arrow pushes the GAGE_ID filter and column selection to the parquet reader,
# so only the relevant columns and rows are loaded into memory.
wldas_raw <- wldas_ds %>%
  dplyr::filter(GAGE_ID %in% all_sites) %>%
  dplyr::select(GAGE_ID, date, rainf_flux, snowf_flux, swe) %>%
  dplyr::collect()

message(sprintf("  Collected: %d site-day rows | %d sites | %s to %s",
                nrow(wldas_raw),
                dplyr::n_distinct(wldas_raw$GAGE_ID),
                as.character(min(wldas_raw$date, na.rm = TRUE)),
                as.character(max(wldas_raw$date, na.rm = TRUE))))

## ── GAGE_ID / site_no alignment check ────────────────────────────────────────
# Check AFTER collecting — verifies the actual values present in the files.
wldas_sites_found <- sort(unique(wldas_raw$GAGE_ID))
in_peaks_not_wldas <- setdiff(all_sites, wldas_sites_found)
in_wldas_not_peaks <- setdiff(wldas_sites_found, all_sites)

if (length(in_peaks_not_wldas) > 0L || length(in_wldas_not_peaks) > 0L) {
  stop(sprintf(
    paste0(
      "\nGAGE_ID / site_no mismatch detected.\n",
      "  In peaks but not WLDAS (%d sites): %s%s\n",
      "  In WLDAS but not peaks (%d sites): %s%s\n",
      "Both should be zero-padded 8-digit USGS IDs (e.g. '09085000').\n",
      "If format differs, standardize in Script 4 before re-running."
    ),
    length(in_peaks_not_wldas),
    paste(head(in_peaks_not_wldas, 5), collapse = ", "),
    if (length(in_peaks_not_wldas) > 5) ", ..." else "",
    length(in_wldas_not_peaks),
    paste(head(in_wldas_not_peaks, 5), collapse = ", "),
    if (length(in_wldas_not_peaks) > 5) ", ..." else ""
  ))
}

message(sprintf("  GAGE_ID / site_no check passed: %d sites in both datasets.",
                length(all_sites)))

## ── WLDAS record start check ──────────────────────────────────────────────────
wldas_start_wy <- wldas_raw %>%
  dplyr::mutate(wy = water_year(date)) %>%
  dplyr::group_by(GAGE_ID) %>%
  dplyr::summarise(first_wy = min(wy, na.rm = TRUE), .groups = "drop")

n_late <- sum(wldas_start_wy$first_wy > 1981L)
pct_late <- 100 * n_late / nrow(wldas_start_wy)

if (pct_late > 10) {
  message(sprintf(
    "\n  *** FLAG: WLDAS record starts after WY1981 for %.0f%% of sites (%d / %d).",
    pct_late, n_late, nrow(wldas_start_wy)
  ))
  message("  ***       This may indicate missing parquet files for early years.")
  message("  ***       Check WLDAS_DIR content and rerun Script 4 if needed.\n")
} else {
  message(sprintf("  Record start check OK (%.0f%% of sites start after WY1981).",
                  pct_late))
}

## ── WLDAS annual day coverage ─────────────────────────────────────────────────
# Count WLDAS days per (site, water_year). Years with fewer than MIN_WY_DAYS_FRAC
# of expected days have their WLDAS-derived sums masked to NA so that partial
# years (e.g. a final year where the parquet ends mid-water-year) are excluded
# from Mann-Kendall testing rather than entering as artificially low values.
message("Computing WLDAS day coverage per site-year...")

wldas_day_coverage <- wldas_raw %>%
  dplyr::mutate(water_year = water_year(date)) %>%
  dplyr::group_by(GAGE_ID, water_year) %>%
  dplyr::summarise(n_days = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(
    expected_days = 365L + dplyr::if_else(
      water_year %% 4L == 0L & (water_year %% 100L != 0L | water_year %% 400L == 0L),
      1L, 0L
    ),
    frac_covered  = n_days / expected_days,
    is_incomplete = frac_covered < MIN_WY_DAYS_FRAC
  ) %>%
  dplyr::rename(site_no = GAGE_ID)

n_incomplete_site_years <- sum(wldas_day_coverage$is_incomplete, na.rm = TRUE)
if (n_incomplete_site_years > 0L) {
  incomplete_summary <- wldas_day_coverage %>%
    dplyr::filter(is_incomplete) %>%
    dplyr::group_by(water_year) %>%
    dplyr::summarise(
      n_sites   = dplyr::n(),
      mean_frac = mean(frac_covered, na.rm = TRUE),
      .groups   = "drop"
    )
  message(sprintf(
    "\n  *** FLAG: %d site-years have WLDAS coverage < %.0f%% of expected days.",
    n_incomplete_site_years, 100 * MIN_WY_DAYS_FRAC
  ))
  for (i in seq_len(nrow(incomplete_summary))) {
    r <- incomplete_summary[i, ]
    message(sprintf("  ***       WY%d: %d sites | mean coverage %.1f%%",
                    r$water_year, r$n_sites, 100 * r$mean_frac))
  }
  message("  ***       total_precip_mm and peak_swe_mm will be NA for these site-years.\n")
} else {
  message(sprintf(
    "  Coverage OK: all site-years have >= %.0f%% of expected WLDAS days.",
    100 * MIN_WY_DAYS_FRAC
  ))
}

## ── Variable 1 & 2: Annual WLDAS summaries ────────────────────────────────────
message("\nDeriving annual WLDAS variables (total precip, peak SWE)...")

wldas_annual <- wldas_raw %>%
  dplyr::mutate(water_year = water_year(date)) %>%
  dplyr::group_by(GAGE_ID, water_year) %>%
  dplyr::summarise(
    # Precip: sum daily rainfall + snowfall (both in mm/day). Days with NA
    # contribute 0 to the sum (na.rm = TRUE); appropriate for a reanalysis
    # product where gaps are rare and not structurally informative.
    total_precip_mm = sum(rainf_flux + snowf_flux, na.rm = TRUE),
    # SWE: state variable — take the peak (max) within the WY, not sum.
    # safe_max returns NA instead of -Inf when all observations are missing.
    peak_swe_mm = safe_max(swe),
    .groups = "drop"
  ) %>%
  dplyr::rename(site_no = GAGE_ID)

# Mask incomplete years: a partial year's sum() is a fraction of the true
# annual total, not a valid climate observation.
wldas_annual <- wldas_annual %>%
  dplyr::left_join(
    dplyr::select(wldas_day_coverage, site_no, water_year, is_incomplete),
    by = c("site_no", "water_year")
  ) %>%
  dplyr::mutate(
    total_precip_mm = dplyr::if_else(is_incomplete, NA_real_, total_precip_mm),
    peak_swe_mm     = dplyr::if_else(is_incomplete, NA_real_, peak_swe_mm)
  ) %>%
  dplyr::select(-is_incomplete)

message(sprintf("  WLDAS annual rows: %d", nrow(wldas_annual)))

## ── Variables 3 & 4: Annual peak max discharge and event count ────────────────
message("Deriving annual peak-based variables (max Q, event rate)...")

# Include ALL peaks (even flood_type_raw == "unclassified") for event rate.
# max() is safe here because every group has at least one row after grouping.
peaks_annual_qr <- peaks %>%
  dplyr::group_by(site_no, water_year) %>%
  dplyr::summarise(
    max_peak_discharge_cfs = max(peak_flow, na.rm = TRUE),
    pot_event_rate         = dplyr::n(),
    .groups = "drop"
  )

## ── Variable 5: Median rain fraction ─────────────────────────────────────────
message("Deriving annual median rain fraction...")

# Restrict to classified (non-"unclassified") peaks with valid rain_frac.
# Years with no classifiable peaks yield NA (not 0) — no imputation.
peaks_annual_rf <- peaks %>%
  dplyr::filter(
    as.character(flood_type_raw) != "unclassified",
    !is.na(rain_frac)
  ) %>%
  dplyr::group_by(site_no, water_year) %>%
  dplyr::summarise(
    median_rain_frac = stats::median(rain_frac, na.rm = TRUE),
    .groups = "drop"
  )

## ── AM peak per site-year (annual maximum of classified peaks) ────────────────
message("Deriving AM peak classification per site-year...")

# The AM peak is the classified peak with the highest peak_flow in each
# site-year. Unclassified peaks are excluded because they lack a valid rain_frac.
peaks_am <- peaks %>%
  dplyr::filter(as.character(flood_type_raw) != "unclassified") %>%
  dplyr::group_by(site_no, water_year) %>%
  dplyr::slice_max(peak_flow, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    site_no,
    water_year,
    am_rain_frac  = rain_frac,
    am_flood_type = as.character(flood_type_raw)
  )

## ── Build unified annual summary ──────────────────────────────────────────────
message("\nJoining annual variables into summary tibble...")

# Full join: retain all site-years present in any source so the coverage
# diagnostic shows the true data availability (NAs preserved, not dropped).
annual_summary <- dplyr::full_join(
    wldas_annual, peaks_annual_qr, by = c("site_no", "water_year")
  ) %>%
  dplyr::full_join(peaks_annual_rf, by = c("site_no", "water_year")) %>%
  dplyr::left_join(peaks_am,        by = c("site_no", "water_year")) %>%
  dplyr::left_join(reg_lookup, by = "site_no") %>%
  dplyr::arrange(site_no, water_year) %>%
  dplyr::select(
    site_no, water_year, regulated,
    total_precip_mm, peak_swe_mm,
    max_peak_discharge_cfs, pot_event_rate, median_rain_frac,
    am_rain_frac, am_flood_type
  )

# Mask WLDAS variables to the gauge-active period: water years with no POT
# peaks are outside the USGS record and should not contribute to trend tests
# or coverage diagnostics, even though WLDAS coverage is continuous.
annual_summary <- annual_summary %>%
  dplyr::mutate(
    total_precip_mm = ifelse(is.na(pot_event_rate), NA_real_, total_precip_mm),
    peak_swe_mm     = ifelse(is.na(pot_event_rate), NA_real_, peak_swe_mm),
    am_rain_frac    = ifelse(is.na(pot_event_rate), NA_real_, am_rain_frac)
  )

message(sprintf("  Summary: %d rows | %d sites | WY %d–%d",
                nrow(annual_summary),
                dplyr::n_distinct(annual_summary$site_no),
                min(annual_summary$water_year, na.rm = TRUE),
                max(annual_summary$water_year, na.rm = TRUE)))

saveRDS(annual_summary, OUT_SUMMARY)
message(sprintf("  Saved: %s", OUT_SUMMARY))

## ── Site eligibility per variable ─────────────────────────────────────────────
message(sprintf(
  "\nApplying eligibility filter: ≥20 non-NA water years per variable..."
))

eligibility <- purrr::map_dfr(VAR_KEYS, function(vk) {
  annual_summary %>%
    dplyr::group_by(site_no) %>%
    dplyr::summarise(
      n_valid  = sum(!is.na(.data[[vk]])),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(variable = vk,
                  eligible = n_valid >= MIN_YEARS_ELIGIBLE)
})

elig_summary <- eligibility %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(
    n_eligible = sum(eligible),
    n_total    = dplyr::n(),
    pct_elig   = round(100 * n_eligible / n_total, 1),
    .groups    = "drop"
  )

message("\n  Eligibility summary (record ≥20 yr):")
for (i in seq_len(nrow(elig_summary))) {
  r <- elig_summary[i, ]
  message(sprintf("    %-28s : %d / %d sites (%4.0f%%)",
                  VAR_LABELS[r$variable], r$n_eligible, r$n_total, r$pct_elig))
  if (r$pct_elig < 80) {
    message(sprintf(
      "\n  *** FLAG: eligibility filter drops >20%% of sites for %s (%.0f%% eligible)\n",
      VAR_LABELS[r$variable], r$pct_elig
    ))
  }
}

eligible_combos <- eligibility %>%
  dplyr::filter(eligible) %>%
  dplyr::select(site_no, variable)

message(sprintf("\n  Total eligible (site \u00d7 variable): %d combinations",
                nrow(eligible_combos)))

## ── Per-site trend tests ──────────────────────────────────────────────────────
message("\nRunning Mann-Kendall + Theil-Sen tests...")

run_mk_sens <- function(site, variable, series_df) {

  site_ts <- series_df %>%
    dplyr::filter(site_no == site) %>%
    dplyr::arrange(water_year)

  vals     <- site_ts[[variable]]
  vals_ok  <- vals[!is.na(vals)]
  n_series <- length(vals_ok)

  if (n_series < MIN_YEARS_ELIGIBLE) return(NULL)

  mk_res <- tryCatch(Kendall::MannKendall(vals_ok), error = function(e) NULL)
  ss_res <- tryCatch(trend::sens.slope(vals_ok),    error = function(e) NULL)

  # Store slopes at full precision; signif() is for display only.
  tibble::tibble(
    site_no       = site,
    variable      = variable,
    n_years       = n_series,
    series_mean   = mean(vals_ok),
    series_median = stats::median(vals_ok),
    series_sd     = stats::sd(vals_ok),
    tau           = if (!is.null(mk_res)) unname(mk_res$tau[1])    else NA_real_,
    p_value       = if (!is.null(mk_res)) unname(mk_res$sl[1])     else NA_real_,
    sen_slope     = if (!is.null(ss_res)) unname(ss_res$estimates) else NA_real_
  )
}

trend_list <- purrr::map2(
  eligible_combos$site_no,
  eligible_combos$variable,
  run_mk_sens,
  series_df = annual_summary
)

trends_raw <- dplyr::bind_rows(trend_list) %>%
  dplyr::left_join(reg_lookup, by = "site_no") %>%
  dplyr::mutate(
    significant  = !is.na(p_value) & p_value < MK_ALPHA,
    sig_positive = significant & !is.na(sen_slope) & sen_slope > 0,
    sig_negative = significant & !is.na(sen_slope) & sen_slope < 0
  ) %>%
  dplyr::arrange(variable, site_no)

message(sprintf("  Tests complete: %d rows", nrow(trends_raw)))

# Save clean version (drop helper flag columns — derivable from p_value / sen_slope)
saveRDS(
  dplyr::select(trends_raw, -significant, -sig_positive, -sig_negative),
  OUT_TRENDS
)
message(sprintf("  Saved: %s", OUT_TRENDS))

## ── Aggregate binomial summary ────────────────────────────────────────────────
message(sprintf("\n%s\n=== AGGREGATE TREND SUMMARY (alpha = %.2f) ===\n%s",
                strrep("=", 60), MK_ALPHA, strrep("=", 60)))

run_binom_agg <- function(df, pop_label) {
  df %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      population      = pop_label,
      n_eligible      = dplyr::n(),
      n_sig           = sum(significant,  na.rm = TRUE),
      n_sig_pos       = sum(sig_positive, na.rm = TRUE),
      n_sig_neg       = sum(sig_negative, na.rm = TRUE),
      binom_p_any_sig = tryCatch(
        binom.test(sum(significant,  na.rm = TRUE), dplyr::n(), p = MK_ALPHA)$p.value,
        error = function(e) NA_real_
      ),
      binom_p_pos_sig = tryCatch(
        binom.test(sum(sig_positive, na.rm = TRUE), dplyr::n(), p = MK_ALPHA)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      pct_sig     = round(100 * n_sig     / n_eligible, 1),
      pct_sig_pos = round(100 * n_sig_pos / n_eligible, 1),
      pct_sig_neg = round(100 * n_sig_neg / n_eligible, 1)
    )
}

aggregate <- run_binom_agg(trends_raw, "all_eligible") %>%
  dplyr::select(
    variable, population, n_eligible,
    n_sig, pct_sig, n_sig_pos, pct_sig_pos, n_sig_neg, pct_sig_neg,
    binom_p_any_sig, binom_p_pos_sig
  ) %>%
  dplyr::arrange(variable)

for (i in seq_len(nrow(aggregate))) {
  r        <- aggregate[i, ]
  sf_flag  <- if (!is.na(r$binom_p_any_sig) && r$binom_p_any_sig < MK_ALPHA) " *" else ""
  pos_flag <- if (!is.na(r$binom_p_pos_sig) && r$binom_p_pos_sig < MK_ALPHA) " *" else ""
  message(sprintf(
    "    %-26s | n=%3d | sig: %3d (%4.0f%%) [p=%.3f%s] | pos: %3d (%4.0f%%) [p=%.3f%s] | neg: %3d (%4.0f%%)",
    VAR_LABELS[r$variable], r$n_eligible,
    r$n_sig,     r$pct_sig,     r$binom_p_any_sig, sf_flag,
    r$n_sig_pos, r$pct_sig_pos, r$binom_p_pos_sig, pos_flag,
    r$n_sig_neg, r$pct_sig_neg
  ))
}
message(strrep("=", 60))
message("  * = binomial p < 0.05 vs. 5% null\n")

## ─────────────────────────────────────────────────────────────────────────────
## FIGURES
## ─────────────────────────────────────────────────────────────────────────────
message("Generating diagnostic figures...")

# Shared: Colorado state boundary
co_map  <- maps::map("state", regions = "colorado", plot = FALSE, fill = TRUE)
co_poly <- data.frame(
  long  = co_map$x,
  lat   = co_map$y,
  group = cumsum(is.na(co_map$x) | is.na(co_map$y))
)

# Shared: site order south-to-north
lat_order <- spatial_pts %>%
  dplyr::arrange(lat) %>%
  dplyr::pull(site_no)

## ── Figure 1: Record coverage heatmap ────────────────────────────────────────
message("  Figure 1: record coverage heatmap...")

coverage_long <- annual_summary %>%
  tidyr::pivot_longer(
    cols      = dplyr::all_of(VAR_KEYS),
    names_to  = "variable",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    has_data  = !is.na(value),
    var_label = factor(VAR_LABELS[variable],
                       levels = unname(VAR_LABELS[VAR_KEYS])),
    site_f    = factor(site_no, levels = lat_order)
  )

cov_panels <- purrr::map(seq_along(VAR_KEYS), function(j) {
  vk       <- VAR_KEYS[j]
  df       <- dplyr::filter(coverage_long, variable == vk)
  is_bot   <- j == length(VAR_KEYS)
  ggplot2::ggplot(df, ggplot2::aes(x = water_year, y = site_f, fill = has_data)) +
    ggplot2::geom_tile(na.rm = TRUE) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = VAR_COLORS[vk], "FALSE" = "grey93"),
      guide  = "none"
    ) +
    ggplot2::scale_x_continuous(breaks = seq(1980, 2025, 5), expand = c(0, 0)) +
    ggplot2::labs(
      title = VAR_LABELS[vk],
      x     = if (is_bot) "Water year" else NULL,
      y     = NULL
    ) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid   = ggplot2::element_blank(),
      plot.title   = ggplot2::element_text(size = 8, face = "bold"),
      axis.text.x  = if (is_bot) ggplot2::element_text()  else ggplot2::element_blank(),
      axis.ticks.x = if (is_bot) ggplot2::element_line()  else ggplot2::element_blank()
    )
})

p_coverage <- patchwork::wrap_plots(cov_panels, ncol = 1) +
  patchwork::plot_annotation(
    title    = "Data coverage: annual values by site and water year",
    subtitle = "Colored = data present; grey = no data. Sites ordered south (bottom) to north (top).",
    theme    = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold")
    )
  )

ggplot2::ggsave(
  file.path(DIAG_DIR, "record_coverage_heatmap.png"),
  p_coverage, width = 14, height = 12, dpi = 150
)
message("    Saved: record_coverage_heatmap.png")

## ── Figure 2: Regional aggregate time series ──────────────────────────────────
message("  Figure 2: regional aggregate time series...")

regional_agg <- purrr::map_dfr(VAR_KEYS, function(vk) {
  annual_summary %>%
    dplyr::group_by(water_year) %>%
    dplyr::summarise(
      n_sites = sum(!is.na(.data[[vk]])),
      med_val = stats::median(.data[[vk]], na.rm = TRUE),
      q25_val = stats::quantile(.data[[vk]], 0.25, na.rm = TRUE),
      q75_val = stats::quantile(.data[[vk]], 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_sites >= MIN_SITES_REGIONAL) %>%
    dplyr::mutate(variable = vk)
})

ts_panels <- purrr::map(seq_along(VAR_KEYS), function(j) {
  vk     <- VAR_KEYS[j]
  df     <- dplyr::filter(regional_agg, variable == vk)
  col    <- VAR_COLORS[vk]
  is_bot <- j == length(VAR_KEYS)
  ggplot2::ggplot(df, ggplot2::aes(x = water_year)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q25_val, ymax = q75_val),
                         fill = col, alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(y = med_val), colour = col, linewidth = 1.1) +
    ggplot2::scale_x_continuous(breaks = seq(1980, 2025, 5)) +
    ggplot2::labs(
      title = VAR_LABELS[vk],
      x     = if (is_bot) "Water year" else NULL,
      y     = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 9, face = "bold"),
      axis.text.x  = if (is_bot) ggplot2::element_text()  else ggplot2::element_blank(),
      axis.ticks.x = if (is_bot) ggplot2::element_line()  else ggplot2::element_blank()
    )
})

p_ts <- patchwork::wrap_plots(ts_panels, ncol = 1) +
  patchwork::plot_annotation(
    title    = "Regional aggregate time series (all Colorado GAGES-II sites)",
    subtitle = sprintf(
      "Thick line: cross-site median; ribbon: 25th–75th percentile. Only WY with ≥10 sites shown."
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold")
    )
  )

ggplot2::ggsave(
  file.path(DIAG_DIR, "regional_aggregate_timeseries.png"),
  p_ts, width = 14, height = 12, dpi = 150
)
message("    Saved: regional_aggregate_timeseries.png")

## ── Figure 3: Five-panel trend map ────────────────────────────────────────────
message("  Figure 3: trend map (five panels)...")

trends_sp <- trends_raw %>%
  dplyr::left_join(spatial_pts, by = "site_no") %>%
  dplyr::filter(!is.na(lon), !is.na(lat), !is.na(sen_slope))

# Legend encoding — defined outside make_trend_map so patchwork's
# guides = "collect" sees identical scale specs across all five panels
# and merges them into one shared legend.
LEGEND_COLOURS <- c(
  pos_sig = "#b2182b",   # matches SLOPE_COLORS["positive"]
  neg_sig = "#2166ac",   # matches SLOPE_COLORS["negative"]
  not_sig = "grey60"
)
LEGEND_SHAPES <- c(pos_sig = 15L, neg_sig = 15L, not_sig = 1L)
LEGEND_LABELS <- c(
  pos_sig = "Positive trend (p < 0.05)",
  neg_sig = "Negative trend (p < 0.05)",
  not_sig = "Not significant"
)
LEGEND_SIZES  <- c(pos_sig = 2.5, neg_sig = 2.5, not_sig = 1.5)

make_trend_map <- function(vk) {
  df <- dplyr::filter(trends_sp, variable == vk) %>%
    dplyr::mutate(
      slope_dir  = dplyr::case_when(
        sen_slope >= 0 ~ "positive",
        sen_slope <  0 ~ "negative"
      ),
      legend_cat = dplyr::case_when(
        significant & slope_dir == "positive" ~ "pos_sig",
        significant & slope_dir == "negative" ~ "neg_sig",
        TRUE                                  ~ "not_sig"
      ),
      legend_cat = factor(legend_cat, levels = c("pos_sig", "neg_sig", "not_sig"))
    )

  ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = co_poly,
      ggplot2::aes(x = long, y = lat, group = group),
      fill = "grey93", colour = "grey55", linewidth = 0.3
    ) +
    # Not-significant sites drawn first (beneath) as open grey circles
    ggplot2::geom_point(
      data = dplyr::filter(df, legend_cat == "not_sig"),
      ggplot2::aes(x = lon, y = lat,
                   colour = legend_cat, shape = legend_cat, size = legend_cat),
      alpha = 0.6
    ) +
    # Significant sites drawn on top as filled colored squares
    ggplot2::geom_point(
      data = dplyr::filter(df, legend_cat != "not_sig"),
      ggplot2::aes(x = lon, y = lat,
                   colour = legend_cat, shape = legend_cat, size = legend_cat),
      alpha = 0.9
    ) +
    ggplot2::scale_colour_manual(
      values = LEGEND_COLOURS, labels = LEGEND_LABELS,
      name = NULL, drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = LEGEND_SHAPES, labels = LEGEND_LABELS,
      name = NULL, drop = FALSE
    ) +
    ggplot2::scale_size_manual(
      values = LEGEND_SIZES, labels = LEGEND_LABELS,
      name = NULL, guide = "none"
    ) +
    ggplot2::coord_fixed(
      ratio = 1.3,
      xlim  = range(co_poly$long, na.rm = TRUE),
      ylim  = range(co_poly$lat,  na.rm = TRUE)
    ) +
    ggplot2::labs(
      title    = VAR_LABELS[vk],
      subtitle = sprintf("n = %d sites", nrow(df)),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text        = ggplot2::element_blank(),
      axis.ticks       = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey88", linewidth = 0.2)
    )
}

map_panels <- purrr::map(VAR_KEYS, make_trend_map)
map_grid <- patchwork::wrap_plots(map_panels, ncol = 3, nrow = 2,
                                  guides = "collect") +
  patchwork::plot_annotation(
    title    = "Site-level Mann-Kendall trend results — Annual climate forcing variables",
    subtitle = "■ = significant trend (p < 0.05), colored by direction; ○ = not significant (grey)",
    theme    = ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "grey40")
    )
  ) &
  ggplot2::theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.key.size  = ggplot2::unit(0.5, "cm"),
    legend.text      = ggplot2::element_text(size = 9)
  )

ggplot2::ggsave(
  file.path(DIAG_DIR, "trend_map_five_panel.png"),
  map_grid, width = 18, height = 12, dpi = 150
)
message("    Saved: trend_map_five_panel.png")

## ── Figure 4: Slope distribution violins — all eligible sites ─────────────────
message("  Figure 4: slope distribution violins...")

# Winsorize sen_slope to 2nd-98th percentile within each variable for display
# only — stored trends_raw output is unchanged.
violin_df <- trends_raw %>%
  dplyr::filter(!is.na(sen_slope)) %>%
  dplyr::group_by(variable) %>%
  dplyr::mutate(
    lo                = stats::quantile(sen_slope, 0.02, na.rm = TRUE),
    hi                = stats::quantile(sen_slope, 0.98, na.rm = TRUE),
    sen_slope_display = pmin(pmax(sen_slope, lo), hi)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    var_label = factor(VAR_LABELS[variable],
                       levels = unname(VAR_LABELS[VAR_KEYS])),
    variable  = factor(variable, levels = VAR_KEYS)   # enforce VAR_KEYS order
  )

violin_annot <- violin_df %>%
  dplyr::group_by(variable, var_label) %>%
  dplyr::summarise(
    med_slope = stats::median(sen_slope,         na.rm = TRUE),
    pct_sig   = round(100 * mean(significant,    na.rm = TRUE), 0),
    y_top     = max(sen_slope_display,           na.rm = TRUE),  # clipped range
    .groups   = "drop"
  ) %>%
  dplyr::mutate(
    label_txt = sprintf("med = %s\n%d%% sig", signif(med_slope, 4), pct_sig)
  )

# Build strip labels with units underneath
facet_labels <- setNames(
  paste0(unname(VAR_LABELS[VAR_KEYS]), "\n", unname(VAR_SLOPE_UNITS[VAR_KEYS])),
  VAR_KEYS
)

p_violin <- ggplot2::ggplot(
  violin_df,
  ggplot2::aes(x = "", y = sen_slope_display)   # winsorized y
) +
  ggplot2::geom_violin(fill = "#4393c3", alpha = 0.6, linewidth = 0.3,
                       trim = FALSE) +
  ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.4,
                        fill = "white", colour = "grey30") +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                      colour = "grey30", linewidth = 0.5) +
  ggplot2::facet_wrap(
    ~variable,
    scales   = "free_y",
    nrow     = 1,
    labeller = ggplot2::labeller(variable = facet_labels)
  ) +
  ggplot2::geom_text(
    data = violin_annot,
    ggplot2::aes(x = "", y = y_top, label = label_txt),
    size = 2.6, colour = "grey20", vjust = -0.15, lineheight = 1.1,
    inherit.aes = FALSE
  ) +
  ggplot2::labs(
    title    = "Distribution of Sen slopes — all eligible sites",
    subtitle = paste0(
      "Inner box: IQR; whiskers: 1.5 \u00d7 IQR; annotation: median slope ",
      "and % significant (p < 0.05). Distributions winsorized to 2nd\u201398th ",
      "percentile within variable for display."
    ),
    x = NULL,
    y = "Sen slope"
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    strip.text   = ggplot2::element_text(size = 8, lineheight = 1.2)
  )

ggplot2::ggsave(
  file.path(DIAG_DIR, "slope_distributions_all_sites.png"),
  p_violin, width = 12, height = 7, dpi = 150
)
message("    Saved: slope_distributions_all_sites.png")

## ── Figure 5: Aggregate summary bar chart ─────────────────────────────────────
message("  Figure 5: aggregate summary bar chart...")

# Reshape to long: one row per (variable × direction)
agg_bar_long <- aggregate %>%
  dplyr::mutate(
    var_label = factor(VAR_LABELS[variable],
                       levels = rev(unname(VAR_LABELS[VAR_KEYS]))),
    null_exp  = round(0.05 * n_eligible)
  ) %>%
  tidyr::pivot_longer(
    cols      = c(n_sig_pos, n_sig_neg),
    names_to  = "direction",
    values_to = "n_sig_dir"
  ) %>%
  dplyr::mutate(
    pct_dir   = dplyr::if_else(direction == "n_sig_pos", pct_sig_pos, pct_sig_neg),
    dir_label = factor(
      dplyr::if_else(direction == "n_sig_pos", "Positive slope", "Negative slope"),
      levels = c("Positive slope", "Negative slope")
    ),
    pct_lbl   = sprintf("%d%%", round(pct_dir))
  )

# Reference data: one row per variable (for null tick and binom p label)
null_ref <- aggregate %>%
  dplyr::mutate(
    var_label = factor(VAR_LABELS[variable],
                       levels = rev(unname(VAR_LABELS[VAR_KEYS]))),
    null_exp  = round(0.05 * n_eligible),
    binom_lbl = sprintf("binom p = %.3f", binom_p_any_sig)
  )

x_max <- max(agg_bar_long$n_sig_dir, na.rm = TRUE)

p_agg_bar <- ggplot2::ggplot(
  agg_bar_long,
  ggplot2::aes(x = var_label, y = n_sig_dir, fill = dir_label)
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.7), width = 0.6
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = pct_lbl),
    position = ggplot2::position_dodge(width = 0.7),
    hjust = -0.1, size = 3.0, colour = "grey20"
  ) +
  # Per-variable null expectation tick (dashed crossbar at 5% of n_eligible)
  ggplot2::geom_errorbar(
    data = null_ref,
    ggplot2::aes(x = var_label, ymin = null_exp, ymax = null_exp),
    width = 0.75, colour = "grey50", linetype = "dashed", linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  # "5% null" label anchored near the bottom-most bar group (first factor level)
  ggplot2::annotate(
    "text",
    x     = 0.62,
    y     = null_ref$null_exp[which.min(as.integer(null_ref$var_label))],
    label = "5% null",
    hjust = 0, size = 2.7, colour = "grey50"
  ) +
  # Binomial p-value at the right margin of each row
  ggplot2::geom_text(
    data = null_ref,
    ggplot2::aes(x = var_label, y = x_max * 1.35, label = binom_lbl),
    hjust = 0.5, size = 2.6, colour = "grey30",
    inherit.aes = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = c("Positive slope" = "#b2182b", "Negative slope" = "#2166ac"),
    name   = NULL
  ) +
  ggplot2::coord_flip(clip = "off") +
  ggplot2::expand_limits(y = x_max * 1.65) +
  ggplot2::labs(
    title    = "Aggregate Mann-Kendall evidence by variable (binomial test vs. 5% null)",
    subtitle = "Bars: # sites with significant positive (red) and negative (blue) Sen slopes; dashed line: 5% null expectation",
    x        = NULL,
    y        = "Number of significant sites"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    legend.position  = "bottom",
    panel.grid.minor = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(DIAG_DIR, "aggregate_trend_summary_table.png"),
  p_agg_bar, width = 10, height = 6, dpi = 150
)
message("    Saved: aggregate_trend_summary_table.png")

## ── Done ─────────────────────────────────────────────────────────────────────
message(sprintf("\n%s\n=== Script 8 (climate forcing trends) complete ===\n%s",
                strrep("=", 60), strrep("=", 60)))
message("  Data outputs:")
message(sprintf("    %s", OUT_SUMMARY))
message(sprintf("    %s", OUT_TRENDS))
message(sprintf("  Figures (%s/):", DIAG_DIR))
for (fn in c("record_coverage_heatmap.png",
             "regional_aggregate_timeseries.png",
             "trend_map_five_panel.png",
             "slope_distributions_all_sites.png",
             "aggregate_trend_summary_table.png")) {
  message(sprintf("    %s", fn))
}
