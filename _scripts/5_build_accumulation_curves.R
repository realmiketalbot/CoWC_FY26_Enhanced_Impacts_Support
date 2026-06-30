#!/usr/bin/env Rscript

# ============================================================================
# Script 5: Build 60-Day Accumulation Curves
# ============================================================================
# Constructs per-peak antecedent accumulation curves from WLDAS daily Parquet
# files, covering the 60-day window preceding each POT peak. Output is
# consumed by the classification scripts (6a, 6b, 6c).
#
# For each peak x lookback_day, the following daily values are extracted:
#   rain_mm_day  : basin-mean rainfall flux (mm/day)
#   melt_mm_day  : basin-mean snowmelt flux (mm/day)
#   swe_mm       : basin-mean snow water equivalent (mm, state variable)
#   cumul_rain_mm: cumulative rainfall from peak day backward
#   cumul_melt_mm: cumulative snowmelt from peak day backward
#
# Note: SWE is a state variable (kg/m² ≡ mm), not a flux, so it is carried
# as a daily value only — not accumulated. Downstream scripts may summarise
# it (mean, min, max) over the antecedent window as needed.
#
# Also produces endmembers_60d.rds: data-driven end-member classification
# from the joint distribution of 60-day rainfall and snowmelt totals.
# End-members are used by diagnostic scripts and are available for use
# in classification validation.
#
# INPUTS:
#   - _data/pot_peaks_all.rds
#   - _data/WLDAS/processed/daily_by_date/YYYYMMDD_wldas.parquet
#   - _data/site_metadata.rds
#
# OUTPUTS:
#   - _data/accum_curves_60d.rds  : long-format accumulation curves,
#                                   one row per (peak_id x lookback_day)
#   - _data/endmembers_60d.rds    : end-member peak IDs with 60-day totals
# ============================================================================

library(tidyverse)
library(arrow)

# ============================================================================
# CONFIGURATION
# ============================================================================
peaks_file     <- "_data/pot_peaks_all.rds"
wldas_dir      <- "_data/WLDAS/processed/daily_by_date"
meta_file      <- "_data/site_metadata.rds"
accum_output   <- "_data/accum_curves_60d.rds"
endmem_output  <- "_data/endmembers_60d.rds"

LOOKBACK_DAYS <- 59   # days before peak; + peak day = 60 days total

# End-member identification thresholds (used for endmembers_60d.rds only;
# do not affect the accumulation curves themselves)
ENDMEMBER_PCTILE <- 0.25  # opposing signal must be in bottom X% to qualify
PRIMARY_PCTILE   <- 0.50  # primary signal must be in top X% to qualify

# ============================================================================
# LOAD PEAKS
# ============================================================================
message("Loading POT peaks...")
peaks <- readRDS(peaks_file) %>%
  mutate(peak_date = as.Date(peak_date),
         GAGE_ID   = site_no) %>%
  distinct(site_no, GAGE_ID, peak_id, peak_date, peak_flow)

message(sprintf("Loaded %d peaks from %d sites",
                nrow(peaks), n_distinct(peaks$site_no)))

# ============================================================================
# CHECK WLDAS COVERAGE
# ============================================================================
parquet_files <- list.files(wldas_dir, pattern = "\\.parquet$", full.names = TRUE)
if (length(parquet_files) == 0) stop(sprintf("No parquet files found in %s", wldas_dir))

wldas_dates <- parquet_files %>%
  basename() %>%
  str_extract("\\d{8}") %>%
  as.Date(format = "%Y%m%d") %>%
  sort()

message(sprintf("WLDAS coverage: %s to %s (%d days)",
                min(wldas_dates), max(wldas_dates), length(wldas_dates)))

# Filter to peaks with full 60-day WLDAS coverage
peaks <- peaks %>%
  mutate(
    window_start    = peak_date - LOOKBACK_DAYS,
    window_end      = peak_date,
    wldas_available = window_start >= min(wldas_dates) &
                      window_end   <= max(wldas_dates)
  ) %>%
  filter(wldas_available) %>%
  select(-wldas_available)

message(sprintf("Peaks with full %d-day WLDAS coverage: %d",
                LOOKBACK_DAYS + 1L, nrow(peaks)))

# ============================================================================
# LOAD WLDAS DAILY DATA
# ============================================================================
required_dates <- peaks %>%
  rowwise() %>%
  reframe(date = seq(window_start, window_end, by = "day")) %>%
  distinct(date) %>%
  filter(date %in% wldas_dates) %>%
  pull(date)

message(sprintf("Loading WLDAS data for %d unique dates...", length(required_dates)))

wldas_data <- map_dfr(required_dates, function(d) {
  f <- file.path(wldas_dir,
                 sprintf("%s_wldas.parquet", format(d, "%Y%m%d")))
  tryCatch(
    read_parquet(f, col_select = c("date", "GAGE_ID",
                                   "rainf_flux", "qsm_flux", "swe")) %>%
      mutate(date    = as.Date(date),
             GAGE_ID = as.character(GAGE_ID)),
    error = function(e) {
      warning(sprintf("Could not read %s: %s", basename(f), e$message))
      NULL
    }
  )
})

message(sprintf("Loaded %d site-day rows", nrow(wldas_data)))

# ============================================================================
# BUILD ACCUMULATION CURVES
# ============================================================================
message("Building accumulation curves...")

peak_windows <- peaks %>%
  rowwise() %>%
  reframe(
    site_no   = site_no,
    GAGE_ID   = GAGE_ID,
    peak_id   = as.character(peak_id),
    peak_date = peak_date,
    peak_flow = peak_flow,
    date      = seq(window_start, window_end, by = "day")
  ) %>%
  mutate(
    # lookback_day: 1 = peak day, 2 = one day before, etc.
    lookback_day = as.integer(peak_date - date) + 1L
  ) %>%
  left_join(wldas_data, by = c("GAGE_ID", "date"))

accum_long <- peak_windows %>%
  group_by(site_no, peak_id, peak_date, peak_flow) %>%
  arrange(lookback_day, .by_group = TRUE) %>%
  mutate(
    cumul_rain_mm = cumsum(replace_na(rainf_flux, 0)),
    cumul_melt_mm = cumsum(replace_na(qsm_flux,   0)),
    n_days_data   = sum(!is.na(rainf_flux))
  ) %>%
  ungroup() %>%
  select(
    site_no, peak_id, peak_date, peak_flow,
    lookback_day, date,
    rain_mm_day   = rainf_flux,
    melt_mm_day   = qsm_flux,
    swe_mm        = swe,
    cumul_rain_mm,
    cumul_melt_mm,
    n_days_data
  )

message(sprintf("Accumulation curves built: %d peak-day rows across %d peaks",
                nrow(accum_long), n_distinct(accum_long$peak_id)))

saveRDS(accum_long, accum_output)
message(sprintf("Saved: %s", accum_output))

# ============================================================================
# END-MEMBER IDENTIFICATION
# ============================================================================
message("\nIdentifying end-members from joint 60-day totals...")

# Require near-complete windows (allow up to 3 missing days)
totals_60d <- accum_long %>%
  filter(lookback_day == 60) %>%
  select(site_no, peak_id, peak_date, peak_flow,
         total_rain_60d = cumul_rain_mm,
         total_melt_60d = cumul_melt_mm,
         n_days_data) %>%
  filter(n_days_data >= 57)

message(sprintf("Peaks with >= 57 days of data: %d", nrow(totals_60d)))

rain_low_cut  <- quantile(totals_60d$total_rain_60d, ENDMEMBER_PCTILE,       na.rm = TRUE)
rain_high_cut <- quantile(totals_60d$total_rain_60d, 1 - PRIMARY_PCTILE,     na.rm = TRUE)
melt_low_cut  <- quantile(totals_60d$total_melt_60d, ENDMEMBER_PCTILE,       na.rm = TRUE)
melt_high_cut <- quantile(totals_60d$total_melt_60d, 1 - PRIMARY_PCTILE,     na.rm = TRUE)

message(sprintf("Rain cutoffs : low = %.1f mm, high = %.1f mm",
                rain_low_cut, rain_high_cut))
message(sprintf("Melt cutoffs : low = %.1f mm, high = %.1f mm",
                melt_low_cut, melt_high_cut))

totals_60d <- totals_60d %>%
  mutate(
    endmember_type = case_when(
      total_rain_60d >= rain_high_cut &
        total_melt_60d <= melt_low_cut  ~ "rain",
      total_melt_60d >= melt_high_cut &
        total_rain_60d <= rain_low_cut  ~ "melt",
      TRUE                              ~ "neither"
    )
  )

message("End-member counts:")
print(count(totals_60d, endmember_type))

saveRDS(totals_60d, endmem_output)
message(sprintf("Saved: %s", endmem_output))

message("\nScript 5 complete.")
message(sprintf("  Columns in accum_curves_60d: %s",
                paste(names(accum_long), collapse = ", ")))