#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Days-to-Accumulation for Snow-Free Rainfall Events
# ============================================================================
# For POT peaks where SWE < 1 mm across a 14-day lookback window (i.e., no
# meaningful snow cover), computes how many days working backward from and
# including the peak date are needed to accumulate >= 5 mm and >= 10 mm of
# basin-mean rainfall.
#
# Purpose: inform the SRF/LRF duration threshold for the Sikorska et al.
# (2015) fuzzy decision tree classification (Script 6b). The distribution of
# days-to-accumulation across qualifying events provides an empirical basis
# for choosing a duration cutoff.
#
# INPUTS:
#   - pot_peaks_all.rds          : POT peaks (output of Script 3)
#   - YYYYMMDD_wldas.parquet     : daily WLDAS basin means (output of Script 4a)
#
# OUTPUTS:
#   - rainfall_duration_diag.rds : one row per qualifying peak with columns:
#       site_no, peak_id, peak_date, peak_flow,
#       days_to_5mm, days_to_10mm
#     NA in either column means the threshold was not reached within the
#     14-day lookback window.
# ============================================================================

library(tidyverse)
library(arrow)

# ============================================================================
# CONFIGURATION
# ============================================================================
pot_peaks_file <- "_data/pot_peaks_all.rds"
wldas_dir      <- "_data/WLDAS/processed/daily_by_date"
output_file    <- "_data/rainfall_duration_diag.rds"

lookback_days  <- 13   # days before peak; combined with peak day = 14 days total
swe_tolerance  <- 1.0  # mm; peaks with any day SWE >= this are excluded
rain_target_lo <- 5.0  # mm accumulation target (lower)
rain_target_hi <- 10.0 # mm accumulation target (upper)

# ============================================================================
# LOAD PEAKS
# ============================================================================
message("Loading POT peaks...")
peaks <- readRDS(pot_peaks_file) %>%
  mutate(peak_date = as.Date(peak_date)) %>%
  distinct(site_no, peak_id, peak_date, peak_flow)

message(sprintf("Loaded %d peaks across %d sites",
                nrow(peaks), n_distinct(peaks$site_no)))

# GAGE_ID == site_no for Colorado GAGES-II sites
peaks <- peaks %>%
  mutate(GAGE_ID = site_no)

# ============================================================================
# CHECK WLDAS COVERAGE
# ============================================================================
parquet_files <- list.files(wldas_dir, pattern = "\\.parquet$", full.names = TRUE)

if (length(parquet_files) == 0) {
  stop(sprintf("No parquet files found in %s", wldas_dir))
}

wldas_dates <- parquet_files %>%
  basename() %>%
  str_extract("\\d{8}") %>%
  as.Date(format = "%Y%m%d") %>%
  sort()

message(sprintf("WLDAS coverage: %s to %s (%d days)",
                min(wldas_dates), max(wldas_dates), length(wldas_dates)))

# Filter to peaks whose full 14-day window is within WLDAS coverage
peaks <- peaks %>%
  mutate(
    window_start    = peak_date - lookback_days,
    window_end      = peak_date,
    wldas_available = window_start >= min(wldas_dates) &
                      window_end   <= max(wldas_dates)
  )

n_outside <- sum(!peaks$wldas_available)
if (n_outside > 0) {
  message(sprintf("Excluding %d peaks outside WLDAS coverage", n_outside))
}

peaks <- peaks %>% filter(wldas_available)
message(sprintf("%d peaks remain after coverage filter", nrow(peaks)))

# ============================================================================
# LOAD REQUIRED WLDAS DATES
# ============================================================================
required_dates <- peaks %>%
  rowwise() %>%
  reframe(date = seq(window_start, window_end, by = "day")) %>%
  distinct(date) %>%
  filter(date %in% wldas_dates) %>%
  pull(date)

message(sprintf("Loading WLDAS data for %d unique dates...", length(required_dates)))

wldas_data <- map_dfr(required_dates, function(d) {
  f <- file.path(wldas_dir, sprintf("%s_wldas.parquet", format(d, "%Y%m%d")))
  tryCatch(
    read_parquet(f, col_select = c("date", "GAGE_ID", "rainf_flux", "swe")) %>%
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
# BUILD PER-PEAK 14-DAY WINDOWS
# ============================================================================
message("Assembling 14-day windows per peak...")

peak_windows <- peaks %>%
  rowwise() %>%
  reframe(
    site_no   = site_no,
    GAGE_ID   = GAGE_ID,
    peak_id   = peak_id,
    peak_date = peak_date,
    peak_flow = peak_flow,
    date      = seq(window_start, window_end, by = "day")
  ) %>%
  left_join(wldas_data, by = c("GAGE_ID", "date"))

# ============================================================================
# APPLY SWE FILTER
# ============================================================================
# Exclude any peak where SWE >= swe_tolerance on ANY day in the window.
# These are events with meaningful snow present and are not candidates for
# the pure rainfall duration diagnostic.
message(sprintf("Applying SWE filter (tolerance = %.1f mm)...", swe_tolerance))

snow_free_peaks <- peak_windows %>%
  group_by(site_no, peak_id) %>%
  summarise(
    max_swe_in_window = max(swe, na.rm = TRUE),
    n_days_data       = sum(!is.na(swe)),
    .groups           = "drop"
  ) %>%
  filter(n_days_data == (lookback_days + 1),   # require complete window
         max_swe_in_window < swe_tolerance)

message(sprintf("%d of %d peaks pass the snow-free filter",
                nrow(snow_free_peaks), n_distinct(peak_windows$peak_id)))

qualifying_windows <- peak_windows %>%
  semi_join(snow_free_peaks, by = c("site_no", "peak_id"))

# ============================================================================
# COMPUTE DAYS-TO-ACCUMULATION
# ============================================================================
# For each peak, work backward from peak_date (day 0) through day -13.
# days_to_Xmm = the number of days (inclusive of peak day) needed to reach
# the target accumulation. NA if the target is not reached within 14 days.
message("Computing days-to-accumulation...")

days_to_accum <- function(rain_vec, dates, peak_d, target_mm) {
  # rain_vec and dates are for the full window; sort descending by date
  # so index 1 = peak day, index 2 = day before, etc.
  ord      <- order(dates, decreasing = TRUE)
  rain_ord <- rain_vec[ord]

  cumrain <- cumsum(replace_na(rain_ord, 0))
  hit     <- which(cumrain >= target_mm)

  if (length(hit) == 0) NA_integer_ else hit[1L]
}

results <- qualifying_windows %>%
  group_by(site_no, peak_id, peak_date, peak_flow) %>%
  summarise(
    days_to_5mm  = days_to_accum(rainf_flux, date, peak_date[1], rain_target_lo),
    days_to_10mm = days_to_accum(rainf_flux, date, peak_date[1], rain_target_hi),
    total_rain_14d = sum(rainf_flux, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================================
# SUMMARY
# ============================================================================
message("\n--- Days-to-5mm summary ---")
print(summary(results$days_to_5mm))
message(sprintf("  Not reached within 14 days: %d (%.1f%%)",
                sum(is.na(results$days_to_5mm)),
                100 * mean(is.na(results$days_to_5mm))))

message("\n--- Days-to-10mm summary ---")
print(summary(results$days_to_10mm))
message(sprintf("  Not reached within 14 days: %d (%.1f%%)",
                sum(is.na(results$days_to_10mm)),
                100 * mean(is.na(results$days_to_10mm))))

message("\n--- Total 14-day rainfall summary (qualifying events) ---")
print(summary(results$total_rain_14d))

# ============================================================================
# SAVE
# ============================================================================
saveRDS(results, output_file)
message(sprintf("\nSaved %d qualifying peak rows to %s", nrow(results), output_file))
message("Columns: site_no, peak_id, peak_date, peak_flow, days_to_5mm, days_to_10mm, total_rain_14d")

site_meta <- readRDS("_data/site_metadata.rds")

elev_breaks <- quantile(
  site_meta$elev_mean_m, 
  probs = seq(0, 1, 0.25), 
  na.rm = TRUE
)

area_breaks <- quantile(
  site_meta$drain_sqkm,
  probs = seq(0, 1, 0.25), 
  na.rm = TRUE
)

plot_data <- inner_join(results, site_meta, by = "site_no") %>%
  mutate(
    elev_quant = cut(elev_mean_m, breaks = elev_breaks, 
                     include.lowest = TRUE, dig.lab = 4)
  ) %>%
  group_by(elev_quant) %>%
  arrange(days_to_10mm, .by_group = TRUE) %>%
  mutate(ecdf_elev = percent_rank(days_to_10mm)) %>%
  ungroup() %>%
  mutate(
    area_quant = cut(drain_sqkm, breaks = area_breaks, 
                     include.lowest = TRUE, dig.lab = 4)
  ) %>%
  group_by(area_quant) %>%
  arrange(days_to_10mm, .by_group = TRUE) %>%
  mutate(ecdf_area = percent_rank(days_to_10mm)) %>%
  ungroup()

ggplot(plot_data) +
  geom_path(aes(x = days_to_10mm, y = ecdf_elev, col = elev_quant)) +
  theme_bw()

ggplot(plot_data) +
  geom_path(aes(x = days_to_10mm, y = ecdf_area, col = area_quant)) +
  theme_bw()
