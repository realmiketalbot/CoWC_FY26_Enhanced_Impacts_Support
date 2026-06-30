#!/usr/bin/env Rscript

# ============================================================================
# Threshold Calibration Diagnostic: End-Member Peak Analysis
# ============================================================================
# Identifies "end-member" flood peaks where classification is unambiguous
# regardless of threshold choice:
#
#   Rain end-members:  total_melt_mm == 0 (strict) AND total_rain_mm > 5 mm
#   Melt end-members:  total_rain_mm == 0 (strict) AND total_melt_mm > 5 mm
#
# For these peaks, compares basin-mean antecedent forcing to:
#   (a) max-cell values, to quantify spatial dilution
#   (b) Davenport et al. thresholds, to assess whether basin-mean values
#       fall below classification thresholds despite unambiguous forcing
#
# The "zero" condition uses a small tolerance (< 0.1 mm cumulative over the
# 8-day window) rather than strict zero to account for floating-point noise
# in WLDAS output.
#
# INPUTS:
#   - pot_peaks_classified.rds   : output of script 5 (includes flood_type,
#                                  total_rain_mm, total_melt_mm, window totals)
#   - YYYYMMDD_wldas.parquet     : daily WLDAS extractions including spatial
#                                  statistics (_max, _sd, _cv, _wet_frac)
#   - site_metadata.rds          : site attributes (elevation, drain area, etc.)
#
# OUTPUTS:
#   - _diagnostics/threshold_calibration/fig01_endmember_counts.png
#   - _diagnostics/threshold_calibration/fig02_rain_dilution.png
#   - _diagnostics/threshold_calibration/fig03_melt_dilution.png
#   - _diagnostics/threshold_calibration/fig04_rain_vs_threshold.png
#   - _diagnostics/threshold_calibration/fig05_melt_vs_threshold.png
#   - _diagnostics/threshold_calibration/fig06_dilution_by_area.png
#   - _diagnostics/threshold_calibration/fig07_dilution_by_elevation.png
#   - _diagnostics/threshold_calibration/endmember_peaks.rds
# ============================================================================

library(tidyverse)
library(arrow)

# ============================================================================
# CONFIGURATION
# ============================================================================
peaks_file   <- "_data/pot_peaks_classified.rds"
wldas_dir    <- "_data/WLDAS/processed/daily_by_date"
meta_file    <- "_data/site_metadata.rds"
out_dir      <- "_diagnostics/threshold_calibration"

# End-member definition
ABSENT_THRESH_MM  <- 0.1   # cumulative mm below which variable is "absent"
PRESENT_THRESH_MM <- 5.0   # cumulative mm above which variable is "present"

# Davenport thresholds for reference lines
THRESH_RAIN_HIGH  <- 10
THRESH_RAIN_LOW   <-  5
THRESH_MELT_HIGH  <- 10
THRESH_MELT_LOW   <-  5
ANTECEDENT_DAYS   <-  7

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading classified peaks...")
peaks <- readRDS(peaks_file) %>%
  mutate(peak_date = as.Date(peak_date))

message(sprintf("Loaded %d peaks from %d sites",
                nrow(peaks), n_distinct(peaks$site_no)))

message("Loading site metadata...")
meta <- readRDS(meta_file)

# ============================================================================
# IDENTIFY END-MEMBER PEAKS
# ============================================================================
# Note: total_rain_mm and total_melt_mm in classified peaks are basin-mean
# cumulative values over the 8-day window (already in mm). We use these
# directly for end-member identification since "zero melt" is a model output
# property independent of threshold choice.

message("Identifying end-member peaks...")

# Rain end-members: negligible melt, meaningful rain
rain_endmembers <- peaks %>%
  filter(
    complete_window,
    total_melt_mm < ABSENT_THRESH_MM,
    total_rain_mm >= PRESENT_THRESH_MM
  ) %>%
  mutate(endmember_type = "rain")

# Melt end-members: negligible rain, meaningful melt
melt_endmembers <- peaks %>%
  filter(
    complete_window,
    total_rain_mm < ABSENT_THRESH_MM,
    total_melt_mm >= PRESENT_THRESH_MM
  ) %>%
  mutate(endmember_type = "melt")

endmembers <- bind_rows(rain_endmembers, melt_endmembers)

message(sprintf("Rain end-members: %d peaks from %d sites",
                nrow(rain_endmembers), n_distinct(rain_endmembers$site_no)))
message(sprintf("Melt end-members: %d peaks from %d sites",
                nrow(melt_endmembers), n_distinct(melt_endmembers$site_no)))
message(sprintf("  Rain end-members already classified 'rainfall':    %d (%.1f%%)",
                sum(rain_endmembers$flood_type == "rainfall"),
                100 * mean(rain_endmembers$flood_type == "rainfall")))
message(sprintf("  Rain end-members classified 'unclassified':        %d (%.1f%%)",
                sum(rain_endmembers$flood_type == "unclassified"),
                100 * mean(rain_endmembers$flood_type == "unclassified")))
message(sprintf("  Melt end-members already classified 'snowmelt':    %d (%.1f%%)",
                sum(melt_endmembers$flood_type == "snowmelt"),
                100 * mean(melt_endmembers$flood_type == "snowmelt")))
message(sprintf("  Melt end-members classified 'unclassified':        %d (%.1f%%)",
                sum(melt_endmembers$flood_type == "unclassified"),
                100 * mean(melt_endmembers$flood_type == "unclassified")))

# ============================================================================
# LOAD SPATIAL STATISTICS FOR END-MEMBER PEAKS
# ============================================================================
# Re-load antecedent windows from Parquet, this time including _max columns
# so we can compare basin-mean to max-cell values.

message("Loading WLDAS spatial statistics for end-member peaks...")

required_dates <- endmembers %>%
  rowwise() %>%
  reframe(date = seq(peak_date - ANTECEDENT_DAYS, peak_date, by = "day")) %>%
  distinct(date) %>%
  pull(date)

message(sprintf("Loading %d unique WLDAS dates...", length(required_dates)))

wldas_spatial <- map_dfr(required_dates, function(d) {
  f <- file.path(wldas_dir, sprintf("%s_wldas.parquet", format(d, "%Y%m%d")))
  tryCatch(
    read_parquet(f, col_select = c(
      "date", "GAGE_ID",
      "rainf_flux", "rainf_flux_max", "rainf_flux_wet_frac",
      "qsm_flux",   "qsm_flux_max",   "qsm_flux_wet_frac"
    )) %>%
      mutate(date = as.Date(date), GAGE_ID = as.character(GAGE_ID)),
    error = function(e) NULL
  )
})

message(sprintf("Loaded %d site-day rows with spatial statistics", nrow(wldas_spatial)))

# ============================================================================
# COMPUTE WINDOW TOTALS FOR END-MEMBER PEAKS (MEAN AND MAX)
# ============================================================================
message("Computing antecedent window totals (mean and max)...")

# Join end-members to basin metadata to get GAGE_ID
endmembers_with_gageid <- endmembers %>%
  mutate(GAGE_ID = site_no)  # GAGE_ID == site_no for CO GAGES-II

# Expand each peak to its full antecedent window
endmember_windows <- endmembers_with_gageid %>%
  select(site_no, GAGE_ID, peak_id, peak_date, endmember_type,
         flood_type, total_rain_mm, total_melt_mm) %>%
  rowwise() %>%
  reframe(
    site_no        = site_no,
    GAGE_ID        = GAGE_ID,
    peak_id        = peak_id,
    peak_date      = peak_date,
    endmember_type = endmember_type,
    flood_type     = flood_type,
    total_rain_mm_mean = total_rain_mm,
    total_melt_mm_mean = total_melt_mm,
    date           = seq(peak_date - ANTECEDENT_DAYS, peak_date, by = "day")
  ) %>%
  left_join(wldas_spatial, by = c("GAGE_ID", "date"))

# Summarise to per-peak window totals for mean and max
endmember_totals <- endmember_windows %>%
  group_by(site_no, peak_id, peak_date, endmember_type,
           flood_type, total_rain_mm_mean, total_melt_mm_mean) %>%
  summarise(
    # Basin-mean cumulative (should match total_rain_mm / total_melt_mm from script 5)
    rain_mean_cumul     = sum(rainf_flux,         na.rm = TRUE),
    melt_mean_cumul     = sum(qsm_flux,           na.rm = TRUE),
    # Max-cell cumulative (sum of daily max-cell values over window)
    rain_max_cumul      = sum(rainf_flux_max,     na.rm = TRUE),
    melt_max_cumul      = sum(qsm_flux_max,       na.rm = TRUE),
    # Mean wetted fraction over window (for days with any rain/melt)
    rain_wet_frac_mean  = mean(rainf_flux_wet_frac[rainf_flux > 0],  na.rm = TRUE),
    melt_wet_frac_mean  = mean(qsm_flux_wet_frac[qsm_flux > 0],     na.rm = TRUE),
    n_days              = sum(!is.na(rainf_flux)),
    .groups = "drop"
  ) %>%
  mutate(
    # Dilution ratio: how much does basin-mean understate the max-cell value?
    rain_dilution_ratio = rain_mean_cumul / rain_max_cumul,
    melt_dilution_ratio = melt_mean_cumul / melt_max_cumul,
    # Classification outcome relative to Davenport thresholds
    rain_above_high     = rain_mean_cumul > THRESH_RAIN_HIGH,
    rain_above_low      = rain_mean_cumul > THRESH_RAIN_LOW,
    melt_above_high     = melt_mean_cumul > THRESH_MELT_HIGH,
    melt_above_low      = melt_mean_cumul > THRESH_MELT_LOW,
    # Would max-cell values change the classification?
    rain_max_above_high = rain_max_cumul  > THRESH_RAIN_HIGH,
    melt_max_above_high = melt_max_cumul  > THRESH_MELT_HIGH
  )

# Join site metadata
endmember_totals <- endmember_totals %>%
  left_join(meta %>% select(site_no, elev_mean_m, drain_sqkm, huc02),
            by = "site_no")

message(sprintf("End-member totals computed for %d peaks", nrow(endmember_totals)))

# Quick sanity check: mean_cumul should match total_*_mm from script 5
rain_check <- cor(endmember_totals$rain_mean_cumul,
                  endmember_totals$total_rain_mm_mean, use = "complete.obs")
melt_check <- cor(endmember_totals$melt_mean_cumul,
                  endmember_totals$total_melt_mm_mean, use = "complete.obs")
message(sprintf("Sanity check correlations (should be ~1.0): rain=%.4f, melt=%.4f",
                rain_check, melt_check))

# ============================================================================
# FIG 01: End-member counts and classification outcome
# ============================================================================
message("Generating figures...")

fig01_data <- endmember_totals %>%
  mutate(
    outcome = case_when(
      endmember_type == "rain" & flood_type == "rainfall"     ~ "Correctly classified",
      endmember_type == "melt" & flood_type == "snowmelt"     ~ "Correctly classified",
      flood_type == "unclassified"                             ~ "Unclassified (sample size\nor threshold failure)",
      TRUE                                                     ~ "Misclassified"
    )
  ) %>%
  count(endmember_type, outcome)

fig01 <- ggplot(fig01_data,
                aes(x = endmember_type, y = n, fill = outcome)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c(
    "Correctly classified"                      = "#2166ac",
    "Unclassified (sample size\nor threshold failure)" = "#d73027",
    "Misclassified"                             = "#f46d43"
  )) +
  labs(title    = "End-member peak classification outcomes",
       subtitle = sprintf("Rain: melt < %.1f mm, rain >= %.1f mm  |  Melt: rain < %.1f mm, melt >= %.1f mm",
                          ABSENT_THRESH_MM, PRESENT_THRESH_MM,
                          ABSENT_THRESH_MM, PRESENT_THRESH_MM),
       x = "End-member type", y = "Number of peaks", fill = NULL) +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "fig01_endmember_counts.png"),
       fig01, width = 7, height = 5, dpi = 150)

# ============================================================================
# FIG 02: Rain dilution — basin-mean vs max-cell cumulative
# ============================================================================
rain_em <- endmember_totals %>% filter(endmember_type == "rain")

fig02 <- ggplot(rain_em, aes(x = rain_max_cumul, y = rain_mean_cumul)) +
  geom_point(alpha = 0.3, size = 1.2, colour = "#2166ac") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = THRESH_RAIN_HIGH, colour = "#d73027",
             linetype = "solid", linewidth = 0.7) +
  geom_hline(yintercept = THRESH_RAIN_LOW, colour = "#f46d43",
             linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = max(rain_em$rain_max_cumul, na.rm = TRUE) * 0.95,
           y = THRESH_RAIN_HIGH + 1, label = "10mm threshold",
           hjust = 1, colour = "#d73027", size = 3) +
  annotate("text", x = max(rain_em$rain_max_cumul, na.rm = TRUE) * 0.95,
           y = THRESH_RAIN_LOW + 1, label = "5mm threshold",
           hjust = 1, colour = "#f46d43", size = 3) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(title    = "Rain end-members: basin-mean vs max-cell cumulative rainfall",
       subtitle = "Points below dashed 1:1 line indicate spatial dilution\nPoints below red line would be unclassified by basin-mean",
       x = "Max-cell cumulative rain (mm, 8-day window)",
       y = "Basin-mean cumulative rain (mm, 8-day window)") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "fig02_rain_dilution.png"),
       fig02, width = 7, height = 6, dpi = 150)

# ============================================================================
# FIG 03: Melt dilution — basin-mean vs max-cell cumulative
# ============================================================================
melt_em <- endmember_totals %>% filter(endmember_type == "melt")

fig03 <- ggplot(melt_em, aes(x = melt_max_cumul, y = melt_mean_cumul)) +
  geom_point(alpha = 0.3, size = 1.2, colour = "#4dac26") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = THRESH_MELT_HIGH, colour = "#d73027",
             linetype = "solid", linewidth = 0.7) +
  geom_hline(yintercept = THRESH_MELT_LOW, colour = "#f46d43",
             linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = max(melt_em$melt_max_cumul, na.rm = TRUE) * 0.95,
           y = THRESH_MELT_HIGH + 1, label = "10mm threshold",
           hjust = 1, colour = "#d73027", size = 3) +
  annotate("text", x = max(melt_em$melt_max_cumul, na.rm = TRUE) * 0.95,
           y = THRESH_MELT_LOW + 1, label = "5mm threshold",
           hjust = 1, colour = "#f46d43", size = 3) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(title    = "Melt end-members: basin-mean vs max-cell cumulative snowmelt",
       subtitle = "Points below dashed 1:1 line indicate spatial dilution\nPoints below red line would be unclassified by basin-mean",
       x = "Max-cell cumulative snowmelt (mm, 8-day window)",
       y = "Basin-mean cumulative snowmelt (mm, 8-day window)") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "fig03_melt_dilution.png"),
       fig03, width = 7, height = 6, dpi = 150)

# ============================================================================
# FIG 04: Rain end-members — distribution of basin-mean relative to thresholds
# ============================================================================
fig04 <- ggplot(rain_em,
                aes(x = rain_mean_cumul, fill = flood_type)) +
  geom_histogram(binwidth = 2, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = THRESH_RAIN_HIGH, colour = "#d73027",
             linewidth = 0.8) +
  geom_vline(xintercept = THRESH_RAIN_LOW, colour = "#f46d43",
             linewidth = 0.8, linetype = "dashed") +
  scale_fill_manual(values = c(
    "rainfall"     = "#2166ac",
    "unclassified" = "#d73027",
    "rain_on_snow" = "#762a83"
  )) +
  labs(title    = "Rain end-members: distribution of basin-mean cumulative rainfall",
       subtitle = "Vertical lines show Davenport 10mm (solid) and 5mm (dashed) thresholds",
       x = "Basin-mean cumulative rain (mm, 8-day window)",
       y = "Number of peaks",
       fill = "Classification") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "fig04_rain_vs_threshold.png"),
       fig04, width = 8, height = 5, dpi = 150)

# ============================================================================
# FIG 05: Melt end-members — distribution of basin-mean relative to thresholds
# ============================================================================
fig05 <- ggplot(melt_em,
                aes(x = melt_mean_cumul, fill = flood_type)) +
  geom_histogram(binwidth = 2, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = THRESH_MELT_HIGH, colour = "#d73027",
             linewidth = 0.8) +
  geom_vline(xintercept = THRESH_MELT_LOW, colour = "#f46d43",
             linewidth = 0.8, linetype = "dashed") +
  scale_fill_manual(values = c(
    "snowmelt"     = "#4dac26",
    "unclassified" = "#d73027",
    "rain_on_snow" = "#762a83"
  )) +
  labs(title    = "Melt end-members: distribution of basin-mean cumulative snowmelt",
       subtitle = "Vertical lines show Davenport 10mm (solid) and 5mm (dashed) thresholds",
       x = "Basin-mean cumulative snowmelt (mm, 8-day window)",
       y = "Number of peaks",
       fill = "Classification") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "fig05_melt_vs_threshold.png"),
       fig05, width = 8, height = 5, dpi = 150)

# ============================================================================
# FIG 06: Dilution ratio by drainage area
# ============================================================================
fig06_data <- endmember_totals %>%
  select(endmember_type, drain_sqkm, elev_mean_m, rain_dilution_ratio, melt_dilution_ratio) %>%
  pivot_longer(cols = c(rain_dilution_ratio, melt_dilution_ratio),
               names_to = "variable", values_to = "dilution_ratio") %>%
  filter(
    (endmember_type == "rain" & variable == "rain_dilution_ratio") |
    (endmember_type == "melt" & variable == "melt_dilution_ratio")
  ) %>%
  filter(!is.na(dilution_ratio), is.finite(dilution_ratio),
         dilution_ratio > 0, dilution_ratio <= 1)

fig06 <- ggplot(fig06_data,
                aes(x = drain_sqkm, y = dilution_ratio,
                    colour = endmember_type)) +
  geom_point(alpha = 0.3, size = 1.2) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_x_log10(labels = scales::comma) +
  scale_colour_manual(values = c("rain" = "#2166ac", "melt" = "#4dac26")) +
  labs(title    = "Spatial dilution ratio vs drainage area",
       subtitle = "Ratio = basin-mean / max-cell cumulative (lower = more dilution)",
       x = "Drainage area (km², log scale)",
       y = "Dilution ratio (basin-mean / max-cell)",
       colour = "End-member type") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "fig06_dilution_by_area.png"),
       fig06, width = 8, height = 5, dpi = 150)

# ============================================================================
# FIG 07: Dilution ratio by elevation
# ============================================================================
fig07 <- ggplot(fig06_data,
                aes(x = elev_mean_m, y = dilution_ratio,
                    colour = endmember_type)) +
  geom_point(alpha = 0.3, size = 1.2) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("rain" = "#2166ac", "melt" = "#4dac26")) +
  labs(title    = "Spatial dilution ratio vs mean basin elevation",
       subtitle = "Ratio = basin-mean / max-cell cumulative (lower = more dilution)",
       x = "Mean basin elevation (m)",
       y = "Dilution ratio (basin-mean / max-cell)",
       colour = "End-member type") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "fig07_dilution_by_elevation.png"),
       fig07, width = 8, height = 5, dpi = 150)

# ============================================================================
# CONSOLE SUMMARY
# ============================================================================
message("\n=== THRESHOLD CALIBRATION SUMMARY ===")

message("\nRain end-members (n = ", nrow(rain_em), "):")
message(sprintf("  Basin-mean cumul rain: median=%.1f mm, mean=%.1f mm",
                median(rain_em$rain_mean_cumul, na.rm = TRUE),
                mean(rain_em$rain_mean_cumul, na.rm = TRUE)))
message(sprintf("  Max-cell cumul rain:   median=%.1f mm, mean=%.1f mm",
                median(rain_em$rain_max_cumul, na.rm = TRUE),
                mean(rain_em$rain_max_cumul, na.rm = TRUE)))
message(sprintf("  Peaks below 10mm threshold (basin-mean): %d (%.1f%%)",
                sum(!rain_em$rain_above_high, na.rm = TRUE),
                100 * mean(!rain_em$rain_above_high, na.rm = TRUE)))
message(sprintf("  Peaks above 10mm threshold (max-cell):   %d (%.1f%%)",
                sum(rain_em$rain_max_above_high, na.rm = TRUE),
                100 * mean(rain_em$rain_max_above_high, na.rm = TRUE)))
message(sprintf("  Median dilution ratio: %.3f",
                median(rain_em$rain_dilution_ratio, na.rm = TRUE)))

message("\nMelt end-members (n = ", nrow(melt_em), "):")
message(sprintf("  Basin-mean cumul melt: median=%.1f mm, mean=%.1f mm",
                median(melt_em$melt_mean_cumul, na.rm = TRUE),
                mean(melt_em$melt_mean_cumul, na.rm = TRUE)))
message(sprintf("  Max-cell cumul melt:   median=%.1f mm, mean=%.1f mm",
                median(melt_em$melt_max_cumul, na.rm = TRUE),
                mean(melt_em$melt_max_cumul, na.rm = TRUE)))
message(sprintf("  Peaks below 10mm threshold (basin-mean): %d (%.1f%%)",
                sum(!melt_em$melt_above_high, na.rm = TRUE),
                100 * mean(!melt_em$melt_above_high, na.rm = TRUE)))
message(sprintf("  Peaks above 10mm threshold (max-cell):   %d (%.1f%%)",
                sum(melt_em$melt_max_above_high, na.rm = TRUE),
                100 * mean(melt_em$melt_max_above_high, na.rm = TRUE)))
message(sprintf("  Median dilution ratio: %.3f",
                median(melt_em$melt_dilution_ratio, na.rm = TRUE)))

# ============================================================================
# SAVE END-MEMBER PEAKS
# ============================================================================
saveRDS(endmember_totals,
        file.path(out_dir, "endmember_peaks.rds"))
message(sprintf("\nSaved %d end-member peaks to %s/endmember_peaks.rds",
                nrow(endmember_totals), out_dir))
message(sprintf("Figures saved to %s/", out_dir))