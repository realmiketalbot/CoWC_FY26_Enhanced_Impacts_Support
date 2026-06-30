#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: 60-Day Accumulation Curves with Data-Driven End-Member ID
# ============================================================================
# Computes day-by-day cumulative rainfall and snowmelt for all POT peaks
# over a 60-day lookback window. End-members are identified statistically
# from the joint distribution of 60-day totals rather than using arbitrary
# thresholds.
#
# Two-stage workflow:
#   Stage 1: Extract 60-day cumulative rain and melt for ALL peaks (long format)
#            Output: one row per peak x lookback_day
#   Stage 2: Identify end-members from joint distribution of 60-day totals,
#            then compare accumulation curves between types
#
# End-member definition:
#   - Computed from the 2D joint distribution of (total_rain_60d, total_melt_60d)
#   - Rain end-members:  high rain / low melt corner
#   - Melt end-members:  high melt / low rain corner
#   - Stringency controlled by ENDMEMBER_PERCENTILE (e.g. 0.25 = bottom quartile
#     of opposing signal among top-half of primary signal)
#
# INPUTS:
#   - _data/pot_peaks_all.rds
#   - _data/WLDAS/processed/daily_by_date/YYYYMMDD_wldas.parquet
#   - _data/site_metadata.rds
#
# OUTPUTS:
#   - _data/accum_curves_60d.rds       : long-format accumulation curves
#                                        (all peaks, rain + melt)
#   - _data/endmembers_60d.rds         : end-member peak IDs with 60d totals
#   - _diagnostics/accum_curves/       : figures
# ============================================================================

library(tidyverse)
library(arrow)
library(patchwork)

# ============================================================================
# CONFIGURATION
# ============================================================================
peaks_file       <- "_data/pot_peaks_all.rds"
wldas_dir        <- "_data/WLDAS/processed/daily_by_date"
site_meta_file   <- "_data/site_metadata.rds"
accum_output     <- "_data/accum_curves_60d.rds"
endmem_output    <- "_data/endmembers_60d.rds"
out_dir          <- "_diagnostics/accum_curves"

LOOKBACK_DAYS       <- 59   # days before peak; + peak day = 60 days total
ENDMEMBER_PCTILE    <- 0.25 # opposing signal must be in bottom X% of its
                             # distribution to qualify as an end-member
PRIMARY_PCTILE      <- 0.50 # primary signal must be in top X% to qualify

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
if (length(parquet_files) == 0) stop(sprintf("No parquet files in %s", wldas_dir))

wldas_dates <- parquet_files %>%
  basename() %>%
  str_extract("\\d{8}") %>%
  as.Date(format = "%Y%m%d") %>%
  sort()

message(sprintf("WLDAS: %s to %s (%d days)",
                min(wldas_dates), max(wldas_dates), length(wldas_dates)))

peaks <- peaks %>%
  mutate(
    window_start    = peak_date - LOOKBACK_DAYS,
    window_end      = peak_date,
    wldas_available = window_start >= min(wldas_dates) &
                      window_end   <= max(wldas_dates)
  )

message(sprintf("Peaks with full 60-day WLDAS coverage: %d of %d",
                sum(peaks$wldas_available), nrow(peaks)))

peaks <- peaks %>% filter(wldas_available)

# ============================================================================
# STAGE 1: LOAD WLDAS AND BUILD ACCUMULATION CURVES
# ============================================================================
# Check for existing output to allow resumption
if (file.exists(accum_output)) {
  message(sprintf("Found existing accumulation output: %s", accum_output))
  message("Loading from cache...")
  accum_long <- readRDS(accum_output)
} else {

  required_dates <- peaks %>%
    rowwise() %>%
    reframe(date = seq(window_start, window_end, by = "day")) %>%
    distinct(date) %>%
    filter(date %in% wldas_dates) %>%
    pull(date)

  message(sprintf("Loading WLDAS data for %d unique dates...",
                  length(required_dates)))

  wldas_data <- map_dfr(required_dates, function(d) {
    f <- file.path(wldas_dir,
                   sprintf("%s_wldas.parquet", format(d, "%Y%m%d")))
    tryCatch(
      read_parquet(f, col_select = c("date", "GAGE_ID",
                                     "rainf_flux", "qsm_flux")) %>%
        mutate(date    = as.Date(date),
               GAGE_ID = as.character(GAGE_ID)),
      error = function(e) {
        warning(sprintf("Could not read %s", basename(f)))
        NULL
      }
    )
  })

  message(sprintf("Loaded %d site-day rows", nrow(wldas_data)))

  # Build per-peak windows with day offset from peak
  message("Building 60-day windows and computing accumulation curves...")

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

  # Compute cumulative sums working backward from peak day
  accum_long <- peak_windows %>%
    group_by(site_no, peak_id, peak_date, peak_flow) %>%
    arrange(lookback_day, .by_group = TRUE) %>%
    mutate(
      cumul_rain_mm = cumsum(replace_na(rainf_flux, 0)),
      cumul_melt_mm = cumsum(replace_na(qsm_flux,   0)),
      n_days_data   = sum(!is.na(rainf_flux))
    ) %>%
    ungroup() %>%
    select(site_no, peak_id, peak_date, peak_flow,
           lookback_day, date,
           rain_mm_day  = rainf_flux,
           melt_mm_day  = qsm_flux,
           cumul_rain_mm, cumul_melt_mm,
           n_days_data)

  message(sprintf("Accumulation curves built: %d peak-day rows", nrow(accum_long)))

  saveRDS(accum_long, accum_output)
  message(sprintf("Saved accumulation curves to %s", accum_output))
}

# ============================================================================
# STAGE 2: DATA-DRIVEN END-MEMBER IDENTIFICATION
# ============================================================================
message("\nIdentifying end-members from joint 60-day totals...")

# Extract 60-day totals per peak
totals_60d <- accum_long %>%
  filter(lookback_day == 60) %>%
  select(site_no, peak_id, peak_date, peak_flow,
         total_rain_60d = cumul_rain_mm,
         total_melt_60d = cumul_melt_mm,
         n_days_data)

# Require complete or near-complete windows (allow up to 3 missing days)
totals_60d <- totals_60d %>%
  filter(n_days_data >= 57)

message(sprintf("Peaks with >= 57 days of data: %d", nrow(totals_60d)))

# Compute percentile cutoffs from the marginal distributions
rain_low_cut  <- quantile(totals_60d$total_rain_60d,
                          ENDMEMBER_PCTILE, na.rm = TRUE)
rain_high_cut <- quantile(totals_60d$total_rain_60d,
                          1 - PRIMARY_PCTILE, na.rm = TRUE)
melt_low_cut  <- quantile(totals_60d$total_melt_60d,
                          ENDMEMBER_PCTILE, na.rm = TRUE)
melt_high_cut <- quantile(totals_60d$total_melt_60d,
                          1 - PRIMARY_PCTILE, na.rm = TRUE)

message(sprintf("\nRain cutoffs:  low (opposing) = %.1f mm, high (primary) = %.1f mm",
                rain_low_cut, rain_high_cut))
message(sprintf("Melt cutoffs:  low (opposing) = %.1f mm, high (primary) = %.1f mm",
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

em_counts <- count(totals_60d, endmember_type)
message("\nEnd-member counts:")
print(em_counts)

# Save end-member table
saveRDS(totals_60d, endmem_output)
message(sprintf("Saved end-member table to %s", endmem_output))

# ============================================================================
# FIGURES
# ============================================================================

# --- Fig 1: Joint distribution of 60-day totals with end-member regions ---
p1 <- ggplot(totals_60d,
             aes(x = total_rain_60d, y = total_melt_60d,
                 color = endmember_type)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_vline(xintercept = c(rain_low_cut, rain_high_cut),
             linetype = c("dashed", "solid"), color = "#2166ac") +
  geom_hline(yintercept = c(melt_low_cut, melt_high_cut),
             linetype = c("dashed", "solid"), color = "#4dac26") +
  annotate("rect",
           xmin = rain_high_cut, xmax = Inf,
           ymin = -Inf,          ymax = melt_low_cut,
           fill = "#2166ac", alpha = 0.08) +
  annotate("rect",
           xmin = -Inf,          xmax = rain_low_cut,
           ymin = melt_high_cut, ymax = Inf,
           fill = "#4dac26", alpha = 0.08) +
  scale_color_manual(values = c(
    rain    = "#2166ac",
    melt    = "#4dac26",
    neither = "grey70"
  )) +
  scale_x_continuous(limits = c(0, quantile(totals_60d$total_rain_60d,
                                            0.99, na.rm = TRUE))) +
  scale_y_continuous(limits = c(0, quantile(totals_60d$total_melt_60d,
                                            0.99, na.rm = TRUE))) +
  labs(
    title    = "Joint distribution of 60-day rainfall and snowmelt totals",
    subtitle = sprintf(
      "Shaded corners = end-member regions (rain top-%d%% / melt bottom-%d%%; melt top-%d%% / rain bottom-%d%%)",
      round(PRIMARY_PCTILE * 100), round(ENDMEMBER_PCTILE * 100),
      round(PRIMARY_PCTILE * 100), round(ENDMEMBER_PCTILE * 100)
    ),
    x     = "Total 60-day rainfall (mm)",
    y     = "Total 60-day snowmelt (mm)",
    color = "End-member type"
  ) +
  theme_bw(base_size = 11)

# --- Fig 2: Median accumulation curves by end-member type ---
em_peaks <- totals_60d %>%
  filter(endmember_type != "neither") %>%
  select(site_no, peak_id, endmember_type)

accum_em <- accum_long %>%
  inner_join(em_peaks, by = c("site_no", "peak_id"))

accum_summary <- accum_em %>%
  pivot_longer(cols = c(cumul_rain_mm, cumul_melt_mm),
               names_to  = "variable",
               values_to = "cumul_mm") %>%
  mutate(variable = if_else(variable == "cumul_rain_mm", "rainfall", "snowmelt")) %>%
  group_by(endmember_type, variable, lookback_day) %>%
  summarise(
    median_cumul = median(cumul_mm, na.rm = TRUE),
    q25_cumul    = quantile(cumul_mm, 0.25, na.rm = TRUE),
    q75_cumul    = quantile(cumul_mm, 0.75, na.rm = TRUE),
    .groups      = "drop"
  )

p2 <- ggplot(accum_summary,
             aes(x = lookback_day, y = median_cumul,
                 color    = interaction(endmember_type, variable),
                 fill     = interaction(endmember_type, variable),
                 linetype = variable)) +
  geom_ribbon(aes(ymin = q25_cumul, ymax = q75_cumul),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 8, linetype = "dashed",
             color = "grey30", linewidth = 0.7) +
  geom_hline(yintercept = 10, linetype = "dotted",
             color = "grey30", linewidth = 0.7) +
  annotate("text", x = 8.3, y = Inf,
           label = "8-day window", hjust = 0, vjust = 1.5,
           color = "grey30", size = 3) +
  annotate("text", x = Inf, y = 10.5,
           label = "10mm threshold", hjust = 1.05,
           color = "grey30", size = 3) +
  scale_color_manual(values = c(
    "rain.rainfall"   = "#2166ac",
    "rain.snowmelt"   = "#92c5de",
    "melt.rainfall"   = "#a6d96a",
    "melt.snowmelt"   = "#4dac26"
  )) +
  scale_fill_manual(values = c(
    "rain.rainfall"   = "#2166ac",
    "rain.snowmelt"   = "#92c5de",
    "melt.rainfall"   = "#a6d96a",
    "melt.snowmelt"   = "#4dac26"
  )) +
  scale_x_continuous(breaks = c(1, 8, 15, 20, 30, 45, 60)) +
  facet_wrap(~ endmember_type, ncol = 1,
             labeller = labeller(endmember_type = c(
               rain = "Rain end-members",
               melt = "Melt end-members"
             ))) +
  labs(
    title    = "60-day accumulation curves by end-member type",
    subtitle = "Median ± IQR. Solid = rainfall, dashed = snowmelt.\nVertical line = current 8-day window; horizontal = 10mm Davenport threshold.",
    x        = "Lookback days (1 = peak day)",
    y        = "Cumulative (mm)",
    color    = "End-member × variable",
    fill     = "End-member × variable",
    linetype = "Variable"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# --- Fig 3: Depth achieved at day 8, 30, and 60 for each end-member type ---
depth_comparison <- accum_long %>%
  filter(lookback_day %in% c(8, 30, 60)) %>%
  inner_join(em_peaks, by = c("site_no", "peak_id")) %>%
  select(site_no, peak_id, endmember_type, lookback_day,
         cumul_rain_mm, cumul_melt_mm) %>%
  pivot_longer(cols = c(cumul_rain_mm, cumul_melt_mm),
               names_to = "variable", values_to = "cumul_mm") %>%
  mutate(
    variable     = if_else(variable == "cumul_rain_mm", "rainfall", "snowmelt"),
    lookback_day = factor(lookback_day,
                          levels = c(8, 30, 60),
                          labels = c("8-day", "30-day", "60-day"))
  )

p3 <- ggplot(depth_comparison,
             aes(x = lookback_day, y = cumul_mm,
                 fill = interaction(endmember_type, variable))) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3, linewidth = 0.4) +
  geom_hline(yintercept = 10, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  annotate("text", x = 0.5, y = 11,
           label = "10mm", hjust = 0, color = "red", size = 3) +
  scale_fill_manual(values = c(
    "rain.rainfall"   = "#2166ac",
    "rain.snowmelt"   = "#92c5de",
    "melt.rainfall"   = "#a6d96a",
    "melt.snowmelt"   = "#4dac26"
  )) +
  facet_grid(variable ~ endmember_type,
             labeller = labeller(
               endmember_type = c(rain = "Rain end-members",
                                  melt = "Melt end-members"),
               variable       = c(rainfall = "Rainfall signal",
                                  snowmelt = "Snowmelt signal")
             )) +
  labs(
    title    = "Signal depth at 8-day, 30-day, and 60-day windows",
    subtitle = "How much additional signal is captured by extending the window?",
    x        = "Window length",
    y        = "Cumulative depth (mm)",
    fill     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# ============================================================================
# COMBINE AND SAVE
# ============================================================================
combined <- p1 / p2 / p3 +
  plot_annotation(
    title   = "60-Day Accumulation Diagnostic",
    caption = sprintf(
      "End-members: primary signal >= %.0fth percentile, opposing signal <= %.0fth percentile of 60-day totals.",
      (1 - PRIMARY_PCTILE) * 100, ENDMEMBER_PCTILE * 100
    )
  )

out_file <- file.path(out_dir, "accum_curves_60d.png")
ggsave(out_file, combined, width = 10, height = 16, dpi = 150)
message(sprintf("\nSaved figure to %s", out_file))

# ============================================================================
# CONSOLE SUMMARY
# ============================================================================
message("\n=== 60-DAY ACCUMULATION SUMMARY ===")

summary_table <- accum_long %>%
  filter(lookback_day %in% c(8, 14, 21, 30, 45, 60)) %>%
  inner_join(em_peaks, by = c("site_no", "peak_id")) %>%
  group_by(endmember_type, lookback_day) %>%
  summarise(
    median_rain = median(cumul_rain_mm, na.rm = TRUE),
    median_melt = median(cumul_melt_mm, na.rm = TRUE),
    pct_rain_above_10mm = round(100 * mean(cumul_rain_mm > 10, na.rm = TRUE), 1),
    pct_melt_above_10mm = round(100 * mean(cumul_melt_mm > 10, na.rm = TRUE), 1),
    .groups = "drop"
  )

message("\nMedian cumulative depths and % above 10mm at key window lengths:")
print(summary_table)