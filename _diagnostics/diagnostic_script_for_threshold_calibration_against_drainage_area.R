#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Davenport Threshold Adequacy vs. Drainage Area
# ============================================================================
# Two-panel scatter plot examining how well the Davenport et al. (2020)
# classification thresholds are calibrated as a function of drainage area:
#
#   Panel 1: Median days-to-10mm vs. drainage area
#     - From the snow-free rainfall duration diagnostic (rainfall_duration_diag.rds)
#     - Addresses window length adequacy: do large basins need a longer lookback?
#
#   Panel 2: Median 8-day total rainfall (rainfall-classified peaks) vs. drainage area
#     - From the basin-mean classified peaks (pot_peaks_classified_mean.rds)
#     - Addresses threshold level adequacy: is the 10mm threshold appropriate
#       across the basin size spectrum?
#
# Both panels use log-scale drainage area and include a loess smoother.
# Points are colored by elevation quartile to reveal whether elevation adds
# explanatory power beyond drainage area.
# ============================================================================

library(tidyverse)

# ============================================================================
# CONFIGURATION
# ============================================================================
duration_diag_file  <- "_data/rainfall_duration_diag.rds"
classified_file     <- "_data/pot_peaks_classified_mean.rds"
site_meta_file      <- "_data/site_metadata.rds"
output_plot         <- "_diagnostics/threshold_calibration/diag_threshold_vs_drainage_area.png"

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading data...")

site_meta   <- readRDS(site_meta_file) %>%
  select(site_no, drain_sqkm, elev_mean_m)

duration    <- readRDS(duration_diag_file)

classified  <- readRDS(classified_file) %>%
  filter(flood_type == "rainfall", complete_window)

# ============================================================================
# PANEL 1: Median days-to-10mm per site vs. drainage area
# ============================================================================
panel1_data <- duration %>%
  group_by(site_no) %>%
  summarise(
    median_days_to_10mm = median(days_to_10mm, na.rm = TRUE),
    n_events            = sum(!is.na(days_to_10mm)),
    pct_unreached       = mean(is.na(days_to_10mm)),
    .groups = "drop"
  ) %>%
  filter(n_events >= 5) %>%       # require minimum events for stable median
  inner_join(site_meta, by = "site_no") %>%
  mutate(
    elev_quartile = cut(elev_mean_m,
                        breaks = quantile(elev_mean_m, probs = seq(0, 1, 0.25),
                                          na.rm = TRUE),
                        include.lowest = TRUE,
                        dig.lab = 4,
                        labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"))
  )

p1 <- ggplot(panel1_data,
             aes(x = drain_sqkm, y = median_days_to_10mm, color = elev_quartile)) +
  geom_point(aes(size = n_events), alpha = 0.6) +
  geom_smooth(aes(group = 1), method = "loess", color = "black",
              linewidth = 0.8, se = TRUE, alpha = 0.15) +
  geom_hline(yintercept = 8, linetype = "dashed", color = "red", linewidth = 0.7) +
  annotate("text", x = Inf, y = 8.3, label = "Current 8-day window",
           hjust = 1.05, color = "red", size = 3) +
  scale_x_log10(labels = scales::comma) +
  scale_size_continuous(range = c(1.5, 5), guide = "none") +
  labs(
    title    = "Panel 1: Window Length Adequacy",
    subtitle = "Median days to accumulate 10mm rainfall (snow-free peaks only)",
    x        = "Drainage area (km², log scale)",
    y        = "Median days to 10mm",
    color    = "Elevation quartile"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

# ============================================================================
# PANEL 2: Median 8-day rainfall total (rainfall peaks) vs. drainage area
# ============================================================================
panel2_data <- classified %>%
  group_by(site_no) %>%
  summarise(
    median_rain_mm = median(total_rain_mm, na.rm = TRUE),
    n_events       = n(),
    .groups        = "drop"
  ) %>%
  filter(n_events >= 5) %>%
  inner_join(site_meta, by = "site_no") %>%
  mutate(
    elev_quartile = cut(elev_mean_m,
                        breaks = quantile(elev_mean_m, probs = seq(0, 1, 0.25),
                                          na.rm = TRUE),
                        include.lowest = TRUE,
                        dig.lab = 4,
                        labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"))
  )

p2 <- ggplot(panel2_data,
             aes(x = drain_sqkm, y = median_rain_mm, color = elev_quartile)) +
  geom_point(aes(size = n_events), alpha = 0.6) +
  geom_smooth(aes(group = 1), method = "loess", color = "black",
              linewidth = 0.8, se = TRUE, alpha = 0.15) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red", linewidth = 0.7) +
  annotate("text", x = Inf, y = 10.6, label = "Current 10mm threshold",
           hjust = 1.05, color = "red", size = 3) +
  scale_x_log10(labels = scales::comma) +
  scale_size_continuous(range = c(1.5, 5), guide = "none") +
  labs(
    title    = "Panel 2: Threshold Level Adequacy",
    subtitle = "Median 8-day rainfall total for rainfall-classified peaks",
    x        = "Drainage area (km², log scale)",
    y        = "Median 8-day rainfall (mm)",
    color    = "Elevation quartile"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

# ============================================================================
# COMBINE AND SAVE
# ============================================================================
library(patchwork)

combined <- p1 / p2 +
  plot_annotation(
    title   = "Davenport Threshold Adequacy vs. Drainage Area",
    caption = "Red dashed lines = current Davenport et al. (2020) parameter values.\nLoess smoother with 95% CI. Point size proportional to number of events. Color = elevation quartile."
  )

ggsave(output_plot, combined, width = 9, height = 10, dpi = 150)
message(sprintf("Saved plot to %s", output_plot))

# ============================================================================
# PRINT CORRELATION SUMMARIES
# ============================================================================
message("\n--- Panel 1: Spearman correlation (days_to_10mm ~ drain_sqkm) ---")
cor_p1 <- cor.test(panel1_data$median_days_to_10mm,
                   log10(panel1_data$drain_sqkm),
                   method = "spearman")
message(sprintf("  rho = %.3f, p = %.4f", cor_p1$estimate, cor_p1$p.value))

message("\n--- Panel 2: Spearman correlation (median_rain_mm ~ drain_sqkm) ---")
cor_p2 <- cor.test(panel2_data$median_rain_mm,
                   log10(panel2_data$drain_sqkm),
                   method = "spearman")
message(sprintf("  rho = %.3f, p = %.4f", cor_p2$estimate, cor_p2$p.value))