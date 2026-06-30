#!/usr/bin/env Rscript
# =============================================================================
# diagnose-unclassified-peaks.R
#
# For a specified site, explains exactly why each unclassified peak failed
# to meet any of the Davenport et al. (2020) classification criteria, and
# shows how close each peak came to qualifying.
#
# Usage:
#   Rscript diagnose-unclassified-peaks.R 09085100
#   source("diagnose-unclassified-peaks.R")   # uses SITE_ID below
# =============================================================================

library(dplyr)
library(ggplot2)
library(scales)

SITE_ID   <- "09085100"
DATA_FILE <- "_data/pot_peaks_classified.rds"
OUT_FILE  <- NULL  # set to e.g. "diag_09085100.png" to save plot

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) SITE_ID  <- args[1]
if (length(args) >= 2) OUT_FILE <- args[2]

# Davenport et al. 2020 thresholds (must match classify script)
thresh_rain_high <- 10
thresh_rain_low  <-  5
thresh_melt_high <- 10
thresh_melt_low  <-  5

# =============================================================================
# LOAD DATA
# =============================================================================
pot <- readRDS(DATA_FILE)

site_peaks <- pot %>%
  filter(site_no == SITE_ID, complete_window) %>%
  arrange(peak_date)

if (nrow(site_peaks) == 0) {
  stop(sprintf("No peaks with complete WLDAS windows found for site %s", SITE_ID))
}

cat(sprintf("\n══ Site: %s ══\n", SITE_ID))
cat(sprintf("Total peaks with complete windows: %d\n", nrow(site_peaks)))
cat(sprintf("Classified:   %d\n", sum(site_peaks$flood_type != "unclassified")))
cat(sprintf("Unclassified: %d\n\n", sum(site_peaks$flood_type == "unclassified")))

# =============================================================================
# EXPLAIN WHY EACH PEAK WAS UNCLASSIFIED
# =============================================================================
# For each unclassified peak, diagnose which criterion it came closest to
# meeting and by how much it missed.

unclassified <- site_peaks %>%
  filter(flood_type == "unclassified") %>%
  mutate(
    # How far each variable is from its qualifying threshold (positive = deficit)
    rain_deficit_for_rainfall  = thresh_rain_high - total_rain_mm,  # need rain > 10
    melt_deficit_for_snowmelt  = thresh_melt_high - total_melt_mm,  # need melt > 10
    rain_deficit_for_ros       = thresh_rain_low  - total_rain_mm,  # need rain > 5
    melt_deficit_for_ros       = thresh_melt_low  - total_melt_mm,  # need melt > 5

    # Would qualify for rainfall if: rain > 10 AND melt < 5
    # Already fails rain threshold; also check melt excess
    melt_excess_for_rainfall   = total_melt_mm - thresh_melt_low,   # need melt < 5

    # Would qualify for snowmelt if: melt > 10 AND rain < 5
    rain_excess_for_snowmelt   = total_rain_mm - thresh_rain_low,   # need rain < 5

    # Reason labels
    reason = case_when(
      # Passes melt threshold but too much rain for snowmelt
      total_melt_mm > thresh_melt_high & total_rain_mm >= thresh_rain_low ~
        sprintf("Near snowmelt: melt=%.1f OK, but rain=%.1f >= %.0f (need < %.0f)",
                total_melt_mm, total_rain_mm, thresh_rain_low, thresh_rain_low),

      # Passes rain threshold but too much melt for rainfall
      total_rain_mm > thresh_rain_high & total_melt_mm >= thresh_melt_low ~
        sprintf("Near rainfall: rain=%.1f OK, but melt=%.1f >= %.0f (need < %.0f)",
                total_rain_mm, total_melt_mm, thresh_melt_low, thresh_melt_low),

      # Both rain and melt are substantial but neither hits the ROS threshold pair
      total_rain_mm >= thresh_rain_low & total_melt_mm >= thresh_melt_low ~
        sprintf("Near ROS: rain=%.1f OK, melt=%.1f OK — but classified as unclassified? Check logic.",
                total_rain_mm, total_melt_mm),

      # Both variables are low — genuinely dry/ambiguous event
      total_rain_mm < thresh_rain_low & total_melt_mm < thresh_melt_low ~
        sprintf("Low signal: rain=%.1f < %.0f, melt=%.1f < %.0f — ambiguous event",
                total_rain_mm, thresh_rain_low, total_melt_mm, thresh_melt_low),

      # Rain nearly enough for rainfall, melt too high
      rain_deficit_for_rainfall < 2 ~
        sprintf("Just missed rainfall: rain=%.1f (need > %.0f), melt=%.1f",
                total_rain_mm, thresh_rain_high, total_melt_mm),

      # Melt nearly enough for snowmelt, rain too high
      melt_deficit_for_snowmelt < 2 ~
        sprintf("Just missed snowmelt: melt=%.1f (need > %.0f), rain=%.1f",
                total_melt_mm, thresh_melt_high, total_rain_mm),

      TRUE ~
        sprintf("Both low: rain=%.1f, melt=%.1f",
                total_rain_mm, total_melt_mm)
    )
  )

cat("── Per-peak diagnosis ──\n")
cat(sprintf("%-12s  %-6s  %8s  %8s  %s\n",
            "Date", "Flow", "Rain(mm)", "Melt(mm)", "Reason"))
cat(strrep("-", 90), "\n")
for (i in seq_len(nrow(unclassified))) {
  r <- unclassified[i, ]
  cat(sprintf("%-12s  %6.0f  %8.1f  %8.1f  %s\n",
              format(r$peak_date), r$peak_flow,
              r$total_rain_mm, r$total_melt_mm, r$reason))
}

# =============================================================================
# SUMMARY OF FAILURE MODES
# =============================================================================
cat("\n── Failure mode summary ──\n")
failure_modes <- unclassified %>%
  mutate(mode = case_when(
    total_melt_mm > thresh_melt_high & total_rain_mm >= thresh_rain_low ~ "Melt OK, rain too high for snowmelt",
    total_rain_mm > thresh_rain_high & total_melt_mm >= thresh_melt_low ~ "Rain OK, melt too high for rainfall",
    total_rain_mm < thresh_rain_low  & total_melt_mm < thresh_melt_low  ~ "Both below thresholds (low signal)",
    TRUE ~ "Other"
  )) %>%
  count(mode) %>%
  mutate(pct = round(100 * n / sum(n), 1))
print(failure_modes)

cat("\n── Distribution of rain and melt values ──\n")
cat(sprintf("  Rain (mm): min=%.1f  median=%.1f  mean=%.1f  max=%.1f\n",
            min(unclassified$total_rain_mm),
            median(unclassified$total_rain_mm),
            mean(unclassified$total_rain_mm),
            max(unclassified$total_rain_mm)))
cat(sprintf("  Melt (mm): min=%.1f  median=%.1f  mean=%.1f  max=%.1f\n",
            min(unclassified$total_melt_mm),
            median(unclassified$total_melt_mm),
            mean(unclassified$total_melt_mm),
            max(unclassified$total_melt_mm)))

# How many peaks are within 2mm of qualifying for any type?
near_misses <- unclassified %>%
  filter(
    (total_rain_mm > thresh_rain_high - 2 & total_melt_mm < thresh_melt_low) |  # near rainfall
    (total_melt_mm > thresh_melt_high - 2 & total_rain_mm < thresh_rain_low) |  # near snowmelt
    (total_rain_mm > thresh_rain_low  - 2 & total_melt_mm > thresh_melt_low - 2) # near ROS
  )
cat(sprintf("\n  Peaks within 2mm of qualifying for any type: %d / %d\n",
            nrow(near_misses), nrow(unclassified)))

# =============================================================================
# SCATTER PLOT: rain vs melt coloured by flood type, threshold lines shown
# =============================================================================
# All classified peaks shown for context alongside the unclassified ones.

plot_data <- site_peaks %>%
  mutate(type_label = case_when(
    flood_type == "rainfall"     ~ "Rainfall",
    flood_type == "snowmelt"     ~ "Snowmelt",
    flood_type == "rain_on_snow" ~ "Rain-on-Snow",
    TRUE                         ~ "Unclassified"
  ))

type_colors <- c(
  "Rainfall"     = "#e41a1c",
  "Snowmelt"     = "#377eb8",
  "Rain-on-Snow" = "#4daf4a",
  "Unclassified" = "#aaaaaa"
)

p <- ggplot(plot_data, aes(x = total_rain_mm, y = total_melt_mm,
                           colour = type_label, size = type_label)) +
  # Threshold lines
  geom_vline(xintercept = thresh_rain_low,  linetype = "dashed", colour = "black", linewidth = 0.4) +
  geom_vline(xintercept = thresh_rain_high, linetype = "solid",  colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = thresh_melt_low,  linetype = "dashed", colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = thresh_melt_high, linetype = "solid",  colour = "black", linewidth = 0.4) +
  # Points
  geom_point(alpha = 0.8) +
  scale_colour_manual(values = type_colors, name = "Flood type") +
  scale_size_manual(values = c("Unclassified" = 2, "Rainfall" = 2.5,
                               "Snowmelt" = 2.5, "Rain-on-Snow" = 2.5),
                    guide = "none") +
  # Annotate threshold zones
  annotate("text", x = thresh_rain_high + 1, y = max(plot_data$total_melt_mm) * 0.95,
           label = "rain > 10", hjust = 0, size = 3, colour = "black") +
  annotate("text", x = max(plot_data$total_rain_mm) * 0.95,
           y = thresh_melt_high + 0.5,
           label = "melt > 10", hjust = 1, size = 3, colour = "black") +
  labs(
    title   = sprintf("%s — Classification space (rain vs. melt)", SITE_ID),
    x       = "Total rain over 8-day window (mm)",
    y       = "Total melt (qsm_flux) over 8-day window (mm)",
    caption = "Solid lines: high thresholds (10mm). Dashed lines: low thresholds (5mm)."
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

if (!is.null(OUT_FILE)) {
  ggsave(OUT_FILE, p, width = 7, height = 6, dpi = 150)
  message(sprintf("Plot saved to %s", OUT_FILE))
} else {
  print(p)
}