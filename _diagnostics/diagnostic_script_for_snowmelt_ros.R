#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Snowmelt Events with High ros_flux_frac
# ============================================================================
# In the adaptive classification (5b), events reach "snowmelt" via:
#   rf_long <= SNOW_EM_THRESH (0.05)  [long window is predominantly melt]
# but ros_flux_frac is computed over the SHORT window only. This means an
# event can be classified as snowmelt while still having substantial
# concurrent rainfall in the days immediately preceding the peak.
#
# This is a potential misclassification pathway: an event with a long
# antecedent melt signal but a rain burst right at the peak would be
# physically more consistent with ROS or mixed than pure snowmelt.
#
# This diagnostic characterises the population of such events and asks:
#   1. How many snowmelt-classified events have ros_flux_frac above various
#      thresholds? What fraction of the snowmelt pool does this represent?
#   2. What do their accumulation curves look like — is the rain arriving
#      as a late burst, or is it distributed throughout the short window?
#   3. Are these events distinguishable from "clean" snowmelt events in
#      terms of peak magnitude, season, or basin characteristics?
#   4. What would their classification be if the decision tree evaluated
#      ros_flux_frac BEFORE the rf_long snowmelt criterion (i.e. a reordered
#      decision tree)?
#   5. Accumulation curve plots for a sample of flagged events vs clean
#      snowmelt events, to visually verify the rain-burst hypothesis.
#
# INPUTS:
#   - _data/pot_peaks_classified_adaptive.rds
#   - _data/accum_curves_60d.rds
#   - _data/site_metadata.rds
#
# OUTPUTS (-> _diagnostics/snowmelt_high_ros/):
#   - 01_ros_frac_within_snowmelt.png   : distribution of ros_flux_frac
#                                         within the snowmelt type
#   - 02_peak_magnitude_comparison.png  : peak flow, flagged vs clean snowmelt
#   - 03_seasonal_comparison.png        : wy_doy, flagged vs clean
#   - 04_reclassification_counts.png    : what would reordered tree produce?
#   - 05_accum_curves_flagged.png       : sample accumulation curves,
#                                         high-ros snowmelt events
#   - 06_accum_curves_clean.png         : sample accumulation curves,
#                                         clean snowmelt events (control)
#   - diag_snowmelt_high_ros.txt        : numeric summary
# ============================================================================

library(tidyverse)
library(patchwork)

# ============================================================================
# CONFIGURATION
# ============================================================================
classified_file <- "_data/pot_peaks_classified_adaptive.rds"
accum_file      <- "_data/accum_curves_60d.rds"
meta_file       <- "_data/site_metadata.rds"
out_dir         <- "_diagnostics/snowmelt_high_ros"

# Thresholds from 5b
RAIN_EM_THRESH    <- 0.95
SNOW_EM_THRESH    <- 0.05
ROS_THRESH        <- 0.50
SNOW_MULT         <- 3.0

# For the "flagged" population: ros_flux_frac above this within snowmelt events
# We'll show results at multiple thresholds but use this as the primary split
FLAG_THRESH <- 0.25   # lower than the ROS criterion — catch events with
                      # *any* meaningful co-occurrence in the short window

# Number of example accumulation curves to plot per group
N_CURVE_EXAMPLES <- 16

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD AND SUBSET
# ============================================================================
message("Loading data...")
peaks <- readRDS(classified_file)
meta  <- readRDS(meta_file)
accum <- readRDS(accum_file)

if (!"regulated" %in% names(peaks)) {
  peaks <- peaks %>%
    left_join(meta %>% select(site_no, regulated), by = "site_no")
}

peaks <- peaks %>%
  mutate(
    peak_date = as.Date(peak_date),
    month     = lubridate::month(peak_date),
    wy_doy    = as.integer(
      peak_date - as.Date(paste0(
        if_else(month >= 10L,
                lubridate::year(peak_date),
                lubridate::year(peak_date) - 1L),
        "-10-01"
      )) + 1L
    )
  )

# Snowmelt-classified events (using flood_type_raw so sample-size demotions
# don't confuse the picture — we want the threshold-based classification)
snowmelt <- peaks %>%
  filter(flood_type_raw == "snowmelt", !is.na(ros_flux_frac))

message(sprintf("Snowmelt events (flood_type_raw): %d", nrow(snowmelt)))
message(sprintf("  With non-missing ros_flux_frac: %d", nrow(snowmelt)))

# Flag events with elevated ros_flux_frac
snowmelt <- snowmelt %>%
  mutate(
    flagged = ros_flux_frac >= FLAG_THRESH,
    flag_label = if_else(flagged,
                         sprintf("ros_flux_frac ≥ %.2f", FLAG_THRESH),
                         sprintf("ros_flux_frac < %.2f", FLAG_THRESH))
  )

n_flagged <- sum(snowmelt$flagged)
n_clean   <- sum(!snowmelt$flagged)

message(sprintf("  Flagged (ros_flux_frac >= %.2f): %d (%.1f%%)",
                FLAG_THRESH, n_flagged, 100 * n_flagged / nrow(snowmelt)))
message(sprintf("  Clean   (ros_flux_frac <  %.2f): %d (%.1f%%)",
                FLAG_THRESH, n_clean,   100 * n_clean   / nrow(snowmelt)))

# ============================================================================
# SHARED THEME
# ============================================================================
theme_diag <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    plot.title       = element_text(size = 12, face = "bold"),
    plot.subtitle    = element_text(size = 10, colour = "grey40")
  )

flag_colors <- c("TRUE" = "#d95f02", "FALSE" = "#377eb8")
flag_labels <- c("TRUE"  = sprintf("Flagged (ros_flux_frac ≥ %.2f)", FLAG_THRESH),
                 "FALSE" = sprintf("Clean   (ros_flux_frac <  %.2f)", FLAG_THRESH))

# ============================================================================
# 1. DISTRIBUTION OF ros_flux_frac WITHIN SNOWMELT
# ============================================================================
message("Plot 1: ros_flux_frac distribution within snowmelt...")

# Counts at several thresholds for the summary table
thresh_summary <- tibble(threshold = c(0.10, 0.20, 0.25, 0.30, 0.40, 0.50)) %>%
  mutate(
    n_above = map_int(threshold, ~sum(snowmelt$ros_flux_frac >= .x)),
    pct     = 100 * n_above / nrow(snowmelt)
  )

p1 <- ggplot(snowmelt, aes(x = ros_flux_frac)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60, fill = "#377eb8", colour = "white", linewidth = 0.2) +
  geom_density(colour = "black", linewidth = 0.8) +
  geom_vline(xintercept = FLAG_THRESH, linetype = "dashed",
             colour = "#d95f02", linewidth = 0.8) +
  geom_vline(xintercept = ROS_THRESH, linetype = "dotted",
             colour = "firebrick", linewidth = 0.8) +
  annotate("text", x = FLAG_THRESH + 0.01, y = Inf,
           label = sprintf("Flag threshold\n(%.2f)", FLAG_THRESH),
           hjust = 0, vjust = 1.4, colour = "#d95f02", size = 3.2) +
  annotate("text", x = ROS_THRESH + 0.01, y = Inf,
           label = sprintf("ROS threshold\n(%.2f)", ROS_THRESH),
           hjust = 0, vjust = 1.4, colour = "firebrick", size = 3.2) +
  labs(
    title    = "ros_flux_frac within snowmelt-classified events",
    subtitle = sprintf(
      "n = %d snowmelt events | %.1f%% have ros_flux_frac ≥ %.2f | %.1f%% have ros_flux_frac ≥ %.2f",
      nrow(snowmelt),
      100 * mean(snowmelt$ros_flux_frac >= FLAG_THRESH), FLAG_THRESH,
      100 * mean(snowmelt$ros_flux_frac >= ROS_THRESH),  ROS_THRESH
    ),
    x = "ros_flux_frac (short window)",
    y = "Density"
  ) +
  theme_diag

ggsave(file.path(out_dir, "01_ros_frac_within_snowmelt.png"),
       p1, width = 9, height = 5, dpi = 150)

# ============================================================================
# 2. PEAK MAGNITUDE: FLAGGED VS CLEAN SNOWMELT
# ============================================================================
message("Plot 2: peak magnitude comparison...")

p2 <- ggplot(snowmelt %>% filter(!regulated),
             aes(x = peak_flow, fill = as.character(flagged),
                 colour = as.character(flagged))) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_x_log10(labels = scales::comma) +
  scale_fill_manual(values = flag_colors, labels = flag_labels, name = NULL) +
  scale_colour_manual(values = flag_colors, labels = flag_labels, name = NULL) +
  labs(
    title    = "Peak flow distribution: flagged vs clean snowmelt (unregulated)",
    subtitle = "Log10 scale. Flagged = snowmelt events with elevated ros_flux_frac.",
    x        = "Peak flow (cfs, log10 scale)",
    y        = "Density"
  ) +
  theme_diag +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "02_peak_magnitude_comparison.png"),
       p2, width = 9, height = 5, dpi = 150)

# ============================================================================
# 3. SEASONAL DISTRIBUTION: FLAGGED VS CLEAN
# ============================================================================
message("Plot 3: seasonal comparison...")

wy_month_breaks <- c(1, 32, 62, 93, 124, 152, 182, 213, 244, 274, 305, 335)
wy_month_labels <- c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep")

p3 <- ggplot(snowmelt, aes(x = wy_doy, fill = as.character(flagged),
                            colour = as.character(flagged))) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_x_continuous(breaks = wy_month_breaks, labels = wy_month_labels) +
  scale_fill_manual(values = flag_colors, labels = flag_labels, name = NULL) +
  scale_colour_manual(values = flag_colors, labels = flag_labels, name = NULL) +
  labs(
    title    = "Seasonal distribution: flagged vs clean snowmelt",
    subtitle = "Do flagged events occur earlier in the season (warm-season rain on residual snowpack)?",
    x        = "Month (water year)",
    y        = "Density"
  ) +
  theme_diag +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "03_seasonal_comparison.png"),
       p3, width = 9, height = 5, dpi = 150)

# ============================================================================
# 4. RECLASSIFICATION UNDER REORDERED DECISION TREE
# ============================================================================
# Current order: rainfall -> snowmelt -> ROS -> mixed
# Alternative:   rainfall -> ROS -> snowmelt -> mixed
# (i.e., check co-occurrence BEFORE committing to snowmelt)
# This reordering would pull high-ros snowmelt events into ROS instead.
message("Plot 4: reclassification under reordered tree...")

reclassified <- peaks %>%
  filter(!is.na(rf_short), !is.na(rf_long), !is.na(ros_flux_frac)) %>%
  mutate(
    flood_type_reordered = case_when(
      flood_type_raw == "unclassified"        ~ "unclassified",
      rf_short >= RAIN_EM_THRESH              ~ "rainfall",
      ros_flux_frac > ROS_THRESH              ~ "ROS",        # moved up
      rf_long  <= SNOW_EM_THRESH              ~ "snowmelt",   # moved down
      TRUE                                    ~ "mixed"
    )
  )

# Compare original vs reordered for events currently in snowmelt
transition_table <- reclassified %>%
  filter(flood_type_raw == "snowmelt") %>%
  count(flood_type_raw, flood_type_reordered) %>%
  mutate(pct = 100 * n / sum(n))

message("\nReclassification of snowmelt events under reordered tree:")
print(transition_table)

# Also look at overall population shifts
overall_transition <- reclassified %>%
  filter(flood_type_raw != "unclassified") %>%
  count(original = flood_type_raw, reordered = flood_type_reordered) %>%
  group_by(original) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

p4a <- ggplot(transition_table,
              aes(x = flood_type_reordered, y = n, fill = flood_type_reordered)) +
  geom_col() +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c(snowmelt = "#377eb8", ROS = "#4daf4a",
                                mixed = "#ff7f00", rainfall = "#e41a1c"),
                    guide = "none") +
  labs(
    title    = "Where do snowmelt events go under a reordered decision tree?",
    subtitle = "Reordered: rainfall → ROS → snowmelt → mixed\n(ros_flux_frac evaluated before rf_long)",
    x        = "Reclassified type",
    y        = "Number of events"
  ) +
  theme_diag

p4b <- ggplot(overall_transition %>% filter(original != reordered),
              aes(x = original, y = n, fill = reordered)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c(snowmelt = "#377eb8", ROS = "#4daf4a",
                                mixed = "#ff7f00", rainfall = "#e41a1c"),
                    name = "Reclassified as") +
  labs(
    title    = "Population-level type shifts from reordering",
    subtitle = "Only events that change classification are shown",
    x        = "Original type",
    y        = "Events reclassified"
  ) +
  theme_diag

ggsave(file.path(out_dir, "04_reclassification_counts.png"),
       p4a / p4b, width = 9, height = 9, dpi = 150)

# ============================================================================
# 5 & 6. ACCUMULATION CURVES: FLAGGED VS CLEAN SNOWMELT
# ============================================================================
message("Plotting accumulation curves...")

# Sample peaks to plot — stratify by flagged/clean, prefer unregulated
set.seed(42)
flagged_sample <- snowmelt %>%
  filter(flagged, !regulated) %>%
  slice_sample(n = min(N_CURVE_EXAMPLES, sum(snowmelt$flagged & !snowmelt$regulated)))

clean_sample <- snowmelt %>%
  filter(!flagged, !regulated) %>%
  slice_sample(n = min(N_CURVE_EXAMPLES, sum(!snowmelt$flagged & !snowmelt$regulated)))

plot_accum_curves <- function(sample_peaks, accum_df, title_str, subtitle_str,
                               short_window_col = "short_window",
                               long_window_col  = "long_window") {
  # Pull accumulation curves for the sampled peaks, up to long window
  curves <- accum_df %>%
    inner_join(sample_peaks %>% select(site_no, peak_id,
                                       short_w = all_of(short_window_col),
                                       long_w  = all_of(long_window_col),
                                       ros_flux_frac, peak_flow),
               by = c("site_no", "peak_id")) %>%
    filter(lookback_day <= long_w) %>%
    mutate(
      panel_label = sprintf("%s\n%s\nros_ff=%.2f",
                            site_no, peak_id, ros_flux_frac)
    )

  # Pivot to long for plotting both rain and melt
  curves_long <- curves %>%
    select(panel_label, lookback_day, short_w,
           Rain = rain_mm_day, Melt = melt_mm_day) %>%
    pivot_longer(c(Rain, Melt), names_to = "flux_type", values_to = "mm_day")

  ggplot(curves_long, aes(x = -lookback_day, y = mm_day,
                           colour = flux_type, fill = flux_type)) +
    geom_area(alpha = 0.3, position = "identity") +
    geom_line(linewidth = 0.5) +
    geom_vline(aes(xintercept = -short_w),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    scale_colour_manual(values = c(Rain = "#e41a1c", Melt = "#377eb8"),
                        name = NULL) +
    scale_fill_manual(values   = c(Rain = "#e41a1c", Melt = "#377eb8"),
                      name = NULL) +
    facet_wrap(~panel_label, ncol = 4, scales = "free_y") +
    labs(
      title    = title_str,
      subtitle = subtitle_str,
      x        = "Days before peak (0 = peak date)",
      y        = "Flux (mm/day)"
    ) +
    theme_diag +
    theme(
      legend.position = "bottom",
      axis.text       = element_text(size = 7),
      strip.text      = element_text(size = 6.5)
    )
}

p5 <- plot_accum_curves(
  flagged_sample, accum,
  title_str    = "Accumulation curves: flagged snowmelt events (high ros_flux_frac)",
  subtitle_str = sprintf(
    "Snowmelt-classified events with ros_flux_frac ≥ %.2f (n = %d shown)\nDashed vertical = short window boundary",
    FLAG_THRESH, nrow(flagged_sample)
  )
)

p6 <- plot_accum_curves(
  clean_sample, accum,
  title_str    = "Accumulation curves: clean snowmelt events (low ros_flux_frac)",
  subtitle_str = sprintf(
    "Snowmelt-classified events with ros_flux_frac < %.2f (n = %d shown)\nDashed vertical = short window boundary",
    FLAG_THRESH, nrow(clean_sample)
  )
)

ggsave(file.path(out_dir, "05_accum_curves_flagged.png"),
       p5, width = 14, height = 10, dpi = 150)

ggsave(file.path(out_dir, "06_accum_curves_clean.png"),
       p6, width = 14, height = 10, dpi = 150)

# ============================================================================
# CONSOLE SUMMARY
# ============================================================================
summary_lines <- c(
  "=== SNOWMELT HIGH ros_flux_frac DIAGNOSTIC SUMMARY ===",
  sprintf("Run date: %s", Sys.Date()),
  "",
  "--- Snowmelt event counts ---",
  sprintf("Total snowmelt (flood_type_raw): %d", nrow(snowmelt)),
  sprintf("Unregulated snowmelt           : %d", sum(!snowmelt$regulated)),
  "",
  "--- ros_flux_frac exceedance within snowmelt ---",
  capture.output(print(thresh_summary, n = Inf)),
  "",
  "--- Reclassification under reordered tree (snowmelt events only) ---",
  capture.output(print(transition_table, n = Inf)),
  "",
  "--- Overall population shifts from reordering ---",
  capture.output(print(overall_transition %>% filter(original != reordered), n = Inf)),
  "",
  sprintf("Outputs written to: %s", out_dir)
)

writeLines(summary_lines,
           file.path(out_dir, "diag_snowmelt_high_ros.txt"))
message(paste(summary_lines, collapse = "\n"))
message("\nDone.")