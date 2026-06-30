#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Rainfall Threshold Sensitivity
# ============================================================================
# Assesses sensitivity of the rainfall end-member classification to the
# rain_frac threshold (currently 0.95 in 6b). The analysis proceeds in
# two stages:
#
#   1. GEOMETRY — Examines the rain_frac distribution in the candidate
#      threshold zone (THRESH_MIN to 1.0) to determine whether the current
#      0.95 value sits in a natural density gap or cuts through a continuous
#      distribution. This is the primary diagnostic: if there is no gap,
#      the threshold is essentially arbitrary and downstream sensitivity
#      matters more.
#
#   2. GPD SENSITIVITY — For each candidate threshold, re-classifies peaks,
#      fits GPD shape parameters (xi) to the rainfall sample at each
#      unregulated site meeting MIN_TYPE_SAMPLE, and summarises the
#      distribution of xi across sites. Stable xi across the threshold range
#      supports the current choice; a shift at a specific threshold value
#      suggests a real process boundary.
#
# The ROS and snowmelt thresholds are held fixed at their 6b values
# throughout — only the rainfall threshold is varied.
#
# INPUTS:
#   - _data/pot_peaks_classified_adaptive.rds
#
# OUTPUTS (in _diagnostics/rainfall_threshold_sensitivity/):
#   - thresh_rain_frac_density.png    : rain_frac density in threshold zone
#   - thresh_rain_frac_gap.png        : density gap test across threshold range
#   - thresh_n_peaks.png              : rainfall peak count vs threshold
#   - thresh_melt_contamination.png   : absolute melt flux in 0.80-0.95 shoulder
#   - thresh_gpd_xi.png               : GPD shape parameter vs threshold
#   - thresh_gpd_xi_by_site.png       : per-site xi trajectories vs threshold
#   - thresh_summary.txt              : written summary
# ============================================================================

library(tidyverse)
library(lmomco)
library(dplyr)

classified_file <- "_data/pot_peaks_classified_adaptive.rds"
outdir          <- "_diagnostics/rainfall_threshold_sensitivity"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Fixed 6b thresholds — only RAIN_EM_THRESH is varied in this diagnostic
SNOW_EM_THRESH     <- 0.05
ROS_COOCCUR_THRESH <- 0.50
CURRENT_THRESH     <- 0.95   # current 6b value; highlighted in all plots

# Threshold scan range and step
THRESH_MIN  <- 0.70
THRESH_MAX  <- 0.99
THRESH_STEP <- 0.01
thresholds  <- seq(THRESH_MIN, THRESH_MAX, by = THRESH_STEP)

# Minimum peaks per site for GPD fitting
MIN_TYPE_SAMPLE <- 20

type_colors <- c(
  rainfall     = "#e41a1c",
  snowmelt     = "#377eb8",
  ROS          = "#4daf4a",
  mixed        = "#ff7f00",
  unclassified = "#cccccc"
)

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading data...")
dat <- readRDS(classified_file)

# Unregulated, data-sufficient peaks with complete features
peaks <- dat %>%
  filter(
    !regulated,
    data_sufficient,
    !is.na(rain_frac),
    !is.na(ros_flux_frac)
  )

message(sprintf("Peaks in scope: %d across %d sites\n",
                nrow(peaks), n_distinct(peaks$site_no)))

# ============================================================================
# STAGE 1: GEOMETRY
# ============================================================================
message("Stage 1: Geometry — rain_frac distribution in threshold zone...")

# --- Plot 1a: Full rain_frac density with threshold zone highlighted ---------
p1a <- ggplot(peaks, aes(x = rain_frac)) +
  annotate("rect", xmin = THRESH_MIN, xmax = 1,
           ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.6) +
  geom_density(fill = "#e41a1c", alpha = 0.4, colour = "#e41a1c",
               linewidth = 0.7) +
  geom_vline(xintercept = CURRENT_THRESH,
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  annotate("text", x = CURRENT_THRESH - 0.01,
           y = Inf, label = sprintf("Current threshold = %.2f", CURRENT_THRESH),
           hjust = 1, vjust = 1.5, size = 3, colour = "grey30") +
  annotate("text", x = (THRESH_MIN + 1) / 2, y = Inf,
           label = "Threshold scan zone", hjust = 0.5, vjust = 1.5,
           size = 3, colour = "grey50") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(
    title    = "rain_frac distribution — all unregulated data-sufficient peaks",
    subtitle = "Shaded zone = candidate threshold range; dashed = current 6b threshold",
    x        = "rain_frac",
    y        = "Density"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(outdir, "thresh_rain_frac_density.png"), p1a,
       width = 8, height = 5, dpi = 150)
message("  Saved: thresh_rain_frac_density.png")

# --- Plot 1b: Zoomed density in threshold zone with gap test ----------------
# Compute KDE over threshold zone to look for local minima (gaps)
rain_zone <- peaks %>% filter(rain_frac >= THRESH_MIN)

kde <- density(rain_zone$rain_frac, from = THRESH_MIN, to = 1.0, n = 512)
kde_df <- tibble(x = kde$x, y = kde$y)

# Find local minima in the KDE (potential gap locations)
kde_df <- kde_df %>%
  mutate(
    is_local_min = y < lag(y, default = Inf) & y < lead(y, default = Inf)
  )

local_mins <- kde_df %>% filter(is_local_min) %>% arrange(y)

p1b <- ggplot(kde_df, aes(x = x, y = y)) +
  geom_area(fill = "#e41a1c", alpha = 0.3) +
  geom_line(colour = "#e41a1c", linewidth = 0.8) +
  geom_vline(xintercept = CURRENT_THRESH,
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  {if (nrow(local_mins) > 0)
    geom_vline(data = local_mins,
               aes(xintercept = x),
               linetype = "dotted", colour = "steelblue", linewidth = 0.6)
  } +
  {if (nrow(local_mins) > 0)
    annotate("text", x = local_mins$x, y = max(kde_df$y) * 0.95,
             label = sprintf("gap\n%.3f", local_mins$x),
             hjust = 0.5, size = 2.8, colour = "steelblue")
  } +
  annotate("text", x = CURRENT_THRESH - 0.002, y = max(kde_df$y) * 0.5,
           label = sprintf("%.2f", CURRENT_THRESH),
           hjust = 1, size = 3, colour = "grey30") +
  scale_x_continuous(breaks = seq(THRESH_MIN, 1.0, by = 0.05)) +
  labs(
    title    = sprintf("rain_frac KDE — threshold zone (%.2f to 1.0)", THRESH_MIN),
    subtitle = "Blue dotted lines = local density minima (potential natural gaps)",
    x        = "rain_frac",
    y        = "Density"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(outdir, "thresh_rain_frac_gap.png"), p1b,
       width = 8, height = 5, dpi = 150)
message("  Saved: thresh_rain_frac_gap.png")

# --- Plot 1c: Peak count entering/leaving rainfall type at each threshold ---
count_df <- map_dfr(thresholds, function(t) {
  tibble(
    threshold    = t,
    n_rainfall   = sum(peaks$rain_frac >= t),
    n_nonrainfall = sum(peaks$rain_frac < t)
  )
}) %>%
  pivot_longer(-threshold, names_to = "type", values_to = "n")

p1c <- ggplot(count_df %>% filter(type == "n_rainfall"),
              aes(x = threshold, y = n)) +
  geom_line(colour = "#e41a1c", linewidth = 0.9) +
  geom_point(colour = "#e41a1c", size = 1.5) +
  geom_vline(xintercept = CURRENT_THRESH,
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  scale_x_continuous(breaks = seq(THRESH_MIN, 1.0, by = 0.05)) +
  labs(
    title    = "Rainfall peak count vs. threshold",
    subtitle = "Steep drop = threshold crossing a dense region; flat = sparse region",
    x        = "rain_frac threshold",
    y        = "Number of peaks classified as rainfall"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(outdir, "thresh_n_peaks.png"), p1c,
       width = 8, height = 5, dpi = 150)
message("  Saved: thresh_n_peaks.png")

# --- Plot 1d: Melt contamination in the shoulder zone -----------------------
# For peaks in the 0.80–0.95 zone, examine absolute melt flux to assess
# whether the non-trivial rain fraction reflects meaningful melt volumes
# or negligible absolute melt that wouldn't materially affect FFA.
SHOULDER_LO <- 0.80
SHOULDER_HI <- CURRENT_THRESH

shoulder <- peaks %>%
  filter(rain_frac >= SHOULDER_LO, rain_frac < SHOULDER_HI) %>%
  mutate(rain_frac_bin = cut(rain_frac,
                             breaks = seq(SHOULDER_LO, SHOULDER_HI, by = 0.05),
                             include.lowest = TRUE, right = FALSE))

# Absolute melt flux distribution within shoulder zone
p1d_abs <- ggplot(shoulder, aes(x = melt_flux_short)) +
  geom_histogram(bins = 40, fill = "#377eb8", colour = "white",
                 linewidth = 0.2) +
  geom_vline(xintercept = 1.0, linetype = "dashed",
             colour = "grey30", linewidth = 0.7) +
  annotate("text", x = 1.05, y = Inf,
           label = "1 mm", hjust = 0, vjust = 1.5,
           size = 3, colour = "grey30") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 5, 10, 25, 50, 100),
                     labels = c("0", "1", "5", "10", "25", "50", "100")) +
  labs(
    title    = sprintf("Absolute melt flux for shoulder zone peaks (%.2f ≤ rain_frac < %.2f)",
                       SHOULDER_LO, SHOULDER_HI),
    subtitle = sprintf("n = %d peaks; x-axis log1p scale", nrow(shoulder)),
    x        = "melt_flux_short (mm, log1p scale)",
    y        = "Count"
  ) +
  theme_bw(base_size = 11)

# Melt flux by rain_frac bin within shoulder
p1d_bin <- ggplot(shoulder %>% filter(!is.na(rain_frac_bin)),
                  aes(x = rain_frac_bin, y = melt_flux_short)) +
  geom_boxplot(fill = "#377eb8", alpha = 0.5, outlier.size = 0.8,
               outlier.alpha = 0.4, linewidth = 0.4) +
  geom_hline(yintercept = 1.0, linetype = "dashed",
             colour = "grey30", linewidth = 0.7) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 5, 10, 25, 50, 100),
                     labels = c("0", "1", "5", "10", "25", "50", "100")) +
  labs(
    title    = "Melt flux by rain_frac bin — shoulder zone",
    subtitle = "y-axis log1p scale; dashed = 1 mm reference",
    x        = "rain_frac bin",
    y        = "melt_flux_short (mm)"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

library(patchwork)
p1d <- p1d_abs / p1d_bin +
  plot_annotation(
    title   = "Melt contamination in the rainfall threshold shoulder zone",
    caption = sprintf(
      "If melt_flux_short is negligible (< 1 mm) throughout, lower threshold is defensible."
    )
  )

ggsave(file.path(outdir, "thresh_melt_contamination.png"), p1d,
       width = 8, height = 9, dpi = 150)
message("  Saved: thresh_melt_contamination.png")

# Console summary of shoulder zone melt flux
message(sprintf("\n  Shoulder zone (%.2f to %.2f): %d peaks",
                SHOULDER_LO, SHOULDER_HI, nrow(shoulder)))
message("  melt_flux_short summary (mm):")
print(summary(shoulder$melt_flux_short))
message(sprintf("  Fraction with melt_flux_short < 1 mm  : %.1f%%",
                100 * mean(shoulder$melt_flux_short < 1.0, na.rm = TRUE)))
message(sprintf("  Fraction with melt_flux_short < 5 mm  : %.1f%%",
                100 * mean(shoulder$melt_flux_short < 5.0, na.rm = TRUE)))
message(sprintf("  Fraction with melt_flux_short < 10 mm : %.1f%%",
                100 * mean(shoulder$melt_flux_short < 10.0, na.rm = TRUE)))
message("\nStage 2: GPD shape parameter sensitivity across threshold range...")

# GPD fit via L-moments using lmomco
# Returns xi (shape), alpha (scale), or NA if fit fails or n < min_n
fit_gpd_xi <- function(x, u, min_n = MIN_TYPE_SAMPLE) {
  # x: vector of peak flows; u: GPD location (pot_threshold_cfs)
  exceedances <- x - u
  exceedances <- exceedances[exceedances > 0]
  if (length(exceedances) < min_n) return(NA_real_)
  tryCatch({
    lmr <- lmoms(exceedances)
    fit <- pargpa(lmr)
    xi  <- fit$para["xi"]
    if (is.null(xi) || length(xi) == 0) NA_real_ else as.numeric(xi)
  }, error = function(e) NA_real_)
}

message(sprintf("  Scanning %d thresholds from %.2f to %.2f...",
                length(thresholds), THRESH_MIN, THRESH_MAX))

gpd_results <- map_dfr(thresholds, function(t) {
  # Re-classify under this threshold — ROS and snowmelt thresholds held fixed
  reclassified <- peaks %>%
    mutate(
      flood_type_t = case_when(
        rain_frac     >= t                  ~ "rainfall",
        ros_flux_frac >  ROS_COOCCUR_THRESH ~ "ROS",
        rain_frac     <= SNOW_EM_THRESH     ~ "snowmelt",
        TRUE                                ~ "mixed"
      )
    )

  # Per-site GPD fits for rainfall type
  site_fits <- reclassified %>%
    filter(flood_type_t == "rainfall") %>%
    group_by(site_no) %>%
    summarise(
      xi = fit_gpd_xi(peak_flow, first(pot_threshold_cfs)),
      n  = sum(flood_type_t == "rainfall"),
      .groups = "drop"
    ) %>%
    filter(!is.na(xi))

  if (nrow(site_fits) == 0) return(NULL)

  tibble(
    threshold    = t,
    n_sites      = nrow(site_fits),
    xi_median    = median(site_fits$xi, na.rm = TRUE),
    xi_q25       = quantile(site_fits$xi, 0.25, na.rm = TRUE),
    xi_q75       = quantile(site_fits$xi, 0.75, na.rm = TRUE),
    xi_mean      = mean(site_fits$xi, na.rm = TRUE),
    xi_sd        = sd(site_fits$xi, na.rm = TRUE)
  )
}, .progress = TRUE)

message(sprintf("  GPD fits complete: %d threshold values", nrow(gpd_results)))

# --- Plot 2a: GPD xi (median + IQR) vs threshold ----------------------------
p2a <- ggplot(gpd_results, aes(x = threshold)) +
  geom_ribbon(aes(ymin = xi_q25, ymax = xi_q75),
              fill = "#e41a1c", alpha = 0.2) +
  geom_line(aes(y = xi_median), colour = "#e41a1c", linewidth = 0.9) +
  geom_point(aes(y = xi_median), colour = "#e41a1c", size = 1.5) +
  geom_vline(xintercept = CURRENT_THRESH,
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  annotate("text", x = CURRENT_THRESH - 0.002,
           y = max(gpd_results$xi_q75, na.rm = TRUE),
           label = sprintf("%.2f", CURRENT_THRESH),
           hjust = 1, size = 3, colour = "grey30") +
  scale_x_continuous(breaks = seq(THRESH_MIN, 1.0, by = 0.05)) +
  labs(
    title    = "GPD shape parameter (xi) for rainfall type vs. threshold",
    subtitle = "Median ± IQR across sites; stable xi = threshold-insensitive",
    x        = "rain_frac threshold",
    y        = "GPD shape parameter (xi)",
    caption  = sprintf("Sites with >= %d rainfall peaks included at each threshold",
                       MIN_TYPE_SAMPLE)
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(outdir, "thresh_gpd_xi.png"), p2a,
       width = 8, height = 5, dpi = 150)
message("  Saved: thresh_gpd_xi.png")

# --- Plot 2b: Per-site xi trajectories ---------------------------------------
# Compute per-site xi at each threshold for spaghetti plot
message("  Computing per-site xi trajectories...")

site_xi_long <- map_dfr(thresholds, function(t) {
  reclassified <- peaks %>%
    mutate(
      flood_type_t = case_when(
        rain_frac     >= t                  ~ "rainfall",
        ros_flux_frac >  ROS_COOCCUR_THRESH ~ "ROS",
        rain_frac     <= SNOW_EM_THRESH     ~ "snowmelt",
        TRUE                                ~ "mixed"
      )
    )

  reclassified %>%
    filter(flood_type_t == "rainfall") %>%
    group_by(site_no) %>%
    summarise(
      xi = fit_gpd_xi(peak_flow, first(pot_threshold_cfs)),
      n  = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(xi)) %>%
    mutate(threshold = t)
}, .progress = TRUE)

p2b <- ggplot(site_xi_long, aes(x = threshold, y = xi, group = site_no)) +
  geom_line(alpha = 0.2, colour = "#e41a1c", linewidth = 0.4) +
  geom_vline(xintercept = CURRENT_THRESH,
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  scale_x_continuous(breaks = seq(THRESH_MIN, 1.0, by = 0.05)) +
  labs(
    title    = "Per-site GPD shape parameter (xi) trajectories vs. threshold",
    subtitle = "Each line = one site; stable trajectories = robust threshold choice",
    x        = "rain_frac threshold",
    y        = "GPD shape parameter (xi)"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(outdir, "thresh_gpd_xi_by_site.png"), p2b,
       width = 8, height = 5, dpi = 150)
message("  Saved: thresh_gpd_xi_by_site.png")

# ============================================================================
# SUMMARY
# ============================================================================
# Identify candidate gap locations from Stage 1
gap_summary <- if (nrow(local_mins) > 0) {
  paste(sprintf("  %.3f (density = %.4f)", local_mins$x, local_mins$y),
        collapse = "\n")
} else {
  "  None detected — distribution appears continuous in threshold zone"
}

# Assess stability: coefficient of variation of xi_median across thresholds
xi_cv <- sd(gpd_results$xi_median, na.rm = TRUE) /
         abs(mean(gpd_results$xi_median, na.rm = TRUE))

summary_lines <- c(
  "=== Rainfall Threshold Sensitivity Summary ===",
  sprintf("Peaks in scope    : %d (%d sites, unregulated, data-sufficient)",
          nrow(peaks), n_distinct(peaks$site_no)),
  sprintf("Threshold range   : %.2f to %.2f (step %.2f)",
          THRESH_MIN, THRESH_MAX, THRESH_STEP),
  sprintf("Current threshold : %.2f", CURRENT_THRESH),
  "",
  "--- Stage 1: Geometry ---",
  sprintf("Peaks with rain_frac >= %.2f (current): %d",
          CURRENT_THRESH, sum(peaks$rain_frac >= CURRENT_THRESH)),
  sprintf("Peaks with rain_frac >= %.2f (min scan): %d",
          THRESH_MIN, sum(peaks$rain_frac >= THRESH_MIN)),
  "Local density minima in threshold zone:",
  gap_summary,
  "",
  sprintf("--- Shoulder zone melt contamination (%.2f to %.2f) ---",
          SHOULDER_LO, SHOULDER_HI),
  sprintf("Peaks in shoulder zone : %d", nrow(shoulder)),
  capture.output(print(summary(shoulder$melt_flux_short))),
  sprintf("Fraction with melt_flux_short < 1 mm  : %.1f%%",
          100 * mean(shoulder$melt_flux_short < 1.0, na.rm = TRUE)),
  sprintf("Fraction with melt_flux_short < 5 mm  : %.1f%%",
          100 * mean(shoulder$melt_flux_short < 5.0, na.rm = TRUE)),
  sprintf("Fraction with melt_flux_short < 10 mm : %.1f%%",
          100 * mean(shoulder$melt_flux_short < 10.0, na.rm = TRUE)),
  "",
  "--- Stage 2: GPD Sensitivity ---",
  capture.output(print(gpd_results %>%
    select(threshold, n_sites, xi_median, xi_q25, xi_q75) %>%
    mutate(across(where(is.numeric), ~ round(., 3))))),
  "",
  sprintf("xi_median CV across threshold range: %.3f", xi_cv),
  sprintf("  (CV < 0.10 suggests low sensitivity; CV > 0.25 suggests meaningful variation)")
)

writeLines(summary_lines, file.path(outdir, "thresh_summary.txt"))
message("  Saved: thresh_summary.txt")

cat(paste(summary_lines, collapse = "\n"), "\n")

message("\nDiagnostic complete. Outputs in: ", outdir)