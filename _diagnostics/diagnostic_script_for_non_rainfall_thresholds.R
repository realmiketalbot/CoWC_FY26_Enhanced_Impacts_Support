#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Non-Rainfall Threshold Sensitivity
# ============================================================================
# Evaluates sensitivity of the ROS and snowmelt end-member thresholds within
# the non-rainfall population (rain_frac < RAIN_EM_THRESH). Two boundaries
# are assessed:
#
#   1. ROS threshold      : ros_flux_frac > ROS_COOCCUR_THRESH (currently 0.50)
#      Peaks above this are classified ROS; below are candidates for snowmelt
#      or mixed. Scan range: 0.20 to 0.80.
#
#   2. Snowmelt threshold : rain_frac <= SNOW_EM_THRESH (currently 0.05)
#      Peaks at or below this are classified snowmelt (after ROS check).
#      Scan range: 0.01 to 0.20.
#
# For each boundary, the analysis proceeds in two stages:
#
#   GEOMETRY — Distribution of the relevant feature in the candidate
#     threshold zone, peak counts vs threshold, and absolute flux
#     contamination in the shoulder zone (rain flux for snowmelt boundary,
#     ros_flux for ROS boundary).
#
#   GPD SENSITIVITY — Per-site GPD shape parameter (xi) fitted to the
#     type sample at each threshold, summarised as median ± IQR across
#     sites and as per-site trajectories. Stable xi = threshold-insensitive.
#
# The rainfall threshold is held fixed at RAIN_EM_THRESH throughout.
# For each ROS threshold scan step, the snowmelt threshold is held at
# its current value, and vice versa.
#
# INPUTS:
#   - _data/pot_peaks_classified_adaptive.rds
#
# OUTPUTS (in _diagnostics/nonrainfall_threshold_sensitivity/):
#   ROS boundary:
#   - ros_flux_density.png           : ros_flux_frac distribution
#   - ros_flux_gap.png               : KDE in scan zone with gap detection
#   - ros_n_peaks.png                : ROS/snowmelt peak count vs threshold
#   - ros_rain_contamination.png     : rain flux in ROS shoulder zone
#   - ros_gpd_xi.png                 : GPD xi (ROS type) vs threshold
#   - ros_gpd_xi_by_site.png         : per-site xi trajectories (ROS)
#   - ros_gpd_xi_snowmelt.png        : GPD xi (snowmelt type) vs ROS threshold
#   Snowmelt boundary:
#   - snow_rain_frac_density.png     : rain_frac distribution (low end)
#   - snow_rain_frac_gap.png         : KDE in scan zone with gap detection
#   - snow_n_peaks.png               : snowmelt peak count vs threshold
#   - snow_melt_contamination.png    : melt flux in snowmelt shoulder zone
#   - snow_gpd_xi.png                : GPD xi (snowmelt type) vs threshold
#   - snow_gpd_xi_by_site.png        : per-site xi trajectories (snowmelt)
#   Summary:
#   - nonrainfall_thresh_summary.txt : written summary
# ============================================================================

library(tidyverse)
library(lmomco)
library(dplyr)
library(patchwork)

classified_file <- "_data/pot_peaks_classified_adaptive.rds"
outdir          <- "_diagnostics/nonrainfall_threshold_sensitivity"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Fixed thresholds
RAIN_EM_THRESH     <- 0.95   # held fixed throughout
CURRENT_ROS_THRESH <- 0.50   # current 6b ROS threshold
CURRENT_SNOW_THRESH <- 0.05  # current 6b snowmelt threshold

# ROS threshold scan
ROS_THRESH_MIN  <- 0.20
ROS_THRESH_MAX  <- 0.80
ROS_THRESH_STEP <- 0.01
ros_thresholds  <- seq(ROS_THRESH_MIN, ROS_THRESH_MAX, by = ROS_THRESH_STEP)

# Snowmelt threshold scan
SNOW_THRESH_MIN  <- 0.01
SNOW_THRESH_MAX  <- 0.20
SNOW_THRESH_STEP <- 0.01
snow_thresholds  <- seq(SNOW_THRESH_MIN, SNOW_THRESH_MAX, by = SNOW_THRESH_STEP)

# Shoulder zones for contamination analysis
ROS_SHOULDER_LO  <- 0.35
ROS_SHOULDER_HI  <- 0.65
SNOW_SHOULDER_LO <- 0.03
SNOW_SHOULDER_HI <- 0.12

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

# Non-rainfall population: unregulated, data-sufficient, rain_frac < 0.95
peaks <- dat %>%
  filter(
    !regulated,
    data_sufficient,
    rain_frac < RAIN_EM_THRESH,
    !is.na(rain_frac),
    !is.na(ros_flux_frac),
    !is.na(melt_flux_short),
    !is.na(rain_flux_short)
  )

message(sprintf("Non-rainfall peaks in scope: %d across %d sites\n",
                nrow(peaks), n_distinct(peaks$site_no)))

# GPD fit helper
fit_gpd_xi <- function(x, u, min_n = MIN_TYPE_SAMPLE) {
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

# KDE gap helper — returns local minima sorted by density
find_gaps <- function(x, from, to, n = 512) {
  kde <- density(x, from = from, to = to, n = n)
  df  <- tibble(x = kde$x, y = kde$y) %>%
    mutate(
      is_local_min = y < lag(y, default = Inf) & y < lead(y, default = Inf)
    )
  df %>% filter(is_local_min) %>% arrange(y)
}

# ============================================================================
# PART 1: ROS THRESHOLD SENSITIVITY
# ============================================================================
message("========================================")
message("Part 1: ROS threshold sensitivity")
message("========================================\n")

# --- Plot R1: ros_flux_frac full distribution --------------------------------
p_r1 <- ggplot(peaks, aes(x = ros_flux_frac)) +
  annotate("rect", xmin = ROS_THRESH_MIN, xmax = ROS_THRESH_MAX,
           ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.6) +
  geom_density(fill = "#4daf4a", alpha = 0.4, colour = "#4daf4a",
               linewidth = 0.7) +
  geom_vline(xintercept = CURRENT_ROS_THRESH,
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  annotate("text", x = CURRENT_ROS_THRESH + 0.01, y = Inf,
           label = sprintf("Current = %.2f", CURRENT_ROS_THRESH),
           hjust = 0, vjust = 1.5, size = 3, colour = "grey30") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(
    title    = "ros_flux_frac distribution — non-rainfall population",
    subtitle = "Shaded zone = candidate threshold range; dashed = current 6b threshold",
    x        = "ros_flux_frac",
    y        = "Density"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(outdir, "ros_flux_density.png"), p_r1,
       width = 8, height = 5, dpi = 150)
message("  Saved: ros_flux_density.png")

# --- Plot R2: KDE in scan zone with gap detection ----------------------------
ros_zone  <- peaks %>% filter(ros_flux_frac >= ROS_THRESH_MIN,
                               ros_flux_frac <= ROS_THRESH_MAX)
ros_gaps  <- find_gaps(ros_zone$ros_flux_frac,
                       from = ROS_THRESH_MIN, to = ROS_THRESH_MAX)
ros_kde   <- density(ros_zone$ros_flux_frac,
                     from = ROS_THRESH_MIN, to = ROS_THRESH_MAX, n = 512)
ros_kde_df <- tibble(x = ros_kde$x, y = ros_kde$y)

p_r2 <- ggplot(ros_kde_df, aes(x = x, y = y)) +
  geom_area(fill = "#4daf4a", alpha = 0.3) +
  geom_line(colour = "#4daf4a", linewidth = 0.8) +
  geom_vline(xintercept = CURRENT_ROS_THRESH,
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  {if (nrow(ros_gaps) > 0)
    geom_vline(data = ros_gaps, aes(xintercept = x),
               linetype = "dotted", colour = "steelblue", linewidth = 0.6)} +
  {if (nrow(ros_gaps) > 0)
    annotate("text", x = ros_gaps$x, y = max(ros_kde_df$y) * 0.95,
             label = sprintf("gap\n%.3f", ros_gaps$x),
             hjust = 0.5, size = 2.8, colour = "steelblue")} +
  annotate("text", x = CURRENT_ROS_THRESH + 0.01,
           y = max(ros_kde_df$y) * 0.5,
           label = sprintf("%.2f", CURRENT_ROS_THRESH),
           hjust = 0, size = 3, colour = "grey30") +
  scale_x_continuous(breaks = seq(ROS_THRESH_MIN, ROS_THRESH_MAX, by = 0.05)) +
  labs(
    title    = sprintf("ros_flux_frac KDE — scan zone (%.2f to %.2f)",
                       ROS_THRESH_MIN, ROS_THRESH_MAX),
    subtitle = "Blue dotted = local density minima; dashed = current threshold",
    x        = "ros_flux_frac",
    y        = "Density"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(outdir, "ros_flux_gap.png"), p_r2,
       width = 8, height = 5, dpi = 150)
message("  Saved: ros_flux_gap.png")

# --- Plot R3: Peak counts vs ROS threshold -----------------------------------
ros_count_df <- map_dfr(ros_thresholds, function(t) {
  tibble(
    threshold  = t,
    n_ROS      = sum(peaks$ros_flux_frac > t),
    n_nonROS   = sum(peaks$ros_flux_frac <= t)
  )
}) %>% pivot_longer(-threshold, names_to = "type", values_to = "n")

p_r3 <- ggplot(ros_count_df, aes(x = threshold, y = n, colour = type)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = CURRENT_ROS_THRESH,
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  scale_colour_manual(values = c(n_ROS = "#4daf4a", n_nonROS = "#377eb8"),
                      labels = c(n_ROS = "ROS", n_nonROS = "Snowmelt + Mixed")) +
  scale_x_continuous(breaks = seq(ROS_THRESH_MIN, ROS_THRESH_MAX, by = 0.05)) +
  labs(
    title    = "Peak counts vs. ROS threshold",
    subtitle = "Steep change = threshold crossing a dense region",
    x        = "ros_flux_frac threshold",
    y        = "Number of peaks",
    colour   = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "ros_n_peaks.png"), p_r3,
       width = 8, height = 5, dpi = 150)
message("  Saved: ros_n_peaks.png")

# --- Plot R4: Rain flux contamination in ROS shoulder zone ------------------
ros_shoulder <- peaks %>%
  filter(ros_flux_frac >= ROS_SHOULDER_LO,
         ros_flux_frac <= ROS_SHOULDER_HI) %>%
  mutate(
    ros_bin = cut(ros_flux_frac,
                  breaks = seq(ROS_SHOULDER_LO, ROS_SHOULDER_HI, by = 0.05),
                  include.lowest = TRUE, right = FALSE),
    side = if_else(ros_flux_frac < CURRENT_ROS_THRESH, "Below threshold\n(→ snowmelt/mixed)",
                   "Above threshold\n(→ ROS)")
  )

p_r4a <- ggplot(ros_shoulder, aes(x = rain_flux_short)) +
  geom_histogram(bins = 40, fill = "#4daf4a", colour = "white",
                 linewidth = 0.2) +
  geom_vline(xintercept = 1.0, linetype = "dashed",
             colour = "grey30", linewidth = 0.7) +
  annotate("text", x = 1.05, y = Inf,
           label = "1 mm", hjust = 0, vjust = 1.5,
           size = 3, colour = "grey30") +
  facet_wrap(~ side, ncol = 2) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 5, 10, 25, 50, 100, 250),
                     labels = c("0", "1", "5", "10", "25", "50", "100", "250")) +
  labs(
    title    = sprintf("Rain flux in ROS shoulder zone (%.2f ≤ ros_flux_frac ≤ %.2f)",
                       ROS_SHOULDER_LO, ROS_SHOULDER_HI),
    subtitle = sprintf("n = %d peaks; x-axis log1p scale", nrow(ros_shoulder)),
    x        = "rain_flux_short (mm, log1p scale)",
    y        = "Count"
  ) +
  theme_bw(base_size = 11)

p_r4b <- ggplot(ros_shoulder %>% filter(!is.na(ros_bin)),
                aes(x = ros_bin, y = rain_flux_short)) +
  geom_boxplot(fill = "#4daf4a", alpha = 0.5, outlier.size = 0.8,
               outlier.alpha = 0.4, linewidth = 0.4) +
  geom_vline(xintercept = which(levels(ros_shoulder$ros_bin) ==
               levels(ros_shoulder$ros_bin)[ceiling(length(levels(
                 ros_shoulder$ros_bin)) / 2)]) + 0.5,
             linetype = "dashed", colour = "grey30", linewidth = 0.6) +
  geom_hline(yintercept = 1.0, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 5, 10, 25, 50, 100, 250),
                     labels = c("0", "1", "5", "10", "25", "50", "100", "250")) +
  labs(
    title    = "Rain flux by ros_flux_frac bin — shoulder zone",
    subtitle = "y-axis log1p scale",
    x        = "ros_flux_frac bin",
    y        = "rain_flux_short (mm)"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_r4 <- p_r4a / p_r4b +
  plot_annotation(
    title   = "Rain flux contamination in ROS/snowmelt boundary zone",
    caption = "If rain flux is similar on both sides of threshold, boundary is arbitrary in flux space."
  )

ggsave(file.path(outdir, "ros_rain_contamination.png"), p_r4,
       width = 8, height = 9, dpi = 150)
message("  Saved: ros_rain_contamination.png")

# --- Plots R5/R6: GPD xi for ROS and snowmelt types vs ROS threshold --------
message("\n  Scanning GPD xi across ROS threshold range...")

ros_gpd_results <- map_dfr(ros_thresholds, function(t) {
  rc <- peaks %>%
    mutate(flood_type_t = case_when(
      ros_flux_frac >  t                   ~ "ROS",
      rain_frac     <= CURRENT_SNOW_THRESH ~ "snowmelt",
      TRUE                                 ~ "mixed"
    ))

  fit_type <- function(type_label) {
    rc %>%
      filter(flood_type_t == type_label) %>%
      group_by(site_no) %>%
      summarise(
        xi = fit_gpd_xi(peak_flow, first(pot_threshold_cfs)),
        .groups = "drop"
      ) %>%
      filter(!is.na(xi)) %>%
      summarise(
        n_sites   = n(),
        xi_median = median(xi),
        xi_q25    = quantile(xi, 0.25),
        xi_q75    = quantile(xi, 0.75)
      ) %>%
      mutate(type = type_label, threshold = t)
  }

  bind_rows(fit_type("ROS"), fit_type("snowmelt"))
}, .progress = TRUE)

make_xi_plot <- function(df, type_label, colour, title_prefix) {
  d <- df %>% filter(type == type_label)
  ggplot(d, aes(x = threshold)) +
    geom_ribbon(aes(ymin = xi_q25, ymax = xi_q75),
                fill = colour, alpha = 0.2) +
    geom_line(aes(y = xi_median), colour = colour, linewidth = 0.9) +
    geom_point(aes(y = xi_median), colour = colour, size = 1.5) +
    geom_vline(xintercept = CURRENT_ROS_THRESH,
               linetype = "dashed", colour = "grey30", linewidth = 0.7) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
    scale_x_continuous(breaks = seq(ROS_THRESH_MIN, ROS_THRESH_MAX, by = 0.05)) +
    labs(
      title    = sprintf("%s GPD xi vs. ros_flux_frac threshold", title_prefix),
      subtitle = "Median ± IQR across sites",
      x        = "ros_flux_frac threshold",
      y        = "GPD shape parameter (xi)",
      caption  = sprintf("Sites with >= %d peaks of this type", MIN_TYPE_SAMPLE)
    ) +
    theme_bw(base_size = 11)
}

p_r5 <- make_xi_plot(ros_gpd_results, "ROS",      "#4daf4a", "ROS type")
p_r6 <- make_xi_plot(ros_gpd_results, "snowmelt", "#377eb8", "Snowmelt type")

ggsave(file.path(outdir, "ros_gpd_xi.png"), p_r5,
       width = 8, height = 5, dpi = 150)
message("  Saved: ros_gpd_xi.png")

ggsave(file.path(outdir, "ros_gpd_xi_snowmelt.png"), p_r6,
       width = 8, height = 5, dpi = 150)
message("  Saved: ros_gpd_xi_snowmelt.png")

# Per-site spaghetti for ROS type
message("  Computing per-site ROS xi trajectories...")
ros_site_xi <- map_dfr(ros_thresholds, function(t) {
  peaks %>%
    mutate(flood_type_t = if_else(ros_flux_frac > t, "ROS", "other")) %>%
    filter(flood_type_t == "ROS") %>%
    group_by(site_no) %>%
    summarise(xi = fit_gpd_xi(peak_flow, first(pot_threshold_cfs)),
              .groups = "drop") %>%
    filter(!is.na(xi)) %>%
    mutate(threshold = t)
}, .progress = TRUE)

p_r7 <- ggplot(ros_site_xi, aes(x = threshold, y = xi, group = site_no)) +
  geom_line(alpha = 0.2, colour = "#4daf4a", linewidth = 0.4) +
  geom_vline(xintercept = CURRENT_ROS_THRESH,
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  scale_x_continuous(breaks = seq(ROS_THRESH_MIN, ROS_THRESH_MAX, by = 0.05)) +
  labs(
    title    = "Per-site GPD xi trajectories (ROS type) vs. ros_flux_frac threshold",
    subtitle = "Each line = one site",
    x        = "ros_flux_frac threshold",
    y        = "GPD shape parameter (xi)"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(outdir, "ros_gpd_xi_by_site.png"), p_r7,
       width = 8, height = 5, dpi = 150)
message("  Saved: ros_gpd_xi_by_site.png")

# ============================================================================
# PART 2: SNOWMELT THRESHOLD SENSITIVITY
# ============================================================================
message("\n========================================")
message("Part 2: Snowmelt threshold sensitivity")
message("========================================\n")

# --- Plot S1: rain_frac distribution (low end) --------------------------------
p_s1 <- ggplot(peaks, aes(x = rain_frac)) +
  annotate("rect", xmin = SNOW_THRESH_MIN, xmax = SNOW_THRESH_MAX,
           ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.6) +
  geom_density(fill = "#377eb8", alpha = 0.4, colour = "#377eb8",
               linewidth = 0.7) +
  geom_vline(xintercept = CURRENT_SNOW_THRESH,
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  annotate("text", x = CURRENT_SNOW_THRESH + 0.002, y = Inf,
           label = sprintf("Current = %.2f", CURRENT_SNOW_THRESH),
           hjust = 0, vjust = 1.5, size = 3, colour = "grey30") +
  scale_x_continuous(limits = c(0, 0.30),
                     breaks = seq(0, 0.30, by = 0.05)) +
  labs(
    title    = "rain_frac distribution (low end) — non-rainfall population",
    subtitle = "Shaded zone = snowmelt threshold scan range; x-axis truncated at 0.30",
    x        = "rain_frac",
    y        = "Density"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(outdir, "snow_rain_frac_density.png"), p_s1,
       width = 8, height = 5, dpi = 150)
message("  Saved: snow_rain_frac_density.png")

# --- Plot S2: KDE in snowmelt scan zone --------------------------------------
snow_zone <- peaks %>%
  filter(rain_frac >= SNOW_THRESH_MIN, rain_frac <= SNOW_THRESH_MAX)
snow_gaps <- find_gaps(snow_zone$rain_frac,
                       from = SNOW_THRESH_MIN, to = SNOW_THRESH_MAX)
snow_kde  <- density(snow_zone$rain_frac,
                     from = SNOW_THRESH_MIN, to = SNOW_THRESH_MAX, n = 512)
snow_kde_df <- tibble(x = snow_kde$x, y = snow_kde$y)

p_s2 <- ggplot(snow_kde_df, aes(x = x, y = y)) +
  geom_area(fill = "#377eb8", alpha = 0.3) +
  geom_line(colour = "#377eb8", linewidth = 0.8) +
  geom_vline(xintercept = CURRENT_SNOW_THRESH,
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  {if (nrow(snow_gaps) > 0)
    geom_vline(data = snow_gaps, aes(xintercept = x),
               linetype = "dotted", colour = "steelblue", linewidth = 0.6)} +
  {if (nrow(snow_gaps) > 0)
    annotate("text", x = snow_gaps$x, y = max(snow_kde_df$y) * 0.95,
             label = sprintf("gap\n%.3f", snow_gaps$x),
             hjust = 0.5, size = 2.8, colour = "steelblue")} +
  scale_x_continuous(breaks = seq(SNOW_THRESH_MIN, SNOW_THRESH_MAX, by = 0.02)) +
  labs(
    title    = sprintf("rain_frac KDE — snowmelt scan zone (%.2f to %.2f)",
                       SNOW_THRESH_MIN, SNOW_THRESH_MAX),
    subtitle = "Blue dotted = local density minima; dashed = current threshold",
    x        = "rain_frac",
    y        = "Density"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(outdir, "snow_rain_frac_gap.png"), p_s2,
       width = 8, height = 5, dpi = 150)
message("  Saved: snow_rain_frac_gap.png")

# --- Plot S3: Peak counts vs snowmelt threshold ------------------------------
snow_count_df <- map_dfr(snow_thresholds, function(t) {
  tibble(
    threshold   = t,
    n_snowmelt  = sum(peaks$rain_frac <= t &
                        peaks$ros_flux_frac <= CURRENT_ROS_THRESH),
    n_other     = sum(peaks$rain_frac > t |
                        peaks$ros_flux_frac > CURRENT_ROS_THRESH)
  )
}) %>% pivot_longer(-threshold, names_to = "type", values_to = "n")

p_s3 <- ggplot(snow_count_df %>% filter(type == "n_snowmelt"),
               aes(x = threshold, y = n)) +
  geom_line(colour = "#377eb8", linewidth = 0.9) +
  geom_point(colour = "#377eb8", size = 1.5) +
  geom_vline(xintercept = CURRENT_SNOW_THRESH,
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  scale_x_continuous(breaks = seq(SNOW_THRESH_MIN, SNOW_THRESH_MAX, by = 0.02)) +
  labs(
    title    = "Snowmelt peak count vs. rain_frac threshold",
    subtitle = "ROS threshold held fixed at 0.50",
    x        = "rain_frac threshold (snowmelt ≤ this)",
    y        = "Number of peaks classified as snowmelt"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(outdir, "snow_n_peaks.png"), p_s3,
       width = 8, height = 5, dpi = 150)
message("  Saved: snow_n_peaks.png")

# --- Plot S4: Melt flux contamination in snowmelt shoulder zone --------------
snow_shoulder <- peaks %>%
  filter(rain_frac >= SNOW_SHOULDER_LO,
         rain_frac <= SNOW_SHOULDER_HI,
         ros_flux_frac <= CURRENT_ROS_THRESH) %>%
  mutate(
    snow_bin = cut(rain_frac,
                   breaks = seq(SNOW_SHOULDER_LO, SNOW_SHOULDER_HI,
                                length.out = 5),
                   include.lowest = TRUE, right = FALSE),
    side = if_else(rain_frac <= CURRENT_SNOW_THRESH,
                   "At or below threshold\n(→ snowmelt)",
                   "Above threshold\n(→ mixed)")
  )

p_s4a <- ggplot(snow_shoulder, aes(x = melt_flux_short)) +
  geom_histogram(bins = 40, fill = "#377eb8", colour = "white",
                 linewidth = 0.2) +
  geom_vline(xintercept = 1.0, linetype = "dashed",
             colour = "grey30", linewidth = 0.7) +
  annotate("text", x = 1.05, y = Inf,
           label = "1 mm", hjust = 0, vjust = 1.5,
           size = 3, colour = "grey30") +
  facet_wrap(~ side, ncol = 2) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 5, 10, 25, 50, 100, 250),
                     labels = c("0", "1", "5", "10", "25", "50", "100", "250")) +
  labs(
    title    = sprintf("Melt flux in snowmelt shoulder zone (%.2f ≤ rain_frac ≤ %.2f)",
                       SNOW_SHOULDER_LO, SNOW_SHOULDER_HI),
    subtitle = sprintf("n = %d peaks; x-axis log1p scale", nrow(snow_shoulder)),
    x        = "melt_flux_short (mm, log1p scale)",
    y        = "Count"
  ) +
  theme_bw(base_size = 11)

p_s4b <- ggplot(snow_shoulder %>% filter(!is.na(snow_bin)),
                aes(x = snow_bin, y = melt_flux_short)) +
  geom_boxplot(fill = "#377eb8", alpha = 0.5, outlier.size = 0.8,
               outlier.alpha = 0.4, linewidth = 0.4) +
  geom_hline(yintercept = 1.0, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 5, 10, 25, 50, 100, 250),
                     labels = c("0", "1", "5", "10", "25", "50", "100", "250")) +
  labs(
    title    = "Melt flux by rain_frac bin — snowmelt shoulder zone",
    subtitle = "y-axis log1p scale",
    x        = "rain_frac bin",
    y        = "melt_flux_short (mm)"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_s4 <- p_s4a / p_s4b +
  plot_annotation(
    title   = "Melt flux contamination in snowmelt/mixed boundary zone",
    caption = "High melt flux above threshold = mixed events are genuinely snow-influenced."
  )

ggsave(file.path(outdir, "snow_melt_contamination.png"), p_s4,
       width = 8, height = 9, dpi = 150)
message("  Saved: snow_melt_contamination.png")

# --- Plots S5/S6: GPD xi for snowmelt type vs snowmelt threshold -------------
message("\n  Scanning GPD xi across snowmelt threshold range...")

snow_gpd_results <- map_dfr(snow_thresholds, function(t) {
  peaks %>%
    mutate(flood_type_t = case_when(
      ros_flux_frac >  CURRENT_ROS_THRESH ~ "ROS",
      rain_frac     <= t                  ~ "snowmelt",
      TRUE                                ~ "mixed"
    )) %>%
    filter(flood_type_t == "snowmelt") %>%
    group_by(site_no) %>%
    summarise(
      xi = fit_gpd_xi(peak_flow, first(pot_threshold_cfs)),
      .groups = "drop"
    ) %>%
    filter(!is.na(xi)) %>%
    summarise(
      n_sites   = n(),
      xi_median = median(xi),
      xi_q25    = quantile(xi, 0.25),
      xi_q75    = quantile(xi, 0.75)
    ) %>%
    mutate(threshold = t)
}, .progress = TRUE)

p_s5 <- ggplot(snow_gpd_results, aes(x = threshold)) +
  geom_ribbon(aes(ymin = xi_q25, ymax = xi_q75),
              fill = "#377eb8", alpha = 0.2) +
  geom_line(aes(y = xi_median), colour = "#377eb8", linewidth = 0.9) +
  geom_point(aes(y = xi_median), colour = "#377eb8", size = 1.5) +
  geom_vline(xintercept = CURRENT_SNOW_THRESH,
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  scale_x_continuous(breaks = seq(SNOW_THRESH_MIN, SNOW_THRESH_MAX, by = 0.02)) +
  labs(
    title    = "Snowmelt GPD xi vs. rain_frac threshold",
    subtitle = "Median ± IQR across sites",
    x        = "rain_frac threshold (snowmelt ≤ this)",
    y        = "GPD shape parameter (xi)",
    caption  = sprintf("Sites with >= %d snowmelt peaks", MIN_TYPE_SAMPLE)
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(outdir, "snow_gpd_xi.png"), p_s5,
       width = 8, height = 5, dpi = 150)
message("  Saved: snow_gpd_xi.png")

# Per-site spaghetti for snowmelt type
message("  Computing per-site snowmelt xi trajectories...")
snow_site_xi <- map_dfr(snow_thresholds, function(t) {
  peaks %>%
    mutate(flood_type_t = case_when(
      ros_flux_frac >  CURRENT_ROS_THRESH ~ "ROS",
      rain_frac     <= t                  ~ "snowmelt",
      TRUE                                ~ "mixed"
    )) %>%
    filter(flood_type_t == "snowmelt") %>%
    group_by(site_no) %>%
    summarise(xi = fit_gpd_xi(peak_flow, first(pot_threshold_cfs)),
              .groups = "drop") %>%
    filter(!is.na(xi)) %>%
    mutate(threshold = t)
}, .progress = TRUE)

p_s6 <- ggplot(snow_site_xi, aes(x = threshold, y = xi, group = site_no)) +
  geom_line(alpha = 0.2, colour = "#377eb8", linewidth = 0.4) +
  geom_vline(xintercept = CURRENT_SNOW_THRESH,
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  scale_x_continuous(breaks = seq(SNOW_THRESH_MIN, SNOW_THRESH_MAX, by = 0.02)) +
  labs(
    title    = "Per-site GPD xi trajectories (snowmelt type) vs. rain_frac threshold",
    subtitle = "Each line = one site",
    x        = "rain_frac threshold",
    y        = "GPD shape parameter (xi)"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(outdir, "snow_gpd_xi_by_site.png"), p_s6,
       width = 8, height = 5, dpi = 150)
message("  Saved: snow_gpd_xi_by_site.png")

# ============================================================================
# SUMMARY
# ============================================================================
ros_xi_cv  <- sd(ros_gpd_results %>% filter(type == "ROS") %>%
                   pull(xi_median), na.rm = TRUE) /
              abs(mean(ros_gpd_results %>% filter(type == "ROS") %>%
                         pull(xi_median), na.rm = TRUE))

snow_xi_cv <- sd(snow_gpd_results$xi_median, na.rm = TRUE) /
              abs(mean(snow_gpd_results$xi_median, na.rm = TRUE))

ros_gap_summary <- if (nrow(ros_gaps) > 0) {
  paste(sprintf("  %.3f (density = %.4f)", ros_gaps$x, ros_gaps$y),
        collapse = "\n")
} else "  None detected"

snow_gap_summary <- if (nrow(snow_gaps) > 0) {
  paste(sprintf("  %.3f (density = %.4f)", snow_gaps$x, snow_gaps$y),
        collapse = "\n")
} else "  None detected"

summary_lines <- c(
  "=== Non-Rainfall Threshold Sensitivity Summary ===",
  sprintf("Non-rainfall peaks in scope : %d (%d sites)",
          nrow(peaks), n_distinct(peaks$site_no)),
  sprintf("Rainfall threshold held fixed at : %.2f", RAIN_EM_THRESH),
  "",
  "============================",
  "Part 1: ROS Threshold",
  "============================",
  sprintf("Scan range    : %.2f to %.2f (step %.2f)",
          ROS_THRESH_MIN, ROS_THRESH_MAX, ROS_THRESH_STEP),
  sprintf("Current value : %.2f", CURRENT_ROS_THRESH),
  sprintf("Peaks with ros_flux_frac > %.2f (current): %d",
          CURRENT_ROS_THRESH, sum(peaks$ros_flux_frac > CURRENT_ROS_THRESH)),
  "Local density minima in scan zone:",
  ros_gap_summary,
  sprintf("ROS xi_median CV across threshold range: %.3f", ros_xi_cv),
  "",
  sprintf("Shoulder zone rain flux (%.2f to %.2f):",
          ROS_SHOULDER_LO, ROS_SHOULDER_HI),
  sprintf("  n = %d peaks", nrow(ros_shoulder)),
  capture.output(print(summary(ros_shoulder$rain_flux_short))),
  sprintf("  Fraction with rain_flux_short < 1mm  : %.1f%%",
          100 * mean(ros_shoulder$rain_flux_short < 1.0, na.rm = TRUE)),
  sprintf("  Fraction with rain_flux_short < 5mm  : %.1f%%",
          100 * mean(ros_shoulder$rain_flux_short < 5.0, na.rm = TRUE)),
  sprintf("  Fraction with rain_flux_short < 10mm : %.1f%%",
          100 * mean(ros_shoulder$rain_flux_short < 10.0, na.rm = TRUE)),
  "",
  "============================",
  "Part 2: Snowmelt Threshold",
  "============================",
  sprintf("Scan range    : %.2f to %.2f (step %.2f)",
          SNOW_THRESH_MIN, SNOW_THRESH_MAX, SNOW_THRESH_STEP),
  sprintf("Current value : %.2f", CURRENT_SNOW_THRESH),
  sprintf("Peaks with rain_frac <= %.2f (current, after ROS): %d",
          CURRENT_SNOW_THRESH,
          sum(peaks$rain_frac <= CURRENT_SNOW_THRESH &
                peaks$ros_flux_frac <= CURRENT_ROS_THRESH)),
  "Local density minima in scan zone:",
  snow_gap_summary,
  sprintf("Snowmelt xi_median CV across threshold range: %.3f", snow_xi_cv),
  "",
  sprintf("Shoulder zone melt flux (%.2f to %.2f):",
          SNOW_SHOULDER_LO, SNOW_SHOULDER_HI),
  sprintf("  n = %d peaks", nrow(snow_shoulder)),
  capture.output(print(summary(snow_shoulder$melt_flux_short))),
  sprintf("  Fraction with melt_flux_short < 1mm  : %.1f%%",
          100 * mean(snow_shoulder$melt_flux_short < 1.0, na.rm = TRUE)),
  sprintf("  Fraction with melt_flux_short < 5mm  : %.1f%%",
          100 * mean(snow_shoulder$melt_flux_short < 5.0, na.rm = TRUE)),
  sprintf("  Fraction with melt_flux_short < 10mm : %.1f%%",
          100 * mean(snow_shoulder$melt_flux_short < 10.0, na.rm = TRUE))
)

writeLines(summary_lines, file.path(outdir, "nonrainfall_thresh_summary.txt"))
message("  Saved: nonrainfall_thresh_summary.txt")

cat(paste(summary_lines, collapse = "\n"), "\n")

message("\nDiagnostic complete. Outputs in: ", outdir)

# ============================================================================
# Diagnostic: 2D Snowmelt Threshold Sensitivity Surface
# ============================================================================
# Varies both the ROS threshold (ros_flux_frac) and the snowmelt end-member
# threshold (rain_frac) jointly, producing a 2D surface of:
#
#   - Snowmelt peak count
#   - GPD shape parameter (xi) for the snowmelt type
#
# This reveals whether there is a stable plateau in xi space that supports
# a particular combination of thresholds, or whether the two thresholds
# are truly independent (i.e. optimizing one doesn't depend on the other).
#
# SCOPE: Unregulated, data-sufficient, non-rainfall peaks
#        (rain_frac < RAIN_EM_THRESH = 0.95)
#
# INPUT:
#   - _data/pot_peaks_classified_adaptive.rds
#
# OUTPUTS (in _diagnostics/nonrainfall_threshold_sensitivity/):
#   - snow_2d_count.png     : heatmap of snowmelt peak count across threshold grid
#   - snow_2d_xi.png        : heatmap of snowmelt GPD xi across threshold grid
#   - snow_2d_xi_stable.png : xi heatmap masked to stable region (n_sites >= N)
# ============================================================================

library(tidyverse)
library(dplyr)
library(lmomco)

classified_file <- "_data/pot_peaks_classified_adaptive.rds"
outdir          <- "_diagnostics/nonrainfall_threshold_sensitivity"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

RAIN_EM_THRESH  <- 0.95
MIN_TYPE_SAMPLE <- 20

# Grid
ROS_THRESHOLDS  <- seq(0.25, 0.75, by = 0.05)
SNOW_THRESHOLDS <- seq(0.02, 0.20, by = 0.02)

# Current values — highlighted on plots
CURRENT_ROS  <- 0.50
CURRENT_SNOW <- 0.05

# ============================================================================
# LOAD
# ============================================================================
message("Loading data...")
dat <- readRDS(classified_file)

peaks <- dat %>%
  filter(
    !regulated,
    data_sufficient,
    rain_frac < RAIN_EM_THRESH,
    !is.na(rain_frac),
    !is.na(ros_flux_frac)
  )

message(sprintf("Non-rainfall peaks in scope: %d across %d sites\n",
                nrow(peaks), n_distinct(peaks$site_no)))

fit_gpd_xi <- function(x, u, min_n = MIN_TYPE_SAMPLE) {
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

# ============================================================================
# COMPUTE GRID
# ============================================================================
message(sprintf("Computing %d x %d threshold grid (%d combinations)...",
                length(ROS_THRESHOLDS), length(SNOW_THRESHOLDS),
                length(ROS_THRESHOLDS) * length(SNOW_THRESHOLDS)))

grid_results <- expand_grid(
  ros_thresh  = ROS_THRESHOLDS,
  snow_thresh = SNOW_THRESHOLDS
) %>%
  mutate(row_id = row_number())

results <- map_dfr(seq_len(nrow(grid_results)), function(i) {
  rt <- grid_results$ros_thresh[i]
  st <- grid_results$snow_thresh[i]

  snowmelt_peaks <- peaks %>%
    filter(
      ros_flux_frac <= rt,
      rain_frac     <= st
    )

  site_fits <- snowmelt_peaks %>%
    group_by(site_no) %>%
    summarise(
      xi = fit_gpd_xi(peak_flow, first(pot_threshold_cfs)),
      n  = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(xi))

  tibble(
    ros_thresh  = rt,
    snow_thresh = st,
    n_peaks     = nrow(snowmelt_peaks),
    n_sites     = nrow(site_fits),
    xi_median   = if (nrow(site_fits) > 0) median(site_fits$xi) else NA_real_,
    xi_q25      = if (nrow(site_fits) > 0) quantile(site_fits$xi, 0.25) else NA_real_,
    xi_q75      = if (nrow(site_fits) > 0) quantile(site_fits$xi, 0.75) else NA_real_
  )
}, .progress = TRUE)

message(sprintf("Grid complete: %d combinations", nrow(results)))

# ============================================================================
# PLOTS
# ============================================================================

# Shared theme elements
thresh_annotation <- function(p) {
  p +
    geom_vline(xintercept = CURRENT_SNOW, linetype = "dashed",
               colour = "white", linewidth = 0.8) +
    geom_hline(yintercept = CURRENT_ROS,  linetype = "dashed",
               colour = "white", linewidth = 0.8) +
    annotate("point", x = CURRENT_SNOW, y = CURRENT_ROS,
             shape = 3, size = 4, stroke = 1.5, colour = "white") +
    annotate("text", x = CURRENT_SNOW + 0.005, y = CURRENT_ROS + 0.025,
             label = sprintf("Current\n(%.2f, %.2f)", CURRENT_SNOW, CURRENT_ROS),
             hjust = 0, size = 2.8, colour = "white") +
    scale_x_continuous(breaks = SNOW_THRESHOLDS,
                       labels = sprintf("%.2f", SNOW_THRESHOLDS)) +
    scale_y_continuous(breaks = ROS_THRESHOLDS,
                       labels = sprintf("%.2f", ROS_THRESHOLDS)) +
    labs(
      x = "Snowmelt threshold (rain_frac ≤ this)",
      y = "ROS threshold (ros_flux_frac > this)"
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid  = element_blank())
}

# --- Plot 1: Peak count surface ----------------------------------------------
p1 <- ggplot(results, aes(x = snow_thresh, y = ros_thresh, fill = n_peaks)) +
  geom_tile() +
  geom_text(aes(label = n_peaks), size = 2.5, colour = "white") +
  scale_fill_viridis_c(option = "plasma", name = "n peaks",
                       direction = 1) +
  labs(
    title    = "Snowmelt peak count across threshold grid",
    subtitle = "Dashed lines and cross = current 6b thresholds (0.05, 0.50)"
  )

p1 <- thresh_annotation(p1)

ggsave(file.path(outdir, "snow_2d_count.png"), p1,
       width = 10, height = 7, dpi = 150)
message("Saved: snow_2d_count.png")

# --- Plot 2: xi median surface -----------------------------------------------
xi_lim <- quantile(abs(results$xi_median), 0.95, na.rm = TRUE)

p2 <- ggplot(results, aes(x = snow_thresh, y = ros_thresh, fill = xi_median)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.0f", xi_median)),
            size = 2.3, colour = "white") +
  scale_fill_viridis_c(option = "magma", name = "xi median",
                       limits = c(-xi_lim, 0),
                       oob = scales::squish) +
  labs(
    title    = "Snowmelt GPD shape parameter (xi) across threshold grid",
    subtitle = "Median across sites; values near 0 = less bounded upper tail"
  )

p2 <- thresh_annotation(p2)

ggsave(file.path(outdir, "snow_2d_xi.png"), p2,
       width = 10, height = 7, dpi = 150)
message("Saved: snow_2d_xi.png")

# --- Plot 3: xi surface masked to cells with n_sites >= threshold -----------
MIN_SITES_FOR_STABLE <- max(5, round(n_distinct(peaks$site_no) * 0.10))
message(sprintf("Masking to cells with n_sites >= %d (10%% of sites)",
                MIN_SITES_FOR_STABLE))

results_stable <- results %>%
  mutate(xi_median_masked = if_else(n_sites >= MIN_SITES_FOR_STABLE,
                                    xi_median, NA_real_))

p3 <- ggplot(results_stable,
             aes(x = snow_thresh, y = ros_thresh, fill = xi_median_masked)) +
  geom_tile() +
  geom_text(aes(label = if_else(!is.na(xi_median_masked),
                                sprintf("%.0f", xi_median_masked), "—")),
            size = 2.3, colour = "white") +
  scale_fill_viridis_c(option = "magma", name = "xi median",
                       na.value = "grey30",
                       limits = c(-xi_lim, 0),
                       oob = scales::squish) +
  labs(
    title    = sprintf("Snowmelt GPD xi — stable cells only (n_sites ≥ %d)",
                       MIN_SITES_FOR_STABLE),
    subtitle = "Grey cells = insufficient sites for reliable xi estimate"
  )

p3 <- thresh_annotation(p3)

ggsave(file.path(outdir, "snow_2d_xi_stable.png"), p3,
       width = 10, height = 7, dpi = 150)
message("Saved: snow_2d_xi_stable.png")

# ============================================================================
# CONSOLE SUMMARY
# ============================================================================
message("\n=== 2D Grid Summary ===")
message(sprintf("Grid: ros_thresh [%.2f, %.2f], snow_thresh [%.2f, %.2f]",
                min(ROS_THRESHOLDS), max(ROS_THRESHOLDS),
                min(SNOW_THRESHOLDS), max(SNOW_THRESHOLDS)))

current_cell <- results %>%
  filter(ros_thresh == CURRENT_ROS, snow_thresh == CURRENT_SNOW)
message(sprintf("\nCurrent thresholds (ros=%.2f, snow=%.2f):",
                CURRENT_ROS, CURRENT_SNOW))
message(sprintf("  n_peaks = %d, n_sites = %d, xi_median = %.1f",
                current_cell$n_peaks, current_cell$n_sites,
                current_cell$xi_median))

# Find most stable region (smallest IQR of xi)
stable_cells <- results %>%
  filter(!is.na(xi_median), n_sites >= MIN_SITES_FOR_STABLE) %>%
  mutate(xi_iqr = xi_q75 - xi_q25) %>%
  arrange(xi_iqr)

message("\nTop 5 most stable threshold combinations (smallest xi IQR):")
print(stable_cells %>%
        select(ros_thresh, snow_thresh, n_peaks, n_sites, xi_median, xi_iqr) %>%
        head(5))

message("\nDiagnostic complete. Outputs in: ", outdir)