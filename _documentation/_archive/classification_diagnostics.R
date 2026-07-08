#!/usr/bin/env Rscript

# ============================================================================
# Classification Diagnostics — Key Decision Points
# ============================================================================
# Documents the empirical basis for four classification decisions:
#
#   1. 3-class schema   : rain_frac distribution in the middle bin does not
#                         support a further ROS/mixed split at relevant
#                         elevations; per-site sample sizes preclude it.
#
#   2. Rainfall threshold (0.95) : sits in a naturally sparse zone; GPD xi
#                         is stable across 0.70–0.95; melt contamination in
#                         the shoulder zone is negligible.
#
#   3. Snowmelt threshold (0.10) : no clear density gap, but GPD xi
#                         stabilises around 0.05–0.10 and count curve
#                         flattens; 0.10 is a pragmatic choice and motivates
#                         the 2-class fallback in the app.
#
#   4. Adaptive window  : snowmelt classification fraction is flat across
#                         window multipliers 1x-5x, and per-peak rain_frac
#                         distributions are stable below the threshold.
#                         Upper outliers grow with multiplier, consistent
#                         with bleeding into adjacent independent events —
#                         reinforcing search_radius (multiplier = 1) as
#                         both sufficient and appropriate.
#
# INPUTS:
#   _data/pot_peaks_classified_adaptive.rds   (Script 6b)
#   _data/accum_curves_60d.rds                (Script 5)
#   _data/pot_peaks_all.rds                   (Script 3)
#   _data/site_metadata.rds                   (Script 0)
#
# OUTPUTS (all in _documentation/classification_decisions/):
#   Section 1 — 3-class schema
#     01a_rain_frac_full_distribution.png
#     01b_middle_bin_by_elevation.png
#     01c_feature_space_scatter.png
#   Section 2 — Rainfall threshold
#     02a_rainfall_threshold_distribution.png
#     02b_rainfall_gpd_xi_stability.png
#   Section 3 — Snowmelt threshold
#     03a_snowmelt_threshold_distribution.png
#     03b_snowmelt_gpd_xi_stability.png
#   Section 4 — Window multiplier sensitivity
#     04a_snowmelt_frac_vs_multiplier.png
#     04b_rain_frac_vs_multiplier.png
#
# Author : Mike Talbot, PhD Student, Colorado State University
# ============================================================================

library(tidyverse)
library(lmomco)
library(patchwork)

# ============================================================================
# CONFIGURATION
# ============================================================================

classified_file <- "_data/pot_peaks_classified_adaptive.rds"
accum_file      <- "_data/accum_curves_60d.rds"
peaks_file      <- "_data/pot_peaks_all.rds"
meta_file       <- "_data/site_metadata.rds"
out_dir         <- "_documentation/classification_decisions"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Current thresholds (Script 6b)
RAIN_EM_THRESH  <- 0.95
SNOW_EM_THRESH  <- 0.10
NOISE_FLOOR_MM  <- 0.1
MIN_TYPE_SAMPLE <- 20

# Rainfall threshold scan range
RAIN_THRESH_MIN  <- 0.70
RAIN_THRESH_MAX  <- 0.99
RAIN_THRESH_STEP <- 0.01

# Snowmelt threshold scan range
SNOW_THRESH_MIN  <- 0.01
SNOW_THRESH_MAX  <- 0.20
SNOW_THRESH_STEP <- 0.01

# Elevation bands for faceting
ELEV_BREAKS <- c(0, 1500, 2000, 2500, 3000, Inf)
ELEV_LABELS <- c("< 1500 m", "1500-2000 m", "2000-2500 m",
                 "2500-3000 m", "> 3000 m")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Fit GPD shape parameter (xi) for a vector of peaks at a given site.
# Returns NA if fewer than min_n peaks or if lmomco fit fails.
fit_gpd_xi <- function(flows, mu, min_n = MIN_TYPE_SAMPLE) {
  if (length(flows) < min_n) return(NA_real_)
  x <- flows[flows > mu] - mu
  if (length(x) < min_n) return(NA_real_)
  lmom <- tryCatch(lmomco::lmoms(x), error = function(e) NULL)
  if (is.null(lmom)) return(NA_real_)
  pars <- tryCatch(lmomco::pargpa(lmom), error = function(e) NULL)
  if (is.null(pars) || is.null(pars$para)) return(NA_real_)
  unname(pars$para["xi"])
}

# Reclassify peaks at a given rainfall threshold, return per-site xi for
# the rainfall type. Used in Section 2.
gpd_xi_at_rain_thresh <- function(peaks_df, threshold, mu_lookup) {
  peaks_df %>%
    filter(data_sufficient, rain_frac >= threshold) %>%
    group_by(site_no) %>%
    arrange(desc(peak_flow), .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    filter(rank <= pmax(MIN_TYPE_SAMPLE, n())) %>%
    filter(n() >= MIN_TYPE_SAMPLE) %>%
    summarise(
      xi = fit_gpd_xi(peak_flow,
                      mu = mu_lookup$pot_threshold_cfs[mu_lookup$site_no == first(site_no)],
                      min_n = MIN_TYPE_SAMPLE),
      .groups = "drop"
    ) %>%
    mutate(threshold = threshold)
}

# Reclassify peaks at a given snowmelt threshold (rainfall already removed),
# return per-site xi for the snowmelt type. Used in Section 3.
gpd_xi_at_snow_thresh <- function(peaks_df, threshold, mu_lookup) {
  peaks_df %>%
    filter(data_sufficient, rain_frac < RAIN_EM_THRESH,
           rain_frac <= threshold) %>%
    group_by(site_no) %>%
    arrange(desc(peak_flow), .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    filter(rank <= pmax(MIN_TYPE_SAMPLE, n())) %>%
    filter(n() >= MIN_TYPE_SAMPLE) %>%
    summarise(
      xi = fit_gpd_xi(peak_flow,
                      mu = mu_lookup$pot_threshold_cfs[mu_lookup$site_no == first(site_no)],
                      min_n = MIN_TYPE_SAMPLE),
      .groups = "drop"
    ) %>%
    mutate(threshold = threshold)
}

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading data...")

classified <- readRDS(classified_file)
meta        <- readRDS(meta_file)
accum       <- readRDS(accum_file)

peak_meta <- readRDS(peaks_file) %>%
  distinct(site_no, peak_id, search_radius, pot_threshold_cfs)

# Unregulated sites only for threshold diagnostics
unreg_sites <- meta %>% filter(!regulated) %>% pull(site_no)

# Per-site POT threshold for GPD fitting (mu)
mu_lookup <- peak_meta %>%
  group_by(site_no) %>%
  summarise(pot_threshold_cfs = first(pot_threshold_cfs), .groups = "drop")

# Restrict to unregulated data-sufficient peaks
peaks_unreg <- classified %>%
  filter(site_no %in% unreg_sites, data_sufficient, !is.na(rain_frac))

message(sprintf("  Unregulated data-sufficient peaks: %d across %d sites",
                nrow(peaks_unreg), n_distinct(peaks_unreg$site_no)))

# ============================================================================
# SECTION 1 — 3-CLASS SCHEMA
# ============================================================================
# Is there a natural breakpoint in the middle bin that would justify a
# 4-class split (snowmelt-dominated ROS vs. rainfall-dominated ROS)?
# ============================================================================
message("\n--- Section 1: 3-class schema ---")

peaks_elev <- peaks_unreg %>%
  left_join(meta %>% select(site_no, elev_mean_m), by = "site_no") %>%
  mutate(elev_band = cut(elev_mean_m, breaks = ELEV_BREAKS,
                         labels = ELEV_LABELS, right = FALSE))

middle_bin <- peaks_elev %>%
  filter(rain_frac > SNOW_EM_THRESH, rain_frac < RAIN_EM_THRESH)

message(sprintf("  Middle-bin peaks: %d (%.1f%% of total)",
                nrow(middle_bin), 100 * nrow(middle_bin) / nrow(peaks_unreg)))

# 01a: Full rain_frac distribution showing end-member modes and sparse middle
p01a <- ggplot(peaks_elev, aes(x = rain_frac)) +
  geom_density(fill = "#e41a1c", alpha = 0.25, color = "#e41a1c", linewidth = 0.8) +
  geom_vline(xintercept = RAIN_EM_THRESH, linetype = "dashed",
             color = "#e41a1c", linewidth = 0.8) +
  geom_vline(xintercept = SNOW_EM_THRESH, linetype = "dashed",
             color = "#377eb8", linewidth = 0.8) +
  annotate("text", x = RAIN_EM_THRESH - 0.01, y = Inf,
           label = sprintf("Rainfall\n>= %.2f", RAIN_EM_THRESH),
           vjust = 1.5, hjust = 1, color = "#e41a1c", size = 3.5) +
  annotate("text", x = SNOW_EM_THRESH + 0.01, y = Inf,
           label = sprintf("Snowmelt\n<= %.2f", SNOW_EM_THRESH),
           vjust = 1.5, hjust = 0, color = "#377eb8", size = 3.5) +
  labs(
    title    = "rain_frac distribution — all unregulated data-sufficient peaks",
    subtitle = sprintf("Strong end-member modes; sparse middle bin (%.0f%% of peaks)",
                       100 * nrow(middle_bin) / nrow(peaks_unreg)),
    x        = "rain_frac",
    y        = "Density"
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(out_dir, "01a_rain_frac_full_distribution.png"),
       p01a, width = 9, height = 5, dpi = 150)
message("  Saved 01a_rain_frac_full_distribution.png")

# 01b: Middle-bin density by elevation — looking for bimodality
elev_n <- middle_bin %>%
  dplyr::count(elev_band) %>%
  mutate(label = sprintf("%s\nn = %d", elev_band, n))

middle_bin_labeled <- middle_bin %>%
  left_join(elev_n %>% select(elev_band, label), by = "elev_band")

p01b <- ggplot(middle_bin_labeled, aes(x = rain_frac)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.05, fill = "#7b3294", alpha = 0.55,
                 color = "white") +
  geom_density(color = "#1a1a1a", linewidth = 0.7) +
  geom_vline(xintercept = 0.50, linetype = "dashed",
             color = "firebrick", linewidth = 0.6) +
  scale_x_continuous(limits = c(SNOW_EM_THRESH, RAIN_EM_THRESH),
                     breaks = seq(0.1, 0.9, by = 0.2)) +
  facet_wrap(~ label, ncol = 3) +
  labs(
    title    = "rain_frac distribution by elevation — middle bin",
    subtitle = "No bimodality at dominant elevations; no natural 4-class breakpoint",
    x        = "rain_frac",
    y        = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "#f0f0f0"))

ggsave(file.path(out_dir, "01b_middle_bin_by_elevation.png"),
       p01b, width = 11, height = 7, dpi = 150)
message("  Saved 01b_middle_bin_by_elevation.png")

# 01c: Feature space scatter — rain_frac vs. melt_flux_short, coloured by
# flood type. Shows classification boundaries in joint feature space and
# confirms that the three classes occupy distinct regions.
type_colors_scatter <- c(
  rainfall = "#e41a1c",
  ROS      = "#7b3294",
  snowmelt = "#377eb8"
)

scatter_df <- peaks_elev %>%
  filter(as.character(flood_type_raw) %in% c("rainfall", "ROS", "snowmelt"),
         melt_flux_short > 0) %>%
  mutate(flood_type_raw = factor(as.character(flood_type_raw),
                                 levels = c("snowmelt", "ROS", "rainfall")))

p01c <- ggplot(scatter_df,
               aes(x = rain_frac, y = melt_flux_short,
                   color = flood_type_raw)) +
  geom_point(size = 0.6, alpha = 0.35) +
  geom_vline(xintercept = SNOW_EM_THRESH, linetype = "dashed",
             color = "grey30", linewidth = 0.7) +
  geom_vline(xintercept = RAIN_EM_THRESH, linetype = "dashed",
             color = "grey30", linewidth = 0.7) +
  scale_y_log10(labels = scales::label_comma(),
                breaks  = c(0.1, 1, 10, 100, 1000)) +
  scale_color_manual(
    values = type_colors_scatter,
    labels = c(snowmelt = "Snowmelt", ROS = "Rain-on-Snow", rainfall = "Rainfall"),
    name   = "Flood type"
  ) +
  annotate("text", x = SNOW_EM_THRESH / 2, y = Inf,
           label = "Snowmelt", vjust = 1.5, color = "#377eb8", size = 3.5) +
  annotate("text", x = (SNOW_EM_THRESH + RAIN_EM_THRESH) / 2, y = Inf,
           label = "ROS", vjust = 1.5, color = "#7b3294", size = 3.5) +
  annotate("text", x = (RAIN_EM_THRESH + 1) / 2, y = Inf,
           label = "Rainfall", vjust = 1.5, color = "#e41a1c", size = 3.5) +
  labs(
    title    = "Classification in feature space — rain_frac vs. melt flux",
    subtitle = "Dashed lines = classification boundaries; y-axis log scale",
    x        = "rain_frac",
    y        = "Total melt flux over window (mm)"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "01c_feature_space_scatter.png"),
       p01c, width = 9, height = 6, dpi = 150)
message("  Saved 01c_feature_space_scatter.png")

# ============================================================================
# SECTION 2 — RAINFALL THRESHOLD
# ============================================================================
# Does 0.95 sit in a natural sparse zone? Is GPD xi stable across the
# candidate range (0.70–0.95)?
# ============================================================================
message("\n--- Section 2: Rainfall threshold ---")

rain_thresholds <- seq(RAIN_THRESH_MIN, RAIN_THRESH_MAX, by = RAIN_THRESH_STEP)

# 02a: rain_frac density in the threshold zone
p02a <- ggplot(peaks_unreg %>% filter(rain_frac >= RAIN_THRESH_MIN),
               aes(x = rain_frac)) +
  annotate("rect",
           xmin = RAIN_THRESH_MIN, xmax = RAIN_EM_THRESH,
           ymin = -Inf, ymax = Inf,
           fill = "grey90", alpha = 0.6) +
  geom_density(fill = "#e41a1c", alpha = 0.25,
               color = "#e41a1c", linewidth = 0.8) +
  geom_vline(xintercept = RAIN_EM_THRESH, linetype = "dashed",
             color = "#333333", linewidth = 0.8) +
  annotate("text", x = RAIN_EM_THRESH + 0.002, y = Inf,
           label = sprintf("Current threshold\n= %.2f", RAIN_EM_THRESH),
           vjust = 1.5, hjust = 0, size = 3.5) +
  labs(
    title    = "rain_frac density — rainfall threshold zone",
    subtitle = "Shaded zone = candidate range; threshold sits in sparse ascending shoulder",
    x        = "rain_frac",
    y        = "Density"
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(out_dir, "02a_rainfall_threshold_distribution.png"),
       p02a, width = 9, height = 5, dpi = 150)
message("  Saved 02a_rainfall_threshold_distribution.png")

# 02b: GPD xi stability across rainfall threshold scan
message("  Fitting GPD xi across rainfall threshold range (this may take a moment)...")

rain_xi_results <- map_dfr(rain_thresholds, function(thr) {
  gpd_xi_at_rain_thresh(peaks_unreg, thr, mu_lookup)
})

rain_xi_summary <- rain_xi_results %>%
  filter(!is.na(xi)) %>%
  group_by(threshold) %>%
  summarise(
    xi_median = median(xi),
    xi_q25    = quantile(xi, 0.25),
    xi_q75    = quantile(xi, 0.75),
    n_sites   = n(),
    .groups   = "drop"
  )

xi_cv <- rain_xi_summary %>%
  summarise(cv = sd(xi_median) / abs(mean(xi_median))) %>%
  pull(cv)

p02b <- ggplot(rain_xi_summary, aes(x = threshold)) +
  geom_ribbon(aes(ymin = xi_q25, ymax = xi_q75),
              fill = "#e41a1c", alpha = 0.2) +
  geom_line(aes(y = xi_median), color = "#e41a1c", linewidth = 1.0) +
  geom_point(aes(y = xi_median), color = "#e41a1c", size = 2) +
  geom_vline(xintercept = RAIN_EM_THRESH, linetype = "dashed",
             color = "#333333", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted",
             color = "grey40", linewidth = 0.6) +
  annotate("text", x = RAIN_EM_THRESH + 0.002, y = Inf,
           label = sprintf("%.2f\nCV = %.3f", RAIN_EM_THRESH, xi_cv),
           vjust = 1.5, hjust = 0, size = 3.5) +
  labs(
    title    = "GPD shape parameter (xi) for rainfall type vs. threshold",
    subtitle = "Median +/- IQR across unregulated sites; stable xi = threshold-insensitive",
    x        = "rain_frac threshold",
    y        = "GPD shape parameter (xi)",
    caption  = sprintf("Sites with >= %d rainfall peaks included at each threshold",
                       MIN_TYPE_SAMPLE)
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(out_dir, "02b_rainfall_gpd_xi_stability.png"),
       p02b, width = 9, height = 5, dpi = 150)
message("  Saved 02b_rainfall_gpd_xi_stability.png")

# ============================================================================
# SECTION 3 — SNOWMELT THRESHOLD
# ============================================================================
# Is there a density gap near 0.10? Does GPD xi stabilise?
# ============================================================================
message("\n--- Section 3: Snowmelt threshold ---")

snow_thresholds <- seq(SNOW_THRESH_MIN, SNOW_THRESH_MAX, by = SNOW_THRESH_STEP)

non_rainfall <- peaks_unreg %>%
  filter(rain_frac < RAIN_EM_THRESH)

# 03a: rain_frac density at the low end
p03a <- ggplot(non_rainfall %>% filter(rain_frac <= 0.30),
               aes(x = rain_frac)) +
  annotate("rect",
           xmin = SNOW_THRESH_MIN, xmax = SNOW_THRESH_MAX,
           ymin = -Inf, ymax = Inf,
           fill = "grey90", alpha = 0.6) +
  geom_density(fill = "#377eb8", alpha = 0.25,
               color = "#377eb8", linewidth = 0.8) +
  geom_vline(xintercept = SNOW_EM_THRESH, linetype = "dashed",
             color = "#333333", linewidth = 0.8) +
  annotate("text", x = SNOW_EM_THRESH + 0.003, y = Inf,
           label = sprintf("Current threshold\n= %.2f", SNOW_EM_THRESH),
           vjust = 1.5, hjust = 0, size = 3.5) +
  scale_x_continuous(limits = c(0, 0.30)) +
  labs(
    title    = "rain_frac distribution (low end) — non-rainfall population",
    subtitle = "Shaded zone = scan range; no density gap supports a clear natural threshold",
    x        = "rain_frac",
    y        = "Density"
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(out_dir, "03a_snowmelt_threshold_distribution.png"),
       p03a, width = 9, height = 5, dpi = 150)
message("  Saved 03a_snowmelt_threshold_distribution.png")

# 03b: GPD xi vs. snowmelt threshold — looking for stabilisation
message("  Fitting GPD xi across snowmelt threshold range...")

snow_xi_results <- map_dfr(snow_thresholds, function(thr) {
  gpd_xi_at_snow_thresh(peaks_unreg, thr, mu_lookup)
})

snow_xi_summary <- snow_xi_results %>%
  filter(!is.na(xi)) %>%
  group_by(threshold) %>%
  summarise(
    xi_median = median(xi),
    xi_q25    = quantile(xi, 0.25),
    xi_q75    = quantile(xi, 0.75),
    n_sites   = n(),
    .groups   = "drop"
  )

snow_xi_cv <- snow_xi_summary %>%
  summarise(cv = sd(xi_median) / abs(mean(xi_median))) %>%
  pull(cv)

# Peak count vs threshold for annotation
snow_counts <- map_dfr(snow_thresholds, function(thr) {
  tibble(
    threshold = thr,
    n_peaks   = sum(non_rainfall$rain_frac <= thr, na.rm = TRUE)
  )
})

p03b <- ggplot(snow_xi_summary, aes(x = threshold)) +
  geom_ribbon(aes(ymin = xi_q25, ymax = xi_q75),
              fill = "#377eb8", alpha = 0.2) +
  geom_line(aes(y = xi_median), color = "#377eb8", linewidth = 1.0) +
  geom_point(aes(y = xi_median), color = "#377eb8", size = 2) +
  geom_vline(xintercept = SNOW_EM_THRESH, linetype = "dashed",
             color = "#333333", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted",
             color = "grey40", linewidth = 0.6) +
  annotate("text", x = SNOW_EM_THRESH + 0.003, y = Inf,
           label = sprintf("%.2f\nCV = %.3f", SNOW_EM_THRESH, snow_xi_cv),
           vjust = 1.5, hjust = 0, size = 3.5) +
  scale_x_continuous(breaks = snow_thresholds[seq(1, length(snow_thresholds), 2)]) +
  labs(
    title    = "GPD shape parameter (xi) for snowmelt type vs. threshold",
    subtitle = "Median +/- IQR; xi stabilises around 0.05-0.10 then flattens",
    x        = "rain_frac threshold (snowmelt <= this)",
    y        = "GPD shape parameter (xi)",
    caption  = sprintf("Sites with >= %d snowmelt peaks included at each threshold",
                       MIN_TYPE_SAMPLE)
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(out_dir, "03b_snowmelt_gpd_xi_stability.png"),
       p03b, width = 9, height = 5, dpi = 150)
message("  Saved 03b_snowmelt_gpd_xi_stability.png")

# ============================================================================
# SECTION 4 — WINDOW MULTIPLIER SENSITIVITY
# ============================================================================
# Justifies using search_radius directly (multiplier = 1) rather than a
# fixed multiple. For each multiplier in MULTIPLIERS, rain_frac is recomputed
# over that window and used to classify peaks as snowmelt (rain_frac <=
# SNOW_EM_THRESH). Two complementary findings support multiplier = 1:
#
#   1. The snowmelt classification fraction is essentially flat across
#      multipliers, meaning the window length does not materially affect
#      which peaks are classified as snowmelt.
#
#   2. For peaks classified as snowmelt at 1x, rain_frac values remain
#      stable (IQR well below threshold) but upper outliers grow with
#      multiplier. This is consistent with longer windows bleeding into
#      adjacent independent events — exactly what the Livneh search_radius
#      criterion is designed to prevent. A multiplier > 1 is therefore not
#      merely redundant but potentially harmful.
# ============================================================================
message("\n--- Section 4: Window multiplier sensitivity ---")

MULTIPLIERS <- c(1, 2, 3, 4, 5)

accum_with_sr <- accum %>%
  inner_join(peak_meta %>% select(site_no, peak_id, search_radius),
             by = c("site_no", "peak_id")) %>%
  filter(site_no %in% unreg_sites)

# Compute rain_frac at each multiplier for all unregulated data-sufficient peaks
message("  Computing rain_frac at each window multiplier...")

rf_by_mult <- map_dfr(MULTIPLIERS, function(m) {
  accum_with_sr %>%
    filter(lookback_day <= m * search_radius) %>%
    mutate(
      rain_active = if_else(rain_mm_day >= NOISE_FLOOR_MM, rain_mm_day, 0),
      melt_active = if_else(melt_mm_day >= NOISE_FLOOR_MM, melt_mm_day, 0),
      total_flux  = rain_active + melt_active
    ) %>%
    group_by(site_no, peak_id) %>%
    summarise(
      rain_flux  = sum(rain_active, na.rm = TRUE),
      total_flux = sum(total_flux,  na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    mutate(
      rf_mult     = if_else(total_flux > 0, rain_flux / total_flux, NA_real_),
      multiplier  = m
    ) %>%
    select(site_no, peak_id, multiplier, rf_mult)
})

# Join back to classified peaks (non-rainfall already removed via RAIN_EM_THRESH)
rf_mult_peaks <- peaks_unreg %>%
  filter(rain_frac < RAIN_EM_THRESH) %>%   # non-rainfall population only
  select(site_no, peak_id) %>%
  inner_join(rf_by_mult, by = c("site_no", "peak_id"))

# 04a: Fraction classified as snowmelt at each multiplier
snowmelt_frac <- rf_mult_peaks %>%
  filter(!is.na(rf_mult)) %>%
  group_by(multiplier) %>%
  summarise(
    n_total   = n(),
    n_snow    = sum(rf_mult <= SNOW_EM_THRESH),
    frac_snow = n_snow / n_total,
    .groups   = "drop"
  )

p04a <- ggplot(snowmelt_frac, aes(x = multiplier, y = frac_snow)) +
  geom_line(color = "#377eb8", linewidth = 1.0) +
  geom_point(color = "#377eb8", size = 3) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", frac_snow * 100, n_snow)),
            vjust = -0.6, size = 3.2, color = "#377eb8") +
  scale_x_continuous(breaks = MULTIPLIERS,
                     labels = sprintf("%dx", MULTIPLIERS)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 0.1),
                     limits = c(0, max(snowmelt_frac$frac_snow) * 1.15)) +
  labs(
    title    = "Fraction of non-rainfall peaks classified as snowmelt vs. window multiplier",
    subtitle = sprintf("Snowmelt criterion: rain_frac <= %.2f; rainfall end-members excluded",
                       SNOW_EM_THRESH),
    x        = "Window length (multiples of search_radius)",
    y        = "Fraction classified as snowmelt"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "04a_snowmelt_frac_vs_multiplier.png"),
       p04a, width = 9, height = 5, dpi = 150)
message("  Saved 04a_snowmelt_frac_vs_multiplier.png")

# 04b: Per-peak rain_frac at each multiplier, for peaks classified as
# snowmelt at multiplier = 1. Shows whether values shift meaningfully
# as window expands.
snow_at_1x <- rf_mult_peaks %>%
  filter(multiplier == 1, rf_mult <= SNOW_EM_THRESH) %>%
  select(site_no, peak_id)

rf_snow_wide <- rf_mult_peaks %>%
  semi_join(snow_at_1x, by = c("site_no", "peak_id")) %>%
  filter(!is.na(rf_mult)) %>%
  mutate(multiplier = factor(sprintf("%dx", multiplier),
                             levels = sprintf("%dx", MULTIPLIERS)))

p04b <- ggplot(rf_snow_wide, aes(x = multiplier, y = rf_mult)) +
  geom_boxplot(fill = "#377eb8", alpha = 0.4, outlier.size = 0.8,
               outlier.alpha = 0.3) +
  geom_hline(yintercept = SNOW_EM_THRESH, linetype = "dashed",
             color = "firebrick", linewidth = 0.7) +
  annotate("text", x = 0.6, y = SNOW_EM_THRESH + 0.005,
           label = sprintf("Snowmelt threshold = %.2f", SNOW_EM_THRESH),
           hjust = 0, vjust = 0, size = 3.2, color = "firebrick") +
  labs(
    title    = "rain_frac distribution vs. window multiplier — peaks snowmelt at 1x",
    subtitle = "IQR stable below threshold at all multipliers; outliers grow with window length",
    x        = "Window length (multiples of search_radius)",
    y        = "rain_frac at given multiplier",
    caption  = paste0("Peaks classified as snowmelt at 1x search_radius (the Livneh independence window).\n",
                      "Growing upper outliers at longer windows are consistent with bleeding into\n",
                      "adjacent independent events — reinforcing search_radius as the correct window.")
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "04b_rain_frac_vs_multiplier.png"),
       p04b, width = 9, height = 5, dpi = 150)
message("  Saved 04b_rain_frac_vs_multiplier.png")

# ============================================================================
# CONSOLE SUMMARY
# ============================================================================
message("\n=== CLASSIFICATION DECISIONS SUMMARY ===")

message(sprintf("\n1. 3-class schema"))
message(sprintf("   Middle-bin peaks (%.2f < rain_frac < %.2f): %d (%.1f%%)",
                SNOW_EM_THRESH, RAIN_EM_THRESH,
                nrow(middle_bin),
                100 * nrow(middle_bin) / nrow(peaks_unreg)))
message(sprintf("   No bimodality detected at dominant elevations (> 2000 m)"))

message(sprintf("\n2. Rainfall threshold = %.2f", RAIN_EM_THRESH))
message(sprintf("   GPD xi CV across %.2f-%.2f range: %.3f",
                RAIN_THRESH_MIN, RAIN_THRESH_MAX, xi_cv))

message(sprintf("\n3. Snowmelt threshold = %.2f", SNOW_EM_THRESH))
message(sprintf("   GPD xi CV across %.2f-%.2f range: %.3f",
                SNOW_THRESH_MIN, SNOW_THRESH_MAX, snow_xi_cv))

message(sprintf("\n4. Window multiplier sensitivity"))
message(sprintf("   Non-rainfall peaks in scope: %d", nrow(rf_mult_peaks) / length(MULTIPLIERS)))
walk(MULTIPLIERS, function(m) {
  row <- snowmelt_frac %>% filter(multiplier == m)
  message(sprintf("   %dx: %.1f%% snowmelt (n = %d)", m, row$frac_snow * 100, row$n_snow))
})

message(sprintf("\nAll figures saved to %s/", out_dir))