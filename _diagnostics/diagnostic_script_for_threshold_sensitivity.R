#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Sensitivity of GPD Shape Parameters to Rain Fraction Thresholds
# ============================================================================
# Evaluates how sensitive the adaptive classifier's GPD shape parameter
# estimates are to the choice of RAIN_EM_THRESH and SNOW_EM_THRESH.
#
# Classification features (rf_short, rf_long, ros_flux_frac) are already
# computed in pot_peaks_classified_adaptive.rds — no need to re-derive them
# from accum_curves_60d.rds. This script re-runs only the decision tree and
# downstream GPD fitting across a two-dimensional threshold grid.
#
# For each (RAIN_EM_THRESH, SNOW_EM_THRESH) combination:
#   1. Re-classify peaks using the adaptive decision tree
#   2. Apply Yu et al. top-n peak selection (MIN_TYPE_SAMPLE = 20)
#   3. Fit GPD per flood type per site
#   4. Extract shape parameter (xi) by type
#
# Outputs:
#   _data/threshold_sensitivity.rds  — long-format results across grid
#   _diagnostics/threshold_sensitivity/  — figures
#
# INPUTS:
#   _data/pot_peaks_classified_adaptive.rds  (pre-computed features)
#   _data/site_metadata.rds
# ============================================================================

library(tidyverse)
library(lmomco)
library(patchwork)
library(furrr)

# ============================================================================
# CONFIGURATION
# ============================================================================
peaks_file  <- "_data/pot_peaks_classified_adaptive.rds"
meta_file   <- "_data/site_metadata.rds"
output_file <- "_data/threshold_sensitivity.rds"
out_dir     <- "_diagnostics/threshold_sensitivity"

# Fixed parameters (not varied)
ROS_COOCCUR_THRESH <- 0.50
MIN_SINGLE_SAMPLE  <- 25
MIN_TYPE_SAMPLE    <- 20
MIN_GPD_N          <- 10   # minimum peaks to attempt GPD fit

# Reference thresholds (current script 5b values)
RAIN_REF <- 0.95
SNOW_REF <- 0.05

# Threshold grid to evaluate
rain_thresholds <- seq(0.80, 0.99, by = 0.01)
snow_thresholds <- seq(0.01, 0.20, by = 0.01)

# Parallelism
N_WORKERS <- max(1L, parallelly::availableCores() - 2L)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# GPD FITTING (inlined from functions.R)
# ============================================================================
fit_gpd <- function(x, dates, mu_threshold, min_n = MIN_GPD_N) {
  keep  <- is.finite(x) & x > 0
  x     <- x[keep]
  dates <- dates[keep]
  if (length(x) < min_n) return(NULL)
  n_years <- as.numeric(diff(range(dates))) / 365.25
  if (n_years <= 0) return(NULL)

  mu <- mu_threshold
  y  <- sort(x - mu); y <- y[y >= 0]
  if (length(y) < min_n) return(NULL)

  pars <- tryCatch(
    pargpa(lmoms(y)),
    error = function(e) NULL
  )
  if (is.null(pars)) return(NULL)

  list(
    xi      = pars$para["kappa"] * -1,   # lmomco: kappa = -xi
    sigma   = pars$para["alpha"],
    mu      = mu,
    n       = length(x),
    n_years = n_years
  )
}

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading data...")
peaks <- readRDS(peaks_file)
meta  <- readRDS(meta_file)

# Unregulated sites only (consistent with snow window derivation)
unreg_sites <- meta %>% filter(!regulated) %>% pull(site_no)

# Base data: one row per peak, with all classification features
# Keep only data-sufficient peaks (this flag is invariant to threshold choice)
base <- peaks %>%
  filter(site_no %in% unreg_sites, data_sufficient) %>%
  select(site_no, peak_id, peak_date, peak_flow, pot_threshold_cfs,
         water_year, rf_short, rf_long, ros_flux_frac)

message(sprintf("Data-sufficient peaks from unregulated sites: %d across %d sites",
                nrow(base), n_distinct(base$site_no)))
message(sprintf("Threshold grid: %d rain × %d snow = %d combinations",
                length(rain_thresholds), length(snow_thresholds),
                length(rain_thresholds) * length(snow_thresholds)))

# ============================================================================
# WATER YEAR HELPER
# ============================================================================
water_year <- function(d) {
  yr <- as.integer(format(d, "%Y"))
  mo <- as.integer(format(d, "%m"))
  ifelse(mo >= 10L, yr + 1L, yr)
}

# ============================================================================
# FUNCTIONS: CLASSIFY -> SELECT -> FIT
# ============================================================================

# Apply the adaptive decision tree for a given threshold pair.
# Returns base with flood_type column added.
classify_peaks <- function(df, rain_thresh, snow_thresh) {
  df %>%
    mutate(
      flood_type = case_when(
        rf_short >= rain_thresh              ~ "rainfall",
        rf_long  <= snow_thresh              ~ "snowmelt",
        ros_flux_frac > ROS_COOCCUR_THRESH   ~ "ROS",
        TRUE                                 ~ "mixed"
      )
    )
}

# Apply Yu et al. top-n peak selection.
# Returns df with in_type_sample column; flood_type demoted to "unclassified"
# where type count < MIN_TYPE_SAMPLE.
select_peaks <- function(df) {
  type_sample_ids <- df %>%
    filter(flood_type != "unclassified") %>%
    group_by(site_no, flood_type, water_year) %>%
    slice_max(peak_flow, n = 1, with_ties = FALSE) %>%
    group_by(site_no, flood_type) %>%
    mutate(
      n    = pmax(MIN_TYPE_SAMPLE, n_distinct(water_year)),
      rank = row_number(desc(peak_flow))
    ) %>%
    filter(rank <= n, n() >= MIN_TYPE_SAMPLE) %>%
    ungroup() %>%
    select(site_no, peak_id, flood_type) %>%
    mutate(in_type_sample = TRUE)

  retained <- type_sample_ids %>% distinct(site_no, flood_type)

  df %>%
    left_join(type_sample_ids %>% select(site_no, peak_id, in_type_sample),
              by = c("site_no", "peak_id")) %>%
    mutate(
      in_type_sample = replace_na(in_type_sample, FALSE),
      flood_type = if_else(
        flood_type != "unclassified" &
          !paste(site_no, flood_type) %in%
            paste(retained$site_no, retained$flood_type),
        "unclassified",
        flood_type
      )
    )
}

# Fit GPD per site per flood type; return long-format xi estimates.
fit_all_gpd <- function(df) {
  df %>%
    filter(in_type_sample, flood_type != "unclassified") %>%
    group_by(site_no, flood_type) %>%
    group_modify(function(d, key) {
      fit <- fit_gpd(
        x            = d$peak_flow,
        dates        = d$peak_date,
        mu_threshold = first(d$pot_threshold_cfs)
      )
      if (is.null(fit)) return(tibble())
      tibble(xi = fit$xi, sigma = fit$sigma, n_peaks = fit$n)
    }) %>%
    ungroup()
}

# ============================================================================
# RUN GRID
# ============================================================================
grid <- expand_grid(
  rain_thresh = rain_thresholds,
  snow_thresh = snow_thresholds
)

message(sprintf("\nRunning sensitivity grid (%d combinations) with %d workers...",
                nrow(grid), N_WORKERS))

plan(multisession, workers = N_WORKERS)

results <- future_map2(
  grid$rain_thresh,
  grid$snow_thresh,
  function(rain_t, snow_t) {
    classified <- classify_peaks(base, rain_t, snow_t)
    selected   <- select_peaks(classified)
    fits       <- fit_all_gpd(selected)

    if (nrow(fits) == 0) return(tibble())

    fits %>%
      mutate(
        rain_thresh = rain_t,
        snow_thresh = snow_t
      )
  },
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

plan(sequential)

sensitivity <- bind_rows(results)

message(sprintf("\nGrid complete. Rows in results: %d", nrow(sensitivity)))

saveRDS(
  list(
    sensitivity    = sensitivity,
    grid           = grid,
    rain_ref       = RAIN_REF,
    snow_ref       = SNOW_REF,
    min_type_sample = MIN_TYPE_SAMPLE
  ),
  output_file
)
message(sprintf("Saved results to %s", output_file))

# ============================================================================
# SUMMARISE
# ============================================================================
# Median xi by type across sites, for each threshold combination
xi_summary <- sensitivity %>%
  group_by(rain_thresh, snow_thresh, flood_type) %>%
  summarise(
    median_xi = median(xi, na.rm = TRUE),
    q25_xi    = quantile(xi, 0.25, na.rm = TRUE),
    q75_xi    = quantile(xi, 0.75, na.rm = TRUE),
    n_sites   = n(),
    .groups   = "drop"
  )

# Count of sites retaining each type at each threshold combination
site_counts <- sensitivity %>%
  group_by(rain_thresh, snow_thresh, flood_type) %>%
  summarise(n_sites = n(), .groups = "drop")

# ============================================================================
# FIGURES
# ============================================================================
types      <- c("rainfall", "snowmelt", "ROS", "mixed")
type_labels <- c(
  rainfall = "Rainfall",
  snowmelt = "Snowmelt",
  ROS      = "Rain-on-snow",
  mixed    = "Mixed"
)

# ── Fig 1: Heatmaps of median xi by type ─────────────────────────────────────
# One facet per flood type; x = rain_thresh, y = snow_thresh,
# fill = median xi. Reference thresholds marked with crosshairs.

xi_heat <- xi_summary %>% filter(flood_type %in% types)

xi_range <- range(xi_heat$median_xi, na.rm = TRUE)
fill_lim <- c(floor(xi_range[1] * 10) / 10,
              ceiling(xi_range[2] * 10) / 10)

p_heat <- ggplot(xi_heat,
                 aes(x = rain_thresh, y = snow_thresh, fill = median_xi)) +
  geom_tile() +
  geom_vline(xintercept = RAIN_REF, color = "white",
             linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = SNOW_REF, color = "white",
             linetype = "dashed", linewidth = 0.6) +
  scale_fill_distiller(
    palette  = "RdBu",
    name     = "Median \u03be",
    limits   = fill_lim,
    oob      = scales::squish
  ) +
  scale_x_continuous(
    name   = "RAIN_EM_THRESH",
    breaks = seq(0.80, 0.99, by = 0.05)
  ) +
  scale_y_continuous(
    name   = "SNOW_EM_THRESH",
    breaks = seq(0.01, 0.20, by = 0.05)
  ) +
  facet_wrap(~ flood_type, ncol = 2,
             labeller = labeller(flood_type = type_labels)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right") +
  labs(
    title    = "Median GPD shape parameter (\u03be) across threshold grid",
    subtitle = sprintf(
      "White dashed lines: reference thresholds (RAIN = %.2f, SNOW = %.2f).",
      RAIN_REF, SNOW_REF
    ),
    caption = sprintf("Unregulated sites only. MIN_TYPE_SAMPLE = %d.",
                      MIN_TYPE_SAMPLE)
  )

ggsave(file.path(out_dir, "fig1_xi_heatmap.png"),
       p_heat, width = 10, height = 8, dpi = 150)
message("Saved fig1_xi_heatmap.png")

# ── Fig 2: Marginal slices through reference thresholds ──────────────────────
# Left panel: median xi vs. rain_thresh at SNOW_EM_THRESH = SNOW_REF
# Right panel: median xi vs. snow_thresh at RAIN_EM_THRESH = RAIN_REF

slice_rain <- xi_summary %>%
  filter(flood_type %in% types,
         abs(snow_thresh - SNOW_REF) < 1e-9)

slice_snow <- xi_summary %>%
  filter(flood_type %in% types,
         abs(rain_thresh - RAIN_REF) < 1e-9)

type_colors <- c(
  rainfall = "#E41A1C",
  snowmelt = "#377EB8",
  ROS      = "#4DAF4A",
  mixed    = "grey50"
)

p_rain_slice <- ggplot(slice_rain,
                       aes(x = rain_thresh, y = median_xi,
                           color = flood_type, fill = flood_type)) +
  geom_ribbon(aes(ymin = q25_xi, ymax = q75_xi),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = RAIN_REF, linetype = "dashed",
             color = "grey30", linewidth = 0.7) +
  scale_color_manual(values = type_colors, name = NULL,
                     labels = type_labels) +
  scale_fill_manual(values  = type_colors, name = NULL,
                    labels  = type_labels) +
  scale_x_continuous(name   = "RAIN_EM_THRESH",
                     breaks = seq(0.80, 0.99, by = 0.05)) +
  scale_y_continuous(name = "Median \u03be (IQR ribbon)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top") +
  labs(
    title    = sprintf("Varying RAIN_EM_THRESH (SNOW fixed at %.2f)", SNOW_REF),
    subtitle = "Median \u00b1 IQR of GPD \u03be across sites."
  )

p_snow_slice <- ggplot(slice_snow,
                       aes(x = snow_thresh, y = median_xi,
                           color = flood_type, fill = flood_type)) +
  geom_ribbon(aes(ymin = q25_xi, ymax = q75_xi),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = SNOW_REF, linetype = "dashed",
             color = "grey30", linewidth = 0.7) +
  scale_color_manual(values = type_colors, name = NULL,
                     labels = type_labels) +
  scale_fill_manual(values  = type_colors, name = NULL,
                    labels  = type_labels) +
  scale_x_continuous(name   = "SNOW_EM_THRESH",
                     breaks = seq(0.01, 0.20, by = 0.05)) +
  scale_y_continuous(name = "Median \u03be (IQR ribbon)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top") +
  labs(
    title    = sprintf("Varying SNOW_EM_THRESH (RAIN fixed at %.2f)", RAIN_REF),
    subtitle = "Median \u00b1 IQR of GPD \u03be across sites."
  )

p_slices <- p_rain_slice + p_snow_slice +
  plot_annotation(
    title   = "Marginal threshold sensitivity — GPD shape parameter",
    caption = sprintf("Unregulated sites only. MIN_TYPE_SAMPLE = %d. ROS_COOCCUR_THRESH fixed at %.2f.",
                      MIN_TYPE_SAMPLE, ROS_COOCCUR_THRESH)
  )

ggsave(file.path(out_dir, "fig2_xi_slices.png"),
       p_slices, width = 12, height = 5, dpi = 150)
message("Saved fig2_xi_slices.png")

# ── Fig 3: Site count retained by type across grid ───────────────────────────
# Heatmap of n_sites retaining each flood type; same layout as fig 1.
# Useful for understanding whether xi shifts are driven by threshold logic
# or by compositional changes (different sites entering/leaving each type).

p_nsites <- ggplot(site_counts %>% filter(flood_type %in% types),
                   aes(x = rain_thresh, y = snow_thresh, fill = n_sites)) +
  geom_tile() +
  geom_vline(xintercept = RAIN_REF, color = "white",
             linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = SNOW_REF, color = "white",
             linetype = "dashed", linewidth = 0.6) +
  scale_fill_distiller(
    palette  = "YlOrRd",
    direction = 1,
    name     = "Sites (n)"
  ) +
  scale_x_continuous(
    name   = "RAIN_EM_THRESH",
    breaks = seq(0.80, 0.99, by = 0.05)
  ) +
  scale_y_continuous(
    name   = "SNOW_EM_THRESH",
    breaks = seq(0.01, 0.20, by = 0.05)
  ) +
  facet_wrap(~ flood_type, ncol = 2,
             labeller = labeller(flood_type = type_labels)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right") +
  labs(
    title    = "Sites retaining each flood type across threshold grid",
    subtitle = sprintf(
      "White dashed lines: reference thresholds (RAIN = %.2f, SNOW = %.2f).",
      RAIN_REF, SNOW_REF
    ),
    caption  = sprintf("Unregulated sites only. MIN_TYPE_SAMPLE = %d.",
                       MIN_TYPE_SAMPLE)
  )

ggsave(file.path(out_dir, "fig3_site_counts.png"),
       p_nsites, width = 10, height = 8, dpi = 150)
message("Saved fig3_site_counts.png")

message("\nDone.")