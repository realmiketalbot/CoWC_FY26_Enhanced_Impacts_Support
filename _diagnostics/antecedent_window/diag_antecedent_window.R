#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Antecedent Window Harmonization
# ============================================================================
# Produces diagnostic outputs to inform a decision on whether to harmonize
# the antecedent window between depth_thresh (fixed 8-day, Davenport et al.
# 2020) and frac_thresh (adaptive search_radius, Livneh et al. equation
# 5 + ln(1.6092 × A)).
#
# Steps:
#   1. Inspect on-disk column schemas
#   2. Characterize search_radius distribution across sites
#   3. Assess rain_frac sensitivity to window length (W1/W2/W3)
#   4. Assess whether 8 days is a reasonable fixed choice
#   5. Document column schema discrepancy between classified outputs
#   6. Write plaintext summary report
#
# READ-ONLY: does not modify any pipeline scripts or _data/ files.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})

# ============================================================================
# PATHS
# ============================================================================
out_dir     <- "_diagnostics/antecedent_window"
peaks_file  <- "_data/pot_peaks_all.rds"
accum_file  <- "_data/accum_curves_60d.rds"
frac_file   <- "_data/pot_peaks_classified_frac_thresh_3class.rds"
depth_file  <- "_data/pot_peaks_classified_depth_thresh.rds"
meta_file   <- "_data/site_metadata.rds"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# CONSTANTS (must match Scripts 6a / 6b exactly)
# ============================================================================
NOISE_FLOOR_MM   <- 0.1    # mm/day
MIN_DAYS_SHORT   <- 0.90   # completeness fraction relative to window length
RAIN_FRAC_THRESH <- 0.99   # rainfall classification boundary
DEPTH_WIN        <- 8L     # Davenport et al. fixed window (days)

# ============================================================================
# STEP 1 — Inspect column schemas
# ============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("STEP 1: Inspect on-disk column schemas\n")
cat(strrep("=", 70), "\n\n", sep = "")

peaks <- readRDS(peaks_file)
cat("Columns in pot_peaks_all.rds:\n")
print(names(peaks))

if (!"search_radius" %in% names(peaks)) {
  stop("MISSING COLUMN: search_radius not found in pot_peaks_all.rds")
}
if (!"antecedent_days_used" %in% names(peaks)) {
  stop("MISSING COLUMN: antecedent_days_used not found in pot_peaks_all.rds")
}

cat("\nsearch_radius summary:\n")
print(summary(peaks$search_radius))
cat("\nantecedent_days_used summary:\n")
print(summary(peaks$antecedent_days_used))

accum <- readRDS(accum_file)
cat("\nColumns in accum_curves_60d.rds:\n")
print(names(accum))

if (!"lookback_day" %in% names(accum)) {
  stop("MISSING COLUMN: lookback_day not found in accum_curves_60d.rds")
}
cat("\nlookback_day range:\n")
print(range(accum$lookback_day, na.rm = TRUE))

frac  <- readRDS(frac_file)
cat("\nColumns in pot_peaks_classified_frac_thresh_3class.rds:\n")
print(names(frac))

depth <- readRDS(depth_file)
cat("\nColumns in pot_peaks_classified_depth_thresh.rds:\n")
print(names(depth))

# ============================================================================
# STEP 2 — Characterize search_radius distribution
# ============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("STEP 2: Characterize search_radius distribution\n")
cat(strrep("=", 70), "\n\n", sep = "")

# One search_radius per site — verify constant within site
peaks_dedup <- peaks %>%
  dplyr::distinct(site_no, peak_id, peak_date, peak_flow,
                  pot_threshold_cfs, search_radius,
                  antecedent_runoff_mm, antecedent_days_used)

site_sr <- peaks_dedup %>%
  dplyr::group_by(site_no) %>%
  dplyr::summarise(
    sr_min = min(search_radius, na.rm = TRUE),
    sr_max = max(search_radius, na.rm = TRUE),
    sr_val = dplyr::first(search_radius),
    .groups = "drop"
  )

n_inconstant <- sum(site_sr$sr_min != site_sr$sr_max)
if (n_inconstant > 0) {
  cat(sprintf(
    "WARNING: search_radius varies within %d site(s). Values below:\n",
    n_inconstant
  ))
  print(dplyr::filter(site_sr, sr_min != sr_max))
} else {
  cat("search_radius is constant within every site (as expected).\n")
}

sr_vals   <- site_sr$sr_val
sr_median <- median(sr_vals, na.rm = TRUE)
W3_window <- as.integer(round(sr_median))

cat(sprintf("\nCross-site search_radius summary:\n"))
cat(sprintf("  Min    : %g\n",    min(sr_vals,  na.rm = TRUE)))
cat(sprintf("  Median : %g\n",    sr_median))
cat(sprintf("  Mean   : %.2f\n",  mean(sr_vals, na.rm = TRUE)))
cat(sprintf("  Max    : %g\n",    max(sr_vals,  na.rm = TRUE)))
cat(sprintf("  Frac <= 8 : %.3f\n", mean(sr_vals <= 8, na.rm = TRUE)))
cat(sprintf("  Frac == 8 : %.3f\n", mean(sr_vals == 8, na.rm = TRUE)))
cat(sprintf("  Frac >  8 : %.3f\n", mean(sr_vals  > 8, na.rm = TRUE)))
cat(sprintf("  W3 (rounded median) : %d days\n\n", W3_window))

# Histogram
p_hist <- ggplot2::ggplot(site_sr, ggplot2::aes(x = sr_val)) +
  ggplot2::geom_histogram(binwidth = 1, fill = "#4472C4", color = "white") +
  ggplot2::geom_vline(xintercept = 8, linetype = "dashed",
                      color = "firebrick", linewidth = 0.9) +
  ggplot2::annotate("text", x = 8.3, y = Inf,
                    label = "8 days\n(depth_thresh)", hjust = 0, vjust = 1.3,
                    color = "firebrick", size = 3.2) +
  ggplot2::labs(
    title    = "Distribution of search_radius across sites",
    subtitle = sprintf("n = %d sites  |  median = %g d  |  W3 (rounded median) = %d d",
                       nrow(site_sr), sr_median, W3_window),
    x        = "search_radius (days)",
    y        = "Number of sites"
  ) +
  ggplot2::theme_bw(base_size = 11)

fig_hist <- file.path(out_dir, "fig_search_radius_distribution.png")
ggplot2::ggsave(fig_hist, p_hist, width = 8, height = 5, dpi = 300)
cat(sprintf("Saved: %s\n", fig_hist))

# Correlation with drainage area: check pot_peaks_all.rds then site_metadata.rds
drain_available <- FALSE

if ("drain_sqkm" %in% names(peaks_dedup)) {
  drain_df <- peaks_dedup %>%
    dplyr::distinct(site_no, drain_sqkm) %>%
    dplyr::left_join(site_sr %>% dplyr::select(site_no, sr_val), by = "site_no") %>%
    dplyr::filter(!is.na(drain_sqkm))
  drain_available <- nrow(drain_df) > 0
} else {
  meta <- readRDS(meta_file)
  if ("drain_sqkm" %in% names(meta)) {
    drain_df <- meta %>%
      dplyr::select(site_no, drain_sqkm) %>%
      dplyr::left_join(site_sr %>% dplyr::select(site_no, sr_val), by = "site_no") %>%
      dplyr::filter(!is.na(drain_sqkm))
    drain_available <- nrow(drain_df) > 0
    cat("drain_sqkm sourced from site_metadata.rds\n")
  } else {
    cat("WARNING: drain_sqkm not found in pot_peaks_all.rds or site_metadata.rds.",
        "Skipping drainage area correlation plot.\n")
  }
}

if (drain_available) {
  drain_df <- drain_df %>%
    dplyr::mutate(
      log_drain   = log(drain_sqkm),
      sr_expected = 5 + log(1.6092 * drain_sqkm)
    ) %>%
    dplyr::arrange(log_drain)

  p_drain <- ggplot2::ggplot(drain_df,
                              ggplot2::aes(x = log_drain, y = sr_val)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.5, color = "#4472C4") +
    ggplot2::geom_line(ggplot2::aes(y = sr_expected),
                       color = "firebrick", linewidth = 0.8, linetype = "solid") +
    ggplot2::annotate("text", x = -Inf, y = Inf,
                      label = "Red line: expected 5 + ln(1.6092 × A)",
                      hjust = -0.05, vjust = 1.4, color = "firebrick", size = 3) +
    ggplot2::labs(
      title    = "search_radius vs. log(drain_sqkm)",
      subtitle = "Points = observed (integer ceiling); red = continuous Livneh formula",
      x        = "log(drainage area, km²)",
      y        = "search_radius (days)"
    ) +
    ggplot2::theme_bw(base_size = 11)

  fig_drain <- file.path(out_dir, "fig_search_radius_vs_drain_area.png")
  ggplot2::ggsave(fig_drain, p_drain, width = 8, height = 5, dpi = 300)
  cat(sprintf("Saved: %s\n", fig_drain))
}

# ============================================================================
# STEP 3 — rain_frac sensitivity to window length
# ============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("STEP 3: Assess rain_frac sensitivity to window length\n")
cat(strrep("=", 70), "\n\n", sep = "")

# Build per-peak window lookup (one row per peak)
peak_windows <- peaks_dedup %>%
  dplyr::distinct(site_no, peak_id, search_radius)

# Helper: compute rain_frac for a given window-length column
compute_rain_frac <- function(accum_df, peak_windows_df, window_col, w3_val) {
  # Attach per-peak window length (or fixed value) to accum rows, then compute
  wl <- peak_windows_df %>%
    dplyr::mutate(win = if (window_col == "search_radius") search_radius
                       else if (window_col == "fixed_8")   DEPTH_WIN
                       else                                 w3_val)

  accum_df %>%
    dplyr::inner_join(wl %>% dplyr::select(site_no, peak_id, win),
                      by = c("site_no", "peak_id")) %>%
    dplyr::filter(lookback_day <= win) %>%
    dplyr::mutate(
      rain_active  = dplyr::if_else(rain_mm_day >= NOISE_FLOOR_MM, rain_mm_day, 0),
      melt_active  = dplyr::if_else(melt_mm_day >= NOISE_FLOOR_MM, melt_mm_day, 0),
      total_active = rain_active + melt_active
    ) %>%
    dplyr::group_by(site_no, peak_id) %>%
    dplyr::summarise(
      total_rain_flux = sum(rain_active,  na.rm = TRUE),
      total_melt_flux = sum(melt_active,  na.rm = TRUE),
      total_flux      = sum(total_active, na.rm = TRUE),
      n_days_data     = dplyr::first(n_days_data),
      window_used     = dplyr::first(win),
      .groups         = "drop"
    ) %>%
    dplyr::mutate(
      rain_frac = dplyr::if_else(total_flux > 0,
                                 total_rain_flux / total_flux,
                                 NA_real_),
      data_sufficient = n_days_data >= ceiling(window_used * MIN_DAYS_SHORT)
    )
}

cat("Computing W1 (fixed 8-day)...\n")
rf_W1 <- compute_rain_frac(accum, peak_windows, "fixed_8",      W3_window) %>%
  dplyr::rename(rain_frac_W1 = rain_frac,
                n_days_W1    = n_days_data,
                suff_W1      = data_sufficient) %>%
  dplyr::select(site_no, peak_id, rain_frac_W1, n_days_W1, suff_W1)

cat("Computing W2 (adaptive search_radius)...\n")
rf_W2 <- compute_rain_frac(accum, peak_windows, "search_radius", W3_window) %>%
  dplyr::rename(rain_frac_W2 = rain_frac,
                n_days_W2    = n_days_data,
                suff_W2      = data_sufficient) %>%
  dplyr::select(site_no, peak_id, rain_frac_W2, n_days_W2, suff_W2,
                window_used)

cat(sprintf("Computing W3 (fixed %d-day, rounded median)...\n", W3_window))
rf_W3 <- compute_rain_frac(accum, peak_windows, "fixed_W3",    W3_window) %>%
  dplyr::rename(rain_frac_W3 = rain_frac,
                n_days_W3    = n_days_data,
                suff_W3      = data_sufficient) %>%
  dplyr::select(site_no, peak_id, rain_frac_W3, n_days_W3, suff_W3)

# Pull flood_type from existing frac_thresh output for colouring plots
frac_type <- frac %>%
  dplyr::distinct(site_no, peak_id, flood_type, regulated)

rf_all <- rf_W1 %>%
  dplyr::inner_join(rf_W2, by = c("site_no", "peak_id")) %>%
  dplyr::inner_join(rf_W3, by = c("site_no", "peak_id")) %>%
  dplyr::left_join(frac_type, by = c("site_no", "peak_id")) %>%
  dplyr::filter(suff_W1, suff_W2, suff_W3)

cat(sprintf("Peaks with all three windows data-sufficient: %d\n\n", nrow(rf_all)))

# Retrieve MIXED_SPLIT from stored attribute
MIXED_SPLIT <- attr(frac, "mixed_split")
cat(sprintf("MIXED_SPLIT (from attr): %.4f\n", MIXED_SPLIT))

# --- Scatter: W1 vs W2 ---
p_w1w2 <- ggplot2::ggplot(
    rf_all %>% dplyr::filter(!is.na(rain_frac_W1), !is.na(rain_frac_W2)),
    ggplot2::aes(x = rain_frac_W1, y = rain_frac_W2, color = flood_type)
  ) +
  ggplot2::geom_point(alpha = 0.35, size = 0.8) +
  ggplot2::geom_abline(slope = 1, intercept = 0,
                       color = "black", linetype = "dashed", linewidth = 0.7) +
  ggplot2::geom_hline(yintercept = RAIN_FRAC_THRESH,
                      linetype = "dotted", color = "grey40") +
  ggplot2::geom_hline(yintercept = MIXED_SPLIT,
                      linetype = "dotted", color = "grey40") +
  ggplot2::geom_vline(xintercept = RAIN_FRAC_THRESH,
                      linetype = "dotted", color = "grey40") +
  ggplot2::geom_vline(xintercept = MIXED_SPLIT,
                      linetype = "dotted", color = "grey40") +
  ggplot2::scale_color_manual(
    values = c(rainfall = "#2166AC", snowmelt = "#92C5DE",
               rain_on_snow = "#F4A582", unclassified = "grey70"),
    na.value = "grey70"
  ) +
  ggplot2::labs(
    title    = "rain_frac: W1 (fixed 8d) vs. W2 (adaptive search_radius)",
    subtitle = sprintf("Data-sufficient peaks (n = %d)  |  dotted lines = RAIN_FRAC_THRESH & MIXED_SPLIT",
                       nrow(dplyr::filter(rf_all, !is.na(rain_frac_W1), !is.na(rain_frac_W2)))),
    x        = "rain_frac  [W1, 8 days]",
    y        = "rain_frac  [W2, adaptive]",
    color    = "flood_type (W2)"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1))

fig_w1w2 <- file.path(out_dir, "fig_rainfrac_W1_vs_W2.png")
ggplot2::ggsave(fig_w1w2, p_w1w2, width = 8, height = 5, dpi = 300)
cat(sprintf("Saved: %s\n", fig_w1w2))

# --- Scatter: W1 vs W3 ---
p_w1w3 <- ggplot2::ggplot(
    rf_all %>% dplyr::filter(!is.na(rain_frac_W1), !is.na(rain_frac_W3)),
    ggplot2::aes(x = rain_frac_W1, y = rain_frac_W3, color = flood_type)
  ) +
  ggplot2::geom_point(alpha = 0.35, size = 0.8) +
  ggplot2::geom_abline(slope = 1, intercept = 0,
                       color = "black", linetype = "dashed", linewidth = 0.7) +
  ggplot2::geom_hline(yintercept = RAIN_FRAC_THRESH,
                      linetype = "dotted", color = "grey40") +
  ggplot2::geom_hline(yintercept = MIXED_SPLIT,
                      linetype = "dotted", color = "grey40") +
  ggplot2::geom_vline(xintercept = RAIN_FRAC_THRESH,
                      linetype = "dotted", color = "grey40") +
  ggplot2::geom_vline(xintercept = MIXED_SPLIT,
                      linetype = "dotted", color = "grey40") +
  ggplot2::scale_color_manual(
    values = c(rainfall = "#2166AC", snowmelt = "#92C5DE",
               rain_on_snow = "#F4A582", unclassified = "grey70"),
    na.value = "grey70"
  ) +
  ggplot2::labs(
    title    = sprintf("rain_frac: W1 (fixed 8d) vs. W3 (fixed %dd, rounded median)", W3_window),
    subtitle = sprintf("Data-sufficient peaks (n = %d)  |  dotted lines = RAIN_FRAC_THRESH & MIXED_SPLIT",
                       nrow(dplyr::filter(rf_all, !is.na(rain_frac_W1), !is.na(rain_frac_W3)))),
    x        = "rain_frac  [W1, 8 days]",
    y        = sprintf("rain_frac  [W3, %d days]", W3_window),
    color    = "flood_type (W2)"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1))

fig_w1w3 <- file.path(out_dir, "fig_rainfrac_W1_vs_W3.png")
ggplot2::ggsave(fig_w1w3, p_w1w3, width = 8, height = 5, dpi = 300)
cat(sprintf("Saved: %s\n", fig_w1w3))

# --- Classification change table: W1 vs W2 ---
classify_rf <- function(rf, mixed_split) {
  dplyr::case_when(
    is.na(rf)            ~ "unclassified",
    rf >= RAIN_FRAC_THRESH ~ "rainfall",
    rf <  mixed_split      ~ "snowmelt",
    TRUE                   ~ "rain_on_snow"
  )
}

change_tbl <- rf_all %>%
  dplyr::filter(!is.na(rain_frac_W1), !is.na(rain_frac_W2)) %>%
  dplyr::mutate(
    type_W1 = classify_rf(rain_frac_W1, MIXED_SPLIT),
    type_W2 = classify_rf(rain_frac_W2, MIXED_SPLIT),
    changes  = type_W1 != type_W2
  ) %>%
  dplyr::group_by(flood_type, regulated) %>%
  dplyr::summarise(
    n_peaks         = dplyr::n(),
    n_changes       = sum(changes),
    frac_changes    = mean(changes),
    .groups         = "drop"
  ) %>%
  dplyr::arrange(flood_type, regulated)

cat("\nClassification change table (W1 → W2, by flood_type and regulation):\n")
print(change_tbl, n = Inf)

overall_frac_changes <- rf_all %>%
  dplyr::filter(!is.na(rain_frac_W1), !is.na(rain_frac_W2)) %>%
  dplyr::mutate(
    type_W1 = classify_rf(rain_frac_W1, MIXED_SPLIT),
    type_W2 = classify_rf(rain_frac_W2, MIXED_SPLIT)
  ) %>%
  dplyr::summarise(
    n_peaks      = dplyr::n(),
    n_changes    = sum(type_W1 != type_W2),
    frac_changes = mean(type_W1 != type_W2)
  )

cat(sprintf("\nOverall fraction of peaks changing class (W1 → W2): %.4f  (%d / %d)\n",
            overall_frac_changes$frac_changes,
            overall_frac_changes$n_changes,
            overall_frac_changes$n_peaks))

# --- Density plot: W1 vs W2 rain_frac distributions ---
density_df <- rf_all %>%
  dplyr::filter(!is.na(rain_frac_W1) | !is.na(rain_frac_W2)) %>%
  tidyr::pivot_longer(
    cols      = c(rain_frac_W1, rain_frac_W2),
    names_to  = "window",
    values_to = "rain_frac"
  ) %>%
  dplyr::filter(!is.na(rain_frac)) %>%
  dplyr::mutate(
    window = dplyr::recode(window,
                           rain_frac_W1 = "W1 (8 days, fixed)",
                           rain_frac_W2 = "W2 (adaptive)")
  )

p_density <- ggplot2::ggplot(density_df,
                              ggplot2::aes(x = rain_frac, fill = window, color = window)) +
  ggplot2::geom_density(alpha = 0.35, linewidth = 0.7) +
  ggplot2::geom_vline(xintercept = RAIN_FRAC_THRESH,
                      linetype = "dashed", color = "black", linewidth = 0.7) +
  ggplot2::geom_vline(xintercept = MIXED_SPLIT,
                      linetype = "dotted", color = "black", linewidth = 0.7) +
  ggplot2::annotate("text", x = RAIN_FRAC_THRESH - 0.01, y = Inf,
                    label = "RAIN_FRAC_THRESH", hjust = 1, vjust = 1.4,
                    size = 3) +
  ggplot2::annotate("text", x = MIXED_SPLIT - 0.01, y = Inf,
                    label = "MIXED_SPLIT", hjust = 1, vjust = 3,
                    size = 3) +
  ggplot2::scale_fill_manual(values  = c("W1 (8 days, fixed)" = "#2166AC",
                                         "W2 (adaptive)"      = "#D6604D")) +
  ggplot2::scale_color_manual(values = c("W1 (8 days, fixed)" = "#2166AC",
                                         "W2 (adaptive)"      = "#D6604D")) +
  ggplot2::labs(
    title = "rain_frac distribution: W1 (fixed 8d) vs. W2 (adaptive)",
    x     = "rain_frac",
    y     = "Density",
    fill  = "Window",
    color = "Window"
  ) +
  ggplot2::theme_bw(base_size = 11)

fig_density <- file.path(out_dir, "fig_rainfrac_density_W1_vs_W2.png")
ggplot2::ggsave(fig_density, p_density, width = 8, height = 5, dpi = 300)
cat(sprintf("Saved: %s\n", fig_density))

# ============================================================================
# STEP 4 — Sensitivity of rain_frac difference vs. search_radius magnitude
# ============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("STEP 4: rain_frac sensitivity vs. window discrepancy magnitude\n")
cat(strrep("=", 70), "\n\n", sep = "")

sensitivity_df <- rf_all %>%
  dplyr::filter(!is.na(rain_frac_W1), !is.na(rain_frac_W2)) %>%
  dplyr::left_join(
    peak_windows %>% dplyr::select(site_no, peak_id, search_radius),
    by = c("site_no", "peak_id")
  ) %>%
  dplyr::mutate(abs_diff = abs(rain_frac_W1 - rain_frac_W2)) %>%
  dplyr::group_by(search_radius) %>%
  dplyr::summarise(
    n_peaks       = dplyr::n(),
    mean_abs_diff = mean(abs_diff, na.rm = TRUE),
    med_abs_diff  = median(abs_diff, na.rm = TRUE),
    .groups       = "drop"
  )

p_sens <- ggplot2::ggplot(sensitivity_df,
                           ggplot2::aes(x = search_radius, y = mean_abs_diff)) +
  ggplot2::geom_point(ggplot2::aes(size = n_peaks), alpha = 0.7, color = "#4472C4") +
  ggplot2::geom_smooth(method = "loess", se = TRUE, color = "firebrick",
                        fill = "firebrick", alpha = 0.15, linewidth = 0.8) +
  ggplot2::geom_vline(xintercept = 8, linetype = "dashed",
                      color = "grey40", linewidth = 0.7) +
  ggplot2::annotate("text", x = 8.3, y = Inf,
                    label = "8 days", hjust = 0, vjust = 1.3,
                    size = 3, color = "grey40") +
  ggplot2::scale_size_continuous(range = c(1, 6), name = "n peaks") +
  ggplot2::labs(
    title    = "Mean |rain_frac_W1 − rain_frac_W2| vs. search_radius",
    subtitle = "Larger points = more peaks; loess smoother with 95% CI",
    x        = "search_radius (days)",
    y        = "Mean absolute rain_frac difference (W1 − W2)"
  ) +
  ggplot2::theme_bw(base_size = 11)

fig_sens <- file.path(out_dir, "fig_rainfrac_sensitivity_vs_window.png")
ggplot2::ggsave(fig_sens, p_sens, width = 8, height = 5, dpi = 300)
cat(sprintf("Saved: %s\n", fig_sens))

# ============================================================================
# STEP 5 — Column schema discrepancy report
# ============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("STEP 5: Column schema discrepancy — depth_thresh vs. frac_thresh_3class\n")
cat(strrep("=", 70), "\n\n", sep = "")

cols_depth <- names(depth)
cols_frac  <- names(frac)
cols_both  <- intersect(cols_depth, cols_frac)
cols_depth_only <- setdiff(cols_depth, cols_frac)
cols_frac_only  <- setdiff(cols_frac,  cols_depth)

# Antecedent-window-relevant columns of interest
aw_keywords <- c("n_days", "complete", "search_radius", "data_sufficient",
                 "antecedent", "window", "short_window", "wldas")
is_aw <- function(x) any(sapply(aw_keywords, function(k) grepl(k, x, ignore.case = TRUE)))

depth_aw <- cols_depth[sapply(cols_depth, is_aw)]
frac_aw  <- cols_frac[ sapply(cols_frac,  is_aw)]

cat("Antecedent-window columns in depth_thresh:\n")
cat(paste("  ", depth_aw, collapse = "\n"), "\n\n")

cat("Antecedent-window columns in frac_thresh_3class:\n")
cat(paste("  ", frac_aw, collapse = "\n"), "\n\n")

cat("Columns in depth_thresh only (absent from frac_thresh_3class):\n")
cat(paste("  ", cols_depth_only, collapse = "\n"), "\n\n")

cat("Columns in frac_thresh_3class only (absent from depth_thresh):\n")
cat(paste("  ", cols_frac_only, collapse = "\n"), "\n\n")

cat("Columns present in both:\n")
cat(paste("  ", cols_both, collapse = "\n"), "\n\n")

# ============================================================================
# STEP 6 — Summary report
# ============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("STEP 6: Writing summary report\n")
cat(strrep("=", 70), "\n\n", sep = "")

summary_path <- file.path(out_dir, "summary.txt")

sr_stats <- list(
  min    = min(sr_vals,  na.rm = TRUE),
  median = median(sr_vals, na.rm = TRUE),
  mean   = mean(sr_vals, na.rm = TRUE),
  max    = max(sr_vals,  na.rm = TRUE),
  frac_lt_8 = mean(sr_vals < 8, na.rm = TRUE),
  frac_eq_8 = mean(sr_vals == 8, na.rm = TRUE),
  frac_gt_8 = mean(sr_vals > 8, na.rm = TRUE),
  n_sites   = length(sr_vals)
)

writeLines(
  c(
    "Antecedent Window Diagnostic — Summary Report",
    sprintf("Generated: %s", Sys.time()),
    "",
    strrep("-", 60),
    "SEARCH_RADIUS DISTRIBUTION (one value per site)",
    strrep("-", 60),
    sprintf("  n sites         : %d", sr_stats$n_sites),
    sprintf("  min             : %g days", sr_stats$min),
    sprintf("  median          : %g days", sr_stats$median),
    sprintf("  mean            : %.2f days", sr_stats$mean),
    sprintf("  max             : %g days", sr_stats$max),
    sprintf("  frac < 8 days   : %.4f", sr_stats$frac_lt_8),
    sprintf("  frac == 8 days  : %.4f", sr_stats$frac_eq_8),
    sprintf("  frac > 8 days   : %.4f", sr_stats$frac_gt_8),
    sprintf("  W3 (rounded median) : %d days", W3_window),
    "",
    strrep("-", 60),
    "MIXED_SPLIT (stored attribute on frac_thresh_3class output)",
    strrep("-", 60),
    sprintf("  MIXED_SPLIT     : %.4f", MIXED_SPLIT),
    sprintf("  RAIN_FRAC_THRESH: %.2f",  RAIN_FRAC_THRESH),
    "",
    strrep("-", 60),
    "CLASSIFICATION CHANGE: W1 (fixed 8d) substituted for W2 (adaptive)",
    strrep("-", 60),
    sprintf("  Overall: %.4f of peaks change class  (%d / %d)",
            overall_frac_changes$frac_changes,
            overall_frac_changes$n_changes,
            overall_frac_changes$n_peaks),
    "",
    "  By flood_type and regulation:",
    paste(
      capture.output(print(change_tbl, n = Inf)),
      collapse = "\n"
    ),
    "",
    strrep("-", 60),
    "COLUMN SCHEMA DISCREPANCIES",
    strrep("-", 60),
    "",
    "Antecedent-window columns in depth_thresh:",
    paste("  ", depth_aw, collapse = "\n"),
    "",
    "Antecedent-window columns in frac_thresh_3class:",
    paste("  ", frac_aw, collapse = "\n"),
    "",
    "Columns in depth_thresh ONLY:",
    paste("  ", if (length(cols_depth_only) > 0) cols_depth_only else "(none)",
          collapse = "\n"),
    "",
    "Columns in frac_thresh_3class ONLY:",
    paste("  ", if (length(cols_frac_only) > 0) cols_frac_only else "(none)",
          collapse = "\n"),
    "",
    strrep("-", 60),
    "NOTE: No methodological recommendation is made in this file.",
    strrep("-", 60)
  ),
  con = summary_path
)

cat(sprintf("Saved: %s\n", summary_path))

# ============================================================================
# COMPLETION SUMMARY
# ============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("DIAGNOSTIC COMPLETE\n")
cat(strrep("=", 70), "\n\n", sep = "")

written_files <- c(
  fig_hist,
  fig_w1w2,
  fig_w1w3,
  fig_density,
  fig_sens,
  summary_path
)
if (drain_available) written_files <- c(fig_drain, written_files)

cat("Files written:\n")
purrr::walk(written_files, ~ cat(sprintf("  %s\n", .x)))
