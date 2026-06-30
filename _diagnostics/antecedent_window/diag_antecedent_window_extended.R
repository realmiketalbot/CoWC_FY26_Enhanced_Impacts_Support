#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Antecedent Window — Extended (Operational Impact of W1)
# ============================================================================
# Extends diag_antecedent_window.R with figures focused on the downstream
# consequences of substituting a fixed 8-day window (W1) for the adaptive
# search_radius window (W2) in frac_thresh classification.
#
# The prior diagnostic established that 6.7% of peaks change class under W1.
# These figures answer: does that matter for NA rates, seasonal coverage,
# and mixture-FFA site-type viability?
#
# Figures produced:
#   fig_na_rainfrac_rate.png               — Step 1: NA rain_frac by type/month
#   fig_reclassified_unclassified_by_season.png — Step 2: W2→unclass by season
#   fig_site_type_viability.png            — Step 3: per-site type viability
#   fig_cumulative_flux_profiles.png       — Step 4: cumul flux profiles
#
# INPUTS (read-only):
#   _data/accum_curves_60d.rds
#   _data/pot_peaks_all.rds
#   _data/pot_peaks_classified_frac_thresh_3class.rds
#   _diagnostics/antecedent_window/summary.txt  (appended at end)
#
# Run from project root. Requires: tidyverse, patchwork, scales.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(patchwork)
  library(scales)
})

# ============================================================================
# CONFIGURATION
# ============================================================================
out_dir    <- "_diagnostics/antecedent_window"
peaks_file <- "_data/pot_peaks_all.rds"
accum_file <- "_data/accum_curves_60d.rds"
frac_file  <- "_data/pot_peaks_classified_frac_thresh_3class.rds"
summary_file <- file.path(out_dir, "summary.txt")

NOISE_FLOOR_MM   <- 0.1
MIN_DAYS_SHORT   <- 0.90
RAIN_FRAC_THRESH <- 0.99
MIN_TYPE_SAMPLE  <- 20L
DEPTH_WIN        <- 8L

TYPE_LEVELS  <- c("rainfall", "snowmelt", "rain_on_snow", "unclassified")
TYPE_COLORS  <- c(rainfall     = "#2166AC",
                  snowmelt     = "#92C5DE",
                  rain_on_snow = "#F4A582",
                  unclassified = "grey70")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading data...")

accum <- readRDS(accum_file)

peaks <- readRDS(peaks_file) %>%
  dplyr::distinct(site_no, peak_id, peak_date, peak_flow, search_radius) %>%
  dplyr::mutate(
    peak_date  = as.Date(peak_date),
    water_year = {
      yr <- as.integer(format(peak_date, "%Y"))
      mo <- as.integer(format(peak_date, "%m"))
      ifelse(mo >= 10L, yr + 1L, yr)
    },
    peak_month = as.integer(format(peak_date, "%m"))
  )

frac <- readRDS(frac_file)
MIXED_SPLIT <- attr(frac, "mixed_split")
stopifnot(!is.na(MIXED_SPLIT))
message(sprintf("MIXED_SPLIT: %.4f", MIXED_SPLIT))

# W3 = rounded median search_radius across sites (for Fig 4 reference line)
W3_window <- peaks %>%
  dplyr::distinct(site_no, search_radius) %>%
  dplyr::pull(search_radius) %>%
  median(na.rm = TRUE) %>%
  round() %>%
  as.integer()
message(sprintf("W3 window (rounded median search_radius): %d days", W3_window))

# Final W2 flood_type from existing classification output (includes demotion)
frac_type <- frac %>%
  dplyr::distinct(site_no, peak_id, flood_type, regulated) %>%
  dplyr::mutate(flood_type = factor(as.character(flood_type), levels = TYPE_LEVELS))

message(sprintf("Peaks: %d across %d sites  |  Accum rows: %d",
                nrow(peaks), dplyr::n_distinct(peaks$site_no), nrow(accum)))

# ============================================================================
# RECOMPUTE rain_frac UNDER W1 AND W2
# ============================================================================
# Uses the same active-flux formula as Script 6b.
# win_val: integer for a fixed window, NA_integer_ to use per-peak search_radius.
compute_rf <- function(accum_df, peak_windows_df, win_val) {
  if (is.na(win_val)) {
    wl <- peak_windows_df %>% dplyr::mutate(win = search_radius)
  } else {
    wl <- peak_windows_df %>% dplyr::mutate(win = win_val)
  }

  accum_df %>%
    dplyr::inner_join(
      wl %>% dplyr::select(site_no, peak_id, win),
      by = c("site_no", "peak_id")
    ) %>%
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
      na_rain_frac    = (total_flux == 0),
      rain_frac       = dplyr::if_else(total_flux > 0,
                                       total_rain_flux / total_flux,
                                       NA_real_),
      data_sufficient = n_days_data >= ceiling(window_used * MIN_DAYS_SHORT)
    )
}

classify_rf <- function(rain_frac, data_sufficient, mixed_split) {
  dplyr::case_when(
    !data_sufficient | is.na(rain_frac) ~ "unclassified",
    rain_frac >= RAIN_FRAC_THRESH       ~ "rainfall",
    rain_frac <  mixed_split            ~ "snowmelt",
    TRUE                                ~ "rain_on_snow"
  )
}

peak_windows <- peaks %>% dplyr::select(site_no, peak_id, search_radius)

message("Computing W1 (fixed 8-day)...")
rf_W1 <- compute_rf(accum, peak_windows, DEPTH_WIN)

message("Computing W2 (adaptive search_radius)...")
rf_W2 <- compute_rf(accum, peak_windows, NA_integer_)

# Combined per-peak table
peaks_all <- peaks %>%
  dplyr::left_join(
    rf_W1 %>% dplyr::select(site_no, peak_id,
                              rain_frac_W1 = rain_frac,
                              suff_W1      = data_sufficient,
                              na_rf_W1     = na_rain_frac),
    by = c("site_no", "peak_id")
  ) %>%
  dplyr::left_join(
    rf_W2 %>% dplyr::select(site_no, peak_id,
                              rain_frac_W2 = rain_frac,
                              suff_W2      = data_sufficient,
                              na_rf_W2     = na_rain_frac),
    by = c("site_no", "peak_id")
  ) %>%
  dplyr::left_join(frac_type, by = c("site_no", "peak_id")) %>%
  dplyr::mutate(
    # Re-derived classifications without sample-size demotion (for change logic)
    type_W1_raw = classify_rf(rain_frac_W1, suff_W1, MIXED_SPLIT),
    type_W2_raw = classify_rf(rain_frac_W2, suff_W2, MIXED_SPLIT),
    # flood_type (from frac) = final W2 classification including demotion
    flood_type_chr = as.character(flood_type)
  )

message(sprintf("Peaks in combined table: %d", nrow(peaks_all)))

# ============================================================================
# FIGURE 1 — NA rain_frac rate: W1 vs. W2, by flood type and month
# ============================================================================
message("\nFigure 1: NA rain_frac rates...")

# Panel A: grouped bar by W2 flood_type (final, from frac)
fig1a_df <- peaks_all %>%
  dplyr::group_by(flood_type) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_na_W1 = sum(na_rf_W1, na.rm = TRUE),
    n_na_W2 = sum(na_rf_W2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    frac_W1 = n_na_W1 / n_total,
    frac_W2 = n_na_W2 / n_total
  ) %>%
  tidyr::pivot_longer(
    cols      = c(frac_W1, frac_W2),
    names_to  = "window",
    values_to = "frac_na"
  ) %>%
  dplyr::mutate(
    n_na   = dplyr::if_else(window == "frac_W1", n_na_W1, n_na_W2),
    window = dplyr::recode(window,
                            frac_W1 = "W1 (8d fixed)",
                            frac_W2 = "W2 (adaptive)")
  )

p1a <- ggplot2::ggplot(
    fig1a_df,
    ggplot2::aes(x = flood_type, y = frac_na, fill = window)
  ) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                    width = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = n_na, y = frac_na + 0.0015),
    position = ggplot2::position_dodge(width = 0.8),
    vjust = 0, size = 2.8
  ) +
  ggplot2::scale_fill_manual(values = c("W1 (8d fixed)" = "#2166AC",
                                         "W2 (adaptive)" = "#D6604D")) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1),
    expand = ggplot2::expansion(mult = c(0, 0.12))
  ) +
  ggplot2::labs(
    title = "A: NA rain_frac rate by W2 flood type",
    x     = "W2 flood type (final, from frac_thresh output)",
    y     = "Fraction of peaks with NA rain_frac",
    fill  = "Window"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 15, hjust = 1))

# Panel B: line by month (1–12)
fig1b_df <- peaks_all %>%
  dplyr::group_by(peak_month) %>%
  dplyr::summarise(
    frac_W1 = mean(na_rf_W1, na.rm = TRUE),
    frac_W2 = mean(na_rf_W2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols      = c(frac_W1, frac_W2),
    names_to  = "window",
    values_to = "frac_na"
  ) %>%
  dplyr::mutate(
    window = dplyr::recode(window,
                            frac_W1 = "W1 (8d fixed)",
                            frac_W2 = "W2 (adaptive)")
  )

p1b <- ggplot2::ggplot(fig1b_df,
                        ggplot2::aes(x = peak_month, y = frac_na, color = window)) +
  ggplot2::annotate("rect", xmin = 2.5, xmax = 6.5,
                    ymin = -Inf, ymax = Inf,
                    fill = "steelblue", alpha = 0.08) +
  ggplot2::annotate("text", x = 4.5, y = Inf,
                    label = "snowmelt\nseason", vjust = 1.3,
                    size = 3, color = "steelblue4") +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_color_manual(values = c("W1 (8d fixed)" = "#2166AC",
                                          "W2 (adaptive)" = "#D6604D")) +
  ggplot2::scale_x_continuous(breaks = 1:12,
                               labels = substr(month.abb, 1, 3)) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  ggplot2::labs(
    title = "B: NA rain_frac rate by month of peak",
    x     = "Month",
    y     = "Fraction of peaks with NA rain_frac",
    color = "Window"
  ) +
  ggplot2::theme_bw(base_size = 11)

fig1_combined <- p1a / p1b +
  patchwork::plot_annotation(
    title   = "NA rain_frac rate: W1 (8d fixed) vs. W2 (adaptive)",
    caption = sprintf("NOISE_FLOOR_MM = %.1f mm/day  |  Shaded band: approx. snowmelt season (Mar–Jun)",
                      NOISE_FLOOR_MM)
  )

fig1_path <- file.path(out_dir, "fig_na_rainfrac_rate.png")
ggplot2::ggsave(fig1_path, fig1_combined, width = 10, height = 6, dpi = 300)
message(sprintf("Saved: %s", fig1_path))

# ============================================================================
# FIGURE 2 — Reclassification to "unclassified" by season
# ============================================================================
message("\nFigure 2: Reclassification to unclassified by season...")

# Peaks that are W2-classified (final flood_type != "unclassified") but
# would become "unclassified" under W1. Reasons are mutually exclusive:
#   data_insufficient: !suff_W1  (data insufficiency, checked first)
#   na_rain_frac:       suff_W1 & na_rf_W1  (sufficient data but zero active flux)
reclass_df <- peaks_all %>%
  dplyr::filter(!is.na(flood_type),
                as.character(flood_type) != "unclassified") %>%
  dplyr::mutate(
    becomes_unclass_W1 = (!suff_W1) | na_rf_W1,
    reason = dplyr::case_when(
      !suff_W1          ~ "data_insufficient",
      na_rf_W1          ~ "na_rain_frac",
      TRUE              ~ "stays_classified"
    )
  )

# Summary counts for appending to summary.txt later
n_w2_classified  <- nrow(reclass_df)
n_reclassified   <- sum(reclass_df$becomes_unclass_W1, na.rm = TRUE)
n_reclass_insuff <- sum(reclass_df$reason == "data_insufficient", na.rm = TRUE)
n_reclass_narf   <- sum(reclass_df$reason == "na_rain_frac",       na.rm = TRUE)
frac_reclassified <- n_reclassified / n_w2_classified

cat(sprintf(
  "\nW2-classified peaks: %d\n  → become unclassified under W1: %d (%.2f%%)\n    - data_insufficient: %d\n    - na_rain_frac: %d\n",
  n_w2_classified, n_reclassified, 100 * frac_reclassified,
  n_reclass_insuff, n_reclass_narf
))

# Monthly counts
monthly_reclass <- reclass_df %>%
  dplyr::filter(becomes_unclass_W1) %>%
  dplyr::count(peak_month, reason) %>%
  tidyr::complete(
    peak_month = 1:12,
    reason     = c("data_insufficient", "na_rain_frac"),
    fill       = list(n = 0)
  ) %>%
  dplyr::mutate(
    reason = factor(reason,
                    levels = c("data_insufficient", "na_rain_frac"),
                    labels = c("Data insufficient", "NA rain_frac (zero active flux)"))
  )

monthly_total_w2 <- reclass_df %>%
  dplyr::count(peak_month, name = "n_w2_month")

monthly_reclass_frac <- monthly_reclass %>%
  dplyr::group_by(peak_month) %>%
  dplyr::summarise(n_reclass_total = sum(n), .groups = "drop") %>%
  dplyr::left_join(monthly_total_w2, by = "peak_month") %>%
  dplyr::mutate(frac_reclass = n_reclass_total / n_w2_month)

# Panel A: stacked bar — absolute counts
p2a <- ggplot2::ggplot(
    monthly_reclass,
    ggplot2::aes(x = peak_month, y = n, fill = reason)
  ) +
  ggplot2::geom_col(width = 0.8) +
  ggplot2::scale_fill_manual(values = c("Data insufficient"             = "#7570B3",
                                         "NA rain_frac (zero active flux)" = "#D95F02")) +
  ggplot2::scale_x_continuous(breaks = 1:12,
                               labels = substr(month.abb, 1, 3)) +
  ggplot2::labs(
    title = "A: Peaks reclassified to 'unclassified' under W1 — count",
    x     = "Month of peak",
    y     = "Count of reclassified peaks",
    fill  = "Reason"
  ) +
  ggplot2::theme_bw(base_size = 11)

# Panel B: fraction of W2-classified peaks in that month
p2b <- ggplot2::ggplot(
    monthly_reclass_frac,
    ggplot2::aes(x = peak_month, y = frac_reclass)
  ) +
  ggplot2::geom_col(fill = "#636363", width = 0.8) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.1f%%", frac_reclass * 100),
                 y     = frac_reclass + 0.001),
    vjust = 0, size = 2.5
  ) +
  ggplot2::scale_x_continuous(breaks = 1:12,
                               labels = substr(month.abb, 1, 3)) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1),
    expand = ggplot2::expansion(mult = c(0, 0.12))
  ) +
  ggplot2::labs(
    title = "B: Peaks reclassified to 'unclassified' under W1 — fraction of W2-classified",
    x     = "Month of peak",
    y     = "Fraction of W2-classified peaks reclassified"
  ) +
  ggplot2::theme_bw(base_size = 11)

fig2_combined <- p2a / p2b +
  patchwork::plot_annotation(
    title   = "Reclassification to 'unclassified' under W1 — seasonal pattern",
    caption = sprintf(
      "Base: %d peaks classified under W2  |  %d (%.1f%%) become unclassified under W1",
      n_w2_classified, n_reclassified, 100 * frac_reclassified
    )
  )

fig2_path <- file.path(out_dir, "fig_reclassified_unclassified_by_season.png")
ggplot2::ggsave(fig2_path, fig2_combined, width = 10, height = 6, dpi = 300)
message(sprintf("Saved: %s", fig2_path))

# ============================================================================
# FIGURE 3 — Per-site type sample viability under W1 vs. W2 (REVISED)
# ============================================================================
message("\nFigure 3: Per-site type sample viability (revised)...")

# --- same data prep as before: viab_W2, viab_W1, viab_combined ---

# Pre-compute per-type summaries for annotations
fig3a_df <- viab_combined %>%
  dplyr::filter(viable_W2) %>%
  dplyr::group_by(flood_type) %>%
  dplyr::summarise(
    n_total  = dplyr::n(),
    n_lost   = sum(!viable_W1),
    n_kept   = sum( viable_W1),
    pct_lost = mean(!viable_W1) * 100,
    .groups  = "drop"
  ) %>%
  dplyr::mutate(flood_type = factor(flood_type, levels = TYPE_LEVELS))

# Reshape to long for stacked bar
fig3a_long <- fig3a_df %>%
  tidyr::pivot_longer(
    cols      = c(n_kept, n_lost),
    names_to  = "status",
    values_to = "n"
  ) %>%
  dplyr::mutate(
    status = dplyr::case_when(
      status == "n_lost" ~ "Loses viability under W1",
      TRUE               ~ "Retains viability under W1"
    ),
    status = factor(status, levels = c(
      "Retains viability under W1",
      "Loses viability under W1"
    ))
  )

p3a <- ggplot2::ggplot(
    fig3a_long,
    ggplot2::aes(x = flood_type, y = n, fill = status)
  ) +
  ggplot2::geom_col(width = 0.55, color = "white", linewidth = 0.3) +
  # Total label at top of each bar
  ggplot2::geom_text(
    data = fig3a_df,
    ggplot2::aes(x = flood_type, y = n_total,
                 label = sprintf("n = %d", n_total), fill = NULL),
    vjust = -0.4, size = 3.4, color = "grey30", inherit.aes = FALSE
  ) +
  # "Lost" annotation in red, inside or just above the lost segment
  ggplot2::geom_text(
    data = dplyr::filter(fig3a_df, n_lost > 0),
    ggplot2::aes(
      x     = flood_type,
      y     = n_kept + n_lost / 2,
      label = sprintf("%d lost\n(%.0f%%)", n_lost, pct_lost),
      fill  = NULL
    ),
    size = 3.5, color = "white", fontface = "bold", inherit.aes = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Retains viability under W1" = "#2C7BB6",
      "Loses viability under W1"   = "#D7191C"
    ),
    name = NULL
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.15))
  ) +
  ggplot2::labs(
    title    = "A: Sites losing mixture-FFA viability if antecedent window fixed at 8 days",
    subtitle = sprintf(
      "Each bar = currently viable sites under W2 (adaptive) | MIN_TYPE_SAMPLE = %d",
      MIN_TYPE_SAMPLE
    ),
    x = "Flood type",
    y = "Number of sites"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    legend.position  = "top",
    legend.key.size  = ggplot2::unit(0.9, "lines"),
    panel.grid.major.x = ggplot2::element_blank()
  )

# --- Panel B: faceted scatter, tight axis limits, red danger zone ---
fig3b_df <- viab_combined %>%
  dplyr::filter(viable_W2) %>%
  dplyr::mutate(
    loses      = viable_W2 & !viable_W1,
    flood_type = factor(flood_type, levels = TYPE_LEVELS)
  )

ax_min <- MIN_TYPE_SAMPLE - 4
ax_max <- max(c(fig3b_df$n_type_W1, fig3b_df$n_type_W2), na.rm = TRUE) + 3

p3b <- ggplot2::ggplot(
    fig3b_df,
    ggplot2::aes(x = n_type_W2, y = n_type_W1)
  ) +
  # Red shaded danger zone below threshold
  ggplot2::annotate("rect",
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = MIN_TYPE_SAMPLE,
    fill = "#D7191C", alpha = 0.07
  ) +
  ggplot2::geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed", color = "grey50", linewidth = 0.5
  ) +
  ggplot2::geom_hline(
    yintercept = MIN_TYPE_SAMPLE,
    linetype   = "solid", color = "#D7191C", linewidth = 0.75
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      color = flood_type,
      shape = loses,
      size  = loses
    ),
    alpha = 0.75
  ) +
  ggplot2::scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 4),
    labels = c("FALSE" = "Retains viability", "TRUE" = "Loses viability"),
    name   = NULL
  ) +
  ggplot2::scale_size_manual(
    values = c("FALSE" = 2.0, "TRUE" = 3.5),
    guide  = "none"
  ) +
  ggplot2::scale_color_manual(values = TYPE_COLORS, name = "Flood type") +
  ggplot2::annotate(
    "text",
    x     = ax_max,
    y     = MIN_TYPE_SAMPLE + 0.6,
    label = sprintf("MIN_TYPE_SAMPLE = %d", MIN_TYPE_SAMPLE),
    hjust = 1, vjust = 0, size = 2.9, color = "#D7191C"
  ) +
  ggplot2::facet_wrap(
    ~ flood_type, nrow = 1,
    labeller = ggplot2::label_value
  ) +
  ggplot2::coord_cartesian(
    xlim = c(ax_min, ax_max),
    ylim = c(ax_min - 1, ax_max)
  ) +
  ggplot2::labs(
    title    = "B: Peak count per site × flood type: W1 vs. W2",
    subtitle = "W2-viable site × flood type pairs only  |  × marks sites falling below viability threshold under W1",
    x        = "n peaks under W2 (adaptive window)",
    y        = "n peaks under W1 (fixed 8-day window)"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    legend.position    = "bottom",
    strip.text         = ggplot2::element_text(face = "bold"),
    panel.grid.minor   = ggplot2::element_blank()
  )

# Combine: give Panel A roughly 40% of height, Panel B 60%
fig3_combined <- p3a / p3b +
  patchwork::plot_layout(heights = c(2, 3)) +
  patchwork::plot_annotation(
    title   = "Impact of fixing antecedent window to 8 days on mixture-FFA site viability",
    caption = sprintf(
      "One-peak-per-water-year selection applied  |  MIN_TYPE_SAMPLE = %d  |  n_type counted without sample-size demotion",
      MIN_TYPE_SAMPLE
    ),
    theme = ggplot2::theme(
      plot.title   = ggplot2::element_text(face = "bold", size = 14),
      plot.caption = ggplot2::element_text(color = "grey40", size = 9)
    )
  )

fig3_path <- file.path(out_dir, "fig_site_type_viability.png")
ggplot2::ggsave(fig3_path, fig3_combined, width = 12, height = 10, dpi = 300)
message(sprintf("Saved: %s", fig3_path))

# ============================================================================
# FIGURE 4 — Cumulative flux profiles by flood type (lookback_day 1–20)
# ============================================================================
message("\nFigure 4: Cumulative flux profiles...")

# Join W2 flood_type (final) onto accum, filter to data-sufficient peaks
peaks_classified <- peaks_all %>%
  dplyr::filter(!is.na(flood_type),
                as.character(flood_type) != "unclassified",
                suff_W2 == TRUE,
                na_rf_W2 == FALSE) %>%
  dplyr::select(site_no, peak_id, flood_type)

accum_classified <- accum %>%
  dplyr::inner_join(peaks_classified, by = c("site_no", "peak_id")) %>%
  dplyr::filter(lookback_day >= 1L, lookback_day <= 20L)

# Median and 25th–75th percentile of cumulative fluxes at each lookback_day
cumul_summary <- accum_classified %>%
  dplyr::group_by(flood_type, lookback_day) %>%
  dplyr::summarise(
    # Rain
    rain_med = median(cumul_rain_mm, na.rm = TRUE),
    rain_p25 = quantile(cumul_rain_mm, 0.25, na.rm = TRUE),
    rain_p75 = quantile(cumul_rain_mm, 0.75, na.rm = TRUE),
    # Melt
    melt_med = median(cumul_melt_mm, na.rm = TRUE),
    melt_p25 = quantile(cumul_melt_mm, 0.25, na.rm = TRUE),
    melt_p75 = quantile(cumul_melt_mm, 0.75, na.rm = TRUE),
    n_peaks  = dplyr::n_distinct(peak_id),
    .groups  = "drop"
  )

# Reshape for plotting: separate rows for rain vs melt
cumul_long <- dplyr::bind_rows(
  cumul_summary %>%
    dplyr::select(flood_type, lookback_day, n_peaks,
                  med = rain_med, p25 = rain_p25, p75 = rain_p75) %>%
    dplyr::mutate(flux = "Rainfall"),
  cumul_summary %>%
    dplyr::select(flood_type, lookback_day, n_peaks,
                  med = melt_med, p25 = melt_p25, p75 = melt_p75) %>%
    dplyr::mutate(flux = "Snowmelt")
) %>%
  dplyr::mutate(
    flux       = factor(flux, levels = c("Rainfall", "Snowmelt")),
    flood_type = factor(flood_type, levels = setdiff(TYPE_LEVELS, "unclassified"))
  )

# Peak counts per flood type for subtitle
n_per_type <- peaks_classified %>%
  dplyr::count(flood_type, name = "n") %>%
  dplyr::mutate(label = sprintf("%s (n=%d)", flood_type, n))
n_label <- paste(n_per_type$label, collapse = "  |  ")

p4 <- ggplot2::ggplot(
    cumul_long,
    ggplot2::aes(x = lookback_day, color = flux, fill = flux)
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = p25, ymax = p75),
    alpha = 0.20, color = NA
  ) +
  ggplot2::geom_line(ggplot2::aes(y = med), linewidth = 0.9) +
  ggplot2::geom_vline(xintercept = DEPTH_WIN,
                      linetype = "dashed", color = "grey30", linewidth = 0.7) +
  ggplot2::geom_vline(xintercept = W3_window,
                      linetype = "dotted", color = "grey50", linewidth = 0.7) +
  ggplot2::annotate("text", x = DEPTH_WIN + 0.3, y = Inf,
                    label = sprintf("W1\n(%dd)", DEPTH_WIN),
                    hjust = 0, vjust = 1.2, size = 2.8, color = "grey30") +
  ggplot2::annotate("text", x = W3_window + 0.3, y = Inf,
                    label = sprintf("W3\n(%dd)", W3_window),
                    hjust = 0, vjust = 1.2, size = 2.8, color = "grey50") +
  ggplot2::scale_color_manual(values = c(Rainfall = "#2166AC", Snowmelt = "#D6604D")) +
  ggplot2::scale_fill_manual(values  = c(Rainfall = "#2166AC", Snowmelt = "#D6604D")) +
  ggplot2::facet_wrap(~ flood_type, scales = "free_y", ncol = 3) +
  ggplot2::labs(
    title    = "Median cumulative flux profiles by W2 flood type (lookback days 1–20)",
    subtitle = sprintf("Ribbon = 25th–75th pctile  |  %s", n_label),
    x        = "Lookback day (1 = peak date)",
    y        = "Cumulative flux (mm)",
    color    = "Flux type",
    fill     = "Flux type"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(legend.position = "top")

fig4_path <- file.path(out_dir, "fig_cumulative_flux_profiles.png")
ggplot2::ggsave(fig4_path, p4, width = 12, height = 8, dpi = 300)
message(sprintf("Saved: %s", fig4_path))

# ============================================================================
# APPEND TO SUMMARY.TXT
# ============================================================================
message("\nAppending to summary.txt...")

# Per-type viability loss counts
viab_by_type <- lost_viability %>%
  dplyr::count(flood_type, name = "n") %>%
  dplyr::arrange(flood_type)

get_viab_count <- function(type_str) {
  val <- viab_by_type %>%
    dplyr::filter(as.character(flood_type) == type_str) %>%
    dplyr::pull(n)
  if (length(val) == 0) 0L else val
}

# NA rate totals
n_total_peaks <- nrow(peaks_all)
n_na_W1_total <- sum(peaks_all$na_rf_W1, na.rm = TRUE)
n_na_W2_total <- sum(peaks_all$na_rf_W2, na.rm = TRUE)

summary_append <- c(
  "",
  strrep("-", 60),
  "EXTENDED DIAGNOSTIC — OPERATIONAL IMPACT OF W1 (fixed 8d)",
  strrep("-", 60),
  sprintf("  NA rain_frac rate under W1          : %.4f  (%d / %d)",
          n_na_W1_total / n_total_peaks, n_na_W1_total, n_total_peaks),
  sprintf("  NA rain_frac rate under W2          : %.4f  (%d / %d)",
          n_na_W2_total / n_total_peaks, n_na_W2_total, n_total_peaks),
  sprintf("  Peaks reclassified to unclassified  : %d / %d (%.1f%%)",
          n_reclassified, n_w2_classified, 100 * frac_reclassified),
  sprintf("    - Due to data insufficiency       : %d", n_reclass_insuff),
  sprintf("    - Due to NA rain_frac             : %d", n_reclass_narf),
  sprintf("  Site-type combinations viable W2    : %d", total_viable_W2),
  sprintf("  Site-type combinations losing viability under W1 : %d",
          total_loses_viab),
  sprintf("    - rainfall    : %d", get_viab_count("rainfall")),
  sprintf("    - snowmelt    : %d", get_viab_count("snowmelt")),
  sprintf("    - rain_on_snow: %d", get_viab_count("rain_on_snow"))
)

cat(summary_append, file = summary_file, sep = "\n", append = TRUE)
message(sprintf("Appended to: %s", summary_file))

# ============================================================================
# COMPLETION
# ============================================================================
written_files <- c(fig1_path, fig2_path, fig3_path, fig4_path, summary_file)

cat("\n", strrep("=", 70), "\n", sep = "")
cat("EXTENDED DIAGNOSTIC COMPLETE\n")
cat(strrep("=", 70), "\n\n", sep = "")
cat("Files written:\n")
purrr::walk(written_files, ~ cat(sprintf("  %s\n", .x)))
