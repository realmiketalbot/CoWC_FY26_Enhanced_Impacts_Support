#!/usr/bin/env Rscript

# ============================================================================
# Script 9: Site-Level Non-Parametric Trend Analysis of Flood Type Frequency
# ============================================================================
#
# PURPOSE:
#   Tests for monotonic trends in the annual frequency and fractional
#   composition of classified flood types at all 206 Colorado GAGES-II sites,
#   across three classification methods:
#
#     depth_thresh       : Rainfall / Snowmelt / Rain-on-Snow (Davenport rules)
#     frac_thresh_3class : Rainfall / Snowmelt / Rain-on-Snow (rain_frac split)
#     frac_thresh_2class : Rainfall / Non-Rainfall (derived from 3-class peaks
#                          by collapsing Snowmelt + Rain-on-Snow to Non-Rainfall)
#
# METHODOLOGICAL CHOICES:
#
#   Why flood_type_raw (not flood_type)?
#     flood_type_raw is the decision-tree label before the MIN_TYPE_SAMPLE
#     demotion applied in Scripts 6a/6b. flood_type collapses rare types to
#     "unclassified" at sites where a type has too few peaks for GPD/LP3
#     fitting. For temporal analysis, we want every classified event even at
#     sites where a type is too rare for FFA, so the pre-demotion label is
#     the correct choice.
#
#   Why all classified peaks (not the post-selection top-n)?
#     The top-n sample in Script 7 is tuned for extreme-value fitting. Temporal
#     composition trends need every classified peak across the full record.
#
#   Why both count and proportion?
#     Count trends detect absolute changes in peak occurrence rate.
#     Proportion trends detect composition shifts independent of total peak rate.
#
#   Why eligibility thresholds (>=20 classified years, >=5 type years)?
#     The Mann-Kendall test has very low power for n < 20 (Yue et al. 2002).
#
#   Why a binomial test against a 5% null?
#     Under a true null, exactly 5% of MK tests are significant by chance.
#     Observing significantly more constitutes regional-scale evidence.
#
# REQUIRED PACKAGES:
#   install.packages(c("Kendall", "trend", "patchwork", "maps"))
#
# INPUTS:
#   _data/pot_peaks_classified_depth_thresh.rds   (Script 6a)
#   _data/pot_peaks_classified_frac_thresh_3class.rds (Script 6b)
#   _data/site_metadata.rds                       (Script 0)
#   _data/station_spatial.rds                     (Script 1)
#   functions.R
#
# OUTPUTS (per method):
#   _data/temporal_trends_{method}.rds            list(site_level, aggregate)
#   _data/annual_type_counts_{method}.rds         one row per site x WY x type
#   _diagnostics/temporal_trends/{method}/        diagnostic figures
#
# OUTPUTS (AM analysis, per method):
#   _data/annual_am_type_{method}.rds         one row per site x WY (AM peak type)
#   _data/am_type_shift_{method}.rds          first- vs. second-half AM type frequency
#   _data/temporal_trends_am_{method}.rds     list(am_type_shift, am_type_indicator_trends)
#                                             [frac_thresh_3class also includes am_rain_frac_trends]
#
# Author : Mike Talbot, PhD Student, Colorado State University
# Date   : 2026-05-21
# ============================================================================

## ── Required-package check ───────────────────────────────────────────────────
required_pkgs <- c("Kendall", "trend", "patchwork", "maps")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace,
                                       quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop(
    "\nThe following packages are required but not installed:\n  ",
    paste(missing_pkgs, collapse = ", "),
    "\n\nInstall with:\n  install.packages(c(",
    paste0('"', missing_pkgs, '"', collapse = ", "), "))\n"
  )
}

## ── Libraries ────────────────────────────────────────────────────────────────
library(tidyverse)
library(patchwork)

source("functions.R")

## ── Shared configuration ─────────────────────────────────────────────────────
MIN_CLASSIFIED_YEARS <- 20L
MIN_TYPE_YEARS       <- 5L
MK_ALPHA             <- 0.05

# Per-method configuration. frac_thresh_2class derives from the 3-class peaks
# by collapsing snowmelt + rain_on_snow to non_rainfall.
METHOD_CONFIGS <- list(

  depth_thresh = list(
    peaks_file   = "_data/pot_peaks_classified_depth_thresh.rds",
    remap_2class = FALSE,
    flood_types  = c("rainfall", "snowmelt", "rain_on_snow"),
    type_colors  = c(rainfall = "#b2182b", snowmelt = "#2166ac",
                     rain_on_snow = "#9970ab"),
    type_labels  = c(rainfall = "Rainfall", snowmelt = "Snowmelt",
                     rain_on_snow = "Rain-on-Snow"),
    method_label = "Depth (3-class)"
  ),

  frac_thresh_3class = list(
    peaks_file   = "_data/pot_peaks_classified_frac_thresh_3class.rds",
    remap_2class = FALSE,
    flood_types  = c("rainfall", "snowmelt", "rain_on_snow"),
    type_colors  = c(rainfall = "#b2182b", snowmelt = "#2166ac",
                     rain_on_snow = "#9970ab"),
    type_labels  = c(rainfall = "Rainfall", snowmelt = "Snowmelt",
                     rain_on_snow = "Rain-on-Snow"),
    method_label = "Rain fraction (3-class)"
  ),

  frac_thresh_2class = list(
    peaks_file   = "_data/pot_peaks_classified_frac_thresh_3class.rds",
    remap_2class = TRUE,
    flood_types  = c("rainfall", "non_rainfall"),
    type_colors  = c(rainfall = "#b2182b", non_rainfall = "#9970ab"),
    type_labels  = c(rainfall = "Rainfall", non_rainfall = "Non-Rainfall"),
    method_label = "Rain fraction (2-class)"
  )
)

## ── Output directory ─────────────────────────────────────────────────────────
diag_dir <- "_diagnostics/temporal_trends"
if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)

## ── Shared spatial data (loaded once) ────────────────────────────────────────
message("Loading shared spatial / metadata...")

meta    <- readRDS("_data/site_metadata.rds")
spatial <- local({
  obj <- readRDS("_data/station_spatial.rds")
  if (inherits(obj, "sf")) sf::st_drop_geometry(obj) else obj
})
spatial <- as.data.frame(spatial)

# Colorado boundary polygon for maps
co_boundary <- maps::map("state", regions = "colorado", plot = FALSE, fill = TRUE)
co_poly <- data.frame(
  long  = co_boundary$x,
  lat   = co_boundary$y,
  group = cumsum(is.na(co_boundary$x) | is.na(co_boundary$y))
)

spatial_pts <- spatial %>% dplyr::select(site_no, lon, lat)

## ── Regulated flag lookup ─────────────────────────────────────────────────────
if ("regulated" %in% names(meta)) {
  reg_lookup <- dplyr::select(meta, site_no, regulated)
  message("  Regulated flag : from site_metadata.rds")
} else {
  stop("'regulated' column not found in site_metadata.rds.")
}

## ── Shared helper: Mann-Kendall + Theil-Sen trend tests ──────────────────────
run_trend_tests <- function(site, ftype, regulated_flag, annual_ts_data) {

  site_yr    <- dplyr::filter(annual_ts_data, site_no == site) %>%
    dplyr::arrange(water_year)
  count_col  <- paste0("n_",    ftype)
  prop_col   <- paste0("prop_", ftype)
  count_series <- site_yr[[count_col]]
  prop_df      <- dplyr::filter(site_yr, classified_total > 0L)
  prop_series  <- prop_df[[prop_col]]

  run_one <- function(series, metric_label) {
    series_clean <- series[!is.na(series)]
    mk_res <- tryCatch(Kendall::MannKendall(series_clean), error = function(e) NULL)
    ss_res <- tryCatch(trend::sens.slope(series_clean),    error = function(e) NULL)
    tibble::tibble(
      site_no    = site,
      flood_type = ftype,
      metric     = metric_label,
      regulated  = regulated_flag,
      n_years    = length(series_clean),
      n_nonzero  = sum(series_clean > 0, na.rm = TRUE),
      mean_val   = mean(series_clean, na.rm = TRUE),
      sd_val     = sd(series_clean,   na.rm = TRUE),
      tau        = if (!is.null(mk_res)) unname(mk_res$tau[1])    else NA_real_,
      mk_pvalue  = if (!is.null(mk_res)) unname(mk_res$sl[1])     else NA_real_,
      sen_slope  = if (!is.null(ss_res)) unname(ss_res$estimates) else NA_real_
    )
  }

  dplyr::bind_rows(run_one(count_series, "count"),
                   run_one(prop_series,  "proportion"))
}

## ── Aggregate helper (shared) ─────────────────────────────────────────────────
run_aggregate <- function(df, pop_label) {
  df %>%
    dplyr::group_by(flood_type, metric) %>%
    dplyr::summarise(
      population      = pop_label,
      n_eligible      = dplyr::n(),
      n_sig           = sum(significant,  na.rm = TRUE),
      n_sig_pos       = sum(sig_positive, na.rm = TRUE),
      n_sig_neg       = sum(sig_negative, na.rm = TRUE),
      binom_p_any_sig = tryCatch(
        binom.test(sum(significant,  na.rm = TRUE), dplyr::n(), p = MK_ALPHA)$p.value,
        error = function(e) NA_real_
      ),
      binom_p_pos_sig = tryCatch(
        binom.test(sum(sig_positive, na.rm = TRUE), dplyr::n(), p = MK_ALPHA)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      pct_sig     = round(100 * n_sig     / n_eligible, 1),
      pct_sig_pos = round(100 * n_sig_pos / n_eligible, 1),
      pct_sig_neg = round(100 * n_sig_neg / n_eligible, 1)
    )
}

## ═══════════════════════════════════════════════════════════════════════════════
## METHOD LOOP
## ═══════════════════════════════════════════════════════════════════════════════

for (method_key in names(METHOD_CONFIGS)) {

  cfg <- METHOD_CONFIGS[[method_key]]

  message(sprintf("\n%s\n=== METHOD: %s (%s) ===\n%s",
                  strrep("=", 60), method_key, cfg$method_label, strrep("=", 60)))

  method_diag_dir <- file.path(diag_dir, method_key)
  if (!dir.exists(method_diag_dir)) dir.create(method_diag_dir, recursive = TRUE)

  ## ── Load and prepare peaks ─────────────────────────────────────────────────
  message("Loading data...")

  peaks_raw <- readRDS(cfg$peaks_file)

  # For frac_thresh_2class, remap snowmelt + rain_on_snow → non_rainfall
  if (isTRUE(cfg$remap_2class)) {
    peaks_raw <- peaks_raw %>%
      dplyr::mutate(
        flood_type_raw = dplyr::case_when(
          as.character(flood_type_raw) == "rainfall"     ~ "rainfall",
          as.character(flood_type_raw) == "unclassified" ~ "unclassified",
          TRUE                                           ~ "non_rainfall"
        )
      )
    message(sprintf(
      "  Remapped: snowmelt + rain_on_snow -> non_rainfall (%s peaks)",
      format(sum(as.character(peaks_raw$flood_type_raw) == "non_rainfall"),
             big.mark = ",")
    ))
  }

  message(sprintf("  Peaks raw     : %d rows | %d sites | WY %d–%d",
                  nrow(peaks_raw),
                  dplyr::n_distinct(peaks_raw$site_no),
                  min(peaks_raw$water_year, na.rm = TRUE),
                  max(peaks_raw$water_year, na.rm = TRUE)))

  # Deduplicate (one row per peak_id)
  peaks <- dplyr::distinct(peaks_raw, site_no, peak_id, .keep_all = TRUE)
  n_dropped <- nrow(peaks_raw) - nrow(peaks)
  if (n_dropped > 0L)
    message(sprintf("  Deduplication : dropped %d duplicate rows", n_dropped))

  peaks <- peaks %>%
    dplyr::select(-dplyr::any_of("regulated")) %>%
    dplyr::left_join(reg_lookup, by = "site_no")

  message(sprintf("  Peaks (deduped): %d rows | %d sites",
                  nrow(peaks), dplyr::n_distinct(peaks$site_no)))

  ## ── Build annual time series (dynamic flood types) ─────────────────────────
  message("\nBuilding annual type counts per site...")

  # Base counts per (site, water_year)
  base_annual <- peaks %>%
    dplyr::group_by(site_no, water_year) %>%
    dplyr::summarise(
      n_total_peaks    = dplyr::n(),
      classified_total = sum(as.character(flood_type_raw) %in% cfg$flood_types,
                             na.rm = TRUE),
      regulated        = dplyr::first(regulated),
      .groups = "drop"
    )

  # Per-type counts via pivot_wider
  type_wide <- peaks %>%
    dplyr::mutate(ftr = as.character(flood_type_raw)) %>%
    dplyr::filter(ftr %in% cfg$flood_types) %>%
    dplyr::count(site_no, water_year, ftr, name = "count") %>%
    tidyr::pivot_wider(names_from = ftr, values_from = count,
                       values_fill = 0L, names_prefix = "n_")

  # Ensure every expected column is present (e.g. a type with zero peaks globally)
  for (ft in cfg$flood_types) {
    col <- paste0("n_", ft)
    if (!col %in% names(type_wide)) type_wide[[col]] <- 0L
  }

  annual_ts <- dplyr::left_join(base_annual, type_wide,
                                by = c("site_no", "water_year")) %>%
    dplyr::mutate(dplyr::across(
      dplyr::all_of(paste0("n_", cfg$flood_types)),
      ~ as.integer(tidyr::replace_na(., 0L))
    ))

  # Add proportion columns
  for (ft in cfg$flood_types) {
    n_col    <- paste0("n_", ft)
    prop_col <- paste0("prop_", ft)
    annual_ts[[prop_col]] <- ifelse(
      annual_ts$classified_total > 0L,
      annual_ts[[n_col]] / annual_ts$classified_total,
      NA_real_
    )
  }

  message(sprintf("  Annual series : %d site-year rows | WY %d–%d",
                  nrow(annual_ts),
                  min(annual_ts$water_year),
                  max(annual_ts$water_year)))

  ## Build Shiny-facing long-form annual counts tibble
  annual_type_counts <- purrr::map_dfr(cfg$flood_types, function(ft) {
    annual_ts %>%
      dplyr::transmute(
        site_no,
        water_year,
        flood_type       = ft,
        count            = .data[[paste0("n_", ft)]],
        classified_total = classified_total,
        proportion       = .data[[paste0("prop_", ft)]],
        regulated        = regulated
      )
  }) %>%
    dplyr::arrange(site_no, water_year, flood_type)

  ## ── Site record stats ──────────────────────────────────────────────────────
  site_record_stats <- annual_ts %>%
    dplyr::group_by(site_no) %>%
    dplyr::summarise(
      n_wy_total       = dplyr::n(),
      n_classified_yrs = sum(classified_total >= 1L),
      regulated        = dplyr::first(regulated),
      dplyr::across(
        dplyr::all_of(paste0("n_", cfg$flood_types)),
        ~ sum(. >= 1L, na.rm = TRUE),
        .names = "{.col}_yrs"
      ),
      .groups = "drop"
    )

  ## ── DIAGNOSTIC FIGURES ────────────────────────────────────────────────────
  message("\nGenerating diagnostic figures...")

  # Colour vector keyed by display label (for ggplot fill scales)
  type_color_vec <- setNames(cfg$type_colors, cfg$type_labels[cfg$flood_types])

  lat_order <- spatial %>%
    dplyr::select(site_no, lat) %>%
    dplyr::arrange(lat) %>%
    dplyr::pull(site_no)

  # Diag 1: Site × water-year heatmap
  p_heatmap <- ggplot2::ggplot(
    dplyr::mutate(annual_ts, site_f = factor(site_no, levels = lat_order)),
    ggplot2::aes(x = water_year, y = site_f, fill = classified_total)
  ) +
    ggplot2::geom_tile(na.rm = TRUE) +
    ggplot2::scale_fill_viridis_c(option = "plasma", name = "Classified\npeaks/yr",
                                  na.value = "white") +
    ggplot2::scale_x_continuous(breaks = seq(1990, 2025, 5), expand = c(0, 0)) +
    ggplot2::labs(
      title    = sprintf("Annual classified peaks per site — %s", cfg$method_label),
      subtitle = "Sites ordered south-to-north; white gaps = no PoT peaks that year",
      x = "Water year", y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())

  ggplot2::ggsave(
    file.path(method_diag_dir, "heatmap_site_year_classified_peaks.png"),
    p_heatmap, width = 14, height = 10, dpi = 150
  )
  message("  Saved: heatmap_site_year_classified_peaks.png")

  # Diag 2: Regional stacked bar
  regional_bar_df <- annual_type_counts %>%
    dplyr::group_by(water_year, flood_type) %>%
    dplyr::summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      flood_type = factor(flood_type, levels = cfg$flood_types),
      type_label = cfg$type_labels[as.character(flood_type)]
    )

  p_regional_bar <- ggplot2::ggplot(
    regional_bar_df,
    ggplot2::aes(x = water_year, y = total_count, fill = type_label)
  ) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::scale_fill_manual(values = type_color_vec, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(1990, 2025, 5)) +
    ggplot2::labs(
      title    = sprintf("Regional peak count by flood type — %s (all sites combined)",
                         cfg$method_label),
      subtitle = "Inspect for regional-scale composition shifts before formal testing",
      x = "Water year", y = "Total classified peaks"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")

  ggplot2::ggsave(
    file.path(method_diag_dir, "regional_stacked_bar_flood_type.png"),
    p_regional_bar, width = 12, height = 6, dpi = 150
  )
  message("  Saved: regional_stacked_bar_flood_type.png")

  # Diag 3: Record-length histograms
  hist_long <- site_record_stats %>%
    tidyr::pivot_longer(
      cols          = dplyr::all_of(paste0("n_", cfg$flood_types, "_yrs")),
      names_to      = "flood_type",
      names_pattern = "n_(.+)_yrs",
      values_to     = "n_type_yrs"
    ) %>%
    dplyr::mutate(
      type_label = factor(cfg$type_labels[flood_type],
                          levels = unname(cfg$type_labels[cfg$flood_types]))
    )

  p_type_hist <- ggplot2::ggplot(
    hist_long, ggplot2::aes(x = n_type_yrs, fill = type_label)
  ) +
    ggplot2::geom_histogram(binwidth = 2, colour = "white", linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = MIN_TYPE_YEARS,
                        linetype = "dashed", colour = "grey20", linewidth = 0.8) +
    ggplot2::scale_fill_manual(values = type_color_vec, name = NULL) +
    ggplot2::facet_wrap(~type_label, ncol = length(cfg$flood_types)) +
    ggplot2::labs(
      title    = sprintf("Water years with >=1 peak of each type — %s (n = %d sites)",
                         cfg$method_label, nrow(site_record_stats)),
      subtitle = sprintf("Dashed line: per-type eligibility minimum (%d yr)", MIN_TYPE_YEARS),
      x = "Water years with >=1 peak of type", y = "Number of sites"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "none")

  p_record_hist <- ggplot2::ggplot(
    site_record_stats, ggplot2::aes(x = n_classified_yrs)
  ) +
    ggplot2::geom_histogram(fill = "#4dac26", colour = "white", binwidth = 2,
                            linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = MIN_CLASSIFIED_YEARS,
                        linetype = "dashed", colour = "grey20", linewidth = 0.8) +
    ggplot2::labs(
      title    = "Water years with >=1 classified peak per site",
      subtitle = sprintf("Dashed line: record eligibility minimum (%d yr)",
                         MIN_CLASSIFIED_YEARS),
      x = "Years with >=1 classified peak", y = "Number of sites"
    ) +
    ggplot2::theme_minimal(base_size = 11)

  p_hist_combined <- (p_record_hist / p_type_hist) +
    patchwork::plot_annotation(
      title = sprintf("Site record-length diagnostics — %s (all sites, pre-eligibility)",
                      cfg$method_label)
    )

  ggplot2::ggsave(
    file.path(method_diag_dir, "histogram_record_length_by_type.png"),
    p_hist_combined, width = 14, height = 9, dpi = 150
  )
  message("  Saved: histogram_record_length_by_type.png")

  ## ── AM analysis: annual maximum classification ────────────────────────────
  message("\nDeriving annual maximum (AM) classification series...")

  ## 1. AM peak per site-year (all methods) ────────────────────────────────────
  has_rain_frac <- "rain_frac" %in% names(peaks)

  annual_am_type <- peaks %>%
    dplyr::filter(as.character(flood_type_raw) %in% cfg$flood_types) %>%
    dplyr::group_by(site_no, water_year) %>%
    dplyr::slice_max(peak_flow, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      site_no,
      water_year,
      peak_flow,
      am_flood_type = as.character(flood_type_raw),
      am_rain_frac  = if (has_rain_frac) rain_frac else NA_real_,
      regulated
    ) %>%
    dplyr::arrange(site_no, water_year)

  message(sprintf("  AM series: %d site-year rows | %d sites",
                  nrow(annual_am_type),
                  dplyr::n_distinct(annual_am_type$site_no)))

  saveRDS(annual_am_type, sprintf("_data/annual_am_type_%s.rds", method_key))
  message(sprintf("  Saved: _data/annual_am_type_%s.rds", method_key))

  ## Eligible sites for AM trend tests (defined here for all methods) ──────────
  am_eligible_sites <- annual_am_type %>%
    dplyr::group_by(site_no) %>%
    dplyr::summarise(
      n_am_yrs  = dplyr::n(),
      regulated = dplyr::first(regulated),
      .groups   = "drop"
    ) %>%
    dplyr::filter(n_am_yrs >= MIN_CLASSIFIED_YEARS)

  message(sprintf("  Eligible sites (>=%d classified AM years): %d",
                  MIN_CLASSIFIED_YEARS, nrow(am_eligible_sites)))

  ## 2. Mann-Kendall trend test on AM rain fraction (frac_thresh_3class only) ──
  if (method_key == "frac_thresh_3class") {

    message("\nRunning Mann-Kendall trend tests on AM rain fraction per site...")

    run_mk_am <- function(site, regulated_flag) {
      site_ts  <- dplyr::filter(annual_am_type, site_no == site) %>%
        dplyr::arrange(water_year)
      vals     <- site_ts$am_rain_frac
      vals_ok  <- vals[!is.na(vals)]
      n_series <- length(vals_ok)
      if (n_series < MIN_CLASSIFIED_YEARS) return(NULL)
      mk_res <- tryCatch(Kendall::MannKendall(vals_ok), error = function(e) NULL)
      ss_res <- tryCatch(trend::sens.slope(vals_ok),    error = function(e) NULL)
      tibble::tibble(
        site_no   = site,
        n_years   = n_series,
        tau       = if (!is.null(mk_res)) unname(mk_res$tau[1])    else NA_real_,
        p_value   = if (!is.null(mk_res)) unname(mk_res$sl[1])     else NA_real_,
        sen_slope = if (!is.null(ss_res)) unname(ss_res$estimates) else NA_real_,
        regulated = regulated_flag
      )
    }

    am_rain_frac_trends <- purrr::map2_dfr(
      am_eligible_sites$site_no,
      am_eligible_sites$regulated,
      run_mk_am
    ) %>%
      dplyr::mutate(significant = !is.na(p_value) & p_value < MK_ALPHA) %>%
      dplyr::arrange(site_no)

    message(sprintf("  AM MK tests complete: %d sites", nrow(am_rain_frac_trends)))

  }  # end frac_thresh_3class rain fraction section

  ## 3. First-half vs. second-half AM type frequency comparison (all methods) ──
  message("\nBuilding first-half vs. second-half AM type frequency comparison...")

  am_type_shift <- annual_am_type %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(
      first_wy = min(water_year),
      last_wy  = max(water_year),
      split_wy = floor(stats::median(water_year))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(half = ifelse(water_year < split_wy, "first", "second")) %>%
    dplyr::group_by(site_no, am_flood_type) %>%
    dplyr::summarise(
      n_first_half  = sum(half == "first",  na.rm = TRUE),
      n_second_half = sum(half == "second", na.rm = TRUE),
      n_total_yrs   = dplyr::n(),
      first_wy      = dplyr::first(first_wy),
      split_wy      = dplyr::first(split_wy),
      last_wy       = dplyr::first(last_wy),
      regulated     = dplyr::first(regulated),
      .groups       = "drop"
    ) %>%
    dplyr::rename(flood_type = am_flood_type) %>%
    dplyr::select(site_no, flood_type, n_first_half, n_second_half, n_total_yrs,
                  first_wy, split_wy, last_wy, regulated) %>%
    dplyr::arrange(site_no, flood_type)

  message(sprintf("  AM type shift: %d site x type rows", nrow(am_type_shift)))

  saveRDS(am_type_shift, sprintf("_data/am_type_shift_%s.rds", method_key))
  message(sprintf("  Saved: _data/am_type_shift_%s.rds", method_key))

  ## 3b. Binary AM type indicator MK tests (all methods) ──────────────────────
  ## Binary indicator = 1 if that type was the AM type that year, 0 otherwise.
  ## NA am_flood_type (unclassified AM peak) is treated as 0, not a type
  ## occurrence. No per-type minimum count threshold is applied — the series is
  ## well-defined (all zeros) even if the type never generated the AM peak.
  message("\nRunning MK trend tests on binary AM type indicator series...")

  run_mk_indicator <- function(site, ftype, regulated_flag) {
    site_ts   <- dplyr::filter(annual_am_type, site_no == site) %>%
      dplyr::arrange(water_year)
    indicator <- as.integer(!is.na(site_ts$am_flood_type) &
                              site_ts$am_flood_type == ftype)
    n_series  <- length(indicator)
    if (n_series < MIN_CLASSIFIED_YEARS) return(NULL)
    mk_res <- tryCatch(Kendall::MannKendall(indicator), error = function(e) NULL)
    ss_res <- tryCatch(trend::sens.slope(indicator),    error = function(e) NULL)
    tibble::tibble(
      site_no      = site,
      flood_type   = ftype,
      regulated    = regulated_flag,
      n_years      = n_series,
      n_type_years = sum(indicator),
      tau          = if (!is.null(mk_res)) unname(mk_res$tau[1])    else NA_real_,
      p_value      = if (!is.null(mk_res)) unname(mk_res$sl[1])     else NA_real_,
      sen_slope    = if (!is.null(ss_res)) unname(ss_res$estimates) else NA_real_
    )
  }

  am_indicator_types <- cfg$flood_types

  am_type_indicator_trends <- purrr::map_dfr(
    am_indicator_types,
    function(ftype) {
      purrr::map2_dfr(
        am_eligible_sites$site_no,
        am_eligible_sites$regulated,
        ~ run_mk_indicator(.x, ftype, .y)
      )
    }
  ) %>%
    dplyr::mutate(
      significant  = !is.na(p_value) & p_value < MK_ALPHA,
      sig_positive = significant & !is.na(sen_slope) & sen_slope > 0,
      sig_negative = significant & !is.na(sen_slope) & sen_slope < 0
    ) %>%
    dplyr::arrange(flood_type, site_no)

  message(sprintf("  AM indicator MK tests complete: %d site x type rows",
                  nrow(am_type_indicator_trends)))

  # Diagnostic: sites with zero occurrences of each type (degenerate MK result)
  for (ftype in am_indicator_types) {
    sub    <- dplyr::filter(am_type_indicator_trends, flood_type == ftype)
    n_zero <- sum(sub$n_type_years == 0L, na.rm = TRUE)
    message(sprintf(
      "  '%s': %d / %d eligible sites with zero AM years of this type (degenerate MK)",
      ftype, n_zero, nrow(sub)
    ))
  }

  # Console summary per type
  message(sprintf("\n  Binary AM type MK summary (alpha = %.2f):", MK_ALPHA))
  for (ftype in am_indicator_types) {
    sub <- dplyr::filter(am_type_indicator_trends, flood_type == ftype)
    message(sprintf(
      "  %-14s : %d eligible | %d sig | %d sig+ | %d sig-",
      ftype,
      nrow(sub),
      sum(sub$significant,  na.rm = TRUE),
      sum(sub$sig_positive, na.rm = TRUE),
      sum(sub$sig_negative, na.rm = TRUE)
    ))
  }

  ## 4. Save combined AM analysis output (all methods) ─────────────────────────
  am_output <- list(
    am_type_shift            = am_type_shift,
    am_type_indicator_trends = am_type_indicator_trends
  )

  if (method_key == "frac_thresh_3class") {
    am_output$am_rain_frac_trends <- am_rain_frac_trends
  }

  saveRDS(am_output, sprintf("_data/temporal_trends_am_%s.rds", method_key))
  message(sprintf("  Saved: _data/temporal_trends_am_%s.rds", method_key))

  ## ── Site eligibility filter ────────────────────────────────────────────────
  message(sprintf(
    "\nApplying eligibility filter: >=%d classified WY AND >=%d WY with that type...",
    MIN_CLASSIFIED_YEARS, MIN_TYPE_YEARS
  ))

  eligibility <- site_record_stats %>%
    dplyr::mutate(eligible_record = n_classified_yrs >= MIN_CLASSIFIED_YEARS) %>%
    tidyr::pivot_longer(
      cols          = dplyr::all_of(paste0("n_", cfg$flood_types, "_yrs")),
      names_to      = "flood_type",
      names_pattern = "n_(.+)_yrs",
      values_to     = "n_type_yrs"
    ) %>%
    dplyr::mutate(eligible = eligible_record & n_type_yrs >= MIN_TYPE_YEARS)

  n_total_sites <- dplyr::n_distinct(peaks$site_no)

  eligibility_summary <- eligibility %>%
    dplyr::group_by(flood_type) %>%
    dplyr::summarise(
      n_eligible   = sum(eligible),
      n_total      = dplyr::n(),
      pct_eligible = round(100 * n_eligible / n_total, 1),
      .groups = "drop"
    )

  message("\n  Eligibility summary (record >=20 yr AND type >=5 yr):")
  for (i in seq_len(nrow(eligibility_summary))) {
    r   <- eligibility_summary[i, ]
    msg <- sprintf("    %-14s : %d / %d sites (%.0f%%) eligible",
                   r$flood_type, r$n_eligible, r$n_total, r$pct_eligible)
    if (r$pct_eligible < 50) {
      message(sprintf("  *** FLAG: %s eligible for only %.0f%% of sites",
                      r$flood_type, r$pct_eligible))
      message(paste0("  ***       ", msg))
    } else {
      message(msg)
    }
  }

  eligible_combos <- eligibility %>%
    dplyr::filter(eligible) %>%
    dplyr::select(site_no, flood_type, regulated) %>%
    dplyr::mutate(flood_type = as.character(flood_type))

  message(sprintf("\n  Total eligible (site x type) combinations : %d",
                  nrow(eligible_combos)))

  ## ── Per-site trend tests ───────────────────────────────────────────────────
  message(sprintf(
    "\nRunning Mann-Kendall + Theil-Sen tests (%d combinations)...",
    nrow(eligible_combos)
  ))

  site_level_list <- purrr::pmap(
    list(eligible_combos$site_no,
         eligible_combos$flood_type,
         eligible_combos$regulated),
    function(s, f, r) run_trend_tests(s, f, r, annual_ts)
  )

  site_level <- dplyr::bind_rows(site_level_list) %>%
    dplyr::mutate(
      significant  = !is.na(mk_pvalue) & mk_pvalue < MK_ALPHA,
      sig_positive = significant & !is.na(sen_slope) & sen_slope > 0,
      sig_negative = significant & !is.na(sen_slope) & sen_slope < 0
    )

  message(sprintf("  Tests complete : %d site x type x metric rows", nrow(site_level)))

  ## ── Aggregate across sites ─────────────────────────────────────────────────
  message("\nAggregating trend results...")

  aggregate <- dplyr::bind_rows(
    run_aggregate(site_level, "all_eligible"),
    run_aggregate(dplyr::filter(site_level, !is.na(regulated) & !regulated),
                  "unregulated")
  ) %>%
    dplyr::select(
      flood_type, metric, population, n_eligible,
      n_sig, pct_sig, n_sig_pos, pct_sig_pos, n_sig_neg, pct_sig_neg,
      binom_p_any_sig, binom_p_pos_sig
    ) %>%
    dplyr::arrange(flood_type, metric, population)

  ## ── Console summary ────────────────────────────────────────────────────────
  message(sprintf("\n%s\n=== AGGREGATE TREND SUMMARY — %s (alpha = %.2f) ===\n%s",
                  strrep("=", 60), cfg$method_label, MK_ALPHA, strrep("=", 60)))
  message(sprintf("  Binomial null: p = %.2f\n", MK_ALPHA))

  for (pop in unique(aggregate$population)) {
    n_pop <- if (pop == "unregulated")
      sum(!reg_lookup$regulated, na.rm = TRUE) else n_total_sites
    message(sprintf("  Population: %s", pop))
    sub <- dplyr::filter(aggregate, population == pop)
    for (i in seq_len(nrow(sub))) {
      r        <- sub[i, ]
      sig_flag <- if (!is.na(r$binom_p_any_sig) && r$binom_p_any_sig < MK_ALPHA) " *" else ""
      pos_flag <- if (!is.na(r$binom_p_pos_sig) && r$binom_p_pos_sig < MK_ALPHA) " *" else ""
      message(sprintf(
        "    %-14s %-12s | n=%3d | sig: %3d (%.0f%%) [binom p=%.3f%s] | pos: %3d (%.0f%%) [p=%.3f%s] | neg: %3d (%.0f%%)",
        r$flood_type, r$metric, r$n_eligible,
        r$n_sig,     r$pct_sig,     r$binom_p_any_sig, sig_flag,
        r$n_sig_pos, r$pct_sig_pos, r$binom_p_pos_sig, pos_flag,
        r$n_sig_neg, r$pct_sig_neg
      ))
    }
    message("")
  }
  message(strrep("=", 60))
  message("  * = binomial p < 0.05 vs. 5% null\n")

  ## ── Save data outputs ──────────────────────────────────────────────────────
  message("Saving output files...")

  trends_out <- sprintf("_data/temporal_trends_%s.rds",     method_key)
  counts_out <- sprintf("_data/annual_type_counts_%s.rds",  method_key)

  saveRDS(list(site_level = site_level, aggregate = aggregate), trends_out)
  message(sprintf("  Saved: %s", trends_out))

  saveRDS(annual_type_counts, counts_out)
  message(sprintf("  Saved: %s", counts_out))

  ## ── VISUALIZATION 1: Six-panel site-level trend map ───────────────────────
  message("\nGenerating: site-level trend maps...")

  site_level_sp <- site_level %>%
    dplyr::left_join(spatial_pts, by = "site_no") %>%
    dplyr::filter(!is.na(lon), !is.na(lat))

  slope_colors <- c(positive = "#b2182b", negative = "#2166ac", zero = "#888888")

  make_trend_map <- function(ftype, metric_val) {
    df <- site_level_sp %>%
      dplyr::filter(flood_type == ftype, metric == metric_val, !is.na(sen_slope)) %>%
      dplyr::mutate(
        slope_dir = dplyr::case_when(
          sen_slope > 0 ~ "positive",
          sen_slope < 0 ~ "negative",
          TRUE          ~ "zero"
        ),
        sig_fill  = ifelse(significant, "significant", "not significant")
      )
    units_str <- if (metric_val == "count") "events yr⁻¹" else "fraction yr⁻¹"
    ggplot2::ggplot() +
      ggplot2::geom_polygon(
        data = co_poly,
        ggplot2::aes(x = long, y = lat, group = group),
        fill = "grey93", colour = "grey55", linewidth = 0.3
      ) +
      ggplot2::geom_point(
        data = dplyr::filter(df, sig_fill == "not significant"),
        ggplot2::aes(x = lon, y = lat, colour = slope_dir),
        shape = 0, size = 2, stroke = 0.9, alpha = 0.75
      ) +
      ggplot2::geom_point(
        data = dplyr::filter(df, sig_fill == "significant"),
        ggplot2::aes(x = lon, y = lat, colour = slope_dir),
        shape = 15, size = 2.5
      ) +
      ggplot2::scale_colour_manual(
        values = slope_colors,
        labels = c(positive = "Positive slope", negative = "Negative slope",
                   zero = "≈Zero"),
        name   = "Sen slope\ndirection",
        drop   = FALSE
      ) +
      ggplot2::coord_fixed(
        ratio = 1.3,
        xlim  = range(co_poly$long, na.rm = TRUE),
        ylim  = range(co_poly$lat,  na.rm = TRUE)
      ) +
      ggplot2::labs(
        title    = sprintf("%s — %s (%s)",
                           cfg$type_labels[ftype], metric_val, units_str),
        subtitle = sprintf("n = %d eligible sites; ■ = p < %.2f, □ = not sig.",
                           nrow(df), MK_ALPHA),
        x = NULL, y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(
        legend.position  = "bottom",
        legend.key.size  = ggplot2::unit(0.4, "cm"),
        axis.text        = ggplot2::element_blank(),
        axis.ticks       = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(colour = "grey88", linewidth = 0.2)
      )
  }

  ncols    <- length(cfg$flood_types)
  n_panels <- ncols * 2  # count + proportion rows
  map_panels <- c(
    purrr::map(cfg$flood_types, function(ft) make_trend_map(ft, "count")),
    purrr::map(cfg$flood_types, function(ft) make_trend_map(ft, "proportion"))
  )

  map_grid <- patchwork::wrap_plots(map_panels, ncol = ncols, nrow = 2,
                                    guides = "collect") +
    patchwork::plot_annotation(
      title    = sprintf("Site-level Mann-Kendall trend results — %s", cfg$method_label),
      subtitle = "Row 1: count (events/yr);  Row 2: proportion (fraction/yr)",
      theme    = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 10, colour = "grey40")
      )
    ) &
    ggplot2::theme(legend.position = "bottom")

  ggplot2::ggsave(
    file.path(method_diag_dir, "trend_map_panels.png"),
    map_grid,
    width  = 6 * ncols, height = 10, dpi = 150
  )
  message("  Saved: trend_map_panels.png")

  ## ── VISUALIZATION 2: Aggregate summary table ──────────────────────────────
  message("Generating: aggregate summary table figure...")

  type_lvls <- rev(unname(cfg$type_labels[cfg$flood_types]))

  agg_tbl_df <- aggregate %>%
    dplyr::mutate(
      type_label   = factor(cfg$type_labels[flood_type], levels = type_lvls),
      metric_label = factor(metric, levels = c("count", "proportion")),
      pop_label    = sprintf("%s\n(n = %d)", population, n_eligible),
      cell_text    = sprintf("%d (%.0f%%)\nbinom p = %.3f",
                             n_sig, pct_sig, binom_p_any_sig),
      neg_log_p    = pmax(-log10(pmax(binom_p_any_sig, 1e-10)), 0)
    )

  sig_threshold <- -log10(MK_ALPHA)

  p_agg_table <- ggplot2::ggplot(
    agg_tbl_df,
    ggplot2::aes(x = metric_label, y = type_label, fill = neg_log_p)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 1.5) +
    ggplot2::geom_text(ggplot2::aes(label = cell_text), size = 3.2, lineheight = 1.1) +
    ggplot2::facet_wrap(~pop_label, ncol = 2) +
    ggplot2::scale_fill_gradientn(
      colours  = c("white", "#ffffcc", "#feb24c", "#f03b20", "#bd0026"),
      values   = scales::rescale(c(0, sig_threshold * 0.5, sig_threshold,
                                   sig_threshold * 2, sig_threshold * 4)),
      na.value = "grey80",
      name     = expression(-log[10](italic(p)))
    ) +
    ggplot2::labs(
      title    = sprintf("Aggregate MK evidence by flood type — %s", cfg$method_label),
      subtitle = sprintf(
        "Fill: −log₁₀(binomial p) vs. %.0f%% null; text: # significant sites (%%)",
        100 * MK_ALPHA
      ),
      x = "Metric", y = "Flood type"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = "right",
                   strip.text = ggplot2::element_text(face = "bold"))

  ggplot2::ggsave(
    file.path(method_diag_dir, "aggregate_trend_summary_table.png"),
    p_agg_table, width = 12, height = 6, dpi = 150
  )
  message("  Saved: aggregate_trend_summary_table.png")

  ## ── VISUALIZATION 3: Top-6 proportion time series (last flood type) ────────
  # For 3-class methods this is rain_on_snow; for 2-class it is non_rainfall.
  focus_type  <- tail(cfg$flood_types, 1)
  focus_label <- cfg$type_labels[focus_type]
  focus_color <- cfg$type_colors[focus_type]

  message(sprintf("Generating: top-6 %s proportion time series...", focus_label))

  top_sites <- site_level %>%
    dplyr::filter(flood_type == focus_type, metric == "proportion",
                  !is.na(sen_slope), !is.na(tau)) %>%
    dplyr::slice_max(abs(sen_slope), n = 6, with_ties = FALSE) %>%
    dplyr::pull(site_no)

  n_top <- length(top_sites)
  if (n_top == 0L) {
    message(sprintf("  NOTE: no eligible %s-proportion sites — skipping Viz 3.", focus_label))
  } else {
    if (n_top < 6L)
      message(sprintf("  NOTE: only %d eligible %s-proportion sites; showing all.",
                      n_top, focus_label))

    trend_meta <- site_level %>%
      dplyr::filter(flood_type == focus_type, metric == "proportion",
                    site_no %in% top_sites) %>%
      dplyr::select(site_no, tau, mk_pvalue, sen_slope)

    ts_df <- annual_type_counts %>%
      dplyr::filter(site_no %in% top_sites, flood_type == focus_type) %>%
      dplyr::left_join(trend_meta, by = "site_no") %>%
      dplyr::mutate(
        panel_label = sprintf("%s\nτ = %.2f, slope = %.4f/yr, p = %.3f",
                              site_no, tau, sen_slope, mk_pvalue)
      )

    sen_lines <- ts_df %>%
      dplyr::filter(!is.na(proportion)) %>%
      dplyr::group_by(site_no, panel_label, sen_slope) %>%
      dplyr::summarise(
        wy_min    = min(water_year),
        wy_max    = max(water_year),
        intercept = stats::median(proportion - sen_slope * water_year, na.rm = TRUE),
        .groups   = "drop"
      ) %>%
      dplyr::mutate(y_start = intercept + sen_slope * wy_min,
                    y_end   = intercept + sen_slope * wy_max)

    p_top_ts <- ggplot2::ggplot(ts_df, ggplot2::aes(x = water_year, y = proportion)) +
      ggplot2::geom_col(fill = focus_color, alpha = 0.65, width = 0.8) +
      ggplot2::geom_segment(
        data = sen_lines,
        ggplot2::aes(x = wy_min, xend = wy_max, y = y_start, yend = y_end),
        colour = "grey15", linewidth = 1.1, inherit.aes = FALSE
      ) +
      ggplot2::facet_wrap(~panel_label, ncol = 3, scales = "free_y") +
      ggplot2::scale_x_continuous(breaks = seq(1990, 2025, 10)) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      ggplot2::labs(
        title    = sprintf("Top %d sites by |Sen slope| — %s proportion (%s)",
                           n_top, focus_label, cfg$method_label),
        subtitle = "Bar: annual proportion; line: Theil-Sen slope (median pairwise)",
        x = "Water year",
        y = sprintf("%s fraction of classified peaks", focus_label)
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(strip.text    = ggplot2::element_text(size = 7),
                     panel.spacing = ggplot2::unit(0.9, "lines"))

    ggplot2::ggsave(
      file.path(method_diag_dir,
                sprintf("%s_proportion_top%d_slopes.png",
                        gsub("_", "", focus_type), n_top)),
      p_top_ts, width = 14, height = 8, dpi = 150
    )
    message(sprintf("  Saved: %s_proportion_top%d_slopes.png", gsub("_", "", focus_type), n_top))
  }

  message(sprintf("\n  Method %s complete.", method_key))

}  # end method loop

## ── Done ─────────────────────────────────────────────────────────────────────
message(sprintf("\n%s\n=== Script 9 complete — all methods processed ===\n%s",
                strrep("=", 60), strrep("=", 60)))
message("  Data outputs (per method):")
for (mk in names(METHOD_CONFIGS)) {
  message(sprintf("    _data/temporal_trends_%s.rds", mk))
  message(sprintf("    _data/annual_type_counts_%s.rds", mk))
}
message("  AM analysis outputs (per method):")
for (mk in names(METHOD_CONFIGS)) {
  message(sprintf("    _data/temporal_trends_am_%s.rds", mk))
  message(sprintf("    _data/annual_am_type_%s.rds", mk))
  message(sprintf("    _data/am_type_shift_%s.rds", mk))
}
message("  Note: am_rain_frac_trends included only in temporal_trends_am_frac_thresh_3class.rds")
message("  Diagnostic figures: _diagnostics/temporal_trends/{method}/")
