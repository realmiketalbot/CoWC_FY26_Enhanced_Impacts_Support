#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Rain-Fraction-Only Classifier vs Wide Geometric
# ============================================================================
# Classifies POT flood peaks using only the marginal distribution of rain_frac,
# with fixed vertical thresholds at RF_SNOW_THRESH and RF_RAIN_THRESH, and
# compares the result against the wide geometric classifier from
# diag_geometric_boundary_comparison.R.
#
# Decision tree (rain_frac-only, 3-class):
#   1. Unclassified : !data_sufficient | is.na(rain_frac)
#   2. Snowmelt     : rain_frac <= RF_SNOW_THRESH
#   3. Rainfall     : rain_frac >= RF_RAIN_THRESH
#   4. ROS          : remainder
#
# Wide geometric boundaries (carried over unchanged):
#   B1 from (0, 0)    to (0.25, swe_max_log)
#   B2 from (0.75, 0) to (1.00, swe_max_log)
#
# Each classifier produces two objects:
#   _raw : pre-demotion labels (all classifiable events) — used for figures
#   _sel : top-n selected per Yu et al. (2022)           — used for statistics
#
# OUTPUTS (to _diagnostics/threshold_evaluation/rainfrac_only/):
#   - decision_space_rainfrac.png       : rain_frac-only decision space
#   - decision_space_comparison.png     : 2-panel comparison (same point cloud)
#   - rainfrac_classifier_stats.csv     : statistics for both classifiers
# ============================================================================

library(tidyverse)
library(patchwork)

# ── Namespace note ────────────────────────────────────────────────────────────
# mclust and MASS can mask purrr::map and dplyr::select. All calls below use
# explicit namespace qualification throughout.

# ============================================================================
# CONFIGURATION
# ============================================================================
peaks_file <- "_data/pot_peaks_classified_adaptive.rds"
out_dir    <- "_diagnostics/threshold_evaluation/rainfrac_only"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Rain-frac-only thresholds
RF_SNOW_THRESH <- 0.05
RF_RAIN_THRESH <- 0.95

# Wide geometric anchor fractions (must match diag_geometric_boundary_comparison.R)
GEO_B1_RF_UPPER <- 0.25
GEO_B2_RF_LOWER <- 0.75

# Carried over from 6b unchanged
MIN_SINGLE_SAMPLE <- 25L
MIN_TYPE_SAMPLE   <- 20L
SWE_LOG_OFFSET    <- 1.0

# ============================================================================
# SHARED AESTHETICS
# ============================================================================
type_colors <- c(
  "rainfall"     = "#b2182b",
  "snowmelt"     = "#2166ac",
  "ROS"          = "#9970ab",
  "unclassified" = "#cccccc"
)

type_labels <- c(
  "rainfall"     = "Rainfall",
  "snowmelt"     = "Snowmelt",
  "ROS"          = "Rain-on-Snow",
  "unclassified" = "Unclassified"
)

base_theme <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "right"
  )

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading peaks_adap...")
peaks_adap <- readRDS(peaks_file)
message(sprintf("  %d peaks across %d sites",
                nrow(peaks_adap), dplyr::n_distinct(peaks_adap$site_no)))

swe_max_log <- log10(max(peaks_adap$swe_mean_mm, na.rm = TRUE) + SWE_LOG_OFFSET)
message(sprintf("  swe_max_log = %.4f  [max SWE = %.1f mm]",
                swe_max_log, max(peaks_adap$swe_mean_mm, na.rm = TRUE)))

# ============================================================================
# YU ET AL. (2022) TOP-N SELECTION
# ============================================================================
water_year <- function(d) {
  yr <- as.integer(format(as.Date(d), "%Y"))
  mo <- as.integer(format(as.Date(d), "%m"))
  ifelse(mo >= 10L, yr + 1L, yr)
}

apply_top_n <- function(peaks, type_col) {
  peaks <- peaks %>%
    dplyr::select(-dplyr::any_of(c("flood_type", "in_single_sample",
                                   "in_type_sample", "water_year"))) %>%
    dplyr::rename(flood_type = !!type_col) %>%
    dplyr::mutate(water_year = water_year(peak_date))

  single_ids <- peaks %>%
    dplyr::filter(as.character(flood_type) != "unclassified") %>%
    dplyr::group_by(site_no, water_year) %>%
    dplyr::slice_max(peak_flow, n = 1, with_ties = FALSE) %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(
      m    = pmax(MIN_SINGLE_SAMPLE, dplyr::n_distinct(water_year)),
      rank = dplyr::row_number(dplyr::desc(peak_flow))
    ) %>%
    dplyr::filter(rank <= m) %>%
    dplyr::ungroup() %>%
    dplyr::select(site_no, peak_id) %>%
    dplyr::mutate(in_single_sample = TRUE)

  type_ids <- peaks %>%
    dplyr::filter(as.character(flood_type) != "unclassified") %>%
    dplyr::group_by(site_no, flood_type, water_year) %>%
    dplyr::slice_max(peak_flow, n = 1, with_ties = FALSE) %>%
    dplyr::group_by(site_no, flood_type) %>%
    dplyr::mutate(
      n    = pmax(MIN_TYPE_SAMPLE, dplyr::n_distinct(water_year)),
      rank = dplyr::row_number(dplyr::desc(peak_flow))
    ) %>%
    dplyr::filter(rank <= n, dplyr::n() >= MIN_TYPE_SAMPLE) %>%
    dplyr::ungroup() %>%
    dplyr::select(site_no, peak_id, flood_type) %>%
    dplyr::mutate(in_type_sample = TRUE)

  retained <- type_ids %>% dplyr::distinct(site_no, flood_type)

  peaks %>%
    dplyr::left_join(single_ids, by = c("site_no", "peak_id")) %>%
    dplyr::left_join(
      type_ids %>% dplyr::select(site_no, peak_id, in_type_sample),
      by = c("site_no", "peak_id")
    ) %>%
    dplyr::mutate(
      in_single_sample = tidyr::replace_na(in_single_sample, FALSE),
      in_type_sample   = tidyr::replace_na(in_type_sample,   FALSE),
      flood_type = dplyr::if_else(
        as.character(flood_type) != "unclassified" &
          !paste(site_no, as.character(flood_type)) %in%
            paste(retained$site_no, as.character(retained$flood_type)),
        "unclassified",
        as.character(flood_type)
      ),
      flood_type = factor(flood_type,
                          levels = c("rainfall", "snowmelt", "ROS", "unclassified"))
    ) %>%
    dplyr::rename(!!type_col := flood_type)
}

# ============================================================================
# STATISTICS SUMMARY (expects _sel objects with in_type_sample column)
# ============================================================================
summarise_classification <- function(peaks, type_col, label) {
  n_total <- nrow(peaks)

  pre_counts <- peaks %>%
    dplyr::count(!!rlang::sym(type_col), name = "n_events") %>%
    dplyr::rename(flood_type = !!type_col)

  pre_viable <- peaks %>%
    dplyr::filter(as.character(!!rlang::sym(type_col)) != "unclassified") %>%
    dplyr::count(site_no, !!rlang::sym(type_col), name = "n_events_site") %>%
    dplyr::group_by(!!rlang::sym(type_col)) %>%
    dplyr::summarise(viable_sites = sum(n_events_site >= MIN_TYPE_SAMPLE),
                     .groups = "drop") %>%
    dplyr::rename(flood_type = !!type_col)

  post_viable <- peaks %>%
    dplyr::filter(in_type_sample) %>%
    dplyr::count(site_no, !!rlang::sym(type_col)) %>%
    dplyr::group_by(!!rlang::sym(type_col)) %>%
    dplyr::summarise(viable_sites_post = dplyr::n_distinct(site_no),
                     .groups = "drop") %>%
    dplyr::rename(flood_type = !!type_col)

  n_single_sites <- peaks %>%
    dplyr::filter(in_single_sample) %>%
    dplyr::count(site_no) %>%
    dplyr::filter(n >= MIN_SINGLE_SAMPLE) %>%
    nrow()

  pre_counts %>%
    dplyr::mutate(pct_events = round(n_events / n_total * 100, 1),
                  config     = label) %>%
    dplyr::left_join(pre_viable,  by = "flood_type") %>%
    dplyr::left_join(post_viable, by = "flood_type") %>%
    dplyr::mutate(
      viable_sites      = tidyr::replace_na(viable_sites,      0L),
      viable_sites_post = tidyr::replace_na(viable_sites_post, 0L),
      n_single_sites    = dplyr::if_else(
        as.character(flood_type) == "rainfall", n_single_sites, NA_integer_
      )
    ) %>%
    dplyr::select(config, flood_type, n_events, pct_events,
                  viable_sites, viable_sites_post, n_single_sites)
}

# ============================================================================
# CLASSIFY: RAIN-FRAC ONLY
# _raw: pre-demotion labels — used for figures
# _sel: top-n selected     — used for statistics only
# ============================================================================
message("\nClassifying: rain_frac-only...")

peaks_rf_raw <- peaks_adap %>%
  dplyr::mutate(
    flood_type_rf = dplyr::case_when(
      !data_sufficient | is.na(rain_frac) ~ "unclassified",
      rain_frac <= RF_SNOW_THRESH          ~ "snowmelt",
      rain_frac >= RF_RAIN_THRESH          ~ "rainfall",
      TRUE                                 ~ "ROS"
    ),
    flood_type_rf = factor(flood_type_rf,
                           levels = c("rainfall", "snowmelt", "ROS", "unclassified"))
  )

peaks_rf_sel <- peaks_rf_raw %>%
  apply_top_n(type_col = "flood_type_rf")

# ============================================================================
# CLASSIFY: WIDE GEOMETRIC
# _raw: pre-demotion labels — used for figures
# _sel: top-n selected     — used for statistics only
# ============================================================================
message("Classifying: wide geometric...")

cfg_wide <- list(
  b1_slope    = swe_max_log / GEO_B1_RF_UPPER,
  b2_rf_lower = GEO_B2_RF_LOWER,
  b2_rf_slope = (1 - GEO_B2_RF_LOWER) / swe_max_log,
  b1_rf_upper = GEO_B1_RF_UPPER,
  swe_max     = swe_max_log
)

peaks_geo_raw <- peaks_adap %>%
  dplyr::mutate(
    flood_type_geo = dplyr::case_when(
      !data_sufficient | is.na(rain_frac)                                   ~ "unclassified",
      log10_swe > cfg_wide$b1_slope * rain_frac                             ~ "snowmelt",
      rain_frac < 0.01                                                       ~ "snowmelt",
      rain_frac > 0.99                                                       ~ "rainfall",
      rain_frac < cfg_wide$b2_rf_lower + cfg_wide$b2_rf_slope * log10_swe   ~ "ROS",
      TRUE                                                                   ~ "rainfall"
    ),
    flood_type_geo = factor(flood_type_geo,
                            levels = c("rainfall", "snowmelt", "ROS", "unclassified"))
  )

peaks_geo_sel <- peaks_geo_raw %>%
  apply_top_n(type_col = "flood_type_geo")

# ============================================================================
# BOUNDARY LINE DATA FRAMES
# ============================================================================
rf_seq <- seq(0, 1, length.out = 300)

# Rain-frac-only: vertical lines
vline_snow <- tibble::tibble(rain_frac = RF_SNOW_THRESH,
                             log10_swe = c(0, swe_max_log))
vline_rain <- tibble::tibble(rain_frac = RF_RAIN_THRESH,
                             log10_swe = c(0, swe_max_log))

# Wide geometric: diagonal lines
b1_geo <- tibble::tibble(
  rain_frac = rf_seq,
  log10_swe = cfg_wide$b1_slope * rf_seq
) %>%
  dplyr::filter(log10_swe >= 0, log10_swe <= swe_max_log,
                rain_frac <= cfg_wide$b1_rf_upper)

b2_geo <- tibble::tibble(
  rain_frac = rf_seq,
  log10_swe = (rf_seq - cfg_wide$b2_rf_lower) / cfg_wide$b2_rf_slope
) %>%
  dplyr::filter(rain_frac >= cfg_wide$b2_rf_lower, rain_frac <= 1,
                log10_swe >= 0, log10_swe <= swe_max_log)

# ============================================================================
# SHARED PLOT SCAFFOLDING
# ============================================================================
shared_scales <- list(
  ggplot2::scale_x_continuous(
    name   = "Antecedent Rainfall Fraction",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ),
  ggplot2::scale_y_continuous(
    name = "Mean Antecedent SWE [log10(mm + 1)]"
  ),
  ggplot2::scale_colour_manual(
    values = type_colors,
    breaks = names(type_labels),
    labels = type_labels,
    name   = "Flood Type"
  ),
  ggplot2::guides(colour = ggplot2::guide_legend(
    override.aes = list(size = 3, alpha = 1, shape = 16)
  ))
)

make_base_layer <- function(peaks, type_col) {
  plot_data <- peaks %>%
    dplyr::filter(!is.na(log10_swe), !is.na(rain_frac)) %>%
    dplyr::mutate(type_chr = as.character(!!rlang::sym(type_col))) %>%
    dplyr::arrange(type_chr == "unclassified")

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = rain_frac, y = log10_swe, colour = type_chr)
  ) +
    ggplot2::geom_point(shape = 16, size = 0.5, alpha = 0.3)
}

boundary_style <- list(colour = "black", linewidth = 0.6, linetype = "dashed")

# ============================================================================
# FIGURE 1: RAIN-FRAC-ONLY DECISION SPACE (pre-demotion)
# ============================================================================
message("\nGenerating Figure 1: rain_frac-only decision space...")

set.seed(42L)
p_rf <- make_base_layer(peaks_rf_raw, "flood_type_rf") +
  ggplot2::geom_line(
    data = vline_snow, ggplot2::aes(x = rain_frac, y = log10_swe),
    colour = boundary_style$colour, linewidth = boundary_style$linewidth,
    linetype = boundary_style$linetype, inherit.aes = FALSE
  ) +
  ggplot2::geom_line(
    data = vline_rain, ggplot2::aes(x = rain_frac, y = log10_swe),
    colour = boundary_style$colour, linewidth = boundary_style$linewidth,
    linetype = boundary_style$linetype, inherit.aes = FALSE
  ) +
  ggplot2::annotate(
    "text", x = RF_SNOW_THRESH + 0.02, y = swe_max_log * 0.97,
    label = sprintf("B1\n(rf=%.2f)", RF_SNOW_THRESH),
    colour = "grey30", size = 3, hjust = 0
  ) +
  ggplot2::annotate(
    "text", x = RF_RAIN_THRESH - 0.02, y = swe_max_log * 0.97,
    label = sprintf("B2\n(rf=%.2f)", RF_RAIN_THRESH),
    colour = "grey30", size = 3, hjust = 1
  ) +
  shared_scales +
  ggplot2::labs(subtitle = sprintf(
    "Rain-fraction only  |  Snowmelt: rf \u2264 %.2f  |  Rainfall: rf \u2265 %.2f  |  ROS: remainder",
    RF_SNOW_THRESH, RF_RAIN_THRESH
  )) +
  base_theme

ggplot2::ggsave(file.path(out_dir, "decision_space_rainfrac.png"),
                plot = p_rf, width = 6.5, height = 4, units = "in", dpi = 150)
message("  Saved: decision_space_rainfrac.png")

# ============================================================================
# FIGURE 2: CONFUSION MATRIX — rain_frac-only vs wide geometric (raw labels)
# Rows = rain_frac-only, Columns = wide geometric
# Cells annotated with event count and row-percent (% of rf-only class)
# ============================================================================
message("\nGenerating Figure 2: classifier agreement heatmap...")

# Join raw labels on peak_id; restrict to events classifiable by both.
# Deduplicate on peak_id first — peaks_adap carries one row per
# (peak, window-day), not one row per peak; joining without deduplication
# produces a many-to-many cross product inflating counts by ~100x.
confusion_df <- peaks_rf_raw %>%
  dplyr::distinct(site_no, peak_id, .keep_all = TRUE) %>%
  dplyr::select(site_no, peak_id, flood_type_rf) %>%
  dplyr::inner_join(
    peaks_geo_raw %>%
      dplyr::distinct(site_no, peak_id, .keep_all = TRUE) %>%
      dplyr::select(site_no, peak_id, flood_type_geo),
    by = c("site_no", "peak_id")
  ) %>%
  dplyr::filter(
    as.character(flood_type_rf)  != "unclassified",
    as.character(flood_type_geo) != "unclassified"
  ) %>%
  dplyr::count(flood_type_rf, flood_type_geo, name = "n") %>%
  dplyr::group_by(flood_type_rf) %>%
  dplyr::mutate(
    row_total  = sum(n),
    row_pct    = round(n / row_total * 100, 1)
  ) %>%
  dplyr::ungroup() %>%
  # Ensure all 9 cells present (fill absent combinations with 0)
  tidyr::complete(
    flood_type_rf  = factor(c("rainfall", "snowmelt", "ROS"),
                            levels = c("snowmelt", "ROS", "rainfall")),
    flood_type_geo = factor(c("rainfall", "snowmelt", "ROS"),
                            levels = c("snowmelt", "ROS", "rainfall")),
    fill = list(n = 0L, row_total = NA_integer_, row_pct = 0)
  ) %>%
  dplyr::mutate(
    # Diagonal cells (agreement) vs off-diagonal (disagreement)
    is_agree   = as.character(flood_type_rf) == as.character(flood_type_geo),
    cell_label = dplyr::if_else(n > 0,
                                sprintf("%d\n(%.1f%%)", n, row_pct),
                                "0")
  )

# Axis label ordering: snowmelt → ROS → rainfall (low to high rain_frac)
type_order <- c("snowmelt", "ROS", "rainfall")
type_display <- c("snowmelt" = "Snowmelt",
                  "ROS"      = "Rain-on-Snow",
                  "rainfall" = "Rainfall")

confusion_df <- confusion_df %>%
  dplyr::mutate(
    flood_type_rf  = factor(as.character(flood_type_rf),  levels = type_order),
    flood_type_geo = factor(as.character(flood_type_geo), levels = type_order)
  )

p_confusion <- ggplot2::ggplot(
  confusion_df,
  ggplot2::aes(x = flood_type_geo, y = flood_type_rf, fill = n)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
  ggplot2::geom_text(
    ggplot2::aes(label = cell_label, colour = is_agree),
    size = 3.8, lineheight = 1.2
  ) +
  ggplot2::scale_fill_gradient(
    low  = "white",
    high = "#4d4d4d",
    name = "Events",
    na.value = "white"
  ) +
  ggplot2::scale_colour_manual(
    values = c("TRUE" = "#1a9641", "FALSE" = "#d7191c"),
    guide  = "none"
  ) +
  ggplot2::scale_x_discrete(
    name   = "Wide geometric classifier",
    labels = type_display
  ) +
  ggplot2::scale_y_discrete(
    name   = "Rain-fraction-only classifier",
    labels = type_display
  ) +
  ggplot2::labs(
    subtitle = sprintf(
      "Agreement (green) / disagreement (red)  |  rf-only: snow\u2264%.2f, rain\u2265%.2f  |  geo: B1\u21920.25, B2\u21900.75",
      RF_SNOW_THRESH, RF_RAIN_THRESH
    ),
    caption  = "Cell %: row proportion (share of rf-only class assigned to each geometric class)"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid   = ggplot2::element_blank(),
    axis.ticks   = ggplot2::element_blank(),
    legend.position = "right"
  )

ggplot2::ggsave(file.path(out_dir, "classifier_agreement.png"),
                plot = p_confusion, width = 6.5, height = 5,
                units = "in", dpi = 150)
message("  Saved: classifier_agreement.png")

# ============================================================================
# STATISTICS (uses _sel objects — top-n applied, post-demotion counts)
# ============================================================================
message("\nComputing summary statistics...")

stats_all <- dplyr::bind_rows(
  summarise_classification(
    peaks_rf_sel,  "flood_type_rf",
    sprintf("Rain-frac only (snow<=%.2f, rain>=%.2f)", RF_SNOW_THRESH, RF_RAIN_THRESH)
  ),
  summarise_classification(
    peaks_geo_sel, "flood_type_geo",
    sprintf("Wide geometric (B1->%.2f, B2<-%.2f)", GEO_B1_RF_UPPER, GEO_B2_RF_LOWER)
  )
)

# ============================================================================
# PRINT STATISTICS
# ============================================================================
message("\n", strrep("=", 72))
message("CLASSIFICATION STATISTICS  (top-n selected, post-demotion)")
message(strrep("=", 72))

for (cfg_label in unique(stats_all$config)) {
  message(sprintf("\n--- %s ---", cfg_label))
  sub <- stats_all %>% dplyr::filter(config == cfg_label)
  message(sprintf("  %-14s  %6s  %6s  %14s  %15s",
                  "Type", "Events", "  Pct", "Pre-sel sites", "Post-sel sites"))
  message(sprintf("  %-14s  %6s  %6s  %14s  %15s",
                  "----", "------", "-----", "  (n>=20 raw)", "  (n>=20 top-n)"))
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    message(sprintf("  %-14s  %6d  %5.1f%%  %14d  %15d",
                    r$flood_type, r$n_events, r$pct_events,
                    r$viable_sites, r$viable_sites_post))
  }
  n_sing <- sub$n_single_sites[!is.na(sub$n_single_sites)]
  if (length(n_sing) > 0)
    message(sprintf("  Sites >= %d single-sample peaks: %d",
                    MIN_SINGLE_SAMPLE, n_sing[1]))
}

message(sprintf("\n%s", strrep("=", 72)))

# ============================================================================
# SAVE STATISTICS CSV
# ============================================================================
readr::write_csv(stats_all, file.path(out_dir, "rainfrac_classifier_stats.csv"))
message(sprintf("\nSaved: %s/rainfrac_classifier_stats.csv", out_dir))
message("Done.")