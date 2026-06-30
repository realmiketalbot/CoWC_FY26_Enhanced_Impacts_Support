#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Geometric Boundary Comparison
# ============================================================================
# Evaluates two alternative geometric boundary configurations for the adaptive
# (rain_frac x log10(SWE+1)) classifier against the existing GMM-derived
# boundaries, producing one decision-space figure and one statistics block
# per configuration (three total).
#
# The two geometric configurations replace the GMM-derived B1/B2 lines with
# straight lines anchored at the data's site-pooled SWE maximum:
#
#   Narrow-ROS  (geo_narrow):  B1 from (0,   0) to (0.50, swe_max_log)
#                              B2 from (0.50, 0) to (1.00, swe_max_log)
#
#   Wide-ROS    (geo_wide):    B1 from (0,   0) to (0.25, swe_max_log)
#                              B2 from (0.75, 0) to (1.00, swe_max_log)
#
# In both cases swe_max_log = log10(max(swe_mean_mm, na.rm=TRUE) + 1) across
# all classified peaks, making the upper anchor data-driven rather than fixed.
#
# Both configurations retain the same sequential decision tree as 6b:
#   1. Unclassified : !data_sufficient | is.na(rain_frac)
#   2. Snowmelt     : log10(SWE+1) above B1
#   3. ROS          : rain_frac below B2  (subject to hard guards)
#   4. Rainfall     : remainder
#
# The hard guards (ROS_RAIN_FRAC_MIN, ROS_RAIN_FRAC_MAX) and Yu et al. (2022)
# top-n selection (MIN_TYPE_SAMPLE, MIN_SINGLE_SAMPLE) are carried over
# unchanged from 6b to keep the comparison apples-to-apples.
#
# INPUTS:
#   - _data/pot_peaks_classified_adaptive.rds : existing adaptive output
#     (provides pre-computed rain_frac, log10_swe, data_sufficient, and the
#      GMM-derived flood_type_raw for the first figure)
#
# OUTPUTS (to _diagnostics/threshold_evaluation/geometric_boundaries/):
#   - decision_space_gmm.png        : existing GMM boundary space (reproduced)
#   - decision_space_geo_narrow.png : narrow-ROS geometric boundaries
#   - decision_space_geo_wide.png   : wide-ROS geometric boundaries
#   - geometric_boundary_stats.csv  : per-configuration summary statistics
#
# Statistics reported per configuration:
#   - Total classified / unclassified peaks
#   - Per-type (rainfall / snowmelt / ROS) pre-selection event counts and %
#   - Sites meeting MIN_TYPE_SAMPLE for each type (pre- and post-selection)
#   - Sites meeting MIN_SINGLE_SAMPLE
# ============================================================================

library(tidyverse)

# ── Namespace note ────────────────────────────────────────────────────────────
# mclust and MASS can mask purrr::map and dplyr::select. All calls below use
# explicit namespace qualification throughout.

# ============================================================================
# CONFIGURATION
# ============================================================================
peaks_file <- "_data/pot_peaks_classified_adaptive.rds"
out_dir    <- "_diagnostics/threshold_evaluation/geometric_boundaries"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Carry over from 6b unchanged
ROS_RAIN_FRAC_MIN <- 0.01
ROS_RAIN_FRAC_MAX <- 0.99
MIN_SINGLE_SAMPLE <- 25L
MIN_TYPE_SAMPLE   <- 20L
SWE_LOG_OFFSET    <- 1.0

# GMM-derived boundary coefficients (from 6b, k=2.0 on B2_INTERCEPT)
GMM_B1_INTERCEPT <- 0.4973
GMM_B1_SLOPE     <- 5.3433
GMM_B2_INTERCEPT <- 0.2156   # 0.1078 * k=2.0
GMM_B2_SLOPE     <- 0.3503

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading peaks_adap...")
peaks_adap <- readRDS(peaks_file)

message(sprintf("  %d peaks across %d sites",
                nrow(peaks_adap), dplyr::n_distinct(peaks_adap$site_no)))

# ============================================================================
# DERIVE swe_max_log FROM EVENT DATA
# ============================================================================
# Upper anchor for geometric boundaries: log10(max event-SWE + 1) pooled
# across all peaks with valid swe_mean_mm. This keeps the anchor consistent
# with the feature space used in classification.
swe_max_log <- log10(max(peaks_adap$swe_mean_mm, na.rm = TRUE) + SWE_LOG_OFFSET)
message(sprintf("  swe_max_log (pooled) = %.4f  [max SWE = %.1f mm]",
                swe_max_log,
                max(peaks_adap$swe_mean_mm, na.rm = TRUE)))

# ============================================================================
# BOUNDARY GEOMETRY HELPERS
# ============================================================================

# Each geometric configuration is defined by two anchor rain_frac values:
#   b1_rf_upper : rain_frac at which B1 reaches swe_max_log  (B1 left bound)
#   b2_rf_lower : rain_frac at which B2 originates at SWE=0  (B2 right bound)
#
# B1: line from (0, 0) to (b1_rf_upper, swe_max_log)
#   log10_swe = (swe_max_log / b1_rf_upper) * rain_frac
#   A point is Snowmelt if log10_swe > B1(rain_frac)
#
# B2: line from (b2_rf_lower, 0) to (1, swe_max_log)
#   log10_swe = swe_max_log / (1 - b2_rf_lower) * (rain_frac - b2_rf_lower)
#   A point is ROS if rain_frac < B2_inv(log10_swe)
#   i.e., rain_frac < b2_rf_lower + (1 - b2_rf_lower) / swe_max_log * log10_swe

make_geo_config <- function(b1_rf_upper, b2_rf_lower, swe_max) {
  list(
    b1_slope      = swe_max / b1_rf_upper,
    b2_rf_lower   = b2_rf_lower,
    b2_rf_slope   = (1 - b2_rf_lower) / swe_max,   # d(rain_frac)/d(log10_swe)
    b1_rf_upper   = b1_rf_upper,
    swe_max       = swe_max
  )
}

# B1 threshold: is the point above B1 (i.e., in snowmelt zone)?
above_b1 <- function(rain_frac, log10_swe, cfg) {
  log10_swe > cfg$b1_slope * rain_frac
}

# B2 threshold: is the point below B2 (i.e., in ROS zone)?
below_b2 <- function(rain_frac, log10_swe, cfg) {
  rain_frac < cfg$b2_rf_lower + cfg$b2_rf_slope * log10_swe
}

# ============================================================================
# CLASSIFY UNDER A GIVEN GEOMETRIC CONFIG
# ============================================================================
classify_geometric <- function(peaks, cfg) {
  peaks %>%
    dplyr::mutate(
      flood_type_geo = dplyr::case_when(
        !data_sufficient | is.na(rain_frac)                          ~ "unclassified",
        above_b1(rain_frac, log10_swe, cfg)                         ~ "snowmelt",
        rain_frac < ROS_RAIN_FRAC_MIN                               ~ "snowmelt",
        rain_frac > ROS_RAIN_FRAC_MAX                               ~ "rainfall",
        below_b2(rain_frac, log10_swe, cfg)                         ~ "ROS",
        TRUE                                                         ~ "rainfall"
      ),
      flood_type_geo = factor(flood_type_geo,
                              levels = c("rainfall", "snowmelt", "ROS", "unclassified"))
    )
}

# ============================================================================
# YU ET AL. (2022) TOP-N SELECTION (mirrored from 6b)
# ============================================================================
water_year <- function(d) {
  yr <- as.integer(format(as.Date(d), "%Y"))
  mo <- as.integer(format(as.Date(d), "%m"))
  ifelse(mo >= 10L, yr + 1L, yr)
}

apply_top_n <- function(peaks, type_col = "flood_type_geo") {
  # Drop any pre-existing flood_type and selection-flag columns carried over
  # from 6b to avoid rename/join collisions, then rename type_col to flood_type.
  peaks <- peaks %>%
    dplyr::select(-dplyr::any_of(c("flood_type", "in_single_sample",
                                   "in_type_sample", "water_year"))) %>%
    dplyr::rename(flood_type = !!type_col)

  peaks <- peaks %>%
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

  peaks <- peaks %>%
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
    )

  # Rename back to the original type_col name
  peaks %>% dplyr::rename(!!type_col := flood_type)
}

# ============================================================================
# STATISTICS SUMMARY
# ============================================================================
summarise_classification <- function(peaks, type_col, label) {
  ft <- dplyr::pull(peaks, !!type_col)

  pre_counts <- table(ft)
  n_total    <- nrow(peaks)

  # Pre-selection viable sites (raw counts before top-n)
  pre_viable <- peaks %>%
    dplyr::filter(as.character(!!rlang::sym(type_col)) != "unclassified") %>%
    dplyr::count(site_no, !!rlang::sym(type_col), name = "n_events") %>%
    dplyr::group_by(!!rlang::sym(type_col)) %>%
    dplyr::summarise(
      viable_sites = sum(n_events >= MIN_TYPE_SAMPLE),
      .groups = "drop"
    ) %>%
    dplyr::rename(flood_type = !!type_col)

  # Post-selection: sites where in_type_sample is TRUE for each type
  post_viable <- peaks %>%
    dplyr::filter(in_type_sample) %>%
    dplyr::count(site_no, !!rlang::sym(type_col), name = "n_sel") %>%
    dplyr::group_by(!!rlang::sym(type_col)) %>%
    dplyr::summarise(
      viable_sites_post = dplyr::n_distinct(site_no),
      .groups = "drop"
    ) %>%
    dplyr::rename(flood_type = !!type_col)

  n_single_sites <- peaks %>%
    dplyr::filter(in_single_sample) %>%
    dplyr::count(site_no) %>%
    dplyr::filter(n >= MIN_SINGLE_SAMPLE) %>%
    nrow()

  tibble::tibble(
    config        = label,
    flood_type    = names(pre_counts),
    n_events      = as.integer(pre_counts),
    pct_events    = round(as.integer(pre_counts) / n_total * 100, 1)
  ) %>%
    dplyr::left_join(pre_viable,  by = "flood_type") %>%
    dplyr::left_join(post_viable, by = "flood_type") %>%
    dplyr::mutate(
      viable_sites      = tidyr::replace_na(viable_sites,      0L),
      viable_sites_post = tidyr::replace_na(viable_sites_post, 0L),
      n_single_sites    = dplyr::if_else(flood_type == "rainfall",
                                         n_single_sites, NA_integer_)
    )
}

# ============================================================================
# FIGURE FUNCTION (matches REPORT.Rmd fig-adaptive-space aesthetic exactly)
# ============================================================================
type_colors <- c(
  "rainfall"     = "#b2182b",
  "snowmelt"     = "#2166ac",
  "ROS"          = "#9970ab",
  "unclassified" = "#cccccc"
)

make_decision_space_plot <- function(peaks, type_col, b1_df, b2_df,
                                     title_label) {
  plot_data <- peaks %>%
    dplyr::mutate(
      flood_type_plot = as.character(!!rlang::sym(type_col))
    ) %>%
    dplyr::filter(!is.na(log10_swe), !is.na(rain_frac)) %>%
    dplyr::group_by(flood_type_plot) %>%
    dplyr::slice_sample(prop = 1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(flood_type_plot != "unclassified")

  ggplot2::ggplot(plot_data,
                  ggplot2::aes(x = rain_frac, y = log10_swe,
                               colour = flood_type_plot)) +
    ggplot2::geom_point(shape = 16, size = 0.5, alpha = 0.3) +
    ggplot2::geom_line(
      data = b1_df,
      ggplot2::aes(x = rain_frac, y = log10_swe),
      colour = "black", linewidth = 0.6, linetype = "dashed",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_line(
      data = b2_df,
      ggplot2::aes(x = rain_frac, y = log10_swe),
      colour = "black", linewidth = 0.6, linetype = "dashed",
      inherit.aes = FALSE
    ) +
    ggplot2::annotate(
      "text",
      x = dplyr::slice_tail(b1_df, n = 1)$rain_frac - 0.03,
      y = dplyr::slice_tail(b1_df, n = 1)$log10_swe,
      label = "B1", colour = "grey30", size = 3, hjust = 1
    ) +
    ggplot2::annotate(
      "text",
      x = dplyr::slice_tail(b2_df, n = 1)$rain_frac - 0.03,
      y = dplyr::slice_tail(b2_df, n = 1)$log10_swe,
      label = "B2", colour = "grey30", size = 3, hjust = 1
    ) +
    ggplot2::scale_x_continuous(
      name   = "Antecedent Rainfall Fraction",
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25)
    ) +
    ggplot2::scale_y_continuous(
      name = "Mean Antecedent SWE [log10(mm + 1)]"
    ) +
    ggplot2::scale_colour_manual(
      values = type_colors,
      breaks = c("rainfall", "snowmelt", "ROS", "unclassified"),
      labels = c("rainfall"     = "Rainfall",
                 "snowmelt"     = "Snowmelt",
                 "ROS"          = "Rain-on-Snow",
                 "unclassified" = "Unclassified"),
      name = "Flood Type"
    ) +
    ggplot2::labs(subtitle = title_label) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "right"
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(size = 3, alpha = 1, shape = 16)
    ))
}

# ============================================================================
# BUILD BOUNDARY LINE DATA FRAMES FOR EACH CONFIGURATION
# ============================================================================
rf_seq    <- seq(0, 1, length.out = 300)
swe_clip  <- stats::quantile(peaks_adap$log10_swe, 0.99, na.rm = TRUE)

# --- GMM boundaries (from 6b) -----------------------------------------------
b1_gmm <- tibble::tibble(
  rain_frac = rf_seq,
  log10_swe = GMM_B1_INTERCEPT + GMM_B1_SLOPE * rf_seq
) %>%
  dplyr::filter(log10_swe >= 0, log10_swe <= swe_clip)

b2_gmm <- tibble::tibble(
  rain_frac = rf_seq,
  log10_swe = (rf_seq - GMM_B2_INTERCEPT) / GMM_B2_SLOPE
) %>%
  dplyr::filter(rain_frac >= GMM_B2_INTERCEPT, rain_frac <= 1,
                log10_swe >= 0, log10_swe <= swe_clip)

# --- Narrow-ROS geometric boundaries ----------------------------------------
cfg_narrow <- make_geo_config(b1_rf_upper = 0.50,
                              b2_rf_lower = 0.50,
                              swe_max     = swe_max_log)

b1_narrow <- tibble::tibble(
  rain_frac = rf_seq,
  log10_swe = cfg_narrow$b1_slope * rf_seq
) %>%
  dplyr::filter(log10_swe >= 0, log10_swe <= swe_max_log,
                rain_frac <= cfg_narrow$b1_rf_upper)

b2_narrow <- tibble::tibble(
  rain_frac = rf_seq,
  log10_swe = (rf_seq - cfg_narrow$b2_rf_lower) / cfg_narrow$b2_rf_slope
) %>%
  dplyr::filter(rain_frac >= cfg_narrow$b2_rf_lower, rain_frac <= 1,
                log10_swe >= 0, log10_swe <= swe_max_log)

# --- Wide-ROS geometric boundaries ------------------------------------------
cfg_wide <- make_geo_config(b1_rf_upper = 0.25,
                            b2_rf_lower = 0.75,
                            swe_max     = swe_max_log)

b1_wide <- tibble::tibble(
  rain_frac = rf_seq,
  log10_swe = cfg_wide$b1_slope * rf_seq
) %>%
  dplyr::filter(log10_swe >= 0, log10_swe <= swe_max_log,
                rain_frac <= cfg_wide$b1_rf_upper)

b2_wide <- tibble::tibble(
  rain_frac = rf_seq,
  log10_swe = (rf_seq - cfg_wide$b2_rf_lower) / cfg_wide$b2_rf_slope
) %>%
  dplyr::filter(rain_frac >= cfg_wide$b2_rf_lower, rain_frac <= 1,
                log10_swe >= 0, log10_swe <= swe_max_log)

# ============================================================================
# FIGURE 1: GMM BOUNDARIES (reproduced from REPORT.Rmd)
# Uses flood_type_raw from peaks_adap (pre-demotion labels from 6b)
# ============================================================================
message("\nGenerating Figure 1: GMM decision space...")

set.seed(42L)
p_gmm <- make_decision_space_plot(
  peaks      = peaks_adap %>% dplyr::mutate(flood_type_raw_chr = as.character(flood_type_raw)),
  type_col   = "flood_type_raw",
  b1_df      = b1_gmm,
  b2_df      = b2_gmm,
  title_label = sprintf(
    "GMM-derived boundaries  |  B1: intercept=%.4f, slope=%.4f  |  B2: intercept=%.4f, slope=%.4f",
    GMM_B1_INTERCEPT, GMM_B1_SLOPE, GMM_B2_INTERCEPT, GMM_B2_SLOPE
  )
)

ggplot2::ggsave(
  file.path(out_dir, "decision_space_gmm.png"),
  plot = p_gmm, width = 6.5, height = 4, units = "in", dpi = 150
)
message("  Saved: decision_space_gmm.png")

# ============================================================================
# CLASSIFY UNDER GEOMETRIC CONFIGURATIONS + APPLY TOP-N
# ============================================================================

# --- Narrow ------------------------------------------------------------------
message("\nClassifying: geo_narrow...")
peaks_narrow <- peaks_adap %>%
  classify_geometric(cfg = cfg_narrow) %>%
  apply_top_n(type_col = "flood_type_geo")

# --- Wide --------------------------------------------------------------------
message("Classifying: geo_wide...")
peaks_wide <- peaks_adap %>%
  classify_geometric(cfg = cfg_wide) %>%
  apply_top_n(type_col = "flood_type_geo")

# ============================================================================
# FIGURE 2: NARROW-ROS GEOMETRIC BOUNDARIES
# ============================================================================
message("\nGenerating Figure 2: geo_narrow decision space...")

set.seed(42L)
p_narrow <- make_decision_space_plot(
  peaks       = peaks_narrow,
  type_col    = "flood_type_geo",
  b1_df       = b1_narrow,
  b2_df       = b2_narrow,
  title_label = sprintf(
    "Geometric (narrow ROS)  |  B1: (0,0)→(0.50, %.2f)  |  B2: (0.50,0)→(1.00, %.2f)",
    swe_max_log, swe_max_log
  )
)

ggplot2::ggsave(
  file.path(out_dir, "decision_space_geo_narrow.png"),
  plot = p_narrow, width = 6.5, height = 4, units = "in", dpi = 150
)
message("  Saved: decision_space_geo_narrow.png")

# ============================================================================
# FIGURE 3: WIDE-ROS GEOMETRIC BOUNDARIES
# ============================================================================
message("\nGenerating Figure 3: geo_wide decision space...")

set.seed(42L)
p_wide <- make_decision_space_plot(
  peaks       = peaks_wide,
  type_col    = "flood_type_geo",
  b1_df       = b1_wide,
  b2_df       = b2_wide,
  title_label = sprintf(
    "Geometric (wide ROS)    |  B1: (0,0)→(0.25, %.2f)  |  B2: (0.75,0)→(1.00, %.2f)",
    swe_max_log, swe_max_log
  )
)

ggplot2::ggsave(
  file.path(out_dir, "decision_space_geo_wide.png"),
  plot = p_wide, width = 6.5, height = 4, units = "in", dpi = 150
)
message("  Saved: decision_space_geo_wide.png")

# ============================================================================
# STATISTICS
# ============================================================================
message("\nComputing summary statistics...")

# GMM stats use flood_type_raw (pre-demotion) so they match the figure.
# We reconstruct in_single_sample and in_type_sample from the existing output
# so that post-selection viable site counts are correct for the GMM config.
peaks_gmm_stats <- peaks_adap %>%
  dplyr::select(-dplyr::any_of(c("flood_type", "in_single_sample",
                                  "in_type_sample", "water_year"))) %>%
  dplyr::rename(flood_type_geo = flood_type_raw) %>%
  dplyr::mutate(flood_type_geo = factor(as.character(flood_type_geo),
                                        levels = c("rainfall", "snowmelt",
                                                   "ROS", "unclassified")))

# Re-apply top-n against flood_type_raw to get selection flags for GMM config
peaks_gmm_stats <- peaks_gmm_stats %>%
  apply_top_n(type_col = "flood_type_geo")

stats_all <- dplyr::bind_rows(
  summarise_classification(peaks_gmm_stats, "flood_type_geo", "GMM (k=2.0)"),
  summarise_classification(peaks_narrow,    "flood_type_geo", "Geometric narrow (B1→0.50, B2←0.50)"),
  summarise_classification(peaks_wide,      "flood_type_geo", "Geometric wide   (B1→0.25, B2←0.75)")
)

# ============================================================================
# PRINT STATISTICS TO CONSOLE
# ============================================================================
message("\n", strrep("=", 72))
message("CLASSIFICATION STATISTICS — PRE-DEMOTION EVENT COUNTS")
message(strrep("=", 72))

for (cfg_label in unique(stats_all$config)) {
  message(sprintf("\n--- %s ---", cfg_label))
  sub <- stats_all %>% dplyr::filter(config == cfg_label)
  message(sprintf("  %-14s  %6s  %6s  %14s  %15s",
                  "Type", "Events", "  Pct",
                  "Pre-sel sites", "Post-sel sites"))
  message(sprintf("  %-14s  %6s  %6s  %14s  %15s",
                  "----", "------", "-----",
                  "  (n>=20 raw)", "  (n>=20 top-n)"))
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    message(sprintf("  %-14s  %6d  %5.1f%%  %14d  %15d",
                    r$flood_type, r$n_events, r$pct_events,
                    r$viable_sites, r$viable_sites_post))
  }
  n_sing <- sub$n_single_sites[!is.na(sub$n_single_sites)]
  if (length(n_sing) > 0) {
    message(sprintf("  Sites >= %d single-sample peaks: %d",
                    MIN_SINGLE_SAMPLE, n_sing[1]))
  }
}

message(sprintf("\n%s", strrep("=", 72)))

# ============================================================================
# SAVE STATISTICS CSV
# ============================================================================
readr::write_csv(stats_all, file.path(out_dir, "geometric_boundary_stats.csv"))
message(sprintf("\nSaved statistics: %s/geometric_boundary_stats.csv", out_dir))
message("Done.")