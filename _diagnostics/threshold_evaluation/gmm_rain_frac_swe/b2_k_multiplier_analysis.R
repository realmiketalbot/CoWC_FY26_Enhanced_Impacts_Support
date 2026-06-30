#!/usr/bin/env Rscript

# ============================================================================
# b2_k_multiplier_analysis.R
#
# Sweeps a scalar multiplier k applied to B2_INTERCEPT in the adaptive
# classification decision tree. The effective B2 boundary becomes:
#
#   rain_frac < (k * B2_INTERCEPT) + B2_SLOPE * log10_swe
#
# A peak is classified as ROS if its rain_frac falls below this boundary
# (having already passed the Snowmelt test). Larger k lowers the boundary,
# widening the ROS region and increasing ROS event counts.
#
# After pipeline fixes (flux floor removal + is.na(rain_frac) guard), viable
# ROS sites dropped from ~37 to ~20. This script finds the smallest k that
# recovers >= 30 viable ROS sites while maintaining >= 50 viable Snowmelt
# sites (MIN_TYPE_SAMPLE = 20 for both).
#
# INPUT:  _data/pot_peaks_classified_adaptive.rds
# OUTPUT: fig_b2_k_viability.png, fig_b2_k_composition.png,
#         fig_b2_k_classification_space.png
# ============================================================================

library(tidyverse)
library(patchwork)
library(scales)

# ============================================================================
# CONFIGURATION
# ============================================================================
PEAKS_FILE      <- "_data/pot_peaks_classified_adaptive.rds"
METADATA_FILE   <- "_data/site_metadata.rds"
OUT_DIR         <- "_diagnostics/threshold_evaluation/gmm_rain_frac_swe"
B1_INTERCEPT    <- 0.4973
B1_SLOPE        <- 5.3433
B2_INTERCEPT    <- 0.1078
B2_SLOPE        <- 0.3503
SWE_LOG_OFFSET  <- 1.0
MIN_TYPE_SAMPLE <- 20L
K_GRID          <- c(1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0)
K_GRID_FINE     <- seq(2.0, 2.5, by = 0.05)

# Viability targets for console flagging
ROS_TARGET      <- 30L
SNOWMELT_TARGET <- 50L

SAVE_WIDTH      <- 10
SAVE_HEIGHT     <- 7
SAVE_DPI        <- 150

# ============================================================================
# SETUP
# ============================================================================
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(p, filename, width = SAVE_WIDTH, height = SAVE_HEIGHT) {
  path <- file.path(OUT_DIR, filename)
  ggsave(path, p, width = width, height = height, dpi = SAVE_DPI)
  message(sprintf("  Saved: %s", path))
  invisible(path)
}

theme_diag <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = NA),
      legend.position  = "bottom"
    )
}

type_colors <- c(
  "ROS"           = "#d7191c",
  "snowmelt"      = "#4dac26",
  "rainfall"      = "#2166ac",
  "unclassified"  = "#d9d9d9"
)

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading classified peaks...")
peaks <- readRDS(PEAKS_FILE)
message(sprintf("  %d peaks across %d sites", nrow(peaks), dplyr::n_distinct(peaks$site_no)))

# Confirm or compute log10_swe
if ("log10_swe" %in% names(peaks)) {
  message("  log10_swe column found — using directly.")
} else {
  message("  log10_swe not found — computing from swe_mean_mm.")
  peaks <- peaks %>%
    dplyr::mutate(
      log10_swe = log10(dplyr::coalesce(swe_mean_mm, 0) + SWE_LOG_OFFSET)
    )
}

# Restrict to peaks with the columns we need
required_cols <- c("site_no", "rain_frac", "log10_swe", "data_sufficient")
missing_cols  <- setdiff(required_cols, names(peaks))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

# Work with a slim table for speed
peaks_slim <- peaks %>%
  dplyr::select(site_no, rain_frac, log10_swe, data_sufficient)

message(sprintf("  Data-sufficient peaks: %d",
                sum(peaks_slim$data_sufficient, na.rm = TRUE)))

# ============================================================================
# K SWEEP
# ============================================================================
message("\nRunning B2 k-multiplier sweep across k = ",
        paste(K_GRID, collapse = ", "), "...")

apply_decision_tree <- function(peaks_df, k) {
  peaks_df %>%
    dplyr::mutate(
      flood_type_k = dplyr::case_when(
        !data_sufficient | is.na(rain_frac)                            ~ "unclassified",
        log10_swe > B1_INTERCEPT + B1_SLOPE * rain_frac               ~ "snowmelt",
        rain_frac < (k * B2_INTERCEPT) + B2_SLOPE * log10_swe         ~ "ROS",
        TRUE                                                           ~ "rainfall"
      )
    )
}

sweep_results <- purrr::map_dfr(K_GRID, function(k) {
  message(sprintf("  k = %.2f ...", k))

  classified <- apply_decision_tree(peaks_slim, k)

  # Per-site viable counts (pre-demotion: count ALL peaks with each label)
  site_counts <- classified %>%
    dplyr::filter(flood_type_k != "unclassified") %>%
    dplyr::count(site_no, flood_type_k) %>%
    tidyr::pivot_wider(names_from = flood_type_k, values_from = n,
                       values_fill = 0L) %>%
    dplyr::mutate(
      ROS      = dplyr::coalesce(ROS,      0L),
      snowmelt = dplyr::coalesce(snowmelt, 0L),
      rainfall = dplyr::coalesce(rainfall, 0L)
    )

  n_data_suff <- sum(classified$data_sufficient, na.rm = TRUE)

  tibble::tibble(
    k                       = k,
    n_ros_sites_viable      = sum(site_counts$ROS      >= MIN_TYPE_SAMPLE),
    n_snowmelt_sites_viable = sum(site_counts$snowmelt >= MIN_TYPE_SAMPLE),
    n_ros_total             = sum(classified$flood_type_k == "ROS"),
    n_snowmelt_total        = sum(classified$flood_type_k == "snowmelt"),
    n_rainfall_total        = sum(classified$flood_type_k == "rainfall"),
    n_unclassified_total    = sum(classified$flood_type_k == "unclassified"),
    median_ros_per_site     = median(
      site_counts$ROS[site_counts$ROS > 0], na.rm = TRUE
    ),
    unclassified_rate       = if (n_data_suff > 0)
      sum(!classified$data_sufficient | is.na(classified$rain_frac)) /
        nrow(classified)
      else NA_real_
  )
})

# ============================================================================
# CONSOLE SUMMARY TABLE
# ============================================================================
message("\n=== B2 K-MULTIPLIER SWEEP RESULTS ===\n")
print(as.data.frame(sweep_results), row.names = FALSE)

# Flag the smallest k meeting both targets
meets_both <- sweep_results %>%
  dplyr::filter(
    n_ros_sites_viable      >= ROS_TARGET,
    n_snowmelt_sites_viable >= SNOWMELT_TARGET
  )

if (nrow(meets_both) > 0) {
  k_min <- min(meets_both$k)
  message(sprintf(
    "\n>>> Smallest k meeting both targets (>= %d ROS sites, >= %d Snowmelt sites): k = %.2f",
    ROS_TARGET, SNOWMELT_TARGET, k_min
  ))
  message(sprintf("    ROS viable sites: %d | Snowmelt viable sites: %d",
                  meets_both$n_ros_sites_viable[meets_both$k == k_min],
                  meets_both$n_snowmelt_sites_viable[meets_both$k == k_min]))
} else {
  message(sprintf(
    "\n>>> No k in the grid meets both targets (>= %d ROS, >= %d Snowmelt).",
    ROS_TARGET, SNOWMELT_TARGET
  ))
  message("    Consider extending K_GRID or revising targets.")
}

# ============================================================================
# FIG 01 — VIABILITY VS. K
# ============================================================================
message("\nBuilding Fig 01 — viability vs. k...")

fig01a <- ggplot(sweep_results, aes(x = k, y = n_ros_sites_viable)) +
  geom_line(colour = type_colors["ROS"], linewidth = 1.0) +
  geom_point(colour = type_colors["ROS"], size = 2.5) +
  geom_hline(yintercept = ROS_TARGET, linetype = "dashed",
             colour = "firebrick", linewidth = 0.7) +
  geom_vline(xintercept = 1.0, linetype = "dashed",
             colour = "grey40", linewidth = 0.7) +
  annotate("text", x = 1.03, y = -Inf, vjust = -0.5, hjust = 0,
           label = "k = 1 (current)", colour = "grey40", size = 3.2) +
  annotate("text", x = max(K_GRID), y = ROS_TARGET + 0.5, hjust = 1, vjust = 0,
           label = sprintf("Target: %d sites", ROS_TARGET),
           colour = "firebrick", size = 3.2) +
  scale_x_continuous(breaks = K_GRID) +
  scale_y_continuous(name = "Viable ROS sites") +
  labs(x = "k (B2_INTERCEPT multiplier)",
       title = "Fig 01a — Viable ROS sites vs. k",
       subtitle = sprintf("Viable = >= %d ROS peaks per site (pre-demotion)",
                          MIN_TYPE_SAMPLE)) +
  theme_diag()

fig01b <- ggplot(sweep_results, aes(x = k, y = n_snowmelt_sites_viable)) +
  geom_line(colour = type_colors["snowmelt"], linewidth = 1.0) +
  geom_point(colour = type_colors["snowmelt"], size = 2.5) +
  geom_hline(yintercept = SNOWMELT_TARGET, linetype = "dashed",
             colour = "darkgreen", linewidth = 0.7) +
  geom_vline(xintercept = 1.0, linetype = "dashed",
             colour = "grey40", linewidth = 0.7) +
  annotate("text", x = 1.03, y = -Inf, vjust = -0.5, hjust = 0,
           label = "k = 1 (current)", colour = "grey40", size = 3.2) +
  annotate("text", x = max(K_GRID), y = SNOWMELT_TARGET + 0.5, hjust = 1, vjust = 0,
           label = sprintf("Target: %d sites", SNOWMELT_TARGET),
           colour = "darkgreen", size = 3.2) +
  scale_x_continuous(breaks = K_GRID) +
  scale_y_continuous(name = "Viable Snowmelt sites") +
  labs(x = "k (B2_INTERCEPT multiplier)",
       title = "Fig 01b — Viable Snowmelt sites vs. k",
       subtitle = sprintf("Viable = >= %d Snowmelt peaks per site (pre-demotion)",
                          MIN_TYPE_SAMPLE)) +
  theme_diag()

fig01 <- fig01a + fig01b
save_fig(fig01, "fig_b2_k_viability.png")

# ============================================================================
# FIG 02 — TYPE COMPOSITION VS. K
# ============================================================================
message("Building Fig 02 — type composition vs. k...")

comp_long <- sweep_results %>%
  dplyr::select(k, ROS = n_ros_total, Snowmelt = n_snowmelt_total,
                Rainfall = n_rainfall_total, Unclassified = n_unclassified_total) %>%
  tidyr::pivot_longer(-k, names_to = "flood_type", values_to = "n_peaks") %>%
  dplyr::mutate(
    flood_type = factor(flood_type,
                        levels = c("Rainfall", "ROS", "Snowmelt", "Unclassified"))
  )

comp_colors <- c(
  "ROS"           = "#d7191c",
  "Snowmelt"      = "#4dac26",
  "Rainfall"      = "#2166ac",
  "Unclassified"  = "#d9d9d9"
)

fig02 <- ggplot(comp_long, aes(x = k, y = n_peaks, colour = flood_type,
                                group = flood_type)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.2) +
  geom_vline(xintercept = 1.0, linetype = "dashed",
             colour = "grey40", linewidth = 0.7) +
  scale_x_continuous(breaks = K_GRID) +
  scale_y_continuous(labels = scales::comma, name = "Total peaks") +
  scale_colour_manual(values = comp_colors, name = "Flood type") +
  labs(x = "k (B2_INTERCEPT multiplier)",
       title = "Fig 02 — Type composition vs. k",
       subtitle = "Total peaks per type across all sites; dashed = current k = 1") +
  theme_diag()

save_fig(fig02, "fig_b2_k_composition.png")

# ============================================================================
# FIG 03 — CLASSIFICATION SPACE AT SELECTED K VALUES
# ============================================================================
message("Building Fig 03 — classification space at selected k values...")

K_PLOT <- c(1.0, 1.5, 2.0, 3.0)

# Subsample for rendering speed
set.seed(42)
N_SUBSAMPLE <- 5000L
peaks_plot  <- if (nrow(peaks_slim) > N_SUBSAMPLE) {
  peaks_slim %>%
    dplyr::filter(data_sufficient, !is.na(rain_frac)) %>%
    dplyr::slice_sample(n = N_SUBSAMPLE)
} else {
  peaks_slim %>% dplyr::filter(data_sufficient, !is.na(rain_frac))
}

# Build classified data for each selected k
space_data <- purrr::map_dfr(K_PLOT, function(k) {
  apply_decision_tree(peaks_plot, k) %>%
    dplyr::mutate(
      k_label = sprintf("k = %.1f", k),
      k_val   = k
    )
}) %>%
  dplyr::mutate(
    k_label    = factor(k_label, levels = sprintf("k = %.1f", K_PLOT)),
    flood_type_k = factor(flood_type_k,
                          levels = c("rainfall", "ROS", "snowmelt", "unclassified"))
  )

# Build boundary lines over the observed rain_frac range
rain_seq <- seq(0, 1, length.out = 300)
b1_line  <- tibble::tibble(
  rain_frac = rain_seq,
  log10_swe = B1_INTERCEPT + B1_SLOPE * rain_frac,
  boundary  = "B1 (Snowmelt)"
)

boundary_lines <- purrr::map_dfr(K_PLOT, function(k) {
  tibble::tibble(
    rain_frac = rain_seq,
    log10_swe = (k * B2_INTERCEPT) + B2_SLOPE * rain_seq,
    boundary  = sprintf("B2 (k = %.1f)", k),
    k_label   = sprintf("k = %.1f", k),
    k_val     = k
  )
}) %>%
  dplyr::mutate(k_label = factor(k_label, levels = sprintf("k = %.1f", K_PLOT)))

fig03 <- ggplot(space_data,
                aes(x = rain_frac, y = log10_swe, colour = flood_type_k)) +
  geom_point(size = 0.5, alpha = 0.3) +
  # B1 boundary (same in every facet)
  geom_line(data = b1_line %>%
              dplyr::cross_join(
                tibble::tibble(k_label = factor(sprintf("k = %.1f", K_PLOT),
                                                levels = sprintf("k = %.1f", K_PLOT)))
              ),
            aes(x = rain_frac, y = log10_swe),
            colour = type_colors["snowmelt"], linewidth = 0.9,
            linetype = "solid", inherit.aes = FALSE) +
  # B2 boundary (varies by facet)
  geom_line(data = boundary_lines,
            aes(x = rain_frac, y = log10_swe),
            colour = type_colors["ROS"], linewidth = 0.9,
            linetype = "solid", inherit.aes = FALSE) +
  scale_colour_manual(values = type_colors, name = "Flood type",
                      labels = c("rainfall" = "Rainfall", "ROS" = "ROS",
                                 "snowmelt" = "Snowmelt",
                                 "unclassified" = "Unclassified")) +
  scale_x_continuous(name = "rain_frac",
                     limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(name = expression(log[10]*"(SWE + 1)")) +
  facet_wrap(~ k_label, ncol = 2) +
  labs(
    title    = "Fig 03 — Classification space at selected k values",
    subtitle = paste0("Green line = B1 (Snowmelt boundary, fixed);",
                      " Red line = effective B2 (varies with k).",
                      "\nSubsample of ", scales::comma(nrow(peaks_plot)),
                      " data-sufficient peaks shown.")
  ) +
  theme_diag() +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

save_fig(fig03, "fig_b2_k_classification_space.png")

message("\n=== B2 K-MULTIPLIER ANALYSIS COMPLETE ===")
message(sprintf("All figures saved to: %s/", OUT_DIR))

# ============================================================================
# FINE-GRID SWEEP (K_GRID_FINE)
# ============================================================================
message("\nRunning fine-grid B2 k-multiplier sweep across k = ",
        paste(sprintf("%.2f", K_GRID_FINE), collapse = ", "), "...")

sweep_results_fine <- purrr::map_dfr(K_GRID_FINE, function(k) {
  message(sprintf("  k = %.2f ...", k))

  classified <- apply_decision_tree(peaks_slim, k)

  site_counts <- classified %>%
    dplyr::filter(flood_type_k != "unclassified") %>%
    dplyr::count(site_no, flood_type_k) %>%
    tidyr::pivot_wider(names_from = flood_type_k, values_from = n,
                       values_fill = 0L) %>%
    dplyr::mutate(
      ROS      = dplyr::coalesce(ROS,      0L),
      snowmelt = dplyr::coalesce(snowmelt, 0L),
      rainfall = dplyr::coalesce(rainfall, 0L)
    )

  n_data_suff <- sum(classified$data_sufficient, na.rm = TRUE)

  tibble::tibble(
    k                       = k,
    n_ros_sites_viable      = sum(site_counts$ROS      >= MIN_TYPE_SAMPLE),
    n_snowmelt_sites_viable = sum(site_counts$snowmelt >= MIN_TYPE_SAMPLE),
    n_ros_total             = sum(classified$flood_type_k == "ROS"),
    n_snowmelt_total        = sum(classified$flood_type_k == "snowmelt"),
    n_rainfall_total        = sum(classified$flood_type_k == "rainfall"),
    n_unclassified_total    = sum(classified$flood_type_k == "unclassified"),
    median_ros_per_site     = median(
      site_counts$ROS[site_counts$ROS > 0], na.rm = TRUE
    ),
    unclassified_rate       = if (n_data_suff > 0)
      sum(!classified$data_sufficient | is.na(classified$rain_frac)) /
        nrow(classified)
      else NA_real_
  )
})

# ============================================================================
# FINE GRID CONSOLE SUMMARY
# ============================================================================
message("\n=== B2 K-MULTIPLIER FINE GRID RESULTS ===\n")

fine_display <- sweep_results_fine %>%
  dplyr::select(k, n_ros_sites_viable, n_snowmelt_sites_viable,
                n_ros_total, n_rainfall_total, median_ros_per_site) %>%
  dplyr::mutate(
    flag = dplyr::if_else(n_ros_sites_viable >= ROS_TARGET, "<<< TARGET MET", "")
  )
print(as.data.frame(fine_display), row.names = FALSE)

meets_fine <- sweep_results_fine %>%
  dplyr::filter(n_ros_sites_viable >= ROS_TARGET)

if (nrow(meets_fine) > 0) {
  k_fine_min <- min(meets_fine$k)
  ros_at_k   <- meets_fine$n_ros_sites_viable[meets_fine$k == k_fine_min]
  snow_at_k  <- meets_fine$n_snowmelt_sites_viable[meets_fine$k == k_fine_min]
  message(sprintf(
    "\n>>> Fine grid: smallest k meeting >= %d ROS sites: k = %.2f",
    ROS_TARGET, k_fine_min
  ))
  message(sprintf("    ROS viable: %d | Snowmelt viable: %d",
                  ros_at_k, snow_at_k))
} else {
  k_fine_min <- NA_real_
  message(sprintf(
    "\n>>> Fine grid: no k in [%.2f, %.2f] meets >= %d ROS sites.",
    min(K_GRID_FINE), max(K_GRID_FINE), ROS_TARGET
  ))
}

# ============================================================================
# FIG 04 — FINE-GRID VIABILITY VS. K
# ============================================================================
message("\nBuilding Fig 04 — fine-grid viability vs. k...")

fine_breaks <- round(K_GRID_FINE, 2)

fig04a <- ggplot(sweep_results_fine, aes(x = k, y = n_ros_sites_viable)) +
  geom_line(colour = type_colors["ROS"], linewidth = 1.0) +
  geom_point(colour = type_colors["ROS"], size = 2.5) +
  geom_hline(yintercept = ROS_TARGET, linetype = "dashed",
             colour = "firebrick", linewidth = 0.7) +
  annotate("text", x = max(K_GRID_FINE), y = ROS_TARGET + 0.3,
           hjust = 1, vjust = 0,
           label = sprintf("Target: %d sites", ROS_TARGET),
           colour = "firebrick", size = 3.2)

if (!is.na(k_fine_min)) {
  ros_cross <- sweep_results_fine$n_ros_sites_viable[sweep_results_fine$k == k_fine_min]
  fig04a <- fig04a +
    geom_vline(xintercept = k_fine_min, linetype = "dashed",
               colour = "grey40", linewidth = 0.7) +
    annotate("point", x = k_fine_min, y = ros_cross,
             colour = "firebrick", fill = "white",
             size = 4.5, shape = 21, stroke = 1.5) +
    annotate("text",
             x = k_fine_min + 0.005, y = ros_cross + 0.3,
             hjust = 0, vjust = 0,
             label = sprintf("k = %.2f\n%d sites", k_fine_min, ros_cross),
             colour = "firebrick", size = 3.2)
}

fig04a <- fig04a +
  scale_x_continuous(breaks = fine_breaks) +
  scale_y_continuous(name = "Viable ROS sites") +
  labs(x = "k (B2_INTERCEPT multiplier)",
       title = "Fig 04a — Viable ROS sites vs. k (fine grid)",
       subtitle = sprintf("Viable = >= %d ROS peaks per site; k in [%.2f, %.2f] by 0.05",
                          MIN_TYPE_SAMPLE, min(K_GRID_FINE), max(K_GRID_FINE))) +
  theme_diag() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig04b <- ggplot(sweep_results_fine, aes(x = k, y = n_snowmelt_sites_viable)) +
  geom_line(colour = type_colors["snowmelt"], linewidth = 1.0) +
  geom_point(colour = type_colors["snowmelt"], size = 2.5) +
  geom_hline(yintercept = SNOWMELT_TARGET, linetype = "dashed",
             colour = "darkgreen", linewidth = 0.7) +
  annotate("text", x = max(K_GRID_FINE), y = SNOWMELT_TARGET + 0.3,
           hjust = 1, vjust = 0,
           label = sprintf("Target: %d sites", SNOWMELT_TARGET),
           colour = "darkgreen", size = 3.2) +
  scale_x_continuous(breaks = fine_breaks) +
  scale_y_continuous(name = "Viable Snowmelt sites") +
  labs(x = "k (B2_INTERCEPT multiplier)",
       title = "Fig 04b — Viable Snowmelt sites vs. k (fine grid)",
       subtitle = sprintf("Viable = >= %d Snowmelt peaks per site",
                          MIN_TYPE_SAMPLE)) +
  theme_diag() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig04 <- fig04a + fig04b
save_fig(fig04, "fig_b2_k_viability_fine.png")

message("\n=== FINE GRID ANALYSIS COMPLETE ===")
message(sprintf("All figures saved to: %s/", OUT_DIR))
