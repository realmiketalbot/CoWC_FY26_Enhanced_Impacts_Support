#!/usr/bin/env Rscript

# ============================================================================
# b2_threshold_comparison.R
#
# Evaluates whether a simple fixed rain_frac threshold for the ROS/Rainfall
# boundary performs comparably to the SWE-dependent GMM B2 boundary:
#
#   Fixed:  rain_frac < thresh          (vertical line in feature space)
#   GMM B2: rain_frac < B2_INTERCEPT + B2_SLOPE * log10_swe  (diagonal)
#
# The B1 Snowmelt boundary is NOT varied. Only the ROS/Rainfall split changes.
# GMM B2 at k=1.0 (current) and k=2.0 are included as reference cases.
#
# Goal: find the smallest fixed threshold meeting >= 30 viable ROS sites while
# maintaining >= 50 viable Snowmelt sites; assess whether the added complexity
# of the SWE-dependent boundary is justified.
#
# INPUT:  _data/pot_peaks_classified_adaptive.rds
# OUTPUT: fig_b2_thresh_viability.png, fig_b2_thresh_composition.png,
#         fig_b2_thresh_classification_space.png,
#         fig_b2_thresh_vs_gmm_agreement.png
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

THRESH_GRID     <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50)

# GMM reference k values included alongside the fixed-threshold sweep
GMM_K_REF       <- c(1.0, 2.0)

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
  "ROS"          = "#d7191c",
  "snowmelt"     = "#4dac26",
  "rainfall"     = "#2166ac",
  "unclassified" = "#d9d9d9"
)

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading classified peaks...")
peaks <- readRDS(PEAKS_FILE)
message(sprintf("  %d peaks across %d sites", nrow(peaks), dplyr::n_distinct(peaks$site_no)))

if ("log10_swe" %in% names(peaks)) {
  message("  log10_swe column found — using directly.")
} else {
  message("  log10_swe not found — computing from swe_mean_mm.")
  peaks <- peaks %>%
    dplyr::mutate(
      log10_swe = log10(dplyr::coalesce(swe_mean_mm, 0) + SWE_LOG_OFFSET)
    )
}

required_cols <- c("site_no", "rain_frac", "log10_swe", "data_sufficient")
missing_cols  <- setdiff(required_cols, names(peaks))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

peaks_slim <- peaks %>%
  dplyr::select(site_no, rain_frac, log10_swe, data_sufficient)

n_data_suff <- sum(peaks_slim$data_sufficient, na.rm = TRUE)
message(sprintf("  Data-sufficient peaks: %d", n_data_suff))

# ============================================================================
# DECISION TREE HELPERS
# ============================================================================

# Fixed rain_frac threshold: ROS if rain_frac < thresh (no SWE dependence)
apply_fixed_thresh <- function(peaks_df, thresh) {
  peaks_df %>%
    dplyr::mutate(
      flood_type_k = dplyr::case_when(
        !data_sufficient | is.na(rain_frac)                    ~ "unclassified",
        log10_swe > B1_INTERCEPT + B1_SLOPE * rain_frac        ~ "snowmelt",
        rain_frac < thresh                                     ~ "ROS",
        TRUE                                                   ~ "rainfall"
      )
    )
}

# GMM SWE-dependent B2 with k multiplier (matches b2_k_multiplier_analysis.R)
apply_gmm_k <- function(peaks_df, k) {
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

# Compute summary metrics from a classified peaks dataframe
compute_metrics <- function(classified_df) {
  site_counts <- classified_df %>%
    dplyr::filter(flood_type_k != "unclassified") %>%
    dplyr::count(site_no, flood_type_k) %>%
    tidyr::pivot_wider(names_from = flood_type_k, values_from = n,
                       values_fill = 0L) %>%
    dplyr::mutate(
      ROS      = dplyr::coalesce(ROS,      0L),
      snowmelt = dplyr::coalesce(snowmelt, 0L),
      rainfall = dplyr::coalesce(rainfall, 0L)
    )

  tibble::tibble(
    n_ros_sites_viable      = sum(site_counts$ROS      >= MIN_TYPE_SAMPLE),
    n_snowmelt_sites_viable = sum(site_counts$snowmelt >= MIN_TYPE_SAMPLE),
    n_ros_total             = sum(classified_df$flood_type_k == "ROS"),
    n_snowmelt_total        = sum(classified_df$flood_type_k == "snowmelt"),
    n_rainfall_total        = sum(classified_df$flood_type_k == "rainfall"),
    n_unclassified_total    = sum(classified_df$flood_type_k == "unclassified"),
    median_ros_per_site     = median(
      site_counts$ROS[site_counts$ROS > 0], na.rm = TRUE
    ),
    unclassified_rate       = sum(classified_df$flood_type_k == "unclassified") /
                              nrow(classified_df)
  )
}

# ============================================================================
# FIXED-THRESHOLD SWEEP
# ============================================================================
message("\nRunning fixed-threshold sweep across thresh = ",
        paste(THRESH_GRID, collapse = ", "), "...")

thresh_results <- purrr::map_dfr(THRESH_GRID, function(thresh) {
  message(sprintf("  thresh = %.2f ...", thresh))
  classified <- apply_fixed_thresh(peaks_slim, thresh)
  dplyr::bind_cols(
    tibble::tibble(boundary = "fixed", thresh = thresh, label = sprintf("%.2f", thresh)),
    compute_metrics(classified)
  )
})

# ============================================================================
# GMM REFERENCE CASES
# ============================================================================
message("\nComputing GMM reference cases (k = ",
        paste(GMM_K_REF, collapse = ", "), ")...")

gmm_results <- purrr::map_dfr(GMM_K_REF, function(k) {
  message(sprintf("  GMM k = %.1f ...", k))
  classified <- apply_gmm_k(peaks_slim, k)
  # Equivalent fixed rain_frac threshold at log10_swe = 0 (SWE = 0 mm):
  # rain_frac < k*B2_INTERCEPT + B2_SLOPE*0 = k*B2_INTERCEPT
  equiv_thresh <- k * B2_INTERCEPT
  dplyr::bind_cols(
    tibble::tibble(
      boundary = "gmm",
      thresh   = equiv_thresh,   # x-axis position: equivalent rain_frac at SWE=0
      label    = sprintf("GMM k=%.1f", k)
    ),
    compute_metrics(classified)
  )
})

all_results <- dplyr::bind_rows(thresh_results, gmm_results) %>%
  dplyr::arrange(boundary, thresh)

# ============================================================================
# CONSOLE SUMMARY
# ============================================================================
message("\n=== FIXED THRESHOLD COMPARISON RESULTS ===\n")
print(as.data.frame(all_results %>% dplyr::select(-boundary)), row.names = FALSE)

# Flag smallest fixed threshold meeting both targets
fixed_meeting <- thresh_results %>%
  dplyr::filter(
    n_ros_sites_viable      >= ROS_TARGET,
    n_snowmelt_sites_viable >= SNOWMELT_TARGET
  )

if (nrow(fixed_meeting) > 0) {
  best_fixed <- fixed_meeting %>% dplyr::slice_min(thresh, n = 1)
  message(sprintf(
    "\n>>> Smallest fixed threshold meeting both targets: rain_frac < %.2f",
    best_fixed$thresh
  ))
  message(sprintf("    ROS viable sites: %d | Snowmelt viable sites: %d",
                  best_fixed$n_ros_sites_viable,
                  best_fixed$n_snowmelt_sites_viable))

  # Compare with GMM k=2.0
  gmm_k2 <- gmm_results %>% dplyr::filter(abs(thresh - 2.0 * B2_INTERCEPT) < 1e-6)
  if (nrow(gmm_k2) > 0) {
    message(sprintf(
      "\n    GMM k=2.0 (equiv. thresh ≈ %.3f): ROS viable = %d, Snowmelt viable = %d",
      2.0 * B2_INTERCEPT, gmm_k2$n_ros_sites_viable, gmm_k2$n_snowmelt_sites_viable
    ))
    rain_diff <- best_fixed$n_rainfall_total - gmm_k2$n_rainfall_total
    message(sprintf(
      "    Rainfall difference (fixed vs GMM k=2): %+d peaks",
      rain_diff
    ))
  }
} else {
  message(sprintf(
    "\n>>> No fixed threshold in THRESH_GRID meets both targets (>= %d ROS, >= %d Snowmelt).",
    ROS_TARGET, SNOWMELT_TARGET
  ))
}

# ============================================================================
# FIG 01 — VIABILITY VS. THRESHOLD
# ============================================================================
message("\nBuilding Fig 01 — viability vs. threshold...")

# GMM reference intercepts for vertical lines (equivalent rain_frac at SWE=0)
gmm_xrefs <- gmm_results %>%
  dplyr::select(thresh, label)

fig01a <- ggplot(thresh_results, aes(x = thresh, y = n_ros_sites_viable)) +
  geom_line(colour = type_colors["ROS"], linewidth = 1.0) +
  geom_point(colour = type_colors["ROS"], size = 2.5) +
  geom_point(data = gmm_results,
             aes(x = thresh, y = n_ros_sites_viable, shape = label),
             colour = "grey20", size = 3.5, fill = "white", stroke = 1.2) +
  geom_vline(data = gmm_xrefs,
             aes(xintercept = thresh, linetype = label),
             colour = "grey40", linewidth = 0.7) +
  geom_hline(yintercept = ROS_TARGET, linetype = "dashed",
             colour = "firebrick", linewidth = 0.7) +
  annotate("text", x = max(THRESH_GRID), y = ROS_TARGET + 0.5,
           hjust = 1, vjust = 0,
           label = sprintf("Target: %d sites", ROS_TARGET),
           colour = "firebrick", size = 3.2) +
  scale_x_continuous(breaks = THRESH_GRID, name = "Fixed rain_frac threshold") +
  scale_y_continuous(name = "Viable ROS sites") +
  scale_shape_manual(values = c("GMM k=1.0" = 21, "GMM k=2.0" = 24),
                     name = "GMM reference") +
  scale_linetype_manual(values = c("GMM k=1.0" = "dotted", "GMM k=2.0" = "dashed"),
                        name = "GMM ref. (SWE=0 equiv.)") +
  labs(title = "Fig 01a — Viable ROS sites vs. fixed threshold",
       subtitle = sprintf("Viable = >= %d ROS peaks per site; GMM refs at rain_frac = k × B2_INTERCEPT",
                          MIN_TYPE_SAMPLE)) +
  theme_diag()

fig01b <- ggplot(thresh_results, aes(x = thresh, y = n_snowmelt_sites_viable)) +
  geom_line(colour = type_colors["snowmelt"], linewidth = 1.0) +
  geom_point(colour = type_colors["snowmelt"], size = 2.5) +
  geom_point(data = gmm_results,
             aes(x = thresh, y = n_snowmelt_sites_viable, shape = label),
             colour = "grey20", size = 3.5, fill = "white", stroke = 1.2) +
  geom_vline(data = gmm_xrefs,
             aes(xintercept = thresh, linetype = label),
             colour = "grey40", linewidth = 0.7) +
  geom_hline(yintercept = SNOWMELT_TARGET, linetype = "dashed",
             colour = "darkgreen", linewidth = 0.7) +
  annotate("text", x = max(THRESH_GRID), y = SNOWMELT_TARGET + 0.5,
           hjust = 1, vjust = 0,
           label = sprintf("Target: %d sites", SNOWMELT_TARGET),
           colour = "darkgreen", size = 3.2) +
  scale_x_continuous(breaks = THRESH_GRID, name = "Fixed rain_frac threshold") +
  scale_y_continuous(name = "Viable Snowmelt sites") +
  scale_shape_manual(values = c("GMM k=1.0" = 21, "GMM k=2.0" = 24),
                     name = "GMM reference") +
  scale_linetype_manual(values = c("GMM k=1.0" = "dotted", "GMM k=2.0" = "dashed"),
                        name = "GMM ref. (SWE=0 equiv.)") +
  labs(title = "Fig 01b — Viable Snowmelt sites vs. fixed threshold") +
  theme_diag()

fig01 <- fig01a + fig01b + patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
save_fig(fig01, "fig_b2_thresh_viability.png")

# ============================================================================
# FIG 02 — TYPE COMPOSITION VS. THRESHOLD
# ============================================================================
message("Building Fig 02 — type composition vs. threshold...")

comp_fixed <- thresh_results %>%
  dplyr::select(thresh, label,
                ROS = n_ros_total, Snowmelt = n_snowmelt_total,
                Rainfall = n_rainfall_total, Unclassified = n_unclassified_total) %>%
  tidyr::pivot_longer(c(ROS, Snowmelt, Rainfall, Unclassified),
                      names_to = "flood_type", values_to = "n_peaks") %>%
  dplyr::mutate(
    flood_type = factor(flood_type,
                        levels = c("Rainfall", "ROS", "Snowmelt", "Unclassified"))
  )

comp_gmm <- gmm_results %>%
  dplyr::select(thresh, label,
                ROS = n_ros_total, Snowmelt = n_snowmelt_total,
                Rainfall = n_rainfall_total, Unclassified = n_unclassified_total) %>%
  tidyr::pivot_longer(c(ROS, Snowmelt, Rainfall, Unclassified),
                      names_to = "flood_type", values_to = "n_peaks") %>%
  dplyr::mutate(
    flood_type = factor(flood_type,
                        levels = c("Rainfall", "ROS", "Snowmelt", "Unclassified"))
  )

comp_colors <- c(
  "ROS"          = "#d7191c",
  "Snowmelt"     = "#4dac26",
  "Rainfall"     = "#2166ac",
  "Unclassified" = "#d9d9d9"
)

fig02 <- ggplot(comp_fixed, aes(x = thresh, y = n_peaks,
                                colour = flood_type, group = flood_type)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.2) +
  geom_point(data = comp_gmm,
             aes(x = thresh, y = n_peaks, shape = label),
             colour = "grey20", size = 3.5, fill = "white", stroke = 1.2) +
  geom_vline(data = gmm_xrefs,
             aes(xintercept = thresh, linetype = label),
             colour = "grey40", linewidth = 0.7) +
  scale_x_continuous(breaks = THRESH_GRID, name = "Fixed rain_frac threshold") +
  scale_y_continuous(labels = scales::comma, name = "Total peaks") +
  scale_colour_manual(values = comp_colors, name = "Flood type") +
  scale_shape_manual(values = c("GMM k=1.0" = 21, "GMM k=2.0" = 24),
                     name = "GMM reference") +
  scale_linetype_manual(values = c("GMM k=1.0" = "dotted", "GMM k=2.0" = "dashed"),
                        name = "GMM ref. (SWE=0 equiv.)") +
  labs(title = "Fig 02 — Type composition vs. fixed threshold",
       subtitle = "Total peaks per type; GMM reference points at SWE=0 equivalent threshold") +
  theme_diag() +
  theme(legend.box = "horizontal")
save_fig(fig02, "fig_b2_thresh_composition.png")

# ============================================================================
# FIG 03 — CLASSIFICATION SPACE AT SELECTED THRESHOLDS
# ============================================================================
message("Building Fig 03 — classification space at selected thresholds...")

THRESH_PLOT <- c(0.10, 0.20, 0.30, 0.40)

set.seed(42)
N_SUBSAMPLE  <- 5000L
peaks_ss     <- peaks_slim %>%
  dplyr::filter(data_sufficient, !is.na(rain_frac))
peaks_plot   <- if (nrow(peaks_ss) > N_SUBSAMPLE) {
  dplyr::slice_sample(peaks_ss, n = N_SUBSAMPLE)
} else {
  peaks_ss
}

# Fixed-threshold panels
space_fixed <- purrr::map_dfr(THRESH_PLOT, function(thresh) {
  apply_fixed_thresh(peaks_plot, thresh) %>%
    dplyr::mutate(
      panel_label = sprintf("Fixed: rain_frac < %.2f", thresh),
      panel_type  = "fixed",
      thresh_val  = thresh
    )
})

# GMM k=1.0 panel
space_gmm <- apply_gmm_k(peaks_plot, k = 1.0) %>%
  dplyr::mutate(
    panel_label = "GMM B2 (k = 1.0)",
    panel_type  = "gmm",
    thresh_val  = NA_real_
  )

space_data <- dplyr::bind_rows(space_fixed, space_gmm) %>%
  dplyr::mutate(
    panel_label  = factor(panel_label,
                          levels = c(sprintf("Fixed: rain_frac < %.2f", THRESH_PLOT),
                                     "GMM B2 (k = 1.0)")),
    flood_type_k = factor(flood_type_k,
                          levels = c("rainfall", "ROS", "snowmelt", "unclassified"))
  )

# B1 boundary line (same across all panels)
rain_seq <- seq(0, 1, length.out = 300)
b1_df    <- tibble::tibble(
  rain_frac = rain_seq,
  log10_swe = B1_INTERCEPT + B1_SLOPE * rain_seq
)

# ROS/Rainfall boundary per panel
all_panels <- c(sprintf("Fixed: rain_frac < %.2f", THRESH_PLOT), "GMM B2 (k = 1.0)")

boundary_df <- purrr::map_dfr(all_panels, function(lbl) {
  if (startsWith(lbl, "Fixed")) {
    thresh <- as.numeric(stringr::str_extract(lbl, "[0-9.]+$"))
    # Vertical line: represented as a short segment spanning the SWE range
    tibble::tibble(
      panel_label = lbl,
      x_start     = thresh,
      x_end       = thresh,
      y_start     = min(peaks_plot$log10_swe, na.rm = TRUE),
      y_end       = max(peaks_plot$log10_swe, na.rm = TRUE),
      is_diagonal = FALSE
    )
  } else {
    # GMM diagonal: rain_frac = B2_INTERCEPT + B2_SLOPE * log10_swe
    swe_seq <- seq(min(peaks_plot$log10_swe, na.rm = TRUE),
                   max(peaks_plot$log10_swe, na.rm = TRUE),
                   length.out = 200)
    tibble::tibble(
      panel_label = lbl,
      x_start     = B2_INTERCEPT + B2_SLOPE * swe_seq,
      x_end       = NA_real_,
      y_start     = swe_seq,
      y_end       = NA_real_,
      is_diagonal = TRUE
    )
  }
}) %>%
  dplyr::mutate(panel_label = factor(panel_label, levels = all_panels))

# Vertical boundary segments (fixed threshold panels)
boundary_vline <- boundary_df %>%
  dplyr::filter(!is_diagonal) %>%
  dplyr::distinct(panel_label, x_start, y_start, y_end)

# Diagonal boundary (GMM panel)
boundary_diag <- boundary_df %>%
  dplyr::filter(is_diagonal) %>%
  dplyr::select(panel_label, rain_frac = x_start, log10_swe = y_start)

# B1 replicated across panels
b1_all <- b1_df %>%
  dplyr::cross_join(tibble::tibble(panel_label = factor(all_panels, levels = all_panels)))

fig03 <- ggplot(space_data,
                aes(x = rain_frac, y = log10_swe, colour = flood_type_k)) +
  geom_point(size = 0.5, alpha = 0.3) +
  # B1 (Snowmelt) boundary — green, consistent across all panels
  geom_line(data = b1_all,
            aes(x = rain_frac, y = log10_swe),
            colour = type_colors["snowmelt"], linewidth = 0.9,
            inherit.aes = FALSE) +
  # ROS/Rainfall boundary — red vertical lines for fixed panels
  geom_segment(data = boundary_vline,
               aes(x = x_start, xend = x_start, y = y_start, yend = y_end),
               colour = type_colors["ROS"], linewidth = 0.9,
               inherit.aes = FALSE) +
  # ROS/Rainfall boundary — red diagonal for GMM panel
  geom_line(data = boundary_diag,
            aes(x = rain_frac, y = log10_swe),
            colour = type_colors["ROS"], linewidth = 0.9,
            inherit.aes = FALSE) +
  scale_colour_manual(values = type_colors, name = "Flood type",
                      labels = c("rainfall"     = "Rainfall",
                                 "ROS"          = "ROS",
                                 "snowmelt"     = "Snowmelt",
                                 "unclassified" = "Unclassified")) +
  scale_x_continuous(name = "rain_frac",
                     limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(name = expression(log[10]*"(SWE + 1)")) +
  facet_wrap(~ panel_label, ncol = 2) +
  labs(
    title    = "Fig 03 — Classification space: fixed thresholds vs. GMM B2",
    subtitle = paste0(
      "Green line = B1 (Snowmelt); Red line = ROS/Rainfall boundary.\n",
      "Fixed: vertical; GMM k=1.0: diagonal. ",
      scales::comma(nrow(peaks_plot)), " data-sufficient peaks."
    )
  ) +
  theme_diag() +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

save_fig(fig03, "fig_b2_thresh_classification_space.png", height = 9)

# ============================================================================
# FIG 04 — AGREEMENT: BEST FIXED THRESHOLD vs. GMM k=1.0
# ============================================================================
message("Building Fig 04 — fixed threshold vs. GMM k=1.0 agreement...")

# Determine best fixed threshold: smallest meeting ROS_TARGET
best_thresh_val <- if (nrow(fixed_meeting) > 0) {
  min(fixed_meeting$thresh)
} else {
  # Fallback: use largest thresh in grid (maximises ROS viability)
  message("  No fixed threshold meets ROS target — using largest THRESH_GRID value for FIG 04.")
  max(THRESH_GRID)
}
message(sprintf("  FIG 04 comparison: fixed thresh = %.2f vs. GMM k=1.0", best_thresh_val))

# Classify full dataset under both schemes (not subsampled)
classified_fixed <- apply_fixed_thresh(peaks_slim, best_thresh_val) %>%
  dplyr::rename(type_fixed = flood_type_k)

classified_gmm <- apply_gmm_k(peaks_slim, k = 1.0) %>%
  dplyr::rename(type_gmm = flood_type_k)

agreement_df <- dplyr::bind_cols(
  classified_fixed %>% dplyr::select(site_no, rain_frac, log10_swe,
                                      data_sufficient, type_fixed),
  classified_gmm   %>% dplyr::select(type_gmm)
) %>%
  dplyr::filter(data_sufficient, !is.na(rain_frac)) %>%
  dplyr::mutate(
    # Agreement status from the fixed-threshold perspective
    agree = type_fixed == type_gmm,
    agreement_label = dplyr::case_when(
      type_fixed == type_gmm                                 ~ "Agree",
      type_fixed == "ROS"      & type_gmm == "rainfall"     ~ "Fixed→ROS, GMM→Rainfall",
      type_fixed == "rainfall" & type_gmm == "ROS"          ~ "Fixed→Rainfall, GMM→ROS",
      type_fixed == "snowmelt" | type_gmm == "snowmelt"     ~ "Snowmelt mismatch",
      TRUE                                                   ~ "Other mismatch"
    ),
    type_gmm_label = factor(
      dplyr::recode(type_gmm,
                    rainfall = "Rainfall", ROS = "ROS",
                    snowmelt = "Snowmelt", unclassified = "Unclassified"),
      levels = c("Rainfall", "ROS", "Snowmelt", "Unclassified")
    )
  )

agree_colors <- c(
  "Agree"                     = "grey75",
  "Fixed→ROS, GMM→Rainfall"   = "#d7191c",
  "Fixed→Rainfall, GMM→ROS"   = "#fdae61",
  "Snowmelt mismatch"         = "#abdda4",
  "Other mismatch"            = "#e0e0e0"
)

fig04 <- ggplot(agreement_df,
                aes(x = type_gmm_label, fill = agreement_label)) +
  geom_bar(position = "stack", colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = agree_colors, name = "Agreement status") +
  scale_y_continuous(labels = scales::comma,
                     name   = "Number of data-sufficient peaks") +
  scale_x_discrete(name = "GMM k=1.0 flood type") +
  labs(
    title    = sprintf("Fig 04 — Agreement: fixed thresh = %.2f vs. GMM k = 1.0", best_thresh_val),
    subtitle = paste0(
      "Red = peaks reclassified Rainfall→ROS by the fixed threshold. ",
      "Orange = peaks reclassified ROS→Rainfall.\n",
      "Key question: are the extra ROS peaks (red) genuinely in the ROS zone of GMM space?"
    ),
    caption  = "Data-sufficient peaks only; no subsampling."
  ) +
  theme_diag()

save_fig(fig04, "fig_b2_thresh_vs_gmm_agreement.png", height = 6)

message("\n=== THRESHOLD COMPARISON COMPLETE ===")
message(sprintf("All figures saved to: %s/", OUT_DIR))
