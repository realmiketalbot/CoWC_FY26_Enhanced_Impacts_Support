#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: ROS vs Snowmelt Clustering
# ============================================================================
# For all candidate melt-season peaks (those with meaningful melt in the
# lookback window), computes two features from the daily rain/melt time series:
#
#   weighted_rain_frac:
#     sum(rain_i) / sum(rain_i + melt_i) across active days
#     Captures overall balance of rain vs melt flux, magnitude-weighted.
#     0 = pure melt, 1 = pure rain.
#
#   weighted_ros_frac:
#     sum(min(rain_i, melt_i)) / sum(rain_i + melt_i) across active days
#     Captures fraction of total flux occurring on days where both rain
#     and melt are simultaneously active. Days with only one process
#     contribute zero to the numerator.
#     0 = no co-occurrence, higher = more ROS-like.
#
# Active days: both rain and melt must individually exceed NOISE_FLOOR_MM
# for a day to contribute to weighted_ros_frac. For weighted_rain_frac,
# only the total flux (rain + melt) must exceed NOISE_FLOOR_MM.
#
# Candidate peaks: all peaks with total melt over the lookback window
# exceeding MELT_PRESENT_MM (i.e., melt is at least a minor contributor).
#
# Clustering: k-means in (weighted_rain_frac, weighted_ros_frac) space,
# k=3 (snowmelt, ROS, rainfall). Silhouette analysis used to evaluate k.
#
# INPUTS:
#   - _data/accum_curves_60d.rds    : long-format accumulation curves
#                                     with daily rain_mm_day, melt_mm_day
#   - _data/pot_peaks_all.rds       : for peak metadata
#   - _data/site_metadata.rds       : for elevation etc.
#
# OUTPUTS:
#   - _data/ros_snowmelt_clusters.rds   : peak-level cluster assignments
#   - _diagnostics/accum_curves/clustering_*.png
# ============================================================================

library(tidyverse)
library(patchwork)
library(cluster)      # silhouette()

# ============================================================================
# CONFIGURATION
# ============================================================================
accum_file    <- "_data/accum_curves_60d.rds"
peaks_file    <- "_data/pot_peaks_all.rds"
meta_file     <- "_data/site_metadata.rds"
cluster_output <- "_data/ros_snowmelt_clusters.rds"
out_dir       <- "_diagnostics/accum_curves"

NOISE_FLOOR_MM  <- 0.1   # mm/day below which a signal is treated as absent
MELT_PRESENT_MM <- 5.0   # total melt over window required to be a candidate peak
LOOKBACK_DAYS   <- 60    # must match accum_curves file

K_MAX           <- 6     # maximum k to evaluate in silhouette analysis
SET_SEED        <- 42

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading accumulation curves...")
accum <- readRDS(accum_file)

message("Loading site metadata...")
meta <- readRDS(meta_file) %>%
  select(site_no, elev_mean_m, drain_sqkm)

# ============================================================================
# COMPUTE PEAK-LEVEL FEATURES
# ============================================================================
message("Computing weighted_rain_frac and weighted_ros_frac...")

features <- accum %>%
  # Treat sub-noise values as zero
  mutate(
    rain_active = if_else(rain_mm_day >= NOISE_FLOOR_MM, rain_mm_day, 0),
    melt_active = if_else(melt_mm_day >= NOISE_FLOOR_MM, melt_mm_day, 0),
    total_flux  = rain_active + melt_active,
    # Co-occurrence flux: lesser of the two on days where both are active
    ros_flux    = if_else(
      rain_active >= NOISE_FLOOR_MM & melt_active >= NOISE_FLOOR_MM,
      pmin(rain_active, melt_active),
      0
    )
  ) %>%
  group_by(site_no, peak_id, peak_date, peak_flow) %>%
  summarise(
    total_rain_flux  = sum(rain_active, na.rm = TRUE),
    total_melt_flux  = sum(melt_active, na.rm = TRUE),
    total_flux       = sum(total_flux,  na.rm = TRUE),
    total_ros_flux   = sum(ros_flux,    na.rm = TRUE),
    n_active_days    = sum(total_flux > 0),
    n_days_data      = first(n_days_data),
    .groups          = "drop"
  ) %>%
  mutate(
    weighted_rain_frac = if_else(
      total_flux > 0,
      total_rain_flux / total_flux,
      NA_real_
    ),
    weighted_ros_frac  = if_else(
      total_flux > 0,
      total_ros_flux / total_flux,
      NA_real_
    )
  )

message(sprintf("Features computed for %d peaks", nrow(features)))

# ============================================================================
# FILTER TO CANDIDATE MELT-SEASON PEAKS
# ============================================================================
candidates <- features %>%
  filter(
    total_melt_flux >= MELT_PRESENT_MM,
    !is.na(weighted_rain_frac),
    !is.na(weighted_ros_frac),
    n_days_data >= 57   # near-complete 60-day window
  )

message(sprintf("Candidate melt-season peaks: %d of %d total",
                nrow(candidates), nrow(features)))

# ============================================================================
# SILHOUETTE ANALYSIS: CHOOSE K
# ============================================================================
message("Running silhouette analysis for k = 2 to ", K_MAX, "...")

set.seed(SET_SEED)
X <- candidates %>%
  select(weighted_rain_frac, weighted_ros_frac) %>%
  as.matrix()

sil_scores <- map_dbl(2:K_MAX, function(k) {
  km  <- kmeans(X, centers = k, nstart = 25, iter.max = 100)
  sil <- silhouette(km$cluster, dist(X))
  mean(sil[, "sil_width"])
})

sil_df <- tibble(
  k         = 2:K_MAX,
  mean_sil  = sil_scores
)

message("\nMean silhouette width by k:")
print(sil_df)

best_k <- sil_df %>% slice_max(mean_sil, n = 1) %>% pull(k)
message(sprintf("Best k by silhouette: %d", best_k))

# ============================================================================
# FIT FINAL CLUSTERING (k=3 forced, plus best_k if different)
# ============================================================================
set.seed(SET_SEED)
km3 <- kmeans(X, centers = 3, nstart = 25, iter.max = 100)

candidates <- candidates %>%
  mutate(cluster_k3 = km3$cluster)

# Label clusters by their centroid position in feature space
# Snowmelt: low rain_frac, low ros_frac
# ROS:      moderate rain_frac, high ros_frac
# Rainfall: high rain_frac, low ros_frac
centroids <- as_tibble(km3$centers) %>%
  mutate(cluster_k3 = row_number()) %>%
  mutate(
    flood_type = case_when(
      weighted_rain_frac < 0.33                          ~ "snowmelt",
      weighted_ros_frac  > 0.15                          ~ "rain_on_snow",
      TRUE                                               ~ "rainfall"
    )
  )

message("\nCluster centroids (k=3):")
print(centroids)

candidates <- candidates %>%
  left_join(centroids %>% select(cluster_k3, flood_type),
            by = "cluster_k3")

# If best_k != 3, also fit and attach
if (best_k != 3) {
  set.seed(SET_SEED)
  km_best <- kmeans(X, centers = best_k, nstart = 25, iter.max = 100)
  candidates <- candidates %>%
    mutate(cluster_best_k = km_best$cluster)
}

# ============================================================================
# SAVE
# ============================================================================
saveRDS(candidates, cluster_output)
message(sprintf("\nSaved cluster assignments to %s", cluster_output))

# ============================================================================
# FIGURES
# ============================================================================

# --- Fig 1: Silhouette scores by k ---
p1 <- ggplot(sil_df, aes(x = k, y = mean_sil)) +
  geom_line(linewidth = 0.9, color = "grey30") +
  geom_point(size = 2.5, color = "grey30") +
  geom_point(data = sil_df %>% filter(k == best_k),
             size = 4, shape = 8, color = "#d73027") +
  geom_point(data = sil_df %>% filter(k == 3),
             size = 4, shape = 1, color = "#2166ac") +
  annotate("text", x = best_k + 0.15, y = sil_df$mean_sil[best_k - 1],
           label = sprintf("best k=%d", best_k),
           hjust = 0, color = "#d73027", size = 3.5) +
  scale_x_continuous(breaks = 2:K_MAX) +
  labs(
    title    = "Silhouette analysis: mean width by k",
    subtitle = "★ = best k, ○ = k=3 (forced). Higher = better-defined clusters.",
    x        = "Number of clusters (k)",
    y        = "Mean silhouette width"
  ) +
  theme_bw(base_size = 11)

# --- Fig 2: 2D feature space coloured by k=3 cluster ---
p2 <- ggplot(candidates,
             aes(x = weighted_rain_frac, y = weighted_ros_frac,
                 color = flood_type)) +
  geom_point(alpha = 0.4, size = 0.9) +
  geom_point(data = centroids,
             aes(x = weighted_rain_frac, y = weighted_ros_frac),
             size = 5, shape = 8, color = "black") +
  scale_color_manual(values = c(
    snowmelt    = "#4dac26",
    rain_on_snow = "#762a83",
    rainfall    = "#2166ac"
  )) +
  labs(
    title    = "Cluster assignments in feature space (k=3)",
    subtitle = "★ = cluster centroids",
    x        = "Weighted rain fraction\n(total rain / total flux)",
    y        = "Weighted ROS fraction\n(co-occurrence flux / total flux)",
    color    = "Flood type"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# --- Fig 3: Marginal distributions of each feature by cluster ---
p3a <- ggplot(candidates,
              aes(x = weighted_rain_frac, fill = flood_type)) +
  geom_histogram(binwidth = 0.05, position = "identity",
                 alpha = 0.5, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c(
    snowmelt     = "#4dac26",
    rain_on_snow = "#762a83",
    rainfall     = "#2166ac"
  )) +
  labs(x = "Weighted rain fraction", y = "Count", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p3b <- ggplot(candidates,
              aes(x = weighted_ros_frac, fill = flood_type)) +
  geom_histogram(binwidth = 0.05, position = "identity",
                 alpha = 0.5, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c(
    snowmelt     = "#4dac26",
    rain_on_snow = "#762a83",
    rainfall     = "#2166ac"
  )) +
  labs(x = "Weighted ROS fraction", y = "Count", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p3 <- p3a / p3b +
  plot_annotation(title = "Marginal distributions by cluster type (k=3)")

# --- Fig 4: Cluster assignments vs elevation and drainage area ---
candidates_meta <- candidates %>%
  left_join(meta, by = "site_no")

p4 <- ggplot(candidates_meta,
             aes(x = elev_mean_m, y = drain_sqkm,
                 color = flood_type)) +
  geom_point(alpha = 0.4, size = 0.9) +
  scale_y_log10(labels = scales::comma) +
  scale_color_manual(values = c(
    snowmelt     = "#4dac26",
    rain_on_snow = "#762a83",
    rainfall     = "#2166ac"
  )) +
  labs(
    title    = "Cluster assignments vs basin characteristics",
    subtitle = "Spatial context: do clusters separate by elevation or area?",
    x        = "Mean basin elevation (m)",
    y        = "Drainage area (km², log scale)",
    color    = "Flood type"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# ============================================================================
# COMBINE AND SAVE
# ============================================================================
combined <- (p1 + p2) / (p3a + p3b) / p4 +
  plot_annotation(
    title   = "ROS vs Snowmelt Clustering",
    caption = sprintf(
      "Candidate peaks: total melt > %.0f mm over %d-day window. Noise floor: %.1f mm/day.",
      MELT_PRESENT_MM, LOOKBACK_DAYS, NOISE_FLOOR_MM
    )
  )

out_file <- file.path(out_dir, "clustering_ros_snowmelt.png")
ggsave(out_file, combined, width = 12, height = 16, dpi = 150)
message(sprintf("Saved figure to %s", out_file))

# ============================================================================
# CONSOLE SUMMARY
# ============================================================================
message("\n=== CLUSTER SUMMARY (k=3) ===")
summary_tbl <- candidates %>%
  group_by(flood_type) %>%
  summarise(
    n_peaks            = n(),
    n_sites            = n_distinct(site_no),
    median_rain_frac   = round(median(weighted_rain_frac), 3),
    median_ros_frac    = round(median(weighted_ros_frac),  3),
    median_total_flux  = round(median(total_flux), 1),
    .groups            = "drop"
  )
print(summary_tbl)