#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Unsupervised Clustering — Non-Rainfall Population
# ============================================================================
# Restricts clustering to the non-rainfall population (rain_frac < 0.95),
# treating the rainfall end-member boundary as already established. This
# avoids the class imbalance problem observed in the full-population analysis,
# where the GMM preferentially subdivided the dominant rainfall population
# along seasonal lines rather than recovering process-coherent clusters.
#
# Feature space (4D):
#   ros_flux_frac : fraction of total flux on rain-melt co-occurrence days
#   swe_log1p     : log1p(swe_mean_mm) — antecedent SWE, log-transformed
#   doy_sin       : sin(2π × DOY / 365) — circular seasonality encoding
#   doy_cos       : cos(2π × DOY / 365) — circular seasonality encoding
#
# rain_frac is excluded from clustering since it has limited variance within
# the non-rainfall population and its boundary with rainfall is treated as
# settled by the 6b decision tree.
#
# The analysis proceeds in three stages:
#
#   1. GEOMETRY — marginal distributions of all 4 features within the
#      non-rainfall population, coloured by 6b type. Primary visual check
#      for separable modes before GMM.
#
#   2. GMM — fit to the 4D feature space for k = 2:6, BIC-based model
#      selection.
#
#   3. COMPARISON — GMM assignments vs. 6b decision tree labels (snowmelt,
#      ROS, mixed, unclassified), with 2D projection plots.
#
# SCOPE: Unregulated, data-sufficient peaks with rain_frac < RAIN_EM_THRESH
# and complete values for all 4 features.
#
# INPUT:
#   - _data/pot_peaks_classified_adaptive.rds
#
# OUTPUTS (in _diagnostics/clustering_nonrainfall/):
#   - clust_marginals.png          : marginals for all 4 features
#   - clust_swe_doy.png            : SWE vs DOY scatter by 6b type
#   - clust_ros_rain.png           : ros_flux_frac vs rain_frac scatter by 6b type
#   - clust_bic.png                : BIC vs k
#   - clust_gmm_ros_rain.png       : GMM assignments: ros_flux_frac vs rain_frac
#   - clust_gmm_swe_doy.png        : GMM assignments: swe_mean_mm vs DOY
#   - clust_uncertainty.png        : posterior uncertainty
#   - clust_vs_decision_tree.png   : confusion heatmap
#   - clust_summary.txt            : written summary
# ============================================================================

library(tidyverse)
library(mclust)
library(dplyr)       # explicit re-attach after mclust to prevent count() masking
library(patchwork)
library(cowplot)

classified_file <- "_data/pot_peaks_classified_adaptive.rds"
outdir          <- "_diagnostics/clustering_nonrainfall"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# 6b decision tree thresholds — used only for comparison overlay
RAIN_EM_THRESH     <- 0.95
SNOW_EM_THRESH     <- 0.05
ROS_COOCCUR_THRESH <- 0.50

K_MIN <- 2
K_MAX <- 6

type_colors <- c(
  rainfall     = "#e41a1c",
  snowmelt     = "#377eb8",
  ROS          = "#4daf4a",
  mixed        = "#ff7f00",
  unclassified = "#cccccc"
)

# ============================================================================
# LOAD AND PREPARE
# ============================================================================
message("Loading data...")
dat <- readRDS(classified_file)

# Unregulated, data-sufficient, NON-RAINFALL peaks with complete features.
# rain_frac < RAIN_EM_THRESH excludes the rainfall end-member population,
# whose boundary is treated as settled. Decision tree labels used only in
# Stage 3 comparison.
peaks <- dat %>%
  filter(
    !regulated,
    data_sufficient,
    rain_frac < RAIN_EM_THRESH,
    !is.na(rain_frac),
    !is.na(ros_flux_frac),
    !is.na(swe_mean_mm)
  ) %>%
  mutate(
    doy     = as.integer(format(peak_date, "%j")),
    doy_sin = sin(2 * pi * doy / 365),
    doy_cos = cos(2 * pi * doy / 365)
  )

message(sprintf("Non-rainfall peaks in scope: %d across %d sites\n",
                nrow(peaks), n_distinct(peaks$site_no)))

# Feature matrix — 4D (rain_frac excluded; limited variance in this population)
peaks <- peaks %>%
  mutate(swe_log1p = log1p(swe_mean_mm))

X <- peaks %>%
  select(ros_flux_frac, swe_log1p, doy_sin, doy_cos) %>%
  as.matrix()

# ============================================================================
# STAGE 1: GEOMETRY — MARGINALS AND FEATURE SCATTERS
# ============================================================================
message("Stage 1: Visualising feature space geometry...")

# --- Plot 1a: Marginals for all 4 features -----------------------------------
p_ros <- ggplot(peaks, aes(x = ros_flux_frac, fill = flood_type)) +
  geom_density(alpha = 0.45, colour = NA) +
  geom_vline(xintercept = ROS_COOCCUR_THRESH,
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_fill_manual(values = type_colors) +
  labs(title = "ros_flux_frac", x = NULL, y = "Density", fill = "6b type") +
  theme_bw(base_size = 10) + theme(legend.position = "none")

p_swe <- ggplot(peaks, aes(x = swe_mean_mm, fill = flood_type)) +
  geom_density(alpha = 0.45, colour = NA) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 10, 50, 100, 250, 500),
                     labels = c("0", "10", "50", "100", "250", "500")) +
  scale_fill_manual(values = type_colors) +
  labs(title = "swe_mean_mm (log1p scale)", x = NULL, y = "Density",
       fill = "6b type") +
  theme_bw(base_size = 10) + theme(legend.position = "none")

p_doy <- ggplot(peaks, aes(x = doy, fill = flood_type)) +
  geom_density(alpha = 0.45, colour = NA) +
  scale_x_continuous(breaks = c(1, 91, 182, 274),
                     labels = c("Jan", "Apr", "Jul", "Oct")) +
  scale_fill_manual(values = type_colors) +
  labs(title = "Peak DOY (seasonality)", x = NULL, y = "Density",
       fill = "6b type") +
  theme_bw(base_size = 10) + theme(legend.position = "none")

p_rain <- ggplot(peaks, aes(x = rain_frac, fill = flood_type)) +
  geom_density(alpha = 0.45, colour = NA) +
  geom_vline(xintercept = SNOW_EM_THRESH,
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_fill_manual(values = type_colors) +
  labs(title = "rain_frac (reference only; not used in GMM)",
       x = NULL, y = "Density", fill = "6b type") +
  theme_bw(base_size = 10) + theme(legend.position = "none")

# Shared legend
p_legend <- ggplot(peaks, aes(x = ros_flux_frac, fill = flood_type)) +
  geom_density(alpha = 0.45, colour = NA) +
  scale_fill_manual(values = type_colors, name = "6b type") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

shared_legend <- cowplot::get_legend(p_legend)

p1a <- (p_ros | p_swe) / (p_doy | p_rain) /
  cowplot::plot_grid(shared_legend) +
  plot_annotation(
    title    = "Marginal feature distributions — non-rainfall population",
    subtitle = sprintf("rain_frac < %.2f; dashed lines = 6b decision tree thresholds",
                       RAIN_EM_THRESH)
  )

ggsave(file.path(outdir, "clust_marginals.png"), p1a,
       width = 10, height = 10, dpi = 150)
message("  Saved: clust_marginals.png")

# --- Plot 1b: SWE vs DOY scatter ---------------------------------------------
p1b <- ggplot(peaks, aes(x = doy, y = swe_mean_mm, colour = flood_type)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_density_2d(colour = "grey20", linewidth = 0.3, alpha = 0.5, bins = 10) +
  scale_colour_manual(values = type_colors) +
  scale_x_continuous(breaks = c(1, 91, 182, 274),
                     labels = c("Jan", "Apr", "Jul", "Oct")) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 10, 50, 100, 250, 500),
                     labels = c("0", "10", "50", "100", "250", "500")) +
  labs(
    title    = "SWE vs peak DOY — non-rainfall population",
    subtitle = "y-axis log1p scale; contours = KDE",
    x        = "Day of year",
    y        = "swe_mean_mm",
    colour   = "6b type"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "clust_swe_doy.png"), p1b,
       width = 7, height = 6, dpi = 150)
message("  Saved: clust_swe_doy.png")

# --- Plot 1c: ros_flux_frac vs rain_frac scatter -----------------------------
p1c <- ggplot(peaks, aes(x = rain_frac, y = ros_flux_frac, colour = flood_type)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_density_2d(colour = "grey20", linewidth = 0.3, alpha = 0.5, bins = 10) +
  geom_vline(xintercept = SNOW_EM_THRESH,
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = ROS_COOCCUR_THRESH,
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_colour_manual(values = type_colors) +
  labs(
    title    = "ros_flux_frac vs rain_frac — non-rainfall population",
    subtitle = "Dashed lines = 6b thresholds; contours = KDE",
    x        = "rain_frac",
    y        = "ros_flux_frac",
    colour   = "6b type"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "clust_ros_rain.png"), p1c,
       width = 7, height = 6, dpi = 150)
message("  Saved: clust_ros_rain.png")

# ============================================================================
# STAGE 2: GMM ON 4D FEATURE SPACE
# ============================================================================
message(sprintf("\nStage 2: Fitting GMMs (k = %d to %d) on 4D feature space...",
                K_MIN, K_MAX))
message("  Features: ros_flux_frac, swe_log1p, doy_sin, doy_cos")

set.seed(42)
gmm_fit <- Mclust(X, G = K_MIN:K_MAX, verbose = FALSE)

if (is.null(gmm_fit)) stop("Mclust failed to fit any model.")

best_k <- gmm_fit$G
message(sprintf("  Best k = %d (model: %s, BIC = %.1f)",
                best_k, gmm_fit$modelName, gmm_fit$bic))

# --- Plot 2a: BIC curve ------------------------------------------------------
bic_mat <- unclass(gmm_fit$BIC)
bic_df <- as.data.frame(bic_mat) %>%
  mutate(k = as.integer(rownames(bic_mat))) %>%
  pivot_longer(-k, names_to = "model", values_to = "BIC") %>%
  filter(!is.na(BIC))

p2a <- ggplot(bic_df, aes(x = k, y = BIC, colour = model, group = model)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_vline(xintercept = best_k, linetype = "dashed", colour = "grey30") +
  annotate("text", x = best_k + 0.1, y = max(bic_df$BIC, na.rm = TRUE),
           label = sprintf("Best k = %d", best_k),
           hjust = 0, size = 3, colour = "grey30") +
  labs(
    title    = "GMM model selection: BIC by k and covariance structure",
    subtitle = "4D feature space (non-rainfall); higher BIC = better fit",
    x        = "Number of components (k)",
    y        = "BIC",
    colour   = "Covariance"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "clust_bic.png"), p2a,
       width = 7, height = 5, dpi = 150)
message("  Saved: clust_bic.png")

# Attach GMM assignments and uncertainty
peaks <- peaks %>%
  mutate(
    gmm_cluster     = factor(gmm_fit$classification),
    gmm_uncertainty = gmm_fit$uncertainty
  )

# --- Plot 2b: GMM assignments — ros_flux_frac vs rain_frac ------------------
# rain_frac retained as reference axis even though excluded from GMM
centroids_2d <- as.data.frame(t(gmm_fit$parameters$mean)) %>%
  select(ros_flux_frac) %>%
  mutate(
    rain_frac   = NA_real_,   # not in GMM; centroid not meaningful on this axis
    gmm_cluster = factor(seq_len(best_k))
  )

p2b <- ggplot(peaks, aes(x = rain_frac, y = ros_flux_frac,
                          colour = gmm_cluster)) +
  geom_point(alpha = 0.35, size = 0.8) +
  geom_vline(xintercept = SNOW_EM_THRESH,
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = ROS_COOCCUR_THRESH,
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    title    = sprintf("GMM assignments (k = %d): ros_flux_frac vs rain_frac",
                       best_k),
    subtitle = "rain_frac reference only (not in GMM); dashed = 6b thresholds",
    x        = "rain_frac",
    y        = "ros_flux_frac",
    colour   = "GMM cluster"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "clust_gmm_ros_rain.png"), p2b,
       width = 7, height = 6, dpi = 150)
message("  Saved: clust_gmm_ros_rain.png")

# --- Plot 2c: GMM assignments — swe_mean_mm vs DOY --------------------------
centroids_swe_doy <- as.data.frame(t(gmm_fit$parameters$mean)) %>%
  select(swe_log1p, doy_sin, doy_cos) %>%
  mutate(
    gmm_cluster = factor(seq_len(best_k)),
    # Recover DOY from sin/cos centroid (atan2 gives angle, convert to DOY)
    doy_centroid = (atan2(doy_sin, doy_cos) %% (2 * pi)) / (2 * pi) * 365,
    swe_centroid = expm1(swe_log1p)
  )

p2c <- ggplot(peaks, aes(x = doy, y = swe_mean_mm, colour = gmm_cluster)) +
  geom_point(alpha = 0.35, size = 0.8) +
  geom_point(data = centroids_swe_doy,
             aes(x = doy_centroid, y = swe_centroid),
             shape = 3, size = 5, stroke = 1.5, colour = "black") +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(1, 91, 182, 274),
                     labels = c("Jan", "Apr", "Jul", "Oct")) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 10, 50, 100, 250, 500),
                     labels = c("0", "10", "50", "100", "250", "500")) +
  labs(
    title    = sprintf("GMM assignments (k = %d): swe_mean_mm vs DOY", best_k),
    subtitle = "Crosses = cluster centroids; y-axis log1p scale",
    x        = "Day of year",
    y        = "swe_mean_mm",
    colour   = "GMM cluster"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "clust_gmm_swe_doy.png"), p2c,
       width = 7, height = 6, dpi = 150)
message("  Saved: clust_gmm_swe_doy.png")

# --- Plot 2d: Posterior uncertainty ------------------------------------------
p2d <- ggplot(peaks, aes(x = rain_frac, y = ros_flux_frac,
                          colour = gmm_uncertainty)) +
  geom_point(alpha = 0.6, size = 0.9) +
  geom_vline(xintercept = SNOW_EM_THRESH,
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = ROS_COOCCUR_THRESH,
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_colour_viridis_c(option = "magma", direction = -1,
                         name = "Uncertainty") +
  labs(
    title    = "GMM posterior assignment uncertainty",
    subtitle = "Projected onto rain_frac × ros_flux_frac (reference axes); high = ambiguous",
    x        = "rain_frac",
    y        = "ros_flux_frac"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(outdir, "clust_uncertainty.png"), p2d,
       width = 7, height = 6, dpi = 150)
message("  Saved: clust_uncertainty.png")

# ============================================================================
# STAGE 3: GMM VS DECISION TREE COMPARISON
# ============================================================================
message("\nStage 3: Comparing GMM assignments to 6b decision tree (non-rainfall types)...")

confusion <- peaks %>%
  dplyr::count(gmm_cluster, flood_type) %>%
  group_by(gmm_cluster) %>%
  mutate(pct = n / sum(n) * 100)

p3 <- ggplot(confusion,
             aes(x = flood_type, y = gmm_cluster, fill = pct)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.0f%%\n(n=%d)", pct, n)),
            size = 3, lineheight = 0.9) +
  scale_fill_distiller(palette = "Blues", direction = 1,
                       name = "% of\nGMM cluster") +
  labs(
    title    = "GMM clusters vs. 6b decision tree — non-rainfall population",
    subtitle = "Row percentages: composition of each GMM cluster by 6b type",
    x        = "6b flood type",
    y        = "GMM cluster"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

ggsave(file.path(outdir, "clust_vs_decision_tree.png"), p3,
       width = 7, height = 4 + best_k * 0.4, dpi = 150)
message("  Saved: clust_vs_decision_tree.png")

# ============================================================================
# SUMMARY TO FILE
# ============================================================================
centroid_df <- as.data.frame(t(gmm_fit$parameters$mean))
names(centroid_df) <- c("ros_flux_frac", "swe_log1p", "doy_sin", "doy_cos")
centroid_df <- centroid_df %>%
  mutate(
    swe_mean_mm  = expm1(swe_log1p),
    doy_approx   = round((atan2(doy_sin, doy_cos) %% (2 * pi)) /
                           (2 * pi) * 365),
    cluster      = seq_len(best_k)
  ) %>%
  select(cluster, ros_flux_frac, swe_mean_mm, doy_approx)

summary_lines <- c(
  "=== Clustering Diagnostic Summary (Non-Rainfall Population) ===",
  sprintf("Peaks in scope : %d (%d sites, unregulated, rain_frac < %.2f)",
          nrow(peaks), n_distinct(peaks$site_no), RAIN_EM_THRESH),
  sprintf("Features       : ros_flux_frac, swe_log1p, doy_sin, doy_cos"),
  sprintf("Best GMM k     : %d  (model: %s, BIC = %.1f)",
          best_k, gmm_fit$modelName, gmm_fit$bic),
  "",
  "--- Cluster centroids (interpretable units) ---",
  capture.output(print(round(centroid_df, 2))),
  "",
  "--- Cluster sizes ---",
  capture.output(print(table(peaks$gmm_cluster))),
  "",
  "--- Uncertainty summary ---",
  capture.output(print(summary(peaks$gmm_uncertainty))),
  sprintf("  Peaks with uncertainty > 0.20 : %d (%.1f%%)",
          sum(peaks$gmm_uncertainty > 0.20),
          100 * mean(peaks$gmm_uncertainty > 0.20)),
  sprintf("  Peaks with uncertainty > 0.40 : %d (%.1f%%)",
          sum(peaks$gmm_uncertainty > 0.40),
          100 * mean(peaks$gmm_uncertainty > 0.40)),
  "",
  "--- GMM vs 6b confusion (counts) ---",
  capture.output(
    print(
      peaks %>%
        dplyr::count(gmm_cluster, flood_type) %>%
        pivot_wider(names_from = flood_type, values_from = n, values_fill = 0)
    )
  )
)

writeLines(summary_lines, file.path(outdir, "clust_summary.txt"))
message("  Saved: clust_summary.txt")

cat(paste(summary_lines, collapse = "\n"), "\n")

message("\nDiagnostic complete. Outputs in: ", outdir)