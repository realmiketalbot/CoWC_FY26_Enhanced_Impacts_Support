#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Unsupervised Clustering in {rain_frac, ros_flux_frac} Space
# ============================================================================
# Investigates whether flood peaks naturally partition into process-coherent
# clusters without relying on the hard thresholds in the 5b decision tree.
# Intended to either validate current threshold choices or suggest revisions.
#
# The analysis proceeds in three stages:
#
#   1. GEOMETRY FIRST — 2D kernel density and marginal distributions to
#      visualise the raw feature space before imposing any model. This is
#      the primary diagnostic: if there are no visible gaps or modes, GMM
#      will not find meaningful structure.
#
#   2. GMM — Gaussian mixture models (mclust) fit to the full 2D feature
#      space for k = 2:5 components, with BIC-based model selection. Cluster
#      centroids and posterior uncertainty are reported.
#
#   3. COMPARISON — Cluster assignments vs. current 5b decision tree labels,
#      to assess whether GMM recovers similar boundaries or suggests
#      different ones.
#
# SCOPE: Unregulated, data-sufficient peaks only. All four current types
# (rainfall, ROS, snowmelt, mixed) are included — the full classification
# problem is treated as a clustering question, not just the snowmelt/ROS
# boundary.
#
# INPUT:
#   - _data/pot_peaks_classified_adaptive.rds
#
# OUTPUTS (in _diagnostics/clustering/):
#   - clust_density_2d.png       : 2D KDE of {rain_frac, ros_flux_frac}
#   - clust_marginals.png        : marginal densities of each feature
#   - clust_bic.png              : BIC vs k for GMM model selection
#   - clust_gmm_assignments.png  : GMM cluster assignments in 2D space
#   - clust_uncertainty.png      : posterior assignment uncertainty
#   - clust_vs_decision_tree.png : GMM vs 5b labels (confusion heatmap)
#   - clust_summary.txt          : console summary written to file
# ============================================================================

library(tidyverse)
library(mclust)
library(dplyr)    # explicit re-attach after mclust to ensure dplyr::count() wins search path

classified_file <- "_data/pot_peaks_classified_adaptive.rds"
outdir          <- "_diagnostics/clustering"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Thresholds from 5b — used only for comparison overlay, not for clustering
RAIN_EM_THRESH   <- 0.95
SNOW_EM_THRESH   <- 0.05
ROS_COOCCUR_THRESH <- 0.50

# Cluster range to evaluate
K_MIN <- 2
K_MAX <- 5

# ============================================================================
# LOAD AND FILTER
# ============================================================================
message("Loading data...")
dat <- readRDS(classified_file)

# Unregulated, data-sufficient peaks with complete features.
# Deliberately NOT filtering on flood_type — the decision tree's threshold-
# based assignments are used only for comparison in Stage 3, not as a gate
# on what enters the clustering. Applying flood_type filters here would
# pre-impose the very boundaries we are trying to evaluate.
peaks <- dat %>%
  filter(
    !regulated,
    data_sufficient,
    !is.na(rain_frac),
    !is.na(ros_flux_frac)
  )

message(sprintf("Peaks in scope: %d across %d sites\n",
                nrow(peaks), n_distinct(peaks$site_no)))

# Feature matrix for clustering
X <- peaks %>% select(rain_frac, ros_flux_frac) %>% as.matrix()

# ============================================================================
# STAGE 1: GEOMETRY — 2D KDE AND MARGINALS
# ============================================================================
message("Stage 1: Visualising feature space geometry...")

type_colors <- c(
  rainfall     = "#e41a1c",
  snowmelt     = "#377eb8",
  ROS          = "#4daf4a",
  mixed        = "#ff7f00",
  unclassified = "#cccccc"
)

# --- Plot 1a: 2D scatter with KDE contours, coloured by current 5b type -----
p1a <- ggplot(peaks, aes(x = rain_frac, y = ros_flux_frac)) +
  geom_point(aes(colour = flood_type), alpha = 0.3, size = 0.8) +
  geom_density_2d(colour = "grey20", linewidth = 0.3, alpha = 0.6,
                  bins = 12) +
  geom_vline(xintercept = RAIN_EM_THRESH,    linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = SNOW_EM_THRESH,    linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = ROS_COOCCUR_THRESH, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  annotate("text", x = RAIN_EM_THRESH - 0.01, y = 0.97,
           label = "rain_frac = 0.95", hjust = 1, size = 2.8, colour = "grey40") +
  annotate("text", x = SNOW_EM_THRESH + 0.01, y = 0.97,
           label = "rain_frac = 0.05", hjust = 0, size = 2.8, colour = "grey40") +
  annotate("text", x = 0.97, y = ROS_COOCCUR_THRESH + 0.02,
           label = "ros_flux_frac = 0.50", hjust = 1, size = 2.8, colour = "grey40") +
  scale_colour_manual(values = type_colors) +
  labs(
    title    = "Feature space: rain_frac vs ros_flux_frac",
    subtitle = "Points coloured by current 5b classification; contours = KDE",
    x        = "rain_frac",
    y        = "ros_flux_frac",
    colour   = "5b type"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "clust_density_2d.png"), p1a,
       width = 7, height = 6, dpi = 150)
message("  Saved: clust_density_2d.png")

# --- Plot 1b: Marginal densities --------------------------------------------
p1b_rain <- ggplot(peaks, aes(x = rain_frac, fill = flood_type)) +
  geom_density(alpha = 0.45, colour = NA) +
  geom_vline(xintercept = c(SNOW_EM_THRESH, RAIN_EM_THRESH),
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_fill_manual(values = type_colors) +
  labs(title = "Marginal density: rain_frac",
       x = "rain_frac", y = "Density", fill = "5b type") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p1b_ros <- ggplot(peaks, aes(x = ros_flux_frac, fill = flood_type)) +
  geom_density(alpha = 0.45, colour = NA) +
  geom_vline(xintercept = ROS_COOCCUR_THRESH,
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_fill_manual(values = type_colors) +
  labs(title = "Marginal density: ros_flux_frac",
       x = "ros_flux_frac", y = "Density", fill = "5b type") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# Stack marginals
library(patchwork)
p1b <- p1b_rain / p1b_ros +
  plot_annotation(title = "Marginal feature distributions by 5b flood type")

ggsave(file.path(outdir, "clust_marginals.png"), p1b,
       width = 7, height = 8, dpi = 150)
message("  Saved: clust_marginals.png")

# ============================================================================
# STAGE 2: GMM
# ============================================================================
message("\nStage 2: Fitting GMMs (k = ", K_MIN, " to ", K_MAX, ")...")

# mclust fits GMM with full BIC-based covariance structure selection
# VVV = unconstrained (variable volume, shape, orientation) — most flexible
set.seed(42)
gmm_fit <- Mclust(X, G = K_MIN:K_MAX, modelNames = "VVV", verbose = FALSE)

if (is.null(gmm_fit)) {
  # Fallback to all covariance types if VVV fails
  message("  VVV model failed, trying all covariance types...")
  gmm_fit <- Mclust(X, G = K_MIN:K_MAX, verbose = FALSE)
}

best_k <- gmm_fit$G
message(sprintf("  Best k = %d (BIC = %.1f)", best_k, gmm_fit$bic))

# --- Plot 2a: BIC curve ------------------------------------------------------
# gmm_fit$BIC is an mclustBIC matrix (k x model); extract as plain matrix first
bic_mat <- unclass(gmm_fit$BIC)  # strip mclustBIC class to get plain matrix
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
    subtitle = "Higher BIC = better fit; dashed line = selected model",
    x        = "Number of components (k)",
    y        = "BIC",
    colour   = "Covariance"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "clust_bic.png"), p2a,
       width = 7, height = 5, dpi = 150)
message("  Saved: clust_bic.png")

# Add GMM assignments and uncertainty to peaks
peaks <- peaks %>%
  mutate(
    gmm_cluster     = factor(gmm_fit$classification),
    gmm_uncertainty = gmm_fit$uncertainty
  )

# --- Plot 2b: GMM assignments in 2D space ------------------------------------
p2b <- ggplot(peaks, aes(x = rain_frac, y = ros_flux_frac, colour = gmm_cluster)) +
  geom_point(alpha = 0.4, size = 0.9) +
  geom_vline(xintercept = c(SNOW_EM_THRESH, RAIN_EM_THRESH),
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = ROS_COOCCUR_THRESH,
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  # Cluster centroids
  geom_point(
    data = as.data.frame(t(gmm_fit$parameters$mean)) %>%
      mutate(gmm_cluster = factor(seq_len(best_k))),
    aes(x = rain_frac, y = ros_flux_frac),
    shape = 3, size = 5, stroke = 1.5, colour = "black"
  ) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    title    = sprintf("GMM cluster assignments (k = %d)", best_k),
    subtitle = "Crosses = cluster centroids; dashed lines = 5b thresholds",
    x        = "rain_frac",
    y        = "ros_flux_frac",
    colour   = "GMM cluster"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "clust_gmm_assignments.png"), p2b,
       width = 7, height = 6, dpi = 150)
message("  Saved: clust_gmm_assignments.png")

# --- Plot 2c: Posterior uncertainty ------------------------------------------
p2c <- ggplot(peaks, aes(x = rain_frac, y = ros_flux_frac, colour = gmm_uncertainty)) +
  geom_point(alpha = 0.6, size = 0.9) +
  geom_vline(xintercept = c(SNOW_EM_THRESH, RAIN_EM_THRESH),
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = ROS_COOCCUR_THRESH,
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_colour_viridis_c(option = "magma", direction = -1,
                         name = "Uncertainty") +
  labs(
    title    = "GMM posterior assignment uncertainty",
    subtitle = "High uncertainty = ambiguous cluster membership",
    x        = "rain_frac",
    y        = "ros_flux_frac"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(outdir, "clust_uncertainty.png"), p2c,
       width = 7, height = 6, dpi = 150)
message("  Saved: clust_uncertainty.png")

# ============================================================================
# STAGE 3: GMM VS DECISION TREE COMPARISON
# ============================================================================
message("\nStage 3: Comparing GMM assignments to 5b decision tree...")

# Confusion heatmap: rows = GMM cluster, cols = 5b flood type
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
    title    = "GMM clusters vs. 5b decision tree classification",
    subtitle = "Row percentages: composition of each GMM cluster by 5b type",
    x        = "5b flood type",
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
summary_lines <- c(
  "=== Clustering Diagnostic Summary ===",
  sprintf("Peaks in scope  : %d (%d sites, unregulated, data-sufficient)",
          nrow(peaks), n_distinct(peaks$site_no)),
  sprintf("Best GMM k      : %d  (model: %s, BIC = %.1f)",
          best_k, gmm_fit$modelName, gmm_fit$bic),
  "",
  "--- Cluster centroids ---",
  capture.output(print(round(t(gmm_fit$parameters$mean), 3))),
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
  "--- GMM vs 5b confusion (counts) ---",
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

# Also echo to console
cat(paste(summary_lines, collapse = "\n"), "\n")

message("\nDiagnostic complete. Outputs in: ", outdir)