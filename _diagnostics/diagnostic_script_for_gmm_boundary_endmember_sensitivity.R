#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: End-Member Sensitivity of GMM-Derived B1/B2 Boundaries
# ============================================================================
# Investigates whether the dense vertical stripes at rain_frac == 0 and
# rain_frac == 1 (events deterministically classified as pure snowmelt or
# pure rainfall under the existing ROS_RAIN_FRAC_MIN / MAX guards) are
# distorting the GMM component covariances and/or the sequential logistic
# regression that produces the linear B1 and B2 boundary coefficients
# currently hardcoded in 6b.
#
# Three configurations are evaluated:
#
#   A. baseline            : GMM on all peaks; logistic on all peaks.
#                            Expected to reproduce the current B1/B2 values
#                            (0.4973, 5.3433, 0.1078, 0.3503) within rounding.
#
#   B. boundary_filter     : GMM on all peaks (so component labels are stable);
#                            logistic on rain_frac in (END_MEMBER_LO, END_MEMBER_HI).
#                            Tests whether the end-member stripes are leveraging
#                            the boundary fit alone.
#
#   C. gmm_and_boundary_filter : GMM and logistic both fit on the interior
#                            subset only. Tests whether the GMM component
#                            covariances themselves are distorted by the
#                            stripes (which would shift labels for interior
#                            events).
#
# For each configuration the script reports:
#   - GMM component count, model name, BIC
#   - Per-component median rain_frac and SWE, and merged 3-class assignment
#   - Sequential logistic regression coefficients in the slope/intercept form
#     that 6b consumes:
#         B1 :  log10(SWE+1) = B1_INTERCEPT + B1_SLOPE * rain_frac
#         B2 :  rain_frac    = B2_INTERCEPT + B2_SLOPE * log10(SWE+1)
#   - Per-class precision and recall against the merged GMM labels
#   - Count of events excluded by the end-member filter
#
# A 3-panel comparison figure overlays the baseline data with the fitted
# B1/B2 lines from each configuration.
#
# INPUTS:
#   _data/pot_peaks_classified_adaptive.rds  (rain_frac, swe_mean_mm, regulated)
#   _data/site_metadata.rds                  (regulated flag)
#
# OUTPUTS (in _diagnostics/threshold_evaluation/gmm_rain_frac_swe/end_member_sensitivity/):
#   coefficients_summary.csv      — B1/B2 across configurations
#   classification_agreement.csv  — per-class precision/recall per config
#   boundary_comparison.png       — 3-panel overlay
#   gmm_summary.txt               — console output
# ============================================================================

# Explicit namespace qualification: mclust masks dplyr::select, purrr::map,
# and dplyr::count. Per project convention, all calls are qualified.
suppressPackageStartupMessages({
  library(tidyverse)
  library(mclust)
  library(patchwork)
})

# ============================================================================
# CONFIGURATION
# ============================================================================
proj_root <- "/Users/mtalbot/Desktop/GitHub/CoWC_FY26_Enhanced_Impacts/"
data_dir  <- file.path(proj_root, "_data")
diag_dir  <- file.path(proj_root, "_diagnostics", "threshold_evaluation",
                       "gmm_rain_frac_swe", "end_member_sensitivity")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

peaks_file <- file.path(data_dir, "pot_peaks_classified_adaptive.rds")
meta_file  <- file.path(data_dir, "site_metadata.rds")

# Match constants used in chunk 1 of the original GMM analysis
SWE_LOG_OFFSET <- 1.0
BIC_MODELS     <- c("EEE", "VVV", "EEV")
GMM_G_RANGE    <- 2:4

# End-member filter thresholds (matches ROS_RAIN_FRAC_MIN / MAX in 6b)
END_MEMBER_LO <- 0.01   # rain_frac <= this is deterministically snowmelt
END_MEMBER_HI <- 0.99   # rain_frac >= this is deterministically rainfall

# Reference values from current 6b (for parity check on baseline)
REF_B1_INT <- 0.4973
REF_B1_SLP <- 5.3433
REF_B2_INT <- 0.1078
REF_B2_SLP <- 0.3503

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading data...")
adaptive <- readRDS(peaks_file)
meta     <- readRDS(meta_file)

unreg_sites <- meta %>%
  dplyr::filter(!regulated) %>%
  dplyr::pull(site_no)

# One row per peak; restrict to unregulated sites with valid features
peaks_full <- adaptive %>%
  dplyr::distinct(site_no, peak_id, .keep_all = TRUE) %>%
  dplyr::filter(site_no %in% unreg_sites,
                !is.na(rain_frac),
                !is.na(swe_mean_mm)) %>%
  dplyr::mutate(log10_swe = log10(swe_mean_mm + SWE_LOG_OFFSET))

n_full <- nrow(peaks_full)
n_low  <- sum(peaks_full$rain_frac <= END_MEMBER_LO)
n_high <- sum(peaks_full$rain_frac >= END_MEMBER_HI)
n_int  <- n_full - n_low - n_high

message(sprintf("Peaks (unregulated, rain_frac & SWE non-NA): %d", n_full))
message(sprintf("  rain_frac <= %.2f : %d (%.1f%%)",
                END_MEMBER_LO, n_low, 100 * n_low / n_full))
message(sprintf("  rain_frac >= %.2f : %d (%.1f%%)",
                END_MEMBER_HI, n_high, 100 * n_high / n_full))
message(sprintf("  interior          : %d (%.1f%%)",
                n_int, 100 * n_int / n_full))

# ============================================================================
# HELPERS
# ============================================================================

# Fit GMM on a given peak subset using BIC selection over GMM_G_RANGE and
# BIC_MODELS. Returns a list with the fit, the scaling parameters used
# during fitting, and the merged 3-class component map.
fit_gmm_with_labels <- function(df) {
  rf_mean <- mean(df$rain_frac)
  rf_sd   <- sd(df$rain_frac)
  sw_mean <- mean(df$log10_swe)
  sw_sd   <- sd(df$log10_swe)

  X <- cbind(
    rain_frac_scaled = (df$rain_frac - rf_mean) / rf_sd,
    swe_log_scaled   = (df$log10_swe - sw_mean) / sw_sd
  )

  fit <- mclust::Mclust(X, G = GMM_G_RANGE,
                        modelNames = BIC_MODELS, verbose = FALSE)

  comp_summary <- tibble::tibble(
    gmm_component = seq_len(fit$G),
    rain_frac_med = vapply(seq_len(fit$G), function(g) {
      median(df$rain_frac[fit$classification == g])
    }, numeric(1)),
    swe_log_med   = vapply(seq_len(fit$G), function(g) {
      median(df$log10_swe[fit$classification == g])
    }, numeric(1)),
    n             = as.integer(table(fit$classification))
  ) %>%
    dplyr::arrange(rain_frac_med)

  # Merged-class mapping: lowest rain_frac = snowmelt, second = ROS,
  # remaining = rainfall (merged). Only valid when fit$G == 4; otherwise
  # we map proportionally and warn.
  if (fit$G == 4L) {
    comp_summary$merged_class <- c("snowmelt", "ROS", "rainfall", "rainfall")
  } else if (fit$G == 3L) {
    comp_summary$merged_class <- c("snowmelt", "ROS", "rainfall")
    warning("GMM selected G=3; no rainfall components to merge.")
  } else if (fit$G == 2L) {
    comp_summary$merged_class <- c("snowmelt", "rainfall")
    warning("GMM selected G=2; ROS class absent — boundary fit will fail.")
  }

  list(
    fit          = fit,
    comp_summary = comp_summary,
    rf_mean = rf_mean, rf_sd = rf_sd,
    sw_mean = sw_mean, sw_sd = sw_sd
  )
}

# Given a fitted GMM and a peak subset, assign merged 3-class labels.
# If the peak subset differs from the GMM training set (config B), the
# subset is projected through the trained GMM via predict().
assign_merged_labels <- function(gmm_obj, df) {
  fit <- gmm_obj$fit
  X_new <- cbind(
    rain_frac_scaled = (df$rain_frac - gmm_obj$rf_mean) / gmm_obj$rf_sd,
    swe_log_scaled   = (df$log10_swe - gmm_obj$sw_mean) / gmm_obj$sw_sd
  )
  pred <- predict(fit, newdata = X_new)
  comp <- pred$classification

  comp_to_class <- setNames(
    gmm_obj$comp_summary$merged_class,
    as.character(gmm_obj$comp_summary$gmm_component)
  )
  factor(comp_to_class[as.character(comp)],
         levels = c("snowmelt", "ROS", "rainfall"))
}

# Fit sequential logistic regression on (rain_frac, log10_swe) to recover
# linear B1 and B2 boundaries.
#
# Decision boundary of glm(y ~ x1 + x2, family=binomial) is the line where
# beta0 + beta1*x1 + beta2*x2 = 0. Solving for the desired form:
#   B1 (snowmelt vs. ROS+rainfall, predict snowmelt=1):
#     log10_swe = -(beta0 + beta1*rain_frac) / beta2
#                = (-beta0/beta2) + (-beta1/beta2) * rain_frac
#     => B1_INT = -beta0/beta2,   B1_SLP = -beta1/beta2
#
#   B2 (ROS vs. rainfall on non-snowmelt subset, predict ROS=1):
#     rain_frac = -(beta0 + beta2*log10_swe) / beta1
#                = (-beta0/beta1) + (-beta2/beta1) * log10_swe
#     => B2_INT = -beta0/beta1,   B2_SLP = -beta2/beta1
#
fit_sequential_logistic <- function(df_with_labels) {
  d <- df_with_labels %>%
    dplyr::mutate(
      y_b1 = as.integer(merged_class == "snowmelt"),
      y_b2 = as.integer(merged_class == "ROS")
    )

  # B1: snowmelt vs. everything else
  m1 <- glm(y_b1 ~ rain_frac + log10_swe,
            data = d, family = binomial(link = "logit"))
  c1 <- coef(m1)
  B1_INT <- unname(-c1["(Intercept)"] / c1["log10_swe"])
  B1_SLP <- unname(-c1["rain_frac"]   / c1["log10_swe"])

  # B2: ROS vs. rainfall on non-snowmelt subset only
  d_b2 <- d %>% dplyr::filter(merged_class != "snowmelt")
  m2 <- glm(y_b2 ~ rain_frac + log10_swe,
            data = d_b2, family = binomial(link = "logit"))
  c2 <- coef(m2)
  B2_INT <- unname(-c2["(Intercept)"] / c2["rain_frac"])
  B2_SLP <- unname(-c2["log10_swe"]   / c2["rain_frac"])

  list(
    B1_INT = B1_INT, B1_SLP = B1_SLP,
    B2_INT = B2_INT, B2_SLP = B2_SLP,
    n_b1 = nrow(d), n_b2 = nrow(d_b2),
    m1 = m1, m2 = m2
  )
}

# Apply fitted boundaries to classify a peak subset using the sequential
# decision tree from 6b.
classify_with_boundaries <- function(df, B1_INT, B1_SLP, B2_INT, B2_SLP) {
  df %>%
    dplyr::mutate(
      pred_class = dplyr::case_when(
        log10_swe > B1_INT + B1_SLP * rain_frac     ~ "snowmelt",
        rain_frac < END_MEMBER_LO                   ~ "snowmelt",
        rain_frac > END_MEMBER_HI                   ~ "rainfall",
        rain_frac < B2_INT + B2_SLP * log10_swe     ~ "ROS",
        TRUE                                        ~ "rainfall"
      ),
      pred_class = factor(pred_class,
                          levels = c("snowmelt", "ROS", "rainfall"))
    )
}

# Per-class precision and recall against merged GMM labels
agreement_metrics <- function(df) {
  classes <- c("snowmelt", "ROS", "rainfall")
  purrr::map_dfr(classes, function(cls) {
    tp <- sum(df$pred_class == cls & df$merged_class == cls)
    fp <- sum(df$pred_class == cls & df$merged_class != cls)
    fn <- sum(df$pred_class != cls & df$merged_class == cls)
    tibble::tibble(
      class     = cls,
      n_true    = sum(df$merged_class == cls),
      precision = if (tp + fp > 0) tp / (tp + fp) else NA_real_,
      recall    = if (tp + fn > 0) tp / (tp + fn) else NA_real_
    )
  })
}

# ============================================================================
# RUN THREE CONFIGURATIONS
# ============================================================================
peaks_interior <- peaks_full %>%
  dplyr::filter(rain_frac > END_MEMBER_LO,
                rain_frac < END_MEMBER_HI)

run_config <- function(label, gmm_data, boundary_data) {
  message(sprintf("\n=== Configuration: %s ===", label))
  message(sprintf("  GMM fit on n = %d", nrow(gmm_data)))
  message(sprintf("  Boundary fit on n = %d", nrow(boundary_data)))

  gmm_obj <- fit_gmm_with_labels(gmm_data)
  message(sprintf("  GMM: G = %d, model = %s",
                  gmm_obj$fit$G, gmm_obj$fit$modelName))

  comp_print <- gmm_obj$comp_summary %>%
    dplyr::mutate(
      rain_frac_med = round(rain_frac_med, 3),
      swe_log_med   = round(swe_log_med,   3)
    )
  message("  Component summary (ascending median rain_frac):")
  message(paste(capture.output(print(comp_print)), collapse = "\n"))

  # Assign labels to the boundary-fit subset using the trained GMM
  bd <- boundary_data %>%
    dplyr::mutate(merged_class = assign_merged_labels(gmm_obj, boundary_data))

  bnd <- fit_sequential_logistic(bd)
  message(sprintf(
    "  B1: log10(SWE+1) = %.4f + %.4f * rain_frac  (n = %d)",
    bnd$B1_INT, bnd$B1_SLP, bnd$n_b1))
  message(sprintf(
    "  B2: rain_frac    = %.4f + %.4f * log10(SWE+1)  (n = %d)",
    bnd$B2_INT, bnd$B2_SLP, bnd$n_b2))

  # Apply boundaries back to the FULL dataset and score against
  # full-dataset merged labels (so configs are comparable on the same
  # population).
  full_labeled <- peaks_full %>%
    dplyr::mutate(merged_class = assign_merged_labels(gmm_obj, peaks_full))
  classified <- classify_with_boundaries(full_labeled,
                                         bnd$B1_INT, bnd$B1_SLP,
                                         bnd$B2_INT, bnd$B2_SLP)
  metrics <- agreement_metrics(classified)
  message("  Per-class agreement (boundary classification vs. merged GMM):")
  message(paste(capture.output(print(metrics, n = Inf)), collapse = "\n"))

  list(
    label    = label,
    gmm      = gmm_obj,
    boundary = bnd,
    metrics  = metrics
  )
}

cfg_A <- run_config("A_baseline",
                    gmm_data = peaks_full,     boundary_data = peaks_full)
cfg_B <- run_config("B_boundary_filter",
                    gmm_data = peaks_full,     boundary_data = peaks_interior)
cfg_C <- run_config("C_gmm_and_boundary_filter",
                    gmm_data = peaks_interior, boundary_data = peaks_interior)

# ============================================================================
# PARITY CHECK ON BASELINE
# ============================================================================
message("\n=== Parity check: configuration A vs. current 6b values ===")
parity <- tibble::tibble(
  coef    = c("B1_INT", "B1_SLP", "B2_INT", "B2_SLP"),
  current = c(REF_B1_INT, REF_B1_SLP, REF_B2_INT, REF_B2_SLP),
  diag_A  = c(cfg_A$boundary$B1_INT, cfg_A$boundary$B1_SLP,
              cfg_A$boundary$B2_INT, cfg_A$boundary$B2_SLP)
) %>%
  dplyr::mutate(
    abs_diff = abs(diag_A - current),
    rel_diff = abs_diff / abs(current)
  )
message(paste(capture.output(print(parity, n = Inf)), collapse = "\n"))

if (max(parity$rel_diff) < 0.05) {
  message("Baseline reproduces current 6b values within 5% — methodology recovered.")
} else {
  message("WARNING: baseline differs from current 6b values by >5%.")
  message("  Possible causes: different filtering, weighted regression, or")
  message("  scaled-feature regression in the original. Inspect before")
  message("  drawing conclusions from configurations B and C.")
}

# ============================================================================
# WRITE OUTPUTS
# ============================================================================
coef_table <- tibble::tibble(
  configuration = c("A_baseline", "B_boundary_filter",
                    "C_gmm_and_boundary_filter", "current_6b_reference"),
  B1_INT = c(cfg_A$boundary$B1_INT, cfg_B$boundary$B1_INT,
             cfg_C$boundary$B1_INT, REF_B1_INT),
  B1_SLP = c(cfg_A$boundary$B1_SLP, cfg_B$boundary$B1_SLP,
             cfg_C$boundary$B1_SLP, REF_B1_SLP),
  B2_INT = c(cfg_A$boundary$B2_INT, cfg_B$boundary$B2_INT,
             cfg_C$boundary$B2_INT, REF_B2_INT),
  B2_SLP = c(cfg_A$boundary$B2_SLP, cfg_B$boundary$B2_SLP,
             cfg_C$boundary$B2_SLP, REF_B2_SLP),
  n_gmm  = c(nrow(peaks_full), nrow(peaks_full),
             nrow(peaks_interior), NA_integer_),
  n_bnd  = c(cfg_A$boundary$n_b1, cfg_B$boundary$n_b1,
             cfg_C$boundary$n_b1, NA_integer_)
)

readr::write_csv(coef_table, file.path(diag_dir, "coefficients_summary.csv"))
message(sprintf("\nWrote: %s", file.path(diag_dir, "coefficients_summary.csv")))

agreement_table <- dplyr::bind_rows(
  cfg_A$metrics %>% dplyr::mutate(configuration = "A_baseline"),
  cfg_B$metrics %>% dplyr::mutate(configuration = "B_boundary_filter"),
  cfg_C$metrics %>% dplyr::mutate(configuration = "C_gmm_and_boundary_filter")
) %>%
  dplyr::select(configuration, class, n_true, precision, recall)

readr::write_csv(agreement_table,
                 file.path(diag_dir, "classification_agreement.csv"))
message(sprintf("Wrote: %s",
                file.path(diag_dir, "classification_agreement.csv")))

# ============================================================================
# COMPARISON FIGURE
# ============================================================================
build_panel <- function(cfg, title_suffix) {
  bnd <- cfg$boundary
  rf_seq <- seq(0, 1, length.out = 400)

  b1_line <- tibble::tibble(
    rain_frac = rf_seq,
    swe_mm    = 10^(bnd$B1_INT + bnd$B1_SLP * rf_seq) - SWE_LOG_OFFSET
  ) %>% dplyr::filter(swe_mm > -SWE_LOG_OFFSET)

  sw_seq <- seq(min(peaks_full$log10_swe),
                max(peaks_full$log10_swe), length.out = 400)
  b2_line <- tibble::tibble(
    log10_swe = sw_seq,
    rain_frac = bnd$B2_INT + bnd$B2_SLP * sw_seq,
    swe_mm    = 10^sw_seq - SWE_LOG_OFFSET
  ) %>% dplyr::filter(rain_frac >= 0, rain_frac <= 1)

  ggplot(peaks_full,
         aes(x = rain_frac, y = swe_mean_mm + SWE_LOG_OFFSET)) +
    geom_point(alpha = 0.15, size = 0.5, color = "grey40") +
    geom_line(data = b1_line,
              aes(x = rain_frac, y = swe_mm + SWE_LOG_OFFSET),
              color = "#5B8DB8", linewidth = 0.9, inherit.aes = FALSE) +
    geom_line(data = b2_line,
              aes(x = rain_frac, y = swe_mm + SWE_LOG_OFFSET),
              color = "#E07B39", linewidth = 0.9, inherit.aes = FALSE) +
    scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
    scale_y_log10() +
    annotation_logticks(sides = "l") +
    labs(
      title    = title_suffix,
      subtitle = sprintf(
        "B1: log10(SWE+1) = %.3f + %.3f * rf\nB2: rf = %.3f + %.3f * log10(SWE+1)",
        bnd$B1_INT, bnd$B1_SLP, bnd$B2_INT, bnd$B2_SLP),
      x = "Rain fraction",
      y = "Mean antecedent SWE + 1 (mm)"
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.subtitle    = element_text(size = 8, color = "grey40"))
}

p_A <- build_panel(cfg_A, "A — Baseline (full data)")
p_B <- build_panel(cfg_B, "B — Boundary filter only")
p_C <- build_panel(cfg_C, "C — GMM + boundary filter")

fig <- (p_A | p_B | p_C) +
  plot_annotation(
    title    = "End-member sensitivity of GMM-derived B1/B2 boundaries",
    subtitle = sprintf(
      "End-member filter: rain_frac in (%.2f, %.2f); excludes %d of %d peaks (%.1f%%)",
      END_MEMBER_LO, END_MEMBER_HI,
      n_low + n_high, n_full, 100 * (n_low + n_high) / n_full),
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9,  color = "grey40")
    )
  )

ggsave(file.path(diag_dir, "boundary_comparison.png"),
       plot = fig, width = 13, height = 5, dpi = 150)
message(sprintf("Wrote: %s",
                file.path(diag_dir, "boundary_comparison.png")))

# ============================================================================
# SUMMARY TXT
# ============================================================================
sink(file.path(diag_dir, "gmm_summary.txt"))
cat("End-member sensitivity diagnostic\n")
cat(format(Sys.time()), "\n\n")
cat(sprintf("Peaks (full):     %d\n", n_full))
cat(sprintf("Peaks (interior): %d  (filter: rain_frac in (%.2f, %.2f))\n",
            nrow(peaks_interior), END_MEMBER_LO, END_MEMBER_HI))
cat(sprintf("Excluded by filter: %d (%.1f%%)\n\n",
            n_low + n_high, 100 * (n_low + n_high) / n_full))

cat("Coefficient summary:\n")
print(coef_table, n = Inf, width = Inf)

cat("\n\nClassification agreement vs. merged GMM labels:\n")
print(agreement_table, n = Inf, width = Inf)

cat("\n\nParity check (configuration A vs. current 6b values):\n")
print(parity, n = Inf, width = Inf)
sink()

message(sprintf("Wrote: %s", file.path(diag_dir, "gmm_summary.txt")))
message("\nDone.")