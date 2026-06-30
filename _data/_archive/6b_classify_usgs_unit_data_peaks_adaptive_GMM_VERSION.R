#!/usr/bin/env Rscript

# ============================================================================
# Script 6b: Flood Type Classification — Adaptive Window
# ============================================================================
# Classifies each POT peak into one of three flood types using a linear SVM
# decision tree applied to WLDAS-derived rain fraction and mean antecedent SWE,
# computed over the basin-adaptive antecedent window (search_radius):
#
#   rain_frac   : weighted rain fraction over search_radius days
#                 = total_rain_flux / (total_rain_flux + total_melt_flux)
#                 on days where either flux exceeds NOISE_FLOOR_MM
#   swe_mean_mm : mean basin SWE (mm) over search_radius days (state variable)
#
# Decision tree (evaluated in order):
#   1. Unclassified : n_days_data < completeness threshold (WLDAS coverage only)
#   2. Snowmelt     : log10(swe_mean_mm + 1) > B1_INTERCEPT + B1_SLOPE * rain_frac
#   3. ROS          : rain_frac < B2_INTERCEPT + B2_SLOPE * log10(swe_mean_mm + 1)
#   4. Rainfall     : all remaining data-sufficient peaks
#
# Boundaries are derived at run time from this pipeline's own feature data:
#   1. Fit a G=4 EEV Gaussian mixture model to (rain_frac, log10_swe) for all
#      data-sufficient peaks with valid rain_frac and swe_mean_mm.
#   2. Merge the four components into three classes by ordering on median
#      rain_frac (ascending: snowmelt, ROS, rainfall [components 3+4 merged]).
#   3. Exclude end-member peaks (rain_frac <= ROS_RAIN_FRAC_MIN or
#      >= ROS_RAIN_FRAC_MAX) from SVM fitting; the well-separated end-member
#      populations would otherwise dominate the margin.
#   4. Fit two linear SVMs on the interior subset (e1071, scale=FALSE):
#        B1: snowmelt vs. ROS      → log10_swe = B1_INTERCEPT + B1_SLOPE * rain_frac
#        B2: ROS vs. rainfall      → rain_frac = B2_INTERCEPT + B2_SLOPE * log10_swe
#   5. Recover slope/intercept by regressing SVM decision values on raw features
#      (see svm_canonical/canonical_boundary_script.R for derivation).
#
# This replaces the previous scalar threshold scheme (RAIN_EM_THRESH = 0.95,
# SNOW_EM_THRESH = 0.10) which used rain_frac alone and ignored SWE entirely.
# The 2D boundaries recover the same process-population structure and are now
# derived directly from the data rather than hardcoded.
#
# INPUTS:
#   - _data/accum_curves_60d.rds     : daily rain/melt/SWE per peak (script 5)
#   - _data/pot_peaks_all.rds        : POT peaks with search_radius per peak
#   - _data/site_metadata.rds        : regulated flag, basin attributes
#
# OUTPUT:
#   - _data/pot_peaks_classified_adaptive.rds : peaks with flood_type and
#                                               classification diagnostic columns
# ============================================================================

library(tidyverse)
library(mclust)
library(e1071)

# ============================================================================
# CONFIGURATION
# ============================================================================
accum_file  <- "_data/accum_curves_60d.rds"
peaks_file  <- "_data/pot_peaks_all.rds"
meta_file   <- "_data/site_metadata.rds"
output_file <- "_data/pot_peaks_classified_adaptive.rds"

# GMM + SVM boundary derivation parameters
# Boundaries are derived at run-time (see DERIVE CLASSIFICATION BOUNDARIES below)
GMM_SEED <- 42L  # reproducibility seed for Mclust
B1_COST  <- 1    # SVM cost, snowmelt vs. ROS boundary
B2_COST  <- 1    # SVM cost, ROS vs. rainfall boundary

# Hard guards on ROS eligibility: prevent boundary artifacts where near-zero
# or near-unity rain_frac values are misclassified as ROS solely because the
# B2 intercept shifts with log10(SWE+1). A peak with rain_frac < ROS_RAIN_FRAC_MIN
# is effectively pure snowmelt and should not reach the ROS branch; one with
# rain_frac > ROS_RAIN_FRAC_MAX is effectively pure rainfall.
ROS_RAIN_FRAC_MIN <- 0.01   # rain_frac below this → cannot be ROS (effectively pure snowmelt)
ROS_RAIN_FRAC_MAX <- 0.99   # rain_frac above this → cannot be ROS (effectively pure rainfall)

# SWE log offset (must match value used during GMM/logistic fitting)
SWE_LOG_OFFSET <- 1.0   # mm

# Noise floor: days below this flux are treated as inactive
NOISE_FLOOR_MM    <- 0.1    # mm/day

# Completeness: minimum active days required for classification
MIN_DAYS_SHORT    <- 0.90   # fraction of search_radius days required
# MIN_TOTAL_FLUX_MM removed: threshold sensitivity sweep showed FFA quantiles
# insensitive across 0.1–10 mm; runoff–discharge mismatch analysis found no
# viable data-driven replacement. All peaks with complete WLDAS coverage are
# now classified regardless of total flux magnitude.

# Yu et al. (2022) top-n sample size thresholds
MIN_SINGLE_SAMPLE <- 25     # minimum peaks for single-sample fit
MIN_TYPE_SAMPLE   <- 20     # minimum peaks per type for per-type fit

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading data...")
accum <- readRDS(accum_file)
peaks <- readRDS(peaks_file) %>%
  mutate(peak_date = as.Date(peak_date))
meta  <- readRDS(meta_file)

# pot_peaks_all.rds is one row per (peak, window day) — deduplicate to one
# row per peak for use as the classification lookup table.
peaks <- peaks %>%
  distinct(site_no, peak_id, peak_date, peak_flow,
           pot_threshold_cfs, search_radius,
           antecedent_runoff_mm, antecedent_days_used)

message(sprintf("Peaks: %d across %d sites",
                nrow(peaks), n_distinct(peaks$site_no)))
message(sprintf("Accum curve rows: %d", nrow(accum)))

# ============================================================================
# COMPUTE FEATURES
# ============================================================================
message("\nComputing classification features...")

# Helper: weighted rain fraction and total fluxes over a given max lookback
# Returns one row per (site_no, peak_id) with summary metrics
compute_features <- function(accum_df, peaks_df, max_lookback_col) {
  peaks_windows <- peaks_df %>%
    select(site_no, peak_id, window_days = {{ max_lookback_col }})

  accum_df %>%
    inner_join(peaks_windows, by = c("site_no", "peak_id")) %>%
    filter(lookback_day <= window_days) %>%
    mutate(
      rain_active  = if_else(rain_mm_day >= NOISE_FLOOR_MM, rain_mm_day, 0),
      melt_active  = if_else(melt_mm_day >= NOISE_FLOOR_MM, melt_mm_day, 0),
      total_active = rain_active + melt_active
    ) %>%
    group_by(site_no, peak_id) %>%
    summarise(
      total_rain_flux = sum(rain_active,  na.rm = TRUE),
      total_melt_flux = sum(melt_active,  na.rm = TRUE),
      total_flux      = sum(total_active, na.rm = TRUE),
      n_days_data     = first(n_days_data),
      .groups = "drop"
    ) %>%
    mutate(
      weighted_rain_frac = if_else(
        total_flux > 0,
        total_rain_flux / total_flux,
        NA_real_
      )
    )
}

# Short window features (search_radius)
peaks_with_windows <- peaks %>%
  mutate(short_window = search_radius)

message("  Computing rain_frac (search_radius window)...")
feat_short <- compute_features(accum, peaks_with_windows, short_window) %>%
  rename(
    rain_frac         = weighted_rain_frac,
    total_flux_short  = total_flux,
    rain_flux_short   = total_rain_flux,
    melt_flux_short   = total_melt_flux
  )

# Mean antecedent SWE over short window (state variable, not flux)
# swe_mm is basin-mean SWE in mm; NA days are excluded from the mean.
message("  Computing mean antecedent SWE (search_radius window)...")
swe_feat <- accum %>%
  inner_join(peaks_with_windows %>% select(site_no, peak_id, short_window),
             by = c("site_no", "peak_id")) %>%
  filter(lookback_day <= short_window) %>%
  group_by(site_no, peak_id) %>%
  summarise(
    swe_mean_mm = mean(swe_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(swe_mean_mm = if_else(is.nan(swe_mean_mm), NA_real_, swe_mean_mm))

# ============================================================================
# DERIVE CLASSIFICATION BOUNDARIES
# ============================================================================
message("\nDeriving classification boundaries (GMM + linear SVM)...")

# Temporary feature table: one row per data-sufficient peak with both
# rain_frac and log10_swe defined — mirrors what canonical_boundary_script.R
# reads from the already-classified output file.
gmm_input <- peaks_with_windows %>%
  dplyr::left_join(feat_short, by = c("site_no", "peak_id")) %>%
  dplyr::left_join(swe_feat,   by = c("site_no", "peak_id")) %>%
  dplyr::mutate(
    data_sufficient = n_days_data >= ceiling(short_window * MIN_DAYS_SHORT),
    log10_swe       = log10(dplyr::coalesce(swe_mean_mm, 0) + SWE_LOG_OFFSET)
  ) %>%
  dplyr::filter(data_sufficient, !is.na(rain_frac), !is.na(log10_swe)) %>%
  dplyr::select(site_no, peak_id, rain_frac, log10_swe)

message(sprintf("  Peaks for GMM: %d", nrow(gmm_input)))

# G=4 EEV GMM on unscaled (rain_frac, log10_swe)
set.seed(GMM_SEED)
mc_fit <- mclust::Mclust(
  data       = cbind(gmm_input$rain_frac, gmm_input$log10_swe),
  G          = 4L,
  modelNames = "EEV",
  verbose    = FALSE
)

gmm_input <- gmm_input %>%
  dplyr::mutate(component = as.integer(mc_fit$classification))

comp_order <- gmm_input %>%
  dplyr::group_by(component) %>%
  dplyr::summarise(rain_frac_med = median(rain_frac), .groups = "drop") %>%
  dplyr::arrange(rain_frac_med)

comp_snowmelt  <- comp_order$component[1]
comp_ros       <- comp_order$component[2]
comp_rain_low  <- comp_order$component[3]
comp_rain_high <- comp_order$component[4]

gmm_input <- gmm_input %>%
  dplyr::mutate(
    gmm_class = dplyr::case_when(
      component %in% c(comp_rain_low, comp_rain_high) ~ "rainfall",
      component == comp_ros                           ~ "ROS",
      component == comp_snowmelt                      ~ "snowmelt"
    ),
    gmm_class = factor(gmm_class, levels = c("snowmelt", "ROS", "rainfall"))
  )

# Interior subset: exclude end-members (same thresholds as ROS guards above)
gmm_interior <- gmm_input %>%
  dplyr::filter(rain_frac > ROS_RAIN_FRAC_MIN,
                rain_frac < ROS_RAIN_FRAC_MAX)

message(sprintf("  Interior subset (%.2f < rain_frac < %.2f): %d of %d peaks",
                ROS_RAIN_FRAC_MIN, ROS_RAIN_FRAC_MAX,
                nrow(gmm_interior), nrow(gmm_input)))

# Linear SVM: extract slope/intercept via linear regression on decision values.
# class_a is placed on the positive side of the hyperplane.
fit_svm_boundary <- function(df, class_a, class_b, cost) {
  d <- df %>%
    dplyr::filter(gmm_class %in% c(class_a, class_b)) %>%
    dplyr::mutate(cls = factor(gmm_class, levels = c(class_a, class_b)))
  model    <- e1071::svm(cls ~ rain_frac + log10_swe, data = d,
                         kernel = "linear", cost = cost, scale = FALSE,
                         decision.values = TRUE)
  pred_obj <- predict(model, d, decision.values = TRUE)
  dv       <- as.numeric(attr(pred_obj, "decision.values"))
  lm_fit   <- lm(dv ~ rain_frac + log10_swe, data = d)
  w        <- coef(lm_fit)[c("rain_frac", "log10_swe")]
  b        <- -coef(lm_fit)["(Intercept)"]
  names(w) <- c("rain_frac", "log10_swe")
  mean_a   <- colMeans(d[d$cls == class_a, c("rain_frac", "log10_swe")])
  if (sum(w * mean_a) - b < 0) { w <- -w; b <- -b }
  list(w = w, b = unname(b))
}

b1_svm <- fit_svm_boundary(gmm_interior, "snowmelt", "ROS",      B1_COST)
b2_svm <- fit_svm_boundary(gmm_interior, "ROS",      "rainfall", B2_COST)

# Convert w/b form to slope/intercept form used by the classification below:
#   B1: log10_swe = B1_INTERCEPT + B1_SLOPE * rain_frac
#   B2: rain_frac = B2_INTERCEPT + B2_SLOPE * log10_swe
B1_INTERCEPT <- unname(b1_svm$b / b1_svm$w["log10_swe"])
B1_SLOPE     <- unname(-b1_svm$w["rain_frac"] / b1_svm$w["log10_swe"])
B2_INTERCEPT <- unname(b2_svm$b / b2_svm$w["rain_frac"])
B2_SLOPE     <- unname(-b2_svm$w["log10_swe"] / b2_svm$w["rain_frac"])

message(sprintf("  B1 (snowmelt/ROS):  log10_swe = %.4f + %.4f * rain_frac",
                B1_INTERCEPT, B1_SLOPE))
message(sprintf("  B2 (ROS/rainfall):  rain_frac = %.4f + %.4f * log10_swe",
                B2_INTERCEPT, B2_SLOPE))

rm(gmm_input, comp_order, comp_snowmelt, comp_ros, comp_rain_low, comp_rain_high,
   gmm_interior, b1_svm, b2_svm, mc_fit)

# ============================================================================
# JOIN FEATURES AND CLASSIFY
# ============================================================================
message("\nClassifying peaks...")

classified <- peaks %>%
  left_join(peaks_with_windows %>% select(site_no, peak_id, short_window),
            by = c("site_no", "peak_id")) %>%
  left_join(feat_short, by = c("site_no", "peak_id")) %>%
  left_join(swe_feat,   by = c("site_no", "peak_id")) %>%
  left_join(meta %>% select(site_no, regulated), by = "site_no") %>%
  mutate(
    # Completeness flag: require MIN_DAYS_SHORT fraction of short window.
    # Flux floor removed — diagnostic analysis found no empirical support for
    # a magnitude threshold (see diagnostic_script_for_unclassified_peaks_all_sites.R).
    data_sufficient = n_days_data >= ceiling(short_window * MIN_DAYS_SHORT),

    # Precompute log10(SWE + 1) for boundary evaluation
    log10_swe = log10(coalesce(swe_mean_mm, 0) + SWE_LOG_OFFSET),

    # Sequential SVM decision tree (B1/B2 derived above):
    #   1. Unclassified : incomplete window OR zero active flux (rain_frac NA)
    #   2. Snowmelt     : log10(SWE+1) above boundary 1
    #   3. ROS          : rain_frac below boundary 2
    #   4. Rainfall     : all remaining data-sufficient peaks
    flood_type = case_when(
      !data_sufficient | is.na(rain_frac)                          ~ "unclassified",
      log10_swe > B1_INTERCEPT + B1_SLOPE * rain_frac             ~ "snowmelt",
      rain_frac < ROS_RAIN_FRAC_MIN                               ~ "snowmelt",
      rain_frac > ROS_RAIN_FRAC_MAX                               ~ "rainfall",
      rain_frac < B2_INTERCEPT + B2_SLOPE * log10_swe             ~ "ROS",
      TRUE                                                         ~ "rainfall"
    ),

    # Preserve pre-sample-size-demotion label for downstream diagnostics
    flood_type_raw = flood_type,
    flood_type     = factor(flood_type,
                            levels = c("rainfall", "snowmelt",
                                       "ROS", "unclassified"))
  )

# ============================================================================
# PRE-SELECTION SUMMARY
# ============================================================================
message("\n=== CLASSIFICATION SUMMARY (pre-selection) ===")
type_counts <- classified %>%
  dplyr::count(flood_type) %>%
  mutate(pct = n / sum(n) * 100)

walk2(type_counts$flood_type, type_counts$n, function(t, n) {
  pct <- type_counts$pct[type_counts$flood_type == t]
  message(sprintf("  %-14s : %5d  (%5.1f%%)", t, n, pct))
})

message(sprintf("\nTotal peaks: %d | Sites: %d",
                nrow(classified), n_distinct(classified$site_no)))

reg_summary <- classified %>%
  dplyr::count(regulated, flood_type) %>%
  mutate(pct = n / sum(n) * 100)

message("\n=== BY REGULATION STATUS (pre-selection) ===")
message(paste(capture.output(print(reg_summary)), collapse = "\n"))

# ============================================================================
# YU ET AL. (2022) TOP-N PEAK SELECTION
# ============================================================================
message("\nApplying Yu et al. (2022) top-n peak selection...")

water_year <- function(d) {
  yr <- as.integer(format(d, "%Y"))
  mo <- as.integer(format(d, "%m"))
  ifelse(mo >= 10L, yr + 1L, yr)
}

classified <- classified %>%
  mutate(water_year = water_year(peak_date))

# Single-sample: largest m peaks, at most one per water year per site
single_sample_ids <- classified %>%
  filter(flood_type != "unclassified") %>%
  group_by(site_no, water_year) %>%
  slice_max(peak_flow, n = 1, with_ties = FALSE) %>%
  group_by(site_no) %>%
  mutate(m    = pmax(MIN_SINGLE_SAMPLE, n_distinct(water_year)),
         rank = row_number(desc(peak_flow))) %>%
  filter(rank <= m) %>%
  ungroup() %>%
  select(site_no, peak_id) %>%
  mutate(in_single_sample = TRUE)

# Per-type: largest n peaks, at most one per water year per site-type;
# types with fewer than MIN_TYPE_SAMPLE events are dropped entirely
type_sample_ids <- classified %>%
  filter(flood_type != "unclassified") %>%
  group_by(site_no, flood_type, water_year) %>%
  slice_max(peak_flow, n = 1, with_ties = FALSE) %>%
  group_by(site_no, flood_type) %>%
  mutate(n    = pmax(MIN_TYPE_SAMPLE, n_distinct(water_year)),
         rank = row_number(desc(peak_flow))) %>%
  filter(rank <= n) %>%
  filter(n() >= MIN_TYPE_SAMPLE) %>%
  ungroup() %>%
  select(site_no, peak_id, flood_type) %>%
  mutate(in_type_sample = TRUE)

classified <- classified %>%
  left_join(single_sample_ids, by = c("site_no", "peak_id")) %>%
  left_join(type_sample_ids %>% select(site_no, peak_id, in_type_sample),
            by = c("site_no", "peak_id")) %>%
  mutate(
    in_single_sample = replace_na(in_single_sample, FALSE),
    in_type_sample   = replace_na(in_type_sample,   FALSE)
  )

# Demote types that did not meet MIN_TYPE_SAMPLE to "unclassified"
retained_site_types <- type_sample_ids %>% distinct(site_no, flood_type)

classified <- classified %>%
  mutate(flood_type = if_else(
    flood_type != "unclassified" &
      !paste(site_no, flood_type) %in%
        paste(retained_site_types$site_no, retained_site_types$flood_type),
    "unclassified",
    as.character(flood_type)
  )) %>%
  mutate(flood_type = factor(flood_type,
                             levels = c("rainfall", "snowmelt", "ROS",
                                        "unclassified")))

# ============================================================================
# POST-SELECTION SUMMARY
# ============================================================================
message("\n=== CLASSIFICATION SUMMARY (post-selection) ===")
type_counts_post <- classified %>%
  dplyr::count(flood_type) %>%
  mutate(pct = n / sum(n) * 100)

walk2(type_counts_post$flood_type, type_counts_post$n, function(t, n) {
  pct <- type_counts_post$pct[type_counts_post$flood_type == t]
  message(sprintf("  %-14s : %5d  (%5.1f%%)", t, n, pct))
})

selection_summary <- classified %>%
  group_by(site_no) %>%
  summarise(
    n_single_sample = sum(in_single_sample),
    n_rainfall      = sum(flood_type == "rainfall" & in_type_sample),
    n_snowmelt      = sum(flood_type == "snowmelt" & in_type_sample),
    n_ROS           = sum(flood_type == "ROS"      & in_type_sample),
    .groups = "drop"
  )

message(sprintf("\nSites with >= %d single-sample peaks : %d",
                MIN_SINGLE_SAMPLE,
                sum(selection_summary$n_single_sample >= MIN_SINGLE_SAMPLE)))
message(sprintf("Sites with >= %d rainfall peaks      : %d",
                MIN_TYPE_SAMPLE, sum(selection_summary$n_rainfall >= MIN_TYPE_SAMPLE)))
message(sprintf("Sites with >= %d snowmelt peaks      : %d",
                MIN_TYPE_SAMPLE, sum(selection_summary$n_snowmelt >= MIN_TYPE_SAMPLE)))
message(sprintf("Sites with >= %d ROS peaks           : %d",
                MIN_TYPE_SAMPLE, sum(selection_summary$n_ROS     >= MIN_TYPE_SAMPLE)))

# ============================================================================
# SAVE
# ============================================================================
saveRDS(classified, output_file)
message(sprintf("\nSaved %d classified peaks to %s", nrow(classified), output_file))
message(sprintf("  MIN_TYPE_SAMPLE  : %d", MIN_TYPE_SAMPLE))
message(sprintf("  Columns: %s", paste(names(classified), collapse = ", ")))