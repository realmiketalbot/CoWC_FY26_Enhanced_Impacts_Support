#!/usr/bin/env Rscript

# ============================================================================
# Script 7: Pre-compute FFA Results Across All Sites and Classification Methods
# ============================================================================
# Runs the mixture-model FFA (Yu et al. 2022) for every combination of:
#   - site   : all 206 Colorado GAGES-II sites
#   - method : depth_thresh, frac_thresh_3class, frac_thresh_2class
#   - family : gpd, lp3
#
# Six family-specific output files are written to _data/:
#
#   ffa_curves_{family}.rds
#     Dense RI grid (Map tab plot). One row per (site_no, method, curve_type,
#     ri). curve_type is one of: a flood type name, "mixture", "single_sample".
#
#   ffa_quantiles_{family}.rds
#     Sparse RI grid (Summary tab scatter). Same structure as ffa_curves but
#     evaluated at RI_GRID only (2, 5, 10, 25, 50, 100, 200, 500 years).
#     Column is `return_period` rather than `ri` to distinguish clearly.
#
#   ffa_params_{family}.rds
#     Fit parameters (Parameters table). One row per (site_no, method,
#     curve_type). GPD: mu, sigma, xi, n. LP3: mu_log, sigma_log, gamma, n.
#
# Plus two family-agnostic files:
#
#   ffa_classification_summary.rds   (unchanged from previous runs)
#   ffa_ros_temporal.rds             (unchanged from previous runs)
#
# Parallelisation: outer loop over sites x methods; inner loop over families
# is cheap and run sequentially within each worker. multicore plan is used
# on Linux/HPC; falls back to sequential if multicore unavailable.
#
# INPUTS:
#   _data/pot_peaks_classified_depth_thresh.rds  (Script 6a)
#   _data/pot_peaks_classified_frac_thresh_3class.rds (Script 6b)
#   _data/site_metadata.rds                   (Script 0)
#   functions.R
#
# OUTPUTS:
#   _data/ffa_curves_gpd.rds
#   _data/ffa_curves_lp3.rds
#   _data/ffa_quantiles_gpd.rds
#   _data/ffa_quantiles_lp3.rds
#   _data/ffa_params_gpd.rds
#   _data/ffa_params_lp3.rds
#   _data/ffa_classification_summary.rds
#   _data/ffa_ros_temporal.rds
#
# Author : Mike Talbot, PhD Student, Colorado State University
# Date   : 2026-03-17
# ============================================================================

library(tidyverse)
library(furrr)
library(lmomco)

source("functions.R")

# ============================================================================
# CONFIGURATION
# ============================================================================

N_WORKERS         <- min(parallelly::availableCores() - 1L, 32L)
MIN_TYPE_SAMPLE   <- 20L
MIN_SINGLE_SAMPLE <- 25L
FAMILIES          <- c("gpd", "lp3")

# Sparse grid for Summary tab quantile scatter
RI_GRID <- c(2, 5, 10, 25, 50, 100, 200, 500)

# Dense grid for Map tab curve rendering
RI_DENSE <- c(seq(1.1, 2, 0.1), seq(2.5, 10, 0.5), seq(11, 500, 1))

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading classified peak files...")

peaks_dav  <- readRDS("_data/pot_peaks_classified_depth_thresh.rds")
peaks_frac <- readRDS("_data/pot_peaks_classified_frac_thresh_3class.rds")
meta       <- readRDS("_data/site_metadata.rds") %>%
  select(site_no, regulated)

message(sprintf("  Depth Threshold  : %d peaks, %d sites",
                nrow(peaks_dav),  n_distinct(peaks_dav$site_no)))
message(sprintf("  Frac Thresh 3cls : %d peaks, %d sites",
                nrow(peaks_frac), n_distinct(peaks_frac$site_no)))

# ============================================================================
# METHOD REGISTRY
# ============================================================================

SNOW_INFLUENCED_TYPES <- c("rain_on_snow", "snowmelt")

remap_2class <- function(peaks) {
  # Remap using flood_type_raw (pre-demotion label) so that peaks individually
  # demoted to "unclassified" by MIN_TYPE_SAMPLE are still eligible for the
  # pooled snow_influenced type. Peaks unclassified due to data insufficiency
  # (flood_type_raw == "unclassified") remain unclassified.
  peaks %>%
    mutate(
      flood_type = case_when(
        as.character(flood_type_raw) %in% SNOW_INFLUENCED_TYPES ~ "snow_influenced",
        as.character(flood_type_raw) == "rainfall"           ~ "rainfall",
        TRUE                                                 ~ "unclassified"
      ),
      flood_type = factor(flood_type,
                          levels = c("rainfall", "snow_influenced", "unclassified"))
    )
}

methods <- list(

  depth_thresh = list(
    peaks        = peaks_dav,
    types        = c("rainfall", "snowmelt", "rain_on_snow"),
    remap        = NULL,
    pooled_types = NULL
  ),

  frac_thresh_3class = list(
    peaks        = peaks_frac,
    types        = c("rainfall", "snowmelt", "rain_on_snow"),
    remap        = NULL,
    pooled_types = NULL
  ),

  frac_thresh_2class = list(
    peaks        = peaks_frac,
    types        = c("rainfall", "snow_influenced"),
    remap        = remap_2class,
    pooled_types = "snow_influenced"
  )

)

all_sites <- sort(unique(c(
  unique(peaks_dav$site_no),
  unique(peaks_frac$site_no)
)))

message(sprintf("\n%d sites x %d methods x %d families = %d total fits\n",
                length(all_sites), length(methods), length(FAMILIES),
                length(all_sites) * length(methods) * length(FAMILIES)))

# ============================================================================
# HELPER: EXTRACT PARAMS FROM A FIT OBJECT
# ============================================================================

extract_params <- function(fit, nm, site, method_name, regulated, n) {
  if (fit$family == "gpd") {
    tibble(
      site_no    = site,
      method     = method_name,
      regulated  = regulated,
      curve_type = nm,
      n          = n,
      mu         = fit$mu,
      sigma      = fit$sigma,
      xi         = fit$xi,
      mu_log     = NA_real_,
      sigma_log  = NA_real_,
      gamma      = NA_real_
    )
  } else {
    tibble(
      site_no    = site,
      method     = method_name,
      regulated  = regulated,
      curve_type = nm,
      n          = n,
      mu         = NA_real_,
      sigma      = NA_real_,
      xi         = NA_real_,
      mu_log     = fit$mu_log,
      sigma_log  = fit$sigma_log,
      gamma      = fit$gamma
    )
  }
}

# ============================================================================
# HELPER: EVALUATE CURVES AT A GIVEN RI GRID
# ============================================================================

evaluate_curves <- function(fits, ri_grid, site, method_name, regulated) {
  if (length(fits) == 0L) return(NULL)
  type_fits <- fits[names(fits) != "single_sample"]

  # Per-type curves
  per_type <- purrr::imap_dfr(type_fits, function(f, typ) {
    q <- tryCatch(
      dist_curve(f, typ, "", ri_grid)$q,
      error = function(e) rep(NA_real_, length(ri_grid))
    )
    tibble(site_no = site, method = method_name, regulated = regulated,
           curve_type = typ, ri = ri_grid, discharge_cfs = pmax(0, q))
  })

  # Mixture
  mix <- if (length(type_fits) >= 2L) {
    q <- tryCatch(
      dist_mixture_curve(type_fits, ri = ri_grid)$q,
      error = function(e) rep(NA_real_, length(ri_grid))
    )
    tibble(site_no = site, method = method_name, regulated = regulated,
           curve_type = "mixture", ri = ri_grid, discharge_cfs = pmax(0, q))
  } else NULL

  # Single-sample
  ss <- if (!is.null(fits[["single_sample"]])) {
    q <- tryCatch(
      dist_curve(fits[["single_sample"]], "single_sample", "", ri_grid)$q,
      error = function(e) rep(NA_real_, length(ri_grid))
    )
    tibble(site_no = site, method = method_name, regulated = regulated,
           curve_type = "single_sample", ri = ri_grid, discharge_cfs = pmax(0, q))
  } else NULL

  bind_rows(per_type, mix, ss)
}

# ============================================================================
# PER-SITE WORKER FUNCTION
# ============================================================================
# Returns a named list:
#   class_rows  : one row per flood_type (family-agnostic)
#   ros_rows    : one row per water_year (family-agnostic)
#   gpd / lp3   : each a list with curve_rows, quantile_rows, param_rows

run_site_method <- function(site, method_name, method_cfg, regulated) {

  # ── Subset and remap ──────────────────────────────────────────────────────
  site_peaks <- method_cfg$peaks %>%
    filter(site_no == site, is.finite(peak_flow), peak_flow > 0)

  if (nrow(site_peaks) == 0L)
    return(list(class_rows = NULL, ros_rows = NULL, gpd = NULL, lp3 = NULL))

  # Remap relabels flood types only — in_type_sample carries through unchanged
  if (!is.null(method_cfg$remap))
    site_peaks <- method_cfg$remap(site_peaks)

  types        <- method_cfg$types
  pooled_types <- method_cfg$pooled_types
  mu_site      <- site_peaks$pot_threshold_cfs[1]

  # ── Classification summary (family-agnostic) ─────────────────────────────
  # Count peaks by flood_type_raw — the decision-tree label before any
  # MIN_TYPE_SAMPLE demotion. Reflects what each method classifies each peak
  # as, independent of whether that type has enough peaks for FFA fitting.
  # For frac_thresh_2class the remap is applied to flood_type_raw so that
  # snowmelt/rain_on_snow -> snow_influenced. Data-insufficient peaks have
  # flood_type_raw == "unclassified" and remain so.
  total_peaks <- nrow(site_peaks)

  # flood_type_raw is pre-remap; for frac_thresh_2class derive composition label
  # from the raw column directly rather than from the remapped flood_type
  composition_type <- if (!is.null(method_cfg$remap)) {
    dplyr::case_when(
      as.character(site_peaks$flood_type_raw) %in% SNOW_INFLUENCED_TYPES ~ "snow_influenced",
      as.character(site_peaks$flood_type_raw) == "rainfall"           ~ "rainfall",
      TRUE                                                            ~ "unclassified"
    )
  } else {
    as.character(site_peaks$flood_type_raw)
  }

  class_rows <- tibble(flood_type = composition_type) %>%
    dplyr::count(flood_type, name = "n_peaks") %>%
    mutate(
      site_no    = site,
      method     = method_name,
      regulated  = regulated,
      frac_peaks = n_peaks / total_peaks
    ) %>%
    select(site_no, method, regulated, flood_type, n_peaks, frac_peaks)

  # ── ROS temporal (family-agnostic) ───────────────────────────────────────
  ros_type <- switch(method_name,
    depth_thresh       = "rain_on_snow",
    frac_thresh_3class = "rain_on_snow",
    frac_thresh_2class = NULL,
    NULL
  )

  ros_rows <- if (!is.null(ros_type) && "water_year" %in% names(site_peaks)) {
    site_peaks %>%
      filter(in_type_sample | in_single_sample) %>%
      group_by(water_year) %>%
      summarise(
        n_peaks_total = n(),
        n_ros      = sum(as.character(flood_type) == ros_type,    na.rm = TRUE),
        n_rainfall = sum(as.character(flood_type) == "rainfall",  na.rm = TRUE),
        n_snowmelt = sum(as.character(flood_type) == "snowmelt",  na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(site_no = site, method = method_name, regulated = regulated) %>%
      select(site_no, method, regulated, water_year, everything())
  } else NULL

  # ── Build per-type fit samples ────────────────────────────────────────────
  # For pooled types, re-derive top-n from merged pool (one per water year).
  # For standard types, use the pre-computed in_type_sample flag.
  get_sample <- function(typ) {
    if (!is.null(pooled_types) && typ %in% pooled_types) {
      site_peaks %>%
        filter(flood_type == typ) %>%
        group_by(water_year) %>%
        slice_max(peak_flow, n = 1, with_ties = FALSE) %>%
        group_by() %>%
        mutate(n    = pmax(MIN_TYPE_SAMPLE, n_distinct(water_year)),
               rank = row_number(desc(peak_flow))) %>%
        filter(rank <= n) %>%
        ungroup()
    } else {
      site_peaks %>% filter(flood_type == typ, in_type_sample)
    }
  }

  peaks_single <- if ("in_single_sample" %in% names(site_peaks)) {
    site_peaks %>% filter(in_single_sample)
  } else site_peaks

  # ── Fit and evaluate for each family ─────────────────────────────────────
  family_results <- lapply(FAMILIES, function(family) {

    fits    <- list()
    n_peaks <- list()

    for (typ in types) {
      sub <- get_sample(typ)
      if (nrow(sub) < MIN_TYPE_SAMPLE) next
      f <- tryCatch(
        fit_dist(sub$peak_flow, sub$peak_date,
                 family = family, mu_threshold = mu_site,
                 min_n  = MIN_TYPE_SAMPLE),
        error = function(e) NULL
      )
      if (!is.null(f)) { fits[[typ]] <- f; n_peaks[[typ]] <- nrow(sub) }
    }

    f_ss <- tryCatch(
      fit_dist(peaks_single$peak_flow, peaks_single$peak_date,
               family = family, mu_threshold = mu_site,
               min_n  = MIN_SINGLE_SAMPLE),
      error = function(e) NULL
    )
    if (!is.null(f_ss)) {
      fits[["single_sample"]] <- f_ss
      n_peaks[["single_sample"]] <- nrow(peaks_single)
    }

    if (length(fits) == 0L) return(NULL)

    # Dense curves
    curve_rows <- evaluate_curves(fits, RI_DENSE, site, method_name, regulated)

    # Sparse quantiles (rename ri -> return_period for clarity)
    quantile_rows <- evaluate_curves(fits, RI_GRID, site, method_name, regulated) %>%
      rename(return_period = ri)

    # Parameters
    param_rows <- purrr::imap_dfr(fits, function(f, nm) {
      extract_params(f, nm, site, method_name, regulated, n_peaks[[nm]])
    })

    list(curve_rows = curve_rows, quantile_rows = quantile_rows, param_rows = param_rows)
  })

  names(family_results) <- FAMILIES

  c(list(class_rows = class_rows, ros_rows = ros_rows), family_results)
}

# ============================================================================
# PARALLELISED EXECUTION
# ============================================================================
message(sprintf("Starting parallel FFA with %d workers...", N_WORKERS))
plan(multicore, workers = N_WORKERS)

jobs <- expand.grid(
  site        = all_sites,
  method_name = names(methods),
  stringsAsFactors = FALSE
)

results <- future_pmap(
  list(jobs$site, jobs$method_name),
  function(site, method_name) {
    cfg       <- methods[[method_name]]
    reg_flag  <- meta$regulated[meta$site_no == site]
    regulated <- if (length(reg_flag) == 1L) reg_flag else NA
    run_site_method(site, method_name, cfg, regulated)
  },
  .progress = TRUE
)

plan(sequential)
message("Parallel run complete. Assembling outputs...")

# ============================================================================
# DIAGNOSTICS
# ============================================================================
message("\n=== PER-METHOD ROW COUNTS (pre-assembly) ===")
for (mn in names(methods)) {
  idx       <- which(jobs$method_name == mn)
  n_nonnull <- sum(!sapply(purrr::map(results[idx], "class_rows"), is.null))
  message(sprintf("  %-20s : %d / %d jobs returned non-NULL class_rows",
                  mn, n_nonnull, length(idx)))
}

# ============================================================================
# ASSEMBLE OUTPUTS
# ============================================================================
safe_bind <- function(lst, arrange_vars) {
  df <- bind_rows(lst)
  if (nrow(df) == 0L) return(df)
  arrange(df, across(all_of(arrange_vars)))
}

ffa_classification_summary <- safe_bind(
  purrr::map(results, "class_rows"),
  c("site_no", "method", "flood_type")
)

ffa_ros_temporal <- safe_bind(
  purrr::map(results, "ros_rows"),
  c("site_no", "method", "water_year")
)

for (fam in FAMILIES) {
  fam_results <- purrr::map(results, fam)

  curves    <- safe_bind(purrr::map(fam_results, "curve_rows"),
                         c("site_no", "method", "curve_type", "ri"))
  quantiles <- safe_bind(purrr::map(fam_results, "quantile_rows"),
                         c("site_no", "method", "curve_type", "return_period"))
  params    <- safe_bind(purrr::map(fam_results, "param_rows"),
                         c("site_no", "method", "curve_type"))

  saveRDS(curves,    sprintf("_data/ffa_curves_%s.rds",    fam))
  saveRDS(quantiles, sprintf("_data/ffa_quantiles_%s.rds", fam))
  saveRDS(params,    sprintf("_data/ffa_params_%s.rds",    fam))

  message(sprintf("\n%s curves    : %d rows | %d sites | %d curve types",
                  toupper(fam), nrow(curves),
                  n_distinct(curves$site_no),
                  n_distinct(curves$curve_type)))
  message(sprintf("%s quantiles : %d rows | %d return periods",
                  toupper(fam), nrow(quantiles),
                  n_distinct(quantiles$return_period)))
  message(sprintf("%s params    : %d rows",
                  toupper(fam), nrow(params)))
}

# ============================================================================
# FAMILY-AGNOSTIC OUTPUTS
# ============================================================================
saveRDS(ffa_classification_summary, "_data/ffa_classification_summary.rds")
saveRDS(ffa_ros_temporal,           "_data/ffa_ros_temporal.rds")

message(sprintf("\nffa_classification_summary : %d rows | %d sites | methods: %s",
                nrow(ffa_classification_summary),
                n_distinct(ffa_classification_summary$site_no),
                paste(sort(unique(ffa_classification_summary$method)), collapse = ", ")))

if (nrow(ffa_ros_temporal) > 0)
  message(sprintf("ffa_ros_temporal           : %d rows | water years %d-%d",
                  nrow(ffa_ros_temporal),
                  min(ffa_ros_temporal$water_year, na.rm = TRUE),
                  max(ffa_ros_temporal$water_year, na.rm = TRUE)))

message("\nScript 7 complete.")