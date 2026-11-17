#!/usr/bin/env Rscript

# ============================================================================
# Script 7: Poisson–GPD flood frequency analysis across all sites and
#           classification methods (milestone M1)
# ============================================================================
# For every site x method (depth_thresh, frac_thresh_3class,
# frac_thresh_2class), fits:
#   - a Poisson–GPD per flood type (Proposal Eq. 1.3–1.4),
#   - their mixture, F_mix = prod_i F_i (Eq. 1.5), when >= 2 types are fitted,
#   - a single-sample Poisson–GPD to all eligible peaks regardless of type.
# Return periods are ANNUAL (T = 1 / (1 - F)), not PoT recurrence intervals.
#
# M1 SCOPE (see _planning/agu2026_interim_plan.md):
#   - Fitting threshold u: INTERIM rate-based rule (D-A), a stand-in for the
#     EQD selection of Murphy et al. (2025). select_threshold_rate() keeps the
#     k = min(n, max(MIN_TYPE_COUNT, n_years)) largest peaks, the closest
#     analogue to the Yu et al. (2022) top-n sample used in the initial
#     results. Recorded in the u_rule column so EQD can replace it.
#   - Fit pool: peaks with fit_eligible == TRUE (Scripts 6a/6b) in the site's
#     QUALIFYING water years (Script 2b). lambda = exceedances / number of
#     qualifying years (D-E).
#   - Single sample: all peaks in the fit pool (D3).
#   - Below-threshold treatment: F_i floored at exp(-lambda_i); quantiles
#     reported only where valid, otherwise NA (D4).
#   - Uncertainty: lo/hi columns exist but are NA until M3 (double bootstrap).
#   - LP3 is retired (D5).
#
# INPUTS:
#   _data/pot_peaks_classified_depth_thresh.rds        (Script 6a)
#   _data/pot_peaks_classified_frac_thresh_3class.rds  (Script 6b)
#   _data/site_qualifying_years.rds                    (Script 2b)
#   _data/site_metadata.rds                            (Script 0)
#   functions.R
#
# OUTPUTS:
#   _data/ffa_curves_pgpd.rds     dense RI grid; one row per
#                                 (site_no, method, curve_type, ri)
#   _data/ffa_quantiles_pgpd.rds  RI_GRID only; column return_period
#   _data/ffa_params_pgpd.rds     one row per (site_no, method, curve_type):
#                                 u, u_rule, k, n_exc, n_years, lambda, sigma, xi
#   _data/ffa_classification_summary.rds  peak composition by method
#                                 (schema unchanged; read by the App and
#                                 report_figures.R)
#   curve_type is a flood type name, "mixture" or "single_sample".
#
# The previous outputs (ffa_*_gpd.rds, ffa_*_lp3.rds, ffa_ros_temporal.rds)
# are no longer produced; they remain at tag initial-results-prelim-proposal.
# ============================================================================

library(tidyverse)
library(furrr)

source("functions.R")

# ============================================================================
# CONFIGURATION
# ============================================================================
N_WORKERS      <- max(1L, min(parallelly::availableCores() - 1L, 32L))
MIN_TYPE_COUNT <- 20L     # minimum exceedances for any fit (D-B)
THRESH_RATE    <- 1       # interim rule: about one exceedance per qualifying year

# Sensitivity runs (_documentation/poster_sensitivity.R) override these
# through environment variables; with none set, the defaults below apply.
THRESH_RATE     <- as.numeric(Sys.getenv("FFA_THRESH_RATE", THRESH_RATE))
PEAKS_FRAC_FILE <- Sys.getenv("FFA_PEAKS_FRAC_FILE", "_data/pot_peaks_classified_frac_thresh_3class.rds")
OUT_DIR         <- Sys.getenv("FFA_OUT_DIR", "_data")

RI_GRID  <- c(2, 5, 10, 25, 50, 100, 200, 500)
RI_DENSE <- c(seq(1.1, 2, 0.1), seq(2.5, 10, 0.5), seq(11, 500, 1))

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading inputs...")
peaks_dav  <- readRDS("_data/pot_peaks_classified_depth_thresh.rds")
peaks_frac <- readRDS(PEAKS_FRAC_FILE)
meta       <- readRDS("_data/site_metadata.rds") %>% select(site_no, regulated)
qual       <- readRDS("_data/site_qualifying_years.rds") %>% filter(selected, qualifies)

qual_years <- split(qual$water_year, qual$site_no)

methods <- list(
  depth_thresh       = list(peaks = peaks_dav,  types = c("rainfall", "snowmelt", "rain_on_snow"), remap = NULL),
  frac_thresh_3class = list(peaks = peaks_frac, types = c("rainfall", "snowmelt", "rain_on_snow"), remap = NULL),
  frac_thresh_2class = list(peaks = peaks_frac, types = c("rainfall", "snow_influenced"),           remap = remap_2class)
)

all_sites <- sort(intersect(unique(c(peaks_dav$site_no, peaks_frac$site_no)), names(qual_years)))
message(sprintf("%d sites x %d methods", length(all_sites), length(methods)))

# ============================================================================
# PER SITE x METHOD
# ============================================================================
run_site_method <- function(site, method_name, cfg, regulated) {

  site_peaks <- cfg$peaks %>% filter(site_no == site, is.finite(peak_flow), peak_flow > 0)
  if (nrow(site_peaks) == 0L) return(NULL)
  if (!is.null(cfg$remap)) site_peaks <- cfg$remap(site_peaks)
  label_col <- if (is.null(cfg$remap)) "flood_type_raw" else "flood_type"

  # Peak composition over all peaks (unchanged from the initial results)
  class_rows <- tibble(flood_type = as.character(site_peaks[[label_col]])) %>%
    count(flood_type, name = "n_peaks") %>%
    mutate(site_no = site, method = method_name, regulated = regulated,
           frac_peaks = n_peaks / sum(n_peaks)) %>%
    select(site_no, method, regulated, flood_type, n_peaks, frac_peaks)

  # Fit pool: eligible peaks in qualifying water years
  wys     <- qual_years[[site]]
  n_years <- length(wys)
  pool    <- site_peaks %>% filter(fit_eligible, water_year %in% wys)
  floor_u <- site_peaks$candidate_threshold_cfs[1]

  fit_one <- function(x) {
    sel <- select_threshold_rate(x, n_years, floor_u, MIN_TYPE_COUNT, THRESH_RATE)
    if (is.null(sel)) return(NULL)
    f <- fit_pgpd(x, sel$u, n_years, min_n = MIN_TYPE_COUNT)
    if (is.null(f)) return(NULL)
    c(f, list(k = sel$k, u_rule = sel$rule))
  }

  fits <- list()
  for (typ in cfg$types) {
    f <- fit_one(pool$peak_flow[as.character(pool[[label_col]]) == typ])
    if (!is.null(f)) fits[[typ]] <- f
  }
  type_fits <- fits
  f_ss <- fit_one(pool$peak_flow)
  if (!is.null(f_ss)) fits[["single_sample"]] <- f_ss

  if (length(fits) == 0L)
    return(list(class_rows = class_rows, curves = NULL, quantiles = NULL, params = NULL))

  curves_at <- function(ri) {
    out <- purrr::imap_dfr(fits, function(f, nm)
      tibble(curve_type = nm, ri = ri, discharge_cfs = pgpd_quantile(ri, f)))
    if (length(type_fits) >= 2L)
      out <- bind_rows(out, tibble(curve_type = "mixture", ri = ri,
                                   discharge_cfs = pgpd_mixture_quantile(ri, type_fits)))
    out %>% mutate(site_no = site, method = method_name, regulated = regulated,
                   lo = NA_real_, hi = NA_real_) %>%
      select(site_no, method, regulated, curve_type, ri, discharge_cfs, lo, hi)
  }

  params <- purrr::imap_dfr(fits, function(f, nm) tibble(
    site_no = site, method = method_name, regulated = regulated, curve_type = nm,
    u = f$u, u_rule = f$u_rule, k = f$k, n_exc = f$n_exc, n_years = f$n_years,
    lambda = f$lambda, sigma = f$sigma, xi = f$xi))

  list(class_rows = class_rows,
       curves     = curves_at(RI_DENSE),
       quantiles  = curves_at(RI_GRID) %>% rename(return_period = ri),
       params     = params)
}

# ============================================================================
# RUN
# ============================================================================
jobs <- expand.grid(site = all_sites, method_name = names(methods), stringsAsFactors = FALSE)
reg_lookup <- setNames(meta$regulated, meta$site_no)

message(sprintf("Fitting with %d worker(s)...", N_WORKERS))
if (N_WORKERS > 1) plan(multisession, workers = N_WORKERS) else plan(sequential)
results <- future_pmap(list(jobs$site, jobs$method_name), function(site, method_name)
  run_site_method(site, method_name, methods[[method_name]], unname(reg_lookup[site])),
  .progress = TRUE, .options = furrr_options(seed = NULL))
plan(sequential)

gather <- function(el, keys) {
  d <- bind_rows(purrr::map(results, el))
  if (nrow(d) == 0) d else arrange(d, across(all_of(keys)))
}
curves    <- gather("curves",     c("site_no", "method", "curve_type", "ri"))
quantiles <- gather("quantiles",  c("site_no", "method", "curve_type", "return_period"))
params    <- gather("params",     c("site_no", "method", "curve_type"))
class_sum <- gather("class_rows", c("site_no", "method", "flood_type"))

run_config <- list(min_type_count = MIN_TYPE_COUNT, thresh_rate = THRESH_RATE,
                   peaks_frac_file = PEAKS_FRAC_FILE,
                   threshold_rule = "select_threshold_rate (interim, D-A)",
                   qualifying_years = attr(readRDS("_data/site_qualifying_years.rds"), "selection_config"))
for (nm in c("curves", "quantiles", "params")) {
  obj <- get(nm); attr(obj, "run_config") <- run_config; assign(nm, obj)
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
saveRDS(curves,    file.path(OUT_DIR, "ffa_curves_pgpd.rds"))
saveRDS(quantiles, file.path(OUT_DIR, "ffa_quantiles_pgpd.rds"))
saveRDS(params,    file.path(OUT_DIR, "ffa_params_pgpd.rds"))
saveRDS(class_sum, file.path(OUT_DIR, "ffa_classification_summary.rds"))

# ============================================================================
# SUMMARY
# ============================================================================
message("\n=== Fitted curves per method ===")
summ <- quantiles %>%
  group_by(method) %>%
  summarise(sites_single  = n_distinct(site_no[curve_type == "single_sample"]),
            sites_mixture = n_distinct(site_no[curve_type == "mixture"]),
            sites_both    = n_distinct(intersect(site_no[curve_type == "single_sample"],
                                                 site_no[curve_type == "mixture"])),
            .groups = "drop")
message(paste(capture.output(print(summ)), collapse = "\n"))
message(sprintf("\nShape xi (single sample): median %.3f, IQR %.3f to %.3f",
                median(params$xi[params$curve_type == "single_sample"]),
                quantile(params$xi[params$curve_type == "single_sample"], .25),
                quantile(params$xi[params$curve_type == "single_sample"], .75)))
message(sprintf("Saved ffa_curves_pgpd (%d rows), ffa_quantiles_pgpd (%d), ffa_params_pgpd (%d), ffa_classification_summary (%d)",
                nrow(curves), nrow(quantiles), nrow(params), nrow(class_sum)))
message("\nScript 7 complete.")
