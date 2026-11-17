#!/usr/bin/env Rscript
# ============================================================================
# Poster sensitivity runs (plan S3 and S4; AGU 2026)
# ============================================================================
# Reruns Script 7 under alternative settings and summarizes how the
# mixture / single-sample ratio responds, by annual return period.
#
#   S3  interim threshold rate (exceedances per qualifying year):
#       0.5, 2, 3 (baseline 1). At 0.5 the 20-exceedance minimum binds at
#       most sites, so it is effectively "the 20 largest peaks".
#   S4  rain-fraction classification:
#       alpha_r (rainfall cutoff) 0.90, 0.95 (baseline 0.99), with alpha_s
#       recomputed as the median of the mixed events (the proposal's rule);
#       alpha_s fixed at 0.20, 0.40 (baseline: network median, ~0.31).
#       Peaks are relabeled from the stored rain_frac, using the same decision
#       tree as Script 6b, so 6b need not be rerun.
#
# Script 7 is called with environment variables (FFA_THRESH_RATE,
# FFA_PEAKS_FRAC_FILE, FFA_OUT_DIR); its defaults are unchanged.
#
# Run from the repo root:  Rscript _documentation/poster_sensitivity.R
# Takes ~90 s per scenario locally.
#
# OUTPUTS
#   _data/sensitivity/<scenario>/   full Script 7 outputs (git-ignored)
#   _data/sensitivity_ratio_summary.rds
#     scenario, family (S3/S4), setting, method, return_period,
#     n_sites, median, q25, q75, share_ge1
# ============================================================================

library(tidyverse)

SENS_DIR  <- "_data/sensitivity"
BASE_FRAC <- "_data/pot_peaks_classified_frac_thresh_3class.rds"
dir.create(SENS_DIR, showWarnings = FALSE)

# Same decision tree as Script 6b.
relabel_frac <- function(peaks, alpha_r, alpha_s = NULL) {
  ok <- peaks$data_sufficient & !is.na(peaks$rain_frac) & peaks$flux_sufficient
  ok <- dplyr::coalesce(ok, FALSE)
  if (is.null(alpha_s))
    alpha_s <- median(peaks$rain_frac[ok & peaks$rain_frac < alpha_r])
  raw <- dplyr::case_when(!ok                      ~ "unclassified",
                          peaks$rain_frac >= alpha_r ~ "rainfall",
                          peaks$rain_frac <  alpha_s ~ "snowmelt",
                          TRUE                       ~ "rain_on_snow")
  peaks$flood_type_raw <- raw
  peaks$flood_type     <- factor(raw, levels = c("rainfall", "snowmelt", "rain_on_snow", "unclassified"))
  attr(peaks, "mixed_split")      <- alpha_s
  attr(peaks, "rain_frac_thresh") <- alpha_r
  peaks
}

scenarios <- tribble(
  ~scenario,   ~family, ~setting,                 ~rate, ~alpha_r, ~alpha_s,
  "rate_0.5",  "S3",    "0.5 per yr (min 20)",    0.5,   0.99,     NA,
  "rate_2",    "S3",    "2 per yr",               2,     0.99,     NA,
  "rate_3",    "S3",    "3 per yr",               3,     0.99,     NA,
  "ar_0.90",   "S4",    "alpha_r = 0.90",         1,     0.90,     NA,
  "ar_0.95",   "S4",    "alpha_r = 0.95",         1,     0.95,     NA,
  "as_0.20",   "S4",    "alpha_s = 0.20",         1,     0.99,     0.20,
  "as_0.40",   "S4",    "alpha_s = 0.40",         1,     0.99,     0.40
)

base_peaks <- readRDS(BASE_FRAC)

for (i in seq_len(nrow(scenarios))) {
  sc  <- scenarios[i, ]
  out <- file.path(SENS_DIR, sc$scenario)
  if (file.exists(file.path(out, "ffa_quantiles_pgpd.rds"))) {
    message("Skipping ", sc$scenario, " (already run)"); next
  }
  dir.create(out, showWarnings = FALSE)
  peaks_file <- BASE_FRAC
  if (sc$alpha_r != 0.99 || !is.na(sc$alpha_s)) {
    relab <- relabel_frac(base_peaks, sc$alpha_r, if (is.na(sc$alpha_s)) NULL else sc$alpha_s)
    message(sprintf("%s: alpha_r = %.2f, alpha_s = %.4f", sc$scenario,
                    attr(relab, "rain_frac_thresh"), attr(relab, "mixed_split")))
    peaks_file <- file.path(out, "pot_peaks_classified_frac_thresh_3class.rds")
    saveRDS(relab, peaks_file)
  }
  message("Running Script 7: ", sc$scenario)
  status <- system2("Rscript", c("--no-init-file", "_scripts/7_compile_ffa_results.R"),
                    env = c(sprintf("FFA_THRESH_RATE=%g", sc$rate),
                            sprintf("FFA_PEAKS_FRAC_FILE=%s", peaks_file),
                            sprintf("FFA_OUT_DIR=%s", out),
                            "R_PROFILE_USER=/dev/null"),
                    stdout = file.path(out, "script7.log"), stderr = file.path(out, "script7.log"))
  if (status != 0) stop("Script 7 failed for ", sc$scenario, "; see ", file.path(out, "script7.log"))
}

# ---- Summarize the ratio by return period --------------------------------
ratio_summary <- function(q) {
  q %>%
    filter(curve_type %in% c("mixture", "single_sample")) %>%
    select(site_no, method, return_period, curve_type, discharge_cfs) %>%
    pivot_wider(names_from = curve_type, values_from = discharge_cfs) %>%
    filter(is.finite(mixture), is.finite(single_sample)) %>%
    mutate(ratio = mixture / single_sample) %>%
    group_by(method, return_period) %>%
    summarise(n_sites = n(), median = median(ratio),
              q25 = quantile(ratio, .25), q75 = quantile(ratio, .75),
              share_ge1 = mean(ratio >= 1), .groups = "drop")
}

baseline <- ratio_summary(readRDS("_data/ffa_quantiles_pgpd.rds")) %>%
  mutate(scenario = "baseline", family = "baseline", setting = "Baseline")

summary_all <- bind_rows(
  baseline,
  purrr::pmap_dfr(scenarios, function(scenario, family, setting, ...) {
    ratio_summary(readRDS(file.path(SENS_DIR, scenario, "ffa_quantiles_pgpd.rds"))) %>%
      mutate(scenario = scenario, family = family, setting = setting)
  })
) %>% select(scenario, family, setting, everything())

saveRDS(summary_all, "_data/sensitivity_ratio_summary.rds")

message("\nMedian ratio, rain fraction (3-class):")
summary_all %>%
  filter(method == "frac_thresh_3class", return_period %in% c(2, 10, 100, 500)) %>%
  select(scenario, return_period, median, share_ge1, n_sites) %>%
  mutate(across(c(median, share_ge1), ~ round(.x, 2))) %>%
  pivot_wider(names_from = return_period, values_from = c(median, share_ge1, n_sites)) %>%
  print(width = Inf)
