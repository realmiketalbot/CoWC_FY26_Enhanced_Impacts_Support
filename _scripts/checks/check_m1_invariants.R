#!/usr/bin/env Rscript

# ============================================================================
# check_m1_invariants.R
#
# Phase U invariant checks for milestone M1 (_planning/agu2026_interim_plan.md,
# item U4). Compares regenerated upstream outputs against the frozen baseline
# at git tag `initial-results-prelim-proposal`.
#
# SITE-AWARE (revised 2026-11-17 for the station re-screen): the site list
# changed (Scripts 2/2b), so outputs are compared on the sites present in both
# the baseline and the new run ("common sites"). Every site is processed
# independently in Scripts 3, 5 and 6a, so common sites must match exactly.
#
# One expected exception: Script 6b's snowmelt / rain-on-snow split
# (mixed_split, alpha_s) is a median over ALL sites, so it moves when the
# site set changes. Peaks whose rain_frac lies between the old and new split
# may legitimately change label. The check verifies that every 3-class label
# change is of that kind. Outputs that depend on the 3-class labels are
# reported as INFO rather than PASS/FAIL. The 2-class labels
# (rainfall vs snow_influenced) do not depend on the split and are checked
# strictly.
#
# USAGE (from the repo root):
#   Rscript _scripts/checks/check_m1_invariants.R [data_dir] [groups]
#     data_dir : directory holding the regenerated .rds files (default _data)
#     groups   : comma-separated subset of script3,script5,classify,trends
#                (default: all)
#
# Baseline files are read with `git show <tag>:_data/<file>`, so this must
# run inside a clone that has the tag. Exit status is 1 if any check FAILs.
# ============================================================================

suppressMessages(library(tidyverse))

args     <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "_data"
groups   <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else
  c("script3", "script5", "classify", "trends")

TAG <- "initial-results-prelim-proposal"

script_path  <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
repo_root    <- normalizePath(file.path(dirname(script_path), "..", ".."))
baseline_dir <- tempfile("m1_baseline_"); dir.create(baseline_dir)

read_baseline <- function(f) {
  out <- file.path(baseline_dir, f)
  status <- system2("git", c("-C", shQuote(repo_root), "show",
                             shQuote(sprintf("%s:_data/%s", TAG, f))),
                    stdout = out, stderr = FALSE)
  if (!identical(status, 0L) || file.size(out) == 0)
    stop(sprintf("Could not read %s from tag %s", f, TAG))
  readRDS(out)
}
read_new <- function(f) readRDS(file.path(data_dir, f))

results <- list()
record <- function(check, pass, detail = "", info = FALSE) {
  status <- if (info) "INFO" else if (pass) "PASS" else "FAIL"
  results[[length(results) + 1]] <<- tibble(check = check, status = status, detail = detail)
  message(sprintf("[%s] %s%s", status, check, if (nzchar(detail)) paste0(" -- ", detail) else ""))
}

# Normalise a data frame for comparison: plain tibble, factors as character,
# names dropped, sorted by keys.
norm <- function(d, keys, cols) {
  d <- as_tibble(d) %>% select(all_of(c(keys, cols)))
  d[] <- lapply(d, function(x) { if (is.factor(x)) x <- as.character(x); unname(x) })
  arrange(d, across(all_of(keys)))
}

# Compare two data frames restricted to `sites` on shared columns.
same_on_sites <- function(a, b, keys, sites, cols = NULL, exclude = character()) {
  if (is.null(cols)) cols <- setdiff(intersect(names(a), names(b)), c(keys, exclude))
  a <- norm(filter(as_tibble(a), site_no %in% sites), keys, cols)
  b <- norm(filter(as_tibble(b), site_no %in% sites), keys, cols)
  res <- all.equal(a, b, check.attributes = FALSE)
  list(pass = isTRUE(res),
       detail = if (isTRUE(res)) sprintf("%d rows", nrow(a)) else paste(head(res, 3), collapse = "; "))
}

# ---------------------------------------------------------------------------
# Site sets
# ---------------------------------------------------------------------------
old_sites <- unique(read_baseline("pot_peaks_all.rds")$site_no)
new_sites <- read_csv(file.path(data_dir, "selected_unit_stations.csv"),
                      show_col_types = FALSE, col_types = cols(site_no = col_character()))$site_no
common    <- intersect(old_sites, new_sites)
message(sprintf("Sites: baseline %d, new selection %d, common %d (dropped %d, added %d)\n",
                length(old_sites), length(new_sites), length(common),
                length(setdiff(old_sites, new_sites)), length(setdiff(new_sites, old_sites))))

# ---------------------------------------------------------------------------
# Script 3: candidate peaks
# ---------------------------------------------------------------------------
if ("script3" %in% groups) {
  old <- read_baseline("pot_peaks_all.rds") %>% rename(candidate_threshold_cfs = pot_threshold_cfs)
  new <- read_new("pot_peaks_all.rds")

  record("Script 3: peaks exist for exactly the selected sites",
         setequal(unique(new$site_no), new_sites),
         sprintf("%d sites with peaks; selected %d", n_distinct(new$site_no), length(new_sites)))

  r <- same_on_sites(old, new, keys = c("site_no", "peak_id", "date"), sites = common)
  record("Script 3: pot_peaks_all identical on common sites (apart from the rename)", r$pass, r$detail)

  q <- read_new("site_qualifying_years.rds")
  record("Script 2b: every selected site has qualifying-year records",
         all(new_sites %in% q$site_no[q$selected]))
}

# ---------------------------------------------------------------------------
# Script 5: accumulation curves and endmembers
# ---------------------------------------------------------------------------
if ("script5" %in% groups) {
  r <- same_on_sites(read_baseline("accum_curves_60d.rds"), read_new("accum_curves_60d.rds"),
                     keys = c("site_no", "peak_id", "lookback_day"), sites = common)
  record("Script 5: accum_curves_60d identical on common sites", r$pass, r$detail)

  eo <- read_baseline("endmembers_60d.rds"); en <- read_new("endmembers_60d.rds")
  r <- same_on_sites(eo, en, keys = intersect(c("site_no", "peak_id"), names(eo)), sites = common)
  record("Script 5: endmembers_60d identical on common sites", r$pass, r$detail)
}

# ---------------------------------------------------------------------------
# Scripts 6a/6b: labels and features
# ---------------------------------------------------------------------------
if ("classify" %in% groups) {
  specs <- list(
    depth_thresh = list(
      file     = "pot_peaks_classified_depth_thresh.rds",
      features = c("peak_date", "peak_flow", "total_rain_mm", "total_melt_mm",
                   "n_days_wldas", "complete_window", "regulated", "water_year"),
      old_elig = function(d) coalesce(d$complete_window, FALSE)
    ),
    frac_thresh_3class = list(
      file     = "pot_peaks_classified_frac_thresh_3class.rds",
      features = c("peak_date", "peak_flow", "total_rain_flux", "total_melt_flux",
                   "total_flux", "n_days_data", "rain_frac", "data_sufficient",
                   "flux_sufficient", "regulated", "water_year"),
      old_elig = function(d) coalesce(d$data_sufficient & !is.na(d$rain_frac) & d$flux_sufficient, FALSE)
    )
  )

  for (m in names(specs)) {
    s   <- specs[[m]]
    old <- read_baseline(s$file) %>% mutate(peak_id = as.character(peak_id)) %>% filter(site_no %in% common)
    new_all <- read_new(s$file)
    new <- new_all %>% mutate(peak_id = as.character(peak_id)) %>% filter(site_no %in% common)

    keys_old <- distinct(old, site_no, peak_id); keys_new <- distinct(new, site_no, peak_id)
    record(sprintf("%s: same peaks on common sites", m),
           nrow(anti_join(keys_old, keys_new, by = c("site_no", "peak_id"))) == 0 &&
             nrow(anti_join(keys_new, keys_old, by = c("site_no", "peak_id"))) == 0,
           sprintf("old %d, new %d", nrow(keys_old), nrow(keys_new)))

    r <- same_on_sites(old, new, keys = c("site_no", "peak_id"), sites = common, cols = s$features)
    record(sprintf("%s: features identical on common sites", m), r$pass, r$detail)

    lab <- inner_join(old %>% select(site_no, peak_id, old_lab = flood_type_raw),
                      new %>% select(site_no, peak_id, new_lab = flood_type_raw, any_of("rain_frac")),
                      by = c("site_no", "peak_id")) %>%
      mutate(across(c(old_lab, new_lab), as.character), changed = old_lab != new_lab)

    if (m == "depth_thresh") {
      record("depth_thresh: flood_type_raw identical on common sites", !any(lab$changed),
             sprintf("%d of %d labels changed", sum(lab$changed), nrow(lab)))
    } else {
      s_old <- attr(read_baseline(s$file), "mixed_split")
      s_new <- attr(new_all, "mixed_split")
      lo <- min(s_old, s_new); hi <- max(s_old, s_new)
      chg <- filter(lab, changed)
      explained <- all(chg$old_lab %in% c("snowmelt", "rain_on_snow") &
                       chg$new_lab %in% c("snowmelt", "rain_on_snow") &
                       chg$rain_frac >= lo & chg$rain_frac <= hi)
      record("frac_thresh_3class: mixed_split moved with the site set", TRUE,
             sprintf("old %.6f, new %.6f", s_old, s_new), info = TRUE)
      record("frac_thresh_3class: every label change is a snowmelt/rain-on-snow swap between the old and new split",
             explained, sprintf("%d of %d labels changed", nrow(chg), nrow(lab)))
      lab2 <- lab %>% mutate(across(c(old_lab, new_lab),
                                    ~ ifelse(.x %in% c("snowmelt", "rain_on_snow"), "snow_influenced", .x)))
      record("frac_thresh_2class: labels identical on common sites (split-independent)",
             all(lab2$old_lab == lab2$new_lab))
    }

    record(sprintf("%s: fit_eligible matches the old single-sample eligibility rule", m),
           identical(new_all$fit_eligible, s$old_elig(new_all)),
           sprintf("%d of %d peaks eligible", sum(new_all$fit_eligible), nrow(new_all)))
    record(sprintf("%s: flood_type is an exact copy of flood_type_raw", m),
           identical(as.character(new_all$flood_type), as.character(new_all$flood_type_raw)))
    gone <- intersect(c("in_type_sample", "in_single_sample", "pot_threshold_cfs"), names(new_all))
    record(sprintf("%s: selection and old threshold columns removed", m),
           length(gone) == 0, paste(gone, collapse = ", "))
    record(sprintf("%s: class_config attribute present", m), !is.null(attr(new_all, "class_config")))
  }
}

# ---------------------------------------------------------------------------
# Scripts 8/9: trend outputs on common sites
# ---------------------------------------------------------------------------
if ("trends" %in% groups) {
  # Strict: independent of the 3-class snowmelt/rain-on-snow split.
  strict <- list(
    list(f = "annual_climate_forcing_summary.rds", keys = c("site_no", "water_year"), exclude = "am_flood_type"),
    list(f = "climate_forcing_trends.rds",         keys = c("site_no", "variable")),
    list(f = "annual_am_peakflow.rds",             keys = c("site_no", "water_year"))
  )
  for (meth in c("depth_thresh", "frac_thresh_2class")) {
    strict <- c(strict, list(
      list(f = sprintf("annual_am_type_%s.rds", meth),      keys = c("site_no", "water_year")),
      list(f = sprintf("am_type_shift_%s.rds", meth),       keys = c("site_no", "flood_type")),
      list(f = sprintf("annual_type_counts_%s.rds", meth),  keys = c("site_no", "water_year", "flood_type")),
      list(f = sprintf("temporal_trends_am_%s.rds", meth),  el = "am_type_indicator_trends", keys = c("site_no", "flood_type")),
      list(f = sprintf("temporal_trends_am_%s.rds", meth),  el = "am_peak_flow_trends",      keys = "site_no"),
      list(f = sprintf("temporal_trends_%s.rds", meth),     el = "site_level", keys = c("site_no", "flood_type", "metric"))
    ))
  }
  strict <- c(strict, list(
    list(f = "temporal_trends_am_frac_thresh_3class.rds", el = "am_peak_flow_trends", keys = "site_no"),
    list(f = "temporal_trends_am_frac_thresh_3class.rds", el = "am_rain_frac_trends", keys = "site_no")
  ))

  get_el <- function(x, el) if (is.null(el)) x else x[[el]]
  for (s in strict) {
    a <- get_el(read_baseline(s$f), s$el); b <- get_el(read_new(s$f), s$el)
    r <- same_on_sites(a, b, keys = s$keys, sites = common,
                       exclude = if (is.null(s$exclude)) character() else s$exclude)
    record(sprintf("Scripts 8/9: %s%s identical on common sites", s$f,
                   if (is.null(s$el)) "" else paste0("$", s$el)), r$pass, r$detail)
  }

  # INFO: depend on 3-class labels, which can shift with mixed_split
  info <- list(
    list(f = "annual_am_type_frac_thresh_3class.rds",     keys = c("site_no", "water_year")),
    list(f = "annual_type_counts_frac_thresh_3class.rds", keys = c("site_no", "water_year", "flood_type")),
    list(f = "temporal_trends_frac_thresh_3class.rds",    el = "site_level", keys = c("site_no", "flood_type", "metric"))
  )
  for (s in info) {
    a <- get_el(read_baseline(s$f), s$el); b <- get_el(read_new(s$f), s$el)
    r <- same_on_sites(a, b, keys = s$keys, sites = common)
    record(sprintf("Scripts 8/9: %s%s on common sites (3-class, split-dependent)", s$f,
                   if (is.null(s$el)) "" else paste0("$", s$el)),
           r$pass, if (r$pass) "identical" else "differs (expected only via mixed_split)", info = TRUE)
  }
}

# ---------------------------------------------------------------------------
summary_tbl <- bind_rows(results)
n_fail <- sum(summary_tbl$status == "FAIL")
message(sprintf("\n%d PASS, %d FAIL, %d INFO.", sum(summary_tbl$status == "PASS"), n_fail,
                sum(summary_tbl$status == "INFO")))
if (n_fail > 0) quit(status = 1)
