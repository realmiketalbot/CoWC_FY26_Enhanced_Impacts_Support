#!/usr/bin/env Rscript

# ============================================================================
# check_m1_invariants.R
#
# Phase U invariant checks for milestone M1 (_planning/agu2026_interim_plan.md,
# item U4). Compares the regenerated upstream outputs against the frozen
# baseline at git tag `initial-results-prelim-proposal`.
#
# Every Phase U change is a rename, an addition, or the removal of
# downstream-only columns. None should change which peaks exist, how they
# are labelled, or any trend output. Any FAIL needs an explanation before
# Script 7 work starts.
#
# USAGE (from the repo root):
#   Rscript _scripts/checks/check_m1_invariants.R [data_dir] [groups]
#     data_dir : directory holding the regenerated .rds files (default _data)
#     groups   : comma-separated subset of script3,script5,classify,trends
#                (default: all)
#
# Baseline files are read with `git show <tag>:_data/<file>`, so this must
# run inside a clone that has the tag.
#
# Exit status is 1 if any check fails.
# ============================================================================

suppressMessages(library(tidyverse))

args     <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "_data"
groups   <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else
  c("script3", "script5", "classify", "trends")

TAG <- "initial-results-prelim-proposal"

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
repo_root   <- normalizePath(file.path(dirname(script_path), "..", ".."))
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
record <- function(check, pass, detail = "") {
  results[[length(results) + 1]] <<- tibble(check = check, pass = pass, detail = detail)
  message(sprintf("[%s] %s%s", if (pass) "PASS" else "FAIL", check,
                  if (nzchar(detail)) paste0(" -- ", detail) else ""))
}

# Compare two data frames on shared columns after sorting by keys.
same_frame <- function(a, b, keys, cols = setdiff(intersect(names(a), names(b)), keys)) {
  a <- as_tibble(a) %>% select(all_of(c(keys, cols))) %>% arrange(across(all_of(keys)))
  b <- as_tibble(b) %>% select(all_of(c(keys, cols))) %>% arrange(across(all_of(keys)))
  a[] <- lapply(a, function(x) { if (is.factor(x)) x <- as.character(x); unname(x) })
  b[] <- lapply(b, function(x) { if (is.factor(x)) x <- as.character(x); unname(x) })
  res <- all.equal(a, b, check.attributes = FALSE)
  list(pass = isTRUE(res), detail = if (isTRUE(res)) "" else paste(head(res, 3), collapse = "; "))
}

# Report key differences between two peak tables.
key_diff <- function(old, new, keys = c("site_no", "peak_id")) {
  old_k <- distinct(as_tibble(old), across(all_of(keys))) %>% mutate(across(everything(), as.character))
  new_k <- distinct(as_tibble(new), across(all_of(keys))) %>% mutate(across(everything(), as.character))
  list(only_old = anti_join(old_k, new_k, by = keys),
       only_new = anti_join(new_k, old_k, by = keys))
}

# ---------------------------------------------------------------------------
# Script 3: candidate peaks (only the threshold column is renamed)
# ---------------------------------------------------------------------------
if ("script3" %in% groups) {
  old <- read_baseline("pot_peaks_all.rds") %>% rename(candidate_threshold_cfs = pot_threshold_cfs)
  new <- read_new("pot_peaks_all.rds")
  r <- same_frame(old, new, keys = c("site_no", "peak_id", "date"))
  record("Script 3: pot_peaks_all identical apart from the rename", r$pass, r$detail)

  rl <- read_new("site_record_length.rds")
  record("Script 3: site_record_length covers every site with peaks",
         setequal(unique(new$site_no), rl$site_no),
         sprintf("%d sites, median eff_years %.1f", nrow(rl), median(rl$eff_years)))
}

# ---------------------------------------------------------------------------
# Script 5: accumulation curves and endmembers (inputs unchanged)
# ---------------------------------------------------------------------------
if ("script5" %in% groups) {
  r <- same_frame(read_baseline("accum_curves_60d.rds"), read_new("accum_curves_60d.rds"),
                  keys = c("site_no", "peak_id", "lookback_day"))
  record("Script 5: accum_curves_60d identical", r$pass, r$detail)

  eo <- read_baseline("endmembers_60d.rds"); en <- read_new("endmembers_60d.rds")
  r <- same_frame(eo, en, keys = intersect(c("site_no", "peak_id"), names(eo)))
  record("Script 5: endmembers_60d identical", r$pass, r$detail)
}

# ---------------------------------------------------------------------------
# Scripts 6a/6b: labels and features unchanged; selection columns gone
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
    old <- read_baseline(s$file)
    new <- read_new(s$file)

    kd <- key_diff(old, new)
    record(sprintf("%s: same peaks (site_no, peak_id)", m),
           nrow(kd$only_old) == 0 && nrow(kd$only_new) == 0,
           sprintf("old %d, new %d, only-old %d, only-new %d",
                   nrow(old), nrow(new), nrow(kd$only_old), nrow(kd$only_new)))
    if (nrow(kd$only_new) > 0) {
      message("  Peaks only in new output:")
      print(semi_join(mutate(new, peak_id = as.character(peak_id)), kd$only_new,
                      by = c("site_no", "peak_id")) %>%
              select(any_of(c("site_no", "peak_id", "peak_date", "peak_flow",
                              "complete_window", "flood_type_raw"))))
    }

    # Compare on the peaks present in both
    old_c <- mutate(old, peak_id = as.character(peak_id))
    new_c <- mutate(new, peak_id = as.character(peak_id))
    common <- inner_join(distinct(old_c, site_no, peak_id), distinct(new_c, site_no, peak_id),
                         by = c("site_no", "peak_id"))
    old_c <- semi_join(old_c, common, by = c("site_no", "peak_id"))
    new_c <- semi_join(new_c, common, by = c("site_no", "peak_id"))

    r <- same_frame(old_c, new_c, keys = c("site_no", "peak_id"), cols = "flood_type_raw")
    record(sprintf("%s: flood_type_raw identical", m), r$pass, r$detail)

    r <- same_frame(old_c, new_c, keys = c("site_no", "peak_id"), cols = s$features)
    record(sprintf("%s: features identical", m), r$pass, r$detail)

    record(sprintf("%s: fit_eligible matches the old single-sample eligibility rule", m),
           identical(new$fit_eligible, s$old_elig(new)),
           sprintf("%d of %d peaks eligible", sum(new$fit_eligible), nrow(new)))

    record(sprintf("%s: flood_type is an exact copy of flood_type_raw", m),
           identical(as.character(new$flood_type), as.character(new$flood_type_raw)))

    gone <- intersect(c("in_type_sample", "in_single_sample", "pot_threshold_cfs"), names(new))
    record(sprintf("%s: selection and old threshold columns removed", m),
           length(gone) == 0, paste(gone, collapse = ", "))

    record(sprintf("%s: class_config attribute present", m),
           !is.null(attr(new, "class_config")))

    if (m == "frac_thresh_3class") {
      record("frac_thresh_3class: mixed_split unchanged",
             isTRUE(all.equal(attr(old, "mixed_split"), attr(new, "mixed_split"))),
             sprintf("old %.6f, new %.6f", attr(old, "mixed_split"), attr(new, "mixed_split")))
    }
  }
}

# ---------------------------------------------------------------------------
# Scripts 8/9: trend outputs (both read flood_type_raw, so unchanged)
# ---------------------------------------------------------------------------
if ("trends" %in% groups) {
  methods <- c("depth_thresh", "frac_thresh_3class", "frac_thresh_2class")
  files <- c("annual_climate_forcing_summary.rds", "climate_forcing_trends.rds",
             "annual_am_peakflow.rds",
             sprintf("annual_am_type_%s.rds", methods),
             sprintf("am_type_shift_%s.rds", methods),
             sprintf("temporal_trends_am_%s.rds", methods),
             sprintf("temporal_trends_%s.rds", methods),
             sprintf("annual_type_counts_%s.rds", methods))
  for (f in files) {
    res <- all.equal(read_baseline(f), read_new(f), check.attributes = FALSE)
    record(sprintf("Scripts 8/9: %s identical", f), isTRUE(res),
           if (isTRUE(res)) "" else paste(head(res, 3), collapse = "; "))
  }
}

# ---------------------------------------------------------------------------
summary_tbl <- bind_rows(results)
message(sprintf("\n%d of %d checks passed.", sum(summary_tbl$pass), nrow(summary_tbl)))
if (!all(summary_tbl$pass)) quit(status = 1)
