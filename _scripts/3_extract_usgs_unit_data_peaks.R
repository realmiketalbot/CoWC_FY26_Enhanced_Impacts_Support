#!/usr/bin/env Rscript

# ============================================================================
# extract_usgs_unit_data_peaks.R
#
# Find all independent candidate peaks for selected unit-data stations.
# Implements independence criteria from USWRC (1976), as reviewed in
# Lang et al. (1999).
# All independent peaks are retained for later flood type classification.
#
# The per-site candidate threshold (q_threshold quantile of daily max flow
# over the full record) is saved as `candidate_threshold_cfs`. It filters
# local maxima into candidate peaks and does nothing else. It is NOT the
# GPD/Poisson-GPD model-fitting threshold: that threshold (u_i) is chosen
# per flood type downstream in Script 7. (Renamed from `pot_threshold_cfs`,
# which the initial results also used as the GPD location; see
# _planning/agu2026_interim_plan.md.)
#
# The effective record length per site is also computed from the daily
# series and saved to _data/site_record_length.rds, as the denominator of
# the Poisson exceedance rate (lambda) and for the dispersion index.
#
# CHANGE FROM PREVIOUS VERSION -----------------------------------------------
# Peak identification is now a three-stage process:
#
#   1. Local maxima are detected across the full daily series using
#      pracma::findpeaks, with NO threshold applied.
#   2. Local maxima are filtered to those exceeding the per-site candidate
#      threshold.
#   3. The USWRC independence criteria are applied by find_independent_peaks().
#
# Previously, stage 1 segmented the record into runs of consecutive days above
# the threshold and kept one maximum per run. At a median-flow threshold those
# runs span a large fraction of the record: a network-wide diagnostic found a
# median of 82% of above-threshold days sitting in runs longer than 2*theta,
# with maximum run lengths exceeding 1000 days. Any additional peak inside a
# long run was therefore discarded before declustering could evaluate it, and
# declustering can only remove peaks, never recover them.
#
# Detecting local maxima before applying the threshold means no candidate is
# eliminated by an arbitrary segmentation; the threshold and the independence
# criteria do all of the filtering, in that order.
#
# find_independent_peaks() has also been optimised to compare each candidate
# only against its immediate date-neighbours among already-accepted peaks.
# This returns an identical peak set (see the note above the function) but
# reduces complexity from O(n*k) to O(n log k), which matters because stage 1
# now presents far more candidates than run segmentation did.
#
# INPUTS:
#   - selected_unit_stations.csv   : site_no column of selected stations
#   - /path/to/usgs/unit/<site_no>/<year>.rds : unit discharge files
#   - conterm_basinid.txt          : GAGES-II metadata (drainage area)
#
# OUTPUTS:
#   - pot_peaks_all.rds : one row per site-day inside each peak window,
#                         with columns:
#                           site_no, date, peak_id, peak_date, peak_flow,
#                           candidate_threshold_cfs, search_radius,
#                           antecedent_runoff_mm, antecedent_days_used
#   - site_record_length.rds : one row per site, with columns:
#       site_no, n_days_obs   : days with a daily-max value
#       eff_years             : n_days_obs / 365.25 (lambda denominator)
#       wy_first, wy_last     : first and last water year with any data
#       n_wy_complete         : water years with >= WY_COMPLETE_FRAC of days
#                               observed (dispersion-index sample size M)
#       complete_wys          : list-column of those water years
# ============================================================================

library(tidyverse)
library(pracma)

# ============================================================================
# CONFIGURATION
# ============================================================================
selected_stations_csv <- "_data/selected_unit_stations.csv"
usgs_unit_base_dir    <- "/davenport-scratch/DATA/USGS/raw/Discharge/Unit"
watershed_meta        <- "/davenport-scratch/DATA/GAGES-II/basinchar_and_report_sept_2011/spreadsheets-in-csv-format/conterm_basinid.txt"
output_file           <- "_data/pot_peaks_all.rds"
record_length_file    <- "_data/site_record_length.rds"

q_threshold <- 0.5   # quantile of daily max flow used as the candidate threshold

# A water year counts as complete for the dispersion index when at least
# this fraction of its days has a daily-max value (decision D-E).
WY_COMPLETE_FRAC <- 0.90

# findpeaks shape parameters. nups/ndowns = 1 is the minimal definition of a
# local maximum (one rise before, one fall after). Raising them requires a
# more sustained rise and recession and yields fewer, cleaner candidates.
NUPS   <- 1
NDOWNS <- 1

# ============================================================================
# FUNCTION: Independence criteria from USWRC (1976), per Lang et al. (1999)
#
# Applies two criteria to decluster candidate peaks for a single watershed:
#   1. Minimum time separation: theta > 5 + ln(0.3861 * A), A in km²
#      (the criterion is stated with A in miles²; 0.3861 converts from km²)
#   2. Inter-peak recession: flow between peaks must drop below 75% of the
#      lesser of the two peak flows
#
# Candidates are considered in descending order of magnitude so that when two
# peaks are dependent, the LARGER is retained.
#
# OPTIMISATION: a candidate is compared only against its immediate date-
# neighbours in the set of already-accepted peaks, rather than against every
# accepted peak. This returns the identical peak set. Let i be the candidate
# and m its nearest accepted neighbour on one side, and let k be any farther
# accepted peak on that same side. If i passes both criteria against m:
#
#   Separation:  |t_k - t_i| > |t_m - t_i| >= theta, so criterion 1 holds.
#
#   Recession:   the minimum flow over a longer interval is no greater than
#                over any sub-interval, so
#                  min(k,i) <= min(m,i) < 0.75 * min(Q_m, Q_i) <= 0.75 * Q_i
#                and, because k and m were both accepted and therefore satisfy
#                the criterion between themselves,
#                  min(k,i) <= min(k,m) < 0.75 * min(Q_k, Q_m) <= 0.75 * Q_k
#                Together: min(k,i) < 0.75 * min(Q_k, Q_i), so criterion 2
#                holds as well.
#
# Independence from the immediate neighbours therefore implies independence
# from every accepted peak beyond them.
# ============================================================================
find_independent_peaks <- function(peak_df, all_flows_df, watershed_drain_sqkm) {

  n <- nrow(peak_df)
  if (n < 2) return(peak_df)

  min_sep_days <- 5 + log(0.3861 * watershed_drain_sqkm)

  # daily flow series, sorted once
  ord        <- order(all_flows_df$date)
  flow_dates <- all_flows_df$date[ord]
  flow_vals  <- all_flows_df$q_daily_max[ord]

  pk_date     <- peak_df$date
  pk_date_num <- as.numeric(pk_date)
  pk_flow     <- peak_df$q_daily_max

  # minimum daily flow strictly between two dates; Inf if no days between
  min_flow_between <- function(d1, d2) {
    lo <- findInterval(d1, flow_dates) + 1L
    hi <- findInterval(d2 - 1, flow_dates)
    if (hi >= lo) min(flow_vals[lo:hi], na.rm = TRUE) else Inf
  }

  # TRUE if peaks i and j fail either independence criterion
  is_dependent <- function(i, j) {

    if (abs(pk_date_num[i] - pk_date_num[j]) < min_sep_days) return(TRUE)

    d1 <- min(pk_date[i], pk_date[j])
    d2 <- max(pk_date[i], pk_date[j])
    mf <- min_flow_between(d1, d2)

    !isTRUE(mf < 0.75 * min(pk_flow[i], pk_flow[j]))
  }

  order_by_size <- order(pk_flow, decreasing = TRUE)

  kept_idx   <- integer(0)   # indices into peak_df, held in DATE order
  kept_dates <- numeric(0)   # matching dates, kept sorted for findInterval

  for (i in order_by_size) {

    di  <- pk_date_num[i]
    pos <- findInterval(di, kept_dates)   # number of accepted dates <= di

    ok <- TRUE

    # nearest accepted peak at or before this date
    if (pos >= 1L) {
      ok <- !is_dependent(i, kept_idx[pos])
    }

    # nearest accepted peak after this date
    if (ok && pos < length(kept_idx)) {
      ok <- !is_dependent(i, kept_idx[pos + 1L])
    }

    if (ok) {
      kept_idx   <- append(kept_idx,   i,  after = pos)
      kept_dates <- append(kept_dates, di, after = pos)
    }
  }

  peak_df[sort(kept_idx), ]
}

# ============================================================================
# FUNCTION: local maxima detection for one site
#
# pracma::findpeaks operates on vector position, not on dates, so the series
# is first split into segments of consecutive calendar days. Without this, a
# gap in the record would allow two observations months apart to be treated
# as neighbours and produce a spurious local maximum at the gap edge.
#
# findpeaks returns a matrix whose columns are: peak height, index of the
# peak, index where the rise begins, index where the fall ends. Only the
# second column is used here. VERIFY this column order against the pracma
# documentation for your installed version before trusting it network-wide.
#
# Segments shorter than nups + ndowns + 1 days cannot contain a local maximum
# by this definition and are skipped.
# ============================================================================
find_local_maxima <- function(flow_df, n_ups = NUPS, n_downs = NDOWNS) {

  flow_df <- flow_df %>%
    filter(!is.na(q_daily_max)) %>%
    arrange(date)

  min_len <- n_ups + n_downs + 1L
  if (nrow(flow_df) < min_len) return(flow_df[0, ])

  gap     <- c(NA_real_, as.numeric(diff(flow_df$date), units = "days"))
  segment <- cumsum(is.na(gap) | gap != 1)

  map_dfr(unique(segment), function(s) {

    seg <- flow_df[segment == s, ]
    if (nrow(seg) < 1L) return(seg[0, ])

    # padding: Pad each segment with a descending ramp on both sides so that a peak at
    # the first or last observation can still be detected. Without padding,
    # findpeaks requires n_downs falling values after a peak, so a peak
    # immediately preceding a data gap is invisible -- a failure mode biased
    # toward the largest floods, since gauges often fail during extreme events.     
    
    # # breaking ties: USGS daily values are rounded, so broad snowmelt crests
    # frequently produce plateaus. findpeaks requires strict rises and falls,
    # so a plateau is not recognised as a peak at all. A monotonic
    # perturbation far below the data's precision resolves each plateau to
    # its first day without altering any real ordering.

    eps <- 1e-6 * max(seg$q_daily_max, na.rm = TRUE)
    q_jit <- seg$q_daily_max - eps * seq_len(nrow(seg))

    padded <- c(min(q_jit) - n_ups:1, q_jit, min(q_jit) - 1:n_downs)

    pk <- pracma::findpeaks(padded, nups = n_ups, ndowns = n_downs)
    if (is.null(pk)) return(seg[0, ])

    # shift indices back onto the unpadded segment and drop anything that
    # falls outside it
    idx <- pk[, 2] - n_ups
    idx <- idx[idx >= 1L & idx <= nrow(seg)]

    seg[idx, ]
  })
}

# ============================================================================
# LOAD STATION LIST AND DRAINAGE AREA METADATA
# ============================================================================
selected_stations <- read_csv(selected_stations_csv, show_col_types = FALSE) %>%
  pull(site_no) %>%
  as.character()

message(sprintf("Loaded %d selected stations", length(selected_stations)))

basin_meta <- read_csv(watershed_meta, show_col_types = FALSE) %>%
  filter(STAID %in% selected_stations) %>%
  transmute(
    site_no    = as.character(STAID),
    drain_sqkm = DRAIN_SQKM
  )

# ============================================================================
# READ UNIT DATA AND AGGREGATE TO DAILY MAX
# ============================================================================
read_daily_max <- function(site_no, base_dir) {

  site_path  <- file.path(base_dir, site_no)
  year_files <- list.files(site_path, pattern = "^\\d{4}\\.rds$", full.names = TRUE)

  if (length(year_files) == 0) {
    warning(sprintf("No year files found for site %s", site_no))
    return(NULL)
  }

  map_dfr(year_files, function(f) {
    tryCatch({
      readRDS(f) %>%
        as_tibble() %>%
        mutate(
          date  = as.Date(dateTime),
          q_cfs = X_00060_00000
        ) %>%
        filter(!is.na(q_cfs), q_cfs >= 0) %>%
        group_by(date) %>%
        summarize(q_daily_max = max(q_cfs, na.rm = TRUE), .groups = "drop")
    }, error = function(e) {
      warning(sprintf("Error reading %s: %s", f, e$message))
      NULL
    })
  })
}

message("Reading unit data and aggregating to daily max...")

all_flows <- map_dfr(selected_stations, function(s) {
  message(sprintf("  Reading site %s", s))
  daily <- read_daily_max(s, usgs_unit_base_dir)
  if (is.null(daily) || nrow(daily) == 0) return(NULL)
  daily %>% mutate(site_no = s)
})

message(sprintf("Loaded daily max data: %d site-day rows across %d sites",
                nrow(all_flows), n_distinct(all_flows$site_no)))

# ============================================================================
# JOIN DRAINAGE AREA — drop sites with missing metadata
# ============================================================================
all_flows <- all_flows %>%
  left_join(basin_meta, by = "site_no")

missing_area <- all_flows %>%
  filter(is.na(drain_sqkm)) %>%
  distinct(site_no)

if (nrow(missing_area) > 0) {
  warning(sprintf(
    "Missing drainage area for %d site(s); dropping: %s",
    nrow(missing_area),
    paste(missing_area$site_no, collapse = ", ")
  ))
  all_flows <- all_flows %>% filter(!is.na(drain_sqkm))
}

# Convert discharge to mm/day over basin area, making it directly comparable
# to WLDAS precipitation depths for mass balance diagnostics.
# q_mm_day = q_cfs * (m3/s per cfs) * (s/day) / (basin area in m2) * (mm/m)
all_flows <- all_flows %>%
  mutate(q_mm_day = q_daily_max * 0.0283168 * 86400 /
                    (drain_sqkm * 1e6) * 1000)

# ============================================================================
# EFFECTIVE RECORD LENGTH PER SITE (decision D-E)
#
# all_flows holds one row per site-day that has a daily-max value, so a row
# count is the number of observed days.
#   eff_years     : observed days / 365.25. Denominator of the Poisson
#                   exceedance rate. Peaks can only be detected on observed
#                   days, so numerator and denominator cover the same period.
#   n_wy_complete : water years (Oct 1 - Sep 30) with >= WY_COMPLETE_FRAC of
#                   their days observed. Used as M in the dispersion index,
#                   where partially observed years would understate annual
#                   counts and inflate the variance.
# ============================================================================
wy_coverage <- all_flows %>%
  distinct(site_no, date) %>%
  mutate(water_year = as.integer(lubridate::year(date)) +
           as.integer(lubridate::month(date) >= 10L)) %>%
  count(site_no, water_year, name = "n_days") %>%
  mutate(
    wy_days  = as.integer(as.Date(sprintf("%d-10-01", water_year)) -
                          as.Date(sprintf("%d-10-01", water_year - 1L))),
    frac_obs = n_days / wy_days
  )

site_record_length <- wy_coverage %>%
  group_by(site_no) %>%
  summarize(
    n_days_obs    = sum(n_days),
    eff_years     = n_days_obs / 365.25,
    wy_first      = min(water_year),
    wy_last       = max(water_year),
    n_wy_complete = sum(frac_obs >= WY_COMPLETE_FRAC),
    complete_wys  = list(sort(water_year[frac_obs >= WY_COMPLETE_FRAC])),
    .groups       = "drop"
  )

attr(site_record_length, "wy_complete_frac") <- WY_COMPLETE_FRAC

message(sprintf(
  "Record length: median %.1f effective years (range %.1f-%.1f); median %d complete water years",
  median(site_record_length$eff_years), min(site_record_length$eff_years),
  max(site_record_length$eff_years), as.integer(median(site_record_length$n_wy_complete))
))

# ============================================================================
# COMPUTE PER-SITE CANDIDATE THRESHOLD
#
# The q_threshold quantile of daily max flow over the full record is the
# minimum flow a local maximum must exceed to become a candidate peak. It is
# not a model-fitting threshold (see header).
#
# Saved separately so it can be joined into the output without needing
# access to the raw daily flow record downstream.
# ============================================================================
site_thresholds <- all_flows %>%
  group_by(site_no) %>%
  summarize(
    candidate_threshold_cfs = unname(quantile(q_daily_max, q_threshold, na.rm = TRUE)),
    .groups = "drop"
  )

message(sprintf(
  "Computed candidate thresholds for %d sites (q=%.2f of daily max flow)",
  nrow(site_thresholds), q_threshold
))

# ============================================================================
# MINIMUM SEPARATION WINDOW
# theta = 5 + ln(0.3861 * A), A in km²
# ============================================================================
search_window <- basin_meta %>%
  filter(site_no %in% unique(all_flows$site_no)) %>%
  transmute(
    site_no,
    search_radius = as.integer(pmax(1L, ceiling(5 + log(0.3861 * drain_sqkm))))
  )

# Declustering boundary check (Eq. 1.1 uses theta > ...; is_dependent() treats
# a gap < theta as dependent, so a gap exactly equal to theta counts as
# independent). Gaps are whole days, so the two only disagree if theta is a
# whole number at some site.
theta_check <- basin_meta %>%
  filter(site_no %in% unique(all_flows$site_no)) %>%
  mutate(theta = 5 + log(0.3861 * drain_sqkm),
         dist_to_integer = abs(theta - round(theta)))

if (any(theta_check$dist_to_integer < 1e-9)) {
  warning(sprintf(
    "theta is a whole number at %d site(s); the '>' vs '>=' boundary in Eq. 1.1 can bind there: %s",
    sum(theta_check$dist_to_integer < 1e-9),
    paste(theta_check$site_no[theta_check$dist_to_integer < 1e-9], collapse = ", ")
  ))
} else {
  message(sprintf(
    "Declustering boundary check: theta is non-integer at all %d sites (closest: %.2e days); '>' vs '>=' cannot bind",
    nrow(theta_check), min(theta_check$dist_to_integer)
  ))
}

# ============================================================================
# STAGE 1: DETECT LOCAL MAXIMA (no threshold applied)
# ============================================================================
message("Stage 1: detecting local maxima with pracma::findpeaks ...")

sites_present <- unique(all_flows$site_no)

local_maxima <- map_dfr(sites_present, function(s) {
  find_local_maxima(filter(all_flows, site_no == s))
})

message(sprintf("  %d local maxima across %d sites (%.1f%% of all site-days)",
                nrow(local_maxima), n_distinct(local_maxima$site_no),
                100 * nrow(local_maxima) / nrow(all_flows)))

# ============================================================================
# STAGE 2: FILTER TO LOCAL MAXIMA ABOVE THE CANDIDATE THRESHOLD
# ============================================================================
candidate_peaks <- local_maxima %>%
  left_join(site_thresholds, by = "site_no") %>%
  filter(q_daily_max > candidate_threshold_cfs) %>%
  select(-candidate_threshold_cfs) %>%
  arrange(site_no, date)

message(sprintf(
  "Stage 2: %d local maxima exceed the per-site candidate threshold (%.1f%% retained)",
  nrow(candidate_peaks),
  100 * nrow(candidate_peaks) / nrow(local_maxima)
))

# ============================================================================
# STAGE 3: DECLUSTER USING USWRC INDEPENDENCE CRITERIA
# ============================================================================
message("Stage 3: applying independence criteria ...")

independent_peaks <- candidate_peaks %>%
  group_by(site_no, drain_sqkm) %>%
  group_modify(~ find_independent_peaks(
    .x,
    all_flows_df         = filter(all_flows, site_no == .y$site_no),
    watershed_drain_sqkm = .y$drain_sqkm
  )) %>%
  ungroup()

message(sprintf(
  "Found %d independent peaks across %d sites (%.1f%% of candidates retained)",
  nrow(independent_peaks), n_distinct(independent_peaks$site_no),
  100 * nrow(independent_peaks) / nrow(candidate_peaks)
))

# ============================================================================
# ANTECEDENT RUNOFF VOLUME PER PEAK (mm over basin area)
#
# For each independent peak, sum q_mm_day over the search_radius days
# PRECEDING the peak date (not including the peak day itself). This is the
# backward-looking antecedent window, symmetric counterpart to the forward
# search window used for peak independence.
#
# Stored as antecedent_runoff_mm for use in mass balance diagnostics against
# WLDAS precipitation totals over the same window.
# ============================================================================
antecedent_runoff <- independent_peaks %>%
  select(site_no, peak_date = date) %>%
  inner_join(search_window, by = "site_no") %>%
  left_join(
    all_flows %>% select(site_no, date, q_mm_day),
    by = "site_no"
  ) %>%
  filter(
    date >= peak_date - search_radius,
    date <  peak_date                   # antecedent only, exclude peak day
  ) %>%
  group_by(site_no, peak_date) %>%
  summarise(
    antecedent_runoff_mm = sum(q_mm_day, na.rm = TRUE),
    antecedent_days_used = sum(!is.na(q_mm_day)),
    .groups = "drop"
  )

message(sprintf("Computed antecedent runoff volumes for %d peaks",
                nrow(antecedent_runoff)))

# ============================================================================
# BUILD BUFFERED SEARCH WINDOWS AROUND EACH PEAK
# ============================================================================
range_df <- all_flows %>%
  group_by(site_no) %>%
  summarize(min_date = min(date), max_date = max(date), .groups = "drop")

buffered_dates <- independent_peaks %>%
  select(site_no, date) %>%
  inner_join(search_window, by = "site_no") %>%
  group_by(site_no) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(
    peak_id   = row_number(),
    peak_date = date,              # keep the peak's own date
    start     = date - search_radius,
    end       = date + search_radius
  ) %>%
  ungroup() %>%
  left_join(range_df, by = "site_no") %>%
  mutate(
    start = pmax(start, min_date),
    end   = pmin(end,   max_date)
  ) %>%
  select(site_no, peak_id, peak_date, start, end) %>%
  mutate(date = map2(start, end, ~ seq(.x, .y, by = "day"))) %>%
  select(site_no, peak_id, peak_date, date) %>%
  unnest(date) %>%
  # resolve overlaps: a day claimed by two windows goes to the nearer peak
  group_by(site_no, date) %>%
  slice_min(abs(as.numeric(date - peak_date)), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(site_no, peak_id, date)

buffered_windows <- buffered_dates %>%
  mutate(peak_id = as.factor(as.character(peak_id)))

# ============================================================================
# ATTACH PEAK DATE AND FLOW
#
# The peak is taken directly from the declustered set rather than recomputed
# as the maximum within the buffered window. Recomputing was wrong: the
# window spans +/- search_radius days, so on a rising limb the window maximum
# can fall on a day belonging to a different, larger event. That relocated
# both the date and the magnitude of the retained peak, inflating the flows
# entering the GPD fits and shifting the antecedent window used for
# classification.
#
# peak_id is assigned exactly as in buffered_dates (row number within site,
# in date order) so the join keys correspond.
# ============================================================================
flows_in_windows <- all_flows %>%
  inner_join(buffered_windows, by = c("site_no", "date"))

peak_info <- independent_peaks %>%
  select(site_no, peak_date = date, peak_flow = q_daily_max) %>%
  group_by(site_no) %>%
  arrange(peak_date, .by_group = TRUE) %>%
  mutate(peak_id = as.factor(as.character(row_number()))) %>%
  ungroup() %>%
  left_join(site_thresholds, by = "site_no")   # attach candidate threshold

# ============================================================================
# TRUNCATE WINDOWS TO EVENT DURATION (75% of peak flow threshold)
# All independent peaks are retained — no top-n filtering applied here.
# ============================================================================
pot_windows <- flows_in_windows %>%
  left_join(peak_info, by = c("site_no", "peak_id")) %>%
  left_join(search_window, by = "site_no") %>%
  left_join(antecedent_runoff, by = c("site_no", "peak_date")) %>%
  group_by(site_no, peak_id) %>%
  mutate(
    original_window_size = n(),
    flow_threshold       = 0.75 * peak_flow,
    is_event             = q_daily_max >= flow_threshold,
    event_block          = cumsum(is_event != lag(is_event, default = first(is_event)))
  ) %>%
  filter(
    original_window_size <= 3 |
      (event_block == event_block[date == peak_date] & is_event)
  ) %>%
  ungroup() %>%
  select(site_no, date, peak_id, peak_date, peak_flow, candidate_threshold_cfs,
         search_radius, antecedent_runoff_mm, antecedent_days_used)

# ============================================================================
# SAVE OUTPUTS AND SUMMARIZE
# ============================================================================
saveRDS(pot_windows, output_file)
message(sprintf("\nSaved POT windows to %s", output_file))

saveRDS(site_record_length, record_length_file)
message(sprintf("Saved record length for %d sites to %s",
                nrow(site_record_length), record_length_file))

site_summary <- pot_windows %>%
  group_by(site_no) %>%
  summarize(
    n_peaks                 = n_distinct(peak_id),
    candidate_threshold_cfs = first(candidate_threshold_cfs),
    .groups                 = "drop"
  ) %>%
  left_join(site_record_length %>% select(site_no, eff_years), by = "site_no") %>%
  mutate(peaks_per_year = n_peaks / eff_years) %>%
  arrange(n_peaks)

print(site_summary, n = Inf)

# mean peaks per year, for comparison against Lang et al.'s mu > 2-3 guidance
cat("\nPeaks per year (mu): median",
    round(median(site_summary$peaks_per_year), 2),
    "| min", round(min(site_summary$peaks_per_year), 2),
    "| max", round(max(site_summary$peaks_per_year), 2), "\n")
