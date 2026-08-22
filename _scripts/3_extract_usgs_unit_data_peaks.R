#!/usr/bin/env Rscript

# ============================================================================
# extract_usgs_unit_data_peaks.R
#
# Find all independent POT peaks for selected unit-data stations.
# Implements independence criteria from Lang et al. 1999.
# All independent peaks are retained for later flood type classification
# and per-type GPD fitting (following Yu et al. 2022).
#
# Key addition vs. previous version: the per-site PoT threshold
# (q_threshold quantile of daily max flow over the full record) is saved
# as `pot_threshold_cfs` in the output. This value is used as the fixed
# GPD location parameter mu in downstream fitting, following Yu et al. 2022.
#
# INPUTS:
#   - selected_unit_stations.csv   : site_no column of selected stations
#   - /path/to/usgs/unit/<site_no>/<year>.rds : unit discharge files
#   - conterm_basinid.txt          : GAGES-II metadata (drainage area)
#
# OUTPUT:
#   - pot_peaks_all.rds : one row per site-day inside each peak window,
#                         with columns:
#                           site_no, date, peak_id, peak_date, peak_flow,
#                           pot_threshold_cfs
# ============================================================================

library(tidyverse)

# ============================================================================
# CONFIGURATION
# ============================================================================
selected_stations_csv <- "_data/selected_unit_stations.csv"
usgs_unit_base_dir    <- "/davenport-scratch/DATA/USGS/raw/Discharge/Unit"
watershed_meta        <- "/davenport-scratch/DATA/GAGES-II/basinchar_and_report_sept_2011/spreadsheets-in-csv-format/conterm_basinid.txt"
output_file           <- "_data/pot_peaks_all.rds"

q_threshold <- 0.5   # quantile of daily max flow used as PoT threshold / GPD mu

# ============================================================================
# FUNCTION: Independence criteria from Lang et al. (1999)
#
# Applies two criteria to decluster candidate peaks for a single watershed:
#   1. Minimum time separation: theta > 5 + ln(0.3861 * A), A in km²
#      (Lang et al. state the criterion with A in miles²; 0.3861 converts km²)
#   2. Inter-peak recession: flow between peaks must drop below 75% of the
#      lesser of the two peak flows
#
# Candidates are considered in descending order of magnitude so that when two
# peaks are dependent, the LARGER is retained. A peak is accepted only if it
# is independent of every peak already accepted.
# ============================================================================
find_independent_peaks <- function(peak_df, all_flows_df, watershed_drain_sqkm) {

  if (nrow(peak_df) < 2) return(peak_df)

  min_sep_days <- 5 + log(0.3861 * watershed_drain_sqkm)

  # pull to plain vectors once; avoid dplyr overhead in the inner loop
  flow_dates <- all_flows_df$date
  flow_vals  <- all_flows_df$q_daily_max
  ord        <- order(flow_dates)
  flow_dates <- flow_dates[ord]
  flow_vals  <- flow_vals[ord]

  pk_date <- peak_df$date
  pk_flow <- peak_df$q_daily_max

  order_by_size <- order(pk_flow, decreasing = TRUE)
  kept <- integer(0)

  for (i in order_by_size) {

    is_independent <- TRUE

    for (j in kept) {

      time_diff <- abs(as.numeric(pk_date[i] - pk_date[j]))

      # criterion 1: cheap, so check first
      if (time_diff < min_sep_days) {
        is_independent <- FALSE
        break
      }

      # criterion 2: index into the sorted vectors, no table scan
      d1 <- min(pk_date[i], pk_date[j])
      d2 <- max(pk_date[i], pk_date[j])
      lo <- findInterval(d1, flow_dates) + 1L
      hi <- findInterval(d2 - 1, flow_dates)

      if (hi >= lo) {
        min_flow_between <- min(flow_vals[lo:hi], na.rm = TRUE)
      } else {
        min_flow_between <- Inf   # no days between; treat as dependent
      }

      if (!isTRUE(min_flow_between < 0.75 * min(pk_flow[i], pk_flow[j]))) {
        is_independent <- FALSE
        break
      }
    }

    if (is_independent) kept <- c(kept, i)
  }

  peak_df[sort(kept), ]
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
# COMPUTE PER-SITE POT THRESHOLD
#
# The q_threshold quantile of daily max flow over the full record becomes:
#   (a) the minimum flow a candidate peak must exceed, and
#   (b) the GPD location parameter mu in downstream fitting (Yu et al. 2022).
#
# Saved separately so it can be joined into the output and used by app.R
# without needing access to the raw daily flow record.
# ============================================================================
site_thresholds <- all_flows %>%
  group_by(site_no) %>%
  summarize(
    pot_threshold_cfs = quantile(q_daily_max, q_threshold, na.rm = TRUE),
    .groups = "drop"
  )

message(sprintf(
  "Computed PoT thresholds for %d sites (q=%.2f of daily max flow)",
  nrow(site_thresholds), q_threshold
))

# ============================================================================
# MINIMUM SEPARATION WINDOW (Lang et al. 1999)
# theta = 5 + ln(0.3861 * A), A in km²
# ============================================================================
search_window <- basin_meta %>%
  filter(site_no %in% unique(all_flows$site_no)) %>%
  transmute(
    site_no,
    search_radius = as.integer(pmax(1L, ceiling(5 + log(0.3861 * drain_sqkm))))
  )

# ============================================================================
# IDENTIFY CANDIDATE PEAKS OVER THRESHOLD
# ============================================================================
candidate_peaks <- all_flows %>%
  left_join(site_thresholds, by = "site_no") %>%
  group_by(site_no) %>%
  filter(q_daily_max > pot_threshold_cfs) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(
    time_gap = as.numeric(date - lag(date), units = "days"),
    event_id = cumsum(is.na(time_gap) | time_gap > 1)
  ) %>%
  group_by(site_no, event_id) %>%
  slice_max(order_by = q_daily_max, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-time_gap, -event_id, -pot_threshold_cfs) %>%
  arrange(site_no, date)

# ============================================================================
# DECLUSTER PEAKS USING LANG ET AL. INDEPENDENCE CRITERIA
# ============================================================================
independent_peaks <- candidate_peaks %>%
  group_by(site_no, drain_sqkm) %>%
  group_modify(~ find_independent_peaks(
    .x,
    all_flows_df         = filter(all_flows, site_no == .y$site_no),
    watershed_drain_sqkm = .y$drain_sqkm
  )) %>%
  ungroup()

message(sprintf("Found %d independent peaks across %d sites",
                nrow(independent_peaks), n_distinct(independent_peaks$site_no)))

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
# IDENTIFY PEAK DATE AND FLOW WITHIN EACH WINDOW
# ============================================================================
flows_in_windows <- all_flows %>%
  inner_join(buffered_windows, by = c("site_no", "date"))

peak_info <- flows_in_windows %>%
  group_by(site_no, peak_id) %>%
  slice_max(order_by = q_daily_max, n = 1, with_ties = FALSE) %>%
  summarize(
    peak_date = date,
    peak_flow = q_daily_max,
    .groups   = "drop"
  ) %>%
  left_join(site_thresholds, by = "site_no")   # attach mu for downstream use

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
  select(site_no, date, peak_id, peak_date, peak_flow, pot_threshold_cfs,
         search_radius, antecedent_runoff_mm, antecedent_days_used)

# ============================================================================
# SAVE OUTPUT AND SUMMARIZE
# ============================================================================
saveRDS(pot_windows, output_file)
message(sprintf("\nSaved POT windows to %s", output_file))

pot_windows %>%
  group_by(site_no) %>%
  summarize(
    n_peaks           = n_distinct(peak_id),
    pot_threshold_cfs = first(pot_threshold_cfs),
    .groups           = "drop"
  ) %>%
  arrange(n_peaks) %>%
  print(n = Inf)