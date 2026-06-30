#!/usr/bin/env Rscript

# ============================================================================
# Classify POT peaks by flood type using WLDAS antecedent climate data
# ============================================================================
# Follows the classification criteria of Davenport et al. (2020):
#   - Rainfall-driven:  total rain > 10 mm AND snowmelt < 5 mm
#   - Snowmelt-driven:  snowmelt > 10 mm AND total rain < 5 mm
#   - Rain-on-snow:     both rain > 5 mm AND snowmelt > 5 mm
#   - Unclassified:     does not meet any of the above criteria
#
# classification_mode controls the spatial support of the forcing variables:
#   "mean"     : basin-mean WLDAS values (Davenport et al. original approach)
#   "max_cell" : maximum cell value within watershed (tests whether spatial
#                aggregation suppresses classification signal, particularly
#                for topographically concentrated snowmelt)
#
# Antecedent window: 8 days ending on peak date (days -7 through 0),
# matching Davenport et al. Fig. 2b.
#
# Snowmelt is computed from qsm_flux (direct WLDAS snowmelt flux, mm/day)
# summed over the window, rather than net SWE decrease. This avoids masking
# of melt by concurrent accumulation mid-window.
#
# Peak selection follows Yu et al. (2022):
#   - Single sample: largest m events where m = max(25, n_water_years),
#                    at most one peak per water year.
#   - Per type:      largest n events where n = max(20, n_water_years),
#                    at most one peak per water year.
#   - Types with fewer than 20 events are excluded (flood_type -> "unclassified").
#   Selection flags are stored in in_single_sample and in_type_sample columns.
#
# INPUTS:
#   - pot_peaks_all.rds             : output of extract_usgs_unit_data_peaks.R
#                                     must contain pot_threshold_cfs column
#   - YYYYMMDD_wldas.parquet        : daily WLDAS basin averages
#   - conterm_basinid.txt           : GAGES-II site identity / location
#   - conterm_topo.txt              : GAGES-II topographic attributes
#   - conterm_hydromod_dams.txt     : GAGES-II dam / hydrologic modification
#   - conterm_bas_classif.txt       : GAGES-II basin classification
#
# OUTPUTS:
#   - pot_peaks_classified_{mode}.rds : peaks with flood_type, flood_type_raw,
#                                 pot_threshold_cfs, water_year,
#                                 in_single_sample, in_type_sample
#                                 flood_type_raw: threshold-based classification
#                                 before Yu et al. sample size filter;
#                                 flood_type: final classification after filter
#                                 (sample size demotions appear as
#                                 flood_type_raw != "unclassified" but
#                                 flood_type == "unclassified")
#   - site_metadata.rds         : one row per site with GAGES-II attributes
# ============================================================================

library(tidyverse)
library(arrow)

# ============================================================================
# CONFIGURATION
# ============================================================================
pot_peaks_file  <- "_data/pot_peaks_all.rds"
wldas_dir       <- "_data/WLDAS/processed/daily_by_date"

gagesii_dir     <- "/davenport-scratch/DATA/GAGES-II/basinchar_and_report_sept_2011/spreadsheets-in-csv-format"
meta_basinid    <- file.path(gagesii_dir, "conterm_basinid.txt")
meta_topo       <- file.path(gagesii_dir, "conterm_topo.txt")
meta_dams       <- file.path(gagesii_dir, "conterm_hydromod_dams.txt")
meta_classif    <- file.path(gagesii_dir, "conterm_bas_classif.txt")

# Classification mode: controls spatial support of forcing variables.
#   "mean"     : basin-mean (original Davenport et al. approach)
#   "max_cell" : maximum cell within watershed (tests spatial aggregation hypothesis)
classification_mode <- "max_cell"

if (!classification_mode %in% c("mean", "max_cell")) {
  stop("classification_mode must be 'mean' or 'max_cell'")
}

# WLDAS columns to use for rain, snow, and melt based on classification_mode
rain_col <- if (classification_mode == "max_cell") "rainf_flux_max" else "rainf_flux"
snow_col <- if (classification_mode == "max_cell") "snowf_flux_max" else "snowf_flux"
melt_col <- if (classification_mode == "max_cell") "qsm_flux_max"   else "qsm_flux"

output_peaks    <- sprintf("_data/pot_peaks_classified_%s.rds", classification_mode)
output_metadata <- "_data/site_metadata.rds"

# Davenport et al. 2020 thresholds (mm over 8-day window)
thresh_rain_high    <- 10
thresh_rain_low     <-  5
thresh_melt_high    <- 10
thresh_melt_low     <-  5
antecedent_days     <-  7

# Yu et al. 2022 sample size thresholds
min_single_sample   <- 25
min_type_sample     <- 20

# ============================================================================
# LOAD PEAKS AND RESOLVE GAGE_ID
# ============================================================================
message("Loading POT peaks...")
peaks <- readRDS(pot_peaks_file) %>%
  mutate(peak_date = as.Date(peak_date))

message(sprintf("Loaded %d peak-day rows covering %d sites",
                nrow(peaks), n_distinct(peaks$site_no)))

if (!"pot_threshold_cfs" %in% names(peaks)) {
  stop("pot_threshold_cfs column not found in pot_peaks_all.rds. ",
       "Re-run extract_usgs_unit_data_peaks.R to regenerate.")
}

peak_events <- peaks %>%
  distinct(site_no, peak_id, peak_date, peak_flow, pot_threshold_cfs)

message(sprintf("Total independent peaks to classify: %d", nrow(peak_events)))

site_ids <- unique(peak_events$site_no)

# GAGE_ID lookup for WLDAS joins
basin_meta <- read_csv(meta_basinid, show_col_types = FALSE) %>%
  filter(STATE == "CO") %>%
  transmute(site_no = as.character(STAID), GAGE_ID = as.character(STAID)) %>%
  distinct()

# ============================================================================
# BUILD AND SAVE SITE METADATA
# ============================================================================
message("Loading GAGES-II metadata tables...")

# Helper: read a GAGES-II table, standardise site_no, filter to our sites
read_gagesii <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(site_no = as.character(STAID)) %>%
    filter(site_no %in% site_ids)
}

tbl_basinid <- read_gagesii(meta_basinid) %>%
  transmute(
    site_no,
    station_name = STANAME,
    state        = STATE,
    lat_deg      = LAT_GAGE,
    lon_deg      = LNG_GAGE,
    drain_sqkm   = DRAIN_SQKM,
    huc02        = HUC02
  )

tbl_topo <- read_gagesii(meta_topo) %>%
  transmute(
    site_no,
    elev_mean_m  = ELEV_MEAN_M_BASIN,
    slope_pct    = SLOPE_PCT
  )

tbl_dams <- read_gagesii(meta_dams) %>%
  transmute(
    site_no,
    ndams_2009        = NDAMS_2009,
    stor_norm_mcm     = STOR_NOR_2009,      # normal storage upstream dams (mcm)
    dist_nearest_dam_km = RAW_DIS_NEAREST_DAM
  )

tbl_classif <- read_gagesii(meta_classif) %>%
  transmute(
    site_no,
    class             = CLASS,              # "Ref" or "Non-ref"
    hydro_disturb_idx = HYDRO_DISTURB_INDX
  )

site_metadata <- tbl_basinid %>%
  left_join(tbl_topo,    by = "site_no") %>%
  left_join(tbl_dams,    by = "site_no") %>%
  left_join(tbl_classif, by = "site_no") %>%
  arrange(site_no)

message(sprintf("Site metadata assembled for %d sites:", nrow(site_metadata)))
message(sprintf("  Reference (CLASS == 'Ref'):     %d",
                sum(site_metadata$class == "Ref",     na.rm = TRUE)))
message(sprintf("  Non-reference (CLASS != 'Ref'): %d",
                sum(site_metadata$class != "Ref",     na.rm = TRUE)))
message(sprintf("  Sites with >= 1 dam upstream:   %d",
                sum(site_metadata$ndams_2009 > 0,     na.rm = TRUE)))

saveRDS(site_metadata, output_metadata)
message(sprintf("Saved site metadata to %s\n", output_metadata))

# ============================================================================
# DETERMINE WLDAS DATE COVERAGE
# ============================================================================
message("Checking WLDAS date coverage...")
parquet_files <- list.files(wldas_dir, pattern = "\\.parquet$", full.names = TRUE)

if (length(parquet_files) == 0) {
  stop(sprintf("No parquet files found in %s", wldas_dir))
}

wldas_dates <- parquet_files %>%
  basename() %>%
  str_extract("\\d{8}") %>%
  as.Date(format = "%Y%m%d") %>%
  sort()

message(sprintf("WLDAS coverage: %s to %s (%d days)",
                min(wldas_dates), max(wldas_dates), length(wldas_dates)))

peak_events <- peak_events %>%
  mutate(
    window_start    = peak_date - antecedent_days,
    window_end      = peak_date,
    wldas_available = window_start >= min(wldas_dates) & window_end <= max(wldas_dates)
  )

n_outside <- sum(!peak_events$wldas_available)
if (n_outside > 0) {
  message(sprintf("Warning: %d peaks have antecedent windows outside WLDAS coverage and will be unclassified",
                  n_outside))
}

# ============================================================================
# IDENTIFY REQUIRED WLDAS DATES AND LOAD THEM
# ============================================================================
classifiable_peaks <- peak_events %>% filter(wldas_available)

required_dates <- classifiable_peaks %>%
  rowwise() %>%
  reframe(date = seq(window_start, window_end, by = "day")) %>%
  distinct(date) %>%
  filter(date %in% wldas_dates) %>%
  pull(date)

message(sprintf("Loading WLDAS data for %d unique dates (mode: %s)...",
                length(required_dates), classification_mode))

wldas_data <- map_dfr(required_dates, function(d) {
  f <- file.path(wldas_dir, sprintf("%s_wldas.parquet", format(d, "%Y%m%d")))
  tryCatch(
    read_parquet(f, col_select = all_of(c("date", "GAGE_ID",
                                          rain_col, snow_col, melt_col))) %>%
      mutate(date = as.Date(date), GAGE_ID = as.character(GAGE_ID)),
    error = function(e) {
      warning(sprintf("Could not read %s: %s", basename(f), e$message))
      NULL
    }
  )
})

message(sprintf("Loaded %d WLDAS site-day rows", nrow(wldas_data)))

# ============================================================================
# COMPUTE ANTECEDENT WINDOW TOTALS PER PEAK
# ============================================================================
message("Computing antecedent window totals...")

peak_windows <- classifiable_peaks %>%
  left_join(basin_meta, by = "site_no") %>%
  rowwise() %>%
  reframe(
    GAGE_ID           = GAGE_ID,
    site_no           = site_no,
    peak_id           = peak_id,
    peak_date         = peak_date,
    peak_flow         = peak_flow,
    pot_threshold_cfs = pot_threshold_cfs,
    date              = seq(window_start, window_end, by = "day")
  ) %>%
  left_join(wldas_data, by = c("GAGE_ID", "date"))

window_totals <- peak_windows %>%
  group_by(site_no, peak_id, peak_date, peak_flow, pot_threshold_cfs) %>%
  summarize(
    total_rain_mm = sum(.data[[rain_col]], na.rm = TRUE),
    total_snow_mm = sum(.data[[snow_col]], na.rm = TRUE),
    total_melt_mm = sum(.data[[melt_col]], na.rm = TRUE),
    n_days_wldas  = sum(!is.na(.data[[rain_col]])),
    .groups = "drop"
  ) %>%
  mutate(complete_window = n_days_wldas == (antecedent_days + 1))

n_incomplete <- sum(!window_totals$complete_window)
if (n_incomplete > 0) {
  message(sprintf("Warning: %d peaks have incomplete WLDAS windows (< 8 days of data)",
                  n_incomplete))
}

# ============================================================================
# APPLY DAVENPORT ET AL. CLASSIFICATION THRESHOLDS
# ============================================================================
message("Classifying flood types...")

window_totals <- window_totals %>%
  mutate(
    flood_type = case_when(
      !complete_window                                                      ~ "unclassified",
      total_rain_mm > thresh_rain_high & total_melt_mm < thresh_melt_low   ~ "rainfall",
      total_melt_mm > thresh_melt_high & total_rain_mm < thresh_rain_low   ~ "snowmelt",
      total_rain_mm > thresh_rain_low  & total_melt_mm > thresh_melt_low   ~ "rain_on_snow",
      TRUE                                                                  ~ "unclassified"
    ),
    # Preserve pre-filter classification to distinguish threshold failures
    # from sample size demotions in downstream diagnostics. flood_type will
    # be set to "unclassified" for types with < min_type_sample events;
    # flood_type_raw retains the original threshold-based assignment.
    flood_type_raw = flood_type
  )

# ============================================================================
# COMBINE WITH UNCLASSIFIABLE PEAKS (OUTSIDE WLDAS COVERAGE)
# ============================================================================
outside_coverage <- peak_events %>%
  filter(!wldas_available) %>%
  select(site_no, peak_id, peak_date, peak_flow, pot_threshold_cfs) %>%
  mutate(
    total_rain_mm   = NA_real_,
    total_snow_mm   = NA_real_,
    total_melt_mm   = NA_real_,
    n_days_wldas    = 0L,
    complete_window = FALSE,
    flood_type      = "unclassified",
    flood_type_raw  = "unclassified"
  )

classified_peaks <- bind_rows(window_totals, outside_coverage) %>%
  arrange(site_no, peak_date)

# ============================================================================
# VERIFY pot_threshold_cfs IS PRESENT AND COMPLETE
# ============================================================================
n_missing_thr <- sum(is.na(classified_peaks$pot_threshold_cfs))
if (n_missing_thr > 0) {
  warning(sprintf("%d rows are missing pot_threshold_cfs", n_missing_thr))
} else {
  message("pot_threshold_cfs verified present for all rows.")
}

# ============================================================================
# YU ET AL. (2022) TOP-N PEAK SELECTION
# ============================================================================
message("Applying Yu et al. (2022) top-n peak selection...")

water_year <- function(d) {
  yr <- as.integer(format(d, "%Y"))
  mo <- as.integer(format(d, "%m"))
  ifelse(mo >= 10L, yr + 1L, yr)
}

classified_peaks <- classified_peaks %>%
  mutate(water_year = water_year(peak_date))

single_sample_ids <- classified_peaks %>%
  filter(flood_type != "unclassified") %>%
  group_by(site_no, water_year) %>%
  slice_max(peak_flow, n = 1, with_ties = FALSE) %>%
  group_by(site_no) %>%
  mutate(m    = pmax(min_single_sample, n_distinct(water_year)),
         rank = row_number(desc(peak_flow))) %>%
  filter(rank <= m) %>%
  ungroup() %>%
  select(site_no, peak_id) %>%
  mutate(in_single_sample = TRUE)

type_sample_ids <- classified_peaks %>%
  filter(flood_type != "unclassified") %>%
  group_by(site_no, flood_type, water_year) %>%
  slice_max(peak_flow, n = 1, with_ties = FALSE) %>%
  group_by(site_no, flood_type) %>%
  mutate(n    = pmax(min_type_sample, n_distinct(water_year)),
         rank = row_number(desc(peak_flow))) %>%
  filter(rank <= n) %>%
  filter(n() >= min_type_sample) %>%
  ungroup() %>%
  select(site_no, peak_id, flood_type) %>%
  mutate(in_type_sample = TRUE)

classified_peaks <- classified_peaks %>%
  left_join(single_sample_ids, by = c("site_no", "peak_id")) %>%
  left_join(type_sample_ids %>% select(site_no, peak_id, in_type_sample),
            by = c("site_no", "peak_id")) %>%
  mutate(
    in_single_sample = replace_na(in_single_sample, FALSE),
    in_type_sample   = replace_na(in_type_sample,   FALSE)
  )

retained_site_types <- type_sample_ids %>% distinct(site_no, flood_type)

classified_peaks <- classified_peaks %>%
  mutate(flood_type = if_else(
    flood_type != "unclassified" &
      !paste(site_no, flood_type) %in%
        paste(retained_site_types$site_no, retained_site_types$flood_type),
    "unclassified",
    flood_type
  ))

# ============================================================================
# SUMMARY
# ============================================================================
type_summary <- classified_peaks %>%
  count(flood_type) %>%
  mutate(pct = round(100 * n / sum(n), 1))

message("\nClassification summary (all peaks, post-selection):")
print(type_summary)

selection_summary <- classified_peaks %>%
  group_by(site_no) %>%
  summarise(
    n_single_sample = sum(in_single_sample),
    n_rainfall      = sum(flood_type == "rainfall"     & in_type_sample),
    n_snowmelt      = sum(flood_type == "snowmelt"     & in_type_sample),
    n_rain_on_snow  = sum(flood_type == "rain_on_snow" & in_type_sample),
    .groups = "drop"
  )

message(sprintf("\nSites with >= %d single-sample peaks: %d",
                min_single_sample,
                sum(selection_summary$n_single_sample >= min_single_sample)))
message(sprintf("Sites with >= %d rainfall peaks:      %d",
                min_type_sample,
                sum(selection_summary$n_rainfall >= min_type_sample)))
message(sprintf("Sites with >= %d snowmelt peaks:      %d",
                min_type_sample,
                sum(selection_summary$n_snowmelt >= min_type_sample)))
message(sprintf("Sites with >= %d rain_on_snow peaks:  %d",
                min_type_sample,
                sum(selection_summary$n_rain_on_snow >= min_type_sample)))

# ============================================================================
# SAVE PEAKS OUTPUT
# ============================================================================
saveRDS(classified_peaks, output_peaks)
message(sprintf("\nSaved %d peaks to %s", nrow(classified_peaks), output_peaks))
message(sprintf("  classification_mode      : %s", classification_mode))
message(sprintf("  in_type_sample == TRUE   -> Yu et al. per-type FFA peaks (n >= %d per type)",
                min_type_sample))
message(sprintf("  in_single_sample == TRUE -> Yu et al. single-sample peaks (m >= %d)",
                min_single_sample))