#!/usr/bin/env Rscript

# ============================================================================
# Script 6b: Peak Classification by Cell-Level WLDAS Fractions
# ============================================================================
# For each POT peak event, loads the per-cell WLDAS Parquets produced by
# script 6a, accumulates 8-day antecedent window totals per cell, and
# computes area-weighted fractions of cells meeting each Davenport et al.
# classification condition.
#
# Cell-level classification (applied to 8-day cumulative per-cell totals):
#   Rainfall    : rain > 10 mm AND melt < 5 mm
#   Snowmelt    : melt > 10 mm AND rain < 5 mm
#   Rain-on-Snow: rain > 5 mm  AND melt > 5 mm
#   Mixed       : none of the above
#
# Area-weighted fractions use exactextractr coverage_frac as the weight.
# Peak-level classification uses a 50% area threshold for dominant type.
#
# Memory strategy: parallelise over sites. Each worker loads only the cell
# Parquets relevant to one site's peak dates, computes fractions, writes a
# per-site Parquet to a scratch directory, and discards all data. A final
# sequential pass binds the per-site Parquets into a single output RDS.
#
# Output schema (one row per peak):
#   site_no         : character
#   peak_id         : factor
#   peak_date       : Date
#   peak_flow       : numeric
#   frac_rainfall   : numeric  area-weighted fraction of rainfall cells
#   frac_snowmelt   : numeric  area-weighted fraction of snowmelt cells
#   frac_ros        : numeric  area-weighted fraction of rain-on-snow cells
#   frac_mixed      : numeric  area-weighted fraction of mixed cells
#   dominant_type   : character  winning type or "mixed" if no type > 50%
#   dominant_frac   : numeric  fraction of dominant type
#   n_days_cell     : integer  days of cell data available (max 8)
#   complete_window : logical  TRUE if all 8 days present
#
# INPUTS:
#   - _data/pot_peaks_all.rds
#   - _data/WLDAS/processed/cell_daily_values/
#       YYYYMMDD_wldas_cells.parquet
#
# OUTPUTS:
#   - _data/pot_peaks_cell_classified.rds
#   - _data/cell_classified_by_site/  (intermediate per-site Parquets,
#                                      safe to delete after successful run)
# ============================================================================

library(tidyverse)
library(arrow)
library(furrr)
library(future)

# ============================================================================
# CONFIGURATION
# ============================================================================
peaks_file       <- "_data/pot_peaks_all.rds"
cell_parquet_dir <- "_data/WLDAS/processed/cell_daily_values"
site_scratch_dir <- "_data/WLDAS/processed/cell_classified_by_site"
output_file      <- "_data/pot_peaks_cell_classified.rds"

# Davenport et al. thresholds (mm, 8-day cumulative)
THRESH_RAIN_HIGH <- 10
THRESH_RAIN_LOW  <-  5
THRESH_MELT_HIGH <- 10
THRESH_MELT_LOW  <-  5
ANTECEDENT_DAYS  <-  7   # window = peak-7 through peak (8 days)

# Minimum area-weighted fraction for dominant type declaration
DOMINANT_THRESH  <- 0.50

n_workers        <- 16

dir.create(site_scratch_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD PEAK INVENTORY
# ============================================================================
message("Loading POT peaks...")
peaks <- readRDS(peaks_file) %>%
  mutate(peak_date = as.Date(peak_date)) %>%
  distinct(site_no, peak_id, peak_date, peak_flow) %>%
  mutate(GAGE_ID = site_no)

message(sprintf("Loaded %d peaks from %d sites",
                nrow(peaks), n_distinct(peaks$site_no)))

# ============================================================================
# IDENTIFY AVAILABLE CELL PARQUET DATES
# ============================================================================
message("Scanning cell parquet directory...")
cell_files <- list.files(cell_parquet_dir, pattern = "\\.parquet$",
                         full.names = TRUE)

if (length(cell_files) == 0)
  stop(sprintf("No cell parquet files found in %s", cell_parquet_dir))

cell_dates <- cell_files %>%
  basename() %>%
  str_extract("\\d{8}") %>%
  as.Date(format = "%Y%m%d")

message(sprintf("Cell data: %s to %s (%d dates)",
                min(cell_dates), max(cell_dates), length(cell_dates)))

# ============================================================================
# FLAG PEAKS OUTSIDE COVERAGE
# ============================================================================
peaks <- peaks %>%
  mutate(
    window_start   = peak_date - ANTECEDENT_DAYS,
    window_end     = peak_date,
    cell_available = window_start >= min(cell_dates) &
                     window_end   <= max(cell_dates)
  )

n_outside <- sum(!peaks$cell_available)
if (n_outside > 0)
  message(sprintf("Warning: %d peaks outside cell data coverage", n_outside))

classifiable <- peaks %>% filter(cell_available)
all_sites    <- unique(classifiable$site_no)

message(sprintf("Sites to process: %d", length(all_sites)))

# ============================================================================
# FUNCTION: Classify peaks for a single site
# ============================================================================
classify_site <- function(site_id, classifiable, cell_parquet_dir,
                          cell_dates, site_scratch_dir,
                          thresh_rain_high, thresh_rain_low,
                          thresh_melt_high, thresh_melt_low,
                          antecedent_days, dominant_thresh) {

  out_file <- file.path(site_scratch_dir,
                        sprintf("%s_cell_classified.parquet", site_id))

  # Skip if already processed (checkpoint behaviour)
  if (file.exists(out_file)) return(invisible(NULL))

  tryCatch({
    site_peaks <- classifiable %>%
      filter(site_no == site_id)

    # Unique dates needed for this site's peaks
    site_dates <- site_peaks %>%
      rowwise() %>%
      reframe(date = seq(window_start, window_end, by = "day")) %>%
      distinct(date) %>%
      filter(date %in% cell_dates) %>%
      pull(date)

    if (length(site_dates) == 0) return(invisible(NULL))

    # Load only this site's cells for the required dates
    cell_data <- map_dfr(site_dates, function(d) {
      f <- file.path(cell_parquet_dir,
                     sprintf("%s_wldas_cells.parquet", format(d, "%Y%m%d")))
      tryCatch(
        read_parquet(f, col_select = c("date", "GAGE_ID", "cell_lon", "cell_lat",
                                       "coverage_frac", "rainf_flux", "qsm_flux")) %>%
          filter(GAGE_ID == site_id) %>%
          mutate(date = as.Date(date)),
        error = function(e) NULL
      )
    })

    if (nrow(cell_data) == 0) return(invisible(NULL))

    # Expand peaks to window dates, join cell data
    peak_cell_windows <- site_peaks %>%
      select(site_no, GAGE_ID, peak_id, peak_date, peak_flow,
             window_start, window_end) %>%
      rowwise() %>%
      reframe(
        site_no   = site_no,
        GAGE_ID   = GAGE_ID,
        peak_id   = as.character(peak_id),
        peak_date = peak_date,
        peak_flow = peak_flow,
        date      = seq(window_start, window_end, by = "day")
      ) %>%
      left_join(cell_data, by = c("GAGE_ID", "date"))

    # 8-day cumulative totals per peak × cell
    cell_totals <- peak_cell_windows %>%
      group_by(site_no, peak_id, peak_date, peak_flow,
               cell_lon, cell_lat, coverage_frac) %>%
      summarise(
        rain_cumul     = sum(rainf_flux, na.rm = TRUE),
        melt_cumul     = sum(qsm_flux,   na.rm = TRUE),
        n_days_present = sum(!is.na(rainf_flux)),
        .groups        = "drop"
      ) %>%
      mutate(
        cell_type = case_when(
          rain_cumul > thresh_rain_high & melt_cumul < thresh_melt_low ~ "rainfall",
          melt_cumul > thresh_melt_high & rain_cumul < thresh_rain_low ~ "snowmelt",
          rain_cumul > thresh_rain_low  & melt_cumul > thresh_melt_low ~ "ros",
          TRUE                                                          ~ "mixed"
        )
      )

    # Area-weighted fractions per peak
    peak_fractions <- cell_totals %>%
      group_by(site_no, peak_id, peak_date, peak_flow) %>%
      summarise(
        total_weight  = sum(coverage_frac, na.rm = TRUE),
        frac_rainfall = sum(coverage_frac[cell_type == "rainfall"], na.rm = TRUE) /
                        total_weight,
        frac_snowmelt = sum(coverage_frac[cell_type == "snowmelt"], na.rm = TRUE) /
                        total_weight,
        frac_ros      = sum(coverage_frac[cell_type == "ros"],      na.rm = TRUE) /
                        total_weight,
        frac_mixed    = sum(coverage_frac[cell_type == "mixed"],    na.rm = TRUE) /
                        total_weight,
        n_days_cell   = as.integer(round(mean(n_days_present))),
        .groups       = "drop"
      ) %>%
      select(-total_weight) %>%
      mutate(
        complete_window = n_days_cell == (antecedent_days + 1L),
        dominant_type   = case_when(
          !complete_window                   ~ NA_character_,
          frac_rainfall > dominant_thresh    ~ "rainfall",
          frac_snowmelt > dominant_thresh    ~ "snowmelt",
          frac_ros      > dominant_thresh    ~ "rain_on_snow",
          TRUE                               ~ "mixed"
        ),
        dominant_frac = pmax(frac_rainfall, frac_snowmelt, frac_ros, frac_mixed,
                             na.rm = TRUE)
      )

    write_parquet(peak_fractions, out_file)
    invisible(NULL)

  }, error = function(e) {
    message(sprintf("ERROR site %s: %s", site_id, e$message))
    invisible(NULL)
  })
}

# ============================================================================
# PARALLEL EXECUTION OVER SITES
# ============================================================================
# Check which sites already have output (checkpoint behaviour)
already_done <- list.files(site_scratch_dir, pattern = "\\.parquet$") %>%
  str_remove("_cell_classified\\.parquet$")

sites_to_run <- setdiff(all_sites, already_done)
message(sprintf("\n%d sites to process (%d already done)",
                length(sites_to_run), length(already_done)))

plan(multicore, workers = n_workers)

future_walk(
  sites_to_run,
  function(site_id) {
    classify_site(
      site_id         = site_id,
      classifiable    = classifiable,
      cell_parquet_dir = cell_parquet_dir,
      cell_dates      = cell_dates,
      site_scratch_dir = site_scratch_dir,
      thresh_rain_high = THRESH_RAIN_HIGH,
      thresh_rain_low  = THRESH_RAIN_LOW,
      thresh_melt_high = THRESH_MELT_HIGH,
      thresh_melt_low  = THRESH_MELT_LOW,
      antecedent_days  = ANTECEDENT_DAYS,
      dominant_thresh  = DOMINANT_THRESH
    )
  },
  .progress = TRUE,
  .options  = furrr_options(seed = TRUE)
)

plan(sequential)

# ============================================================================
# BIND PER-SITE PARQUETS INTO FINAL OUTPUT
# ============================================================================
message("\nBinding per-site results...")

site_files <- list.files(site_scratch_dir, pattern = "\\.parquet$",
                         full.names = TRUE)

message(sprintf("Reading %d site Parquets...", length(site_files)))

classified <- map_dfr(site_files, read_parquet)

# Append peaks outside cell coverage
outside_coverage <- peaks %>%
  filter(!cell_available) %>%
  select(site_no, peak_id, peak_date, peak_flow) %>%
  mutate(
    peak_id         = as.character(peak_id),
    frac_rainfall   = NA_real_,
    frac_snowmelt   = NA_real_,
    frac_ros        = NA_real_,
    frac_mixed      = NA_real_,
    n_days_cell     = 0L,
    complete_window = FALSE,
    dominant_type   = NA_character_,
    dominant_frac   = NA_real_
  )

result <- bind_rows(classified, outside_coverage) %>%
  arrange(site_no, peak_date)

# ============================================================================
# SUMMARY
# ============================================================================
message("\n=== CELL CLASSIFICATION SUMMARY ===")
message(sprintf("Total peaks: %d", nrow(result)))
message(sprintf("  Complete windows: %d (%.1f%%)",
                sum(result$complete_window, na.rm = TRUE),
                100 * mean(result$complete_window, na.rm = TRUE)))

message("\nDominant type breakdown (complete windows):")
result %>%
  filter(complete_window) %>%
  count(dominant_type) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(desc(n)) %>%
  print()

message("\nMean area-weighted fractions (complete windows):")
result %>%
  filter(complete_window) %>%
  summarise(across(starts_with("frac_"),
                   ~ round(mean(.x, na.rm = TRUE), 3))) %>%
  print()

# ============================================================================
# SAVE
# ============================================================================
saveRDS(result, output_file)
message(sprintf("\nSaved %d peaks to %s", nrow(result), output_file))
message(sprintf("Per-site scratch files retained in %s", site_scratch_dir))
message("These can be deleted once the output RDS is verified.")