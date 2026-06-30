#!/usr/bin/env Rscript

# ============================================================================
# WLDAS Daily to GAGES-II Watershed Extraction
# ============================================================================
# Extract basin-averaged WLDAS daily variables for Colorado GAGES-II watersheds
# Output: One Parquet file per date with all watersheds
#
# KEY FEATURES:
#   - Loops through NetCDF files (reads each file ONCE)
#   - Extracts to ALL watersheds at once using exactextractr
#   - Saves by date as Parquet (efficient storage & querying)
#   - Proper NetCDF coordinate system handling
#   - Loads multiple shapefiles (reference + non-reference by ecoregion)
#   - Detailed variable metadata with unit conversions
#
# SPATIAL STATISTICS (added):
#   Spatial statistics are computed only for rainfall, snowfall, and snowmelt
#   (SPATIAL_STATS_VARS), where spatial heterogeneity is scientifically
#   relevant for flood classification. All other variables get basin mean only.
#
#   _max     : maximum cell value within the watershed (post unit-conversion)
#   _wet_frac: fraction of cells exceeding WETTED_THRESH_RAW (in native units,
#              i.e. before conversion). Only computed for flux variables
#              (rainf, snowf, qsm).
#   _sd      : population standard deviation of cell values (post conversion)
#   _cv      : coefficient of variation = sd / mean. Set to NA when mean < 
#              CV_MEAN_MIN_MM to avoid instability near zero.
#
#   These are appended as additional columns: e.g. rainf_flux_max,
#   rainf_flux_wet_frac, rainf_flux_sd, rainf_flux_cv.
#
#   The wetted fraction threshold is defined in native flux units (kg/m2/s)
#   as WETTED_THRESH_RAW = 1mm/day / 86400 ≈ 1.157e-5 kg/m2/s, so the 
#   comparison is made before unit conversion. This is equivalent to asking
#   "what fraction of cells received > 1 mm/day?".
# ============================================================================

library(tidyverse)
library(sf)
library(terra)
library(exactextractr)
library(arrow)
library(lubridate)
library(furrr)
library(future)

# ============================================================================
# FILE PATHS AND CONFIGURATION
# ============================================================================
watershed_dataset <- "GAGES-II"
watershed_file <- "/davenport-scratch/DATA/GAGES-II/boundaries_shapefiles_by_aggeco/boundaries-shapefiles-by-aggeco/bas_ref_all.shp"
watershed_meta <- "/davenport-scratch/DATA/GAGES-II/basinchar_and_report_sept_2011/spreadsheets-in-csv-format/conterm_basinid.txt"
wldas_dir <- "_data/WLDAS/raw"
base_output_dir <- "_data/WLDAS/processed"
parquet_output_dir <- file.path(base_output_dir, "daily_by_date")

# ============================================================================
# SPATIAL STATISTICS CONFIGURATION
# ============================================================================
# Wetted fraction threshold in native WLDAS flux units (kg/m2/s).
# 1 mm/day = 1 kg/m2/day / 86400 s/day ≈ 1.157e-5 kg/m2/s.
# The comparison is made on raw (pre-conversion) values so the threshold
# is physically interpretable regardless of which variable is being tested.
WETTED_THRESH_RAW <- 1 / 86400   # kg/m2/s  (≈ 1 mm/day equivalent)

# Variables for which spatial statistics (max, sd, cv, wet_frac) are computed.
# Restricted to rainfall, snowfall, and snowmelt — the variables relevant to
# the spatial aggregation / flood classification hypothesis. All other variables
# get basin mean only, which reduces extraction time significantly.
SPATIAL_STATS_VARS <- c("rainf_flux", "snowf_flux", "qsm_flux")

# Variables for which wetted fraction is meaningful (must be subset of
# SPATIAL_STATS_VARS; flux variables only).
WETTED_FRAC_VARS <- c("rainf_flux", "snowf_flux", "qsm_flux")

# CV is set to NA when the basin-mean (post-conversion, mm) is below this
# threshold, to avoid numerically unstable CV values near zero.
CV_MEAN_MIN_MM <- 0.5  # mm (post unit-conversion)

# Command line arguments for chunking (optional)
# Usage: Rscript script.R [chunk_index] [num_chunks]
args <- commandArgs(trailingOnly = TRUE)
chunk_index <- ifelse(length(args) >= 1, as.integer(args[1]), 1)
num_chunks  <- ifelse(length(args) >= 2, as.integer(args[2]), 1)

# ============================================================================
# SETUP LOGGING
# ============================================================================
log_file <- file.path(base_output_dir, "extraction_LOG.txt")
dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(parquet_output_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  msg <- paste(..., sep = " ")
  cat(sprintf("[%s] %s\n", Sys.time(), msg))
  cat(sprintf("[%s] %s\n", Sys.time(), msg), file = log_file, append = TRUE)
}

log_msg("=== WLDAS Extraction Started ===")
log_msg(sprintf("Chunk: %d of %d", chunk_index, num_chunks))
log_msg(sprintf("Output directory: %s", parquet_output_dir))
log_msg(sprintf("Wetted fraction threshold: %.2e kg/m2/s (= 1 mm/day) for vars: %s",
                WETTED_THRESH_RAW, paste(WETTED_FRAC_VARS, collapse = ", ")))
log_msg(sprintf("Spatial statistics (max, sd, cv) computed for: %s",
                paste(SPATIAL_STATS_VARS, collapse = ", ")))
log_msg(sprintf("CV suppressed when basin mean < %.1f mm", CV_MEAN_MIN_MM))

# ============================================================================
# LOAD WATERSHEDS AND FILTER TO COLORADO
# ============================================================================
log_msg("Loading watershed metadata...")
basin_meta <- read_csv(watershed_meta, show_col_types = FALSE) %>%
  filter(STATE == "CO") %>%
  select(GAGE_ID = STAID) %>%
  distinct()

log_msg(sprintf("Found %d Colorado watersheds in metadata", nrow(basin_meta)))

log_msg("Loading watershed boundaries from all GAGES-II shapefiles...")

watershed_dir <- dirname(watershed_file)
shapefile_list <- c(
  "bas_ref_all.shp",
  "bas_nonref_WestMnts.shp",
  "bas_nonref_CntlPlains.shp",
  "bas_nonref_WestPlains.shp",
  "bas_nonref_WestXeric.shp"
)

all_watersheds <- map_df(
  shapefile_list,
  function(shp_name) {
    shp_path <- file.path(watershed_dir, shp_name)
    if (!file.exists(shp_path)) {
      log_msg(sprintf("  Warning: %s not found, skipping", shp_name))
      return(NULL)
    }
    log_msg(sprintf("  Loading %s...", shp_name))
    st_read(shp_path, quiet = TRUE) %>%
      st_transform("EPSG:4326")
  }
)

log_msg(sprintf("Loaded %d total watersheds from all shapefiles", nrow(all_watersheds)))

watersheds <- all_watersheds %>%
  inner_join(basin_meta, by = "GAGE_ID") %>%
  arrange(GAGE_ID) %>%
  mutate(chunk = ntile(row_number(), num_chunks)) %>%
  filter(chunk == chunk_index)

log_msg(sprintf("Processing %d Colorado watersheds (chunk %d of %d)",
                nrow(watersheds), chunk_index, num_chunks))

# ============================================================================
# DEFINE WLDAS VARIABLES
# ============================================================================
wldas_vars <- tibble(
  wldas_name = c(
    "Rainf_tavg", "Snowf_tavg", "SWE_tavg", "Qsm_tavg", "Tair_f_tavg",
    "SoilMoi00_10cm_tavg", "SoilMoi10_40cm_tavg", "SoilMoi40_100cm_tavg",
    "Qsb_tavg", "Qs_tavg", "Evap_tavg"
  ),
  output_name = c(
    "rainf_flux", "snowf_flux", "swe", "qsm_flux", "tair",
    "soilm_0_10cm", "soilm_10_40cm", "soilm_40_100cm",
    "qsb_flux", "qs_flux", "evap_flux"
  ),
  units = c(
    "kg/m2/s", "kg/m2/s", "kg/m2", "kg/m2/s", "K",
    "m3/m3", "m3/m3", "m3/m3",
    "kg/m2/s", "kg/m2/s", "kg/m2/s"
  ),
  conversion_factor = c(
    86400, 86400, 1, 86400, NA,
    1, 1, 1,
    86400, 86400, 86400
  )
)

log_msg(sprintf("Extracting %d variables from WLDAS:", nrow(wldas_vars)))
log_msg(paste(wldas_vars$output_name, collapse = ", "))

# ============================================================================
# CALCULATE WLDAS EXTENT FROM KNOWN PARAMETERS
# ============================================================================
get_wldas_extent <- function() {
  log_msg("Calculating WLDAS extent from known parameters...")
  resolution <- 0.01
  xmin <- -125.0
  ymin <- 25.0
  xmax <- xmin + (3591 * resolution)
  ymax <- ymin + (2787 * resolution)
  extent <- ext(xmin, xmax, ymin, ymax)
  log_msg(sprintf("WLDAS extent: lon=[%.2f, %.2f], lat=[%.2f, %.2f]",
                  xmin, xmax, ymin, ymax))
  return(extent)
}

wldas_extent <- get_wldas_extent()

# ============================================================================
# FUNCTION: Load and process single NetCDF file
# ============================================================================
load_wldas_nc <- function(nc_file, wldas_extent) {
  tryCatch({
    # Load raster
    r <- terra::rast(nc_file)

    # Force known extent and CRS — WLDAS NetCDF metadata can be unreliable
    terra::ext(r) <- wldas_extent
    terra::crs(r) <- "EPSG:4326"

    # Subset to requested variables only
    available_vars <- names(r)
    vars_to_use    <- wldas_vars$wldas_name[wldas_vars$wldas_name %in% available_vars]

    if (length(vars_to_use) == 0) {
      message(sprintf("WARNING: No matching variables in %s", basename(nc_file)))
      return(NULL)
    }

    r_subset <- terra::subset(r, vars_to_use)

    # Parse date from filename: WLDAS_NOAHMP001_DA1_YYYYMMDD.D10.nc
    date_str <- str_extract(basename(nc_file), "(?<=_)\\d{8}(?=\\.)")
    date_val <- as.Date(date_str, format = "%Y%m%d")

    list(raster = r_subset, date = date_val)

  }, error = function(e) {
    message(sprintf("ERROR loading %s: %s", basename(nc_file), e$message))
    return(NULL)
  })
}

# ============================================================================
# GET LIST OF WLDAS FILES
# ============================================================================
log_msg("\nFinding WLDAS files...")
nc_files <- list.files(wldas_dir, pattern = "\\.nc$", full.names = TRUE) %>% sort()

if (length(nc_files) == 0) stop(sprintf("No NetCDF files found in %s", wldas_dir))

log_msg(sprintf("Found %d WLDAS files", length(nc_files)))

get_date_from_filename <- function(filepath) {
  date_str <- str_extract(basename(filepath), "(?<=_)\\d{8}(?=\\.)")
  if (is.na(date_str)) return(NA)
  as.Date(date_str, format = "%Y%m%d")
}

file_dates       <- map_vec(nc_files, get_date_from_filename)
file_dates_valid <- file_dates[!is.na(file_dates)]

if (length(file_dates_valid) == 0) stop("Could not extract dates from any WLDAS files")

date_range <- range(file_dates_valid, na.rm = TRUE)
log_msg(sprintf("Date range: %s to %s (%d days)",
                date_range[1], date_range[2], length(file_dates_valid)))

# ============================================================================
# CREATE CHECKPOINT
# ============================================================================
checkpoint_file <- file.path(base_output_dir,
                             sprintf("extract_wldas.chunk%d_of_%d.checkpoint",
                                     chunk_index, num_chunks))

completed_files <- if (file.exists(checkpoint_file)) read_lines(checkpoint_file) else character(0)
if (length(completed_files) > 0) log_msg(sprintf("Found %d previously completed files", length(completed_files)))

nc_files_to_process <- nc_files[!(basename(nc_files) %in% completed_files)]
log_msg(sprintf("Files to process: %d (out of %d)", length(nc_files_to_process), length(nc_files)))

# ============================================================================
# PARALLEL SETUP
# ============================================================================
n_workers <- 16
plan(multicore, workers = n_workers)
log_msg(sprintf("Using %d parallel workers for NetCDF extraction\n", n_workers))

# ============================================================================
# FUNCTION: Extract spatial statistics for one variable, all watersheds
# ============================================================================
# Returns a named list of vectors (one element per stat), each of length
# nrow(ws_sf). The caller assembles these into the output tibble.
#
# Arguments:
#   var_raster   : single-layer SpatRaster in NATIVE (pre-conversion) units
#   ws_sf        : sf watersheds
#   output_name  : character, used to name the output columns
#   conversion   : numeric or NA; multiplied after extraction
#   compute_wet  : logical; whether to compute wetted fraction for this var
#   wet_thresh   : threshold in native units for wetted fraction

extract_var_stats <- function(var_raster, ws_sf, output_name, conversion,
                              compute_wet, wet_thresh, cv_min) {

  has_conversion <- !is.na(conversion)

  # ---- 1. Native-unit statistics needed before conversion ------------------
  # Wetted fraction requires comparison in native units, so compute it now.
  # We also extract mean and stdev in native units, then convert.

  # exactextractr accepts a character vector of built-in summary functions;
  # this reads each cell's values only ONCE and computes all stats together,
  # which is meaningfully faster than multiple exact_extract() calls.
  native_stats <- exactextractr::exact_extract(
    var_raster, ws_sf,
    fun = c("mean", "max", "stdev"),
    progress = FALSE
  )
  # exact_extract returns a data.frame when fun has length > 1;
  # columns are named "mean", "max", "stdev"

  mean_native <- native_stats[["mean"]]
  max_native  <- native_stats[["max"]]
  sd_native   <- native_stats[["stdev"]]

  # ---- 2. Apply unit conversion --------------------------------------------
  if (has_conversion) {
    mean_conv <- mean_native * conversion
    max_conv  <- max_native  * conversion
    sd_conv   <- sd_native   * conversion   # SD scales linearly with conversion
  } else {
    mean_conv <- mean_native
    max_conv  <- max_native
    sd_conv   <- sd_native
  }

  # ---- 3. Coefficient of variation -----------------------------------------
  # Suppress CV when the converted mean is below cv_min to avoid instability.
  cv_conv <- if_else(
    !is.na(mean_conv) & mean_conv >= cv_min,
    sd_conv / mean_conv,
    NA_real_
  )

  # ---- 4. Wetted fraction (flux variables only) ----------------------------
  # Uses a custom summary function operating on raw cell values (pre-conversion).
  # exact_extract passes a data.frame with columns 'value' and 'coverage_fraction'
  # to custom functions; we ignore coverage_fraction here (treat all cells equally)
  # since we want a simple cell count ratio. If area-weighted wetted fraction is
  # preferred in future, weight by coverage_fraction instead.
  if (compute_wet) {
    wet_frac <- exactextractr::exact_extract(
      var_raster, ws_sf,
      fun = function(values, coverage_fractions) {
        valid <- values[!is.na(values)]
        if (length(valid) == 0L) return(NA_real_)
        sum(valid > wet_thresh) / length(valid)
      },
      progress = FALSE
    )
  } else {
    wet_frac <- rep(NA_real_, nrow(ws_sf))
  }

  # ---- 5. Return named list of vectors -------------------------------------
  out <- list(
    mean     = as.numeric(mean_conv),
    max      = as.numeric(max_conv),
    sd       = as.numeric(sd_conv),
    cv       = as.numeric(cv_conv)
  )
  if (compute_wet) out[["wet_frac"]] <- as.numeric(wet_frac)

  out
}

# ============================================================================
# FUNCTION: Process and save results for a single NetCDF file
# ============================================================================
process_and_save_file <- function(nc_file, ws_sf, wldas_extent, wldas_vars,
                                  parquet_output_dir, checkpoint_file,
                                  wetted_thresh, wetted_frac_vars,
                                  spatial_stats_vars, cv_min) {

  tryCatch({
    data <- load_wldas_nc(nc_file, wldas_extent)
    if (is.null(data)) return(NULL)

    date_val <- data$date
    if (is.na(date_val)) return(NULL)

    r <- data$raster

    # Start output with identifiers
    extracted_vars <- tibble(
      date    = date_val,
      GAGE_ID = ws_sf$GAGE_ID
    )

    for (var_idx in seq_len(nrow(wldas_vars))) {
      wldas_name          <- wldas_vars$wldas_name[var_idx]
      output_name         <- wldas_vars$output_name[var_idx]
      conversion          <- wldas_vars$conversion_factor[var_idx]
      compute_wet         <- output_name %in% wetted_frac_vars
      compute_spatial     <- output_name %in% spatial_stats_vars

      if (!wldas_name %in% names(r)) {
        extracted_vars[[output_name]] <- NA_real_
        if (compute_spatial) {
          extracted_vars[[paste0(output_name, "_max")]] <- NA_real_
          extracted_vars[[paste0(output_name, "_sd")]]  <- NA_real_
          extracted_vars[[paste0(output_name, "_cv")]]  <- NA_real_
        }
        if (compute_wet)
          extracted_vars[[paste0(output_name, "_wet_frac")]] <- NA_real_
        next
      }

      var_raster <- terra::subset(r, wldas_name)

      if (compute_spatial) {
        # Full spatial statistics for rainfall/snowfall/snowmelt
        stats <- extract_var_stats(
          var_raster  = var_raster,
          ws_sf       = ws_sf,
          output_name = output_name,
          conversion  = conversion,
          compute_wet = compute_wet,
          wet_thresh  = wetted_thresh,
          cv_min      = cv_min
        )
        extracted_vars[[output_name]]                  <- stats$mean
        extracted_vars[[paste0(output_name, "_max")]]  <- stats$max
        extracted_vars[[paste0(output_name, "_sd")]]   <- stats$sd
        extracted_vars[[paste0(output_name, "_cv")]]   <- stats$cv
        if (compute_wet)
          extracted_vars[[paste0(output_name, "_wet_frac")]] <- stats$wet_frac
      } else {
        # Basin mean only for all other variables
        vals <- exactextractr::exact_extract(
          var_raster, ws_sf,
          fun      = "mean",
          progress = FALSE
        )
        if (!is.na(conversion)) vals <- vals * conversion
        extracted_vars[[output_name]] <- as.numeric(vals)
      }
    }

    list(
      nc_file   = basename(nc_file),
      file_date = date_val,
      data      = extracted_vars,
      n_rows    = nrow(extracted_vars)
    )

  }, error = function(e) {
    return(NULL)
  })
}

# ============================================================================
# MAIN EXTRACTION LOOP: Parallel across NetCDF Files
# ============================================================================
# Each worker processes one NetCDF file, writes its own Parquet file, and
# discards its result. Since each input file maps to a unique output date,
# there are no write conflicts between workers. The checkpoint is rebuilt
# from the output directory after the run rather than being written
# incrementally, which avoids any race condition on that file.
# ============================================================================
nc_files_to_process <- nc_files_to_process[!str_detect(nc_files_to_process, "test")]
#nc_files_to_process <- nc_files_to_process[1:100]
n_files    <- length(nc_files_to_process)
start_time <- now()

log_msg(sprintf("\nStarting parallel extraction of %d files...\n", n_files))

future_walk(
  nc_files_to_process,
  function(nc_file) {
    tryCatch({
      result <- process_and_save_file(
        nc_file, watersheds, wldas_extent, wldas_vars,
        parquet_output_dir, checkpoint_file,
        wetted_thresh      = WETTED_THRESH_RAW,
        wetted_frac_vars   = WETTED_FRAC_VARS,
        spatial_stats_vars = SPATIAL_STATS_VARS,
        cv_min             = CV_MEAN_MIN_MM
      )
      if (!is.null(result)) {
        parquet_file <- file.path(parquet_output_dir,
                                  sprintf("%s_wldas.parquet",
                                          format(result$file_date, "%Y%m%d")))
        write_parquet(result$data, parquet_file)
      }
    }, error = function(e) NULL)
  },
  .progress = TRUE,
  .options  = furrr_options(seed = TRUE)
)

# ============================================================================
# CLEANUP AND SUMMARY
# ============================================================================
plan(sequential)

elapsed_time <- now() - start_time
log_msg(sprintf("\n\n=== EXTRACTION COMPLETE ==="))
log_msg(sprintf("Elapsed time: %s", format(elapsed_time)))
log_msg(sprintf("Output directory: %s", parquet_output_dir))

# Rebuild checkpoint from completed Parquet files so resumption works
# correctly on re-runs. Derives original NC filename from parquet filename.
parquet_files <- list.files(parquet_output_dir, pattern = "\\.parquet$")
completed_nc  <- parquet_files %>%
  str_remove("_wldas\\.parquet$") %>%
  str_c("WLDAS_NOAHMP001_DA1_", ., ".D10.nc")
write_lines(completed_nc, checkpoint_file)
log_msg(sprintf("Checkpoint rebuilt: %d files recorded", length(completed_nc)))
log_msg(sprintf("Checkpoint file: %s", checkpoint_file))

log_msg(sprintf("Created/updated %d Parquet files", length(parquet_files)))

# Report the new column structure for one file as a sanity check
if (length(parquet_files) > 0) {
  sample_file <- file.path(parquet_output_dir, parquet_files[[1]])
  sample_cols <- names(read_parquet(sample_file, as_data_frame = FALSE)$schema)
  log_msg(sprintf("Output columns (%d): %s", length(sample_cols),
                  paste(sample_cols, collapse = ", ")))
}