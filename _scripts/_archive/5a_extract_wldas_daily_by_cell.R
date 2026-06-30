#!/usr/bin/env Rscript

# ============================================================================
# Script 6a: WLDAS Per-Cell Extraction to GAGES-II Watersheds
# ============================================================================
# Extracts per-cell WLDAS values for all Colorado GAGES-II watersheds,
# retaining cell identity (lon/lat centroid) and exactextractr coverage
# fraction for area-weighted aggregation in script 6b.
#
# Only rainf_flux and qsm_flux are extracted — the minimal set needed for
# Davenport et al. cell-level flood classification in script 6b.
#
# Parallelisation strategy mirrors script 4: future_walk over NetCDF files,
# each worker writes its own Parquet and discards results immediately.
# Checkpoint is rebuilt from the output directory after the run.
#
# Output schema (one Parquet per date):
#   date          : Date
#   GAGE_ID       : character  watershed identifier
#   cell_lon      : numeric    cell centroid longitude (degrees, EPSG:4326)
#   cell_lat      : numeric    cell centroid latitude  (degrees, EPSG:4326)
#   coverage_frac : numeric    fraction of cell area within watershed [0, 1]
#   rainf_flux    : numeric    rainfall flux (mm/day, converted from kg/m2/s)
#   qsm_flux      : numeric    snowmelt flux (mm/day, converted from kg/m2/s)
#
# INPUTS:
#   - WLDAS NetCDF files (raw, same directory as script 4)
#   - GAGES-II watershed boundary shapefiles
#   - GAGES-II basin metadata (for Colorado filtering)
#
# OUTPUTS:
#   - _data/WLDAS/processed/cell_daily_values/
#       YYYYMMDD_wldas_cells.parquet  (one per date)
#   - _data/WLDAS/processed/
#       extract_wldas_cells.chunk{i}_of_{n}.checkpoint
#       extraction_cells_LOG.txt
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
watershed_file     <- "/davenport-scratch/DATA/GAGES-II/boundaries_shapefiles_by_aggeco/boundaries-shapefiles-by-aggeco/bas_ref_all.shp"
watershed_meta     <- "/davenport-scratch/DATA/GAGES-II/basinchar_and_report_sept_2011/spreadsheets-in-csv-format/conterm_basinid.txt"
wldas_dir          <- "_data/WLDAS/raw"
base_output_dir    <- "_data/WLDAS/processed"
parquet_output_dir <- file.path(base_output_dir, "cell_daily_values")

# Conversion factors from native WLDAS flux units (kg/m2/s) to mm/day
RAINF_CONVERSION <- 86400
QSM_CONVERSION   <- 86400

# Command line arguments for chunking (optional)
# Usage: Rscript 6a_extract_wldas_cells_to_gages-ii.R [chunk_index] [num_chunks]
args        <- commandArgs(trailingOnly = TRUE)
chunk_index <- ifelse(length(args) >= 1, as.integer(args[1]), 1)
num_chunks  <- ifelse(length(args) >= 2, as.integer(args[2]), 1)

# ============================================================================
# SETUP LOGGING
# ============================================================================
log_file <- file.path(base_output_dir, "extraction_cells_LOG.txt")
dir.create(base_output_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(parquet_output_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  msg <- paste(..., sep = " ")
  cat(sprintf("[%s] %s\n", Sys.time(), msg))
  cat(sprintf("[%s] %s\n", Sys.time(), msg), file = log_file, append = TRUE)
}

log_msg("=== WLDAS Cell Extraction Started ===")
log_msg(sprintf("Chunk: %d of %d", chunk_index, num_chunks))
log_msg(sprintf("Output directory: %s", parquet_output_dir))
log_msg(sprintf("Variables: rainf_flux (x%d), qsm_flux (x%d)",
                RAINF_CONVERSION, QSM_CONVERSION))

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

watershed_dir  <- dirname(watershed_file)
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
# WLDAS EXTENT (mirrors script 4)
# ============================================================================
get_wldas_extent <- function() {
  resolution <- 0.01
  xmin <- -125.0
  ymin <-   25.0
  xmax <- xmin + (3591 * resolution)
  ymax <- ymin + (2787 * resolution)
  ext(xmin, xmax, ymin, ymax)
}

wldas_extent <- get_wldas_extent()
log_msg(sprintf("WLDAS extent set: lon=[%.2f, %.2f], lat=[%.2f, %.2f]",
                xmin(wldas_extent), xmax(wldas_extent),
                ymin(wldas_extent), ymax(wldas_extent)))

# ============================================================================
# FUNCTION: Load single WLDAS NetCDF file (mirrors script 4)
# ============================================================================
load_wldas_nc <- function(nc_file, extent) {
  tryCatch({
    r <- terra::rast(nc_file)
    terra::ext(r) <- extent
    terra::crs(r) <- "EPSG:4326"

    date_str <- str_extract(basename(nc_file), "(?<=_)\\d{8}(?=\\.)")
    date_val <- as.Date(date_str, format = "%Y%m%d")

    list(raster = r, date = date_val)
  }, error = function(e) {
    message(sprintf("ERROR loading %s: %s", basename(nc_file), e$message))
    NULL
  })
}

# ============================================================================
# FUNCTION: Extract per-cell values for one variable, all watersheds
# ============================================================================
# exact_extract with fun = NULL returns a list of data.frames (one per
# watershed feature), each with columns: value, coverage_fraction, x, y.
# We bind these into a single tibble, attaching GAGE_ID from ws_sf.
# ============================================================================
extract_cells_one_var <- function(var_raster, ws_sf, output_name, conversion) {

  cell_list <- exactextractr::exact_extract(
    var_raster, ws_sf,
    fun        = NULL,   # return raw per-cell data.frames
    include_xy = TRUE,   # append cell centroid lon (x) and lat (y)
    progress   = FALSE
  )

  map2_dfr(cell_list, ws_sf$GAGE_ID, function(df, gid) {
    tibble(
      GAGE_ID       = gid,
      cell_lon      = df$x,
      cell_lat      = df$y,
      coverage_frac = df$coverage_fraction,
      val           = if (!is.na(conversion)) df$value * conversion else df$value
    ) %>%
      rename(!!output_name := val)
  })
}

# ============================================================================
# FUNCTION: Process and save one NetCDF file
# ============================================================================
process_and_save_file <- function(nc_file, ws_sf, wldas_extent,
                                  parquet_output_dir) {
  tryCatch({
    data <- load_wldas_nc(nc_file, wldas_extent)
    if (is.null(data))    return(invisible(NULL))
    if (is.na(data$date)) return(invisible(NULL))

    date_val <- data$date
    r        <- data$raster

    # ---- Extract rainf_flux -------------------------------------------------
    if (!"Rainf_tavg" %in% names(r)) {
      message(sprintf("WARNING: Rainf_tavg missing in %s", basename(nc_file)))
      return(invisible(NULL))
    }

    rain_cells <- extract_cells_one_var(
      var_raster  = terra::subset(r, "Rainf_tavg"),
      ws_sf       = ws_sf,
      output_name = "rainf_flux",
      conversion  = RAINF_CONVERSION
    )

    # ---- Extract qsm_flux ---------------------------------------------------
    if (!"Qsm_tavg" %in% names(r)) {
      message(sprintf("WARNING: Qsm_tavg missing in %s", basename(nc_file)))
      return(invisible(NULL))
    }

    qsm_cells <- extract_cells_one_var(
      var_raster  = terra::subset(r, "Qsm_tavg"),
      ws_sf       = ws_sf,
      output_name = "qsm_flux",
      conversion  = QSM_CONVERSION
    )

    # ---- Join on GAGE_ID + cell centroid ------------------------------------
    # Both variables are extracted from the same raster grid, so cell_lon and
    # cell_lat are identical between the two extractions. Join on all three
    # identity columns. coverage_frac is taken from rain_cells (same for both).
    result <- rain_cells %>%
      inner_join(
        qsm_cells %>% select(GAGE_ID, cell_lon, cell_lat, qsm_flux),
        by = c("GAGE_ID", "cell_lon", "cell_lat")
      ) %>%
      mutate(date = date_val) %>%
      select(date, GAGE_ID, cell_lon, cell_lat, coverage_frac,
             rainf_flux, qsm_flux)

    # ---- Write Parquet -------------------------------------------------------
    parquet_file <- file.path(
      parquet_output_dir,
      sprintf("%s_wldas_cells.parquet", format(date_val, "%Y%m%d"))
    )
    write_parquet(result, parquet_file)

    invisible(NULL)

  }, error = function(e) {
    message(sprintf("ERROR processing %s: %s", basename(nc_file), e$message))
    invisible(NULL)
  })
}

# ============================================================================
# GET LIST OF WLDAS FILES AND APPLY CHECKPOINT
# ============================================================================
log_msg("\nFinding WLDAS files...")
nc_files <- list.files(wldas_dir, pattern = "\\.nc$", full.names = TRUE) %>% sort()

if (length(nc_files) == 0) stop(sprintf("No NetCDF files found in %s", wldas_dir))
log_msg(sprintf("Found %d WLDAS files", length(nc_files)))

checkpoint_file <- file.path(
  base_output_dir,
  sprintf("extract_wldas_cells.chunk%d_of_%d.checkpoint", chunk_index, num_chunks)
)

completed_files <- if (file.exists(checkpoint_file)) read_lines(checkpoint_file) else character(0)
if (length(completed_files) > 0)
  log_msg(sprintf("Found %d previously completed files in checkpoint",
                  length(completed_files)))

nc_files_to_process <- nc_files[!(basename(nc_files) %in% completed_files)]
log_msg(sprintf("Files to process: %d (out of %d)",
                length(nc_files_to_process), length(nc_files)))

# ============================================================================
# PARALLEL SETUP
# ============================================================================
n_workers <- 16
plan(multicore, workers = n_workers)
log_msg(sprintf("Using %d parallel workers\n", n_workers))

# ============================================================================
# MAIN EXTRACTION LOOP
# ============================================================================
n_files    <- length(nc_files_to_process)
start_time <- now()

log_msg(sprintf("Starting parallel cell extraction of %d files...\n", n_files))

future_walk(
  nc_files_to_process,
  function(nc_file) {
    process_and_save_file(
      nc_file, watersheds, wldas_extent, parquet_output_dir
    )
  },
  .progress = TRUE,
  .options  = furrr_options(seed = TRUE)
)

# ============================================================================
# CLEANUP AND SUMMARY
# ============================================================================
plan(sequential)

elapsed_time <- now() - start_time
log_msg(sprintf("\n\n=== CELL EXTRACTION COMPLETE ==="))
log_msg(sprintf("Elapsed time: %s", format(elapsed_time)))
log_msg(sprintf("Output directory: %s", parquet_output_dir))

# Rebuild checkpoint from completed Parquet files
parquet_files <- list.files(parquet_output_dir, pattern = "\\.parquet$")
completed_nc  <- parquet_files %>%
  str_remove("_wldas_cells\\.parquet$") %>%
  str_c("WLDAS_NOAHMP001_DA1_", ., ".D10.nc")
write_lines(completed_nc, checkpoint_file)
log_msg(sprintf("Checkpoint rebuilt: %d files recorded", length(completed_nc)))
log_msg(sprintf("Created %d Parquet files", length(parquet_files)))

# Sanity check on one output file
if (length(parquet_files) > 0) {
  sample_file <- file.path(parquet_output_dir, parquet_files[[1]])
  sample      <- read_parquet(sample_file)
  log_msg(sprintf("Sample file: %s", parquet_files[[1]]))
  log_msg(sprintf("  Rows:      %d", nrow(sample)))
  log_msg(sprintf("  Watersheds: %d", n_distinct(sample$GAGE_ID)))
  log_msg(sprintf("  Columns:   %s", paste(names(sample), collapse = ", ")))
  log_msg(sprintf("  Cells per watershed (median): %.0f",
                  median(table(sample$GAGE_ID))))
}