#!/usr/bin/env Rscript

# ============================================================================
# Build station spatial dataset for Shiny map tab
# ============================================================================
# Reads selected stations, joins GAGES-II metadata (location, basin
# characteristics, topography), and saves an sf object RDS for use
# in the dashboard map tab.
#
# Output: _data/station_spatial.rds
#   An sf object (points, WGS84) with one row per station and
#   columns for all dashboard-relevant metadata.
# ============================================================================

library(tidyverse)
library(sf)

# ============================================================================
# FILE PATHS
# ============================================================================
gagesii_dir       <- "/davenport-scratch/DATA/GAGES-II/basinchar_and_report_sept_2011/spreadsheets-in-csv-format"
stations_file     <- "_data/selected_unit_stations.csv"
output_file       <- "_data/station_spatial.rds"

# GAGES-II metadata files
basinid_file      <- file.path(gagesii_dir, "conterm_basinid.txt")
topo_file         <- file.path(gagesii_dir, "conterm_topo.txt")
classif_file      <- file.path(gagesii_dir, "conterm_bas_classif.txt")

# ============================================================================
# LOAD SELECTED STATIONS
# ============================================================================
message("Loading selected stations...")

stations <- read_csv(stations_file, show_col_types = FALSE) %>%
  mutate(site_no = str_pad(as.character(site_no), 8, "left", "0"))

message(sprintf("  %d stations loaded", nrow(stations)))

# ============================================================================
# LOAD GAGES-II METADATA
# ============================================================================
message("Loading GAGES-II metadata...")

# conterm_basinid: location, drainage area, HUC, state
basinid <- read_csv(basinid_file, show_col_types = FALSE) %>%
  transmute(
    site_no      = str_pad(as.character(STAID), 8, "left", "0"),
    station_name = STANAME,
    lon          = LNG_GAGE,
    lat          = LAT_GAGE,
    drain_sqkm   = DRAIN_SQKM,
    state        = STATE,
  )

# conterm_topo: elevation, relief, slope
topo <- read_csv(topo_file, show_col_types = FALSE) %>%
  transmute(
    site_no      = str_pad(as.character(STAID), 8, "left", "0"),
    elev_mean_m  = ELEV_MEAN_M_BASIN,
    elev_max_m   = ELEV_MAX_M_BASIN,
    slope_pct    = SLOPE_PCT
  )

# conterm_bas_classif: Ref vs Non-Ref classification
classif <- read_csv(classif_file, show_col_types = FALSE) %>%
  transmute(
    site_no            = str_pad(as.character(STAID), 8, "left", "0"),
    class              = CLASS,       # "Ref" or "Non-Ref"
    hydro_disturb_indx = HYDRO_DISTURB_INDX
  )

# ============================================================================
# JOIN
# ============================================================================
message("Joining metadata...")

station_meta <- stations %>%
  left_join(basinid,  by = "site_no") %>%
  left_join(topo,     by = "site_no") %>%
  left_join(classif,  by = "site_no")

# Report any stations missing coordinates (shouldn't happen for GAGES-II sites)
n_missing_coords <- sum(is.na(station_meta$lat) | is.na(station_meta$lon))
if (n_missing_coords > 0) {
  warning(sprintf(
    "%d stations have missing coordinates and will be dropped from spatial dataset",
    n_missing_coords
  ))
  missing <- station_meta %>% filter(is.na(lat) | is.na(lon)) %>% pull(site_no)
  message("  Missing: ", paste(missing, collapse = ", "))
  station_meta <- station_meta %>% filter(!is.na(lat) & !is.na(lon))
}

message(sprintf("  %d stations with complete metadata", nrow(station_meta)))

# ============================================================================
# BUILD SF OBJECT
# ============================================================================
message("Building sf object...")

sv <- st_as_sf(
  station_meta,
  coords = c("lon", "lat"),
  crs    = 4326,
  remove = FALSE   # retain lon/lat columns for convenience
)

message(sprintf("  sf object: %d points, %d attributes, CRS = %s",
                nrow(sv), ncol(sv) - 1, st_crs(sv)$epsg))

# ============================================================================
# SAVE
# ============================================================================
saveRDS(sv, output_file)
message(sprintf("Saved to %s", output_file))

# Quick summary for sanity check
message("\n--- Attribute summary ---")
df <- st_drop_geometry(sv)
message(sprintf("  Columns: %s", paste(names(df), collapse = ", ")))
message(sprintf("  States:  %s", paste(sort(unique(df$state)), collapse = ", ")))
message(sprintf("  Class:   %d Ref / %d Non-Ref",
                sum(df$class == "Ref",     na.rm = TRUE),
                sum(df$class == "Non-Ref", na.rm = TRUE)))
message(sprintf("  Elev range (m): %.0f – %.0f",
                min(df$elev_mean_m, na.rm = TRUE),
                max(df$elev_mean_m, na.rm = TRUE)))