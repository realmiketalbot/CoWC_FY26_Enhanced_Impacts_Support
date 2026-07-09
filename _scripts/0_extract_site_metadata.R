#!/usr/bin/env Rscript

# ============================================================================
# Build Site Metadata
# ============================================================================
# Assembles site_metadata.rds from GAGES-II attribute tables. This is a
# standalone script so that metadata can be updated independently of the
# classification pipeline, and so that all downstream scripts have a single
# authoritative source for site attributes.
#
# Regulation screening flag:
#   regulated = STOR_NOR_2009 >= 1.233 mcm (= 1,000 acre-feet)
#   This matches the upstream storage threshold used by Davenport et al. (2020)
#   to exclude regulated sites. Sites below this threshold are treated as
#   minimally regulated for the purposes of mass balance diagnostics and
#   end-member selection.
#
# INPUTS (GAGES-II spreadsheets):
#   - conterm_basinid.txt        : location, drainage area, HUC
#   - conterm_topo.txt           : elevation, slope
#   - conterm_hydromod_dams.txt  : dam count, storage, distance to nearest dam
#   - conterm_bas_classif.txt    : reference/non-reference, disturbance index
#
# OUTPUT:
#   - _data/site_metadata.rds    : one row per site
#
# COLUMNS IN OUTPUT:
#   site_no             : USGS gage ID (character)
#   station_name        : gage name
#   state               : state abbreviation
#   lat_deg, lon_deg    : gage coordinates
#   drain_sqkm          : drainage area (km²)
#   huc02               : 2-digit HUC
#   elev_mean_m         : mean basin elevation (m)
#   slope_pct           : mean basin slope (%)
#   ndams_2009          : number of upstream dams as of 2009
#   stor_norm_mcm       : total normal upstream storage (million m³)
#   stor_norm_acft      : total normal upstream storage (acre-feet)
#   dist_nearest_dam_km : distance to nearest dam (km); -999 = no dam recorded
#   class               : GAGES-II classification ("Ref" or "Non-ref")
#   hydro_disturb_idx   : hydrologic disturbance index
#   regulated           : logical; TRUE if stor_norm_mcm >= 1.233 mcm
#                         (Davenport et al. 2020 threshold = 1,000 acre-feet)
# ============================================================================

library(tidyverse)

# ============================================================================
# CONFIGURATION
# ============================================================================
gagesii_dir  <- "/davenport-scratch/DATA/GAGES-II/basinchar_and_report_sept_2011/spreadsheets-in-csv-format"
meta_basinid <- file.path(gagesii_dir, "conterm_basinid.txt")
meta_topo    <- file.path(gagesii_dir, "conterm_topo.txt")
meta_dams    <- file.path(gagesii_dir, "conterm_hydromod_dams.txt")
meta_classif <- file.path(gagesii_dir, "conterm_bas_classif.txt")

peaks_file   <- "_data/pot_peaks_all.rds"
output_file  <- "_data/site_metadata.rds"

# Davenport et al. (2020) regulation threshold: 1,000 acre-feet = 1.2335 mcm
REGULATION_THRESH_MCM  <- 1.2335   # million m³
MCM_PER_ACFT           <- 0.0012335

# ============================================================================
# LOAD SITE LIST FROM PEAKS
# ============================================================================
message("Loading site list from peaks file...")
peaks    <- readRDS(peaks_file)
site_ids <- unique(as.character(peaks$site_no))
message(sprintf("Found %d unique sites", length(site_ids)))

# ============================================================================
# HELPER: Read a GAGES-II table, standardise site_no, filter to our sites
# ============================================================================
read_gagesii <- function(path, site_ids) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(site_no = as.character(STAID)) %>%
    filter(site_no %in% site_ids)
}

# ============================================================================
# LOAD AND JOIN GAGES-II TABLES
# ============================================================================
message("Loading GAGES-II attribute tables...")

tbl_basinid <- read_gagesii(meta_basinid, site_ids) %>%
  transmute(
    site_no,
    station_name = STANAME,
    state        = STATE,
    lat_deg      = LAT_GAGE,
    lon_deg      = LNG_GAGE,
    drain_sqkm   = DRAIN_SQKM,
    huc02        = HUC02
  )

tbl_topo <- read_gagesii(meta_topo, site_ids) %>%
  transmute(
    site_no,
    elev_mean_m = ELEV_MEAN_M_BASIN,
    slope_pct   = SLOPE_PCT
  )

tbl_dams <- read_gagesii(meta_dams, site_ids) %>%
  transmute(
    site_no,
    ndams_2009          = NDAMS_2009,
    stor_norm_mcm       = STOR_NOR_2009,
    stor_norm_acft      = STOR_NOR_2009 / MCM_PER_ACFT,
    dist_nearest_dam_km = RAW_DIS_NEAREST_DAM
  )

tbl_classif <- read_gagesii(meta_classif, site_ids) %>%
  transmute(
    site_no,
    class             = CLASS,
    hydro_disturb_idx = HYDRO_DISTURB_INDX
  )

# ============================================================================
# JOIN AND DERIVE REGULATION FLAG
# ============================================================================
site_metadata <- tbl_basinid %>%
  left_join(tbl_topo,    by = "site_no") %>%
  left_join(tbl_dams,    by = "site_no") %>%
  left_join(tbl_classif, by = "site_no") %>%
  mutate(
    # Regulated if upstream storage exceeds Davenport et al. threshold.
    # NA storage (no dam recorded) treated as unregulated.
    regulated = !is.na(stor_norm_mcm) & stor_norm_mcm >= REGULATION_THRESH_MCM
  ) %>%
  arrange(site_no)

# ============================================================================
# SUMMARY
# ============================================================================
message(sprintf("\nSite metadata assembled for %d sites:", nrow(site_metadata)))
message(sprintf("  Reference (CLASS == 'Ref'):           %d",
                sum(site_metadata$class == "Ref",     na.rm = TRUE)))
message(sprintf("  Non-reference (CLASS != 'Ref'):       %d",
                sum(site_metadata$class != "Ref",     na.rm = TRUE)))
message(sprintf("  Sites with >= 1 dam upstream:         %d",
                sum(site_metadata$ndams_2009 > 0,     na.rm = TRUE)))
message(sprintf("  Regulated (storage >= %.3f mcm):      %d",
                REGULATION_THRESH_MCM,
                sum(site_metadata$regulated,          na.rm = TRUE)))
message(sprintf("  Unregulated (storage < %.3f mcm):     %d",
                REGULATION_THRESH_MCM,
                sum(!site_metadata$regulated,         na.rm = TRUE)))
message(sprintf("  Storage NA (no dam recorded):         %d",
                sum(is.na(site_metadata$stor_norm_mcm))))

# ============================================================================
# SAVE
# ============================================================================
saveRDS(site_metadata, output_file)
message(sprintf("\nSaved site metadata to %s", output_file))