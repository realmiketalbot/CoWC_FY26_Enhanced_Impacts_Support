#!/usr/bin/env Rscript
# ============================================================================
# Script 1b: Basin boundary polygons for the selected stations (App Data tab)
# ============================================================================
# Extracts the GAGES-II basin boundary for every station in
# _data/station_spatial.rds (Script 1) and writes them, unsimplified, in
# EPSG:4326 to _data/boundaries_spatial.rds (columns: site_no, geometry).
# The App's Data tab overlays the selected site's boundary.
#
# This file was previously created ad hoc for the 206-site release; this
# script makes it reproducible for the current station list.
#
# Usage: Rscript _scripts/1b_build_boundaries_spatial.R [gagesii_boundaries_dir]
#   default dir: /davenport-scratch/DATA/GAGES-II/boundaries_shapefiles_by_aggeco/
#                boundaries-shapefiles-by-aggeco
# ============================================================================

suppressPackageStartupMessages(library(sf))

args <- commandArgs(trailingOnly = TRUE)
bnd_dir <- if (length(args)) args[1] else
  "/davenport-scratch/DATA/GAGES-II/boundaries_shapefiles_by_aggeco/boundaries-shapefiles-by-aggeco"
output_file <- "_data/boundaries_spatial.rds"

stations <- readRDS("_data/station_spatial.rds")
shps <- list.files(path.expand(bnd_dir), pattern = "[.]shp$", full.names = TRUE)
if (!length(shps)) stop("No shapefiles found in ", bnd_dir)

b <- do.call(rbind, lapply(shps, function(f) {
  x <- st_read(f, quiet = TRUE)
  x <- x[x$GAGE_ID %in% stations$site_no, "GAGE_ID"]
  st_transform(x, 4326)
}))
names(b)[names(b) == "GAGE_ID"] <- "site_no"
b <- b[!duplicated(b$site_no), ]
b <- b[order(b$site_no), ]

missing <- setdiff(stations$site_no, b$site_no)
if (length(missing)) warning("No boundary for: ", paste(missing, collapse = ", "))
message(sprintf("%d of %d stations have boundaries", nrow(b), nrow(stations)))

saveRDS(b, output_file)
message("Saved ", output_file)
