#!/usr/bin/env Rscript

# ============================================================================
# Identify sites viable for full three-component mixture model
# (Yu et al. 2022 framework: rainfall + snowmelt + rain-on-snow)
# Minimum threshold: >= 20 peaks per flood type (Yu et al. criterion)
# ============================================================================

library(tidyverse)

# ============================================================================
# CONFIGURATION
# ============================================================================
classified_file <- "_data/pot_peaks_classified.rds"
watershed_meta  <- "/davenport-scratch/DATA/GAGES-II/basinchar_and_report_sept_2011/spreadsheets-in-csv-format/conterm_basinid.txt"
topo_meta_file  <- "/davenport-scratch/DATA/GAGES-II/basinchar_and_report_sept_2011/spreadsheets-in-csv-format/conterm_topo.txt"
output_file     <- "_data/mixture_viable_sites.csv"

min_peaks_per_type <- 20

# ============================================================================
# LOAD DATA
# ============================================================================
peaks <- readRDS(classified_file) %>%
  filter(flood_type != "unclassified")

basin_meta <- read_csv(watershed_meta, show_col_types = FALSE) %>%
  filter(STATE == "CO") %>%
  transmute(site_no = as.character(STAID), drain_sqkm = DRAIN_SQKM)

topo_meta <- read_csv(topo_meta_file, show_col_types = FALSE) %>%
  transmute(site_no = as.character(STAID), elev_m = ELEV_MEAN_M_BASIN)

site_meta <- left_join(basin_meta, topo_meta, by = "site_no")

# ============================================================================
# COUNT PEAKS PER SITE PER TYPE
# ============================================================================
type_counts <- peaks %>%
  count(site_no, flood_type) %>%
  pivot_wider(names_from = flood_type, values_from = n, values_fill = 0) %>%
  # Ensure all three columns exist even if no sites have them
  { if (!"rainfall"    %in% names(.)) mutate(., rainfall    = 0L) else . } %>%
  { if (!"snowmelt"    %in% names(.)) mutate(., snowmelt    = 0L) else . } %>%
  { if (!"rain_on_snow" %in% names(.)) mutate(., rain_on_snow = 0L) else . } %>%
  mutate(
    total_classified   = rainfall + snowmelt + rain_on_snow,
    viable_rainfall    = rainfall    >= min_peaks_per_type,
    viable_snowmelt    = snowmelt    >= min_peaks_per_type,
    viable_ros         = rain_on_snow >= min_peaks_per_type,
    viable_full_mixture = viable_rainfall & viable_snowmelt & viable_ros,
    # Partial mixture viability
    viable_rain_snow   = viable_rainfall & viable_snowmelt,
    viable_rain_ros    = viable_rainfall & viable_ros,
    viable_snow_ros    = viable_snowmelt & viable_ros,
    n_viable_types     = viable_rainfall + viable_snowmelt + viable_ros
  ) %>%
  left_join(site_meta, by = "site_no") %>%
  arrange(desc(viable_full_mixture), desc(n_viable_types), desc(total_classified))

# ============================================================================
# SUMMARY
# ============================================================================
message("=== Mixture model viability summary ===\n")

message(sprintf("Total sites with any classified peaks: %d", nrow(type_counts)))
message(sprintf("Sites viable for full mixture (all 3 types >= %d): %d",
                min_peaks_per_type,
                sum(type_counts$viable_full_mixture)))
message(sprintf("Sites viable for rainfall + snowmelt only:         %d",
                sum(type_counts$viable_rain_snow & !type_counts$viable_ros)))
message(sprintf("Sites viable for rainfall + ROS only:              %d",
                sum(type_counts$viable_rain_ros & !type_counts$viable_snowmelt)))
message(sprintf("Sites viable for rainfall only:                    %d",
                sum(type_counts$viable_rainfall & type_counts$n_viable_types == 1)))
message(sprintf("Sites with only 1 viable type total:               %d",
                sum(type_counts$n_viable_types == 1)))
message(sprintf("Sites with no viable types (all < %d peaks):       %d",
                min_peaks_per_type,
                sum(type_counts$n_viable_types == 0)))

message("\n--- Full mixture sites ---")
full_mixture_sites <- type_counts %>%
  filter(viable_full_mixture) %>%
  select(site_no, rainfall, snowmelt, rain_on_snow, total_classified,
         elev_m, drain_sqkm)

print(full_mixture_sites, n = Inf)

message(sprintf("\nElevation of full mixture sites (m):"))
print(summary(full_mixture_sites$elev_m))

message("\n--- Distribution of viable type counts across all sites ---")
print(count(type_counts, n_viable_types))

# ============================================================================
# SAVE
# ============================================================================
write_csv(type_counts, output_file)
message(sprintf("\nFull site viability table saved to %s", output_file))