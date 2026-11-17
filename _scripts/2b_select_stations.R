#!/usr/bin/env Rscript

# ============================================================================
# Script 2b: Station selection from flood-season record coverage
# ============================================================================
#
# Applies the station screen to the coverage table from Script 2 and defines
# which water years count as observed at each site. Runs locally (no access
# to raw unit data needed).
#
# RULE (decided 2026-11-17; see _planning/agu2026_interim_plan.md, D-E)
#   A water year QUALIFIES at a site when at least min_season_frac of the
#   days in its core flood season (flood_season_months, May-Sep) have an
#   observed daily value, AND at least min_approved_frac of its flood-season
#   unit values are USGS-approved ("A"). A site is SELECTED when it has at
#   least min_qual_years qualifying water years within [wy_min, wy_max].
#   Qualifying years need not be consecutive. Peaks outside the core season
#   (e.g., April or October) still count whenever they were observed; the
#   core season only decides whether a year counts toward record length.
#
# WHY A FLOOD SEASON
#   Unit records are routinely missing in winter (ice-affected or seasonal
#   operation): network-mean coverage is ~90% May-Sep but ~50% Dec-Jan. At
#   site-years where winter was observed, 97% of annual maxima fell in
#   Apr-Oct, and Nov-Mar annual maxima were typically half the site's median
#   annual maximum. Winter gaps therefore hide almost no relevant floods.
#
# WHY MAY-SEP RATHER THAN APR-OCT
#   21 gages operate mostly May-Sep (Apr ~57%, Oct ~59% observed). They are
#   high-elevation (median 3,210 m vs 2,750 m), snowmelt-dominated sites. In
#   years when their April or October WAS observed, only 2.4% (Apr) and 0.4%
#   (Oct) of their large peaks (top ~1 per year) fell in those months, versus
#   4.6% and 1.2% at other sites. Requiring Apr-Oct would drop these good
#   flood records; requiring May-Sep keeps 16 of them. The cost is a slight
#   undercount of lambda where an April/October flood went unrecorded in a
#   qualifying year, estimated at under ~1% for these gages and less elsewhere.
#   Site counts (>= 90%, >= 20 years, approval condition on): May-Sep 215,
#   Apr-Oct 194.
#
# WHY THE APPROVAL CONDITION
#   It removes a handful of mostly provisional WY2022-2023 site-years and no
#   sites, so qualifying years rest on reviewed data only.
#
# DOWNSTREAM USE OF QUALIFYING YEARS
#   - Poisson exceedance rate: lambda = exceedances in qualifying years /
#     number of qualifying years. Peaks in non-qualifying years are excluded
#     from fitting (Script 7), so numerator and denominator match.
#   - Dispersion index: M = number of qualifying years.
#
# WY_MAX = 2023 because WLDAS forcing ends 2023-12-31, so WY2024 peaks
# cannot be classified.
#
# INPUTS:
#   - _data/usgs_unit_coverage_co.rds (Script 2)
#
# OUTPUTS:
#   - _data/selected_unit_stations.csv : site_no of selected stations
#   - _data/site_qualifying_years.rds  : one row per candidate site x water
#       year in [wy_min, wy_max] with season_frac, approved_frac and
#       qualifies (TRUE/FALSE),
#       plus n_qual_years and selected per site. Attribute "selection_config".
# ============================================================================

library(tidyverse)

# ============================================================================
# CONFIGURATION
# ============================================================================
coverage_file   <- "_data/usgs_unit_coverage_co.rds"
stations_file   <- "_data/selected_unit_stations.csv"
qual_years_file <- "_data/site_qualifying_years.rds"

SELECTION_CONFIG <- list(
  flood_season_months = 5:9,    # May-September (core flood season)
  min_season_frac     = 0.90,
  min_approved_frac   = 0.95,  # share of flood-season unit values flagged "A"
  min_qual_years      = 20L,
  wy_min              = 1986L,  # start of unit records in this network
  wy_max              = 2023L   # WLDAS forcing ends 2023-12-31
)

cfg <- SELECTION_CONFIG

# ============================================================================
# QUALIFYING YEARS
# ============================================================================
coverage <- readRDS(coverage_file)

qual <- coverage %>%
  filter(water_year >= cfg$wy_min, water_year <= cfg$wy_max,
         month %in% cfg$flood_season_months) %>%
  group_by(site_no, water_year) %>%
  summarize(
    season_days_obs = sum(n_days_obs),
    season_days     = sum(days_in_month),
    season_frac     = season_days_obs / season_days,
    approved_frac   = if (sum(n_obs) > 0) sum(n_obs_approved) / sum(n_obs) else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(qualifies = season_frac >= cfg$min_season_frac &
                     coalesce(approved_frac >= cfg$min_approved_frac, FALSE)) %>%
  group_by(site_no) %>%
  mutate(
    n_qual_years = sum(qualifies),
    selected     = n_qual_years >= cfg$min_qual_years
  ) %>%
  ungroup() %>%
  arrange(site_no, water_year)

attr(qual, "selection_config") <- SELECTION_CONFIG

selected <- qual %>% filter(selected) %>% distinct(site_no) %>% pull(site_no) %>% sort()

# ============================================================================
# SUMMARY
# ============================================================================
site_tbl <- qual %>% distinct(site_no, n_qual_years, selected)

message(sprintf("Candidate sites (Colorado GAGES-II with unit data): %d",
                n_distinct(coverage$site_no)))
message(sprintf("Rule: %s-%s >= %.0f%% observed and >= %.0f%% approved; >= %d qualifying water years in WY%d-%d",
                month.abb[min(cfg$flood_season_months)], month.abb[max(cfg$flood_season_months)],
                100 * cfg$min_season_frac, 100 * cfg$min_approved_frac,
                cfg$min_qual_years, cfg$wy_min, cfg$wy_max))
message(sprintf("Selected: %d sites | qualifying years per selected site: median %d (range %d-%d)",
                length(selected),
                as.integer(median(site_tbl$n_qual_years[site_tbl$selected])),
                min(site_tbl$n_qual_years[site_tbl$selected]),
                max(site_tbl$n_qual_years[site_tbl$selected])))

# Robustness: site counts under neighbouring minimum-year settings
for (n in c(20L, 25L, 30L)) {
  message(sprintf("  sites with >= %d qualifying years: %d", n, sum(site_tbl$n_qual_years >= n)))
}

# Change relative to the previous station list, if present
if (file.exists(stations_file)) {
  prev <- read_csv(stations_file, show_col_types = FALSE, col_types = cols(site_no = col_character()))$site_no
  message(sprintf("Versus previous list (%d sites): kept %d, dropped %d, added %d",
                  length(prev), length(intersect(prev, selected)),
                  length(setdiff(prev, selected)), length(setdiff(selected, prev))))
}

# ============================================================================
# SAVE
# ============================================================================
write_csv(tibble(site_no = selected), stations_file)
saveRDS(qual, qual_years_file)
message(sprintf("Wrote %s and %s", stations_file, qual_years_file))
