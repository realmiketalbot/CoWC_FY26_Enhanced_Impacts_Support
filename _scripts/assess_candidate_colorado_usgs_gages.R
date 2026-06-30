#!/usr/bin/env Rscript

# ============================================================================
# Identify Colorado USGS gages outside GAGES-II with adequate record length
# ============================================================================
# Uses the new dataRetrieval Water Data API (read_waterdata_* functions),
# replacing the deprecated NWIS services (whatNWISsites, readNWISstat).
#
# Workflow:
#   1. Pull all CO streamgage time-series metadata for discharge (00060)
#      directly via read_waterdata_ts_meta() — gives begin/end per site
#   2. Filter to daily mean series with adequate record length and recency
#   3. Separately check for instantaneous (IV) series at the same sites
#   4. Cross-reference against GAGES-II Colorado sites
#   5. Report what's gained by expanding beyond GAGES-II
#
# Note: requires internet access. Runtime ~1-2 min.
# Rate limit: 1000 requests/hour per IP. If hitting limits, register a token
# at https://api.waterdata.usgs.gov/signup/ and set WATERDATA_TOKEN in .Renviron
# ============================================================================

library(tidyverse)
library(dataRetrieval)
library(sf)

# ============================================================================
# CONFIGURATION
# ============================================================================
min_years_short <- 20
min_years_long  <- 30
record_end_year <- 2022   # require data through at least this year
param_cd        <- "00060"  # discharge

gagesii_meta <- "/davenport-scratch/DATA/GAGES-II/basinchar_and_report_sept_2011/spreadsheets-in-csv-format/conterm_basinid.txt"
output_file  <- "_scripts/non_gagesii_colorado_gages.csv"

# ============================================================================
# STEP 1: PULL ALL COLORADO DISCHARGE TIME-SERIES METADATA
# ============================================================================
# read_waterdata_ts_meta() returns one row per time series (a site can have
# multiple: daily mean, instantaneous, max, min, etc.)
# We request only the fields we need to keep the response lean.

message("Querying time-series metadata for all Colorado discharge gages...")

ts_meta_raw <- read_waterdata_ts_meta(
  state_name     = "Colorado",
  parameter_code = param_cd,
  skipGeometry   = TRUE,
  properties     = c(
    "monitoring_location_id",
    "parameter_code",
    "computation_period_identifier",
    "computation_identifier",
    "statistic_id",
    "begin",
    "end"
  )
)

message(sprintf("Total discharge time-series rows returned: %d", nrow(ts_meta_raw)))

# ============================================================================
# STEP 2: FILTER TO DAILY MEAN SERIES
# ============================================================================
# computation_period_identifier == "Daily" and statistic_id == "00003" (mean)
# This excludes instantaneous, max, min, and other aggregate series.

daily_mean <- ts_meta_raw %>%
  filter(
    computation_period_identifier == "Daily",
    statistic_id == "00003"
  ) %>%
  mutate(
    site_no    = str_remove(monitoring_location_id, "^USGS-"),
    begin_date = as.Date(begin),
    end_date   = as.Date(end),
    # Calendar span — note this is not the same as years with data,
    # but begin/end from this API reflect actual data extent, not just
    # gage establishment, so it's a reasonable proxy.
    record_span_yr = as.numeric(end_date - begin_date) / 365.25
  ) %>%
  # Keep one row per site (take the longest-running series if duplicates)
  group_by(site_no) %>%
  slice_max(record_span_yr, n = 1, with_ties = FALSE) %>%
  ungroup()

message(sprintf("Sites with daily mean discharge data: %d", nrow(daily_mean)))

# ============================================================================
# STEP 3: IDENTIFY SITES WITH INSTANTANEOUS DATA
# ============================================================================
# computation_identifier == "Instantaneous" flags IV series

iv_sites <- ts_meta_raw %>%
  filter(computation_identifier == "Instantaneous") %>%
  mutate(site_no = str_remove(monitoring_location_id, "^USGS-")) %>%
  distinct(site_no) %>%
  mutate(has_iv = TRUE)

message(sprintf("Sites with any instantaneous discharge data: %d", nrow(iv_sites)))

# ============================================================================
# STEP 4: GET SITE METADATA (NAME, COORDINATES, DRAINAGE AREA)
# ============================================================================
# Query by state + site_type — a supported server-side filter.
# Passing a large vector of monitoring_location_ids triggers HTTP 403
# (query complexity limit), so we fetch all CO stream sites and join instead.
message("Fetching Colorado stream gage site metadata...")

site_meta_sf <- read_waterdata_monitoring_location(
  state_name   = "Colorado",
  site_type    = "Stream",
  skipGeometry = FALSE,
  properties   = c(
    "monitoring_location_id",
    "monitoring_location_name",
    "drainage_area",
    "hydrologic_unit_code",
    "altitude"
  )
)

# Extract coordinates from sf geometry, then drop geometry for plain data frame
site_meta <- site_meta_sf %>%
  mutate(
    longitude = sf::st_coordinates(.)[, 1],
    latitude  = sf::st_coordinates(.)[, 2]
  ) %>%
  sf::st_drop_geometry() %>%
  mutate(site_no = str_remove(monitoring_location_id, "^USGS-"))

message(sprintf("CO stream gages with site metadata: %d", nrow(site_meta)))

# ============================================================================
# STEP 5: JOIN AND APPLY FILTERS
# ============================================================================
co_gages <- daily_mean %>%
  inner_join(site_meta, by = "site_no") %>%
  left_join(iv_sites, by = "site_no") %>%
  mutate(
    has_iv     = replace_na(has_iv, FALSE),
    has_recent = year(end_date) >= record_end_year,
    has_20yr   = record_span_yr >= min_years_short & has_recent,
    has_30yr   = record_span_yr >= min_years_long  & has_recent
  )

message(sprintf("\nAll CO stream gages with daily Q data:       %d", nrow(co_gages)))
message(sprintf("With data through %d:                    %d",
                record_end_year, sum(co_gages$has_recent)))
message(sprintf("With >= %d-yr span through %d:           %d",
                min_years_short, record_end_year, sum(co_gages$has_20yr)))
message(sprintf("With >= %d-yr span through %d:           %d",
                min_years_long,  record_end_year, sum(co_gages$has_30yr)))

# ============================================================================
# STEP 6: CROSS-REFERENCE AGAINST GAGES-II
# ============================================================================
message("\nLoading GAGES-II Colorado sites...")

gagesii_co <- read_csv(gagesii_meta, show_col_types = FALSE) %>%
  filter(STATE == "CO") %>%
  transmute(
    site_no    = str_pad(as.character(STAID), 8, "left", "0"),
    in_gagesii = TRUE
  )

message(sprintf("GAGES-II Colorado sites: %d", nrow(gagesii_co)))

co_gages <- co_gages %>%
  mutate(site_no = str_pad(site_no, 8, "left", "0")) %>%
  left_join(gagesii_co, by = "site_no") %>%
  mutate(
    in_gagesii = replace_na(in_gagesii, FALSE)
  )

# ============================================================================
# STEP 7: SUMMARY REPORT
# ============================================================================
message("\n=====================================================")
message("SUMMARY: GAGES-II vs Full Colorado USGS Inventory")
message("=====================================================\n")

summary_tbl <- co_gages %>%
  group_by(in_gagesii) %>%
  summarize(
    total_sites   = n(),
    recent_data   = sum(has_recent),
    ge_20yr       = sum(has_20yr),
    ge_30yr       = sum(has_30yr),
    with_iv_data  = sum(has_iv),
    ge_20yr_w_iv  = sum(has_20yr & has_iv),
    .groups = "drop"
  ) %>%
  mutate(in_gagesii = ifelse(in_gagesii, "In GAGES-II", "Not in GAGES-II"))

print(summary_tbl, width = Inf)

message("\n--- Non-GAGES-II sites with >= 20yr record through 2022 ---")
non_gagesii_viable <- co_gages %>%
  filter(!in_gagesii, has_20yr) %>%
  arrange(desc(record_span_yr)) %>%
  transmute(
    site_no,
    name            = monitoring_location_name,
    drainage_area_sqmi = drainage_area,
    record_start    = year(begin_date),
    record_end      = year(end_date),
    span_yr         = round(record_span_yr, 1),
    has_iv,
    huc             = hydrologic_unit_code
  )

print(non_gagesii_viable, n = Inf)

message(sprintf("\nTotal non-GAGES-II sites viable (>=20yr, through %d): %d",
                record_end_year, nrow(non_gagesii_viable)))
message(sprintf("Of which also have instantaneous data:                %d",
                sum(non_gagesii_viable$has_iv)))

# ============================================================================
# STEP 8: SAVE
# ============================================================================
write_csv(
  co_gages %>% filter(has_20yr) %>%
    select(site_no, monitoring_location_name, drainage_area,
           begin_date, end_date, record_span_yr,
           has_iv, has_20yr, has_30yr, in_gagesii,
           latitude, longitude, hydrologic_unit_code),
  output_file
)
message(sprintf("\nFull viability table saved to %s", output_file))