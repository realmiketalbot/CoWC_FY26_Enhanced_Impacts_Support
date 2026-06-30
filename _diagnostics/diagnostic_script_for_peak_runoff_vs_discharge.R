#!/usr/bin/env Rscript

# ============================================================================
# diagnostic_peak_runoff_vs_discharge.R
#
# Compares peak WLDAS total runoff (surface + subsurface) against peak USGS
# discharge depth for each classified POT event. Large mismatches — where
# observed discharge far exceeds what WLDAS produced — indicate events where
# WLDAS failed to capture the hydrometeorological forcing. This provides a
# data-driven alternative to the MIN_TOTAL_FLUX_MM noise floor in script 6b.
#
# Mismatch metric:
#   log10_ratio = log10( (usgs_peak_depth_mm + 1e-3) /
#                        (wldas_peak_runoff_mm + 1e-3) )
#   Positive:     discharge exceeds WLDAS runoff (expected: most peaks)
#   Large (> 2):  discharge ≥ 100× WLDAS runoff → strong signal miss
#
# WLDAS parquet organisation: one file per date, all sites.
#   Path: _data/WLDAS/processed/daily_by_date/YYYYMMDD_wldas.parquet
#   (Not per-site files — this is a common point of confusion.)
#
# INPUTS:
#   _data/WLDAS/processed/daily_by_date/  — per-date WLDAS parquet files
#   _data/pot_peaks_classified_adaptive.rds
#   _data/site_metadata.rds
#
# OUTPUTS (written to OUT_DIR):
#   fig01_log10_ratio_distribution.png
#   fig02_peak_scatter_log_log.png
#   fig03_site_outlier_rate_distribution.png
#   fig04_outlier_rate_vs_basin_chars.png
#   fig05_overlap_floor_vs_mismatch.png
# ============================================================================

library(tidyverse)
library(arrow)
library(scales)
library(patchwork)

# ============================================================================
# CONFIGURATION
# ============================================================================
# NOTE: WLDAS files are per-date, not per-site. Path corrected from the
# originally specified "_data/wldas" to the actual extraction output directory.
WLDAS_DIR         <- "_data/WLDAS/processed/daily_by_date"
PEAKS_FILE        <- "_data/pot_peaks_classified_adaptive.rds"
METADATA_FILE     <- "_data/site_metadata.rds"
OUT_DIR           <- "_diagnostics/peak_runoff_vs_discharge"

# Confirmed from 0_extract_site_metadata.R: DRAIN_SQKM → drain_sqkm (km²)
AREA_COL          <- "drain_sqkm"

# WLDAS runoff: no pre-combined column exists. Confirmed from
# 4_extract_wldas_daily_to_gages-ii_basin_mean.R variable table:
#   Qs_tavg  → qs_flux  (surface runoff,    mm/day after conversion)
#   Qsb_tavg → qsb_flux (subsurface runoff, mm/day after conversion)
# Combined total runoff = qs_flux + qsb_flux.
RUNOFF_SURFACE_COL    <- "qs_flux"
RUNOFF_SUBSURFACE_COL <- "qsb_flux"

# Unit conversion: cfs → mm/day over basin area.
# Formula (matches script 3_extract_usgs_unit_data_peaks.R exactly):
#   usgs_peak_depth_mm = peak_flow * (m³/s per cfs) * (s/day) /
#                        (drain_sqkm * km²→m²) * (mm/m)
#   = peak_flow * 0.0283168 * 86400 / (drain_sqkm * 1e6) * 1000
#   = peak_flow * CFS_TO_MM_PER_KM2 / drain_sqkm
# NOTE: the user-supplied constant (0.0283168 * 86400 / 1e6) was missing the
# ×1000 mm→m conversion; the correct scalar is 2.447, not 0.002447.
CFS_TO_MM_PER_KM2 <- 0.0283168 * 86400 * 1000 / 1e6  # ≈ 2.447

# Outlier criterion: discharge exceeds WLDAS runoff by a factor of 100
OUTLIER_LOG10     <- 2.0

# Flux floor from script 6b (for FIG 05 comparison)
MIN_TOTAL_FLUX_MM <- 1.0

SAVE_WIDTH        <- 10
SAVE_HEIGHT       <- 7
SAVE_DPI          <- 150

# ============================================================================
# SETUP
# ============================================================================
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(p, filename, width = SAVE_WIDTH, height = SAVE_HEIGHT) {
  path <- file.path(OUT_DIR, filename)
  ggsave(path, p, width = width, height = height, dpi = SAVE_DPI)
  message(sprintf("  Saved: %s", path))
  invisible(path)
}

theme_diag <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = NA),
      legend.position  = "bottom"
    )
}

# ============================================================================
# STEP 0 — SCHEMA INSPECTION
# ============================================================================
# This section inspects all three input sources at runtime and confirms that
# the column names hardcoded above exist. It prints the full schema to console
# so the output log captures any drift between expected and actual schemas.

message("\n=== STEP 0: SCHEMA INSPECTION ===\n")

## 0a — WLDAS parquet schema
parquet_files <- list.files(WLDAS_DIR, pattern = "\\.parquet$", full.names = TRUE)
if (length(parquet_files) == 0) {
  stop(sprintf(
    "No .parquet files found in '%s'. Check WLDAS_DIR path.", WLDAS_DIR
  ))
}

example_file <- parquet_files[1]
message(sprintf("Reading schema from: %s", basename(example_file)))

wldas_ds     <- arrow::open_dataset(example_file)
wldas_schema <- wldas_ds$schema
wldas_cols   <- names(wldas_schema)
wldas_types <- purrr::map_chr(wldas_schema$fields, ~ .x$type$ToString())

message(sprintf("WLDAS parquet columns (%d total):", length(wldas_cols)))
for (i in seq_along(wldas_cols)) {
  message(sprintf("  %-35s %s", wldas_cols[i], wldas_types[i]))
}

missing_runoff <- setdiff(c(RUNOFF_SURFACE_COL, RUNOFF_SUBSURFACE_COL), wldas_cols)
if (length(missing_runoff) > 0) {
  stop(sprintf(
    "Runoff columns not found in WLDAS parquet: %s\n",
    paste(missing_runoff, collapse = ", ")
  ))
}
message(sprintf(
  "\n  Runoff columns confirmed: %s (surface) + %s (subsurface)",
  RUNOFF_SURFACE_COL, RUNOFF_SUBSURFACE_COL
))

## 0b — POT peaks schema
message("\nReading POT peaks RDS...")
peaks_raw <- readRDS(PEAKS_FILE)
message(sprintf("  %d rows x %d columns, %d sites",
                nrow(peaks_raw), ncol(peaks_raw),
                dplyr::n_distinct(peaks_raw$site_no)))
message("  Columns:")
for (nm in names(peaks_raw)) {
  message(sprintf("    %-40s %s", nm, class(peaks_raw[[nm]])[1]))
}

required_peak_cols <- c(
  "site_no", "peak_id", "peak_date", "peak_flow", "search_radius",
  "flood_type", "flood_type_raw", "in_single_sample", "in_type_sample",
  "total_flux_short", "data_sufficient"
)
missing_peak_cols <- setdiff(required_peak_cols, names(peaks_raw))
if (length(missing_peak_cols) > 0) {
  warning(sprintf(
    "Missing expected columns in peaks RDS: %s",
    paste(missing_peak_cols, collapse = ", ")
  ))
} else {
  message("\n  All required peak columns confirmed.")
}

## 0c — Site metadata schema
message("\nReading site metadata RDS...")
meta_raw <- readRDS(METADATA_FILE)
message(sprintf("  %d sites, %d columns", nrow(meta_raw), ncol(meta_raw)))
message("  Columns:")
for (nm in names(meta_raw)) {
  message(sprintf("    %-40s %s", nm, class(meta_raw[[nm]])[1]))
}

if (!AREA_COL %in% names(meta_raw)) {
  stop(sprintf("AREA_COL '%s' not found in site metadata.", AREA_COL))
}
message(sprintf("\n  Drainage area column confirmed: %s (km²)", AREA_COL))

has_elevation <- "elev_mean_m" %in% names(meta_raw)
if (has_elevation) {
  message("  Elevation column confirmed: elev_mean_m (m)")
} else {
  message("  WARNING: elev_mean_m not in metadata — FIG 02/04 will use area-quartile fallback")
}

message("\n=== STEP 0 COMPLETE ===\n")

# ============================================================================
# LOAD AND PREPARE BASE DATA
# ============================================================================
message("Preparing peaks and site metadata...")

peaks <- peaks_raw %>%
  dplyr::mutate(
    peak_date = as.Date(peak_date),
    peak_id   = as.character(peak_id)
  )

meta <- meta_raw %>%
  dplyr::select(
    site_no,
    dplyr::all_of(AREA_COL),
    dplyr::any_of(c("elev_mean_m", "ndams_2009", "class"))
  )

peaks <- peaks %>%
  dplyr::left_join(meta, by = "site_no")

n_missing_area <- sum(is.na(peaks[[AREA_COL]]))
if (n_missing_area > 0) {
  warning(sprintf(
    "%d peaks have no drainage area — usgs_peak_depth_mm will be NA for those peaks",
    n_missing_area
  ))
}

# ============================================================================
# WLDAS DATE RANGE AND COVERAGE CHECK
# ============================================================================
message("Checking WLDAS archive coverage...")

wldas_dates_available <- parquet_files %>%
  basename() %>%
  stringr::str_extract("\\d{8}") %>%
  as.Date(format = "%Y%m%d") %>%
  sort()

message(sprintf("  WLDAS archive: %s → %s (%d daily files)",
                min(wldas_dates_available),
                max(wldas_dates_available),
                length(wldas_dates_available)))

# Antecedent window per peak: [peak_date - (search_radius - 1), peak_date]
peaks <- peaks %>%
  dplyr::mutate(
    window_start    = peak_date - (as.integer(search_radius) - 1L),
    window_end      = peak_date,
    wldas_available = (window_start >= min(wldas_dates_available)) &
                      (window_end   <= max(wldas_dates_available))
  )

n_no_wldas  <- sum(!peaks$wldas_available)
n_has_wldas <- sum( peaks$wldas_available)
message(sprintf("  Peaks with full WLDAS window:    %d", n_has_wldas))
message(sprintf("  Peaks outside WLDAS archive:     %d", n_no_wldas))

# ============================================================================
# LOAD WLDAS RUNOFF DATA
# ============================================================================
# Strategy: collect all unique dates needed, load each file once, bind rows.
# Files are per-date (not per-site), so we load by date and join by site.

needed_dates <- peaks %>%
  dplyr::filter(wldas_available) %>%
  dplyr::rowwise() %>%
  dplyr::reframe(date = seq(window_start, window_end, by = "day")) %>%
  dplyr::distinct(date) %>%
  dplyr::filter(date %in% wldas_dates_available) %>%
  dplyr::pull(date)

message(sprintf("\nUnique WLDAS dates to load: %d", length(needed_dates)))
message("Loading WLDAS runoff data (may take several minutes)...")

needed_files <- file.path(
  WLDAS_DIR,
  sprintf("%s_wldas.parquet", format(needed_dates, "%Y%m%d"))
)

wldas_runoff <- purrr::map_dfr(
  seq_along(needed_files),
  function(i) {
    if (i %% 200 == 0 || i == length(needed_files)) {
      message(sprintf("  Loading file %d / %d...", i, length(needed_files)))
    }
    f <- needed_files[i]
    tryCatch(
      arrow::read_parquet(
        f,
        col_select = c("date", "GAGE_ID", RUNOFF_SURFACE_COL, RUNOFF_SUBSURFACE_COL)
      ) %>%
        dplyr::mutate(
          date         = as.Date(date),
          site_no      = as.character(GAGE_ID),
          total_runoff = .data[[RUNOFF_SURFACE_COL]] + .data[[RUNOFF_SUBSURFACE_COL]]
        ) %>%
        dplyr::select(date, site_no, total_runoff),
      error = function(e) {
        warning(sprintf("Could not read %s: %s", basename(f), e$message))
        NULL
      }
    )
  }
)

message(sprintf("  Loaded %d site-day rows across %d dates",
                nrow(wldas_runoff),
                dplyr::n_distinct(wldas_runoff$date)))

# ============================================================================
# COMPUTE wldas_peak_runoff_mm PER PEAK
# ============================================================================
# Expand each peak to its antecedent window, join WLDAS, take window maximum.

message("\nComputing wldas_peak_runoff_mm per peak...")

peak_windows <- peaks %>%
  dplyr::filter(wldas_available) %>%
  dplyr::select(site_no, peak_id, window_start, window_end) %>%
  dplyr::rowwise() %>%
  dplyr::reframe(
    site_no  = site_no,
    peak_id  = peak_id,
    date     = seq(window_start, window_end, by = "day")
  )

wldas_peak_max <- peak_windows %>%
  dplyr::left_join(wldas_runoff, by = c("site_no", "date")) %>%
  dplyr::group_by(site_no, peak_id) %>%
  dplyr::summarise(
    wldas_peak_runoff_mm = max(total_runoff, na.rm = TRUE),
    n_wldas_days         = sum(!is.na(total_runoff)),
    .groups              = "drop"
  ) %>%
  dplyr::mutate(
    # max() of all-NA returns -Inf; convert to NA
    wldas_peak_runoff_mm = dplyr::if_else(
      is.infinite(wldas_peak_runoff_mm), NA_real_, wldas_peak_runoff_mm
    )
  )

message(sprintf("  Peak runoff extracted for %d peaks", nrow(wldas_peak_max)))

# ============================================================================
# ASSEMBLE ANALYSIS TABLE
# ============================================================================
message("\nAssembling analysis table...")

analysis <- peaks %>%
  dplyr::left_join(wldas_peak_max, by = c("site_no", "peak_id")) %>%
  dplyr::mutate(
    usgs_peak_depth_mm   = peak_flow * CFS_TO_MM_PER_KM2 / .data[[AREA_COL]],
    log10_ratio          = log10(
      (usgs_peak_depth_mm  + 1e-3) /
      (wldas_peak_runoff_mm + 1e-3)
    ),
    is_outlier           = !is.na(log10_ratio) & log10_ratio > OUTLIER_LOG10,
    flux_floor_fail      = !is.na(total_flux_short) &
                           total_flux_short < MIN_TOTAL_FLUX_MM,
    flood_type_raw       = dplyr::coalesce(as.character(flood_type_raw), "unclassified"),
    flood_type_raw       = factor(
      flood_type_raw,
      levels = c("rainfall", "ROS", "snowmelt", "unclassified")
    )
  )

n_complete_wldas <- sum(!is.na(analysis$wldas_peak_runoff_mm))
n_missing_wldas  <- sum( is.na(analysis$wldas_peak_runoff_mm))

# ============================================================================
# CONSOLE SUMMARY
# ============================================================================
message("\n==============================")
message("  SUMMARY")
message("==============================")
message(sprintf("Total peaks:                    %d", nrow(analysis)))
message(sprintf("Peaks with WLDAS coverage:      %d (%.1f%%)",
                n_complete_wldas, 100 * n_complete_wldas / nrow(analysis)))
message(sprintf("Peaks missing WLDAS data:       %d (%.1f%%)",
                n_missing_wldas,  100 * n_missing_wldas  / nrow(analysis)))

analysis_complete <- analysis %>% dplyr::filter(!is.na(log10_ratio))

q_overall <- quantile(analysis_complete$log10_ratio, c(0.25, 0.5, 0.75), na.rm = TRUE)
message(sprintf(
  "\nlog10_ratio overall: median = %.2f  |  IQR = [%.2f, %.2f]",
  q_overall[2], q_overall[1], q_overall[3]
))

by_type_summary <- analysis_complete %>%
  dplyr::group_by(flood_type_raw) %>%
  dplyr::summarise(
    n             = dplyr::n(),
    median_ratio  = median(log10_ratio),
    q25_ratio     = quantile(log10_ratio, 0.25),
    q75_ratio     = quantile(log10_ratio, 0.75),
    n_outlier     = sum(log10_ratio > OUTLIER_LOG10),
    pct_outlier   = 100 * mean(log10_ratio > OUTLIER_LOG10),
    .groups       = "drop"
  )

message("\nBy flood_type_raw:")
print(as.data.frame(by_type_summary), row.names = FALSE)

n_outlier_total  <- sum(analysis_complete$log10_ratio > OUTLIER_LOG10)
pct_outlier_total <- 100 * n_outlier_total / nrow(analysis_complete)
message(sprintf(
  "\nTotal peaks with log10_ratio > %.0f: %d (%.1f%%)",
  OUTLIER_LOG10, n_outlier_total, pct_outlier_total
))

site_outlier_rates <- analysis_complete %>%
  dplyr::group_by(site_no) %>%
  dplyr::summarise(
    outlier_rate = mean(log10_ratio > OUTLIER_LOG10),
    n_peaks      = dplyr::n(),
    .groups      = "drop"
  )

message(sprintf("Sites with outlier rate > 10%%:  %d", sum(site_outlier_rates$outlier_rate > 0.10)))
message(sprintf("Sites with outlier rate > 25%%:  %d", sum(site_outlier_rates$outlier_rate > 0.25)))
message(sprintf("Sites with outlier rate > 50%%:  %d", sum(site_outlier_rates$outlier_rate > 0.50)))

# ============================================================================
# SHARED STYLE ELEMENTS
# ============================================================================
type_colors <- c(
  "rainfall"     = "#E69F00",
  "ROS"          = "#56B4E9",
  "snowmelt"     = "#009E73",
  "unclassified" = "grey60"
)

# ============================================================================
# FIG 01 — LOG10_RATIO DISTRIBUTION
# ECDF + histogram by flood_type_raw; vertical lines at 0 and OUTLIER_LOG10.
# ============================================================================
message("\nBuilding Fig 01 — log10_ratio distribution...")

fig01_hist <- ggplot(analysis_complete,
                     aes(x = log10_ratio, fill = flood_type_raw)) +
  geom_histogram(binwidth = 0.25, position = "identity", alpha = 0.5, colour = NA) +
  geom_vline(xintercept = 0,              linetype = "dashed", colour = "grey30",
             linewidth = 0.7) +
  geom_vline(xintercept = OUTLIER_LOG10,  linetype = "dashed", colour = "firebrick",
             linewidth = 0.7) +
  scale_fill_manual(values = type_colors, name = "Flood type") +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~ flood_type_raw, scales = "free_y", ncol = 2) +
  labs(
    x       = expression(log[10]*"(USGS discharge depth / WLDAS runoff)"),
    y       = "Count",
    title   = "Fig 01a — Discharge–runoff mismatch distribution by flood type",
    subtitle = sprintf(
      "Dashed grey = perfect match (0); dashed red = outlier threshold (log₁₀ > %.0f)",
      OUTLIER_LOG10
    )
  ) +
  theme_diag() +
  theme(legend.position = "none")

fig01_ecdf <- ggplot(analysis_complete,
                     aes(x = log10_ratio, colour = flood_type_raw)) +
  stat_ecdf(geom = "step", linewidth = 0.9) +
  geom_vline(xintercept = 0,              linetype = "dashed", colour = "grey30",
             linewidth = 0.7) +
  geom_vline(xintercept = OUTLIER_LOG10,  linetype = "dashed", colour = "firebrick",
             linewidth = 0.7) +
  scale_colour_manual(values = type_colors, name = "Flood type") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x       = expression(log[10]*"(USGS discharge depth / WLDAS runoff)"),
    y       = "Cumulative fraction",
    title   = "Fig 01b — ECDF of mismatch metric by flood type"
  ) +
  theme_diag()

fig01 <- fig01_hist / fig01_ecdf + patchwork::plot_layout(heights = c(2.5, 1))
save_fig(fig01, "fig01_log10_ratio_distribution.png", height = 11)

# ============================================================================
# FIG 02 — SCATTER: usgs_peak_depth_mm vs wldas_peak_runoff_mm (log-log)
# Faceted by mean basin elevation bin (or drainage area quartile if unavailable).
# ============================================================================
message("Building Fig 02 — log-log scatter...")

fig02_base <- analysis_complete %>%
  dplyr::filter(usgs_peak_depth_mm > 0,
                !is.na(wldas_peak_runoff_mm),
                wldas_peak_runoff_mm >= 0)

if (has_elevation) {
  fig02_data <- fig02_base %>%
    dplyr::mutate(
      facet_var = cut(
        elev_mean_m,
        breaks = c(-Inf, 1500, 2500, 3500, Inf),
        labels = c("< 1500 m", "1500–2500 m", "2500–3500 m", "> 3500 m"),
        right  = TRUE
      )
    )
  facet_label <- "Mean basin elevation"
} else {
  fig02_data <- fig02_base %>%
    dplyr::mutate(
      facet_var = cut(
        .data[[AREA_COL]],
        breaks         = quantile(.data[[AREA_COL]], c(0, .25, .5, .75, 1), na.rm = TRUE),
        labels         = c("Q1 (smallest)", "Q2", "Q3", "Q4 (largest)"),
        include.lowest = TRUE
      )
    )
  facet_label <- "Drainage area quartile"
}

xy_lims <- range(
  c(fig02_data$usgs_peak_depth_mm,
    fig02_data$wldas_peak_runoff_mm + 1e-3),
  na.rm = TRUE
)

fig02 <- ggplot(
  fig02_data,
  aes(x = wldas_peak_runoff_mm + 1e-3,
      y = usgs_peak_depth_mm,
      colour = flood_type_raw)
) +
  geom_point(alpha = 0.25, size = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey30", linewidth = 0.7) +
  scale_x_log10(labels = scales::comma, limits = xy_lims) +
  scale_y_log10(labels = scales::comma, limits = xy_lims) +
  scale_colour_manual(values = type_colors, name = "Flood type") +
  facet_wrap(~ facet_var, ncol = 2) +
  labs(
    x        = "WLDAS peak total runoff (mm/day, log scale)",
    y        = "USGS peak discharge depth (mm/day, log scale)",
    title    = "Fig 02 — Peak discharge depth vs. WLDAS peak runoff",
    subtitle = sprintf("1:1 line dashed; faceted by %s", facet_label),
    caption  = "Points above the 1:1 line: discharge exceeds WLDAS runoff"
  ) +
  theme_diag()

save_fig(fig02, "fig02_peak_scatter_log_log.png")

# ============================================================================
# FIG 03 — SITE-LEVEL OUTLIER RATE: histogram + ECDF
# Overlay with flux-floor failure rate for direct comparison.
# ============================================================================
message("Building Fig 03 — site outlier rate distribution...")

site_flux_fail_rate <- analysis_complete %>%
  dplyr::group_by(site_no) %>%
  dplyr::summarise(
    flux_fail_rate = mean(flux_floor_fail, na.rm = TRUE),
    .groups        = "drop"
  )

site_diag <- site_outlier_rates %>%
  dplyr::left_join(site_flux_fail_rate, by = "site_no")

criterion_labels <- c(
  outlier_rate   = sprintf("Runoff mismatch  (log₁₀ ratio > %.0f)", OUTLIER_LOG10),
  flux_fail_rate = sprintf("Flux floor fail  (total flux < %.1f mm)", MIN_TOTAL_FLUX_MM)
)
criterion_colors <- c(
  "outlier_rate"   = "#D55E00",
  "flux_fail_rate" = "#0072B2"
)

long_rates <- site_diag %>%
  tidyr::pivot_longer(
    c(outlier_rate, flux_fail_rate),
    names_to  = "criterion",
    values_to = "rate"
  ) %>%
  dplyr::mutate(criterion = criterion_labels[criterion])

fig03_hist <- ggplot(long_rates, aes(x = rate * 100, fill = criterion)) +
  geom_histogram(binwidth = 5, position = "identity", alpha = 0.6, colour = NA) +
  scale_fill_manual(values = unname(criterion_colors), name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x     = "Peak failure / outlier rate per site (%)",
    y     = "Number of sites",
    title = "Fig 03a — Site-level outlier rate: runoff mismatch vs. flux floor"
  ) +
  theme_diag()

fig03_ecdf <- ggplot(long_rates, aes(x = rate, colour = criterion)) +
  stat_ecdf(geom = "step", linewidth = 0.9) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  scale_colour_manual(values = unname(criterion_colors), name = NULL) +
  labs(
    x     = "Per-site failure rate",
    y     = "Cumulative fraction of sites",
    title = "Fig 03b — ECDF of per-site rates"
  ) +
  theme_diag()

fig03 <- fig03_hist / fig03_ecdf
save_fig(fig03, "fig03_site_outlier_rate_distribution.png", height = 9)

# ============================================================================
# FIG 04 — OUTLIER RATE VS. BASIN CHARACTERISTICS
# Per-site outlier rate against mean basin elevation and drainage area.
# ============================================================================
message("Building Fig 04 — outlier rate vs. basin characteristics...")

site_chars <- site_outlier_rates %>%
  dplyr::left_join(meta, by = "site_no")

if (has_elevation && any(!is.na(site_chars$elev_mean_m))) {
  fig04a <- ggplot(
    site_chars %>% dplyr::filter(!is.na(elev_mean_m)),
    aes(x = elev_mean_m, y = outlier_rate)
  ) +
    geom_point(alpha = 0.5, colour = "#D55E00") +
    geom_smooth(method  = "loess", se = TRUE,
                colour  = "#D55E00", fill = "#D55E00",
                alpha   = 0.2, linewidth = 0.8) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      x     = "Mean basin elevation (m)",
      y     = sprintf("Outlier rate (log₁₀ ratio > %.0f)", OUTLIER_LOG10),
      title = "Fig 04a — Outlier rate vs. mean basin elevation"
    ) +
    theme_diag()
} else {
  fig04a <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "elev_mean_m not available in site metadata",
             colour = "grey50", size = 5) +
    theme_void()
}

fig04b <- ggplot(
  site_chars %>% dplyr::filter(!is.na(.data[[AREA_COL]])),
  aes(x = .data[[AREA_COL]], y = outlier_rate)
) +
  geom_point(alpha = 0.5, colour = "#0072B2") +
  geom_smooth(method  = "loess", se = TRUE,
              colour  = "#0072B2", fill = "#0072B2",
              alpha   = 0.2, linewidth = 0.8) +
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x     = "Drainage area (km², log scale)",
    y     = sprintf("Outlier rate (log₁₀ ratio > %.0f)", OUTLIER_LOG10),
    title = "Fig 04b — Outlier rate vs. drainage area"
  ) +
  theme_diag()

fig04 <- fig04a + fig04b
save_fig(fig04, "fig04_outlier_rate_vs_basin_chars.png", height = 5)

# ============================================================================
# FIG 05 — OVERLAP: FLUX FLOOR vs. RUNOFF MISMATCH
# Stacked bar by flood_type_raw showing which peaks are flagged by each
# criterion alone, both, or neither.
# ============================================================================
message("Building Fig 05 — overlap between flux floor and runoff mismatch flags...")

overlap_levels <- c(
  "Neither (passes both)",
  "Flux floor only",
  "Runoff mismatch only",
  "Both flags"
)
overlap_colors <- c(
  "Neither (passes both)" = "grey70",
  "Flux floor only"       = "#0072B2",
  "Runoff mismatch only"  = "#D55E00",
  "Both flags"            = "#CC79A7"
)

fig05_data <- analysis_complete %>%
  dplyr::filter(!is.na(flux_floor_fail)) %>%
  dplyr::mutate(
    overlap_group = dplyr::case_when(
      !flux_floor_fail & !is_outlier ~ "Neither (passes both)",
      flux_floor_fail  & !is_outlier ~ "Flux floor only",
      !flux_floor_fail &  is_outlier ~ "Runoff mismatch only",
      flux_floor_fail  &  is_outlier ~ "Both flags"
    ),
    overlap_group = factor(overlap_group, levels = overlap_levels)
  )

fig05 <- ggplot(fig05_data,
                aes(x = flood_type_raw, fill = overlap_group)) +
  geom_bar(position = "stack", colour = NA) +
  scale_fill_manual(values = overlap_colors, name = "Flag status") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x        = "Flood type",
    y        = "Number of peaks",
    title    = "Fig 05 — Overlap: flux floor failures vs. runoff–discharge mismatch",
    subtitle = sprintf(
      "Flux floor: total_flux_short < %.1f mm  |  Mismatch: log₁₀(USGS/WLDAS) > %.0f",
      MIN_TOTAL_FLUX_MM, OUTLIER_LOG10
    ),
    caption  = "Restricted to peaks with complete WLDAS coverage and non-NA log10_ratio"
  ) +
  theme_diag()

save_fig(fig05, "fig05_overlap_floor_vs_mismatch.png", height = 6)

# ============================================================================
# FINAL OVERLAP SUMMARY TABLE
# ============================================================================
message("\n=== OVERLAP SUMMARY TABLE (FIG 05) ===")

overlap_wide <- fig05_data %>%
  dplyr::count(flood_type_raw, overlap_group) %>%
  tidyr::pivot_wider(
    names_from  = overlap_group,
    values_from = n,
    values_fill = 0L
  )
print(as.data.frame(overlap_wide), row.names = FALSE)

margin_totals <- fig05_data %>%
  dplyr::count(overlap_group) %>%
  tidyr::pivot_wider(
    names_from  = overlap_group,
    values_from = n,
    values_fill = 0L
  ) %>%
  dplyr::mutate(flood_type_raw = factor("TOTAL"))

message("\nMargin totals across all flood types:")
print(as.data.frame(margin_totals %>%
  dplyr::select(flood_type_raw, dplyr::everything())), row.names = FALSE)

message("\n=== DIAGNOSTIC COMPLETE ===")
message(sprintf("All figures written to: %s/", OUT_DIR))
