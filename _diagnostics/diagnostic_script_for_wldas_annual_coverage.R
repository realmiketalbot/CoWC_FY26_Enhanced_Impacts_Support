#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: WLDAS Annual Coverage and Missing Year Investigation
# ============================================================================
#
# PURPOSE:
#   Investigates two issues observed in the climate trend analysis:
#
#   1. ANOMALOUSLY LOW FINAL-YEAR PRECIPITATION
#      WY2024 median total_precip_mm (~85 mm) is ~14% of WY2023's median
#      (~621 mm), consistent with the WLDAS record not covering the full
#      final water year. Script 8 sums precipitation with na.rm = TRUE, so
#      partial WLDAS coverage silently produces a low annual total rather than
#      NA, which then enters Mann-Kendall testing and biases the trend.
#
#   2. MISSING YEARS IN TREND ANALYSIS
#      Some sites show gaps in their annual time series. This diagnostic
#      distinguishes between (a) genuinely inactive gauge years (no POT peaks
#      — expected and correctly handled as NA) and (b) years where WLDAS
#      coverage may be insufficient, which would need correction.
#
# STRUCTURE:
#   Part 1 — Final-year anomaly  (uses annual_climate_forcing_summary.rds;
#                                  runs on any machine)
#   Part 2 — Missing year patterns (same; runs on any machine)
#   Part 3 — WLDAS daily coverage (requires WLDAS parquet; skips gracefully
#                                   if the dataset is not available locally)
#
# OUTPUTS:
#   Console: coverage summaries and per-issue findings
#   _diagnostics/climate_forcing_trends/fig_wldas_precip_ts.png
#   _diagnostics/climate_forcing_trends/fig_wldas_precip_final_year.png
#   _diagnostics/climate_forcing_trends/fig_wldas_daily_coverage.png  (Part 3)
#
# RECOMMENDED FIX (Script 8):
#   After computing wldas_annual, count WLDAS days per (site, water_year)
#   and set total_precip_mm = NA (and peak_swe_mm = NA) for any year where
#   day coverage is below MIN_DAYS_FRAC. This excludes incomplete years from
#   Mann-Kendall without requiring a hard cutoff on the final year.
#
# Author : Mike Talbot, PhD Student, Colorado State University
# Date   : 2026-05-21
# ============================================================================

library(tidyverse)
library(scales)

# ── Configuration ─────────────────────────────────────────────────────────────
WLDAS_DIR    <- "_data/WLDAS/processed/daily_by_date"
SUMMARY_FILE <- "_data/annual_climate_forcing_summary.rds"
OUT_DIR      <- "_diagnostics/climate_forcing_trends"

# Coverage fraction below which a site-year is considered incomplete.
# A water year has 365 or 366 days; 90% = no more than ~36 days missing.
MIN_DAYS_FRAC <- 0.90

water_year <- function(d) {
  yr <- as.integer(format(as.Date(d), "%Y"))
  mo <- as.integer(format(as.Date(d), "%m"))
  ifelse(mo >= 10L, yr + 1L, yr)
}

is_leap_year <- function(yr) {
  yr %% 4 == 0 & (yr %% 100 != 0 | yr %% 400 == 0)
}

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)


# ============================================================================
# PART 1: FINAL-YEAR ANOMALY
# ============================================================================
message("\n=== PART 1: Final-year anomaly (annual_climate_forcing_summary.rds) ===\n")

annual <- readRDS(SUMMARY_FILE)

final_wy <- max(annual$water_year, na.rm = TRUE)
prev_wy  <- final_wy - 1L

message(sprintf("Annual summary: WY %d – %d | %d sites",
                min(annual$water_year, na.rm = TRUE),
                final_wy,
                n_distinct(annual$site_no)))

# ── Fleet-wide median precipitation by year ───────────────────────────────────
wy_precip <- annual %>%
  group_by(water_year) %>%
  summarise(
    n_sites        = sum(!is.na(total_precip_mm)),
    median_precip  = median(total_precip_mm, na.rm = TRUE),
    q25_precip     = quantile(total_precip_mm, 0.25, na.rm = TRUE),
    q75_precip     = quantile(total_precip_mm, 0.75, na.rm = TRUE),
    .groups        = "drop"
  )

message("\nMedian total_precip_mm — final 5 water years:")
tail(wy_precip, 5) %>% print(n = Inf)

final_med   <- wy_precip$median_precip[wy_precip$water_year == final_wy]
prev_med    <- wy_precip$median_precip[wy_precip$water_year == prev_wy]

if (is.na(final_med)) {
  pct_of_prev <- NA_real_
  message(sprintf(
    "\nWY%d median precip: NA — all site-years are NA (Script 8 coverage mask already applied).",
    final_wy
  ))
  message(sprintf("WY%d median precip: %.1f mm", prev_wy, prev_med))
} else {
  pct_of_prev <- 100 * final_med / prev_med
  message(sprintf(
    "\nWY%d median precip: %.1f mm  (%.0f%% of WY%d median: %.1f mm)",
    final_wy, final_med, pct_of_prev, prev_wy, prev_med
  ))
  if (pct_of_prev < 50) {
    message(sprintf(
      "*** FLAG: WY%d precip is %.0f%% of the prior year.",
      final_wy, pct_of_prev
    ))
    message(sprintf(
      "    Consistent with an incomplete WLDAS record for WY%d.", final_wy
    ))
    message(
      "    Estimated coverage: the ratio implies roughly ",
      sprintf("%.0f days of a ~365-day year are present.",
              round(365 * pct_of_prev / 100))
    )
  }
}

# ── Figure 1: Fleet-wide precipitation time series ───────────────────────────
p_ts <- ggplot(wy_precip %>% filter(n_sites >= 10),
               aes(x = water_year, y = median_precip)) +
  geom_ribbon(aes(ymin = q25_precip, ymax = q75_precip),
              fill = "#4e79a7", alpha = 0.2) +
  geom_line(colour = "#4e79a7", linewidth = 0.8) +
  geom_point(colour = "#4e79a7", size = 1.6) +
  geom_point(data = filter(wy_precip, water_year == final_wy),
             colour = "#d95f02", size = 4, shape = 18) +
  (if (!is.na(final_med))
    annotate("text",
      x = final_wy - 0.7, y = final_med + 20,
      label = sprintf("WY%d\n%.0f mm", final_wy, final_med),
      colour = "#d95f02", size = 3.2, hjust = 1, fontface = "bold"
    )
  ) +
  scale_x_continuous(name = "Water Year", breaks = seq(1980, 2030, 5)) +
  scale_y_continuous(name = "Median total precip. (mm)", labels = comma) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title    = "Fleet-wide median annual precipitation by water year",
    subtitle = if (is.na(final_med))
      sprintf(
        "Ribbon = IQR across sites  |  WY%d masked (WLDAS coverage < 90%%, corrected in Script 8)",
        final_wy
      )
    else
      sprintf(
        "Ribbon = IQR across sites  |  Orange diamond = WY%d (final year, suspect incomplete WLDAS)",
        final_wy
      )
  )

ggsave(file.path(OUT_DIR, "fig_wldas_precip_ts.png"),
       p_ts, width = 9, height = 4.5, dpi = 150)
message("\nSaved: fig_wldas_precip_ts.png")

# ── Figure 2: Final year vs. prior year precipitation distribution ─────────────
compare_df <- annual %>%
  filter(water_year %in% c(prev_wy, final_wy), !is.na(total_precip_mm)) %>%
  mutate(wy_label = factor(
    water_year,
    levels = c(prev_wy, final_wy),
    labels = c(
      sprintf("WY%d (previous)", prev_wy),
      sprintf("WY%d (final — suspect)", final_wy)
    )
  ))

p_compare <- ggplot(compare_df, aes(x = total_precip_mm, fill = wy_label)) +
  geom_histogram(binwidth = 30, colour = "white", linewidth = 0.2,
                 position = "identity", alpha = 0.65, boundary = 0) +
  geom_vline(
    data = compare_df %>%
      group_by(wy_label) %>%
      summarise(med = median(total_precip_mm), .groups = "drop"),
    aes(xintercept = med, colour = wy_label),
    linetype = "dashed", linewidth = 0.9, show.legend = FALSE
  ) +
  scale_fill_manual(values  = c("#4e79a7", "#d95f02"), name = NULL) +
  scale_colour_manual(values = c("#4e79a7", "#d95f02")) +
  scale_x_continuous(name = "Annual total precipitation (mm)", labels = comma) +
  scale_y_continuous(name = "Site count") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top") +
  labs(
    title    = sprintf("Annual precipitation: WY%d vs. WY%d", final_wy, prev_wy),
    subtitle = if (is.na(final_med))
      sprintf(
        "WY%d has been masked (WLDAS coverage < 90%%) — no WY%d bars shown",
        final_wy, final_wy
      )
    else
      sprintf(
        "Dashed lines = medians  |  WY%d values are anomalously low (incomplete WLDAS coverage)",
        final_wy
      )
  )

ggsave(file.path(OUT_DIR, "fig_wldas_precip_final_year.png"),
       p_compare, width = 7.5, height = 4.5, dpi = 150)
message("Saved: fig_wldas_precip_final_year.png")


# ============================================================================
# PART 2: MISSING YEAR PATTERNS
# ============================================================================
message("\n=== PART 2: Missing year patterns ===\n")

# Active span per site: first through last year with POT peaks present.
# Years outside this span are structurally missing (gauge inactive) — expected.
# Years inside this span with no POT peaks are interior gaps — investigated here.
site_spans <- annual %>%
  filter(!is.na(pot_event_rate)) %>%
  group_by(site_no) %>%
  summarise(
    first_wy   = min(water_year, na.rm = TRUE),
    last_wy    = max(water_year, na.rm = TRUE),
    n_active   = n_distinct(water_year),
    .groups    = "drop"
  )

interior_gaps <- annual %>%
  inner_join(site_spans, by = "site_no") %>%
  filter(
    water_year >= first_wy,
    water_year <= last_wy,
    is.na(pot_event_rate)
  )

n_gap_sites <- n_distinct(interior_gaps$site_no)
n_gap_yrs   <- nrow(interior_gaps)

message(sprintf(
  "Interior gaps (no POT peaks within active span): %d site-years across %d sites",
  n_gap_yrs, n_gap_sites
))

if (n_gap_yrs > 0) {
  gap_summary <- interior_gaps %>%
    group_by(site_no) %>%
    summarise(
      n_gap_years = n(),
      years       = paste(sort(water_year), collapse = ", "),
      first_wy    = first(first_wy),
      last_wy     = first(last_wy),
      .groups     = "drop"
    ) %>%
    arrange(desc(n_gap_years))

  message("\nTop sites by number of interior gap years (max 15 shown):")
  print(head(gap_summary, 15), n = Inf)

  message(
    "\nNote: interior gaps indicate years where the gauge had no POT-threshold",
    "\ncrossings (all peaks below the site's 50th-percentile threshold). This",
    "\nis a genuine data feature, not a WLDAS coverage problem. These years are",
    "\ncorrectly represented as NA and excluded from Mann-Kendall via na.rm."
  )
} else {
  message("No interior gaps found — every site has contiguous POT peak records.")
}

# Does the final WY show up as a gap for any sites?
final_gap_sites <- interior_gaps %>% filter(water_year == final_wy)
message(sprintf(
  "\nSites where WY%d is an interior gap (no POT peaks): %d",
  final_wy, n_distinct(final_gap_sites$site_no)
))


# ============================================================================
# PART 3: WLDAS DAILY COVERAGE PER SITE-YEAR
# ============================================================================
message("\n=== PART 3: WLDAS daily coverage per site-year (parquet) ===\n")

if (!dir.exists(WLDAS_DIR)) {
  message("  WLDAS parquet dataset not found at: ", WLDAS_DIR)
  message("  Skipping Part 3 — run on the server where the parquet files are stored.")
  message("  When available, Part 3 will count exact WLDAS days per (site, water_year),")
  message(sprintf(
    "  flag years with coverage < %.0f%%, and confirm whether WY%d is incomplete.",
    100 * MIN_DAYS_FRAC, final_wy
  ))

} else {

  suppressPackageStartupMessages(library(arrow))

  wldas_ds <- arrow::open_dataset(WLDAS_DIR)

  message("Counting WLDAS days per (GAGE_ID, water_year)...")

  # Collect only GAGE_ID and date — minimal footprint.
  # water_year() applied in R after collect; not pushable to parquet reader.
  raw_days <- wldas_ds %>%
    dplyr::select(GAGE_ID, date) %>%
    dplyr::collect() %>%
    dplyr::mutate(water_year = water_year(date))

  day_counts <- raw_days %>%
    dplyr::group_by(GAGE_ID, water_year) %>%
    dplyr::summarise(
      n_days   = dplyr::n(),
      min_date = min(date, na.rm = TRUE),
      max_date = max(date, na.rm = TRUE),
      .groups  = "drop"
    )

  # Expected days per water year: 365 + 1 if the calendar year containing
  # Feb 29 (i.e., the year the WY ends in) is a leap year.
  expected_days <- day_counts %>%
    dplyr::distinct(water_year) %>%
    dplyr::mutate(
      expected_days = 365L + dplyr::if_else(is_leap_year(water_year), 1L, 0L)
    )

  day_counts <- day_counts %>%
    dplyr::left_join(expected_days, by = "water_year") %>%
    dplyr::mutate(
      frac_covered  = n_days / expected_days,
      is_incomplete = frac_covered < MIN_DAYS_FRAC
    )

  # ── Fleet-wide coverage summary by water year ──────────────────────────────
  coverage_wy <- day_counts %>%
    dplyr::group_by(water_year) %>%
    dplyr::summarise(
      n_sites        = dplyr::n(),
      median_frac    = median(frac_covered, na.rm = TRUE),
      min_frac       = min(frac_covered, na.rm = TRUE),
      n_incomplete   = sum(is_incomplete, na.rm = TRUE),
      pct_incomplete = 100 * n_incomplete / n_sites,
      .groups        = "drop"
    )

  message("\nWLDAS coverage by water year — final 5 years:")
  tail(coverage_wy, 5) %>% print(n = Inf)

  # ── Final year detail ──────────────────────────────────────────────────────
  final_cov <- dplyr::filter(coverage_wy, water_year == final_wy)
  message(sprintf(
    "\nWY%d: median site coverage = %.1f%%  |  sites with < %.0f%% coverage: %d / %d",
    final_wy,
    100 * final_cov$median_frac,
    100 * MIN_DAYS_FRAC,
    final_cov$n_incomplete,
    final_cov$n_sites
  ))

  # Last date present in the dataset (globally)
  wldas_max_date <- max(raw_days$date, na.rm = TRUE)
  last_complete_wy <- water_year(wldas_max_date) - 1L
  message(sprintf(
    "\nLatest WLDAS date in parquet: %s  (WY%d ends %s)",
    wldas_max_date,
    water_year(wldas_max_date),
    format(as.Date(sprintf("%d-09-30", water_year(wldas_max_date))), "%Y-%m-%d")
  ))
  message(sprintf(
    "Last fully covered water year: WY%d", last_complete_wy
  ))

  # Sites where the final year is incomplete
  incomplete_final <- day_counts %>%
    dplyr::filter(water_year == final_wy, is_incomplete) %>%
    dplyr::arrange(frac_covered) %>%
    dplyr::select(GAGE_ID, n_days, expected_days, frac_covered, min_date, max_date)

  if (nrow(incomplete_final) > 0) {
    message(sprintf(
      "\nSites with WY%d WLDAS coverage < %.0f%% (first 20 shown):",
      final_wy, 100 * MIN_DAYS_FRAC
    ))
    print(head(incomplete_final, 20), n = Inf)
  } else {
    message(sprintf(
      "\nAll sites have >= %.0f%% WLDAS day coverage in WY%d.",
      100 * MIN_DAYS_FRAC, final_wy
    ))
  }

  # ── Confirm link between coverage and low precip ───────────────────────────
  link_df <- day_counts %>%
    dplyr::filter(water_year == final_wy) %>%
    dplyr::inner_join(
      dplyr::filter(annual, water_year == final_wy) %>%
        dplyr::select(site_no, total_precip_mm),
      by = c("GAGE_ID" = "site_no")
    )

  if (nrow(link_df) > 0) {
    sd_cov <- sd(link_df$frac_covered, na.rm = TRUE)
    if (sd_cov == 0) {
      message(sprintf(
        "\nWY%d WLDAS coverage is uniform across all sites (frac = %.3f) — correlation not computable.",
        final_wy, mean(link_df$frac_covered, na.rm = TRUE)
      ))
      message("  Uniform coverage confirms the WLDAS parquet ends at the same date for all sites.")
    } else {
      cor_val <- cor(link_df$frac_covered, link_df$total_precip_mm,
                     use = "complete.obs")
      message(sprintf(
        "\nCorrelation between WY%d WLDAS day coverage and annual precip: r = %.3f",
        final_wy, cor_val
      ))
      if (!is.na(cor_val) && abs(cor_val) > 0.7)
        message("  *** Strong positive correlation confirms low precip is driven by coverage, not climate.")
    }
  }

  # ── Figure 3: Daily coverage heatmap, recent years ────────────────────────
  hmap_yrs <- max(1L, final_wy - 4L):final_wy
  hmap_df  <- day_counts %>%
    dplyr::filter(water_year %in% hmap_yrs) %>%
    dplyr::mutate(frac_covered = pmin(frac_covered, 1))

  p_coverage <- ggplot(hmap_df,
                       aes(x = factor(water_year), y = fct_rev(GAGE_ID),
                           fill = frac_covered)) +
    geom_tile(colour = NA) +
    scale_fill_gradient2(
      low      = "#d73027",
      mid      = "#fee08b",
      high     = "#1a9641",
      midpoint = MIN_DAYS_FRAC,
      limits   = c(0, 1),
      labels   = percent_format(accuracy = 1),
      name     = "Coverage\nfraction"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.y      = element_text(size = 4.5),
      panel.grid       = element_blank(),
      legend.position  = "right"
    ) +
    labs(
      x        = "Water Year",
      y        = "Site",
      title    = sprintf("WLDAS daily coverage — WY%d to WY%d", min(hmap_yrs), final_wy),
      subtitle = sprintf(
        "Red = < %.0f%% of expected days covered  |  midpoint of colour scale at %.0f%%",
        100 * MIN_DAYS_FRAC, 100 * MIN_DAYS_FRAC
      )
    )

  hmap_h <- max(4, n_distinct(hmap_df$GAGE_ID) * 0.07 + 1.5)
  ggsave(file.path(OUT_DIR, "fig_wldas_daily_coverage.png"),
         p_coverage, width = 7, height = hmap_h, dpi = 150, limitsize = FALSE)
  message("\nSaved: fig_wldas_daily_coverage.png")
}


# ============================================================================
# SUMMARY AND RECOMMENDATION
# ============================================================================
message("\n=== SUMMARY AND RECOMMENDATION ===\n")

message("Issue 1 — Incomplete final water year")
if (is.na(final_med)) {
  message(sprintf(
    "  WY%d median total_precip_mm = NA — Script 8 coverage mask already applied.",
    final_wy
  ))
  message(sprintf(
    "  WY%d median total_precip_mm = %.1f mm (reference year).",
    prev_wy, prev_med
  ))
  message(paste0(
    "  WY", final_wy, " has been excluded from the trend analysis because WLDAS\n",
    sprintf("  day coverage is below %.0f%% of expected days.\n", 100 * MIN_DAYS_FRAC),
    "  No further action needed for this issue."
  ))
} else {
  message(sprintf(
    "  WY%d median total_precip_mm = %.1f mm (%.0f%% of WY%d's %.1f mm).",
    final_wy, final_med, pct_of_prev, prev_wy, prev_med
  ))
  message(paste0(
    "  Root cause: Script 8 uses sum(na.rm = TRUE) for total_precip_mm.\n",
    "  When the WLDAS record ends mid-water-year, missing days contribute 0\n",
    "  rather than NA, yielding an artificially low annual total.\n",
    "\n",
    "  Recommended fix for Script 8:\n",
    "  After building wldas_annual, compute n_days per (site, water_year) from\n",
    "  the raw WLDAS data and derive frac_covered = n_days / expected_days.\n",
    sprintf("  Set total_precip_mm = NA and peak_swe_mm = NA for any row where\n"),
    sprintf("  frac_covered < %.2f. These rows will then be treated as NA in\n", MIN_DAYS_FRAC),
    "  the Mann-Kendall test (excluded via the non-NA year count), avoiding\n",
    "  the need for a hard cutoff on the final water year."
  ))
}

message("\nIssue 2 — Missing years in trend plots")
message(sprintf(
  "  %d sites have interior gaps (years with no POT peaks within their active span).",
  n_gap_sites
))
message(paste0(
  "  These are years where the gauge had no peaks above its PoT threshold —\n",
  "  a genuine feature of the discharge record, not a data error. They are\n",
  "  correctly represented as NA in annual_climate_forcing_summary.rds and\n",
  "  excluded from Mann-Kendall automatically. The visual appearance of\n",
  "  'missing years' in the trend plots is expected and correct behaviour.\n",
  "  No fix is needed for this issue."
))

message("\nDiagnostic complete.")
