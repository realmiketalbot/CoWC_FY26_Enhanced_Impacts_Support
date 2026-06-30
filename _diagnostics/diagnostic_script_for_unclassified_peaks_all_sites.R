#!/usr/bin/env Rscript

# =============================================================================
# diagnostic_unclassified_by_site_metadata.R
#
# Fleet-level diagnostic: explores whether unclassified peak rates vary
# systematically with site characteristics (elevation, drainage area, dam
# presence, reference class, hydrologic disturbance, HUC02 / geographic
# region).
#
# Two failure modes are tracked separately:
#   (1) WLDAS-coverage failures   — peak predates or postdates WLDAS archive
#   (2) Threshold failures         — complete WLDAS window but no type matched
#
# Unclassified rate denominator: peaks with a complete WLDAS window
# (i.e., n_days_wldas == 8; coverage failures excluded from rate calculation
#  but shown separately).
#
# INPUTS:
#   _data/pot_peaks_classified_adaptive.rds  — output of script 6b
#   _data/site_metadata.rds         — output of script 5
#
# OUTPUTS (all written to OUT_DIR):
#   fig01_coverage_vs_threshold_failures.png
#   fig02_unclassified_rate_by_elevation_bin.png
#   fig03_unclassified_rate_by_area_bin.png
#   fig04_unclassified_rate_by_dam_presence.png
#   fig05_unclassified_rate_by_ref_class.png
#   fig06_unclassified_rate_by_huc02.png
#   fig07_rain_melt_scatter_all_sites.png
#   fig08_site_rate_distribution.png
# =============================================================================

library(tidyverse)
library(ggplot2)
library(scales)
library(patchwork)   # for combining panels; install if missing

OUT_DIR        <- "_diagnostics/unclassified_metadata"
PEAKS_FILE     <- "_data/pot_peaks_classified_adaptive.rds"
METADATA_FILE  <- "_data/site_metadata.rds"

# Davenport et al. 2020 thresholds (must match script 5)
thresh_rain_high <- 10
thresh_rain_low  <-  5
thresh_melt_high <- 10
thresh_melt_low  <-  5

SAVE_WIDTH  <- 9    # inches
SAVE_HEIGHT <- 6
SAVE_DPI    <- 150

# =============================================================================
# SETUP
# =============================================================================
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(p, filename, width = SAVE_WIDTH, height = SAVE_HEIGHT) {
  path <- file.path(OUT_DIR, filename)
  ggsave(path, p, width = width, height = height, dpi = SAVE_DPI)
  message(sprintf("  Saved: %s", path))
}

# Shared theme
theme_diag <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "grey92", colour = NA),
      legend.position   = "bottom"
    )
}

# =============================================================================
# LOAD DATA
# =============================================================================
message("Loading classified peaks and site metadata...")
peaks <- readRDS(PEAKS_FILE)
meta  <- readRDS(METADATA_FILE)

# Reconstruct a pure WLDAS-coverage flag (n_days_data >= 90% of search_radius).
# data_sufficient in the adaptive RDS bundles the MIN_TOTAL_FLUX_MM floor in with
# the coverage check, which would collapse the coverage-vs-threshold distinction
# fig01 is designed to show. complete_window is rebuilt here without the flux
# floor so the original failure-mode separation is preserved.
peaks <- peaks %>%
  dplyr::mutate(complete_window = n_days_data >= ceiling(short_window * 0.90))

message(sprintf("  Peaks:  %d rows, %d sites", nrow(peaks), n_distinct(peaks$site_no)))
message(sprintf("  Meta:   %d sites", nrow(meta)))

# =============================================================================
# CLASSIFY EACH PEAK'S FAILURE MODE
# =============================================================================
# complete_window == FALSE can mean two things:
#   (a) outside WLDAS temporal coverage  → n_days_wldas == 0
#   (b) partial window (rare edge cases) → 0 < n_days_wldas < 8
# We treat (a) and (b) together as "coverage failure" and separate them
# from "threshold failure" (complete_window == TRUE, flood_type == "unclassified").

peaks <- peaks %>%
  mutate(
    failure_mode = case_when(
      flood_type != "unclassified"          ~ "classified",
      !complete_window                      ~ "coverage_failure",
      TRUE                                  ~ "threshold_failure"
    )
  )

message("\nFleet-wide failure mode summary:")
peaks %>% count(failure_mode) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

# =============================================================================
# SITE-LEVEL SUMMARY
# =============================================================================
# Unclassified RATE = threshold failures / (classified + threshold failures)
# Coverage failures are tallied separately and not included in the denominator.

site_summary <- peaks %>%
  group_by(site_no) %>%
  summarise(
    n_total             = n(),
    n_classified        = sum(failure_mode == "classified"),
    n_threshold_fail    = sum(failure_mode == "threshold_failure"),
    n_coverage_fail     = sum(failure_mode == "coverage_failure"),
    n_complete_window   = n_classified + n_threshold_fail,
    unclassified_rate   = if_else(
      n_complete_window > 0,
      n_threshold_fail / n_complete_window,
      NA_real_
    ),
    coverage_fail_rate  = n_coverage_fail / n_total,
    .groups = "drop"
  ) %>%
  left_join(meta, by = "site_no")

message(sprintf("\nSite summary: %d sites, %d with at least one complete-window peak",
                nrow(site_summary), sum(!is.na(site_summary$unclassified_rate))))

# Bin continuous predictors
site_summary <- site_summary %>%
  mutate(
    elev_bin  = cut(elev_mean_m,
                    breaks = c(-Inf, 1500, 2500, 3500, Inf),
                    labels = c("< 1500 m", "1500–2500 m", "2500–3500 m", "> 3500 m")),
    area_bin  = cut(drain_sqkm,
                    breaks = c(0, 100, 500, 2000, Inf),
                    labels = c("< 100 km²", "100–500 km²", "500–2000 km²", "> 2000 km²")),
    dam_group = if_else(ndams_2009 > 0, "≥ 1 dam upstream", "No dams upstream"),
    ref_group = case_when(
      class == "Ref"     ~ "Reference",
      class == "Non-ref" ~ "Non-Reference",
      TRUE               ~ "Unknown"
    )
  )

# =============================================================================
# FIG 01: FLEET-WIDE FAILURE MODE STACKED BAR
# Summary of classified / threshold-fail / coverage-fail across the whole fleet,
# sorted by unclassified rate so patterns by site are visible.
# =============================================================================
message("\nBuilding Fig 01 — coverage vs threshold failures...")

fig01_data <- site_summary %>%
  filter(!is.na(unclassified_rate)) %>%
  arrange(desc(unclassified_rate)) %>%
  mutate(site_rank = row_number()) %>%
  select(site_rank, site_no, n_classified, n_threshold_fail, n_coverage_fail) %>%
  pivot_longer(c(n_classified, n_threshold_fail, n_coverage_fail),
               names_to = "mode", values_to = "count") %>%
  mutate(mode = recode(mode,
    n_classified     = "Classified",
    n_threshold_fail = "Threshold failure",
    n_coverage_fail  = "Coverage failure"
  ))

mode_colors <- c(
  "Classified"        = "#4dac26",
  "Threshold failure" = "#d7191c",
  "Coverage failure"  = "#fdae61"
)

p01 <- ggplot(fig01_data, aes(x = site_rank, y = count, fill = mode)) +
  geom_col(position = "stack", width = 1) +
  scale_fill_manual(values = mode_colors, name = NULL) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title   = "Failure mode breakdown by site (sorted by unclassified rate)",
    x       = "Site rank (highest → lowest unclassified rate)",
    y       = "Peak count",
    caption = "Coverage failure: peak outside WLDAS temporal coverage. Threshold failure: complete window but no type matched."
  ) +
  theme_diag()

save_fig(p01, "fig01_coverage_vs_threshold_failures.png")

# =============================================================================
# FIG 02: UNCLASSIFIED RATE BY ELEVATION BIN
# =============================================================================
message("Building Fig 02 — unclassified rate by elevation...")

elev_summary <- site_summary %>%
  filter(!is.na(unclassified_rate), !is.na(elev_bin)) %>%
  group_by(elev_bin) %>%
  summarise(
    n_sites       = n(),
    median_rate   = median(unclassified_rate),
    q25           = quantile(unclassified_rate, 0.25),
    q75           = quantile(unclassified_rate, 0.75),
    .groups = "drop"
  )

p02a <- ggplot(site_summary %>% filter(!is.na(unclassified_rate), !is.na(elev_bin)),
               aes(x = elev_bin, y = unclassified_rate)) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, colour = "steelblue") +
  geom_pointrange(data = elev_summary,
                  aes(y = median_rate, ymin = q25, ymax = q75),
                  colour = "black", size = 0.6, linewidth = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title    = "Threshold-failure unclassified rate by mean basin elevation",
    subtitle = "Black point = median; bar = IQR. Points = individual sites.",
    x        = "Mean basin elevation",
    y        = "Unclassified rate (complete windows)"
  ) +
  theme_diag()

p02b <- ggplot(site_summary %>% filter(!is.na(unclassified_rate), !is.na(elev_mean_m)),
               aes(x = elev_mean_m, y = unclassified_rate)) +
  geom_point(alpha = 0.5, size = 1.5, colour = "steelblue") +
  geom_smooth(method = "loess", span = 0.6, se = TRUE, colour = "black",
              linewidth = 0.8) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    x = "Mean basin elevation (m)",
    y = "Unclassified rate"
  ) +
  theme_diag()

p02 <- p02a / p02b +
  plot_annotation(tag_levels = "A")

save_fig(p02, "fig02_unclassified_rate_by_elevation_bin.png",
         height = 10)

# =============================================================================
# FIG 03: UNCLASSIFIED RATE BY DRAINAGE AREA BIN
# =============================================================================
message("Building Fig 03 — unclassified rate by drainage area...")

area_summary <- site_summary %>%
  filter(!is.na(unclassified_rate), !is.na(area_bin)) %>%
  group_by(area_bin) %>%
  summarise(
    n_sites       = n(),
    median_rate   = median(unclassified_rate),
    q25           = quantile(unclassified_rate, 0.25),
    q75           = quantile(unclassified_rate, 0.75),
    .groups = "drop"
  )

p03a <- ggplot(site_summary %>% filter(!is.na(unclassified_rate), !is.na(area_bin)),
               aes(x = area_bin, y = unclassified_rate)) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, colour = "#e66101") +
  geom_pointrange(data = area_summary,
                  aes(y = median_rate, ymin = q25, ymax = q75),
                  colour = "black", size = 0.6, linewidth = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title    = "Threshold-failure unclassified rate by drainage area",
    subtitle = "Black point = median; bar = IQR.",
    x        = "Drainage area",
    y        = "Unclassified rate (complete windows)"
  ) +
  theme_diag()

p03b <- ggplot(site_summary %>% filter(!is.na(unclassified_rate), !is.na(drain_sqkm)),
               aes(x = drain_sqkm, y = unclassified_rate)) +
  geom_point(alpha = 0.5, size = 1.5, colour = "#e66101") +
  geom_smooth(method = "loess", span = 0.6, se = TRUE, colour = "black",
              linewidth = 0.8) +
  scale_x_log10(labels = comma) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    x = "Drainage area (km², log scale)",
    y = "Unclassified rate"
  ) +
  theme_diag()

p03 <- p03a / p03b +
  plot_annotation(tag_levels = "A")

save_fig(p03, "fig03_unclassified_rate_by_area_bin.png",
         height = 10)

# =============================================================================
# FIG 04: UNCLASSIFIED RATE BY DAM PRESENCE
# =============================================================================
message("Building Fig 04 — unclassified rate by dam presence...")

dam_summary <- site_summary %>%
  filter(!is.na(unclassified_rate), !is.na(dam_group)) %>%
  group_by(dam_group) %>%
  summarise(
    n_sites     = n(),
    median_rate = median(unclassified_rate),
    q25         = quantile(unclassified_rate, 0.25),
    q75         = quantile(unclassified_rate, 0.75),
    .groups = "drop"
  )

p04 <- ggplot(site_summary %>% filter(!is.na(unclassified_rate), !is.na(dam_group)),
              aes(x = dam_group, y = unclassified_rate, colour = dam_group)) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2, show.legend = FALSE) +
  geom_pointrange(data = dam_summary,
                  aes(y = median_rate, ymin = q25, ymax = q75),
                  colour = "black", size = 0.7, linewidth = 1.2) +
  geom_text(data = dam_summary,
            aes(y = -0.04, label = sprintf("n = %d", n_sites)),
            colour = "grey30", size = 3.5) +
  scale_colour_manual(values = c("No dams upstream" = "steelblue",
                                 "≥ 1 dam upstream" = "#d7191c")) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(-0.06, 1)) +
  labs(
    title    = "Threshold-failure unclassified rate by dam presence",
    subtitle = "Black point = median; bar = IQR.",
    x        = NULL,
    y        = "Unclassified rate (complete windows)"
  ) +
  theme_diag()

save_fig(p04, "fig04_unclassified_rate_by_dam_presence.png",
         width = 6, height = 6)

# =============================================================================
# FIG 05: UNCLASSIFIED RATE BY REFERENCE CLASS AND DISTURBANCE INDEX
# =============================================================================
message("Building Fig 05 — unclassified rate by reference class...")

ref_summary <- site_summary %>%
  filter(!is.na(unclassified_rate), ref_group != "Unknown") %>%
  group_by(ref_group) %>%
  summarise(
    n_sites     = n(),
    median_rate = median(unclassified_rate),
    q25         = quantile(unclassified_rate, 0.25),
    q75         = quantile(unclassified_rate, 0.75),
    .groups = "drop"
  )

p05a <- ggplot(site_summary %>% filter(!is.na(unclassified_rate), ref_group != "Unknown"),
               aes(x = ref_group, y = unclassified_rate, colour = ref_group)) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2, show.legend = FALSE) +
  geom_pointrange(data = ref_summary,
                  aes(y = median_rate, ymin = q25, ymax = q75),
                  colour = "black", size = 0.7, linewidth = 1.2) +
  geom_text(data = ref_summary,
            aes(y = -0.04, label = sprintf("n = %d", n_sites)),
            colour = "grey30", size = 3.5) +
  scale_colour_manual(values = c("Reference" = "#2c7bb6",
                                 "Non-Reference" = "#d7191c")) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(-0.06, 1)) +
  labs(
    title    = "Unclassified rate by GAGES-II reference class",
    x        = NULL,
    y        = "Unclassified rate (complete windows)"
  ) +
  theme_diag()

p05b <- ggplot(
  site_summary %>% filter(!is.na(unclassified_rate), !is.na(hydro_disturb_idx)),
  aes(x = hydro_disturb_idx, y = unclassified_rate)
) +
  geom_point(alpha = 0.5, size = 1.5, colour = "grey40") +
  geom_smooth(method = "loess", span = 0.8, se = TRUE, colour = "#d7191c",
              linewidth = 0.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = "Unclassified rate vs. hydrologic disturbance index",
    x     = "Hydrologic disturbance index (GAGES-II)",
    y     = "Unclassified rate"
  ) +
  theme_diag()

p05 <- p05a / p05b +
  plot_annotation(tag_levels = "A")

save_fig(p05, "fig05_unclassified_rate_by_ref_class.png", height = 10)

# =============================================================================
# FIG 06: UNCLASSIFIED RATE BY HUC02
# =============================================================================
message("Building Fig 06 — unclassified rate by HUC02...")

huc_summary <- site_summary %>%
  filter(!is.na(unclassified_rate), !is.na(huc02)) %>%
  group_by(huc02) %>%
  summarise(
    n_sites     = n(),
    median_rate = median(unclassified_rate),
    q25         = quantile(unclassified_rate, 0.25),
    q75         = quantile(unclassified_rate, 0.75),
    .groups = "drop"
  ) %>%
  arrange(median_rate) %>%
  mutate(huc02 = factor(huc02, levels = huc02))

p06 <- ggplot(huc_summary, aes(x = huc02, y = median_rate)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.3, colour = "grey30") +
  geom_text(aes(y = -0.02, label = sprintf("n=%d", n_sites)),
            colour = "grey30", size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(-0.04, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Threshold-failure unclassified rate by HUC02 region",
    subtitle = "Bar = median; whiskers = IQR. n = number of sites.",
    x        = "HUC02 code",
    y        = "Unclassified rate (complete windows)"
  ) +
  theme_diag() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p06, "fig06_unclassified_rate_by_huc02.png")

# =============================================================================
# FIG 07: RAIN VS MELT SCATTER — ALL THRESHOLD-FAILURE PEAKS
# Shows where unclassified peaks sit in classification space relative to
# classified ones; faceted by elevation bin. Adaptive flux variables are used
# (rain_flux_short, melt_flux_short — search_radius window with NOISE_FLOOR_MM
# applied), so the axes are not directly comparable to davenport's 8-day totals.
# =============================================================================
MIN_TOTAL_FLUX_MM_FIG07 <- 1.0  # mirrors MIN_TOTAL_FLUX_MM in 6b

message("Building Fig 07 — rain vs melt scatter (all threshold-failure peaks)...")

scatter_data <- peaks %>%
  filter(complete_window, !is.na(rain_flux_short), !is.na(melt_flux_short)) %>%
  left_join(site_summary %>% select(site_no, elev_bin), by = "site_no") %>%
  filter(!is.na(elev_bin)) %>%
  mutate(
    type_label = case_when(
      flood_type == "rainfall" ~ "Rainfall",
      flood_type == "snowmelt" ~ "Snowmelt",
      flood_type == "ROS"      ~ "Rain-on-Snow",
      TRUE                     ~ "Unclassified"
    ),
    # Cap axes for legibility; the vast majority of peaks are in low range
    rain_plot = pmin(rain_flux_short, 60),
    melt_plot = pmin(melt_flux_short, 60)
  )

type_colors <- c(
  "Rainfall"     = "#e41a1c",
  "Snowmelt"     = "#377eb8",
  "Rain-on-Snow" = "#4daf4a",
  "Unclassified" = "#aaaaaa"
)
type_alpha <- c(
  "Rainfall" = 0.25, "Snowmelt" = 0.25,
  "Rain-on-Snow" = 0.25, "Unclassified" = 0.7
)
type_size <- c(
  "Rainfall" = 0.8, "Snowmelt" = 0.8,
  "Rain-on-Snow" = 0.8, "Unclassified" = 1.2
)

# Unclassified plotted last (on top)
scatter_data <- scatter_data %>%
  arrange(type_label == "Unclassified")

p07 <- ggplot(scatter_data,
              aes(x = rain_plot, y = melt_plot, colour = type_label)) +
  # MIN_TOTAL_FLUX_MM floor: rain_flux_short + melt_flux_short = 1.0 mm.
  # Peaks below this diagonal are unclassified by the flux floor.
  geom_abline(slope = -1, intercept = MIN_TOTAL_FLUX_MM_FIG07,
              linetype = "dashed", colour = "black", linewidth = 0.4) +
  geom_point(aes(alpha = type_label, size = type_label)) +
  scale_colour_manual(values = type_colors, name = "Flood type") +
  scale_alpha_manual(values = type_alpha,   guide = "none") +
  scale_size_manual( values = type_size,    guide = "none") +
  facet_wrap(~ elev_bin, ncol = 2) +
  labs(
    title    = "Adaptive classification space: rain vs. melt flux over search_radius window",
    subtitle = "Unclassified peaks highlighted in grey. Axes capped at 60 mm.",
    x        = "rain_flux_short (mm, capped at 60)",
    y        = "melt_flux_short (mm, capped at 60)",
    caption  = sprintf(
      "Dashed diagonal: MIN_TOTAL_FLUX_MM floor (rain + melt = %.1f mm). Peaks below the line fail the flux-sufficiency check.",
      MIN_TOTAL_FLUX_MM_FIG07
    )
  ) +
  theme_diag() +
  theme(legend.position = "bottom")

save_fig(p07, "fig07_rain_melt_scatter_all_sites.png",
         width = 10, height = 8)

# =============================================================================
# FIG 08: DISTRIBUTION OF SITE-LEVEL UNCLASSIFIED RATES
# Histogram + density, split by reference class
# =============================================================================
message("Building Fig 08 — distribution of site-level unclassified rates...")

p08 <- ggplot(
  site_summary %>% filter(!is.na(unclassified_rate), ref_group != "Unknown"),
  aes(x = unclassified_rate, fill = ref_group)
) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 position = "identity", alpha = 0.55) +
  geom_density(aes(colour = ref_group), linewidth = 0.8, fill = NA) +
  geom_vline(
    data = site_summary %>%
      filter(!is.na(unclassified_rate), ref_group != "Unknown") %>%
      group_by(ref_group) %>%
      summarise(med = median(unclassified_rate), .groups = "drop"),
    aes(xintercept = med, colour = ref_group),
    linetype = "dashed", linewidth = 0.8
  ) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values  = c("Reference" = "#2c7bb6", "Non-Reference" = "#d7191c"),
                    name = NULL) +
  scale_colour_manual(values = c("Reference" = "#2c7bb6", "Non-Reference" = "#d7191c"),
                      name = NULL) +
  labs(
    title    = "Distribution of site-level threshold-failure rates",
    subtitle = "Dashed lines = group medians.",
    x        = "Unclassified rate (complete windows)",
    y        = "Density"
  ) +
  theme_diag()

save_fig(p08, "fig08_site_rate_distribution.png", width = 7, height = 5)

# =============================================================================
# THRESHOLD JUSTIFICATION DIAGNOSTICS
# =============================================================================
# Figures 09–11 evaluate whether MIN_TOTAL_FLUX_MM (the data-sufficiency floor
# used by 6b) is set sensibly. 09 and 10 examine the empirical distribution of
# total_flux_short and the stability of rain_frac across the threshold; 11 runs
# a full pipeline sweep, regenerating 6b/7 outputs at a range of thresholds and
# reporting how downstream metrics change.
# =============================================================================

ADAPTIVE_PEAKS_FILE <- "_data/pot_peaks_classified_adaptive.rds"
MIN_TOTAL_FLUX_MM_OPERATIONAL <- 1.0   # current pipeline value

message("\nLoading adaptive classified peaks for threshold diagnostics...")
peaks_adap <- readRDS(ADAPTIVE_PEAKS_FILE)
message(sprintf("  %d adaptive peaks across %d sites",
                nrow(peaks_adap), n_distinct(peaks_adap$site_no)))

# =============================================================================
# FIG 09: TOTAL FLUX DISTRIBUTION
# Empirical distribution of total_flux_short across all peaks. Two panels: a
# log10 ECDF spanning the full range (where the threshold sits in the tail
# matters), and a linear-scale histogram zoomed to 0–5 mm (where the threshold
# discriminates inactive from active peaks).
# =============================================================================
message("Building Fig 09 — total_flux_short distribution...")

flux_data <- peaks_adap %>%
  filter(!is.na(total_flux_short), is.finite(total_flux_short))

p09a <- ggplot(flux_data, aes(x = pmax(total_flux_short, 1e-3))) +
  stat_ecdf(geom = "step", colour = "steelblue", linewidth = 0.8) +
  geom_vline(xintercept = MIN_TOTAL_FLUX_MM_OPERATIONAL,
             linetype = "dashed", colour = "#d7191c", linewidth = 0.7) +
  annotate("text", x = MIN_TOTAL_FLUX_MM_OPERATIONAL, y = 0.05,
           label = sprintf(" MIN_TOTAL_FLUX_MM = %.1f", MIN_TOTAL_FLUX_MM_OPERATIONAL),
           hjust = 0, colour = "#d7191c", size = 3.5) +
  scale_x_log10(labels = comma,
                breaks = c(0.01, 0.1, 1, 10, 100, 1000)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Total flux ECDF (log10) — full range",
    x     = "total_flux_short (mm, log10; 0 floored to 1e-3)",
    y     = "Cumulative fraction of peaks"
  ) +
  theme_diag()

p09b <- ggplot(flux_data %>% filter(total_flux_short <= 5),
               aes(x = total_flux_short)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", alpha = 0.8,
                 boundary = 0) +
  geom_vline(xintercept = MIN_TOTAL_FLUX_MM_OPERATIONAL,
             linetype = "dashed", colour = "#d7191c", linewidth = 0.7) +
  scale_x_continuous(breaks = seq(0, 5, 0.5)) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Near-zero detail (0–5 mm, linear)",
    x     = "total_flux_short (mm)",
    y     = "Peak count"
  ) +
  theme_diag()

p09 <- p09a / p09b + plot_annotation(tag_levels = "A")

save_fig(p09, "fig09_total_flux_distribution.png", height = 9)

# =============================================================================
# FIG 10: RAIN FRACTION STABILITY VS. TOTAL FLUX
# Below the threshold, rain_frac is dominated by 0/1 splits driven by single
# noisy days. Above the threshold, the pool of contributing days grows and
# rain_frac should populate the interior (0,1) interval. Panel A: hexbin scatter
# 0–10 mm. Panel B: rolling IQR of rain_frac as a function of total_flux_short.
# =============================================================================
message("Building Fig 10 — rain_frac stability vs total_flux_short...")

stability_data <- peaks_adap %>%
  filter(!is.na(rain_frac), !is.na(total_flux_short),
         total_flux_short > 0, total_flux_short <= 10)

p10a <- ggplot(stability_data, aes(x = total_flux_short, y = rain_frac)) +
  geom_hex(bins = 50) +
  geom_vline(xintercept = MIN_TOTAL_FLUX_MM_OPERATIONAL,
             linetype = "dashed", colour = "#d7191c", linewidth = 0.7) +
  scale_fill_viridis_c(option = "magma", direction = -1, trans = "log10",
                       labels = comma, name = "Peaks\n(log10)") +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  labs(
    title = "rain_frac vs. total_flux_short (0–10 mm)",
    x     = "total_flux_short (mm)",
    y     = "rain_frac"
  ) +
  theme_diag()

# Rolling IQR in 0.2-mm flux bins (centred)
flux_bins <- seq(0, 10, 0.2)
bin_summary <- stability_data %>%
  mutate(flux_bin = cut(total_flux_short, breaks = flux_bins,
                        include.lowest = TRUE)) %>%
  group_by(flux_bin) %>%
  summarise(
    n            = n(),
    flux_mid     = mean(total_flux_short),
    rain_frac_q25 = quantile(rain_frac, 0.25),
    rain_frac_q75 = quantile(rain_frac, 0.75),
    rain_frac_iqr = rain_frac_q75 - rain_frac_q25,
    rain_frac_sd  = sd(rain_frac),
    .groups = "drop"
  ) %>%
  filter(n >= 20)

p10b <- ggplot(bin_summary, aes(x = flux_mid)) +
  geom_line(aes(y = rain_frac_iqr, colour = "IQR"), linewidth = 0.9) +
  geom_line(aes(y = rain_frac_sd,  colour = "SD"),  linewidth = 0.9) +
  geom_vline(xintercept = MIN_TOTAL_FLUX_MM_OPERATIONAL,
             linetype = "dashed", colour = "#d7191c", linewidth = 0.7) +
  scale_colour_manual(values = c("IQR" = "#1b7837", "SD" = "#762a83"),
                      name = NULL) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    title    = "Within-bin spread of rain_frac (0.2-mm bins, n ≥ 20)",
    subtitle = "Lower spread above threshold ⇒ rain_frac is more stable as data sufficiency improves.",
    x        = "total_flux_short (mm)",
    y        = "Spread of rain_frac"
  ) +
  theme_diag()

p10 <- p10a / p10b + plot_annotation(tag_levels = "A")

save_fig(p10, "fig10_rain_frac_stability_vs_flux.png", height = 10)

# =============================================================================
# FIG 11: THRESHOLD SENSITIVITY SWEEP
# Re-runs 6b and 7 across a sweep of MIN_TOTAL_FLUX_MM values and plots how
# the unclassified rate, per-type counts, and 100-year mixture quantile change.
# Operational RDS files are backed up and restored regardless of sweep outcome.
# =============================================================================
message("\nBuilding Fig 11 — threshold sensitivity sweep...")
message("  WARNING: this section re-runs scripts 6b and 7 multiple times.")
message("  Expect a substantial wall-clock cost (parallel FFA × n_thresholds).")

THRESH_GRID    <- c(0.1, 0.5, 1.0, 2.0, 5.0, 10.0)
SWEEP_DIR      <- file.path(OUT_DIR, "_threshold_sweep")
SWEEP_BACKUP   <- file.path(OUT_DIR, "_threshold_sweep_backup")
SCRIPT_6B      <- "_scripts/6b_classify_usgs_unit_data_peaks_adaptive.R"
SCRIPT_7       <- "_scripts/7_compile_ffa_results.R"

dir.create(SWEEP_DIR,    recursive = TRUE, showWarnings = FALSE)
dir.create(SWEEP_BACKUP, recursive = TRUE, showWarnings = FALSE)

OPERATIONAL_FILES <- c(
  "pot_peaks_classified_adaptive.rds",
  "ffa_curves_gpd.rds",    "ffa_curves_lp3.rds",
  "ffa_quantiles_gpd.rds", "ffa_quantiles_lp3.rds",
  "ffa_params_gpd.rds",    "ffa_params_lp3.rds",
  "ffa_classification_summary.rds",
  "ffa_ros_temporal.rds"
)

backup_operational <- function() {
  for (f in OPERATIONAL_FILES) {
    src <- file.path("_data", f)
    if (file.exists(src))
      file.copy(src, file.path(SWEEP_BACKUP, f), overwrite = TRUE)
  }
}

restore_operational <- function() {
  for (f in OPERATIONAL_FILES) {
    bkp <- file.path(SWEEP_BACKUP, f)
    if (file.exists(bkp))
      file.copy(bkp, file.path("_data", f), overwrite = TRUE)
  }
  message("  Operational RDS files restored from backup.")
}

stash_iteration <- function(thresh) {
  tag <- sprintf("thresh_%g", thresh)
  for (f in OPERATIONAL_FILES) {
    src <- file.path("_data", f)
    if (file.exists(src)) {
      dst <- file.path(SWEEP_DIR, sub("\\.rds$", paste0("_", tag, ".rds"), f))
      file.copy(src, dst, overwrite = TRUE)
    }
  }
}

run_iteration <- function(thresh) {
  env <- new.env(parent = globalenv())
  assign("MIN_TOTAL_FLUX_MM", thresh, envir = env)
  sys.source(SCRIPT_6B, envir = env)
  sys.source(SCRIPT_7,  envir = env)
  invisible(NULL)
}

collect_metrics <- function(thresh) {
  tag         <- sprintf("thresh_%g", thresh)
  peaks_path  <- file.path(SWEEP_DIR, paste0("pot_peaks_classified_adaptive_", tag, ".rds"))
  qgpd_path   <- file.path(SWEEP_DIR, paste0("ffa_quantiles_gpd_", tag, ".rds"))

  pk  <- readRDS(peaks_path)
  qg  <- readRDS(qgpd_path)

  site_rates <- pk %>%
    group_by(site_no) %>%
    summarise(
      n_total          = n(),
      n_unclass_raw    = sum(as.character(flood_type_raw) == "unclassified"),
      unclass_rate     = n_unclass_raw / n_total,
      .groups = "drop"
    )

  type_counts <- pk %>%
    mutate(ftr = as.character(flood_type_raw)) %>%
    filter(ftr %in% c("rainfall", "snowmelt", "ROS")) %>%
    count(ftr, name = "n") %>%
    rename(flood_type_raw = ftr)

  q100 <- qg %>%
    filter(method == "adaptive",
           curve_type == "mixture",
           return_period == 100) %>%
    select(site_no, q100 = discharge_cfs)

  list(
    threshold     = thresh,
    site_rates    = site_rates %>% mutate(threshold = thresh),
    type_counts   = type_counts %>% mutate(threshold = thresh),
    q100          = q100 %>% mutate(threshold = thresh),
    median_uncl   = median(site_rates$unclass_rate, na.rm = TRUE),
    iqr_uncl_q25  = quantile(site_rates$unclass_rate, 0.25, na.rm = TRUE),
    iqr_uncl_q75  = quantile(site_rates$unclass_rate, 0.75, na.rm = TRUE),
    n_rainfall    = sum(type_counts$n[type_counts$flood_type_raw == "rainfall"]),
    n_snowmelt    = sum(type_counts$n[type_counts$flood_type_raw == "snowmelt"]),
    n_ros         = sum(type_counts$n[type_counts$flood_type_raw == "ROS"])
  )
}

options(future.globals.maxSize = 2 * 1024^3)  # 2 GB; required for sweep loop's nested future_pmap

message(sprintf("  Backing up operational RDS files to %s/", SWEEP_BACKUP))
backup_operational()

sweep_results <- tryCatch({
  acc <- list()
  for (i in seq_along(THRESH_GRID)) {
    thresh <- THRESH_GRID[i]
    message(sprintf("\n── Iteration %d/%d: MIN_TOTAL_FLUX_MM = %g ──",
                    i, length(THRESH_GRID), thresh))

    run_iteration(thresh)
    stash_iteration(thresh)
    m <- collect_metrics(thresh)
    acc[[length(acc) + 1L]] <- m

    message(sprintf(
      "  thresh=%-5g | median unclass rate=%.3f (IQR %.3f–%.3f) | rainfall=%d snowmelt=%d ROS=%d",
      thresh, m$median_uncl, m$iqr_uncl_q25, m$iqr_uncl_q75,
      m$n_rainfall, m$n_snowmelt, m$n_ros
    ))
  }
  acc
}, finally = restore_operational())

# Assemble long-format frames for plotting
sweep_site_rates  <- bind_rows(lapply(sweep_results, `[[`, "site_rates"))
sweep_type_counts <- bind_rows(lapply(sweep_results, `[[`, "type_counts"))
sweep_q100        <- bind_rows(lapply(sweep_results, `[[`, "q100"))

# Per-site % change in 100-yr mixture quantile relative to thresh = 1.0
baseline_q100 <- sweep_q100 %>%
  filter(threshold == 1.0) %>%
  select(site_no, q100_baseline = q100)

q100_pct_change <- sweep_q100 %>%
  inner_join(baseline_q100, by = "site_no") %>%
  filter(!is.na(q100), !is.na(q100_baseline), q100_baseline > 0) %>%
  mutate(pct_change = 100 * (q100 - q100_baseline) / q100_baseline)

# ── Panel A: median unclassified rate ± IQR vs. threshold ────────────────────
panel_a_summary <- sweep_site_rates %>%
  group_by(threshold) %>%
  summarise(
    median_rate = median(unclass_rate, na.rm = TRUE),
    q25         = quantile(unclass_rate, 0.25, na.rm = TRUE),
    q75         = quantile(unclass_rate, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

pA <- ggplot(panel_a_summary, aes(x = threshold, y = median_rate)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "steelblue", alpha = 0.25) +
  geom_line(colour = "steelblue", linewidth = 0.9) +
  geom_point(colour = "steelblue", size = 2.5) +
  geom_vline(xintercept = MIN_TOTAL_FLUX_MM_OPERATIONAL,
             linetype = "dashed", colour = "#d7191c", linewidth = 0.7) +
  scale_x_log10(breaks = THRESH_GRID) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = "A. Site-level unclassified rate (median ± IQR)",
    x     = "MIN_TOTAL_FLUX_MM (mm, log10)",
    y     = "Per-site unclassified rate (flood_type_raw)"
  ) +
  theme_diag()

# ── Panel B: per-type total counts vs. threshold ─────────────────────────────
panel_b_data <- sweep_type_counts %>%
  mutate(flood_type_raw = factor(flood_type_raw,
                                 levels = c("rainfall", "snowmelt", "ROS")))

pB <- ggplot(panel_b_data,
             aes(x = threshold, y = n, colour = flood_type_raw,
                 group = flood_type_raw)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = MIN_TOTAL_FLUX_MM_OPERATIONAL,
             linetype = "dashed", colour = "#d7191c", linewidth = 0.7) +
  scale_x_log10(breaks = THRESH_GRID) +
  scale_y_continuous(labels = comma) +
  scale_colour_manual(values = c("rainfall" = "#e41a1c",
                                 "snowmelt" = "#377eb8",
                                 "ROS"      = "#4daf4a"),
                      name = NULL) +
  labs(
    title = "B. Per-type peak counts (flood_type_raw, all sites)",
    x     = "MIN_TOTAL_FLUX_MM (mm, log10)",
    y     = "Total peaks classified"
  ) +
  theme_diag()

# ── Panel C: % change in 100-yr mixture quantile vs. baseline (thresh = 1.0) ─
pC <- ggplot(q100_pct_change,
             aes(x = factor(threshold, levels = THRESH_GRID), y = pct_change)) +
  geom_hline(yintercept = 0, linetype = "solid",
             colour = "grey50", linewidth = 0.4) +
  geom_boxplot(fill = "steelblue", alpha = 0.6, outlier.size = 0.8) +
  scale_y_continuous(limits = c(-10, 10)) +
  labs(
    title    = "C. 100-yr adaptive mixture quantile vs. threshold = 1.0 baseline",
    subtitle = "Per-site % change. Box = IQR; whiskers = 1.5×IQR; horizontal line = no change.",
    x        = "MIN_TOTAL_FLUX_MM (mm)",
    y        = "% change in Q100 (mixture, GPD, adaptive)"
  ) +
  theme_diag()

p11 <- (pA | pB) / pC + plot_annotation(tag_levels = "A")

save_fig(p11, "fig11_threshold_sensitivity_sweep.png",
         width = 12, height = 8)

# ── Console summary table for the sweep ──────────────────────────────────────
sweep_summary <- tibble(
  threshold       = vapply(sweep_results, `[[`, numeric(1), "threshold"),
  median_uncl     = vapply(sweep_results, `[[`, numeric(1), "median_uncl"),
  iqr_uncl_q25    = vapply(sweep_results, `[[`, numeric(1), "iqr_uncl_q25"),
  iqr_uncl_q75    = vapply(sweep_results, `[[`, numeric(1), "iqr_uncl_q75"),
  n_rainfall      = vapply(sweep_results, `[[`, integer(1), "n_rainfall"),
  n_snowmelt      = vapply(sweep_results, `[[`, integer(1), "n_snowmelt"),
  n_ros           = vapply(sweep_results, `[[`, integer(1), "n_ros")
)

q100_summary <- q100_pct_change %>%
  group_by(threshold) %>%
  summarise(
    n_sites          = n(),
    median_pct_change = round(median(pct_change), 2),
    q25_pct_change    = round(quantile(pct_change, 0.25), 2),
    q75_pct_change    = round(quantile(pct_change, 0.75), 2),
    .groups = "drop"
  )

sweep_summary <- sweep_summary %>%
  left_join(q100_summary, by = "threshold")

message("\n── Threshold sensitivity sweep summary ──")
print(sweep_summary)

# =============================================================================
# CONSOLE SUMMARY TABLE
# =============================================================================
message("\n\n══ Summary tables ══\n")

message("── By elevation bin ──")
site_summary %>%
  filter(!is.na(unclassified_rate), !is.na(elev_bin)) %>%
  group_by(elev_bin) %>%
  summarise(
    n_sites      = n(),
    median_rate  = round(median(unclassified_rate), 3),
    mean_rate    = round(mean(unclassified_rate),   3),
    q25          = round(quantile(unclassified_rate, 0.25), 3),
    q75          = round(quantile(unclassified_rate, 0.75), 3),
    .groups = "drop"
  ) %>% print()

message("\n── By drainage area bin ──")
site_summary %>%
  filter(!is.na(unclassified_rate), !is.na(area_bin)) %>%
  group_by(area_bin) %>%
  summarise(
    n_sites      = n(),
    median_rate  = round(median(unclassified_rate), 3),
    mean_rate    = round(mean(unclassified_rate),   3),
    q25          = round(quantile(unclassified_rate, 0.25), 3),
    q75          = round(quantile(unclassified_rate, 0.75), 3),
    .groups = "drop"
  ) %>% print()

message("\n── By dam presence ──")
site_summary %>%
  filter(!is.na(unclassified_rate), !is.na(dam_group)) %>%
  group_by(dam_group) %>%
  summarise(
    n_sites      = n(),
    median_rate  = round(median(unclassified_rate), 3),
    mean_rate    = round(mean(unclassified_rate),   3),
    .groups = "drop"
  ) %>% print()

message("\n── By reference class ──")
site_summary %>%
  filter(!is.na(unclassified_rate), ref_group != "Unknown") %>%
  group_by(ref_group) %>%
  summarise(
    n_sites      = n(),
    median_rate  = round(median(unclassified_rate), 3),
    mean_rate    = round(mean(unclassified_rate),   3),
    .groups = "drop"
  ) %>% print()

message("\n── By HUC02 ──")
site_summary %>%
  filter(!is.na(unclassified_rate), !is.na(huc02)) %>%
  group_by(huc02) %>%
  summarise(
    n_sites      = n(),
    median_rate  = round(median(unclassified_rate), 3),
    mean_rate    = round(mean(unclassified_rate),   3),
    .groups = "drop"
  ) %>%
  arrange(desc(median_rate)) %>%
  print()

message(sprintf("\nAll figures written to: %s/", OUT_DIR))