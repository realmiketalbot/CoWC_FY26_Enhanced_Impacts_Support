#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Antecedent Window Mass Balance
# ============================================================================
# For each POT peak, compares the WLDAS water budget over the antecedent
# search_radius window against observed streamflow runoff depth over the
# same window.
#
# Water balance over the antecedent window:
#
#   Inputs  : rain + snow (total precip) + ΔSWE (snowpack release, if negative)
#   Outputs : ET + Qs (surface runoff) + Qsb (subsurface runoff)
#   Storage : ΔSoilMoisture + ΔSWE (if positive, snowpack is building)
#
#   Balance equation (all terms in mm):
#     P + max(0, -ΔSWE) ≈ ET + Qs + Qsb + ΔSoilMoi + max(0, ΔSWE) + Q_obs
#
#   Rearranged as residual:
#     residual = (P - ΔSWE - ET - Qs - Qsb - ΔSoilMoi) - Q_obs
#
#   A residual near zero means the WLDAS budget closes against observed
#   runoff over the antecedent window. Systematic bias indicates either:
#     - The antecedent window is too short (missing upstream forcing)
#     - Regulated flow at the gauge
#     - ET or storage terms are poorly estimated in WLDAS
#
# NOTE: Soil moisture in WLDAS is volumetric (m3/m3). To convert to mm,
# multiply by layer depth in mm:
#   0-10 cm  = 100 mm
#   10-40 cm = 300 mm
#   40-100cm = 600 mm
#   Total    = 1000 mm equivalent depth
#
# INPUTS:
#   - _data/pot_peaks_all.rds
#   - _data/WLDAS/processed/daily_by_date/*.parquet
#   - _data/site_metadata.rds
#
# OUTPUTS:
#   - _data/mass_balance.rds
#   - _diagnostics/mass_balance/mass_balance.png
# ============================================================================

library(tidyverse)
library(arrow)
library(patchwork)

# ============================================================================
# CONFIGURATION
# ============================================================================
peaks_file     <- "_data/pot_peaks_all.rds"
parquet_dir    <- "_data/WLDAS/processed/daily_by_date"
meta_file      <- "_data/site_metadata.rds"
output_file    <- "_data/mass_balance.rds"
out_dir        <- "_diagnostics/mass_balance"

# Soil moisture layer depths in mm (for converting m3/m3 to mm)
SOILM_DEPTH_MM <- c(
  soilm_0_10cm   = 100,
  soilm_10_40cm  = 300,
  soilm_40_100cm = 600
)

# Minimum antecedent days required (fraction of search_radius)
# Peaks where antecedent_days_used < MIN_DAYS_FRAC * search_radius are flagged
MIN_DAYS_FRAC <- 0.75

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD PEAKS AND METADATA
# ============================================================================
message("Loading peaks and metadata...")
peaks <- readRDS(peaks_file)
meta  <- readRDS(meta_file)

# One row per peak (not per peak-day)
peak_index <- peaks %>%
  distinct(site_no, peak_id, peak_date, search_radius,
           antecedent_runoff_mm, antecedent_days_used) %>%
  mutate(
    window_start = peak_date - search_radius,
    window_end   = peak_date - 1L,          # antecedent only, exclude peak day
    days_complete = antecedent_days_used >= MIN_DAYS_FRAC * search_radius
  )

message(sprintf("Peaks: %d total, %d with complete antecedent windows",
                nrow(peak_index), sum(peak_index$days_complete)))

# ============================================================================
# IDENTIFY PARQUET FILES NEEDED
# ============================================================================
# All dates spanned by any antecedent window
all_dates_needed <- peak_index %>%
  rowwise() %>%
  reframe(date = seq(window_start, window_end, by = "day")) %>%
  distinct(date)

parquet_files_all <- tibble(
  path = list.files(parquet_dir, pattern = "\\.parquet$", full.names = TRUE)
) %>%
  mutate(date = as.Date(str_extract(basename(path), "\\d{8}"), "%Y%m%d")) %>%
  filter(date %in% all_dates_needed$date)

message(sprintf("Loading %d Parquet files covering antecedent windows...",
                nrow(parquet_files_all)))

# ============================================================================
# LOAD AND STACK WLDAS DATA
# ============================================================================
wldas_daily <- map_dfr(parquet_files_all$path, function(f) {
  tryCatch(
    read_parquet(f) %>%
      rename(site_no = GAGE_ID) %>%
      mutate(site_no = as.character(site_no)),
    error = function(e) {
      warning(sprintf("Failed to read %s: %s", basename(f), e$message))
      NULL
    }
  )
})

message(sprintf("Loaded %d site-day rows from WLDAS", nrow(wldas_daily)))

# ============================================================================
# COMPUTE ANTECEDENT FLUX INTEGRALS PER PEAK
# ============================================================================
# Join each WLDAS day to its peak(s) via date range, then sum fluxes.
# A given WLDAS day may contribute to multiple peaks if windows overlap
# (rare given the independence criterion, but possible at large basins).
message("Integrating WLDAS fluxes over antecedent windows...")

wldas_site_dates <- wldas_daily %>%
  select(site_no, date,
         rainf_flux, snowf_flux, qsm_flux,
         qs_flux, qsb_flux, evap_flux,
         swe,
         soilm_0_10cm, soilm_10_40cm, soilm_40_100cm)

# For each peak, pull matching site+date rows and summarise
flux_integrals <- peak_index %>%
  select(site_no, peak_id, peak_date, window_start, window_end) %>%
  left_join(wldas_site_dates, by = "site_no") %>%
  filter(date >= window_start, date <= window_end) %>%
  group_by(site_no, peak_id, peak_date, window_start, window_end) %>%
  summarise(
    # Integrated fluxes (mm): sum of daily mm values over window
    rain_mm      = sum(rainf_flux,  na.rm = TRUE),
    snow_mm      = sum(snowf_flux,  na.rm = TRUE),
    melt_mm      = sum(qsm_flux,    na.rm = TRUE),
    qs_mm        = sum(qs_flux,     na.rm = TRUE),
    qsb_mm       = sum(qsb_flux,    na.rm = TRUE),
    et_mm        = sum(evap_flux,   na.rm = TRUE),
    # SWE endpoints for storage change
    swe_start    = first(swe[order(date)]),
    swe_end      = last(swe[order(date)]),
    # Soil moisture endpoints (m3/m3) — convert to mm below
    sm1_start    = first(soilm_0_10cm[order(date)]),
    sm1_end      = last(soilm_0_10cm[order(date)]),
    sm2_start    = first(soilm_10_40cm[order(date)]),
    sm2_end      = last(soilm_10_40cm[order(date)]),
    sm3_start    = first(soilm_40_100cm[order(date)]),
    sm3_end      = last(soilm_40_100cm[order(date)]),
    n_wldas_days = n(),
    .groups      = "drop"
  )

# ============================================================================
# COMPUTE STORAGE CHANGES AND MASS BALANCE
# ============================================================================
message("Computing storage changes and mass balance residuals...")

mass_balance <- flux_integrals %>%
  mutate(
    # ΔSWE (mm): positive = snowpack building, negative = releasing
    delta_swe_mm  = swe_end - swe_start,

    # ΔSoil moisture (mm): convert volumetric to depth then difference
    delta_soilm_mm = (sm1_end - sm1_start) * SOILM_DEPTH_MM["soilm_0_10cm"]  +
                     (sm2_end - sm2_start) * SOILM_DEPTH_MM["soilm_10_40cm"] +
                     (sm3_end - sm3_start) * SOILM_DEPTH_MM["soilm_40_100cm"],

    # Total precipitation input (rain + snow, both phases)
    precip_mm     = rain_mm + snow_mm,

    # WLDAS-modelled total output (ET + surface + subsurface runoff)
    wldas_q_mm    = qs_mm + qsb_mm
  ) %>%
  left_join(
    peak_index %>% select(site_no, peak_id, antecedent_runoff_mm,
                          antecedent_days_used, search_radius, days_complete),
    by = c("site_no", "peak_id")
  ) %>%
  mutate(
    # Full balance residual:
    # (precip - delta_swe - ET - WLDAS_Q - delta_soilm) - observed_Q
    # Near zero = budget closes. Positive = WLDAS overestimates supply.
    # Negative = observed runoff exceeds WLDAS supply (common for regulated sites
    # or where snowpack release precedes the window).
    residual_mm = (precip_mm - delta_swe_mm - et_mm -
                   wldas_q_mm - delta_soilm_mm) - antecedent_runoff_mm,

    # Runoff ratio: observed runoff / total precip input
    runoff_ratio = antecedent_runoff_mm / pmax(precip_mm, 1),

    # Season for faceting
    month        = month(peak_date),
    season       = case_when(
      month %in% 3:5  ~ "Spring (MAM)",
      month %in% 6:8  ~ "Summer (JJA)",
      month %in% 9:11 ~ "Fall (SON)",
      TRUE            ~ "Winter (DJF)"
    )
  ) %>%
  left_join(meta %>% select(site_no, drain_sqkm, elev_mean_m), by = "site_no")

saveRDS(mass_balance, output_file)
message(sprintf("Saved mass balance to %s (%d peaks)", output_file, nrow(mass_balance)))

# Summary statistics
mb_complete <- mass_balance %>% filter(days_complete)
message(sprintf("\n=== MASS BALANCE SUMMARY (complete windows only, n=%d) ===",
                nrow(mb_complete)))
message(sprintf("  Median residual:      %+.1f mm", median(mb_complete$residual_mm, na.rm = TRUE)))
message(sprintf("  Mean residual:        %+.1f mm", mean(mb_complete$residual_mm,   na.rm = TRUE)))
message(sprintf("  SD of residual:       %.1f mm",  sd(mb_complete$residual_mm,     na.rm = TRUE)))
message(sprintf("  Median runoff ratio:  %.2f",     median(mb_complete$runoff_ratio, na.rm = TRUE)))

# ============================================================================
# FIGURES
# ============================================================================

# Colour consistent with rest of project
season_colors <- c(
  "Spring (MAM)" = "#4dac26",
  "Summer (JJA)" = "#d73027",
  "Fall (SON)"   = "#f1a340",
  "Winter (DJF)" = "#2166ac"
)

# --- Fig 1: Residual distribution by season ---
p1 <- mb_complete %>%
  ggplot(aes(x = residual_mm, fill = season)) +
  geom_histogram(bins = 60, color = NA, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey20") +
  facet_wrap(~ season, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = season_colors) +
  labs(
    title    = "Mass balance residual by season",
    subtitle = "Residual = (WLDAS supply - storage change) - observed runoff. Near zero = budget closes.",
    x        = "Residual (mm)",
    y        = "Number of peaks",
    fill     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# --- Fig 2: Observed vs WLDAS-implied runoff ---
p2 <- mb_complete %>%
  ggplot(aes(x = wldas_q_mm, y = antecedent_runoff_mm, color = season)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey30") +
  scale_color_manual(values = season_colors) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  labs(
    title    = "Observed vs WLDAS-modelled runoff (log scale)",
    subtitle = "Dashed line = 1:1. Points above line = observed exceeds WLDAS.",
    x        = "WLDAS Qs + Qsb (mm)",
    y        = "Observed antecedent runoff (mm)",
    color    = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# --- Fig 3: Residual vs search_radius (window length) ---
p3 <- mb_complete %>%
  ggplot(aes(x = search_radius, y = residual_mm)) +
  geom_point(alpha = 0.2, size = 0.7, color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20") +
  geom_smooth(method = "loess", span = 0.5,
              color = "#762a83", fill = "#762a83", alpha = 0.2) +
  labs(
    title    = "Mass balance residual vs antecedent window length",
    subtitle = "Systematic trend would suggest window is too short or too long.",
    x        = "Search radius (days)",
    y        = "Residual (mm)"
  ) +
  theme_bw(base_size = 11)

# --- Fig 4: Runoff ratio vs drainage area ---
p4 <- mb_complete %>%
  filter(runoff_ratio > 0, runoff_ratio < 5) %>%   # trim extreme outliers
  ggplot(aes(x = drain_sqkm, y = runoff_ratio, color = season)) +
  geom_point(alpha = 0.25, size = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey30") +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = season_colors) +
  labs(
    title    = "Runoff ratio vs drainage area",
    subtitle = "Ratio > 1 suggests regulated flow or snowpack release outside window.",
    x        = "Drainage area (km²)",
    y        = "Runoff ratio (observed Q / precip)",
    color    = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# ============================================================================
# COMBINE AND SAVE
# ============================================================================
combined <- (p1 / (p2 + p3) / p4) +
  plot_annotation(
    title   = "Antecedent Window Mass Balance Diagnostic",
    caption = sprintf(
      paste0("Antecedent window = search_radius days prior to peak (Livneh et al. independence criterion).\n",
             "Complete windows only (>= %.0f%% of search_radius days with WLDAS data). n = %d peaks."),
      MIN_DAYS_FRAC * 100, nrow(mb_complete)
    )
  )

out_file <- file.path(out_dir, "mass_balance.png")
ggsave(out_file, combined, width = 12, height = 16, dpi = 150)
message(sprintf("Saved figure to %s", out_file))