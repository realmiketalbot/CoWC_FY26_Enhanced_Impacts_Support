#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Modeled Forcing vs. Observed Runoff Mismatch
# ============================================================================
# PURPOSE:
#   Inspects whether there are peaks with a suspected mismatch between
#   modeled hydrometeorological forcing (WLDAS rain + melt, total_flux) and
#   the observed streamflow response (antecedent_runoff_mm) over the same
#   antecedent window. This is exploratory only — it does not change the
#   classification pipeline, does not touch NOISE_FLOOR_MM, and does not
#   exclude any peaks from the classified output. The goal is to give Mike
#   real numbers to decide on an exclusion rule afterward, if warranted.
#
# Two mismatch directions are flagged, using site-relative medians (baseline
# flux/runoff magnitudes vary a lot across the 206 sites, so fixed global
# cutoffs are not meaningful):
#   - "Model-missed"    : total_flux ~ 0 but antecedent_runoff_mm is large
#                         relative to the site's own peaks.
#   - "Model-overshoot" : total_flux is large relative to the site's own
#                         peaks but antecedent_runoff_mm ~ 0.
#
# INPUTS:
#   _data/pot_peaks_classified_frac_thresh_3class.rds  (Script 6b)
#
# OUTPUTS (diagnostics only — not consumed by any pipeline script):
#   _diagnostics/forcing_response_mismatch/scatter_total_flux_vs_antecedent_runoff.png
#   _diagnostics/forcing_response_mismatch/flagged_model_missed.csv
#   _diagnostics/forcing_response_mismatch/flagged_model_overshoot.csv
#
# Author : Mike Talbot, PhD Student, Colorado State University
# Date   : 2026-07-08
# ============================================================================

library(tidyverse)

# ============================================================================
# CONFIGURATION
# ============================================================================
peaks_file <- "_data/pot_peaks_classified_frac_thresh_3class.rds"
out_dir    <- "_diagnostics/forcing_response_mismatch"

# A flux/runoff value at or below this (mm) is treated as "at or near zero".
NEAR_ZERO_MM <- 1

# "Large" is defined relative to the site's own median of the other axis
# (computed only across that site's data-sufficient peaks), not a fixed
# global cutoff — baseline flux/runoff magnitudes vary widely across sites.
LARGE_REL_TO <- "site median"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

type_colors <- c(rainfall = "#b2182b", snowmelt = "#2166ac",
                 rain_on_snow = "#9970ab", unclassified = "grey60")

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading classified peaks...")
peaks <- readRDS(peaks_file) %>%
  filter(data_sufficient == TRUE)

message(sprintf("Peaks evaluated (data_sufficient): %d across %d sites",
                nrow(peaks), n_distinct(peaks$site_no)))

# ============================================================================
# SCATTER: total_flux vs. antecedent_runoff_mm
# ============================================================================
message("\nBuilding scatter plot...")

axis_max <- max(log1p(peaks$total_flux), log1p(peaks$antecedent_runoff_mm),
                na.rm = TRUE)

p_scatter <- ggplot(peaks,
                    aes(x = log1p(antecedent_runoff_mm), y = log1p(total_flux),
                        color = flood_type_raw)) +
  geom_abline(intercept = 0, slope = 1,
             linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_point(alpha = 0.4, size = 1.2) +
  scale_color_manual(values = type_colors, name = "Flood type\n(pre-demotion)") +
  coord_equal(xlim = c(0, axis_max), ylim = c(0, axis_max)) +
  labs(
    title    = "Modeled forcing vs. observed antecedent runoff, per peak",
    subtitle = "Dashed line: 1:1 (proportional) reference. Points below/right diverge toward model-missed; above/left toward model-overshoot.",
    x = "log1p(antecedent_runoff_mm)  [observed]",
    y = "log1p(total_flux)  [modeled rain + melt]",
    caption = sprintf("data_sufficient == TRUE peaks only, n = %d", nrow(peaks))
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "scatter_total_flux_vs_antecedent_runoff.png"),
      p_scatter, width = 9, height = 8, dpi = 150)
message("  Saved: scatter_total_flux_vs_antecedent_runoff.png")

# ============================================================================
# FLAGGED MISMATCH TABLES (site-relative medians)
# ============================================================================
message("\nFlagging mismatch candidates (site-relative medians)...")

site_medians <- peaks %>%
  group_by(site_no) %>%
  summarise(
    site_median_total_flux           = median(total_flux, na.rm = TRUE),
    site_median_antecedent_runoff_mm = median(antecedent_runoff_mm, na.rm = TRUE),
    .groups = "drop"
  )

peaks_annotated <- peaks %>%
  left_join(site_medians, by = "site_no")

model_missed <- peaks_annotated %>%
  filter(total_flux < NEAR_ZERO_MM,
        antecedent_runoff_mm > site_median_antecedent_runoff_mm) %>%
  mutate(divergence = antecedent_runoff_mm - total_flux) %>%
  arrange(desc(divergence)) %>%
  select(site_no, peak_id, peak_date, peak_flow, flood_type_raw,
        total_flux, antecedent_runoff_mm,
        site_median_antecedent_runoff_mm, divergence)

model_overshoot <- peaks_annotated %>%
  filter(antecedent_runoff_mm < NEAR_ZERO_MM,
        total_flux > site_median_total_flux) %>%
  mutate(divergence = total_flux - antecedent_runoff_mm) %>%
  arrange(desc(divergence)) %>%
  select(site_no, peak_id, peak_date, peak_flow, flood_type_raw,
        total_flux, antecedent_runoff_mm,
        site_median_total_flux, divergence)

write_csv(model_missed, file.path(out_dir, "flagged_model_missed.csv"))
write_csv(model_overshoot, file.path(out_dir, "flagged_model_overshoot.csv"))
message("  Saved: flagged_model_missed.csv")
message("  Saved: flagged_model_overshoot.csv")

# ============================================================================
# SUMMARY COUNTS
# ============================================================================
message("\n=== SUMMARY ===")

n_total <- nrow(peaks)
message(sprintf("Total peaks evaluated (data_sufficient == TRUE): %d", n_total))

message(sprintf(
  "Model-missed candidates    : %d (%.2f%%), %d distinct sites",
  nrow(model_missed),
  100 * nrow(model_missed) / n_total,
  n_distinct(model_missed$site_no)
))
message(sprintf(
  "Model-overshoot candidates : %d (%.2f%%), %d distinct sites",
  nrow(model_overshoot),
  100 * nrow(model_overshoot) / n_total,
  n_distinct(model_overshoot$site_no)
))

message("\nModel-missed by flood_type_raw:")
model_missed %>% count(flood_type_raw, name = "n_flagged") %>% print(n = Inf)

message("\nModel-overshoot by flood_type_raw:")
model_overshoot %>% count(flood_type_raw, name = "n_flagged") %>% print(n = Inf)

message("\nDone.")
