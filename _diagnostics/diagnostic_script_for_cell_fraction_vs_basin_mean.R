#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Basin-Mean vs Cell-Fraction Classification Comparison
# ============================================================================
# Evaluates whether the cell-fraction approach recovers peaks that are
# unclassified under the basin-mean Davenport scheme, and whether the
# type fractions show better signal separation than basin-mean totals.
#
# Key questions:
#   1. How many basin-mean unclassified peaks get a dominant type under
#      the cell-fraction scheme, and what types do they receive?
#   2. How do the frac_* distributions look — is there clear separation
#      between types, or is frac_mixed dominating?
#   3. Does recovery rate vary with drainage area or elevation?
#   4. For peaks classified by both schemes, do they agree?
#
# INPUTS:
#   - _data/pot_peaks_classified_mean.rds   : basin-mean Davenport output
#   - _data/pot_peaks_cell_classified.rds   : cell-fraction output
#   - _data/site_metadata.rds              : site attributes
#
# OUTPUTS:
#   - Console summaries
#   - _data/diag_cellfrac_vs_basin_mean.png : multi-panel figure
# ============================================================================

library(tidyverse)
library(patchwork)

# ============================================================================
# CONFIGURATION
# ============================================================================
basin_mean_file  <- "_data/pot_peaks_classified_mean.rds"
cell_frac_file   <- "_data/pot_peaks_cell_classified.rds"
site_meta_file   <- "_data/site_metadata.rds"
output_plot      <- "_diagnostics/threshold_calibration/diag_cellfrac_vs_basin_mean.png"

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading data...")

site_meta  <- readRDS(site_meta_file) %>%
  select(site_no, drain_sqkm, elev_mean_m) %>%
  mutate(
    area_quartile = cut(drain_sqkm,
                        breaks = quantile(drain_sqkm, probs = seq(0, 1, 0.25),
                                          na.rm = TRUE),
                        include.lowest = TRUE, dig.lab = 4,
                        labels = c("Q1 (small)", "Q2", "Q3", "Q4 (large)")),
    elev_quartile = cut(elev_mean_m,
                        breaks = quantile(elev_mean_m, probs = seq(0, 1, 0.25),
                                          na.rm = TRUE),
                        include.lowest = TRUE, dig.lab = 4,
                        labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"))
  )

basin_mean <- readRDS(basin_mean_file) %>%
  select(site_no, peak_id, peak_date, peak_flow,
         bm_type      = flood_type,
         bm_type_raw  = flood_type_raw,
         total_rain_mm, total_melt_mm,
         complete_window) %>%
  mutate(peak_id = as.character(peak_id))

cell_frac  <- readRDS(cell_frac_file) %>%
  select(site_no, peak_id, peak_date, peak_flow,
         frac_rainfall, frac_snowmelt, frac_ros, frac_mixed,
         dominant_type, dominant_frac,
         cf_complete = complete_window) %>%
  mutate(peak_id = as.character(peak_id))

# ============================================================================
# JOIN
# ============================================================================
combined <- inner_join(basin_mean, cell_frac,
                       by = c("site_no", "peak_id", "peak_date", "peak_flow")) %>%
  inner_join(site_meta, by = "site_no") %>%
  filter(complete_window, cf_complete)   # require complete windows in both

message(sprintf("Peaks with complete windows in both schemes: %d", nrow(combined)))

# ============================================================================
# QUESTION 1: Recovery of basin-mean unclassified peaks
# ============================================================================
message("\n--- Q1: Cell-fraction recovery of basin-mean unclassified peaks ---")

recovery <- combined %>%
  filter(bm_type == "unclassified") %>%
  count(dominant_type) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(desc(n))

message("Basin-mean UNCLASSIFIED peaks -> cell-fraction dominant type:")
print(recovery)

# Recovery rate by area quartile
recovery_by_area <- combined %>%
  filter(bm_type == "unclassified") %>%
  group_by(area_quartile) %>%
  summarise(
    n_total     = n(),
    n_recovered = sum(!is.na(dominant_type) & dominant_type != "mixed"),
    pct_recovered = round(100 * n_recovered / n_total, 1),
    .groups = "drop"
  )

message("\nRecovery rate by drainage area quartile:")
print(recovery_by_area)

# ============================================================================
# QUESTION 2: Distribution of frac_* values by basin-mean type
# ============================================================================
# Reshape to long for plotting
frac_long <- combined %>%
  select(site_no, peak_id, bm_type, area_quartile,
         frac_rainfall, frac_snowmelt, frac_ros, frac_mixed) %>%
  pivot_longer(cols = starts_with("frac_"),
               names_to  = "frac_type",
               values_to = "fraction") %>%
  mutate(
    frac_type = str_remove(frac_type, "frac_"),
    frac_type = factor(frac_type,
                       levels = c("rainfall", "snowmelt", "ros", "mixed"))
  )

# ============================================================================
# QUESTION 3: Agreement between schemes for classified peaks
# ============================================================================
message("\n--- Q3: Agreement between schemes (classified peaks only) ---")

# Normalise cell-fraction dominant type labels to match basin-mean
agreement <- combined %>%
  filter(bm_type != "unclassified") %>%
  mutate(
    cf_type_norm = case_when(
      dominant_type == "rainfall"    ~ "rainfall",
      dominant_type == "snowmelt"    ~ "snowmelt",
      dominant_type == "rain_on_snow"~ "rain_on_snow",
      TRUE                           ~ "mixed/unclassified"
    ),
    bm_type_norm = bm_type,
    agree = cf_type_norm == bm_type_norm
  )

message(sprintf("Agreement rate (both classified, same type): %.1f%%",
                100 * mean(agreement$agree, na.rm = TRUE)))

confusion <- agreement %>%
  count(bm_type_norm, cf_type_norm) %>%
  group_by(bm_type_norm) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()

message("\nConfusion matrix (basin-mean type vs cell-fraction type):")
print(confusion)

# ============================================================================
# FIGURES
# ============================================================================

# --- Fig 1: Recovery counts by type ---
p1 <- combined %>%
  filter(bm_type == "unclassified", !is.na(dominant_type)) %>%
  mutate(dominant_type = replace_na(dominant_type, "no dominant type")) %>%
  count(dominant_type) %>%
  mutate(dominant_type = fct_reorder(dominant_type, n)) %>%
  ggplot(aes(x = dominant_type, y = n, fill = dominant_type)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c(
    rainfall     = "#4575b4",
    snowmelt     = "#74add1",
    rain_on_snow = "#f46d43",
    mixed        = "#d9d9d9"
  )) +
  coord_flip() +
  labs(
    title    = "Cell-fraction type for basin-mean unclassified peaks",
    subtitle = "How many unclassified peaks get a dominant type?",
    x        = NULL,
    y        = "Number of peaks"
  ) +
  theme_bw(base_size = 11)

# --- Fig 2: Recovery rate by drainage area quartile ---
p2 <- recovery_by_area %>%
  ggplot(aes(x = area_quartile, y = pct_recovered, fill = area_quartile)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = sprintf("n=%d\n%.1f%%", n_total, pct_recovered)),
            vjust = -0.3, size = 3) +
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    title    = "Recovery rate by drainage area quartile",
    subtitle = "% of basin-mean unclassified peaks with non-mixed cell-fraction type",
    x        = "Drainage area quartile",
    y        = "Recovery rate (%)"
  ) +
  theme_bw(base_size = 11)

# --- Fig 3: frac_* distributions by basin-mean type ---
p3 <- frac_long %>%
  filter(bm_type != "unclassified") %>%
  ggplot(aes(x = fraction, fill = frac_type)) +
  geom_histogram(bins = 40, alpha = 0.8, position = "identity") +
  facet_grid(bm_type ~ frac_type, scales = "free_y") +
  scale_fill_manual(values = c(
    rainfall = "#4575b4",
    snowmelt = "#74add1",
    ros      = "#f46d43",
    mixed    = "#d9d9d9"
  )) +
  labs(
    title    = "Cell-fraction distributions by basin-mean classified type",
    subtitle = "Rows = basin-mean type; columns = cell fraction component",
    x        = "Area-weighted fraction",
    y        = "Count",
    fill     = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

# --- Fig 4: dominant_frac distribution, classified vs unclassified ---
p4 <- combined %>%
  mutate(bm_classified = if_else(bm_type == "unclassified",
                                 "basin-mean unclassified",
                                 "basin-mean classified")) %>%
  filter(!is.na(dominant_frac)) %>%
  ggplot(aes(x = dominant_frac, fill = bm_classified)) +
  geom_histogram(bins = 40, alpha = 0.7, position = "dodge") +
  geom_vline(xintercept = 0.50, linetype = "dashed", color = "red") +
  annotate("text", x = 0.52, y = Inf, label = "50% threshold",
           hjust = 0, vjust = 1.5, color = "red", size = 3) +
  scale_fill_manual(values = c(
    "basin-mean classified"   = "#4575b4",
    "basin-mean unclassified" = "#d73027"
  )) +
  labs(
    title    = "Cell-fraction dominant type confidence",
    subtitle = "Distribution of dominant area-weighted fraction",
    x        = "Dominant type fraction",
    y        = "Count",
    fill     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# ============================================================================
# COMBINE AND SAVE
# ============================================================================
combined_plot <- (p1 + p2) / p3 / p4 +
  plot_annotation(
    title   = "Cell-Fraction vs Basin-Mean Classification Comparison",
    caption = "Only peaks with complete 8-day windows in both schemes included."
  )

ggsave(output_plot, combined_plot, width = 12, height = 14, dpi = 150)
message(sprintf("\nSaved plot to %s", output_plot))