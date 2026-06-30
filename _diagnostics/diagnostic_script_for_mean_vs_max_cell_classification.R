#!/usr/bin/env Rscript

# ============================================================================
# Classification Comparison: Basin-Mean vs Max-Cell WLDAS Values
# ============================================================================
# Compares flood peak classification outcomes between two modes of script 5:
#   "mean"     : basin-mean antecedent forcing (Davenport et al. approach)
#   "max_cell" : max-cell antecedent forcing (spatial aggregation test)
#
# Key questions:
#   1. Overall: how many peaks move from unclassified to classified?
#   2. Are there regressions — peaks that lose classification under max_cell?
#   3. At the site level: which sites gain flood types, which lose them?
#   4. Do classification gains concentrate in particular flood types,
#      elevations, or drainage areas?
#   5. What is the net change in sites qualifying for per-type FFA
#      (i.e., reaching the Yu et al. 20-event minimum)?
#
# INPUTS:
#   - _data/pot_peaks_classified_mean.rds
#   - _data/pot_peaks_classified_max_cell.rds
#   - _data/site_metadata.rds
#
# OUTPUTS:
#   - _diagnostics/mean_vs_maxcell/fig01_overall_classification.png
#   - _diagnostics/mean_vs_maxcell/fig02_peak_transitions.png
#   - _diagnostics/mean_vs_maxcell/fig03_site_type_gains.png
#   - _diagnostics/mean_vs_maxcell/fig04_gains_by_elevation.png
#   - _diagnostics/mean_vs_maxcell/fig05_gains_by_area.png
#   - _diagnostics/mean_vs_maxcell/fig06_ffa_eligible_sites.png
#   - _diagnostics/mean_vs_maxcell/comparison_summary.rds
# ============================================================================

library(tidyverse)

# ============================================================================
# CONFIGURATION
# ============================================================================
peaks_mean_file    <- "_data/pot_peaks_classified_mean.rds"
peaks_maxcell_file <- "_data/pot_peaks_classified_max_cell.rds"
meta_file          <- "_data/site_metadata.rds"
out_dir            <- "_diagnostics/mean_vs_maxcell"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

FLOOD_TYPES <- c("rainfall", "snowmelt", "rain_on_snow")

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading classified peaks...")
mean_peaks    <- readRDS(peaks_mean_file)    %>% mutate(peak_date = as.Date(peak_date))
maxcell_peaks <- readRDS(peaks_maxcell_file) %>% mutate(peak_date = as.Date(peak_date))
meta          <- readRDS(meta_file)

message(sprintf("Mean mode:     %d peaks from %d sites",
                nrow(mean_peaks),    n_distinct(mean_peaks$site_no)))
message(sprintf("Max-cell mode: %d peaks from %d sites",
                nrow(maxcell_peaks), n_distinct(maxcell_peaks$site_no)))

# ============================================================================
# PEAK-LEVEL COMPARISON
# ============================================================================
# Join on site_no + peak_id, which uniquely identifies each peak event.
# flood_type_raw preserves the threshold result before sample size filter.

message("Joining peak-level classifications...")

peak_compare <- mean_peaks %>%
  select(site_no, peak_id, peak_date, peak_flow,
         flood_type_mean     = flood_type,
         flood_type_raw_mean = flood_type_raw) %>%
  inner_join(
    maxcell_peaks %>%
      select(site_no, peak_id,
             flood_type_maxcell     = flood_type,
             flood_type_raw_maxcell = flood_type_raw),
    by = c("site_no", "peak_id")
  )

message(sprintf("Joined %d peaks present in both datasets", nrow(peak_compare)))

# Classify each peak's transition
peak_compare <- peak_compare %>%
  mutate(
    transition = case_when(
      flood_type_mean == flood_type_maxcell                           ~ "no change",
      flood_type_mean == "unclassified" & flood_type_maxcell != "unclassified" ~ "gained classification",
      flood_type_mean != "unclassified" & flood_type_maxcell == "unclassified" ~ "lost classification",
      flood_type_mean != flood_type_maxcell &
        flood_type_mean != "unclassified" &
        flood_type_maxcell != "unclassified"                          ~ "changed type",
      TRUE                                                            ~ "other"
    ),
    # For gained peaks, what type did they gain?
    gained_type = if_else(transition == "gained classification",
                          flood_type_maxcell, NA_character_),
    # For lost peaks, what type did they lose?
    lost_type   = if_else(transition == "lost classification",
                          flood_type_mean, NA_character_)
  )

# ============================================================================
# CONSOLE SUMMARY: PEAK-LEVEL
# ============================================================================
message("\n=== PEAK-LEVEL TRANSITION SUMMARY ===")

transition_summary <- peak_compare %>%
  count(transition) %>%
  mutate(pct = round(100 * n / sum(n), 1))
print(transition_summary)

message("\nGained classifications by type:")
peak_compare %>%
  filter(transition == "gained classification") %>%
  count(gained_type) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

message("\nLost classifications by type:")
peak_compare %>%
  filter(transition == "lost classification") %>%
  count(lost_type) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

message("\nChanged type (mean -> maxcell):")
peak_compare %>%
  filter(transition == "changed type") %>%
  count(flood_type_mean, flood_type_maxcell) %>%
  arrange(desc(n)) %>%
  print()

# ============================================================================
# SITE-LEVEL COMPARISON
# ============================================================================
message("\nComputing site-level summaries...")

# Per-site, per-type counts under each mode
site_type_counts <- function(peaks, mode_label) {
  peaks %>%
    filter(flood_type %in% FLOOD_TYPES) %>%
    count(site_no, flood_type, name = "n_peaks") %>%
    mutate(mode = mode_label)
}

site_counts_mean    <- site_type_counts(mean_peaks,    "mean")
site_counts_maxcell <- site_type_counts(maxcell_peaks, "max_cell")

# Full grid: all sites x all flood types x both modes
all_sites <- union(mean_peaks$site_no, maxcell_peaks$site_no)
full_grid <- expand_grid(
  site_no    = all_sites,
  flood_type = FLOOD_TYPES,
  mode       = c("mean", "max_cell")
)

site_counts_full <- full_grid %>%
  left_join(bind_rows(site_counts_mean, site_counts_maxcell),
            by = c("site_no", "flood_type", "mode")) %>%
  replace_na(list(n_peaks = 0L))

# Pivot wide for comparison
site_wide <- site_counts_full %>%
  pivot_wider(names_from = mode, values_from = n_peaks,
              names_prefix = "n_") %>%
  mutate(
    delta          = n_max_cell - n_mean,
    qualifies_mean    = n_mean    >= 20,
    qualifies_maxcell = n_max_cell >= 20,
    ffa_change = case_when(
      !qualifies_mean & qualifies_maxcell ~ "gained FFA eligibility",
      qualifies_mean & !qualifies_maxcell ~ "lost FFA eligibility",
      qualifies_mean & qualifies_maxcell  ~ "eligible in both",
      TRUE                                ~ "ineligible in both"
    )
  ) %>%
  left_join(meta %>% select(site_no, elev_mean_m, drain_sqkm, huc02),
            by = "site_no")

# ============================================================================
# CONSOLE SUMMARY: SITE-LEVEL FFA ELIGIBILITY
# ============================================================================
message("\n=== SITE-LEVEL FFA ELIGIBILITY (>= 20 events per type) ===")

ffa_summary <- site_wide %>%
  count(flood_type, ffa_change) %>%
  pivot_wider(names_from = ffa_change, values_from = n, values_fill = 0L)
print(ffa_summary)

message("\nNet change in FFA-eligible sites per flood type:")
site_wide %>%
  group_by(flood_type) %>%
  summarise(
    eligible_mean    = sum(qualifies_mean),
    eligible_maxcell = sum(qualifies_maxcell),
    net_change       = eligible_maxcell - eligible_mean,
    .groups = "drop"
  ) %>%
  print()

# ============================================================================
# FIG 01: Overall classification breakdown — stacked bar, both modes
# ============================================================================
message("\nGenerating figures...")

fig01_data <- bind_rows(
  mean_peaks    %>% mutate(mode = "Basin mean"),
  maxcell_peaks %>% mutate(mode = "Max cell")
) %>%
  mutate(
    flood_type = factor(flood_type,
                        levels = c("rainfall", "snowmelt", "rain_on_snow", "unclassified"),
                        labels = c("Rainfall", "Snowmelt", "Rain-on-Snow", "Unclassified"))
  ) %>%
  count(mode, flood_type) %>%
  group_by(mode) %>%
  mutate(pct = 100 * n / sum(n))

fig01 <- ggplot(fig01_data,
                aes(x = mode, y = pct, fill = flood_type)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_stack(vjust = 0.5),
            size = 3.5, colour = "white", fontface = "bold") +
  scale_fill_manual(values = c(
    "Rainfall"      = "#e41a1c",
    "Snowmelt"      = "#377eb8",
    "Rain-on-Snow"  = "#4daf4a",
    "Unclassified"  = "#cccccc"
  )) +
  labs(title    = "Overall classification breakdown: basin mean vs max cell",
       x        = NULL,
       y        = "Percentage of peaks (%)",
       fill     = "Flood type") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "fig01_overall_classification.png"),
       fig01, width = 6, height = 5, dpi = 150)

# ============================================================================
# FIG 02: Peak-level transition counts
# ============================================================================
transition_levels <- c("no change", "gained classification",
                       "lost classification", "changed type", "other")
transition_colors <- c(
  "no change"             = "#aaaaaa",
  "gained classification" = "#2166ac",
  "lost classification"   = "#d73027",
  "changed type"          = "#f46d43",
  "other"                 = "#cccccc"
)

fig02_data <- peak_compare %>%
  count(transition) %>%
  mutate(transition = factor(transition, levels = transition_levels))

fig02 <- ggplot(fig02_data, aes(x = transition, y = n, fill = transition)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = scales::comma(n)), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = transition_colors) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +
  labs(title    = "Peak-level classification transitions: mean -> max cell",
       x        = NULL,
       y        = "Number of peaks") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(out_dir, "fig02_peak_transitions.png"),
       fig02, width = 7, height = 5, dpi = 150)

# ============================================================================
# FIG 03: Site-level FFA eligibility changes by flood type
# ============================================================================
fig03_data <- site_wide %>%
  filter(ffa_change %in% c("gained FFA eligibility", "lost FFA eligibility",
                           "eligible in both")) %>%
  count(flood_type, ffa_change) %>%
  mutate(
    flood_type = factor(flood_type,
                        levels = c("rainfall", "snowmelt", "rain_on_snow"),
                        labels = c("Rainfall", "Snowmelt", "Rain-on-Snow")),
    ffa_change = factor(ffa_change,
                        levels = c("eligible in both", "gained FFA eligibility",
                                   "lost FFA eligibility"))
  )

fig03 <- ggplot(fig03_data,
                aes(x = flood_type, y = n, fill = ffa_change)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c(
    "eligible in both"       = "#aaaaaa",
    "gained FFA eligibility" = "#2166ac",
    "lost FFA eligibility"   = "#d73027"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title    = "Site-level FFA eligibility changes by flood type",
       subtitle = "FFA eligible = >= 20 classified events of that type (Yu et al.)",
       x        = NULL,
       y        = "Number of sites",
       fill     = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "fig03_site_type_gains.png"),
       fig03, width = 7, height = 5, dpi = 150)

# ============================================================================
# FIG 04: Classification gains by elevation
# ============================================================================
fig04_data <- peak_compare %>%
  left_join(meta %>% select(site_no, elev_mean_m), by = "site_no") %>%
  filter(transition %in% c("gained classification", "lost classification",
                           "no change")) %>%
  filter(!is.na(elev_mean_m))

fig04 <- ggplot(fig04_data,
                aes(x = elev_mean_m, fill = transition)) +
  geom_histogram(binwidth = 100, position = "stack", colour = "white",
                 linewidth = 0.2) +
  scale_fill_manual(values = transition_colors) +
  labs(title    = "Classification transitions by mean basin elevation",
       x        = "Mean basin elevation (m)",
       y        = "Number of peaks",
       fill     = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "fig04_gains_by_elevation.png"),
       fig04, width = 8, height = 5, dpi = 150)

# ============================================================================
# FIG 05: Classification gains by drainage area
# ============================================================================
fig05_data <- peak_compare %>%
  left_join(meta %>% select(site_no, drain_sqkm), by = "site_no") %>%
  filter(transition %in% c("gained classification", "lost classification",
                           "no change")) %>%
  filter(!is.na(drain_sqkm), drain_sqkm > 0)

fig05 <- ggplot(fig05_data,
                aes(x = drain_sqkm, fill = transition)) +
  geom_histogram(bins = 40, position = "stack", colour = "white",
                 linewidth = 0.2) +
  scale_x_log10(labels = scales::comma) +
  scale_fill_manual(values = transition_colors) +
  labs(title    = "Classification transitions by drainage area",
       x        = "Drainage area (km², log scale)",
       y        = "Number of peaks",
       fill     = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "fig05_gains_by_area.png"),
       fig05, width = 8, height = 5, dpi = 150)

# ============================================================================
# FIG 06: FFA-eligible site counts — before and after, by flood type
# ============================================================================
fig06_data <- site_wide %>%
  group_by(flood_type) %>%
  summarise(
    `Basin mean` = sum(qualifies_mean),
    `Max cell`   = sum(qualifies_maxcell),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c("Basin mean", "Max cell"),
               names_to = "mode", values_to = "n_sites") %>%
  mutate(
    flood_type = factor(flood_type,
                        levels = c("rainfall", "snowmelt", "rain_on_snow"),
                        labels = c("Rainfall", "Snowmelt", "Rain-on-Snow")),
    mode = factor(mode, levels = c("Basin mean", "Max cell"))
  )

fig06 <- ggplot(fig06_data,
                aes(x = flood_type, y = n_sites, fill = mode)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = n_sites),
            position = position_dodge(width = 0.6),
            vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("Basin mean" = "#969696", "Max cell" = "#2166ac")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title    = "Sites with >= 20 events per flood type (FFA eligible)",
       subtitle = "Basin mean vs max cell classification",
       x        = NULL,
       y        = "Number of sites",
       fill     = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "fig06_ffa_eligible_sites.png"),
       fig06, width = 6, height = 5, dpi = 150)

# ============================================================================
# FIG 07: Changed-type peaks — distribution of max-cell melt values
# ============================================================================
# Focuses on rainfall -> rain_on_snow transitions, which are the most
# consequential. If max-cell melt in these peaks is clustered just above the
# 5mm threshold, the reclassification is likely threshold-sensitive noise.
# If it is broadly distributed, the reclassification is more physically credible.

changed_type_peaks <- peak_compare %>%
  filter(transition == "changed type") %>%
  left_join(
    maxcell_peaks %>%
      select(site_no, peak_id, total_rain_mm, total_melt_mm),
    by = c("site_no", "peak_id")
  ) %>%
  mutate(
    type_transition = paste0(flood_type_mean, " -> ", flood_type_maxcell)
  )

# Console summary of changed-type transitions
message("\n=== CHANGED-TYPE TRANSITION BREAKDOWN ===")
changed_type_peaks %>%
  count(type_transition) %>%
  arrange(desc(n)) %>%
  print()

# Focus on rainfall -> rain_on_snow, the most numerous and consequential
rain_to_ros <- changed_type_peaks %>%
  filter(flood_type_mean == "rainfall", flood_type_maxcell == "rain_on_snow")

message(sprintf("\nRainfall -> Rain-on-Snow transitions: %d peaks", nrow(rain_to_ros)))
message(sprintf("  Max-cell melt: median=%.1f mm, mean=%.1f mm",
                median(rain_to_ros$total_melt_mm, na.rm = TRUE),
                mean(rain_to_ros$total_melt_mm,   na.rm = TRUE)))
message(sprintf("  Max-cell melt 5-15mm (near-threshold): %d (%.1f%%)",
                sum(rain_to_ros$total_melt_mm >= 5  &
                    rain_to_ros$total_melt_mm <  15, na.rm = TRUE),
                100 * mean(rain_to_ros$total_melt_mm >= 5 &
                           rain_to_ros$total_melt_mm <  15, na.rm = TRUE)))
message(sprintf("  Max-cell melt >= 15mm (substantial):   %d (%.1f%%)",
                sum(rain_to_ros$total_melt_mm >= 15, na.rm = TRUE),
                100 * mean(rain_to_ros$total_melt_mm >= 15, na.rm = TRUE)))

# All changed-type transitions: scatter of max-cell rain vs melt,
# coloured by transition direction
fig07_data <- changed_type_peaks %>%
  filter(!is.na(total_rain_mm), !is.na(total_melt_mm)) %>%
  mutate(
    type_transition = factor(type_transition),
    # Flag near-threshold melt (5-15mm) for rainfall->ROS transitions
    near_threshold = flood_type_mean == "rainfall" &
                     flood_type_maxcell == "rain_on_snow" &
                     total_melt_mm >= 5 & total_melt_mm < 15
  )

fig07 <- ggplot(fig07_data,
                aes(x = total_rain_mm, y = total_melt_mm,
                    colour = type_transition)) +
  geom_point(alpha = 0.35, size = 1.4) +
  # Davenport threshold lines for reference
  geom_hline(yintercept =  5, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = 10, linetype = "solid",  colour = "grey40", linewidth = 0.5) +
  geom_vline(xintercept =  5, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 10, linetype = "solid",  colour = "grey40", linewidth = 0.5) +
  annotate("text", x = 195, y =  6, label = "melt 5mm",  hjust = 1,
           size = 3, colour = "grey40") +
  annotate("text", x = 195, y = 11, label = "melt 10mm", hjust = 1,
           size = 3, colour = "grey40") +
  scale_x_continuous(limits = c(0, 200)) +
  scale_y_continuous(limits = c(0, 200)) +
  labs(title    = "Changed-type peaks: max-cell rain vs melt (8-day cumulative)",
       subtitle = "Threshold lines at 5mm (dashed) and 10mm (solid)",
       x        = "Max-cell cumulative rain (mm)",
       y        = "Max-cell cumulative melt (mm)",
       colour   = "Transition") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "fig07_changed_type_melt_distribution.png"),
       fig07, width = 8, height = 6, dpi = 150)

# Fig 07b: histogram of max-cell melt for rainfall -> ROS transitions only,
# to make the near-threshold concentration visible
fig07b <- ggplot(rain_to_ros, aes(x = total_melt_mm)) +
  geom_histogram(binwidth = 2, fill = "#f46d43", colour = "white",
                 linewidth = 0.2) +
  geom_vline(xintercept =  5, linetype = "dashed", colour = "grey20",
             linewidth = 0.8) +
  geom_vline(xintercept = 15, linetype = "dotted", colour = "grey20",
             linewidth = 0.8) +
  annotate("text", x =  5.5, y = Inf, label = "5mm threshold",
           hjust = 0, vjust = 1.5, size = 3.2, colour = "grey20") +
  annotate("text", x = 15.5, y = Inf, label = "15mm (arbitrary\n'substantial' cutoff)",
           hjust = 0, vjust = 1.5, size = 3.2, colour = "grey20") +
  labs(title    = "Rainfall -> Rain-on-Snow transitions: max-cell melt distribution",
       subtitle = sprintf("n = %d peaks; near-threshold (5-15mm) = %.1f%%",
                          nrow(rain_to_ros),
                          100 * mean(rain_to_ros$total_melt_mm >= 5 &
                                     rain_to_ros$total_melt_mm < 15,
                                     na.rm = TRUE)),
       x = "Max-cell cumulative melt (mm, 8-day window)",
       y = "Number of peaks") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "fig07b_rainfall_to_ros_melt_histogram.png"),
       fig07b, width = 7, height = 5, dpi = 150)

# ============================================================================
# SAVE COMPARISON SUMMARY
# ============================================================================
comparison_summary <- list(
  peak_transitions  = transition_summary,
  gained_by_type    = peak_compare %>%
    filter(transition == "gained classification") %>%
    count(gained_type),
  lost_by_type      = peak_compare %>%
    filter(transition == "lost classification") %>%
    count(lost_type),
  changed_type_breakdown = changed_type_peaks %>%
    count(type_transition) %>%
    arrange(desc(n)),
  ffa_eligibility   = site_wide %>%
    group_by(flood_type) %>%
    summarise(
      eligible_mean    = sum(qualifies_mean),
      eligible_maxcell = sum(qualifies_maxcell),
      net_change       = eligible_maxcell - eligible_mean,
      .groups = "drop"
    ),
  peak_compare      = peak_compare,
  site_wide         = site_wide
)

saveRDS(comparison_summary, file.path(out_dir, "comparison_summary.rds"))
message(sprintf("\nSaved comparison summary to %s/comparison_summary.rds", out_dir))
message(sprintf("Figures saved to %s/", out_dir))