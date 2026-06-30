#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Three-Way Classification Comparison
# ============================================================================
# Visualises key findings from comparing basin-mean, max-cell, and
# cell-fraction classification approaches.
#
# Figures:
#   fig01  Three-way overall classification breakdown (stacked bar)
#   fig02  Cell fraction distributions by dominant type
#   fig03  Spatial map of dominant type by site (cell-fraction approach)
#          with elevation and drainage area panels
#   fig04  ROS story: max-cell melt vs area-weighted ROS fraction
#   fig05  FFA-eligible site counts across all three approaches
#
# INPUTS:
#   - _data/pot_peaks_classified_mean.rds
#   - _data/pot_peaks_classified_max_cell.rds
#   - _data/pot_peaks_cell_classified.rds
#   - _data/site_metadata.rds
#   - _data/station_spatial.rds  (for map)
# ============================================================================

library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

out_dir <- "_diagnostics/three_way_comparison"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MIN_TYPE_SAMPLE <- 20   # Yu et al. FFA eligibility threshold
FLOOD_TYPES     <- c("rainfall", "snowmelt", "rain_on_snow")

TYPE_COLORS <- c(
  "Rainfall"     = "#e41a1c",
  "Snowmelt"     = "#377eb8",
  "Rain-on-Snow" = "#4daf4a",
  "Mixed"        = "#ff7f00",
  "Unclassified" = "#cccccc"
)

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading data...")
mean_peaks    <- readRDS("_data/pot_peaks_classified_mean.rds")    %>%
  mutate(peak_date = as.Date(peak_date))
maxcell_peaks <- readRDS("_data/pot_peaks_classified_max_cell.rds") %>%
  mutate(peak_date = as.Date(peak_date))
cell_peaks    <- readRDS("_data/pot_peaks_cell_classified.rds")     %>%
  mutate(peak_date = as.Date(peak_date))
meta          <- readRDS("_data/site_metadata.rds")
stations      <- readRDS("_data/station_spatial.rds")

# Join max_cell melt values to cell_peaks for fig04
# (total_melt_mm in maxcell_peaks is the max-cell 8-day cumulative)
maxcell_melt <- maxcell_peaks %>%
  select(site_no, peak_id, melt_mm_maxcell = total_melt_mm)

cell_peaks <- cell_peaks %>%
  left_join(maxcell_melt, by = c("site_no", "peak_id"))

# ============================================================================
# FIG 01: Three-way overall classification breakdown
# ============================================================================
message("Fig 01: Three-way classification breakdown...")

fig01_data <- bind_rows(
  mean_peaks %>%
    mutate(mode = "Basin mean",
           type = flood_type,
           type = recode(type,
             rainfall     = "Rainfall",
             snowmelt     = "Snowmelt",
             rain_on_snow = "Rain-on-Snow",
             unclassified = "Unclassified"
           )),
  maxcell_peaks %>%
    mutate(mode = "Max cell",
           type = flood_type,
           type = recode(type,
             rainfall     = "Rainfall",
             snowmelt     = "Snowmelt",
             rain_on_snow = "Rain-on-Snow",
             unclassified = "Unclassified"
           )),
  cell_peaks %>%
    filter(complete_window) %>%
    mutate(mode = "Cell fraction",
           type = if_else(is.na(dominant_type), "Unclassified",
                   recode(dominant_type,
                     rainfall     = "Rainfall",
                     snowmelt     = "Snowmelt",
                     rain_on_snow = "Rain-on-Snow",
                     mixed        = "Mixed"
                   )))
) %>%
  count(mode, type) %>%
  group_by(mode) %>%
  mutate(
    pct  = 100 * n / sum(n),
    mode = factor(mode, levels = c("Basin mean", "Max cell", "Cell fraction")),
    type = factor(type, levels = c("Rainfall", "Snowmelt", "Rain-on-Snow",
                                   "Mixed", "Unclassified"))
  )

fig01 <- ggplot(fig01_data, aes(x = mode, y = pct, fill = type)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_stack(vjust = 0.5),
            size = 3.2, colour = "white", fontface = "bold") +
  scale_fill_manual(values = TYPE_COLORS) +
  labs(title    = "Peak classification breakdown across three approaches",
       x        = NULL,
       y        = "Percentage of peaks (%)",
       fill     = "Flood type") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "fig01_three_way_breakdown.png"),
       fig01, width = 7, height = 5, dpi = 150)

# ============================================================================
# FIG 02: Cell fraction distributions by dominant type
# ============================================================================
message("Fig 02: Cell fraction distributions...")

# For each peak, show the four fractions as a stacked bar sorted by
# dominant type. Summarise as violin/boxplot distributions per dominant type
# to show whether "mixed" peaks are truly ambiguous or lean toward one type.

fig02_data <- cell_peaks %>%
  filter(complete_window, !is.na(dominant_type)) %>%
  mutate(
    dominant_label = recode(dominant_type,
      rainfall     = "Rainfall",
      snowmelt     = "Snowmelt",
      rain_on_snow = "Rain-on-Snow",
      mixed        = "Mixed"
    ),
    dominant_label = factor(dominant_label,
                            levels = c("Rainfall", "Snowmelt",
                                       "Rain-on-Snow", "Mixed"))
  ) %>%
  select(dominant_label, frac_rainfall, frac_snowmelt, frac_ros, frac_mixed) %>%
  pivot_longer(cols = starts_with("frac_"),
               names_to  = "fraction_type",
               values_to = "fraction") %>%
  mutate(
    fraction_type = recode(fraction_type,
      frac_rainfall = "Rainfall",
      frac_snowmelt = "Snowmelt",
      frac_ros      = "Rain-on-Snow",
      frac_mixed    = "Mixed"
    ),
    fraction_type = factor(fraction_type,
                           levels = c("Rainfall", "Snowmelt",
                                      "Rain-on-Snow", "Mixed"))
  )

fig02 <- ggplot(fig02_data,
                aes(x = fraction_type, y = fraction, fill = fraction_type)) +
  geom_violin(scale = "width", alpha = 0.7, colour = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, outlier.alpha = 0.3,
               colour = "grey30", fill = "white") +
  facet_wrap(~ dominant_label, nrow = 1) +
  scale_fill_manual(values = TYPE_COLORS, guide = "none") +
  scale_y_continuous(labels = scales::percent) +
  labs(title    = "Cell fraction distributions by dominant type (cell-fraction approach)",
       subtitle = "Each panel shows the distribution of all four cell-type fractions\nfor peaks assigned to that dominant type",
       x        = "Cell type fraction",
       y        = "Area-weighted fraction") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "grey92"))

ggsave(file.path(out_dir, "fig02_fraction_distributions.png"),
       fig02, width = 10, height = 5, dpi = 150)

# ============================================================================
# FIG 03: Spatial map + elevation/area panels (cell-fraction approach)
# ============================================================================
message("Fig 03: Spatial patterns...")

# Site-level dominant type: plurality type across all peaks at the site
site_dominant <- cell_peaks %>%
  filter(complete_window, !is.na(dominant_type)) %>%
  count(site_no, dominant_type) %>%
  group_by(site_no) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(site_dominant_type = dominant_type)

# Join to spatial data
stations_df <- stations %>%
  st_drop_geometry() %>%
  left_join(site_dominant, by = "site_no") %>%
  filter(!is.na(site_dominant_type)) %>%
  mutate(
    type_label = recode(site_dominant_type,
      rainfall     = "Rainfall",
      snowmelt     = "Snowmelt",
      rain_on_snow = "Rain-on-Snow",
      mixed        = "Mixed"
    ),
    type_label = factor(type_label,
                        levels = c("Rainfall", "Snowmelt",
                                   "Rain-on-Snow", "Mixed"))
  )

# Natural Earth boundaries
us_states <- ne_states(country = "united states of america", returnclass = "sf")
co_border <- us_states %>% filter(name == "Colorado")
neighbors <- us_states %>%
  filter(name %in% c("Wyoming", "Nebraska", "Kansas",
                     "Oklahoma", "New Mexico", "Utah"))

# Map
fig03a <- ggplot() +
  geom_sf(data = neighbors,  fill = "grey96", colour = "grey70", linewidth = 0.3) +
  geom_sf(data = co_border,  fill = "grey92", colour = "grey40", linewidth = 0.5) +
  geom_point(data = stations_df,
             aes(x = lon, y = lat, colour = type_label),
             size = 2.2, alpha = 0.8) +
  scale_colour_manual(values = TYPE_COLORS) +
  coord_sf(xlim = c(-109.2, -102), ylim = c(36.8, 41.1), expand = FALSE) +
  labs(title  = "Dominant flood type by site (cell-fraction approach)",
       colour = "Dominant type", x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text = element_blank(), axis.ticks = element_blank())

# Elevation panel
fig03b <- ggplot(stations_df %>% filter(!is.na(elev_mean_m)),
                 aes(x = type_label, y = elev_mean_m, fill = type_label)) +
  geom_violin(scale = "width", alpha = 0.7, colour = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, fill = "white",
               colour = "grey30") +
  scale_fill_manual(values = TYPE_COLORS, guide = "none") +
  labs(x = NULL, y = "Mean basin elevation (m)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# Drainage area panel
fig03c <- ggplot(stations_df %>% filter(!is.na(drain_sqkm)),
                 aes(x = type_label, y = drain_sqkm, fill = type_label)) +
  geom_violin(scale = "width", alpha = 0.7, colour = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, fill = "white",
               colour = "grey30") +
  scale_fill_manual(values = TYPE_COLORS, guide = "none") +
  scale_y_log10(labels = scales::comma) +
  labs(x = NULL, y = "Drainage area (km², log scale)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# Combine with patchwork-style layout using cowplot if available,
# otherwise save separately
if (requireNamespace("cowplot", quietly = TRUE)) {
  library(cowplot)
  bottom_row <- plot_grid(fig03b, fig03c, nrow = 1, labels = c("B", "C"))
  fig03_combined <- plot_grid(fig03a, bottom_row, ncol = 1,
                              rel_heights = c(1.4, 1), labels = c("A", ""))
  ggsave(file.path(out_dir, "fig03_spatial_patterns.png"),
         fig03_combined, width = 10, height = 10, dpi = 150)
} else {
  ggsave(file.path(out_dir, "fig03a_map.png"),
         fig03a, width = 7, height = 5, dpi = 150)
  ggsave(file.path(out_dir, "fig03b_elevation.png"),
         fig03b, width = 5, height = 4, dpi = 150)
  ggsave(file.path(out_dir, "fig03c_area.png"),
         fig03c, width = 5, height = 4, dpi = 150)
}

# ============================================================================
# FIG 04: ROS story — max-cell melt vs area-weighted ROS fraction
# ============================================================================
message("Fig 04: ROS story...")

fig04_data <- cell_peaks %>%
  filter(complete_window, !is.na(melt_mm_maxcell)) %>%
  mutate(
    dominant_label = if_else(is.na(dominant_type), "Unclassified",
                      recode(dominant_type,
                        rainfall     = "Rainfall",
                        snowmelt     = "Snowmelt",
                        rain_on_snow = "Rain-on-Snow",
                        mixed        = "Mixed"
                      )),
    dominant_label = factor(dominant_label,
                            levels = c("Rainfall", "Snowmelt",
                                       "Rain-on-Snow", "Mixed"))
  )

fig04 <- ggplot(fig04_data,
                aes(x = melt_mm_maxcell, y = frac_ros,
                    colour = dominant_label)) +
  geom_point(alpha = 0.25, size = 1.2) +
  geom_hline(yintercept = 0.50, linetype = "dashed",
             colour = "grey30", linewidth = 0.7) +
  annotate("text", x = max(fig04_data$melt_mm_maxcell, na.rm = TRUE) * 0.95,
           y = 0.52, label = "50% threshold", hjust = 1,
           size = 3.2, colour = "grey30") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::percent) +
  scale_colour_manual(values = TYPE_COLORS) +
  labs(title    = "Max-cell melt vs area-weighted ROS cell fraction",
       subtitle = "Points at frac_ros > 50% are classified as Rain-on-Snow by cell-fraction approach\nMany high max-cell melt peaks have low ROS area fraction — spatially isolated melt",
       x        = "Max-cell cumulative melt (mm, 8-day window)",
       y        = "Area-weighted ROS fraction",
       colour   = "Dominant type") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "fig04_ros_maxcell_vs_fraction.png"),
       fig04, width = 8, height = 6, dpi = 150)

# ============================================================================
# FIG 05: FFA-eligible site counts across all three approaches
# ============================================================================
message("Fig 05: FFA-eligible site counts...")

# Helper: count FFA-eligible sites (>= 20 events) per flood type per approach
ffa_eligible <- function(peaks, type_col, mode_label) {
  peaks %>%
    filter(!!sym(type_col) %in% FLOOD_TYPES) %>%
    count(site_no, !!sym(type_col), name = "n_peaks") %>%
    filter(n_peaks >= MIN_TYPE_SAMPLE) %>%
    count(!!sym(type_col), name = "n_sites") %>%
    rename(flood_type = !!sym(type_col)) %>%
    mutate(mode = mode_label)
}

# Cell-fraction approach uses dominant_type with rain_on_snow label
cell_ffa <- cell_peaks %>%
  filter(complete_window, !is.na(dominant_type),
         dominant_type != "mixed") %>%
  mutate(flood_type = recode(dominant_type,
    rain_on_snow = "rain_on_snow",
    rainfall     = "rainfall",
    snowmelt     = "snowmelt"
  )) %>%
  count(site_no, flood_type) %>%
  filter(n >= MIN_TYPE_SAMPLE) %>%
  count(flood_type, name = "n_sites") %>%
  mutate(mode = "Cell fraction")

fig05_data <- bind_rows(
  ffa_eligible(mean_peaks,    "flood_type", "Basin mean"),
  ffa_eligible(maxcell_peaks, "flood_type", "Max cell"),
  cell_ffa
) %>%
  mutate(
    flood_type = recode(flood_type,
      rainfall     = "Rainfall",
      snowmelt     = "Snowmelt",
      rain_on_snow = "Rain-on-Snow"
    ),
    flood_type = factor(flood_type,
                        levels = c("Rainfall", "Snowmelt", "Rain-on-Snow")),
    mode = factor(mode, levels = c("Basin mean", "Max cell", "Cell fraction"))
  )

fig05 <- ggplot(fig05_data,
                aes(x = flood_type, y = n_sites, fill = mode)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = n_sites),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.8) +
  scale_fill_manual(values = c(
    "Basin mean"    = "#969696",
    "Max cell"      = "#2166ac",
    "Cell fraction" = "#4dac26"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title    = "FFA-eligible sites (>= 20 events per flood type)",
       subtitle = "Across all three classification approaches",
       x        = NULL,
       y        = "Number of sites",
       fill     = "Approach") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "fig05_ffa_eligible_three_way.png"),
       fig05, width = 7, height = 5, dpi = 150)

message(sprintf("\nAll figures saved to %s/", out_dir))