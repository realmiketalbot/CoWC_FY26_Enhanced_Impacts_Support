# Setup ----
library(tidyverse)
library(ggpattern)
library(tigris)
library(rnaturalearth)
library(rnaturalearthdata)

peaks_dt   <- readRDS("_data/pot_peaks_classified_depth_thresh.rds")
peaks_frac <- readRDS("_data/pot_peaks_classified_frac_thresh_3class.rds")
ffa_q      <- readRDS("_data/ffa_quantiles_gpd.rds")
ffa_p      <- readRDS("_data/ffa_params_gpd.rds")
meta       <- readRDS("_data/site_metadata.rds")

climate_trends <- readRDS("_data/climate_forcing_trends.rds")
fft_trends     <- readRDS("_data/temporal_trends_frac_thresh_3class.rds")
fig_dir <- "_documentation/report_figures"

station_spatial     <- readRDS("_data/station_spatial.rds")
ffa_q_lp3           <- readRDS("_data/ffa_quantiles_lp3.rds")
am_rain_frac_trends <- readRDS("_data/temporal_trends_am_frac_thresh_3class.rds")[["am_rain_frac_trends"]]

# Shared recode vector for curve_type -> display label (mirrors type_labels in app.R)
flood_type_recode <- c(
  rainfall     = "Rainfall",
  snowmelt     = "Snowmelt",
  rain_on_snow = "Rain-on-Snow"
)

# Section: Results ----

## Sub-section: Frequency of flood types across Colorado ----

### Figure: Map of study area (sites colored based on regulated/unregulated) ----

co_boundary_study <- tigris::states(cb = TRUE, class = "sf", progress_bar = FALSE) %>%
  dplyr::filter(STUSPS == "CO") %>%
  sf::st_transform(4326)

neighbor_states_study <- tigris::states(cb = TRUE, class = "sf", progress_bar = FALSE) %>%
  dplyr::filter(STUSPS %in% c("WY", "NE", "KS", "OK", "NM", "UT")) %>%
  sf::st_transform(4326)

study_area_data <- station_spatial %>%
  dplyr::left_join(
    sf::st_drop_geometry(meta) %>% dplyr::select(site_no, regulated),
    by = "site_no"
  ) %>%
  dplyr::mutate(
    reg_label = dplyr::if_else(regulated, "Regulated", "Unregulated"),
    reg_label = dplyr::if_else(is.na(reg_label), "Unknown", reg_label)
  )

ggplot() +
  ggplot2::geom_sf(data = neighbor_states_study, fill = "white", color = "grey60", linewidth = 0.3) +
  ggplot2::geom_sf(data = co_boundary_study, fill = "white", color = "grey20", linewidth = 0.5) +
  ggplot2::geom_sf(data = study_area_data, aes(color = reg_label, shape = reg_label), size = 2.2, alpha = 0.9) +
  ggplot2::scale_color_manual(
    values = c("Unregulated" = "#1E4D2B", "Regulated" = "#C8C372", "Unknown" = "grey50"),
    name   = "Regulation status"
  ) +
  ggplot2::scale_shape_manual(
    values = c("Unregulated" = 16, "Regulated" = 17, "Unknown" = 4),
    name   = "Regulation status"
  ) +
  ggplot2::coord_sf(
    xlim = sf::st_bbox(co_boundary_study)[c("xmin", "xmax")] + c(-0.5, 0.5),
    ylim = sf::st_bbox(co_boundary_study)[c("ymin", "ymax")] + c(-0.5, 0.5)
  ) +
  ggplot2::labs(
    title    = "Study Area and Streamgage Network",
    subtitle = "206 USGS gages across Colorado"
  ) +
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(legend.position = "right")

ggsave(file.path(fig_dir, "study_area_map.png"), width = 6.5, height = 4.5, units = "in")

### Figure: Histogram of rain fraction (population-level) ----
plot_data <- peaks_frac %>%
  dplyr::filter(data_sufficient, !is.na(rain_frac)) %>%
  dplyr::select(rain_frac)

ggplot(plot_data, aes(x = rain_frac)) +
  geom_histogram(binwidth = 0.04) +
  theme_bw() +
  theme(
    plot.subtitle       = element_text(size = 11, face = "italic"),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor.x  = element_blank()
  )

### Figure: Classification Composition (population-level) ----

# Build tidy labels for each method from the classified peaks directly.
# flood_type_raw captures the threshold-only classification before any
# sample-size demotion, allowing a like-for-like composition comparison.
label_dt <- peaks_dt %>%
  dplyr::mutate(
    method     = "Depth threshold",
    type_label = dplyr::case_when(
      as.character(flood_type_raw) == "rainfall"     ~ "Rainfall",
      as.character(flood_type_raw) == "snowmelt"     ~ "Snowmelt",
      as.character(flood_type_raw) == "rain_on_snow" ~ "Rain-on-Snow",
      TRUE                                           ~ "Unclassified"
    )
  ) %>%
  dplyr::select(method, type_label)

label_frac <- peaks_frac %>%
  dplyr::mutate(
    method     = "Rain fraction (3-class)",
    type_label = dplyr::case_when(
      as.character(flood_type_raw) == "rainfall"     ~ "Rainfall",
      as.character(flood_type_raw) == "snowmelt"     ~ "Snowmelt",
      as.character(flood_type_raw) == "rain_on_snow" ~ "Rain-on-Snow",
      TRUE                                           ~ "Unclassified"
    )
  ) %>%
  dplyr::select(method, type_label)

# 2-class: rainfall stays as-is, snowmelt + rain-on-snow collapse to
# Snow-Influenced (mirrors remap_2class() in Script 7 / SNOW_INFLUENCED_TYPES)
label_frac2 <- peaks_frac %>%
  dplyr::mutate(
    method     = "Rain fraction (2-class)",
    type_label = dplyr::case_when(
      as.character(flood_type_raw) == "rainfall"                      ~ "Rainfall",
      as.character(flood_type_raw) %in% c("snowmelt", "rain_on_snow") ~ "Snow-Influenced",
      TRUE                                                             ~ "Unclassified"
    )
  ) %>%
  dplyr::select(method, type_label)

comp_data <- dplyr::bind_rows(label_dt, label_frac, label_frac2) %>%
  dplyr::mutate(
    type_label = factor(
      type_label,
      levels = c("Unclassified", "Rainfall", "Rain-on-Snow", "Snow-Influenced", "Snowmelt")
    )
  ) %>%
  dplyr::count(method, type_label) %>%
  dplyr::group_by(method) %>%
  dplyr::mutate(
    pct    = n / sum(n),
    method = factor(
      method,
      levels = c("Depth threshold", "Rain fraction (3-class)", "Rain fraction (2-class)")
    )
  ) %>%
  dplyr::ungroup()

# Solid fill colors -- "Snow-Influenced" gets the same blue as Snowmelt as its
# background; the stripe foreground color is set separately below via
# pattern_fill, since that's a fixed parameter rather than a mapped aesthetic.
comp_colors <- c(
  "Rainfall"        = "#b2182b",
  "Snowmelt"        = "#2166ac",
  "Rain-on-Snow"    = "#9970ab",
  "Snow-Influenced" = "#2166ac",
  "Unclassified"    = "#d9d9d9"
)

# Only "Snow-Influenced" gets a pattern; everything else is a solid fill
comp_patterns <- c(
  "Rainfall"        = "none",
  "Snowmelt"        = "none",
  "Rain-on-Snow"    = "none",
  "Snow-Influenced" = "stripe",
  "Unclassified"    = "none"
)

ggplot(comp_data, aes(x = method, y = pct, fill = type_label, pattern = type_label)) +
  ggpattern::geom_col_pattern(
    position                 = "stack",
    linewidth                = 0,
    color                    = NA,
    pattern_fill             = "#9970ab",  # purple stripe foreground
    pattern_colour           = NA,
    pattern_density          = 0.5,
    pattern_spacing          = 0.04,
    pattern_angle            = 45,
    pattern_key_scale_factor = 0.4
  ) +
  scale_fill_manual(values = comp_colors, name = "Flood type") +
  scale_pattern_manual(values = comp_patterns, name = "Flood type") +
  scale_y_continuous(labels = scales::percent, name = "Fraction of all POT peaks", expand = c(0, 0)) +
  scale_x_discrete(name = "Classification method", expand = c(0, 0)) +
  labs(
    title    = "Flood Type Composition Across Classification Methods",
    subtitle = "Share of all peaks-over-threshold events attributed to each flood-generating process"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position     = "right",
    plot.subtitle       = element_text(size = 11, face = "italic"),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor.x  = element_blank()
  )

ggsave(file.path(fig_dir, "classification_composition_population_level.png"), width = 6.5, height = 3.5, units = "in")

### Figure: Frequency of flood types across Colorado ----
plot_data <- peaks_frac %>%
  dplyr::select(site_no, flood_type_raw) %>%
  dplyr::filter(flood_type_raw != "unclassified") %>%
  dplyr::group_by(site_no) %>%
  dplyr::mutate(total = n()) %>%
  dplyr::group_by(site_no, flood_type_raw) %>%
  dplyr::summarize(
    fraction = n() / mean(total),
    .groups = "drop"
  )

station_order <- plot_data %>%
  tidyr::pivot_wider(
    names_from = flood_type_raw, 
    values_from = fraction
  ) %>%
  dplyr::mutate(raininess = rainfall / (rainfall + rain_on_snow + snowmelt)) %>%
  dplyr::arrange(raininess) %>%
  dplyr::pull(site_no)

plot_data <- plot_data %>% 
  dplyr::mutate(site_no = factor(site_no, levels = station_order)) %>%
  dplyr::mutate(
    type_label = dplyr::case_when(
      as.character(flood_type_raw) == "rainfall"     ~ "Rainfall",
      as.character(flood_type_raw) == "snowmelt"     ~ "Snowmelt",
      as.character(flood_type_raw) == "rain_on_snow" ~ "Rain-on-Snow",
      TRUE                                           ~ "Unclassified"
    ),
    type_label = factor(type_label, levels = c("Unclassified", "Rainfall", "Rain-on-Snow", "Snowmelt")) 
  )

comp_colors <- c(
  "Rainfall"     = "#b2182b",
  "Snowmelt"     = "#2166ac",
  "Rain-on-Snow" = "#9970ab",
  "Unclassified" = "#d9d9d9"
)

ggplot(plot_data, aes(x = site_no, y = fraction, fill = type_label)) +
  geom_histogram(stat = "identity", position = "stack", color = "white", linewidth = 0) +
  scale_fill_manual(values = comp_colors, name = "Flood type") +
  scale_y_continuous(labels = scales::percent, name = "Fraction of classified POT peaks", expand = c(0, 0)) +
  scale_x_discrete(name = "Site (ordered by fraction of rainfall peaks)", expand = c(0, 0)) +
  labs(
    title    = "Site-Level Variation in Flood Type Composition",
    subtitle = "Based on the rain fraction (3-class) method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    plot.subtitle = element_text(size = 11, face = "italic"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
  
ggsave(file.path(fig_dir, "classification_composition_site_level.png"), width = 6.5, height = 3.5, units = "in")

### Figure: Fraction of rainfall peaks across Colorado ----

site_raininess <- peaks_frac %>%
  dplyr::filter(flood_type_raw != "unclassified") %>%
  dplyr::count(site_no, flood_type_raw) %>%
  tidyr::pivot_wider(
    names_from  = flood_type_raw,
    values_from = n,
    values_fill = 0
  ) %>%
  dplyr::mutate(raininess = rainfall / (rainfall + rain_on_snow + snowmelt)) %>%
  dplyr::select(site_no, raininess)

map_data <- site_raininess %>%
  dplyr::inner_join(meta, by = "site_no") %>%
  sf::st_as_sf(coords = c("lon_deg", "lat_deg"), crs = 4326)

co_boundary <- tigris::states(cb = TRUE, class = "sf", progress_bar = FALSE) %>%
  dplyr::filter(STUSPS == "CO") %>%
  sf::st_transform(4326)

neighbor_states <- tigris::states(cb = TRUE, class = "sf", progress_bar = FALSE) %>%
  dplyr::filter(STUSPS %in% c("WY", "NE", "KS", "OK", "NM", "UT")) %>%
  sf::st_transform(4326)

ggplot() +
  ggplot2::geom_sf(data = neighbor_states, fill = "white", color = "grey60", linewidth = 0.3) +
  ggplot2::geom_sf(data = co_boundary, fill = "white", color = "grey20", linewidth = 0.5) +
  ggplot2::geom_sf(data = map_data, aes(color = raininess), size = 2.2, alpha = 0.9) +
  ggplot2::scale_color_gradientn(
    colours = raininess_pal,
    limits  = c(0, 1),
    name    = "Fraction of \nrainfall peaks"
  ) +
  ggplot2::coord_sf(
    xlim = sf::st_bbox(co_boundary)[c("xmin", "xmax")] + c(-0.5, 0.5),
    ylim = sf::st_bbox(co_boundary)[c("ymin", "ymax")] + c(-0.5, 0.5)
  ) +
  ggplot2::labs(
    title    = "Spatial Patterns of Flood Classifications",
    subtitle = "Based on the rain fraction (3-class) method"
  ) +
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(
    legend.position = "right",
    plot.subtitle       = element_text(size = 11, face = "italic"),
  )

ggsave(file.path(fig_dir, "fraction_of_rainfall_peaks_map.png"), width = 6.5, height = 5, units = "in")

### Figure: fraction of rainfall peaks vs elevation ----

# site_raininess (computed above for the raininess map) is one row per site:
# fraction of that site's classified POT peaks that are rainfall-driven.
elev_raininess <- site_raininess %>%
  dplyr::inner_join(
    sf::st_drop_geometry(station_spatial) %>% dplyr::select(site_no, elev_mean_m),
    by = "site_no"
  )

ggplot(elev_raininess, aes(x = elev_mean_m, y = raininess)) +
  geom_point() +
  geom_smooth(method = "loess", se = F) +
  scale_y_continuous(labels = scales::percent, name = "Fraction of rainfall peaks", limits = c(0, 1)) +
  scale_x_continuous(name = "Basin mean elevation (m)") +
  labs(
    title    = "Fraction of Rainfall Peaks vs Basin Elevation",
    subtitle = "Based on the rain fraction (3-class) method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.subtitle       = element_text(size = 11, face = "italic"),
  )

ggsave(file.path(fig_dir, "fraction_of_rainfall_peaks_vs_elevation_site_level.png"), width = 6.5, height = 4.5, units = "in")

### Table: number of watersheds with sufficient data by flood type ----

sufficient_dt <- peaks_dt %>%
  dplyr::filter(flood_type != "unclassified") %>%
  dplyr::distinct(site_no, flood_type) %>%
  dplyr::count(flood_type, name = "n_sites") %>%
  dplyr::mutate(method = "Depth threshold")

sufficient_frac <- peaks_frac %>%
  dplyr::filter(flood_type != "unclassified") %>%
  dplyr::distinct(site_no, flood_type) %>%
  dplyr::count(flood_type, name = "n_sites") %>%
  dplyr::mutate(method = "Rain fraction (3-class)")

sufficient_table <- dplyr::bind_rows(sufficient_dt, sufficient_frac) %>%
  dplyr::mutate(
    flood_type = dplyr::case_when(
      flood_type == "rainfall"     ~ "Rainfall",
      flood_type == "snowmelt"     ~ "Snowmelt",
      flood_type == "rain_on_snow" ~ "Rain-on-Snow",
      TRUE                         ~ flood_type
    )
  ) %>%
  tidyr::pivot_wider(names_from = method, values_from = n_sites)

knitr::kable(
  sufficient_table,
  col.names = c("Flood type", "Depth threshold", "Rain fraction (3-class)"),
  caption   = "Number of watersheds (of 206) with sufficient peaks (\u2265 20) to support a per-type distribution fit."
)

### Figure: scatterplot of fraction of rainfall peaks vs elevation

# This looks like a duplicate of "Figure: fraction of rainfall peaks vs
# elevation" above (raininess_vs_elevation_site_level.png) -- same variables,
# same subsection. I haven't regenerated it here; if you meant something
# different (e.g. faceted by classification method, or scoped to a subset of
# sites), let me know and I'll build that version instead.

## Sub-section: Mixed vs. single sample distribution FFA results ----

### Figure: map of dominant flood type for the 100-yr return period, GPD, Rain fraction (3-class) ----

# Mirrors the "dominant" symbology logic in app.R: at each site, the
# per-type curve (excluding mixture/single_sample) with the highest Q100
# is the "dominant" flood type driving that return period.
dominant_100yr <- ffa_q %>%
  dplyr::filter(
    method        == "frac_thresh_3class",
    return_period == 100,
    !curve_type   %in% c("mixture", "single_sample")
  ) %>%
  dplyr::filter(is.finite(discharge_cfs), discharge_cfs > 0) %>%
  dplyr::group_by(site_no) %>%
  dplyr::slice_max(discharge_cfs, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(dominant_label = dplyr::recode(as.character(curve_type), !!!flood_type_recode)) %>%
  dplyr::select(site_no, dominant_label)

dominant_map_data <- station_spatial %>%
  dplyr::inner_join(dominant_100yr, by = "site_no")

ggplot() +
  ggplot2::geom_sf(data = neighbor_states, fill = "white", color = "grey60", linewidth = 0.3) +
  ggplot2::geom_sf(data = co_boundary, fill = "white", color = "grey20", linewidth = 0.5) +
  ggplot2::geom_sf(data = dominant_map_data, aes(color = dominant_label), size = 2.2, alpha = 0.9) +
  ggplot2::scale_color_manual(values = comp_colors, name = "Dominant flood type\n(Q100, GPD)") +
  ggplot2::coord_sf(
    xlim = sf::st_bbox(co_boundary)[c("xmin", "xmax")] + c(-0.5, 0.5),
    ylim = sf::st_bbox(co_boundary)[c("ymin", "ymax")] + c(-0.5, 0.5)
  ) +
  ggplot2::labs(
    title    = "Dominant Flood Type at the 100-Year Return Period",
    subtitle = "Based on the rain fraction (3-class) method and the Generalized Pareto distribution"
  ) +
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(
    legend.position = "right",
    plot.subtitle = element_text(size = 11, face = "italic")
  )

ggsave(file.path(fig_dir, "dominant_flood_type_100yr_gpd.png"), width = 6.5, height = 4.5, units = "in")

# Shared helper for the next two plots: build a site x ratio table (mixture / single_sample) at a
# given return period, for a given ffa_quantiles family object.
compute_ratio_data <- function(qdata, rp = 100, method_name = "frac_thresh_3class") {
  qdata %>%
    dplyr::filter(
      method        == method_name,
      return_period == rp,
      curve_type    %in% c("mixture", "single_sample")
    ) %>%
    dplyr::select(site_no, curve_type, discharge_cfs) %>%
    tidyr::pivot_wider(names_from = curve_type, values_from = discharge_cfs) %>%
    dplyr::filter(!is.na(mixture), !is.na(single_sample), mixture > 0, single_sample > 0) %>%
    dplyr::mutate(ratio = mixture / single_sample, log_ratio = log(ratio))
}

# Shared helper: render the ratio map for one family and save it.
plot_ratio_map <- function(ratio_df, family_label_txt, out_name) {
  map_df   <- station_spatial %>% dplyr::inner_join(ratio_df, by = "site_no")
  max_abs  <- max(abs(map_df$log_ratio), na.rm = TRUE)

  p <- ggplot() +
    ggplot2::geom_sf(data = neighbor_states, fill = "white", color = "grey60", linewidth = 0.3) +
    ggplot2::geom_sf(data = co_boundary, fill = "white", color = "grey20", linewidth = 0.5) +
    ggplot2::geom_sf(data = map_df, aes(color = log_ratio), size = 2.2, alpha = 0.9) +
    ggplot2::scale_color_gradient2(
      low      = "#2166ac", mid = "grey90", high = "#b2182b", midpoint = 0,
      limits   = c(-max_abs, max_abs),
      labels   = function(x) round(exp(x), 2),
      name     = "Q100 ratio\n(mixture / single sample)"
    ) +
    ggplot2::coord_sf(
      xlim = sf::st_bbox(co_boundary)[c("xmin", "xmax")] + c(-0.5, 0.5),
      ylim = sf::st_bbox(co_boundary)[c("ymin", "ymax")] + c(-0.5, 0.5)
    ) +
    ggplot2::labs(
      title    = paste0("Mixture vs. Single-Sample Q100 Ratio (", family_label_txt, ")"),
      subtitle = "Based on the rain fraction (3-class) method"
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      plot.subtitle = element_text(size = 11, face = "italic")
    )

  ggsave(file.path(fig_dir, out_name), plot = p, width = 6.5, height = 4.5, units = "in")
  p
}

ratio_gpd_100 <- compute_ratio_data(ffa_q,     rp = 100)
ratio_lp3_100 <- compute_ratio_data(ffa_q_lp3, rp = 100)

### Figure: map of flood ratio of 100-year return period, GPD, Rain fraction (3-class) (mixed over single sample) ----
plot_ratio_map(ratio_gpd_100, "GPD", "ratio_map_100yr_gpd.png")

### Figure: map of flood ratio of 100-year return period, LP3, Rain fraction (3-class) (mixed over single sample) ----
plot_ratio_map(ratio_lp3_100, "LP3", "ratio_map_100yr_lp3.png")

### Figure: scatterplot of mixed vs single sample for the 100-year, GPD ----

# Reuses ratio_gpd_100 computed just above; joins back to `regulated` for
# point styling, matching the convention used in app.R's ffa_q_plot.
scatter_100_gpd <- ratio_gpd_100 %>%
  dplyr::left_join(
    sf::st_drop_geometry(meta) %>% dplyr::select(site_no, regulated),
    by = "site_no"
  ) %>%
  dplyr::mutate(reg_label = dplyr::if_else(regulated, "Regulated", "Unregulated"))

ggplot(scatter_100_gpd, aes(x = single_sample, y = mixture, color = reg_label, shape = reg_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_point(alpha = 0.7, size = 2.2) +
  scale_x_log10(labels = scales::label_comma(), name = "Single-sample Q100 (cfs)") +
  scale_y_log10(labels = scales::label_comma(), name = "Mixture Q100 (cfs)") +
  scale_color_manual(values = c("Unregulated" = "#1E4D2B", "Regulated" = "#C8C372"), name = "Regulation") +
  scale_shape_manual(values = c("Unregulated" = 16, "Regulated" = 17), name = "Regulation") +
  labs(
    title    = "Mixture vs. Single-Sample Discharge Estimates at Q100",
    subtitle = "Based on the rain fraction (3-class) method and the Generalized Pareto distribution"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom", 
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(size = 11, face = "italic")
  )

ggsave(file.path(fig_dir, "scatter_mixture_vs_single_100yr_gpd.png"), width = 6.5, height = 5, units = "in")

### Figure: box plot of ratios comparing across all return periods for GPD and Rain fraction (3-class) (mixed over single sample) ----

ratio_all_rp_gpd <- ffa_q %>%
  dplyr::filter(method == "frac_thresh_3class", curve_type %in% c("mixture", "single_sample")) %>%
  dplyr::select(site_no, return_period, curve_type, discharge_cfs) %>%
  tidyr::pivot_wider(names_from = curve_type, values_from = discharge_cfs) %>%
  dplyr::filter(!is.na(mixture), !is.na(single_sample), mixture > 0, single_sample > 0) %>%
  dplyr::mutate(
    ratio         = mixture / single_sample,
    return_period = factor(return_period, levels = sort(unique(return_period)))
  )

ggplot(ratio_all_rp_gpd, aes(x = return_period, y = ratio)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6) +
  scale_y_log10(labels = scales::label_comma(), name = "Ratio (mixture / single sample)") +
  scale_x_discrete(name = "Return period (years)") +
  labs(
    title    = "Mixture-to-Single-Sample Ratio Across Return Periods",
    subtitle = "Based on the rain fraction (3-class) method and the Generalized Pareto distribution"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(size = 11, face = "italic")
  )

ggsave(file.path(fig_dir, "ratio_boxplot_all_rp_gpd.png"), width = 6.5, height = 4.5, units = "in")

### Figure: column plot of binary (above/below 1.0 ratio) fraction for all return periods for GPD and Rain fraction (3-class) (mixed over single sample) ----

# Reuses ratio_all_rp_gpd computed just above.
ratio_binary_counts <- ratio_all_rp_gpd %>%
  dplyr::mutate(
    above_below = dplyr::if_else(ratio >= 1, "Mixture \u2265 Single Sample", "Mixture < Single Sample")
  ) %>%
  dplyr::count(return_period, above_below) %>%
  dplyr::group_by(return_period) %>%
  dplyr::mutate(
    frac        = n / sum(n),
    signed_frac = dplyr::if_else(above_below == "Mixture \u2265 Single Sample", frac, -frac)
  ) %>%
  dplyr::ungroup()

ggplot(ratio_binary_counts, aes(x = return_period, y = signed_frac, fill = above_below)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
  scale_fill_manual(
    values = c("Mixture \u2265 Single Sample" = "#b2182b", "Mixture < Single Sample" = "#2166ac"),
    name   = NULL
  ) +
  scale_y_continuous(
    name   = "Fraction of sites",
    labels = function(x) scales::percent(abs(x)),
    limits = c(-1, 1)
  ) +
  scale_x_discrete(name = "Return period (years)") +
  labs(
    title    = "Sites Where Mixture Exceeds Single-Sample Estimate",
    subtitle = "Based on the rain fraction (3-class) method and the Generalized Pareto distribution"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position     = "bottom",
    panel.grid.major.x  = element_blank(),
    panel.grid.minor.x  = element_blank(),
    plot.subtitle       = element_text(size = 11, face = "italic")
  )

ggsave(file.path(fig_dir, "ratio_binary_count_all_rp_gpd.png"), width = 6.5, height = 4.5, units = "in")

## Sub-section: Trends in flood drivers ----
 
### Figure: map of annual maxima rain fraction ----

# Interpreting this placeholder as the Mann-Kendall trend in AM rain
# fraction (Sen slope), matching the "Annual Maxima Rain Fraction" layer on
# the app's Trends tab -- as opposed to a single-year snapshot, which
# wouldn't have an obvious "map" to show. If you actually wanted a snapshot
# (e.g. mean AM rain fraction over the record, or for one water year), let
# me know and I'll swap this out.
rf_trend_map_data <- station_spatial %>%
  dplyr::inner_join(
    am_rain_frac_trends %>% dplyr::select(site_no, sen_slope, significant),
    by = "site_no"
  ) %>%
  dplyr::filter(!is.na(sen_slope))

max_abs_slope <- max(abs(rf_trend_map_data$sen_slope), na.rm = TRUE)

ggplot() +
  ggplot2::geom_sf(data = neighbor_states, fill = "white", color = "grey60", linewidth = 0.3) +
  ggplot2::geom_sf(data = co_boundary, fill = "white", color = "grey20", linewidth = 0.5) +
  ggplot2::geom_sf(data = rf_trend_map_data, aes(color = sen_slope, shape = significant), size = 2.2, alpha = 0.9) +
  ggplot2::scale_color_gradient2(
    low = "#2166ac", mid = "grey90", high = "#b2182b", midpoint = 0,
    limits = c(-max_abs_slope, max_abs_slope),
    name   = "Sen slope\n(rain fraction/yr)"
  ) +
  ggplot2::scale_shape_manual(
    values = c(`TRUE` = 17, `FALSE` = 16),
    labels = c(`TRUE` = "p < 0.05", `FALSE` = "Not significant"),
    name   = "Mann-Kendall test"
  ) +
  ggplot2::coord_sf(
    xlim = sf::st_bbox(co_boundary)[c("xmin", "xmax")] + c(-0.5, 0.5),
    ylim = sf::st_bbox(co_boundary)[c("ymin", "ymax")] + c(-0.5, 0.5)
  ) +
  ggplot2::labs(
    title    = "Trend in Annual Maxima Rain Fraction",
    subtitle = "Direction and magnitude of long-term change in rainfall's share of the annual peak flow, by site"
  ) +
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(
    legend.position = "right",
    plot.subtitle = element_text(size = 11, face = "italic")
  )

ggsave(file.path(fig_dir, "am_rain_fraction_trend_map.png"), width = 6.5, height = 4.5, units = "in")

### Figure: scatterplot of rain fraction vs elevation (all PoT events) ----

# Peak-level (not site-aggregated): every classified POT event's individual
# rain_frac value against its site's elevation, so this shows the underlying
# classification variable directly rather than a per-site summary.
peak_elev_data <- peaks_frac %>%
  dplyr::filter(!is.na(rain_frac)) %>%
  dplyr::select(site_no, rain_frac) %>%
  dplyr::inner_join(
    sf::st_drop_geometry(station_spatial) %>% dplyr::select(site_no, elev_mean_m),
    by = "site_no"
  )

ggplot(peak_elev_data, aes(x = elev_mean_m, y = rain_frac)) +
  geom_point(alpha = 0.15, size = 1) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8) +
  scale_y_continuous(name = "Rain fraction (individual POT events)", limits = c(0, 1)) +
  scale_x_continuous(name = "Basin mean elevation (m)") +
  labs(
    title    = "Rain Fraction of Individual Flood Events by Elevation",
    subtitle = "Every classified peak-over-threshold event across all sites; line shows the overall trend"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.subtitle       = element_text(size = 11, face = "italic")
  )

ggsave(file.path(fig_dir, "rain_fraction_vs_elevation_all_events.png"), width = 6.5, height = 4.5, units = "in")

### Figure: trend in annual maxima streamflow ----
 
# Mirrors the rain fraction trend map above, using the peak-flow layer that's
# already in the app's Trends tab (input$trends_map_layer == "peak_flow").
# Same Sen slope / significance encoding, just a different variable and units
# (cfs/yr instead of fraction/yr) -- kept as a separate figure rather than a
# facet since the color scale domain differs by orders of magnitude across
# variables and would wash out one or the other on a shared scale.
pf_trend_map_data <- station_spatial %>%
  dplyr::inner_join(
    am_peak_flow_trends %>% dplyr::select(site_no, sen_slope, significant),
    by = "site_no"
  ) %>%
  dplyr::filter(!is.na(sen_slope))
 
max_abs_slope_pf <- max(abs(pf_trend_map_data$sen_slope), na.rm = TRUE)
 
ggplot() +
  ggplot2::geom_sf(data = neighbor_states, fill = "white", color = "grey60", linewidth = 0.3) +
  ggplot2::geom_sf(data = co_boundary, fill = "white", color = "grey20", linewidth = 0.5) +
  ggplot2::geom_sf(data = pf_trend_map_data, aes(color = sen_slope, shape = significant), size = 2.2, alpha = 0.9) +
  ggplot2::scale_color_gradient2(
    low = "#2166ac", mid = "grey90", high = "#b2182b", midpoint = 0,
    limits = c(-max_abs_slope_pf, max_abs_slope_pf),
    labels = scales::label_comma(),
    name   = "Sen slope\n(cfs/yr)"
  ) +
  ggplot2::scale_shape_manual(
    values = c(`TRUE` = 17, `FALSE` = 16),
    labels = c(`TRUE` = "p < 0.05", `FALSE` = "Not significant"),
    name   = "Mann-Kendall test"
  ) +
  ggplot2::coord_sf(
    xlim = sf::st_bbox(co_boundary)[c("xmin", "xmax")] + c(-0.5, 0.5),
    ylim = sf::st_bbox(co_boundary)[c("ymin", "ymax")] + c(-0.5, 0.5)
  ) +
  ggplot2::labs(
    title    = "Trend in Annual Maxima Streamflow",
    subtitle = "Direction and magnitude of long-term change in the annual peak flow, by gauge"
  ) +
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(
    legend.position = "right",
    plot.subtitle       = element_text(size = 11, face = "italic")
  )
 
ggsave(file.path(fig_dir, "am_peak_flow_trend_map.png"), width = 6.5, height = 4.5, units = "in")
