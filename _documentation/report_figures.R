# Setup ----
library(tidyverse)
library(ggpattern)
library(legendry)
library(patchwork)
library(colorspace)
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

co_counties_study <- tigris::counties(state = "CO", cb = TRUE, class = "sf", progress_bar = FALSE) %>%
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
  ggplot2::geom_sf(data = co_boundary_study, fill = "white", color = NA) +
  ggplot2::geom_sf(data = co_counties_study, fill = NA, color = "grey85", linewidth = 0.2) +
  ggplot2::geom_sf(data = co_boundary_study, fill = NA, color = "grey20", linewidth = 0.5) +
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
    xlim = sf::st_bbox(co_boundary_study)[c("xmin", "xmax")] + c(0, 0),
    ylim = sf::st_bbox(co_boundary_study)[c("ymin", "ymax")] + c(0, 0)
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
  geom_histogram(binwidth = 0.05) +
  labs(
    title = "Histogram of Rain Fraction",
    subtitle = "Population-level histogram across all sites; bin width of 0.05",
    x = "Rain Fraction", 
    y = "# of peaks"
  ) +
  theme_bw() +
  theme(
    plot.subtitle       = element_text(size = 11, face = "italic"),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor.x  = element_blank()
  )

ggsave(file.path(fig_dir, "rain_fraction_histogram.png"), width = 6.5, height = 3.5, units = "in")

### Figure: Classification Composition (population-level) ----

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

# Collapsing mirrors remap_2class() in Script 7 / SNOW_INFLUENCED_TYPES.
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

# The stripe foreground for "Snow-Influenced" is set separately below via
# pattern_fill, since that's a fixed ggpattern parameter, not a mapped aesthetic.
comp_colors <- c(
  "Rainfall"        = "#b2182b",
  "Snowmelt"        = "#2166ac",
  "Rain-on-Snow"    = "#9970ab",
  "Snow-Influenced" = "#2166ac",
  "Unclassified"    = "#d9d9d9"
)

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

co_counties <- tigris::counties(state = "CO", cb = TRUE, class = "sf", progress_bar = FALSE) %>%
  sf::st_transform(4326)

raininess_pal <- c("#2166ac", "#9970ab", "#b2182b")

ggplot() +
  ggplot2::geom_sf(data = co_boundary, fill = "white", color = NA) +
  ggplot2::geom_sf(data = co_counties, fill = NA, color = "grey85", linewidth = 0.2) +
  ggplot2::geom_sf(data = co_boundary, fill = NA, color = "grey20", linewidth = 0.5) +
  ggplot2::geom_sf(data = map_data, aes(color = raininess), size = 2.2, alpha = 0.9) +
  ggplot2::scale_color_gradientn(
    colours = raininess_pal,
    limits  = c(0, 1),
    name    = "Fraction of \nrainfall peaks"
  ) +
  ggplot2::coord_sf(
    xlim = sf::st_bbox(co_boundary)[c("xmin", "xmax")] + c(0, 0),
    ylim = sf::st_bbox(co_boundary)[c("ymin", "ymax")] + c(0, 0)
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
  geom_point(aes(color = raininess), size = 2.2, alpha = 0.9) +
  geom_smooth(method = "loess", se = FALSE, color = "grey20", linewidth = 0.8) +
  scale_color_gradientn(
    colours = raininess_pal,
    limits  = c(0, 1),
    name    = "Fraction of \nrainfall peaks"
  ) +
  scale_y_continuous(labels = scales::percent, name = "Fraction of rainfall peaks", limits = c(0, 1)) +
  scale_x_continuous(name = "Basin mean elevation (m)") +
  labs(
    title    = "Fraction of Rainfall Peaks vs Basin Elevation",
    subtitle = "Based on the rain fraction (3-class) method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(size = 11, face = "italic"),
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

### Table: number of watersheds by number of flood type curves available ----

site_type_counts <- function(peaks) {
  peaks %>%
    dplyr::filter(flood_type != "unclassified") %>%
    dplyr::distinct(site_no, flood_type) %>%
    dplyr::count(site_no, name = "n_types")
}

bucket_labels <- c(
  "0" = "None (single-sample only)",
  "1" = "One flood type curve",
  "2" = "Two flood type curves",
  "3" = "Three flood type curves"
)

curve_count_dt <- tibble::tibble(site_no = unique(peaks_dt$site_no)) %>%
  dplyr::left_join(site_type_counts(peaks_dt), by = "site_no") %>%
  dplyr::mutate(n_types = tidyr::replace_na(n_types, 0L)) %>%
  dplyr::count(n_types, name = "n_sites") %>%
  dplyr::mutate(method = "Depth threshold")

curve_count_frac <- tibble::tibble(site_no = unique(peaks_frac$site_no)) %>%
  dplyr::left_join(site_type_counts(peaks_frac), by = "site_no") %>%
  dplyr::mutate(n_types = tidyr::replace_na(n_types, 0L)) %>%
  dplyr::count(n_types, name = "n_sites") %>%
  dplyr::mutate(method = "Rain fraction (3-class)")

curve_count_table <- dplyr::bind_rows(curve_count_dt, curve_count_frac) %>%
  dplyr::mutate(
    curves_available = dplyr::recode(as.character(n_types), !!!bucket_labels),
    curves_available = factor(curves_available, levels = unname(bucket_labels))
  ) %>%
  tidyr::complete(curves_available, method, fill = list(n_sites = 0L)) %>%
  dplyr::select(curves_available, method, n_sites) %>%
  tidyr::pivot_wider(names_from = method, values_from = n_sites)

knitr::kable(
  curve_count_table,
  col.names = c("Curves available", "Depth threshold", "Rain fraction (3-class)"),
  caption   = "Number of watersheds (of 206) by number of flood-type curves available for the mixture distribution."
)

## Sub-section: Mixed vs. single sample distribution FFA results ----

### Figure: map of dominant flood type for the 100-yr return period, GPD, Rain fraction (3-class) ----

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
  ggplot2::geom_sf(data = co_boundary, fill = "white", color = NA) +
  ggplot2::geom_sf(data = co_counties, fill = NA, color = "grey85", linewidth = 0.2) +
  ggplot2::geom_sf(data = co_boundary, fill = NA, color = "grey20", linewidth = 0.5) +
  ggplot2::geom_sf(data = dominant_map_data, aes(color = dominant_label), size = 2.2, alpha = 0.9) +
  ggplot2::scale_color_manual(values = comp_colors, name = "Dominant flood type\n(Q100, GPD)") +
  ggplot2::coord_sf(
    xlim = sf::st_bbox(co_boundary)[c("xmin", "xmax")] + c(0, 0),
    ylim = sf::st_bbox(co_boundary)[c("ymin", "ymax")] + c(0, 0)
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
    ggplot2::geom_sf(data = co_boundary, fill = "white", color = NA) +
    ggplot2::geom_sf(data = co_counties, fill = NA, color = "grey85", linewidth = 0.2) +
    ggplot2::geom_sf(data = co_boundary, fill = NA, color = "grey20", linewidth = 0.5) +
    ggplot2::geom_sf(data = map_df, aes(color = log_ratio), size = 2.2, alpha = 0.9) +
    ggplot2::scale_color_gradient2(
      low      = "#2166ac", mid = "grey90", high = "#b2182b", midpoint = 0,
      limits   = c(-max_abs, max_abs),
      labels   = function(x) round(exp(x), 2),
      name     = "Q100 ratio\n(mixture / single sample)"
    ) +
    ggplot2::coord_sf(
      xlim = sf::st_bbox(co_boundary)[c("xmin", "xmax")] + c(0, 0),
      ylim = sf::st_bbox(co_boundary)[c("ymin", "ymax")] + c(0, 0)
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
 
### Figure: map of trend in annual maxima rain fraction ----
am_trends           <- readRDS("_data/temporal_trends_am_frac_thresh_3class.rds")
am_rain_frac_trends <- am_trends[["am_rain_frac_trends"]]
am_peak_flow_trends <- am_trends[["am_peak_flow_trends"]]

rf_trend_map_data <- station_spatial %>%
  dplyr::inner_join(
    am_rain_frac_trends %>% dplyr::select(site_no, sen_slope, significant),
    by = "site_no"
  ) %>%
  dplyr::filter(!is.na(sen_slope))

# color_limit is the 90th percentile of |slope|; values beyond it are squished
# to the endpoint color, and the guide's triangular end caps mark those colors
# as "beyond this value," not an exact reading.
color_limit <- stats::quantile(abs(rf_trend_map_data$sen_slope), 0.90, na.rm = TRUE)

ggplot() +
  ggplot2::geom_sf(data = co_boundary, fill = "white", color = NA) +
  ggplot2::geom_sf(data = co_counties, fill = NA, color = "grey85", linewidth = 0.2) +
  ggplot2::geom_sf(data = co_boundary, fill = NA, color = "grey20", linewidth = 0.5) +
  ggplot2::geom_sf(
    data = rf_trend_map_data,
    aes(color = sen_slope, shape = significant), size = 2.4, alpha = 1
  ) +
  ggplot2::scale_color_gradient2(
    low = "#2166ac", mid = "grey85", high = "#b2182b", midpoint = 0,
    limits = c(-color_limit, color_limit),
    oob    = scales::oob_squish,
    name   = "Sen slope\n(rain fraction/yr)",
    guide  = legendry::guide_colbar(
      shape = "triangle",
      show  = NA,
      oob   = "squish"
    )
  ) +
  ggplot2::scale_shape_manual(
    values = c(`TRUE` = 17, `FALSE` = 16),
    labels = c(`TRUE` = "p < 0.05", `FALSE` = "Not significant"),
    name   = "Mann-Kendall test"
  ) +
  ggplot2::coord_sf(
    xlim = sf::st_bbox(co_boundary)[c("xmin", "xmax")] + c(0, 0),
    ylim = sf::st_bbox(co_boundary)[c("ymin", "ymax")] + c(0, 0)
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

### Figure: scatterplot of rain fraction trend vs elevation (site-level) ----

trend_elev_data <- am_rain_frac_trends %>%
  dplyr::select(site_no, sen_slope, significant) %>%
  dplyr::filter(!is.na(sen_slope)) %>%
  dplyr::inner_join(
    sf::st_drop_geometry(station_spatial) %>% dplyr::select(site_no, elev_mean_m),
    by = "site_no"
  )

ggplot(trend_elev_data, aes(x = elev_mean_m, y = sen_slope)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4, linetype = "dashed") +
  geom_point(aes(shape = significant, fill = significant), color = "grey20", size = 2.4, stroke = 0.4) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, color = "black") +
  scale_shape_manual(
    values = c(`TRUE` = 23, `FALSE` = 21),
    labels = c(`TRUE` = "p < 0.05", `FALSE` = "Not significant"),
    name   = "Mann-Kendall test"
  ) +
  scale_fill_manual(
    values = c(`TRUE` = "#b2182b", `FALSE` = "grey70"),
    labels = c(`TRUE` = "p < 0.05", `FALSE` = "Not significant"),
    name   = "Mann-Kendall test"
  ) +
  scale_y_continuous(name = "Sen's slope, annual maximum rain fraction (per yr)") +
  scale_x_continuous(name = "Basin mean elevation (m)") +
  labs(
    title    = "Trend in Annual Maximum Rain Fraction by Elevation",
    subtitle = "Each point is one site's long-term trend"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(size = 11, face = "italic")
  )

ggsave(file.path(fig_dir, "rain_fraction_trend_vs_elevation.png"), width = 6.5, height = 4.5, units = "in")

### Figure: trend in annual maxima streamflow ----
 
pf_trend_map_data <- station_spatial %>%
  dplyr::inner_join(
    am_peak_flow_trends %>% dplyr::select(site_no, sen_slope, significant),
    by = "site_no"
  ) %>%
  dplyr::filter(!is.na(sen_slope))
 
max_abs_slope_pf <- max(abs(pf_trend_map_data$sen_slope), na.rm = TRUE)
 
ggplot() +
  ggplot2::geom_sf(data = co_boundary, fill = "white", color = NA) +
  ggplot2::geom_sf(data = co_counties, fill = NA, color = "grey85", linewidth = 0.2) +
  ggplot2::geom_sf(data = co_boundary, fill = NA, color = "grey20", linewidth = 0.5) +
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
    xlim = sf::st_bbox(co_boundary)[c("xmin", "xmax")] + c(0, 0),
    ylim = sf::st_bbox(co_boundary)[c("ymin", "ymax")] + c(0, 0)
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

### Example station figures: site 09081600 (Crystal River AB Avalanche Creek, near Redstone, CO) ----

source("functions.R")   # needed for empirical_ri()

ffa_curves <- list(
  gpd = readRDS("_data/ffa_curves_gpd.rds"),
  lp3 = readRDS("_data/ffa_curves_lp3.rds")
)

am_annual <- readRDS("_data/annual_am_type_frac_thresh_3class.rds")

am_annual_list <- setNames(
  lapply(c("depth_thresh", "frac_thresh_3class", "frac_thresh_2class"), function(m) {
    readRDS(sprintf("_data/annual_am_type_%s.rds", m))
  }),
  c("depth_thresh", "frac_thresh_3class", "frac_thresh_2class")
)

EXAMPLE_SITE   <- "09081600"
EXAMPLE_METHOD <- "frac_thresh_3class"
EXAMPLE_FAMILY <- "gpd"

# Local color/label maps matching app.R's global type_colors/type_labels —
# defined here explicitly so this section doesn't depend on app.R's globals
# being loaded into this script's environment.
site_type_labels <- c(
  rainfall        = "Rainfall",
  snowmelt        = "Snowmelt",
  rain_on_snow    = "Rain-on-Snow",
  snow_influenced = "Snow-Influenced",
  single_sample   = "Single Sample",
  unclassified    = "Unclassified"
)

site_type_colors <- c(
  Rainfall         = "#b2182b",
  Snowmelt         = "#2166ac",
  "Rain-on-Snow"   = "#9970ab",
  "Snow-Influenced" = "#9970ab",
  "Single Sample"  = "#888888",
  Unclassified     = "#cccccc"
)

# ---- Load Script 7/8 outputs if not already in this session ----
if (!exists("ffa_classification_summary"))
  ffa_classification_summary <- readRDS("_data/ffa_classification_summary.rds")

method_labels <- c(
  depth_thresh       = "Depth threshold",
  frac_thresh_3class = "Rain fraction (3-class)",
  frac_thresh_2class = "Rain fraction (2-class)"
)

#### 1. Classification bar chart (cross-method composition, one site) ----

stack_order <- c("Unclassified", "Snowmelt", "Rain-on-Snow", "Snow-Influenced", "Rainfall")

class_plot_df <- ffa_classification_summary %>%
  dplyr::filter(site_no == EXAMPLE_SITE) %>%
  dplyr::group_by(method) %>%
  dplyr::mutate(frac = n_peaks / sum(n_peaks)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    flood_label  = dplyr::recode(flood_type, !!!site_type_labels),
    method_label = factor(
      dplyr::recode(method, !!!method_labels),
      levels = unname(method_labels[c("depth_thresh", "frac_thresh_3class", "frac_thresh_2class")])
    ),
    flood_label  = factor(flood_label, levels = stack_order[stack_order %in% unique(flood_label)])
  )

color_vec_class <- site_type_colors[names(site_type_colors) %in% levels(class_plot_df$flood_label)]
if ("Snow-Influenced" %in% names(color_vec_class)) color_vec_class["Snow-Influenced"] <- "#2166ac"

pattern_vec_class <- setNames(
  ifelse(names(color_vec_class) == "Snow-Influenced", "stripe", "none"),
  names(color_vec_class)
)

ggplot(class_plot_df, aes(x = method_label, y = frac, fill = flood_label, pattern = flood_label)) +
  geom_col_pattern(
    position = position_stack(reverse = FALSE),
    width = 0.65, colour = "white", linewidth = 0.3,
    pattern_fill             = "#9970ab",
    pattern_colour           = NA,
    pattern_density          = 0.5,
    pattern_spacing          = 0.04,
    pattern_angle            = 45,
    pattern_key_scale_factor = 0.4
  ) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.03))) +
  scale_fill_manual(values = color_vec_class, name = NULL) +
  scale_pattern_manual(values = pattern_vec_class, name = NULL) +
  labs(
    title    = "Flood Type Composition by Classification Method",
    subtitle = paste0("Station ", EXAMPLE_SITE, " — example of dashboard's Peak Classification tab"),
    x = NULL, y = "Fraction of peaks"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position    = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.subtitle      = element_text(size = 11, face = "italic"),
    axis.text.x        = element_text(angle = 25, hjust = 1)
  )

ggsave(file.path(fig_dir, "example_site_classification_composition.png"), width = 6.5, height = 4.5, units = "in")

#### 2. Return period (frequency curve) plot ----

curves_df <- ffa_curves[[EXAMPLE_FAMILY]] %>%
  dplyr::filter(site_no == EXAMPLE_SITE, method == EXAMPLE_METHOD) %>%
  dplyr::mutate(
    label = dplyr::case_when(
      curve_type == "mixture"       ~ "Mixture",
      curve_type == "single_sample" ~ "Single Sample",
      TRUE                          ~ site_type_labels[curve_type]
    ),
    color = dplyr::case_when(
      curve_type == "mixture"       ~ "black",
      curve_type == "single_sample" ~ site_type_colors["Single Sample"],
      TRUE                          ~ site_type_colors[site_type_labels[curve_type]]
    ),
    q = discharge_cfs
  )

present_types <- unique(curves_df$label)
lty_vec   <- setNames(ifelse(present_types %in% c("Mixture", "Single Sample"), "solid", "dashed"), present_types)
color_vec_ffa <- setNames(
  curves_df %>% dplyr::distinct(label, color) %>% dplyr::pull(color),
  curves_df %>% dplyr::distinct(label, color) %>% dplyr::pull(label)
)

type_curve_types <- unique(curves_df$curve_type[!curves_df$curve_type %in% c("mixture", "single_sample")])

emp_peaks <- peaks_frac %>%
  dplyr::filter(site_no == EXAMPLE_SITE, in_type_sample, as.character(flood_type) %in% type_curve_types)

emp <- if (nrow(emp_peaks) > 0) {
  empirical_ri(emp_peaks %>% dplyr::mutate(flood_type = factor(flood_type))) %>%
    dplyr::mutate(flood_label = site_type_labels[as.character(flood_type)])
} else data.frame()

ss_ids <- peaks_frac %>%
  dplyr::filter(site_no == EXAMPLE_SITE, in_single_sample) %>%
  dplyr::pull(peak_flow)

emp <- if (nrow(emp) > 0) dplyr::mutate(emp, in_ss = peak_flow %in% ss_ids) else emp
emp_ss <- if (nrow(emp) > 0) dplyr::filter(emp, in_ss) else data.frame()

ss_label <- "In Single Sample"
color_vec_p <- if (nrow(emp_ss) > 0) c(color_vec_ffa, setNames("#888888", ss_label)) else color_vec_ffa
legend_breaks <- if (nrow(emp_ss) > 0) c(names(color_vec_ffa), ss_label) else names(color_vec_ffa)
line_only <- c("Mixture", "Single Sample")
shape_overrides <- dplyr::case_when(
  legend_breaks == ss_label    ~ 1L,
  legend_breaks %in% line_only ~ NA_integer_,
  TRUE                         ~ 16L
)

ggplot() +
  geom_line(data = curves_df, aes(x = ri, y = q, colour = label, linetype = label), linewidth = 1.1, na.rm = TRUE) +
  { if (nrow(emp_ss) > 0)
      geom_point(data = emp_ss, aes(x = ri, y = peak_flow, colour = ss_label),
                 shape = 1, size = 3.5, stroke = 0.9, alpha = 0.8) } +
  { if (nrow(emp) > 0)
      geom_point(data = emp, aes(x = ri, y = peak_flow, colour = flood_label),
                 shape = 16, size = 1.8, alpha = 0.7) } +
  scale_x_log10(breaks = c(2, 5, 10, 25, 50, 100, 200, 500)) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = color_vec_p, breaks = legend_breaks, name = NULL) +
  scale_linetype_manual(values = lty_vec, guide = "none") +
  guides(colour = guide_legend(nrow = 2, override.aes = list(shape = shape_overrides))) +
  labs(
    title    = "Flood Discharge Estimates by Return Period",
    subtitle = paste0("Station ", EXAMPLE_SITE, " — ", method_labels[EXAMPLE_METHOD],
                      ", ", toupper(EXAMPLE_FAMILY), " — example of dashboard's FFA tab"),
    x = "Return Period (years)", y = "Peak Discharge (cfs)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(size = 11, face = "italic")
  )

ggsave(file.path(fig_dir, "example_site_return_period_curve.png"), width = 6.5, height = 5, units = "in")

#### 3. Peak time series ----

timeline_data <- peaks_frac %>%
  dplyr::filter(site_no == EXAMPLE_SITE) %>%
  dplyr::mutate(flood_label = site_type_labels[as.character(flood_type_raw)])

ggplot(timeline_data, aes(x = peak_date, y = peak_flow, colour = flood_label)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(values = site_type_colors, name = NULL) +
  scale_x_date(name = "Peak Date") +
  scale_y_continuous(name = "Peak Flow (cfs)", labels = scales::label_comma()) +
  labs(
    title    = "Peak Flows Over Time by Flood Type",
    subtitle = paste0("Station ", EXAMPLE_SITE, " — ", method_labels[EXAMPLE_METHOD],
                      " — example of dashboard's FFA tab")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.subtitle   = element_text(size = 11, face = "italic")
  )

ggsave(file.path(fig_dir, "example_site_peak_timeline.png"), width = 6.5, height = 4.5, units = "in")

#### 4. Rain fraction trend plot ----

rf_site_data <- am_annual %>%
  dplyr::filter(site_no == EXAMPLE_SITE) %>%
  dplyr::select(site_no, water_year, am_rain_frac)

type_site_data <- am_annual_list[[EXAMPLE_METHOD]] %>%
  dplyr::filter(site_no == EXAMPLE_SITE) %>%
  dplyr::select(site_no, water_year, am_flood_type)

site_am <- dplyr::inner_join(rf_site_data, type_site_data, by = c("site_no", "water_year")) %>%
  dplyr::mutate(flood_label = site_type_labels[as.character(am_flood_type)])

site_trend <- am_rain_frac_trends %>% dplyr::filter(site_no == EXAMPLE_SITE)

sub_txt <- if (nrow(site_trend) > 0) {
  t_row    <- site_trend[1, ]
  sig_flag <- dplyr::case_when(
    is.na(t_row$p_value)  ~ "",
    t_row$p_value < 0.001  ~ " ***",
    t_row$p_value < 0.01   ~ " **",
    t_row$p_value < 0.05   ~ " *",
    TRUE                   ~ ""
  )
  paste0("\u03c4 = ", round(t_row$tau, 3),
         "  |  Sen slope = ", signif(t_row$sen_slope, 4), " fraction/yr",
         "  |  p = ", round(t_row$p_value, 4), sig_flag,
         "  |  n = ", t_row$n_years, " years")
} else {
  "Insufficient record for trend test (\u226520 years required)"
}

am_clean <- dplyr::filter(site_am, !is.na(am_rain_frac))
sen_line <- if (nrow(site_trend) > 0 && nrow(am_clean) > 0) {
  t_row  <- site_trend[1, ]
  intcpt <- stats::median(am_clean$am_rain_frac - t_row$sen_slope * am_clean$water_year, na.rm = TRUE)
  wy_min <- min(am_clean$water_year)
  wy_max <- max(am_clean$water_year)
  data.frame(wy_min = wy_min, wy_max = wy_max,
             y_min = intcpt + t_row$sen_slope * wy_min,
             y_max = intcpt + t_row$sen_slope * wy_max)
} else NULL

present_labels_rf <- unique(site_am$flood_label[!is.na(site_am$flood_label)])
color_vec_rf <- site_type_colors[names(site_type_colors) %in% present_labels_rf]

ggplot(site_am, aes(x = water_year, y = am_rain_frac, colour = flood_label)) +
  geom_point(size = 2.5, alpha = 0.85, na.rm = TRUE) +
  { if (!is.null(sen_line))
      geom_segment(data = sen_line, aes(x = wy_min, xend = wy_max, y = y_min, yend = y_max),
                   colour = "grey20", linewidth = 1.1, inherit.aes = FALSE) } +
  scale_colour_manual(values = color_vec_rf, name = NULL) +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1)) +
  labs(
    title    = "Annual Maximum Rain Fraction by Water Year",
    subtitle = paste0(sub_txt, "\nStation ", EXAMPLE_SITE, " \u2014 example of dashboard's Trends tab"),
    x = "Water Year", y = "Annual Maxima Rain Fraction"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(size = 11, face = "italic")
  )

ggsave(file.path(fig_dir, "example_site_rain_fraction_trend.png"), width = 6.5, height = 4.5, units = "in")

