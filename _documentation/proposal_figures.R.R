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
fig_dir <- "_documentation/proposal_figures"

station_spatial     <- readRDS("_data/station_spatial.rds")
ffa_q_lp3           <- readRDS("_data/ffa_quantiles_lp3.rds")
am_rain_frac_trends <- readRDS("_data/temporal_trends_am_frac_thresh_3class.rds")[["am_rain_frac_trends"]]

flood_type_recode <- c(
  rainfall     = "Rainfall",
  snowmelt     = "Snowmelt",
  rain_on_snow = "Rain-on-Snow"
)

## Map of study area (sites colored based on regulated/unregulated) but smaller ----

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
  geom_sf(data = co_boundary_study, fill = "white", color = NA) +
  geom_sf(data = co_counties_study, fill = NA, color = "grey85", linewidth = 0.2) +
  geom_sf(data = co_boundary_study, fill = NA, color = "grey20", linewidth = 0.5) +
  geom_sf(data = study_area_data, aes(color = reg_label, shape = reg_label), size = 2, alpha = 0.8) +
  scale_color_manual(
    values = c("Unregulated" = "#1E4D2B", "Regulated" = "#C8C372", "Unknown" = "grey50"),
    name   = "Regulation status"
  ) +
  scale_shape_manual(
    values = c("Unregulated" = 16, "Regulated" = 17, "Unknown" = 4),
    name   = "Regulation status"
  ) +
  coord_sf(
    xlim = sf::st_bbox(co_boundary_study)[c("xmin", "xmax")] + c(0, 0),
    ylim = sf::st_bbox(co_boundary_study)[c("ymin", "ymax")] + c(0, 0)
  ) +
  theme_void() +
  theme(
    text = element_text(size = 8),
    legend.position    = "right",
    legend.box.just = "top",
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.spacing = unit(0, "mm"),
    legend.key.height = unit(2.5, "mm"),
    legend.key.width  = unit(3.5, "mm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    plot.margin = margin(0, 0, 0, 0, unit = "mm")
  ) 

ggsave(file.path(fig_dir, "study_area_map.png"), width = 6.5, height = 4.5, units = "in")

## Depth threshold decision space ----
OFFSET <- 0.1
type_labels <- c(
  "rainfall"      = "Rainfall",
  "snowmelt"      = "Snowmelt",
  "rain_on_snow"  = "Rain-on-snow",
  "unclassified"  = "Unclassified"
)

plot_data <- peaks_dt %>%
  filter(!is.na(total_rain_mm), !is.na(total_melt_mm)) %>%
  mutate(
    rain_plot = total_rain_mm + OFFSET,
    melt_plot = total_melt_mm + OFFSET,
    flood_type_raw = factor(flood_type_raw,
                            levels = names(type_labels),
                            labels = names(comp_colors))
  ) %>%
  # plot unclassified first so classified points render on top
  arrange(desc(flood_type_raw))

p_thresh_space <- ggplot(plot_data, aes(x = rain_plot, y = melt_plot, color = flood_type_raw)) +
  geom_point(size = 0.1, alpha = 0.4) +
  # Davenport threshold boundaries
  geom_vline(xintercept = c(5, 10) + OFFSET, linetype = "dashed",
             color = "grey30", linewidth = 0.4) +
  geom_hline(yintercept = c(5, 10) + OFFSET, linetype = "dashed",
             color = "grey30", linewidth = 0.4) +
  scale_x_log10(
    name   = expression(F[rainfall] ~ "(mm, log space)"),
    breaks = c(0.1, 1, 5, 10, 50, 100, 200) + OFFSET,
    labels = c(0, 1, 5, 10, 50, 100, 200)
  ) +
  scale_y_log10(
    name   = expression(F[snowmelt] ~ "(mm, log space)"),
    breaks = c(0.1, 1, 5, 10, 50, 100, 200) + OFFSET,
    labels = c(0, 1, 5, 10, 50, 100, 200)
  ) +
  scale_color_manual(values = comp_colors,
                     name = "Flood type") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_bw() +
  theme(
    text = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    panel.grid = element_blank(),
    panel.spacing.x = unit(1.75, "lines"),
    panel.spacing.y = unit(0.4, "lines"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    legend.box.just = "top",
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.spacing = unit(0, "mm"),
    legend.key.height = unit(2.5, "mm"),
    legend.key.width  = unit(3.5, "mm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    plot.margin = margin(1, 1, 1, 1, unit = "mm")
  ) 

p_thresh_space

## Rain fraction decision space
OFFSET <- 0

ALPHA_R <- 0.99
ALPHA_S <- round(0.312, digits = 2)

plot_data <- peaks_frac %>%
  filter(!is.na(total_rain_flux), !is.na(total_melt_flux)) %>%
  mutate(
    frac_class = case_when(
      rain_frac >= ALPHA_R ~ "Rainfall",
      rain_frac >= ALPHA_S ~ "Rain-on-Snow",
      TRUE                 ~ "Snowmelt"
    ),
    frac_class = factor(frac_class,
                        levels = c("Rainfall", "Rain-on-Snow", "Snowmelt"))
  )

p_frac_space <- ggplot(plot_data, aes(x = rain_frac, fill = frac_class)) +
  geom_histogram(binwidth = 0.01, boundary = 0)  +
  scale_fill_manual(values = comp_colors, name = "Flood type") +
  # Proposed threshold boundaries
  scale_x_continuous(name = expression(alpha)) +
  scale_y_sqrt(name = "sqrt(count)") +
  geom_vline(xintercept = c(ALPHA_S, ALPHA_R) + OFFSET, linetype = "dashed",
             color = "grey30", linewidth = 0.4) +
  theme_bw() +
  theme(
    text = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    panel.grid = element_blank(),
    panel.spacing.x = unit(1.75, "lines"),
    panel.spacing.y = unit(0.4, "lines"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    legend.box.just = "top",
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.spacing = unit(0, "mm"),
    legend.key.height = unit(2.5, "mm"),
    legend.key.width  = unit(3.5, "mm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    plot.margin = margin(1, 1, 1, 1, unit = "mm")
  )

p_frac_space

combined_plot <- (p_thresh_space + p_frac_space) + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')

ggsave(file.path(fig_dir, "classification_decision_spaces.png"),
       combined_plot, width = 6.5, height = 2.8, dpi = 150)

## Classification bar chart (cross-method composition, all sites) ----

WEIGHTING <- "site"   # "pooled" = weight by peak count; "site" = equal weight per site

stack_order <- c("Unclassified", "Snowmelt", "Rain-on-Snow", "Snow-Influenced", "Rainfall")

n_sites <- dplyr::n_distinct(ffa_classification_summary$site_no)

if (WEIGHTING == "pooled") {

  class_plot_df <- ffa_classification_summary %>%
    dplyr::group_by(method, flood_type) %>%
    dplyr::summarise(n_peaks = sum(n_peaks), .groups = "drop_last") %>%
    dplyr::group_by(method) %>%
    dplyr::mutate(frac = n_peaks / sum(n_peaks)) %>%
    dplyr::ungroup()

  subtitle_txt <- paste0("All peaks pooled across ", n_sites, " sites")

} else {

  class_plot_df <- ffa_classification_summary %>%
    dplyr::group_by(method, site_no) %>%
    dplyr::mutate(site_frac = n_peaks / sum(n_peaks)) %>%
    dplyr::group_by(method, flood_type) %>%
    dplyr::summarise(frac = mean(site_frac), .groups = "drop")

  subtitle_txt <- paste0("Mean of per-site fractions, ", n_sites, " sites")
}

class_plot_df <- class_plot_df %>%
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

p_composition <- ggplot(class_plot_df,
                        aes(x = method_label, y = frac,
                            fill = flood_label, pattern = flood_label)) +
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
  scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                     expand = expansion(mult = c(0, 0.03))) +
  scale_fill_manual(values = color_vec_class, name = NULL) +
  scale_pattern_manual(values = pattern_vec_class, name = NULL) +
  labs(
    x = NULL, y = "Fraction of peaks"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 15, hjust = 1),
    axis.title = element_text(size = 7),
    panel.grid = element_blank(),
    panel.spacing.x = unit(1.75, "lines"),
    panel.spacing.y = unit(0.4, "lines"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    legend.box.just = "top",
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.spacing = unit(0, "mm"),
    legend.key.height = unit(2.5, "mm"),
    legend.key.width  = unit(2.5, "mm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    plot.margin = margin(1, 1, 1, 1, unit = "mm")
  )

p_composition

## Figure: Map of fraction of rainfall peaks across Colorado ----

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

p_fraction <- ggplot() +
  geom_sf(data = co_boundary, fill = "white", color = NA) +
  geom_sf(data = co_counties, fill = NA, color = "grey85", linewidth = 0.2) +
  geom_sf(data = co_boundary, fill = NA, color = "grey20", linewidth = 0.5) +
  geom_sf(data = map_data, aes(color = raininess), size = 1, alpha = 0.8) +
  scale_color_gradientn(
    colours = raininess_pal,
    limits  = c(0, 1),
    name    = "Fraction"
  ) +
  coord_sf(
    xlim = sf::st_bbox(co_boundary)[c("xmin", "xmax")] + c(0, 0),
    ylim = sf::st_bbox(co_boundary)[c("ymin", "ymax")] + c(0, 0)
  ) +
  theme_void() +
  theme(
    text = element_text(size = 6),
    legend.position    = "bottom",
    #legend.box.just = "top",
    legend.box.spacing = unit(0, "mm"),
    legend.margin      = margin(0, 0, 0, 0, unit = "mm"),
    #legend.spacing = unit(0, "mm"),
    legend.key.height = unit(2.5, "mm"),
    legend.key.width  = unit(5, "mm"),
    legend.title = element_text(size = 7, margin = margin(r = 10, unit = "pt")),
    legend.text = element_text(size = 6),
    plot.margin = margin(1, 1, 1, 1, unit = "mm")
  ) 

p_combined <- (p_composition + patchwork::free(p_fraction, type = "space", side = "b")) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') +
  plot_layout(widths = c(1, 2))
ggsave(file.path(fig_dir, "composition_and_fraction_combined.png"), width = 6.5, height = 3, units = "in")
p_combined

## Figure: map of dominant flood type for the 100-yr return period, GPD, Rain fraction (3-class) ----

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
  geom_sf(data = co_boundary, fill = "white", color = NA) +
  geom_sf(data = co_counties, fill = NA, color = "grey85", linewidth = 0.2) +
  geom_sf(data = co_boundary, fill = NA, color = "grey20", linewidth = 0.5) +
  geom_sf(data = dominant_map_data, aes(color = dominant_label), size = 1, alpha = 0.8) +
  scale_color_manual(values = comp_colors, name = "Dominant flood type\n(Q100, GPD)") +
  coord_sf(
    xlim = sf::st_bbox(co_boundary)[c("xmin", "xmax")] + c(0, 0),
    ylim = sf::st_bbox(co_boundary)[c("ymin", "ymax")] + c(0, 0)
  ) +
  theme_void() +
  theme(
    text = element_text(size = 6),
    legend.position    = "right",
    legend.box.just = "top",
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.spacing = unit(0, "mm"),
    legend.key.height = unit(2.5, "mm"),
    legend.key.width  = unit(3.5, "mm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    plot.margin = margin(1, 1, 1, 1, unit = "mm")
  ) 

ggsave(file.path(fig_dir, "dominant_flood_type_100yr_gpd.png"), width = 4.6, height = 2.8, units = "in")

## Figure: scatterplot of mixed vs single sample for the 100-year, GPD ----

scatter_100_gpd <- ratio_gpd_100 %>%
  dplyr::left_join(
    sf::st_drop_geometry(meta) %>% dplyr::select(site_no, regulated),
    by = "site_no"
  ) %>%
  dplyr::mutate(reg_label = dplyr::if_else(regulated, "Regulated", "Unregulated"))

p_scatter <- ggplot(scatter_100_gpd, aes(x = single_sample, y = mixture, color = reg_label, shape = reg_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_point(alpha = 0.8, size = 1) +
  scale_x_log10(labels = scales::label_comma(), name = "Single-sample Q100 (cfs)") +
  scale_y_log10(labels = scales::label_comma(), name = "Mixture Q100 (cfs)") +
  scale_color_manual(values = c("Unregulated" = "#1E4D2B", "Regulated" = "#C8C372"), name = "Regulation") +
  scale_shape_manual(values = c("Unregulated" = 16, "Regulated" = 17), name = "Regulation") +
  theme_bw() +
  theme(
    text = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    panel.grid = element_blank(),
    panel.spacing.x = unit(1.75, "lines"),
    panel.spacing.y = unit(0.4, "lines"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    legend.box.just = "top",
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.spacing = unit(0, "mm"),
    legend.key.height = unit(2.5, "mm"),
    legend.key.width  = unit(3.5, "mm"),
    legend.text = element_text(size = 6),
    legend.title = element_blank(),
    plot.margin = margin(1, 1, 1, 1, unit = "mm")
  )

ggsave(file.path(fig_dir, "scatter_mixture_vs_single_100yr_gpd.png"), width = 3.25, height = 3, units = "in")

## Figure: column plot of binary (above/below 1.0 ratio) fraction for all return periods for GPD and Rain fraction (3-class) (mixed over single sample) ----

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

p_ratio <- ggplot(ratio_binary_counts, aes(x = return_period, y = signed_frac, fill = above_below)) +
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
  theme_bw() +
  theme(
    text = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    panel.grid = element_blank(),
    panel.spacing.x = unit(1.75, "lines"),
    panel.spacing.y = unit(0.4, "lines"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    legend.box.just = "top",
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.spacing = unit(0, "mm"),
    legend.key.height = unit(2, "mm"),
    legend.key.width  = unit(2, "mm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    plot.margin = margin(1, 1, 1, 1, unit = "mm")
  )

ggsave(file.path(fig_dir, "ratio_binary_count_all_rp_gpd.png"), width = 3.25, height = 3, units = "in")

p_combined <- (p_scatter + p_ratio) + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')
ggsave(file.path(fig_dir, "scatter_and_ratio_combined.png"), width = 6, height = 3, units = "in")
p_combined
