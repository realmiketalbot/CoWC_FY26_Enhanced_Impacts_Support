# Setup ----
library(tidyverse)
library(ggpattern)
library(legendry)
library(patchwork)
library(colorspace)
library(tigris)
library(rnaturalearth)
library(rnaturalearthdata)
library(extRemes)

peaks_dt   <- readRDS("_data/pot_peaks_classified_depth_thresh.rds")
peaks_frac <- readRDS("_data/pot_peaks_classified_frac_thresh_3class.rds")
ffa_q      <- readRDS("_data/ffa_quantiles_gpd.rds")
ffa_p      <- readRDS("_data/ffa_params_gpd.rds")
meta       <- readRDS("_data/site_metadata.rds")

climate_trends <- readRDS("_data/climate_forcing_trends.rds")
fft_trends     <- readRDS("_data/temporal_trends_frac_thresh_3class.rds")
fig_dir <- "_documentation/presentation_figures"

station_spatial     <- readRDS("_data/station_spatial.rds")
ffa_q_lp3           <- readRDS("_data/ffa_quantiles_lp3.rds")
am_rain_frac_trends <- readRDS("_data/temporal_trends_am_frac_thresh_3class.rds")[["am_rain_frac_trends"]]

flood_type_recode <- c(
  rainfall     = "Rainfall",
  snowmelt     = "Snowmelt",
  rain_on_snow = "Rain-on-Snow"
)

comp_colors <- c(
  "Rainfall"        = "#b2182b",
  "Snowmelt"        = "#2166ac",
  "Rain-on-Snow"    = "#9970ab"
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
  geom_sf(data = study_area_data, aes(color = reg_label, shape = reg_label), alpha = 0.8) +
  scale_color_manual(
    values = c("Unregulated" = "#1E4D2B", "Regulated" = "#C8C372", "Unknown" = "grey50"),
    name   = NULL
  ) +
  scale_shape_manual(
    values = c("Unregulated" = 16, "Regulated" = 17, "Unknown" = 4),
    name   = NULL
  ) +
  coord_sf(
    xlim = sf::st_bbox(co_boundary_study)[c("xmin", "xmax")] + c(0, 0),
    ylim = sf::st_bbox(co_boundary_study)[c("ymin", "ymax")] + c(0, 0)
  ) +
  theme_void() +
  theme(
    legend.position    = "bottom"
  ) 

ggsave(file.path(fig_dir, "prelim_study_area_map.png"), width = 4, height = 3, units = "in")

## Depth threshold decision space ----
OFFSET <- 0.1
type_labels2 <- c(
  "rainfall"      = "Rainfall",
  "snowmelt"      = "Snowmelt",
  "rain_on_snow"  = "Rain-on-snow",
  "unclassified"  = "Unclassified"
)

comp_colors2 <- c(
  "Rainfall"        = "#b2182b",
  "Snowmelt"        = "#2166ac",
  "Rain-on-Snow"    = "#9970ab",
  "Unclassified"    = "grey"
)

plot_data <- peaks_dt %>%
  filter(!is.na(total_rain_mm), !is.na(total_melt_mm)) %>%
  mutate(
    rain_plot = total_rain_mm + OFFSET,
    melt_plot = total_melt_mm + OFFSET,
    flood_type_raw = factor(flood_type_raw,
                            levels = names(type_labels2),
                            labels = names(comp_colors2))
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
  scale_x_continuous(name = "Rain fraction") +
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

p_frac_space

combined_plot <- (p_thresh_space + p_frac_space)

ggsave(file.path(fig_dir, "classification_decision_spaces.png"),
       combined_plot, width = 4, height = 2, dpi = 150)

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
  geom_sf(data = dominant_map_data, aes(color = dominant_label), alpha = 0.8) +
  scale_color_manual(values = comp_colors, name = NULL) +
  coord_sf(
    xlim = sf::st_bbox(co_boundary)[c("xmin", "xmax")] + c(0, 0),
    ylim = sf::st_bbox(co_boundary)[c("ymin", "ymax")] + c(0, 0)
  ) +
  theme_void() +
  theme(
    legend.position    = "bottom"
  ) 

ggsave(file.path(fig_dir, "dominant_flood_type_100yr_gpd.png"), width = 4, height = 4, units = "in")

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
    legend.position    = "bottom",
    panel.grid = element_blank()
  )

ggsave(file.path(fig_dir, "scatter_mixture_vs_single_100yr_gpd.png"), width = 4, height = 4, units = "in")

## Figure: column plot of binary (above/below 1.0 ratio) fraction for all return periods for GPD and Rain fraction (3-class) (mixed over single sample) ----

ratio_binary_counts <- ratio_all_rp_gpd %>%
  dplyr::mutate(
    above_below = dplyr::if_else(ratio >= 1, "Mixture \u2265 Sing. Samp.", "Mixture < Sing. Samp.")
  ) %>%
  dplyr::count(return_period, above_below) %>%
  dplyr::group_by(return_period) %>%
  dplyr::mutate(
    frac        = n / sum(n),
    signed_frac = dplyr::if_else(above_below == "Mixture \u2265 Sing. Samp.", frac, -frac)
  ) %>%
  dplyr::ungroup()

p_ratio <- ggplot(ratio_binary_counts, aes(x = return_period, y = signed_frac, fill = above_below)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
  scale_fill_manual(
    values = c("Mixture \u2265 Sing. Samp." = "#b2182b", "Mixture < Sing. Samp." = "#2166ac"),
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
    panel.grid   = element_blank(),
    legend.position    = "bottom"
  )

ggsave(file.path(fig_dir, "ratio_binary_count_all_rp_gpd.png"), width = 4, height = 4, units = "in")

p_combined <- (p_scatter + p_ratio) + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')
ggsave(file.path(fig_dir, "scatter_and_ratio_combined.png"), width = 6, height = 3, units = "in")
p_combined

## Figure: example streamflow plot with all local maxima with independent peaks highlighted ----
library(tidyverse)
library(pracma)

site        <- "09081600"
plot_year   <- NULL     # NULL -> pick the year with the most independent peaks
q_threshold <- 0.5
NUPS        <- 1
NDOWNS      <- 1

unit_dir  <- "/davenport-scratch/DATA/USGS/raw/Discharge/Unit"
meta_file  <- paste0("/davenport-scratch/DATA/GAGES-II/",
                     "basinchar_and_report_sept_2011/",
                     "spreadsheets-in-csv-format/conterm_basinid.txt")

find_local_maxima <- function(flow_df, n_ups = NUPS, n_downs = NDOWNS) {

  flow_df <- flow_df %>% filter(!is.na(q_daily_max)) %>% arrange(date)
  if (nrow(flow_df) < n_ups + n_downs + 1L) return(flow_df[0, ])

  gap     <- c(NA_real_, as.numeric(diff(flow_df$date), units = "days"))
  segment <- cumsum(is.na(gap) | gap != 1)

  map_dfr(unique(segment), function(s) {
    seg <- flow_df[segment == s, ]
    if (nrow(seg) < 1L) return(seg[0, ])

    eps    <- 1e-6 * max(seg$q_daily_max)
    q_jit  <- seg$q_daily_max - eps * seq_len(nrow(seg))
    padded <- c(min(q_jit) - n_ups:1, q_jit, min(q_jit) - 1:n_downs)

    pk <- pracma::findpeaks(padded, nups = n_ups, ndowns = n_downs)
    if (is.null(pk)) return(seg[0, ])

    idx <- pk[, 2] - n_ups
    seg[idx[idx >= 1L & idx <= nrow(seg)], ]
  })
}

find_independent_peaks <- function(peak_df, all_flows_df, watershed_drain_sqkm) {

  n <- nrow(peak_df)
  if (n < 2) return(peak_df)

  min_sep_days <- 5 + log(0.3861 * watershed_drain_sqkm)

  ord        <- order(all_flows_df$date)
  flow_dates <- all_flows_df$date[ord]
  flow_vals  <- all_flows_df$q_daily_max[ord]

  pk_date     <- peak_df$date
  pk_date_num <- as.numeric(pk_date)
  pk_flow     <- peak_df$q_daily_max

  min_flow_between <- function(d1, d2) {
    lo <- findInterval(d1, flow_dates) + 1L
    hi <- findInterval(d2 - 1, flow_dates)
    if (hi >= lo) min(flow_vals[lo:hi], na.rm = TRUE) else Inf
  }

  is_dependent <- function(i, j) {
    if (abs(pk_date_num[i] - pk_date_num[j]) < min_sep_days) return(TRUE)
    mf <- min_flow_between(min(pk_date[i], pk_date[j]),
                           max(pk_date[i], pk_date[j]))
    !isTRUE(mf < 0.75 * min(pk_flow[i], pk_flow[j]))
  }

  kept_idx   <- integer(0)
  kept_dates <- numeric(0)

  for (i in order(pk_flow, decreasing = TRUE)) {
    di  <- pk_date_num[i]
    pos <- findInterval(di, kept_dates)
    ok  <- TRUE
    if (pos >= 1L)                       ok <- !is_dependent(i, kept_idx[pos])
    if (ok && pos < length(kept_idx))    ok <- !is_dependent(i, kept_idx[pos + 1L])
    if (ok) {
      kept_idx   <- append(kept_idx,   i,  after = pos)
      kept_dates <- append(kept_dates, di, after = pos)
    }
  }

  peak_df[sort(kept_idx), ]
}

read_year <- function(f) {
  x <- readRDS(f)
  if (!is.data.frame(x) || nrow(x) == 0)        return(NULL)
  q_col <- grep("00060_00000$", names(x), value = TRUE)
  if (!"dateTime" %in% names(x) || !length(q_col)) return(NULL)
  if (length(q_col) > 1)
    warning(sprintf("%s: %d discharge columns (%s); using the first",
                    basename(f), length(q_col), paste(q_col, collapse = ", ")))
  tibble(date = as.Date(x$dateTime), q_cfs = as.numeric(x[[q_col[1]]])) %>%
    filter(!is.na(q_cfs), q_cfs >= 0)
}

year_files <- list.files(file.path(unit_dir, site),
                         pattern = "^\\d{4}\\.rds$", full.names = TRUE)

parsed <- map(year_files, read_year)
message(sprintf("%s: %d of %d year files contain usable unit data",
                site, sum(!map_lgl(parsed, is.null)), length(year_files)))

daily <- bind_rows(parsed) %>%
  group_by(date) %>%
  summarize(q_daily_max = max(q_cfs), .groups = "drop") %>%
  arrange(date)

drain_sqkm <- read_csv(meta_file,
                       col_types = cols(STAID = col_character(),
                                        .default = col_guess())) %>%
  filter(STAID == site) %>%
  pull(DRAIN_SQKM)

stopifnot(length(drain_sqkm) == 1, nrow(daily) > 0)

thr    <- unname(quantile(daily$q_daily_max, q_threshold, na.rm = TRUE))
theta  <- 5 + log(0.3861 * drain_sqkm)

maxima <- find_local_maxima(daily)
cand   <- filter(maxima, q_daily_max > thr)
indep  <- find_independent_peaks(cand, daily, drain_sqkm)

message(sprintf("%s: %d local maxima -> %d above threshold -> %d independent",
                site, nrow(maxima), nrow(cand), nrow(indep)))

yr <- function(d) as.integer(format(d, "%Y"))

if (is.null(plot_year)) {
  plot_year <- indep %>%
    count(y = yr(date)) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    pull(y)
}

dat <- filter(daily, yr(date) == plot_year)
pts <- maxima %>%
  filter(yr(date) == plot_year) %>%
  mutate(status = factor(
    case_when(date %in% indep$date ~ "Independent peak",
              q_daily_max > thr    ~ "Local maximum"),
    levels = c("Independent peak",
               "Local maximum")))

p <- ggplot(dat, aes(date, q_daily_max)) +
  geom_line(colour = "lightgrey", linewidth = 0.4) +
  geom_point(data = pts %>% filter(!is.na(status)),
             aes(colour = status, shape = status), size = 2) +
  scale_colour_manual(values = c("cyan", "darkgrey")) +
  scale_shape_manual(values  = c(19, 1, 3)) +
  scale_x_date(
    date_breaks = "1 month", 
    labels = function(x) substr(format(x, "%B"), 1, 1),
    limits = c(as_date("1997-04-01"), as_date("1997-09-01"))
  ) +
  labs(x = NULL, y = "Discharge (cfs)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        panel.grid = element_blank())

ggsave(file.path(fig_dir, "prelim_independent_peaks.png"), width = 4, height = 2.5, dpi = 150)

## Figure: example streamflow peaks colored by event type
plot_data <- peaks_frac %>%
  filter(site_no == site, year(peak_date) == 1997) %>%
  select(peak_date, peak_flow, flood_type_raw) %>%
  mutate(flood_type = factor(unname(flood_type_recode[flood_type_raw]),
                             levels = names(comp_colors)))

p <- ggplot() +
  geom_line(data = dat,
            aes(date, q_daily_max), colour = "lightgrey", linewidth = 0.4) +
  geom_point(data = plot_data,
             aes(x = peak_date, y = peak_flow, colour = flood_type), size = 2) +
  scale_colour_manual(values = comp_colors, na.value = "grey60", drop = FALSE) +
  scale_x_date(
    date_breaks = "1 month",
    labels = function(x) substr(format(x, "%B"), 1, 1),
    limits = c(as_date("1997-04-01"), as_date("1997-09-01"))
  ) +
  labs(x = NULL, y = "Discharge (cfs)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        panel.grid = element_blank())

ggsave(file.path(fig_dir, "prelim_classified_peaks.png"), width = 4, height = 2.5, dpi = 150)

#--------------------------------------------------------
# AM fit with GEV (2 figures)
#--------------------------------------------------------
site <- "09025000" # mixed conditions
am_data <- readRDS("_data/annual_am_type_frac_thresh_3class.rds") %>%
  filter(site_no == site) %>%
  select(water_year, peak_flow, am_flood_type)

# --- GEV fit (MLE) to the full annual maximum series ---
am_clean <- am_data %>% filter(!is.na(peak_flow))

gev_fit <- fevd(am_clean$peak_flow, type = "GEV", method = "MLE")

# --- fitted return-level curve ---
rp_grid <- exp(seq(log(1.01), log(200), length.out = 400))

rl_curve <- tibble(
  return_period = rp_grid,
  return_level  = as.numeric(return.level(gev_fit, return.period = rp_grid))
)

# --- empirical plotting positions (Cunnane) ---
am_points <- am_clean %>%
  arrange(desc(peak_flow)) %>%
  mutate(
    rank          = row_number(),
    p_exceed      = (rank - 0.4) / (n() + 0.2),
    return_period = 1 / p_exceed,
    flood_type    = factor(flood_type_recode[am_flood_type],
                           levels = names(comp_colors))
  )

rp_breaks <- c(1, 2, 5, 10, 25, 50, 100, 200)

# --- shared theme ---
theme_rl <- theme_bw(base_size = 10) +
  theme(
    panel.grid      = element_blank(),
    panel.border    = element_rect(color = "grey30", linewidth = 0.4),
    axis.ticks      = element_line(linewidth = 0.3),
    legend.title    = element_blank(),
    legend.key      = element_blank(),
    legend.background = element_blank()
  )

# --- Figure 1: black points ---
p1 <- ggplot() +
  geom_line(data = rl_curve,
            aes(return_period, return_level),
            linewidth = 0.5) +
  geom_point(data = am_points,
             aes(return_period, peak_flow),
             size = 1.5, shape = 16) +
  scale_x_log10(breaks = rp_breaks, labels = rp_breaks,
                expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Return period (yr)", y = "Peak flow (cfs)") +
  theme_rl

ggsave(filename = file.path(fig_dir, "prelim_gev_curve.png"),
       width = 4, height = 3, dpi = 150)

# --- Figure 2: black points with vertical dashed line ---
p2 <- ggplot() +
  geom_line(data = rl_curve,
            aes(return_period, return_level),
            linewidth = 0.5) +
  geom_point(data = am_points,
             aes(return_period, peak_flow),
             size = 1.5, shape = 16) +
  geom_vline(xintercept = 100, color = "blue", linetype = "dashed") +
  scale_x_log10(breaks = rp_breaks, labels = rp_breaks,
                expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Return period (yr)", y = "Peak flow (cfs)") +
  theme_rl

ggsave(filename = file.path(fig_dir, "prelim_gev_curve_with_dash.png"),
       width = 4, height = 3, dpi = 150)

# --- Figure 3: points colored by flood type ---
p3 <- ggplot() +
  geom_line(data = rl_curve,
            aes(return_period, return_level),
            linewidth = 0.5, color = "grey20") +
  geom_point(data = am_points,
             aes(return_period, peak_flow, color = flood_type),
             size = 1.5, shape = 16) +
  scale_x_log10(breaks = rp_breaks, labels = rp_breaks,
                expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(values = comp_colors, na.translate = FALSE) +
  labs(x = "Return period (yr)", y = "Peak flow (cfs)") +
  theme_rl +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.margin = margin(0, 0, 0, 0))

ggsave(filename = file.path(fig_dir, "prelim_gev_curve_classified.png"),
       width = 4, height = 3, dpi = 150)

## Figure: proposed mixture model expected quantile discrepancy diagram ----
fig_w <- 4.00
fig_h <- 4.00

## ---- palette -------------------------------------------------------------
comp_colors <- c(
  "Rainfall" = "#b2182b",
  "Snowmelt" = "#2166ac"
)

## Light fill derived from each component colour, so the palette stays at
## two hues plus grey.
tint <- function(col, p) {
  rgb(colorRamp(c("white", col))(p), maxColorValue = 255)
}
comp_fills <- vapply(comp_colors, tint, character(1), p = 0.12)

grey_bg <- "#F0F0F0"; grey_st <- "#808080"; grey_tx <- "#333333"
neut_bg <- "#FFFFFF"; neut_st <- "#BFBFBF"

## ---- type sizes ----------------------------------------------------------
cex_title <- 0.95   # node titles and the outer-loop heading
cex_sub   <- 0.80   # subtitles, captions, edge label

## ---- helpers -------------------------------------------------------------

## Rounded rectangle, specified by its top-left corner.
rrect <- function(x, y, w, h, r = 7, fill = NA, border = "black",
                  lty = 1, lwd = 1) {
  a  <- seq(0, pi / 2, length.out = 12)
  px <- c((x + w - r) + r * cos(a - pi / 2),
          (x + w - r) + r * cos(a),
          (x + r)     + r * cos(a + pi / 2),
          (x + r)     + r * cos(a + pi))
  py <- c((y + r)     + r * sin(a - pi / 2),
          (y + h - r) + r * sin(a),
          (y + h - r) + r * sin(a + pi / 2),
          (y + r)     + r * sin(a + pi))
  polygon(px, py, col = fill, border = border, lty = lty, lwd = lwd)
}

## Two-line node: bold title, muted subtitle.
node2 <- function(x, y, w, h, title, sub, fill, border, tcol) {
  rrect(x, y, w, h, fill = fill, border = border, lwd = 1)
  text(x + w / 2, y + 17, title, col = tcol, font = 2, cex = cex_title)
  text(x + w / 2, y + 36, sub,   col = tcol, font = 1, cex = cex_sub)
}

## Single-line node.
node1 <- function(x, y, w, h, title, fill, border, tcol) {
  rrect(x, y, w, h, fill = fill, border = border, lwd = 1)
  text(x + w / 2, y + h / 2, title, col = tcol, font = 2, cex = cex_title)
}

arr <- function(x1, y1, x2, y2, col = grey_st) {
  arrows(x1, y1, x2, y2, length = 0.06, angle = 22, col = col, lwd = 1)
}

## Column of nodes for one flood type.
component_column <- function(x, label) {
  col  <- comp_colors[[label]]
  fill <- comp_fills[[label]]
  cx   <- x + 77.5
  text(cx, 100, label, col = col, font = 2, cex = cex_title)
  node2(x, 114, 155, 52, "Run EQD",
        expression(paste(B, " fits, picks threshold")), fill, col, col)
  arr(cx, 168, cx, 190)
  node2(x, 194, 155, 52, "Algorithm 1",
        expression(paste(B[1], " parameter draws")), fill, col, col)
}

## --------------------------------------------------------------------------
## ORIGINAL 
## ---- canvas --------------------------------------------------------------
out_file <- file.path(fig_dir, "mixture_eqd_bootstrap_diagram_original.png")
png(out_file, width = fig_w, height = fig_h, units = "in", res = 300,
    pointsize = 10, bg = "white")
par(mar = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(0, 400), ylim = c(448, 0), xaxs = "i", yaxs = "i",
            asp = 1)

## ---- outer loop container ------------------------------------------------
rrect(24, 24, 352, 324, r = 16, fill = grey_bg, border = grey_st)
text(38, 50, expression(bold(paste("Outer loop, repeated ", B[2], " times"))),
     adj = c(0, 0.5), col = grey_tx, cex = cex_title)
text(38, 70, "resample years, keep flood-type labels", adj = c(0, 0.5),
     col = grey_st, cex = cex_sub)

## ---- components ----------------------------------------------------------
component_column(38,  "Rainfall")
component_column(207, "Snowmelt")

## ---- combine (inside the outer loop) -------------------------------------
arr(115.5, 248, 150, 272)
arr(284.5, 248, 250, 272)
rrect(60, 276, 280, 52, fill = neut_bg, border = neut_st)
text(200, 293, "Combine the two fits", col = grey_tx, font = 2,
     cex = cex_title)
text(200, 312, "solve product of CDFs for return level",
     col = grey_st, cex = cex_sub)

## ---- collected across replicates (outside the loop) ----------------------
arr(200, 348, 200, 378)
text(210, 363, expression(paste(B[2] %*% B[1], " values")),
     adj = c(0, 0.5), col = grey_st, cex = cex_sub)
node1(75, 382, 250, 42, "Return level CIs",
      neut_bg, neut_st, grey_tx)

invisible(dev.off())
cat("written:", out_file, "\n")

## Figure: mixture EQD bootstrap nested inside a sensitivity loop.
## Base R only. If you are appending this to the file that already defines
## tint(), rrect(), node2(), node1(), arr() and component_column(), delete
## the helper section below and keep only the canvas block.
##
## Canvas is 400 x 486 units at 100 units per inch.

fig_w <- 4.00
fig_h <- 4.86

## ---- palette -------------------------------------------------------------
comp_colors <- c(
  "Rainfall" = "#b2182b",
  "Snowmelt" = "#2166ac"
)

tint <- function(col, p) {
  rgb(colorRamp(c("white", col))(p), maxColorValue = 255)
}
comp_fills <- vapply(comp_colors, tint, character(1), p = 0.12)

outer_bg <- "#FAFAFA"   # sensitivity loop
grey_bg  <- "#EBEBEB"   # bootstrap loop, one step darker so nesting reads
grey_st  <- "#808080"; grey_tx <- "#333333"
neut_bg  <- "#FFFFFF"; neut_st <- "#BFBFBF"

## ---- type sizes ----------------------------------------------------------
cex_title <- 0.95
cex_sub   <- 0.80

## ---- helpers -------------------------------------------------------------

rrect <- function(x, y, w, h, r = 7, fill = NA, border = "black",
                  lty = 1, lwd = 1) {
  a  <- seq(0, pi / 2, length.out = 12)
  px <- c((x + w - r) + r * cos(a - pi / 2),
          (x + w - r) + r * cos(a),
          (x + r)     + r * cos(a + pi / 2),
          (x + r)     + r * cos(a + pi))
  py <- c((y + r)     + r * sin(a - pi / 2),
          (y + h - r) + r * sin(a),
          (y + h - r) + r * sin(a + pi / 2),
          (y + r)     + r * sin(a + pi))
  polygon(px, py, col = fill, border = border, lty = lty, lwd = lwd)
}

## Two-line node. Text offsets are proportional to h so the same helper
## works at any node height.
node2 <- function(x, y, w, h, title, sub, fill, border, tcol) {
  rrect(x, y, w, h, fill = fill, border = border, lwd = 1)
  text(x + w / 2, y + 0.35 * h, title, col = tcol, font = 2, cex = cex_title)
  text(x + w / 2, y + 0.72 * h, sub,   col = tcol, font = 1, cex = cex_sub)
}

node1 <- function(x, y, w, h, title, fill, border, tcol) {
  rrect(x, y, w, h, fill = fill, border = border, lwd = 1)
  text(x + w / 2, y + h / 2, title, col = tcol, font = 2, cex = cex_title)
}

arr <- function(x1, y1, x2, y2, col = grey_st) {
  arrows(x1, y1, x2, y2, length = 0.06, angle = 22, col = col, lwd = 1)
}

## Column of nodes for one flood type.
component_column <- function(x, label) {
  col  <- comp_colors[[label]]
  fill <- comp_fills[[label]]
  cx   <- x + 76
  text(cx, 125, label, col = col, font = 2, cex = cex_title)
  node2(x, 134, 152, 46, "Run EQD",
        expression(paste(B, " fits, picks threshold")), fill, col, col)
  arr(cx, 182, cx, 198)
  node2(x, 202, 152, 46, "Algorithm 1",
        expression(paste(B[1], " parameter draws")), fill, col, col)
}

## ---- canvas --------------------------------------------------------------
out_file <- file.path(fig_dir, "mixture_eqd_bootstrap_diagram_sensitivity.png")
png(out_file, width = fig_w, height = fig_h, units = "in", res = 300,
    pointsize = 10, bg = "white")
par(mar = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(0, 400), ylim = c(486, 0), xaxs = "i", yaxs = "i",
            asp = 1)

## ---- sensitivity loop (outermost) ----------------------------------------
rrect(16, 16, 368, 392, r = 16, fill = outer_bg, border = grey_st)
text(30, 37, expression(bold("Sensitivity loop, over settings")),
     adj = c(0, 0.5), col = grey_tx, cex = cex_title)
text(30, 53, "classification threshold x min. sample size", adj = c(0, 0.5),
     col = grey_st, cex = cex_sub)

## ---- bootstrap loop ------------------------------------------------------
rrect(28, 66, 344, 267, r = 14, fill = grey_bg, border = grey_st)
text(40, 87, expression(bold(paste("Bootstrap loop, repeated ", B[2],
                                   " times"))),
     adj = c(0, 0.5), col = grey_tx, cex = cex_title)
text(40, 103, "resample years, keep flood-type labels", adj = c(0, 0.5),
     col = grey_st, cex = cex_sub)

## ---- components ----------------------------------------------------------
component_column(40,  "Rainfall")
component_column(208, "Snowmelt")

## ---- combine (inside the bootstrap loop) ---------------------------------
arr(116, 250, 150, 270)
arr(284, 250, 250, 270)
rrect(58, 274, 284, 46, fill = neut_bg, border = neut_st)
text(200, 290, "Combine the two fits", col = grey_tx, font = 2,
     cex = cex_title)
text(200, 307, "solve product of CDFs for return level",
     col = grey_st, cex = cex_sub)

## ---- one CI per setting (inside the sensitivity loop) --------------------
arr(200, 333, 200, 355)
text(210, 344, expression(paste(B[2] %*% B[1], " values")),
     adj = c(0, 0.5), col = grey_st, cex = cex_sub)
node1(85, 359, 230, 36, "Return level CIs", neut_bg, neut_st, grey_tx)

## ---- collected across settings (outside everything) ----------------------
arr(200, 408, 200, 430)
text(210, 419, "one CI per setting", adj = c(0, 0.5), col = grey_st,
     cex = cex_sub)
node1(60, 434, 280, 36, "Compare CIs across settings",
      neut_bg, neut_st, grey_tx)

invisible(dev.off())
cat("written:", out_file, "\n")