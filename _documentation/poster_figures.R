#!/usr/bin/env Rscript
# ============================================================================
# AGU 2026 poster figures (brainstorm set)
# ============================================================================
# Every candidate panel for the poster, built from the M1 Poisson–GPD outputs
# (Script 7). Only the new results are shown; there is no comparison with the
# initial results. Return periods are ANNUAL. Results beyond 100 yr are
# marked indicative until the M3 bootstrap intervals exist (plan D-F).
#
# Run from the repo root:  Rscript _documentation/poster_figures.R
# Sensitivity panels (S3, S4) need _data/sensitivity_ratio_summary.rds from
# _documentation/poster_sensitivity.R and are skipped if it is missing.
#
# Output: _documentation/poster_figures/*.png
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(patchwork)
})
source("functions.R")

fig_dir <- "_documentation/poster_figures"
dir.create(fig_dir, showWarnings = FALSE)

# ---- Data -------------------------------------------------------------------
ffa_q   <- readRDS("_data/ffa_quantiles_pgpd.rds")
ffa_c   <- readRDS("_data/ffa_curves_pgpd.rds")
ffa_p   <- readRDS("_data/ffa_params_pgpd.rds")
peaks   <- readRDS("_data/pot_peaks_classified_frac_thresh_3class.rds")
meta    <- readRDS("_data/site_metadata.rds")
stations <- readRDS("_data/station_spatial.rds")
qual    <- readRDS("_data/site_qualifying_years.rds") %>% filter(selected, qualifies)
qual_years <- split(qual$water_year, qual$site_no)
n_qual  <- qual %>% count(site_no, name = "n_qual_years")

MAIN_METHOD <- "frac_thresh_3class"
SHOWCASE_SITE <- "09115500"   # Tomichi Creek at Sargents (see poster_showcase.R)
RP_BREAKS   <- c(2, 5, 10, 25, 50, 100, 200, 500)

# Labels and colours follow the Hydrology Days 2026 talk (Set1 type colours;
# "Fixed" = depth threshold, "Adaptive" = rain fraction). Method colours avoid
# the type colours.
method_labels <- c(depth_thresh       = "Fixed",
                   frac_thresh_3class = "Adaptive (3-class)",
                   frac_thresh_2class = "Adaptive (2-class)")
method_colors <- c("Adaptive (3-class)" = "black",
                   "Adaptive (2-class)" = "#FF7F00",
                   "Fixed"              = "#A65628")

type_labels <- c(rainfall = "Rainfall", snowmelt = "Snowmelt", rain_on_snow = "Rain-on-Snow",
                 snow_influenced = "Non-Rainfall", unclassified = "Unclassified",
                 mixture = "Mixture", single_sample = "Single Sample")
type_colors <- c("Rainfall" = "#E41A1C", "Snowmelt" = "#377EB8", "Rain-on-Snow" = "#4DAF4A",
                 "Non-Rainfall" = "#984EA3", "Unclassified" = "#CCCCCC",
                 "Mixture" = "black", "Single Sample" = "grey55")

theme_poster <- theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 11, colour = "grey40"),
        legend.position = "bottom")

# Shading for the indicative range (> 100 yr)
indicative_layer <- function(ymin = -Inf, ymax = Inf) list(
  annotate("rect", xmin = 100, xmax = Inf, ymin = ymin, ymax = ymax, fill = "grey92"),
  geom_vline(xintercept = 100, linetype = "dashed"),
  annotate("text", x = 220, y = Inf, label = "indicative", vjust = 1.5, size = 3.5,
           colour = "grey35", fontface = "italic"))
rp_axis <- scale_x_log10(breaks = RP_BREAKS, name = "Annual return period (years)")

save_fig <- function(name, p, w = 7, h = 5) {
  ggsave(file.path(fig_dir, name), p, width = w, height = h, units = "in", dpi = 200)
  message("  ", name)
}

# Colorado boundary: tigris when available (as in report_figures.R), else maps.
co_layers <- function() {
  if (requireNamespace("tigris", quietly = TRUE)) {
    st  <- tigris::states(cb = TRUE, class = "sf", progress_bar = FALSE) %>% filter(STUSPS == "CO")
    cty <- tigris::counties(state = "CO", cb = TRUE, class = "sf", progress_bar = FALSE)
  } else {
    st  <- st_as_sf(maps::map("state", "colorado", plot = FALSE, fill = TRUE))
    cty <- st_as_sf(maps::map("county", "colorado", plot = FALSE, fill = TRUE))
  }
  list(state = st_transform(st, 4326), county = st_transform(cty, 4326))
}
co <- co_layers()
co_map <- function(data, ...) {
  bb <- st_bbox(co$state)
  ggplot() +
    geom_sf(data = co$state, fill = "white", colour = NA) +
    geom_sf(data = co$county, fill = NA, colour = "grey85", linewidth = 0.2) +
    geom_sf(data = co$state, fill = NA, colour = "grey20", linewidth = 0.5) +
    geom_sf(data = data, ...) +
    coord_sf(xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")]) +
    theme_void(base_size = 13) +
    theme(legend.position = "right", plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 11, colour = "grey40"))
}

# Mixture / single-sample ratio per site x method x return period
ratios <- ffa_q %>%
  filter(curve_type %in% c("mixture", "single_sample")) %>%
  select(site_no, method, return_period, curve_type, discharge_cfs) %>%
  pivot_wider(names_from = curve_type, values_from = discharge_cfs) %>%
  filter(is.finite(mixture), is.finite(single_sample)) %>%
  mutate(ratio = mixture / single_sample,
         method_label = factor(method_labels[method], levels = method_labels))

# Fit pool (eligible peaks in qualifying water years), 3-class labels
pool <- peaks %>%
  filter(fit_eligible, is.finite(peak_flow), peak_flow > 0) %>%
  semi_join(qual, by = c("site_no", "water_year"))

message("Writing figures to ", fig_dir)

# ============================================================================
# CONTEXT
# ============================================================================

## 01 Study area: sites coloured by qualifying years ----
p <- co_map(stations %>% left_join(n_qual, by = "site_no"),
            aes(colour = n_qual_years), size = 2.2) +
  scale_colour_viridis_c(name = "Qualifying\nwater years", option = "D", direction = -1) +
  labs(title = "Study area and streamgage network",
       subtitle = sprintf("%d USGS gages, >= 20 qualifying water years (WY1986–2023)", nrow(stations)))
save_fig("01_study_area_qualifying_years.png", p, 7, 5)

## 02 Network composition of fitted peaks by method ----
comp <- bind_rows(
  readRDS("_data/pot_peaks_classified_depth_thresh.rds") %>%
    mutate(method = "depth_thresh", lab = as.character(flood_type_raw)),
  peaks %>% mutate(method = "frac_thresh_3class", lab = as.character(flood_type_raw)),
  remap_2class(peaks) %>% mutate(method = "frac_thresh_2class", lab = as.character(flood_type))
) %>%
  mutate(lab = if_else(fit_eligible, lab, "unclassified")) %>%
  count(method, lab) %>% group_by(method) %>% mutate(pct = n / sum(n)) %>% ungroup() %>%
  mutate(method = factor(method_labels[method], levels = method_labels),
         lab = factor(type_labels[lab], levels = c("Unclassified", "Snowmelt", "Non-Rainfall",
                                                   "Rain-on-Snow", "Rainfall")))
p <- ggplot(comp, aes(method, pct, fill = lab)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(pct > 0.04, scales::percent(pct, 1), "")),
            position = position_stack(vjust = 0.5), colour = "white", size = 3.6) +
  scale_fill_manual(values = type_colors, name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "Share of candidate peaks", title = "Flood-type composition by classification method",
       subtitle = "All candidate peaks. Unclassified: Fixed criteria not met, or insufficient WLDAS data") +
  theme_poster
save_fig("02_composition_by_method.png", p, 7, 5)

## 03 Share of rainfall peaks per site (map) ----
rain_share <- pool %>% group_by(site_no) %>%
  summarise(rain_share = mean(flood_type_raw == "rainfall"), .groups = "drop")
p <- co_map(stations %>% inner_join(rain_share, by = "site_no"),
            aes(colour = rain_share), size = 2.2) +
  scale_colour_distiller(palette = "RdBu", direction = -1, labels = scales::percent,
                         name = "Rainfall\npeaks") +
  labs(title = "Share of flood peaks driven by rainfall",
       subtitle = "Adaptive (3-class); eligible peaks in qualifying water years")
save_fig("03_rainfall_share_map.png", p, 7, 5)

# ============================================================================
# MAIN RESULT: mixture vs single sample
# ============================================================================

## 04 Ratio by return period (median and IQR), by method ----
ratio_sum <- ratios %>% group_by(method_label, return_period) %>%
  summarise(median = median(ratio), q25 = quantile(ratio, .25), q75 = quantile(ratio, .75),
            share_ge1 = mean(ratio >= 1), n = n(), .groups = "drop")
p <- ggplot(ratio_sum, aes(return_period, median, colour = method_label, fill = method_label)) +
  indicative_layer() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  rp_axis +
  scale_colour_manual(values = method_colors, name = NULL, aesthetics = c("colour", "fill")) +
  labs(y = "Mixture / single-sample discharge", title = "Where do mixture and single-sample estimates diverge?",
       subtitle = "Median across sites with interquartile range") +
  theme_poster
save_fig("04_ratio_by_return_period.png", p, 7, 5)

## 04b Same, main method only, with the site-level spread ----
p <- ratios %>% filter(method == MAIN_METHOD) %>%
  ggplot(aes(factor(return_period), ratio)) +
  annotate("rect", xmin = 6.5, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_violin(fill = "#1b7837", alpha = 0.25, colour = NA, scale = "width") +
  geom_boxplot(width = 0.18, outlier.size = 0.8, fill = "white") +
  scale_y_log10() +
  labs(x = "Annual return period (years)", y = "Mixture / single-sample discharge (log scale)",
       title = "Mixture vs single-sample ratio across sites",
       subtitle = "Adaptive (3-class); shaded = indicative (> 100 yr)") +
  theme_poster
save_fig("04b_ratio_distribution_3class.png", p, 7, 5)

## 05 Share of sites where the mixture is higher ----
p <- ggplot(ratio_sum, aes(return_period, share_ge1, colour = method_label)) +
  indicative_layer() +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  rp_axis +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_colour_manual(values = method_colors, name = NULL) +
  labs(y = "Sites where mixture >= single sample", title = "The mixture is higher mainly in the tail") +
  theme_poster
save_fig("05_share_mixture_higher.png", p, 7, 5)

## 06 Ratio map at 100 yr ----
r100 <- ratios %>% filter(method == MAIN_METHOD, return_period == 100)
lim  <- max(abs(log(r100$ratio)))
p <- co_map(stations %>% inner_join(r100, by = "site_no"), aes(colour = log(ratio)), size = 2.4) +
  scale_colour_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b", midpoint = 0,
                         limits = c(-lim, lim), labels = function(x) round(exp(x), 2),
                         name = "Q100 ratio\n(mixture /\nsingle sample)") +
  labs(title = "Mixture vs single-sample 100-year flood",
       subtitle = "Adaptive (3-class), Poisson–GPD. Circled: showcase site (Tomichi Creek)") +
  geom_sf(data = stations %>% filter(site_no == SHOWCASE_SITE), shape = 21, size = 5.5,
          stroke = 1.1, colour = "black", fill = NA)
save_fig("06_ratio_map_q100.png", p, 7, 5)

## 07 Mixture vs single sample, Q100 scatter (by elevation) ----
p <- r100 %>% left_join(meta %>% select(site_no, elev_mean_m), by = "site_no") %>%
  ggplot(aes(single_sample, mixture, colour = elev_mean_m)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_x_log10(labels = scales::label_comma(), name = "Single-sample Q100 (cfs)") +
  scale_y_log10(labels = scales::label_comma(), name = "Mixture Q100 (cfs)") +
  scale_colour_viridis_c(option = "C", name = "Mean basin\nelevation (m)") +
  labs(title = "Mixture vs single-sample 100-year flood") +
  theme_poster + theme(legend.position = "right")
save_fig("07_scatter_q100.png", p, 7, 5.5)

## 08 Ratio vs elevation, coloured by rainfall share ----
p <- r100 %>% left_join(meta %>% select(site_no, elev_mean_m), by = "site_no") %>%
  left_join(rain_share, by = "site_no") %>%
  ggplot(aes(elev_mean_m, ratio, colour = rain_share)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "loess", formula = y ~ x, colour = "black", se = FALSE, linewidth = 0.7) +
  scale_y_log10() +
  scale_colour_distiller(palette = "RdBu", direction = -1, labels = scales::percent, name = "Rainfall\npeaks") +
  labs(x = "Mean basin elevation (m)", y = "Q100 ratio (mixture / single sample)",
       title = "Ratio vs elevation") +
  theme_poster + theme(legend.position = "right")
save_fig("08_ratio_vs_elevation.png", p, 7, 5)

# ============================================================================
# WHICH FLOOD TYPE CONTROLS THE TAIL
# ============================================================================
type_q <- ffa_q %>%
  filter(method == MAIN_METHOD, !curve_type %in% c("mixture", "single_sample"), is.finite(discharge_cfs))
dominant <- type_q %>% group_by(site_no, return_period) %>%
  slice_max(discharge_cfs, n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(lab = factor(type_labels[curve_type], levels = c("Snowmelt", "Rain-on-Snow", "Rainfall")))

## 09 Dominant type by return period ----
p <- dominant %>% count(return_period, lab) %>% group_by(return_period) %>% mutate(pct = n / sum(n)) %>%
  ggplot(aes(factor(return_period), pct, fill = lab)) +
  geom_col(width = 0.75) +
  annotate("rect", xmin = 6.5, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "white", alpha = 0.35) +
  scale_fill_manual(values = type_colors, name = "Largest flood type") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = "Annual return period (years)", y = "Share of sites",
       title = "Which flood type produces the largest floods?",
       subtitle = "Type with the highest quantile at each site; Adaptive (3-class)") +
  theme_poster
save_fig("09_dominant_type_by_rp.png", p, 7, 5)

## 10 Dominant type map at Q2 and Q100 ----
dom_map <- stations %>% inner_join(dominant %>% filter(return_period %in% c(2, 100)), by = "site_no") %>%
  mutate(rp = factor(paste0(return_period, "-year"), levels = c("2-year", "100-year")))
p <- co_map(dom_map, aes(colour = lab), size = 1.9) + facet_wrap(~ rp) +
  scale_colour_manual(values = type_colors, name = "Largest\nflood type") +
  labs(title = "Dominant flood type: common vs rare floods")
save_fig("10_dominant_type_map_q2_q100.png", p, 10, 4.5)

# ============================================================================
# EXAMPLE SITES
# ============================================================================
r3 <- ratios %>% filter(method == MAIN_METHOD, return_period == 100)
three_type <- ffa_p %>% filter(method == MAIN_METHOD, curve_type %in% c("rainfall", "snowmelt", "rain_on_snow")) %>%
  count(site_no) %>% filter(n == 3) %>% pull(site_no)
r3t <- r3 %>% filter(site_no %in% three_type)
example_sites <- unique(c(
  if ("09081600" %in% r3$site_no) "09081600",                          # report example (Crystal R.)
  r3t$site_no[which.min(abs(r3t$ratio - median(r3$ratio)))],          # typical ratio
  r3$site_no[which.max(r3$ratio)]                                       # largest divergence
))

site_name <- function(s) {
  nm <- meta$station_name[meta$site_no == s]
  str_trunc(str_to_title(nm), 42)
}

plot_example <- function(s) {
  cur <- ffa_c %>% filter(site_no == s, method == MAIN_METHOD, is.finite(discharge_cfs)) %>%
    mutate(lab = factor(type_labels[curve_type], levels = names(type_colors)))
  prm <- ffa_p %>% filter(site_no == s, method == MAIN_METHOD)
  sp  <- pool %>% filter(site_no == s)
  pts <- prm %>% filter(curve_type != "single_sample") %>%
    pmap_dfr(function(curve_type, u, n_years, ...) {
      empirical_annual_rp(sp$peak_flow[sp$flood_type_raw == curve_type], u, n_years) %>%
        mutate(lab = type_labels[curve_type])
    }) %>% mutate(lab = factor(lab, levels = names(type_colors)))
  ratio <- r3$ratio[r3$site_no == s]
  ggplot() +
    indicative_layer() +
    geom_point(data = pts, aes(ri, peak_flow, colour = lab), size = 1.4, alpha = 0.6) +
    geom_line(data = cur, aes(ri, discharge_cfs, colour = lab, linetype = lab), linewidth = 0.9) +
    scale_x_log10(breaks = RP_BREAKS, limits = c(1.1, 500), name = "Annual return period (years)") +
    scale_y_log10(labels = scales::label_comma(), name = "Peak discharge (cfs)") +
    scale_colour_manual(values = type_colors, name = NULL, drop = TRUE) +
    scale_linetype_manual(values = c("Rainfall" = "solid", "Snowmelt" = "solid", "Rain-on-Snow" = "solid",
                                     "Mixture" = "solid", "Single Sample" = "dashed"), name = NULL) +
    labs(title = sprintf("%s (%s)", site_name(s), s),
         subtitle = sprintf("Q100 ratio = %.2f", ratio)) +
    theme_poster + theme(plot.title = element_text(size = 11))
}
message("  example sites: ", paste(example_sites, collapse = ", "))
for (s in example_sites) save_fig(sprintf("11_example_curves_%s.png", s), plot_example(s), 6.5, 5)
p <- wrap_plots(lapply(example_sites, plot_example), nrow = 1, guides = "collect") &
  theme(legend.position = "bottom")
save_fig("11_example_curves_panel.png", p, 15, 5.5)

# ============================================================================
# DIAGNOSTICS
# ============================================================================

## 12 Shape parameter distribution ----
xi_df <- ffa_p %>% filter(method == MAIN_METHOD) %>%
  mutate(lab = factor(type_labels[curve_type], levels = c("Rainfall", "Snowmelt", "Rain-on-Snow", "Single Sample")))
p <- ggplot(xi_df, aes(lab, xi, fill = lab)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = Inf, fill = "#fddbc7", alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_violin(alpha = 0.4, colour = NA, scale = "width") +
  geom_boxplot(width = 0.18, fill = "white", outlier.size = 0.8) +
  scale_fill_manual(values = type_colors, guide = "none") +
  labs(x = NULL, y = expression("GPD shape " * xi),
       title = "Fitted tail shape by flood type",
       subtitle = sprintf("Shaded: xi > 0.5 (%.0f%% of type fits); rain fraction (3-class)",
                          100 * mean(xi_df$xi[xi_df$curve_type != "single_sample"] > 0.5))) +
  theme_poster
save_fig("12_shape_parameter.png", p, 7, 5)

## 13 Dispersion index of annual exceedance counts (S1; Proposal App. A.2) ----
disp <- ffa_p %>% filter(method == MAIN_METHOD) %>%
  pmap_dfr(function(site_no, curve_type, u, ...) {
    wys <- qual_years[[site_no]]
    sp  <- pool[pool$site_no == site_no, ]
    x   <- if (curve_type == "single_sample") sp else sp[sp$flood_type_raw == curve_type, ]
    cnt <- tabulate(match(x$water_year[x$peak_flow > u], wys), nbins = length(wys))
    M   <- length(wys)
    tibble(site_no = site_no, curve_type = curve_type, M = M, D = var(cnt) / mean(cnt),
           lo = qchisq(0.025, M - 1) / (M - 1), hi = qchisq(0.975, M - 1) / (M - 1))
  }) %>%
  mutate(status = case_when(D < lo ~ "Under-dispersed", D > hi ~ "Over-dispersed", TRUE ~ "Consistent with Poisson"),
         lab = factor(type_labels[curve_type], levels = c("Rainfall", "Snowmelt", "Rain-on-Snow", "Single Sample")))
disp_share <- disp %>% group_by(lab) %>%
  summarise(txt = sprintf("%.0f%% within", 100 * mean(status == "Consistent with Poisson")), .groups = "drop")
p <- ggplot(disp, aes(lab, D)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_jitter(aes(colour = status), width = 0.25, height = 0, size = 1.3, alpha = 0.7) +
  geom_text(data = disp_share, aes(lab, Inf, label = txt), vjust = 1.4, size = 3.6) +
  scale_colour_manual(values = c("Consistent with Poisson" = "grey55", "Over-dispersed" = "#b2182b",
                                 "Under-dispersed" = "#2166ac"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL, y = "Dispersion index (variance / mean)",
       title = "Are annual exceedance counts Poisson?",
       subtitle = "Counts above u per qualifying water year; 95% chi-square bounds per site") +
  theme_poster
save_fig("13_dispersion_index.png", p, 7, 5)
saveRDS(disp, file.path(fig_dir, "dispersion_index.rds"))

# ============================================================================
# SENSITIVITY (S3, S4)
# ============================================================================
sens_file <- "_data/sensitivity_ratio_summary.rds"
if (file.exists(sens_file)) {
  sens <- readRDS(sens_file) %>% filter(method == MAIN_METHOD)
  sens_plot <- function(fam, title, subtitle) {
    d <- sens %>% filter(family %in% c("baseline", fam)) %>%
      mutate(setting = fct_relevel(factor(setting), "Baseline"))
    ggplot(d, aes(return_period, median, colour = setting)) +
      indicative_layer() +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
      geom_line(aes(linewidth = setting == "Baseline")) + geom_point(size = 1.8) +
      scale_linewidth_manual(values = c(`TRUE` = 1.4, `FALSE` = 0.8), guide = "none") +
      scale_colour_manual(values = c("Baseline" = "black", setNames(
        scales::brewer_pal(palette = "Set2")(max(3, n_distinct(d$setting) - 1))[seq_len(n_distinct(d$setting) - 1)],
        setdiff(levels(d$setting), "Baseline"))), name = NULL) +
      rp_axis +
      labs(y = "Median mixture / single sample", title = title, subtitle = subtitle) +
      theme_poster
  }
  save_fig("14_sensitivity_threshold.png",
           sens_plot("S3", "Sensitivity to the fitting threshold",
                     "Interim rate-based threshold (baseline: ~1 exceedance per year); EQD replaces it later"), 7, 5)
  save_fig("15_sensitivity_rain_fraction.png",
           sens_plot("S4", "Sensitivity to the rain-fraction cutoffs",
                     "Baseline: alpha_r = 0.99, alpha_s = network median of mixed events"), 7, 5)
} else {
  message("  (skipping 14-15: run _documentation/poster_sensitivity.R first)")
}

message("Done.")
