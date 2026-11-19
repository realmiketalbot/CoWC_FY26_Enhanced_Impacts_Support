#!/usr/bin/env Rscript
# ============================================================================
# Old vs new FFA results: figures for INTERNAL review (not for the poster)
# ============================================================================
# Compares the initial results in the Preliminary Proposal (tag
# initial-results-prelim-proposal) with the M1 Poisson–GPD results.
#
# What changed between them (see _planning/agu2026_interim_plan.md):
#   OLD  206 sites (calendar-year approval share). Per type, GPD fit by
#        L-moments to the top-n peaks, with pargpa() estimating a location
#        that was then discarded. Mixture S_mix = 1 - prod_i F_i on per-EVENT
#        CDFs, read as RI = 1 / S_mix: a peaks-over-threshold recurrence
#        interval, not an annual return period. LP3 as an alternative.
#   NEW  215 sites (>= 20 qualifying May–Sep water years). Poisson–GPD per type
#        above an interim rate-based threshold (lower bound fixed at 0);
#        F_mix = prod_i exp(-lambda_i (1 - H_i)); ANNUAL T = 1 / (1 - F).
# So the two x-axes are not the same quantity: "old RI" vs "annual T".
# Figures use the rain-fraction (3-class) method unless noted.
#
# Run from the repo root:  Rscript _documentation/advisor_old_vs_new.R
# Output: _documentation/advisor_old_vs_new/*.png
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(patchwork)
})
source("functions.R")

TAG     <- "initial-results-prelim-proposal"
METHOD  <- "frac_thresh_3class"
fig_dir <- "_documentation/advisor_old_vs_new"
dir.create(fig_dir, showWarnings = FALSE)

# ---- Data -------------------------------------------------------------------
read_tag <- function(path) {
  tmp <- tempfile(fileext = ".rds")
  status <- system2("git", c("show", sprintf("%s:%s", TAG, path)), stdout = tmp)
  if (status != 0) stop("Could not read ", path, " at tag ", TAG)
  readRDS(tmp)
}
old_q   <- read_tag("_data/ffa_quantiles_gpd.rds")
old_lp3 <- read_tag("_data/ffa_quantiles_lp3.rds")
old_c   <- read_tag("_data/ffa_curves_gpd.rds")
old_sites <- unique(read_tag("_data/site_metadata.rds")$site_no)

new_q <- readRDS("_data/ffa_quantiles_pgpd.rds")
new_c <- readRDS("_data/ffa_curves_pgpd.rds")
new_p <- readRDS("_data/ffa_params_pgpd.rds")
meta  <- readRDS("_data/site_metadata.rds")
stations <- readRDS("_data/station_spatial.rds")
peaks <- readRDS("_data/pot_peaks_classified_frac_thresh_3class.rds")
qual  <- readRDS("_data/site_qualifying_years.rds") %>% filter(selected, qualifies)
pool  <- peaks %>% filter(fit_eligible, is.finite(peak_flow), peak_flow > 0) %>%
  semi_join(qual, by = c("site_no", "water_year"))

RP_BREAKS <- c(2, 5, 10, 25, 50, 100, 200, 500)
version_colors <- c("Old: GPD (top-n, PoT RI)" = "#d95f02",
                    "Old: LP3 (top-n, PoT RI)" = "#fdb863",
                    "New: Poisson–GPD (annual T)" = "#1b7837")
type_labels <- c(rainfall = "Rainfall", snowmelt = "Snowmelt", rain_on_snow = "Rain-on-Snow")
type_colors <- c("Rainfall" = "#b2182b", "Snowmelt" = "#2166ac", "Rain-on-Snow" = "#9970ab")

theme_review <- theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.subtitle = element_text(size = 10.5, face = "italic", colour = "grey30"),
        plot.caption = element_text(size = 9, colour = "grey40", hjust = 0))
axis_note <- "Old x-axis: peaks-over-threshold recurrence interval. New x-axis: annual return period."

save_fig <- function(name, p, w = 7, h = 5) {
  ggsave(file.path(fig_dir, name), p, width = w, height = h, units = "in", dpi = 200)
  message("  ", name)
}

ratio_tbl <- function(q, version) {
  q %>% filter(method == METHOD, curve_type %in% c("mixture", "single_sample")) %>%
    select(site_no, return_period, curve_type, discharge_cfs) %>%
    pivot_wider(names_from = curve_type, values_from = discharge_cfs) %>%
    filter(is.finite(mixture), is.finite(single_sample), mixture > 0, single_sample > 0) %>%
    mutate(ratio = mixture / single_sample, version = version)
}
ratios <- bind_rows(ratio_tbl(old_q,   names(version_colors)[1]),
                    ratio_tbl(old_lp3, names(version_colors)[2]),
                    ratio_tbl(new_q,   names(version_colors)[3])) %>%
  mutate(version = factor(version, levels = names(version_colors)))

common <- intersect(old_sites, unique(new_q$site_no))
message(sprintf("Sites: old %d, new %d, common %d", length(old_sites),
                n_distinct(new_q$site_no), length(common)))
message("Writing figures to ", fig_dir)

# ============================================================================
# 01 Site set: kept, dropped, added
# ============================================================================
gages <- readRDS("_data/station_spatial.rds")
all_pts <- meta %>% select(site_no, lon_deg, lat_deg) %>%
  bind_rows(read_tag("_data/site_metadata.rds") %>% select(site_no, lon_deg, lat_deg)) %>%
  distinct(site_no, .keep_all = TRUE) %>%
  filter(site_no %in% union(old_sites, gages$site_no)) %>%
  mutate(status = case_when(site_no %in% old_sites & site_no %in% gages$site_no ~ "Kept",
                            site_no %in% old_sites ~ "Dropped (new screen)",
                            TRUE ~ "Added (new screen)")) %>%
  st_as_sf(coords = c("lon_deg", "lat_deg"), crs = 4326)
co <- if (requireNamespace("tigris", quietly = TRUE)) {
  tigris::states(cb = TRUE, class = "sf", progress_bar = FALSE) %>% filter(STUSPS == "CO")
} else st_as_sf(maps::map("state", "colorado", plot = FALSE, fill = TRUE))
co <- st_transform(co, 4326)
cnt <- count(st_drop_geometry(all_pts), status)
p <- ggplot() +
  geom_sf(data = co, fill = "white", colour = "grey20") +
  geom_sf(data = all_pts %>% arrange(status == "Kept"), aes(colour = status, shape = status), size = 2.3) +
  scale_colour_manual(values = c("Kept" = "grey60", "Dropped (new screen)" = "#d95f02",
                                 "Added (new screen)" = "#1b7837"), name = NULL) +
  scale_shape_manual(values = c("Kept" = 16, "Dropped (new screen)" = 4, "Added (new screen)" = 17), name = NULL) +
  labs(title = "Station set: old (206) vs new (215)",
       subtitle = paste(sprintf("%s: %d", cnt$status, cnt$n), collapse = "; ")) +
  theme_void(base_size = 13) + theme(legend.position = "bottom",
                                     plot.subtitle = element_text(size = 10.5, face = "italic"))
save_fig("01_site_set_changes.png", p, 7, 5.5)

# ============================================================================
# 02 Mixture / single-sample ratio by return period, old vs new
# ============================================================================
rs <- ratios %>% group_by(version, return_period) %>%
  summarise(median = median(ratio), q25 = quantile(ratio, .25), q75 = quantile(ratio, .75),
            share_ge1 = mean(ratio >= 1), n = n(), .groups = "drop")
p <- ggplot(rs, aes(return_period, median, colour = version, fill = version)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_x_log10(breaks = RP_BREAKS, name = "Return period / recurrence interval (years)") +
  scale_colour_manual(values = version_colors, name = NULL, aesthetics = c("colour", "fill")) +
  guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2)) +
  labs(y = "Mixture / single-sample discharge", title = "Mixture / single-sample ratio: old vs new",
       subtitle = "Median across sites with interquartile range; rain fraction (3-class)",
       caption = axis_note) +
  theme_review
save_fig("02_ratio_by_rp_old_vs_new.png", p, 7.5, 5.5)

# ============================================================================
# 03 Share of sites where mixture >= single sample, old vs new
# ============================================================================
p <- ggplot(rs, aes(return_period, share_ge1, colour = version)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_x_log10(breaks = RP_BREAKS, name = "Return period / recurrence interval (years)") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_colour_manual(values = version_colors, name = NULL) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(y = "Sites where mixture >= single sample",
       title = "How often the mixture exceeds the single sample", caption = axis_note) +
  theme_review
save_fig("03_share_ge1_old_vs_new.png", p, 7.5, 5.5)

# ============================================================================
# 04 Site-by-site Q100 ratio, old vs new (common sites)
# ============================================================================
r100 <- ratios %>% filter(return_period == 100, site_no %in% common,
                          version != names(version_colors)[2]) %>%
  mutate(v = if_else(str_starts(as.character(version), "Old"), "old", "new")) %>%
  select(site_no, v, ratio) %>% pivot_wider(names_from = v, values_from = ratio) %>%
  filter(is.finite(old), is.finite(new)) %>%
  left_join(meta %>% select(site_no, elev_mean_m), by = "site_no")
quad <- r100 %>% summarise(both_up = mean(old >= 1 & new >= 1), up_down = mean(old >= 1 & new < 1),
                           down_up = mean(old < 1 & new >= 1), both_down = mean(old < 1 & new < 1))
p <- ggplot(r100, aes(old, new, colour = elev_mean_m)) +
  geom_hline(yintercept = 1, colour = "grey60") + geom_vline(xintercept = 1, colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_viridis_c(option = "C", name = "Mean basin\nelevation (m)") +
  labs(x = "Old Q100 ratio (GPD, top-n)", y = "New Q100 ratio (Poisson–GPD)",
       title = "Site-level Q100 mixture / single-sample ratio",
       subtitle = sprintf("%d common sites. Both >= 1: %.0f%%; old >= 1 but new < 1: %.0f%%; old < 1 but new >= 1: %.0f%%",
                          nrow(r100), 100 * quad$both_up, 100 * quad$up_down, 100 * quad$down_up)) +
  theme_review + theme(legend.position = "right")
save_fig("04_q100_ratio_site_scatter.png", p, 7.5, 5.5)

# ============================================================================
# 05 Q100 discharge, old vs new, for single sample and mixture
# ============================================================================
q100 <- bind_rows(old_q %>% mutate(v = "old"), new_q %>% mutate(v = "new")) %>%
  filter(method == METHOD, return_period == 100, curve_type %in% c("mixture", "single_sample"),
         site_no %in% common) %>%
  select(site_no, curve_type, v, discharge_cfs) %>%
  pivot_wider(names_from = v, values_from = discharge_cfs) %>%
  filter(is.finite(old), is.finite(new), old > 0, new > 0) %>%
  mutate(curve_type = recode(curve_type, mixture = "Mixture", single_sample = "Single sample"))
lab05 <- q100 %>% group_by(curve_type) %>%
  summarise(txt = sprintf("median new/old = %.2f\n(n = %d)", median(new / old), n()), .groups = "drop") %>%
  mutate(x = min(q100$old), y = max(q100$new))
p <- ggplot(q100, aes(old, new)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.6, size = 1.8, colour = "#1b7837") +
  geom_text(data = lab05, aes(x = x, y = y, label = txt), hjust = 0, vjust = 1, size = 3.6,
            inherit.aes = FALSE) +
  facet_wrap(~ curve_type) +
  scale_x_log10(labels = scales::label_comma(), name = "Old Q100 (cfs)") +
  scale_y_log10(labels = scales::label_comma(), name = "New Q100 (cfs)") +
  coord_cartesian(clip = "off") +
  labs(title = "100-year discharge: old vs new", caption = axis_note) +
  theme_review
save_fig("05_q100_discharge_old_vs_new.png", p, 9, 5)

# ============================================================================
# 06 New / old discharge by return period (median, IQR)
# ============================================================================
nvo <- bind_rows(old_q %>% mutate(v = "old"), new_q %>% mutate(v = "new")) %>%
  filter(method == METHOD, curve_type %in% c("mixture", "single_sample"), site_no %in% common) %>%
  select(site_no, curve_type, return_period, v, discharge_cfs) %>%
  pivot_wider(names_from = v, values_from = discharge_cfs) %>%
  filter(is.finite(old), is.finite(new), old > 0, new > 0) %>%
  mutate(curve_type = recode(curve_type, mixture = "Mixture", single_sample = "Single sample")) %>%
  group_by(curve_type, return_period) %>%
  summarise(median = median(new / old), q25 = quantile(new / old, .25), q75 = quantile(new / old, .75),
            .groups = "drop")
p <- ggplot(nvo, aes(return_period, median, colour = curve_type, fill = curve_type)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_x_log10(breaks = RP_BREAKS, name = "Return period / recurrence interval (years)") +
  scale_colour_manual(values = c("Mixture" = "black", "Single sample" = "grey50"), name = NULL,
                      aesthetics = c("colour", "fill")) +
  labs(y = "New / old discharge", title = "How much each estimate moved",
       subtitle = "Median across common sites with IQR; rain fraction (3-class)", caption = axis_note) +
  theme_review
save_fig("06_new_over_old_by_rp.png", p, 7.5, 5.5)

# ============================================================================
# 07 Dominant flood type at Q100: old vs new
# ============================================================================
dom <- function(q) q %>%
  filter(method == METHOD, return_period == 100, curve_type %in% names(type_labels),
         is.finite(discharge_cfs), discharge_cfs > 0) %>%
  group_by(site_no) %>% slice_max(discharge_cfs, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(site_no, type = type_labels[curve_type])
dd <- inner_join(dom(old_q) %>% rename(old = type), dom(new_q) %>% rename(new = type), by = "site_no") %>%
  count(old, new) %>%
  mutate(old = factor(old, levels = type_labels), new = factor(new, levels = type_labels))
p <- ggplot(dd, aes(new, old, fill = n)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n), size = 5) +
  scale_fill_gradient(low = "grey95", high = "#1b7837", guide = "none") +
  labs(x = "New: dominant type at Q100", y = "Old: dominant type at Q100",
       title = "Which flood type has the largest 100-year flood?",
       subtitle = sprintf("%d common sites with type fits in both; %.0f%% unchanged",
                          sum(dd$n), 100 * sum(dd$n[dd$old == dd$new]) / sum(dd$n))) +
  theme_review
save_fig("07_dominant_type_q100_old_vs_new.png", p, 6.5, 5.5)

# ============================================================================
# 08 Example sites: old vs new curves
# ============================================================================
examples <- c("09081600", "09146020", "09025000")
examples <- examples[examples %in% common]
plot_example <- function(s) {
  cv <- bind_rows(
    old_c %>% filter(site_no == s, method == METHOD) %>% mutate(v = "Old (PoT RI)"),
    new_c %>% filter(site_no == s, method == METHOD) %>% mutate(v = "New (annual T)")) %>%
    filter(curve_type %in% c("mixture", "single_sample"), is.finite(discharge_cfs), discharge_cfs > 0) %>%
    mutate(curve = recode(curve_type, mixture = "Mixture", single_sample = "Single sample"),
           v = factor(v, levels = c("Old (PoT RI)", "New (annual T)")))
  u_ss <- new_p %>% filter(site_no == s, method == METHOD, curve_type == "single_sample")
  pts <- empirical_annual_rp(pool$peak_flow[pool$site_no == s], u_ss$u, u_ss$n_years)
  nm  <- str_trunc(str_to_title(meta$station_name[meta$site_no == s]), 40)
  ggplot() +
    geom_point(data = pts, aes(ri, peak_flow), colour = "grey55", size = 1.2, alpha = 0.6) +
    geom_line(data = cv, aes(ri, discharge_cfs, colour = curve, linetype = v), linewidth = 0.9) +
    scale_x_log10(breaks = RP_BREAKS, limits = c(1.1, 500), name = "Return period / RI (years)") +
    scale_y_log10(labels = scales::label_comma(), name = "Peak discharge (cfs)") +
    scale_colour_manual(values = c("Mixture" = "black", "Single sample" = "#e7298a"), name = NULL) +
    scale_linetype_manual(values = c("Old (PoT RI)" = "dotted", "New (annual T)" = "solid"), name = NULL) +
    labs(title = sprintf("%s (%s)", nm, s)) +
    theme_review + theme(plot.title = element_text(size = 11))
}
p <- wrap_plots(lapply(examples, plot_example), nrow = 1, guides = "collect") +
  plot_annotation(caption = paste("Points: observed peaks above the new single-sample threshold, annual plotting positions.",
                                  axis_note),
                  theme = theme(plot.caption = element_text(size = 9, colour = "grey40", hjust = 0))) &
  theme(legend.position = "bottom")
save_fig("08_example_sites_old_vs_new.png", p, 15, 5.5)

# ============================================================================
# 09 Side-by-side scatter: mixture vs single sample, old vs new, at 2/10/100
# ============================================================================
sc <- bind_rows(
  ratio_tbl(old_q, "Old: GPD (top-n, PoT RI)"),
  ratio_tbl(new_q, "New: Poisson–GPD (annual T)")) %>%
  filter(site_no %in% common, return_period %in% c(2, 10, 100)) %>%
  mutate(version = factor(version, levels = c("Old: GPD (top-n, PoT RI)", "New: Poisson–GPD (annual T)")),
         rp = factor(sprintf("%d-year", as.integer(return_period)), levels = c("2-year", "10-year", "100-year")))
lim <- range(c(sc$mixture, sc$single_sample))
lab09 <- sc %>% group_by(version, rp) %>%
  summarise(txt = sprintf("mixture > single: %.0f%%\nmedian ratio: %.2f", 100 * mean(ratio > 1), median(ratio)),
            .groups = "drop")
p <- ggplot(sc, aes(single_sample, mixture)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(aes(colour = version), alpha = 0.6, size = 1.6) +
  geom_text(data = lab09, aes(x = lim[1], y = lim[2], label = txt), hjust = 0, vjust = 1, size = 3.4) +
  facet_grid(version ~ rp) +
  scale_x_log10(labels = scales::label_comma(), limits = lim, name = "Single-sample discharge (cfs)") +
  scale_y_log10(labels = scales::label_comma(), limits = lim, name = "Mixture discharge (cfs)") +
  scale_colour_manual(values = version_colors, guide = "none") +
  coord_equal() +
  labs(title = "Mixture vs single sample, site by site: old vs new",
       subtitle = sprintf("%d common sites; Adaptive (3-class). Dashed: 1:1", length(common)),
       caption = axis_note) +
  theme_review + theme(strip.text = element_text(face = "bold"))
save_fig("09_scatter_grid_old_vs_new.png", p, 11, 8)

# 09b The same at Q100 only, two panels side by side (proposal-figure style)
p <- sc %>% filter(rp == "100-year") %>%
  ggplot(aes(single_sample, mixture)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(aes(colour = version), alpha = 0.7, size = 2) +
  geom_text(data = lab09 %>% filter(rp == "100-year"), aes(x = lim[1], y = lim[2], label = txt),
            hjust = 0, vjust = 1, size = 3.8) +
  facet_wrap(~ version) +
  scale_x_log10(labels = scales::label_comma(), limits = lim, name = "Single-sample Q100 (cfs)") +
  scale_y_log10(labels = scales::label_comma(), limits = lim, name = "Mixture Q100 (cfs)") +
  scale_colour_manual(values = version_colors, guide = "none") +
  coord_equal() +
  labs(title = "Mixture vs single-sample Q100: old vs new", caption = axis_note) +
  theme_review + theme(strip.text = element_text(face = "bold"))
save_fig("09b_scatter_q100_old_vs_new.png", p, 10, 5.5)

message("Done.")
