#!/usr/bin/env Rscript
# ============================================================================
# Poster showcase-site panels, styled after the Hydrology Days 2026 talk
# ============================================================================
# Paired panels per candidate site:
#   A  Fixed (depth threshold) vs Adaptive (3-class rain fraction)
#   B  Adaptive (3-class) vs Adaptive (2-class)
# Style follows the talk: Set1 type colours, dashed per-type curves, black
# mixture, grey single sample, dashed line at 100 yr. New Poisson–GPD fits;
# annual return periods to 500 yr (plan D-F, revised 2026-11-19).
#
# Candidates were screened (2026-11-19) for: unregulated; all three types
# fitted (3-class); Q100 mixture/single ratio > 1 under every method and under
# the S3 threshold settings; GPD shape xi < 0.55 for every type.
#
# Run from the repo root:  Rscript _documentation/poster_showcase.R [site ...]
# Output: _documentation/poster_figures/showcase/
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
source("functions.R")

args <- commandArgs(trailingOnly = TRUE)
SITES <- if (length(args)) args else c("09113980", "09115500", "09112500", "09342500")

out_dir <- "_documentation/poster_figures/showcase"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

curves <- readRDS("_data/ffa_curves_pgpd.rds")
params <- readRDS("_data/ffa_params_pgpd.rds")
qnt    <- readRDS("_data/ffa_quantiles_pgpd.rds")
meta   <- readRDS("_data/site_metadata.rds")
qual   <- readRDS("_data/site_qualifying_years.rds") %>% filter(selected, qualifies)
peaks  <- list(
  depth_thresh       = readRDS("_data/pot_peaks_classified_depth_thresh.rds"),
  frac_thresh_3class = readRDS("_data/pot_peaks_classified_frac_thresh_3class.rds"))
peaks$frac_thresh_2class <- remap_2class(peaks$frac_thresh_3class)

method_titles <- c(depth_thresh = "Fixed", frac_thresh_3class = "Adaptive (3-class)",
                   frac_thresh_2class = "Adaptive (2-class)")
labels <- c(rainfall = "Rainfall", snowmelt = "Snowmelt", rain_on_snow = "Rain-on-Snow",
            snow_influenced = "Non-Rainfall", mixture = "Mixture", single_sample = "Single Sample")
colors <- c("Rainfall" = "#E41A1C", "Snowmelt" = "#377EB8", "Rain-on-Snow" = "#4DAF4A",
            "Non-Rainfall" = "#984EA3", "Mixture" = "black", "Single Sample" = "grey55")
ltys   <- c("Rainfall" = "dashed", "Snowmelt" = "dashed", "Rain-on-Snow" = "dashed",
            "Non-Rainfall" = "dashed", "Mixture" = "solid", "Single Sample" = "solid")
RP_BREAKS <- c(2, 5, 10, 25, 50, 100, 200, 500)

panel <- function(s, m, ylim, brks) {
  label_col <- if (m == "frac_thresh_2class") "flood_type" else "flood_type_raw"
  pool <- peaks[[m]] %>%
    filter(site_no == s, fit_eligible, is.finite(peak_flow), peak_flow > 0) %>%
    semi_join(qual, by = c("site_no", "water_year"))
  prm <- params %>% filter(site_no == s, method == m)
  pts <- prm %>% filter(curve_type != "single_sample") %>%
    pmap_dfr(function(curve_type, u, n_years, ...)
      empirical_annual_rp(pool$peak_flow[as.character(pool[[label_col]]) == curve_type], u, n_years) %>%
        mutate(lab = labels[curve_type]))
  ss  <- prm %>% filter(curve_type == "single_sample")
  ss_pts <- empirical_annual_rp(pool$peak_flow, ss$u, ss$n_years)
  cv <- curves %>% filter(site_no == s, method == m, is.finite(discharge_cfs)) %>%
    mutate(lab = labels[curve_type])
  ggplot() +
    geom_vline(xintercept = 100, linetype = "dashed") +
    geom_point(data = ss_pts, aes(ri, peak_flow, shape = "Observed (single sample)"),
               colour = "grey45", size = 3.6, stroke = 0.9) +
    geom_point(data = pts, aes(ri, peak_flow, colour = lab), size = 2.4, alpha = 0.75) +
    geom_line(data = cv, aes(ri, discharge_cfs, colour = lab, linetype = lab), linewidth = 1.4) +
    scale_colour_manual(values = colors[brks], breaks = brks, limits = brks, name = NULL) +
    scale_linetype_manual(values = ltys[brks], breaks = brks, limits = brks, name = NULL) +
    scale_shape_manual(values = c("Observed (single sample)" = 1), name = NULL) +
    scale_x_log10(breaks = RP_BREAKS, limits = c(1.05, 500), name = "Annual Return Period (years)") +
    scale_y_log10(labels = scales::label_comma(), limits = ylim, name = "Peak Discharge (cfs)") +
    labs(title = method_titles[[m]]) +
    theme_bw(base_size = 18) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          plot.title = element_text(face = "bold", hjust = 0.5, size = 24),
          axis.title = element_text(size = 16), axis.text = element_text(size = 14),
          legend.text = element_text(size = 16), legend.key.width = unit(1.6, "cm")) +
    guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
}

pair <- function(s, m1, m2, tag) {
  nm <- str_remove(str_to_title(meta$station_name[meta$site_no == s]), "\\.$")
  rng <- curves %>% filter(site_no == s, method %in% c(m1, m2), is.finite(discharge_cfs)) %>%
    pull(discharge_cfs) %>% range()
  ratio <- function(m) {
    q <- qnt %>% filter(site_no == s, method == m, return_period == 100)
    q$discharge_cfs[q$curve_type == "mixture"] / q$discharge_cfs[q$curve_type == "single_sample"]
  }
  brks <- curves %>% filter(site_no == s, method %in% c(m1, m2)) %>% distinct(curve_type) %>%
    mutate(lab = labels[curve_type]) %>% pull(lab)
  brks <- names(colors)[names(colors) %in% brks]
  p <- (panel(s, m1, rng, brks) | panel(s, m2, rng, brks)) + plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("Station %s — Poisson–GPD Flood Frequency Analysis", s),
      subtitle = sprintf("%s. %s vs. %s. Q100 mixture / single sample: %.2f vs. %.2f",
                         nm, method_titles[[m1]], method_titles[[m2]], ratio(m1), ratio(m2)),
      theme = theme(plot.title = element_text(face = "bold", size = 24),
                    plot.subtitle = element_text(colour = "grey40", size = 18))) &
    theme(legend.position = "bottom")
  f <- file.path(out_dir, sprintf("%s_%s.png", s, tag))
  # Exported at its placed size on the poster (two-column block, 20.15 in), 300 dpi,
  # with the AGU 2025 text sizes (title 24, subtitle 18, axes 16/14, legend 16).
  ggsave(f, p, width = 20.15, height = 10, dpi = 300, bg = "white")
  message("  ", f)
}

for (s in SITES) {
  pair(s, "depth_thresh", "frac_thresh_3class", "A_fixed_vs_adaptive3")
  pair(s, "frac_thresh_3class", "frac_thresh_2class", "B_adaptive3_vs_adaptive2")
}
