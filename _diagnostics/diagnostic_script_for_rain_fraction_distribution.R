#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: rain_frac Distribution in the Middle Bin
# ============================================================================
# Examines the distribution of rain_frac for events that fall between the
# rainfall and snowmelt end-member thresholds (0.10 < rain_frac < 0.95).
# These are the events that would be classified as ROS under the simplified
# 3-class scheme, or split into rainfall-dominated ROS / snowmelt-dominated
# ROS under a 4-class scheme.
#
# Produces three figures:
#   1. Density of rain_frac across all middle-bin events (pooled)
#   2. Same, faceted by elevation band
#   3. Per-site sample sizes in the middle bin
#
# INPUTS:
#   - _data/pot_peaks_classified_adaptive.rds : output of script 6b
#   - _data/site_metadata.rds                : site attributes (elevation)
#
# OUTPUT:
#   - _diagnostics/rain_frac_dist/rain_frac_pooled_density.png
#   - _diagnostics/rain_frac_dist/rain_frac_density_by_elevation.png
#   - _diagnostics/rain_frac_dist/rain_frac_middle_bin_n_per_site.png
# ============================================================================

library(tidyverse)

# ============================================================================
# CONFIGURATION
# ============================================================================
classified_file <- "_data/pot_peaks_classified_adaptive.rds"
meta_file       <- "_data/site_metadata.rds"
out_dir         <- "_diagnostics/rain_frac_dist"

RAIN_EM_THRESH  <- 0.95
SNOW_EM_THRESH  <- 0.10

# Elevation breaks (m) for faceting
ELEV_BREAKS     <- c(0, 1500, 2000, 2500, 3000, Inf)
ELEV_LABELS     <- c("< 1500 m", "1500–2000 m", "2000–2500 m",
                     "2500–3000 m", "> 3000 m")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD AND PREPARE DATA
# ============================================================================
message("Loading data...")
classified <- readRDS(classified_file)
meta       <- readRDS(meta_file) %>%
  select(site_no, elev_mean_m)

# Use flood_type_raw (pre-demotion) so sample-size demotions don't suppress
# events that were legitimately classified. Join elevation.
df <- classified %>%
  left_join(meta, by = "site_no") %>%
  filter(data_sufficient) %>%
  mutate(
    elev_band = cut(elev_mean_m,
                    breaks = ELEV_BREAKS,
                    labels = ELEV_LABELS,
                    right  = FALSE)
  )

# Middle bin: events not captured by either end-member threshold
middle <- df %>%
  filter(rain_frac > SNOW_EM_THRESH, rain_frac < RAIN_EM_THRESH)

message(sprintf("Total data-sufficient peaks : %d", nrow(df)))
message(sprintf("Middle-bin peaks            : %d (%.1f%%)",
                nrow(middle), 100 * nrow(middle) / nrow(df)))

# ============================================================================
# FIGURE 1: Pooled density of rain_frac in middle bin
# ============================================================================
message("Plotting Figure 1: pooled density...")

p1 <- ggplot(middle, aes(x = rain_frac)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.025, fill = "#7b3294", alpha = 0.6, color = "white") +
  geom_density(color = "#1a1a1a", linewidth = 0.8) +
  geom_vline(xintercept = 0.50, linetype = "dashed", color = "firebrick", linewidth = 0.7) +
  annotate("text", x = 0.52, y = Inf, label = "0.50", vjust = 1.5,
           hjust = 0, color = "firebrick", size = 3.5) +
  scale_x_continuous(limits = c(SNOW_EM_THRESH, RAIN_EM_THRESH),
                     breaks = seq(0.1, 0.95, by = 0.1)) +
  labs(
    title    = "rain_frac distribution — middle bin (pooled)",
    subtitle = sprintf("Events with %.2f < rain_frac < %.2f  |  n = %d",
                       SNOW_EM_THRESH, RAIN_EM_THRESH, nrow(middle)),
    x        = "rain_frac",
    y        = "Density"
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(out_dir, "rain_frac_pooled_density.png"),
       p1, width = 8, height = 5, dpi = 150)
message("  Saved rain_frac_pooled_density.png")

# ============================================================================
# FIGURE 2: Density faceted by elevation band
# ============================================================================
message("Plotting Figure 2: density by elevation band...")

elev_n <- middle %>%
  dplyr::count(elev_band) %>%
  mutate(label = sprintf("%s\nn = %d", elev_band, n))

middle_elev <- middle %>%
  left_join(elev_n %>% select(elev_band, label), by = "elev_band")

p2 <- ggplot(middle_elev, aes(x = rain_frac)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.05, fill = "#7b3294", alpha = 0.6, color = "white") +
  geom_density(color = "#1a1a1a", linewidth = 0.7) +
  geom_vline(xintercept = 0.50, linetype = "dashed",
             color = "firebrick", linewidth = 0.6) +
  scale_x_continuous(limits = c(SNOW_EM_THRESH, RAIN_EM_THRESH),
                     breaks = seq(0.1, 0.9, by = 0.2)) +
  facet_wrap(~ label, ncol = 3) +
  labs(
    title    = "rain_frac distribution by elevation band — middle bin",
    subtitle = sprintf("Events with %.2f < rain_frac < %.2f",
                       SNOW_EM_THRESH, RAIN_EM_THRESH),
    x        = "rain_frac",
    y        = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "#f0f0f0"))

ggsave(file.path(out_dir, "rain_frac_density_by_elevation.png"),
       p2, width = 11, height = 7, dpi = 150)
message("  Saved rain_frac_density_by_elevation.png")

# ============================================================================
# FIGURE 3: Per-site sample sizes in the middle bin
# ============================================================================
message("Plotting Figure 3: per-site n in middle bin...")

site_n <- middle %>%
  dplyr::count(site_no) %>%
  left_join(meta, by = "site_no") %>%
  arrange(elev_mean_m) %>%
  mutate(site_no = factor(site_no, levels = site_no))

p3 <- ggplot(site_n, aes(x = n, y = elev_mean_m)) +
  geom_point(alpha = 0.5, size = 1.5, color = "#7b3294") +
  geom_vline(xintercept = 20, linetype = "dashed",
             color = "firebrick", linewidth = 0.7) +
  annotate("text", x = 21, y = Inf, label = "MIN_TYPE_SAMPLE = 20",
           vjust = 1.5, hjust = 0, color = "firebrick", size = 3.5) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  labs(
    title    = "Middle-bin peak count per site vs. mean basin elevation",
    subtitle = "Each point is one site; dashed line = MIN_TYPE_SAMPLE threshold",
    x        = "Number of middle-bin peaks",
    y        = "Mean basin elevation (m)"
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(out_dir, "rain_frac_middle_bin_n_per_site.png"),
       p3, width = 8, height = 6, dpi = 150)
message("  Saved rain_frac_middle_bin_n_per_site.png")

# ============================================================================
# CONSOLE SUMMARY
# ============================================================================
message("\n=== MIDDLE BIN SUMMARY ===")
message(sprintf("  Median rain_frac : %.3f", median(middle$rain_frac, na.rm = TRUE)))
message(sprintf("  Mean rain_frac   : %.3f", mean(middle$rain_frac,   na.rm = TRUE)))
message(sprintf("  SD rain_frac     : %.3f", sd(middle$rain_frac,     na.rm = TRUE)))

q_summary <- quantile(middle$rain_frac, probs = c(0.10, 0.25, 0.50, 0.75, 0.90),
                      na.rm = TRUE)
message("  Quantiles:")
walk2(names(q_summary), q_summary,
      ~ message(sprintf("    %s : %.3f", .x, .y)))

sites_above_min <- sum(site_n$n >= 20)
message(sprintf("\n  Sites with >= 20 middle-bin peaks : %d / %d (%.0f%%)",
                sites_above_min, nrow(site_n),
                100 * sites_above_min / nrow(site_n)))

message(sprintf("\nDone. Outputs in %s/", out_dir))