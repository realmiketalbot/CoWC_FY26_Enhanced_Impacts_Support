#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: ROS Feature Space Structure
# ============================================================================
# Examines whether the snowmelt/ROS boundary reflects a genuine discontinuity
# in the data or an arbitrary cut through a continuous gradient.
#
# Specifically investigates:
#   1. Marginal distribution of ros_flux_frac (population and per-site)
#   2. Joint distribution of ros_flux_frac x rf_long
#   3. Joint distribution of ros_flux_frac x rf_short
#   4. All three features in the "interior" — events that are neither pure
#      rainfall nor pure snowmelt end-members (i.e., the space where the
#      snowmelt/ROS/mixed boundary is actually being drawn)
#   5. Seasonal structure of interior events (day of water year)
#   6. Within-site bimodality of ros_flux_frac using Hartigan's dip test
#
# The key question: does ros_flux_frac show bimodality (consistent with
# a real process discontinuity) or unimodality/uniformity (consistent
# with an arbitrary threshold on a gradient)?
#
# INPUTS:
#   - _data/pot_peaks_classified_adaptive.rds : classified peaks with
#     rf_short, rf_long, ros_flux_frac columns
#   - _data/site_metadata.rds                 : regulated flag, basin attrs
#
# OUTPUTS (written to _diagnostics/ros_feature_space/):
#   - 01_ros_flux_frac_marginal.png    : population-level density + histogram
#   - 02_ros_flux_frac_by_type.png     : marginal by current flood_type
#   - 03_joint_ros_vs_rf_long.png      : 2D joint density, interior events
#   - 04_joint_ros_vs_rf_short.png     : 2D joint density, interior events
#   - 05_feature_triangle.png          : rf_short x rf_long x ros_flux_frac
#                                        (color = ros_flux_frac)
#   - 06_seasonal_interior.png         : day-of-water-year distribution,
#                                        interior events split by current type
#   - 07_dip_test_summary.png          : site-level Hartigan dip p-values
#   - 08_ros_frac_site_densities.png   : small-multiple per-site densities
#                                        (unregulated sites only, n >= 20)
#   - diag_ros_feature_space.txt       : console summary of dip test results
# ============================================================================

library(tidyverse)
library(patchwork)
library(diptest)   # Hartigan's dip test for unimodality

# ============================================================================
# CONFIGURATION
# ============================================================================
classified_file <- "_data/pot_peaks_classified_adaptive.rds"
meta_file       <- "_data/site_metadata.rds"
out_dir         <- "_diagnostics/ros_feature_space"

# Thresholds from 5b — used only for labeling and for defining "interior"
RAIN_EM_THRESH   <- 0.95
SNOW_EM_THRESH   <- 0.05
ROS_THRESH       <- 0.50

# Minimum peaks per site for within-site density plots
MIN_SITE_N       <- 20

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading classified peaks...")
peaks <- readRDS(classified_file)
meta  <- readRDS(meta_file)

message(sprintf("Loaded %d peaks across %d sites", nrow(peaks), n_distinct(peaks$site_no)))

# Attach regulated flag if not already present
if (!"regulated" %in% names(peaks)) {
  peaks <- peaks %>%
    left_join(meta %>% select(site_no, regulated), by = "site_no")
}

# Water year day-of-year: water year starts Oct 1
# peak_date should already exist; compute wy_doy
peaks <- peaks %>%
  mutate(
    peak_date = as.Date(peak_date),
    month     = lubridate::month(peak_date),
    # Water year day: Oct 1 = day 1
    wy_doy    = as.integer(
      peak_date - as.Date(paste0(
        if_else(month >= 10,
                lubridate::year(peak_date),
                lubridate::year(peak_date) - 1L),
        "-10-01"
      )) + 1L
    )
  )

# Convenience subsets
unreg <- peaks %>% filter(!regulated)

# "Interior" events: not pure rainfall end-members, not pure snowmelt
# end-members — the space where the snowmelt/ROS/mixed boundary operates
interior <- peaks %>%
  filter(
    !is.na(rf_short), !is.na(rf_long), !is.na(ros_flux_frac),
    rf_short < RAIN_EM_THRESH,
    rf_long  > SNOW_EM_THRESH
  )

message(sprintf(
  "Interior events (neither rainfall nor snowmelt EM): %d of %d (%.1f%%)",
  nrow(interior), nrow(peaks), 100 * nrow(interior) / nrow(peaks)
))

# ============================================================================
# SHARED THEME
# ============================================================================
theme_diag <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    plot.title       = element_text(size = 12, face = "bold"),
    plot.subtitle    = element_text(size = 10, colour = "grey40")
  )

type_colors <- c(
  rainfall    = "#e41a1c",
  snowmelt    = "#377eb8",
  ROS         = "#4daf4a",
  mixed       = "#ff7f00",
  unclassified = "grey60"
)

# ============================================================================
# 1. MARGINAL DISTRIBUTION OF ros_flux_frac (POPULATION)
# ============================================================================
message("Plot 1: marginal distribution of ros_flux_frac...")

dat1 <- peaks %>% filter(!is.na(ros_flux_frac))

p1a <- ggplot(dat1, aes(x = ros_flux_frac)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 80, fill = "steelblue", colour = "white", linewidth = 0.2) +
  geom_density(colour = "black", linewidth = 0.8) +
  geom_vline(xintercept = ROS_THRESH, linetype = "dashed", colour = "firebrick", linewidth = 0.8) +
  annotate("text", x = ROS_THRESH + 0.02, y = Inf, label = "ROS threshold",
           hjust = 0, vjust = 1.5, colour = "firebrick", size = 3.5) +
  labs(
    title    = "Marginal distribution of ros_flux_frac (all events)",
    subtitle = sprintf("n = %d peaks with non-missing ros_flux_frac", nrow(dat1)),
    x        = "ros_flux_frac (co-occurrence flux fraction, short window)",
    y        = "Density"
  ) +
  theme_diag

# Separate for unregulated only
p1b <- ggplot(dat1 %>% filter(!regulated), aes(x = ros_flux_frac)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 80, fill = "steelblue", colour = "white", linewidth = 0.2) +
  geom_density(colour = "black", linewidth = 0.8) +
  geom_vline(xintercept = ROS_THRESH, linetype = "dashed", colour = "firebrick", linewidth = 0.8) +
  labs(
    title    = "Marginal distribution of ros_flux_frac (unregulated only)",
    subtitle = sprintf("n = %d unregulated peaks", sum(!dat1$regulated)),
    x        = "ros_flux_frac",
    y        = "Density"
  ) +
  theme_diag

ggsave(file.path(out_dir, "01_ros_flux_frac_marginal.png"),
       p1a / p1b, width = 9, height = 8, dpi = 150)

# ============================================================================
# 2. ros_flux_frac BY CURRENT FLOOD TYPE
# ============================================================================
message("Plot 2: ros_flux_frac by flood_type...")

dat2 <- peaks %>%
  filter(!is.na(ros_flux_frac), flood_type != "unclassified") %>%
  mutate(flood_type = factor(flood_type,
                             levels = c("rainfall","snowmelt","ROS","mixed")))

p2 <- ggplot(dat2, aes(x = ros_flux_frac, fill = flood_type, colour = flood_type)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  geom_vline(xintercept = ROS_THRESH, linetype = "dashed",
             colour = "black", linewidth = 0.8) +
  scale_fill_manual(values = type_colors) +
  scale_colour_manual(values = type_colors) +
  facet_wrap(~flood_type, ncol = 2, scales = "free_y") +
  labs(
    title    = "ros_flux_frac distribution by flood type",
    subtitle = "Vertical dashed line = current ROS threshold (0.50)",
    x        = "ros_flux_frac",
    y        = "Density",
    fill     = NULL, colour = NULL
  ) +
  theme_diag +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "02_ros_flux_frac_by_type.png"),
       p2, width = 9, height = 7, dpi = 150)

# ============================================================================
# 3. JOINT DISTRIBUTION: ros_flux_frac x rf_long (INTERIOR EVENTS)
# ============================================================================
message("Plot 3: joint ros_flux_frac x rf_long (interior)...")

p3 <- ggplot(interior, aes(x = rf_long, y = ros_flux_frac)) +
  geom_hex(bins = 50) +
  geom_hline(yintercept = ROS_THRESH, linetype = "dashed",
             colour = "firebrick", linewidth = 0.8) +
  scale_fill_viridis_c(option = "magma", trans = "log10",
                       name = "Count\n(log10)") +
  labs(
    title    = "Joint distribution: rf_long × ros_flux_frac",
    subtitle = sprintf("Interior events only (rf_short < %.2f, rf_long > %.2f), n = %d",
                       RAIN_EM_THRESH, SNOW_EM_THRESH, nrow(interior)),
    x        = "rf_long (weighted rain fraction, long window)",
    y        = "ros_flux_frac (co-occurrence fraction)"
  ) +
  theme_diag

# Overlay current classification as colour
p3b <- ggplot(interior %>% filter(flood_type != "unclassified"),
              aes(x = rf_long, y = ros_flux_frac, colour = flood_type)) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_hline(yintercept = ROS_THRESH, linetype = "dashed",
             colour = "black", linewidth = 0.8) +
  scale_colour_manual(values = type_colors, name = "Flood type") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(
    title    = "rf_long × ros_flux_frac, coloured by current classification",
    subtitle = "Interior events only",
    x        = "rf_long",
    y        = "ros_flux_frac"
  ) +
  theme_diag

ggsave(file.path(out_dir, "03_joint_ros_vs_rf_long.png"),
       p3 / p3b, width = 9, height = 10, dpi = 150)

# ============================================================================
# 4. JOINT DISTRIBUTION: ros_flux_frac x rf_short (INTERIOR EVENTS)
# ============================================================================
message("Plot 4: joint ros_flux_frac x rf_short (interior)...")

p4 <- ggplot(interior, aes(x = rf_short, y = ros_flux_frac)) +
  geom_hex(bins = 50) +
  geom_hline(yintercept = ROS_THRESH, linetype = "dashed",
             colour = "firebrick", linewidth = 0.8) +
  scale_fill_viridis_c(option = "magma", trans = "log10",
                       name = "Count\n(log10)") +
  labs(
    title    = "Joint distribution: rf_short × ros_flux_frac",
    subtitle = sprintf("Interior events only, n = %d", nrow(interior)),
    x        = "rf_short (weighted rain fraction, short window)",
    y        = "ros_flux_frac"
  ) +
  theme_diag

p4b <- ggplot(interior %>% filter(flood_type != "unclassified"),
              aes(x = rf_short, y = ros_flux_frac, colour = flood_type)) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_hline(yintercept = ROS_THRESH, linetype = "dashed",
             colour = "black", linewidth = 0.8) +
  scale_colour_manual(values = type_colors, name = "Flood type") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(
    title    = "rf_short × ros_flux_frac, coloured by current classification",
    x        = "rf_short",
    y        = "ros_flux_frac"
  ) +
  theme_diag

ggsave(file.path(out_dir, "04_joint_ros_vs_rf_short.png"),
       p4 / p4b, width = 9, height = 10, dpi = 150)

# ============================================================================
# 5. THREE-WAY: rf_short x rf_long, colour = ros_flux_frac (INTERIOR)
# ============================================================================
message("Plot 5: feature triangle rf_short x rf_long, colour = ros_flux_frac...")

p5 <- ggplot(interior, aes(x = rf_long, y = rf_short, colour = ros_flux_frac)) +
  geom_point(alpha = 0.35, size = 0.9) +
  scale_colour_viridis_c(option = "plasma", name = "ros_flux_frac",
                         limits = c(0, 1)) +
  labs(
    title    = "rf_long × rf_short, coloured by ros_flux_frac",
    subtitle = "Interior events — the decision space for snowmelt vs ROS vs mixed",
    x        = "rf_long (long-window rain fraction)",
    y        = "rf_short (short-window rain fraction)"
  ) +
  theme_diag

ggsave(file.path(out_dir, "05_feature_triangle.png"),
       p5, width = 9, height = 7, dpi = 150)

# ============================================================================
# 6. SEASONAL STRUCTURE OF INTERIOR EVENTS
# ============================================================================
message("Plot 6: seasonal structure of interior events by current type...")

# Month labels for water year axis
wy_month_breaks <- c(1, 32, 62, 93, 124, 152, 182, 213, 244, 274, 305, 335)
wy_month_labels <- c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep")

interior_typed <- interior %>%
  filter(flood_type %in% c("snowmelt","ROS","mixed"))

p6a <- ggplot(interior_typed, aes(x = wy_doy, fill = flood_type, colour = flood_type)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  scale_x_continuous(breaks = wy_month_breaks, labels = wy_month_labels) +
  scale_fill_manual(values = type_colors) +
  scale_colour_manual(values = type_colors) +
  labs(
    title    = "Seasonal distribution of interior events by current type",
    subtitle = "Snowmelt, ROS, and Mixed only",
    x        = "Month (water year)",
    y        = "Density",
    fill     = "Flood type", colour = "Flood type"
  ) +
  theme_diag

# ros_flux_frac as a function of wy_doy — does it vary seasonally?
p6b <- ggplot(interior_typed, aes(x = wy_doy, y = ros_flux_frac, colour = flood_type)) +
  geom_point(alpha = 0.15, size = 0.7) +
  geom_smooth(method = "loess", span = 0.4, se = TRUE, linewidth = 1.0) +
  geom_hline(yintercept = ROS_THRESH, linetype = "dashed",
             colour = "black", linewidth = 0.8) +
  scale_x_continuous(breaks = wy_month_breaks, labels = wy_month_labels) +
  scale_colour_manual(values = type_colors) +
  labs(
    title    = "ros_flux_frac across the water year (interior events)",
    subtitle = "Loess smoother per type; dashed line = ROS threshold",
    x        = "Month (water year)",
    y        = "ros_flux_frac",
    colour   = "Flood type"
  ) +
  theme_diag

ggsave(file.path(out_dir, "06_seasonal_interior.png"),
       p6a / p6b, width = 10, height = 9, dpi = 150)

# ============================================================================
# 7. HARTIGAN'S DIP TEST: WITHIN-SITE BIMODALITY OF ros_flux_frac
# ============================================================================
message("Running Hartigan's dip test per site...")

# Test within each unregulated site with at least MIN_SITE_N interior events
dip_results <- interior %>%
  filter(!regulated) %>%
  group_by(site_no) %>%
  filter(n() >= MIN_SITE_N, !is.na(ros_flux_frac)) %>%
  summarise(
    n       = n(),
    dip     = dip.test(ros_flux_frac)$statistic,
    p_value = dip.test(ros_flux_frac)$p.value,
    .groups = "drop"
  ) %>%
  arrange(p_value)

message(sprintf(
  "\nDip test: %d of %d eligible sites show p < 0.05 (bimodal ros_flux_frac)",
  sum(dip_results$p_value < 0.05), nrow(dip_results)
))
message(sprintf(
  "  p < 0.10: %d sites | p < 0.05: %d sites | p < 0.01: %d sites",
  sum(dip_results$p_value < 0.10),
  sum(dip_results$p_value < 0.05),
  sum(dip_results$p_value < 0.01)
))

p7 <- ggplot(dip_results, aes(x = p_value)) +
  geom_histogram(bins = 30, fill = "steelblue", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.05, linetype = "dashed",
             colour = "firebrick", linewidth = 0.8) +
  annotate("text", x = 0.07, y = Inf, label = "p = 0.05",
           hjust = 0, vjust = 1.5, colour = "firebrick", size = 3.5) +
  labs(
    title    = "Hartigan's dip test: within-site bimodality of ros_flux_frac",
    subtitle = sprintf(
      "Unregulated sites with ≥ %d interior events (n = %d sites)\n%d of %d sites show p < 0.05",
      MIN_SITE_N, nrow(dip_results),
      sum(dip_results$p_value < 0.05), nrow(dip_results)
    ),
    x        = "Dip test p-value (low = bimodal)",
    y        = "Number of sites"
  ) +
  theme_diag

ggsave(file.path(out_dir, "07_dip_test_summary.png"),
       p7, width = 8, height = 5, dpi = 150)

# ============================================================================
# 8. PER-SITE DENSITY PLOTS (small multiples, unregulated, n >= MIN_SITE_N)
# ============================================================================
message("Plot 8: per-site ros_flux_frac densities...")

sites_for_density <- interior %>%
  filter(!regulated) %>%
  group_by(site_no) %>%
  filter(n() >= MIN_SITE_N) %>%
  ungroup() %>%
  left_join(dip_results %>% select(site_no, p_value), by = "site_no") %>%
  mutate(
    site_label = sprintf("%s\n(p=%.2f, n=%d)",
                         site_no,
                         p_value,
                         n()),
    bimodal    = p_value < 0.05
  )

n_sites <- n_distinct(sites_for_density$site_no)
message(sprintf("  Plotting %d sites", n_sites))

p8 <- ggplot(sites_for_density,
             aes(x = ros_flux_frac, fill = bimodal)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, colour = "white", linewidth = 0.15) +
  geom_density(colour = "black", linewidth = 0.5) +
  geom_vline(xintercept = ROS_THRESH, linetype = "dashed",
             colour = "firebrick", linewidth = 0.5) +
  scale_fill_manual(values = c("TRUE" = "#4daf4a", "FALSE" = "steelblue"),
                    labels = c("TRUE" = "p < 0.05 (bimodal)",
                               "FALSE" = "p ≥ 0.05 (unimodal)"),
                    name = "Dip test") +
  facet_wrap(~site_no, ncol = 6) +
  labs(
    title    = "Per-site ros_flux_frac distributions (interior events, unregulated)",
    subtitle = sprintf("Sites with ≥ %d interior events | red dashed = ROS threshold (0.50)", MIN_SITE_N),
    x        = "ros_flux_frac",
    y        = "Density"
  ) +
  theme_diag +
  theme(
    axis.text  = element_text(size = 6),
    strip.text = element_text(size = 6)
  )

# Height scales with number of sites
plot_height <- max(8, ceiling(n_sites / 6) * 2.5)
ggsave(file.path(out_dir, "08_ros_frac_site_densities.png"),
       p8, width = 15, height = plot_height, dpi = 150)

# ============================================================================
# CONSOLE SUMMARY
# ============================================================================
summary_lines <- c(
  "=== ROS FEATURE SPACE DIAGNOSTIC SUMMARY ===",
  sprintf("Run date       : %s", Sys.Date()),
  sprintf("Classified file: %s", classified_file),
  "",
  "--- Population counts ---",
  sprintf("Total peaks                   : %d", nrow(peaks)),
  sprintf("With non-missing ros_flux_frac: %d", nrow(dat1)),
  sprintf("Interior events               : %d (%.1f%%)",
          nrow(interior), 100 * nrow(interior) / nrow(peaks)),
  "",
  "--- Current type breakdown (interior) ---"
)

type_summary <- interior %>%
  count(flood_type) %>%
  mutate(pct = 100 * n / sum(n))

summary_lines <- c(summary_lines,
  capture.output(print(type_summary, n = Inf)),
  "",
  "--- Dip test (within-site bimodality of ros_flux_frac) ---",
  sprintf("Eligible sites (unregulated, ≥ %d interior peaks): %d",
          MIN_SITE_N, nrow(dip_results)),
  sprintf("  p < 0.10 : %d sites (%.1f%%)",
          sum(dip_results$p_value < 0.10),
          100 * mean(dip_results$p_value < 0.10)),
  sprintf("  p < 0.05 : %d sites (%.1f%%)",
          sum(dip_results$p_value < 0.05),
          100 * mean(dip_results$p_value < 0.05)),
  sprintf("  p < 0.01 : %d sites (%.1f%%)",
          sum(dip_results$p_value < 0.01),
          100 * mean(dip_results$p_value < 0.01)),
  "",
  "Most bimodal sites:",
  capture.output(print(head(dip_results, 10), n = 10)),
  "",
  sprintf("Outputs written to: %s", out_dir)
)

writeLines(summary_lines,
           file.path(out_dir, "diag_ros_feature_space.txt"))
message(paste(summary_lines, collapse = "\n"))
message("\nDone. All outputs written to ", out_dir)