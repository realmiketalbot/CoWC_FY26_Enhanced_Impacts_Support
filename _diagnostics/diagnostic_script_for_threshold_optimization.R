#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Joint Optimization of Window Length and Threshold
# ============================================================================
# For each candidate (window_length, threshold) combination, computes
# sensitivity and specificity for rainfall and snowmelt classification
# using the data-driven end-members from the 60-day accumulation diagnostic.
#
# Sensitivity: fraction of true end-members exceeding the threshold
#              within the window
# Specificity: fraction of opposing end-members staying below the threshold
#              within the window
# Youden's J:  sensitivity + specificity - 1  (combined score, 0-1)
#
# Evaluated separately for rainfall and snowmelt, over a grid of:
#   - Window lengths: 1-14 days (rainfall), 5-60 days (snowmelt)
#   - Thresholds:     2-30 mm
#
# End-member populations are weighted equally by type (not by count) so
# the optimization is not dominated by the larger rainfall end-member pool.
#
# INPUTS:
#   - _data/accum_curves_60d.rds    : long-format accumulation curves
#   - _data/endmembers_60d.rds      : end-member peak IDs with 60d totals
#
# OUTPUTS:
#   - _data/window_threshold_optim.rds   : full grid results
#   - _diagnostics/accum_curves/optim_*  : figures
# ============================================================================

library(tidyverse)
library(patchwork)

# ============================================================================
# CONFIGURATION
# ============================================================================
accum_file   <- "_data/accum_curves_60d.rds"
endmem_file  <- "_data/endmembers_60d.rds"
optim_output <- "_data/window_threshold_optim.rds"
out_dir      <- "_diagnostics/accum_curves"

# Candidate grid
RAIN_WINDOWS  <- 1:14
MELT_WINDOWS  <- 5:60
THRESHOLDS    <- seq(2, 30, by = 1)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading accumulation curves and end-members...")

accum   <- readRDS(accum_file)
endmems <- readRDS(endmem_file) %>%
  filter(endmember_type != "neither") %>%
  select(site_no, peak_id, endmember_type)

message(sprintf("End-members: %d rain, %d melt",
                sum(endmems$endmember_type == "rain"),
                sum(endmems$endmember_type == "melt")))

# Join accumulation curves to end-member classification
accum_em <- accum %>%
  inner_join(endmems, by = c("site_no", "peak_id"))

# ============================================================================
# BUILD EVALUATION TABLE
# ============================================================================
message("Extracting cumulative values at each candidate window length...")

# Rain end-members: primary signal = cumul_rain_mm, opposing = cumul_melt_mm
rain_em_vals <- accum_em %>%
  filter(
    endmember_type %in% c("rain", "melt"),
    lookback_day %in% RAIN_WINDOWS
  ) %>%
  select(site_no, peak_id, endmember_type, lookback_day,
         primary   = cumul_rain_mm,
         opposing  = cumul_melt_mm) %>%
  mutate(variable = "rainfall")

# Melt end-members: primary signal = cumul_melt_mm, opposing = cumul_rain_mm
melt_em_vals <- accum_em %>%
  filter(
    endmember_type %in% c("rain", "melt"),
    lookback_day %in% MELT_WINDOWS
  ) %>%
  select(site_no, peak_id, endmember_type, lookback_day,
         primary   = cumul_melt_mm,
         opposing  = cumul_rain_mm) %>%
  mutate(variable = "snowmelt")

# ============================================================================
# COMPUTE SENSITIVITY, SPECIFICITY, YOUDEN'S J OVER GRID
# ============================================================================
message("Computing sensitivity/specificity grid...")

compute_grid <- function(em_vals, variable_name, windows, thresholds) {

  true_pos  <- em_vals %>%
    filter(endmember_type == if_else(variable_name == "rainfall",
                                     "rain", "melt"))
  true_neg  <- em_vals %>%
    filter(endmember_type == if_else(variable_name == "rainfall",
                                     "melt", "rain"))

  n_pos <- n_distinct(true_pos$peak_id)
  n_neg <- n_distinct(true_neg$peak_id)

  message(sprintf("  %s: %d positive end-members, %d negative end-members",
                  variable_name, n_pos, n_neg))

  grid <- expand_grid(
    window_days = windows,
    threshold   = thresholds
  )

  grid %>%
    rowwise() %>%
    mutate(
      sensitivity = {
        tp_at_window <- true_pos %>% filter(lookback_day == window_days)
        if (nrow(tp_at_window) == 0) NA_real_
        else mean(tp_at_window$primary >= threshold, na.rm = TRUE)
      },
      specificity = {
        tn_at_window <- true_neg %>% filter(lookback_day == window_days)
        if (nrow(tn_at_window) == 0) NA_real_
        else mean(tn_at_window$primary < threshold, na.rm = TRUE)
      }
    ) %>%
    ungroup() %>%
    mutate(
      youdens_j = sensitivity + specificity - 1,
      variable  = variable_name
    )
}

message("Computing rainfall grid...")
rain_grid <- compute_grid(rain_em_vals, "rainfall", RAIN_WINDOWS, THRESHOLDS)

message("Computing snowmelt grid...")
melt_grid <- compute_grid(melt_em_vals, "snowmelt", MELT_WINDOWS, THRESHOLDS)

full_grid <- bind_rows(rain_grid, melt_grid)

saveRDS(full_grid, optim_output)
message(sprintf("Saved optimization grid to %s", optim_output))

# ============================================================================
# FIND OPTIMAL COMBINATIONS
# ============================================================================
optimal <- full_grid %>%
  filter(!is.na(youdens_j)) %>%
  group_by(variable) %>%
  slice_max(youdens_j, n = 1, with_ties = FALSE) %>%
  ungroup()

message("\n=== OPTIMAL WINDOW/THRESHOLD COMBINATIONS ===")
print(optimal %>% select(variable, window_days, threshold,
                          sensitivity, specificity, youdens_j))

rain_opt  <- optimal %>% filter(variable == "rainfall")
melt_opt  <- optimal %>% filter(variable == "snowmelt")

# ============================================================================
# FIGURES
# ============================================================================

# --- Fig 1: Youden's J heatmap for rainfall ---
p1 <- rain_grid %>%
  filter(!is.na(youdens_j)) %>%
  ggplot(aes(x = window_days, y = threshold, fill = youdens_j)) +
  geom_tile() +
  geom_point(data = rain_opt,
             aes(x = window_days, y = threshold),
             color = "white", size = 3, shape = 8) +
  scale_fill_distiller(palette = "RdYlGn", direction = 1,
                       limits = c(0, 1), name = "Youden's J") +
  scale_x_continuous(breaks = c(1, 2, 4, 6, 8, 10, 12, 14)) +
  scale_y_continuous(breaks = seq(2, 30, by = 4)) +
  annotate("segment",
           x = rain_opt$window_days, xend = rain_opt$window_days,
           y = -Inf, yend = rain_opt$threshold,
           linetype = "dashed", color = "white", linewidth = 0.5) +
  annotate("segment",
           x = -Inf, xend = rain_opt$window_days,
           y = rain_opt$threshold, yend = rain_opt$threshold,
           linetype = "dashed", color = "white", linewidth = 0.5) +
  labs(
    title    = "Rainfall: Youden's J across window × threshold grid",
    subtitle = sprintf("★ = optimal: %d-day window, %.0f mm threshold (J = %.3f)",
                       rain_opt$window_days, rain_opt$threshold,
                       rain_opt$youdens_j),
    x = "Window length (days)",
    y = "Threshold (mm)"
  ) +
  theme_bw(base_size = 11)

# --- Fig 2: Youden's J heatmap for snowmelt ---
p2 <- melt_grid %>%
  filter(!is.na(youdens_j)) %>%
  ggplot(aes(x = window_days, y = threshold, fill = youdens_j)) +
  geom_tile() +
  geom_point(data = melt_opt,
             aes(x = window_days, y = threshold),
             color = "white", size = 3, shape = 8) +
  scale_fill_distiller(palette = "RdYlGn", direction = 1,
                       limits = c(0, 1), name = "Youden's J") +
  scale_x_continuous(breaks = c(5, 10, 15, 20, 30, 45, 60)) +
  scale_y_continuous(breaks = seq(2, 30, by = 4)) +
  annotate("segment",
           x = melt_opt$window_days, xend = melt_opt$window_days,
           y = -Inf, yend = melt_opt$threshold,
           linetype = "dashed", color = "white", linewidth = 0.5) +
  annotate("segment",
           x = -Inf, xend = melt_opt$window_days,
           y = melt_opt$threshold, yend = melt_opt$threshold,
           linetype = "dashed", color = "white", linewidth = 0.5) +
  labs(
    title    = "Snowmelt: Youden's J across window × threshold grid",
    subtitle = sprintf("★ = optimal: %d-day window, %.0f mm threshold (J = %.3f)",
                       melt_opt$window_days, melt_opt$threshold,
                       melt_opt$youdens_j),
    x = "Window length (days)",
    y = "Threshold (mm)"
  ) +
  theme_bw(base_size = 11)

# --- Fig 3: Sensitivity and specificity curves at optimal window ---
sens_spec_curves <- full_grid %>%
  filter(
    (variable == "rainfall" & window_days == rain_opt$window_days) |
    (variable == "snowmelt" & window_days == melt_opt$window_days)
  ) %>%
  pivot_longer(cols = c(sensitivity, specificity),
               names_to = "metric", values_to = "value")

opt_lines <- bind_rows(
  rain_opt %>% mutate(variable = "rainfall"),
  melt_opt %>% mutate(variable = "snowmelt")
)

p3 <- ggplot(sens_spec_curves,
             aes(x = threshold, y = value,
                 color = metric, linetype = metric)) +
  geom_line(linewidth = 0.9) +
  geom_vline(data = opt_lines,
             aes(xintercept = threshold),
             linetype = "dashed", color = "grey30", linewidth = 0.7) +
  geom_text(data = opt_lines,
            aes(x = threshold + 0.5, y = 0.05,
                label = sprintf("%.0f mm", threshold)),
            hjust = 0, color = "grey30", size = 3, inherit.aes = FALSE) +
  scale_color_manual(values = c(sensitivity = "#d73027",
                                specificity = "#4575b4")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  facet_wrap(~ variable, ncol = 2,
             labeller = labeller(variable = c(
               rainfall = "Rainfall (optimal window)",
               snowmelt = "Snowmelt (optimal window)"
             ))) +
  labs(
    title    = "Sensitivity and specificity vs threshold at optimal window length",
    subtitle = "Dashed vertical line = optimal threshold (max Youden's J)",
    x        = "Threshold (mm)",
    y        = "Value",
    color    = NULL, linetype = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# --- Fig 4: Youden's J vs window length at optimal threshold ---
j_vs_window <- full_grid %>%
  filter(
    (variable == "rainfall" & threshold == rain_opt$threshold) |
    (variable == "snowmelt" & threshold == melt_opt$threshold)
  )

p4 <- ggplot(j_vs_window,
             aes(x = window_days, y = youdens_j, color = variable)) +
  geom_line(linewidth = 0.9) +
  geom_point(data = opt_lines,
             aes(x = window_days, y = youdens_j),
             size = 3, shape = 8) +
  geom_vline(xintercept = 8, linetype = "dashed",
             color = "grey50", linewidth = 0.7) +
  annotate("text", x = 8.3, y = 0.05,
           label = "Current\n8-day window",
           hjust = 0, color = "grey50", size = 3) +
  scale_color_manual(values = c(rainfall = "#2166ac",
                                snowmelt = "#4dac26")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title    = "Youden's J vs window length at optimal threshold",
    subtitle = "★ = optimal window. Shows marginal gain from extending window.",
    x        = "Window length (days)",
    y        = "Youden's J",
    color    = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# ============================================================================
# COMBINE AND SAVE
# ============================================================================
combined <- (p1 + p2) / p3 / p4 +
  plot_annotation(
    title   = "Joint Window Length × Threshold Optimization",
    caption = "Youden's J = sensitivity + specificity - 1. End-member types weighted equally."
  )

out_file <- file.path(out_dir, "optim_window_threshold.png")
ggsave(out_file, combined, width = 12, height = 14, dpi = 150)
message(sprintf("\nSaved figure to %s", out_file))

# ============================================================================
# SUMMARY TABLE
# ============================================================================
message("\n=== COMPARISON: DAVENPORT VS OPTIMAL ===")
message(sprintf("Rainfall:  Davenport = 8-day / 10mm"))
message(sprintf("           Optimal   = %d-day / %.0f mm  (J = %.3f)",
                rain_opt$window_days, rain_opt$threshold,
                rain_opt$youdens_j))

davenport_rain_j <- rain_grid %>%
  filter(window_days == 8, threshold == 10) %>%
  pull(youdens_j)
message(sprintf("           Davenport J = %.3f", davenport_rain_j))

message(sprintf("\nSnowmelt:  Davenport = 8-day / 10mm"))
message(sprintf("           Optimal   = %d-day / %.0f mm  (J = %.3f)",
                melt_opt$window_days, melt_opt$threshold,
                melt_opt$youdens_j))

davenport_melt_j <- melt_grid %>%
  filter(window_days == 8, threshold == 10) %>%
  pull(youdens_j)
message(sprintf("           Davenport J = %.3f", davenport_melt_j))