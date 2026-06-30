#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Snowmelt Lookback Window from Marginal Melt Flux
# ============================================================================
# For snowmelt end-member peaks, characterises how cumulative melt flux
# accumulates over the 60-day lookback window. The optimal window is expressed
# as a multiplier of each peak's site-specific search_radius (Livneh et al.
# independence criterion), making the result directly applicable to the
# per-peak variable window approach.
#
# Per-peak optimal window (days) is found under two criteria:
#   Marginal : smallest d where melt_mm_day < MARGINAL_THRESH_PCT of 60d total
#   Cumulative: smallest d where cumul_melt >= CUMUL_TARGET_PCT of 60d total
#
# Per-peak multiplier = optimal_window_days / search_radius
#
# The distribution of per-peak multipliers across all snowmelt end-members
# gives a data-driven, site-adaptive scaling factor.
#
# INPUTS:
#   - _data/accum_curves_60d.rds
#   - _data/pot_peaks_all.rds      (for search_radius per peak)
#   - _data/site_metadata.rds
#
# OUTPUTS:
#   - _data/snowmelt_window.rds
#   - _diagnostics/accum_curves/snowmelt_window.png
# ============================================================================

library(tidyverse)
library(patchwork)

# ============================================================================
# CONFIGURATION
# ============================================================================
accum_file   <- "_data/accum_curves_60d.rds"
peaks_file   <- "_data/pot_peaks_all.rds"
meta_file    <- "_data/site_metadata.rds"
output_file  <- "_data/snowmelt_window.rds"
out_dir      <- "_diagnostics/accum_curves"

# End-member threshold
SNOW_EM_THRESH    <- 0.05   # weighted rain frac < this => snowmelt end-member
NOISE_FLOOR_MM    <- 0.1    # mm/day below which flux treated as absent

# Optimal window criteria
MARGINAL_THRESH_PCT <- 1.0  # % of 60-day total below which marginal is negligible
CUMUL_TARGET_PCT    <- 0.90 # window where cumulative reaches this % of 60-day total

# Regulation threshold (Davenport et al. 2020: 1,000 acre-feet = 1.233 mcm)
REGULATION_THRESH_MCM <- 1.2335

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD AND FILTER
# ============================================================================
message("Loading data...")
accum <- readRDS(accum_file)
meta  <- readRDS(meta_file)

# Pull search_radius per peak from pot_peaks_all
peak_search_radius <- readRDS(peaks_file) %>%
  distinct(site_no, peak_id, search_radius)

unregulated_sites <- meta %>%
  filter(!regulated) %>%
  pull(site_no)

message(sprintf("Unregulated sites: %d of %d total",
                length(unregulated_sites), nrow(meta)))

accum_unreg <- accum %>%
  filter(site_no %in% unregulated_sites) %>%
  left_join(peak_search_radius, by = c("site_no", "peak_id"))

# ============================================================================
# IDENTIFY SNOWMELT END-MEMBERS
# ============================================================================
message("Identifying snowmelt end-members...")

# Compute weighted rain fraction over full 60-day window per peak
rain_frac_60d <- accum_unreg %>%
  mutate(
    rain_active = if_else(rain_mm_day >= NOISE_FLOOR_MM, rain_mm_day, 0),
    melt_active = if_else(melt_mm_day >= NOISE_FLOOR_MM, melt_mm_day, 0),
    total_flux  = rain_active + melt_active
  ) %>%
  group_by(site_no, peak_id) %>%
  summarise(
    total_rain_flux    = sum(rain_active, na.rm = TRUE),
    total_melt_flux    = sum(melt_active, na.rm = TRUE),
    total_flux         = sum(total_flux,  na.rm = TRUE),
    n_days_data        = first(n_days_data),
    .groups            = "drop"
  ) %>%
  mutate(
    weighted_rain_frac = if_else(total_flux > 0,
                                 total_rain_flux / total_flux,
                                 NA_real_)
  )

snow_ems <- rain_frac_60d %>%
  filter(!is.na(weighted_rain_frac),
         weighted_rain_frac < SNOW_EM_THRESH,
         n_days_data >= 57)     # require near-complete 60-day window

message(sprintf("Snowmelt end-members: %d peaks from %d sites",
                nrow(snow_ems), n_distinct(snow_ems$site_no)))

# ============================================================================
# COMPUTE MARGINAL AND CUMULATIVE MELT FLUX PER LOOKBACK DAY
# ============================================================================
message("Computing accumulation curves for snowmelt end-members...")

# Join end-members back to daily accumulation curves
# lookback_day = 1 is the peak day, increasing values go further back
accum_snow <- accum_unreg %>%
  semi_join(snow_ems, by = c("site_no", "peak_id")) %>%
  left_join(snow_ems %>% select(site_no, peak_id, total_melt_flux),
            by = c("site_no", "peak_id")) %>%
  arrange(site_no, peak_id, lookback_day) %>%
  # Marginal contribution = melt on that day (accum curves already store
  # daily values; cumul_melt_mm is the running sum from day 1 outward)
  group_by(site_no, peak_id) %>%
  mutate(
    # Normalise marginal by each peak's 60-day total melt
    marginal_pct = if_else(
      total_melt_flux > 0,
      melt_mm_day / total_melt_flux * 100,
      NA_real_
    ),
    # Fraction of 60-day total accumulated by this lookback day
    cumul_pct = if_else(
      total_melt_flux > 0,
      cumul_melt_mm / total_melt_flux * 100,
      NA_real_
    )
  ) %>%
  ungroup()

# ============================================================================
# SUMMARISE ACROSS END-MEMBERS
# ============================================================================
curve_summary <- accum_snow %>%
  group_by(lookback_day) %>%
  summarise(
    # Cumulative melt (absolute mm)
    median_cumul_mm  = median(cumul_melt_mm,  na.rm = TRUE),
    q10_cumul_mm     = quantile(cumul_melt_mm, 0.10, na.rm = TRUE),
    q90_cumul_mm     = quantile(cumul_melt_mm, 0.90, na.rm = TRUE),
    # Cumulative melt (% of 60-day total)
    median_cumul_pct = median(cumul_pct,  na.rm = TRUE),
    q10_cumul_pct    = quantile(cumul_pct, 0.10, na.rm = TRUE),
    q90_cumul_pct    = quantile(cumul_pct, 0.90, na.rm = TRUE),
    # Marginal contribution (% of 60-day total)
    median_marginal_pct = median(marginal_pct, na.rm = TRUE),
    q10_marginal_pct    = quantile(marginal_pct, 0.10, na.rm = TRUE),
    q90_marginal_pct    = quantile(marginal_pct, 0.90, na.rm = TRUE),
    n_peaks          = n(),
    .groups          = "drop"
  )

# ============================================================================
# PER-PEAK OPTIMAL WINDOW AND MULTIPLIER
# ============================================================================
# For each end-member peak, find the smallest lookback day satisfying each
# criterion, then express as a ratio to that peak's search_radius.

per_peak_windows <- accum_snow %>%
  group_by(site_no, peak_id, search_radius) %>%
  summarise(
    # Marginal criterion: first day where marginal_pct < threshold
    opt_marginal = {
      d <- lookback_day[marginal_pct < MARGINAL_THRESH_PCT]
      if (length(d) == 0) NA_integer_ else min(d)
    },
    # Cumulative criterion: first day where cumul_pct >= target
    opt_cumul = {
      d <- lookback_day[!is.na(cumul_pct) & cumul_pct >= CUMUL_TARGET_PCT * 100]
      if (length(d) == 0) NA_integer_ else min(d)
    },
    .groups = "drop"
  ) %>%
  mutate(
    multiplier_marginal = opt_marginal / search_radius,
    multiplier_cumul    = opt_cumul    / search_radius
  )

# Summary statistics on multipliers
mult_summary <- per_peak_windows %>%
  summarise(
    across(
      c(multiplier_marginal, multiplier_cumul),
      list(
        median = ~ median(.x, na.rm = TRUE),
        q25    = ~ quantile(.x, 0.25, na.rm = TRUE),
        q75    = ~ quantile(.x, 0.75, na.rm = TRUE),
        n_valid = ~ sum(!is.na(.x))
      ),
      .names = "{.col}__{.fn}"
    )
  )

message(sprintf("\n=== PER-PEAK MULTIPLIER DISTRIBUTION ==="))
message(sprintf("  Marginal criterion  : median = %.2f (IQR: %.2f – %.2f), n = %d",
                mult_summary$multiplier_marginal__median,
                mult_summary$multiplier_marginal__q25,
                mult_summary$multiplier_marginal__q75,
                mult_summary$multiplier_marginal__n_valid))
message(sprintf("  Cumulative criterion: median = %.2f (IQR: %.2f – %.2f), n = %d",
                mult_summary$multiplier_cumul__median,
                mult_summary$multiplier_cumul__q25,
                mult_summary$multiplier_cumul__q75,
                mult_summary$multiplier_cumul__n_valid))

# ============================================================================
# SAVE
# ============================================================================
results <- list(
  per_peak_windows     = per_peak_windows,
  mult_summary         = mult_summary,
  curve_summary        = curve_summary,
  marginal_thresh_pct  = MARGINAL_THRESH_PCT,
  cumul_target_pct     = CUMUL_TARGET_PCT,
  snow_em_thresh       = SNOW_EM_THRESH,
  n_snow_em            = nrow(snow_ems)
)
saveRDS(results, output_file)
message(sprintf("Saved results to %s", output_file))

# ============================================================================
# FIGURES
# ============================================================================

med_mult_marg  <- mult_summary$multiplier_marginal__median
med_mult_cumul <- mult_summary$multiplier_cumul__median
q25_marg       <- mult_summary$multiplier_marginal__q25
q75_marg       <- mult_summary$multiplier_marginal__q75
q25_cumul      <- mult_summary$multiplier_cumul__q25
q75_cumul      <- mult_summary$multiplier_cumul__q75

# --- Fig 1: Marginal melt contribution vs lookback day ---
p1 <- ggplot(curve_summary, aes(x = lookback_day)) +
  geom_ribbon(aes(ymin = q10_marginal_pct, ymax = q90_marginal_pct),
              fill = "#2166ac", alpha = 0.15) +
  geom_line(aes(y = median_marginal_pct), color = "#2166ac", linewidth = 1.0) +
  geom_hline(yintercept = MARGINAL_THRESH_PCT,
             linetype = "dotted", color = "#d73027", linewidth = 0.8) +
  annotate("text", x = 2, y = MARGINAL_THRESH_PCT + 0.15,
           label = sprintf("%.1f%% threshold", MARGINAL_THRESH_PCT),
           hjust = 0, color = "#d73027", size = 3.5) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 30, 45, 60)) +
  labs(
    title    = "Marginal melt contribution per lookback day — snowmelt end-members",
    subtitle = sprintf(
      "Median ± 10th/90th percentile. Threshold at %.1f%% of 60-day total.",
      MARGINAL_THRESH_PCT),
    x = "Lookback days (1 = peak day)",
    y = "Marginal melt (% of 60-day total)"
  ) +
  theme_bw(base_size = 11)

# --- Fig 2: Cumulative melt (% of 60-day total) ---
p2 <- ggplot(curve_summary, aes(x = lookback_day)) +
  geom_ribbon(aes(ymin = q10_cumul_pct, ymax = q90_cumul_pct),
              fill = "#2166ac", alpha = 0.15) +
  geom_line(aes(y = median_cumul_pct), color = "#2166ac", linewidth = 1.0) +
  geom_hline(yintercept = CUMUL_TARGET_PCT * 100,
             linetype = "dashed", color = "#4dac26", linewidth = 0.8) +
  annotate("text", x = 2, y = CUMUL_TARGET_PCT * 100 + 2,
           label = sprintf("%.0f%% target", CUMUL_TARGET_PCT * 100),
           hjust = 0, color = "#4dac26", size = 3.5) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 30, 45, 60)) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Cumulative melt as % of 60-day total — snowmelt end-members",
    subtitle = sprintf("Median ± 10th/90th percentile. Green dashed = %.0f%% target.",
                       CUMUL_TARGET_PCT * 100),
    x = "Lookback days (1 = peak day)",
    y = "% of 60-day cumulative melt"
  ) +
  theme_bw(base_size = 11)

# --- Fig 3: Per-peak multiplier distributions ---
mult_long <- per_peak_windows %>%
  select(site_no, peak_id, search_radius,
         Marginal   = multiplier_marginal,
         Cumulative = multiplier_cumul) %>%
  pivot_longer(c(Marginal, Cumulative),
               names_to = "criterion", values_to = "multiplier") %>%
  filter(!is.na(multiplier), multiplier <= 10)   # clip extreme outliers for display

vlines <- tibble(
  criterion  = c("Marginal", "Cumulative"),
  median_val = c(med_mult_marg, med_mult_cumul),
  color      = c("#d73027", "#4dac26")
)

p3 <- ggplot(mult_long, aes(x = multiplier, fill = criterion)) +
  geom_histogram(binwidth = 0.25, color = NA, alpha = 0.55,
                 position = "identity") +
  geom_vline(data = vlines,
             aes(xintercept = median_val, color = color),
             linetype = "dashed", linewidth = 0.9, show.legend = FALSE) +
  scale_fill_manual(
    values = c("Marginal" = "#d73027", "Cumulative" = "#4dac26"),
    name   = "Criterion"
  ) +
  scale_color_identity() +
  annotate("text", x = med_mult_marg  + 0.1, y = Inf, vjust = 1.5, hjust = 0,
           label = sprintf("Marginal median: %.2f", med_mult_marg),
           color = "#d73027", size = 3.5) +
  annotate("text", x = med_mult_cumul + 0.1, y = Inf, vjust = 3.5, hjust = 0,
           label = sprintf("Cumulative median: %.2f", med_mult_cumul),
           color = "#4dac26", size = 3.5) +
  labs(
    title    = "Per-peak optimal window as a multiple of search_radius",
    subtitle = sprintf(
      paste0("Marginal criterion: median = %.2f (IQR: %.2f–%.2f). ",
             "Cumulative criterion: median = %.2f (IQR: %.2f–%.2f)."),
      med_mult_marg,  q25_marg,  q75_marg,
      med_mult_cumul, q25_cumul, q75_cumul),
    x = "Optimal window / search_radius",
    y = "Count of end-member peaks"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.75))

# ============================================================================
# COMBINE AND SAVE
# ============================================================================
combined <- p1 / p2 / p3 +
  plot_annotation(
    title   = "Snowmelt Lookback Window — Ratio to search_radius",
    caption = sprintf(
      paste0("Snowmelt end-members: rain frac < %.2f over 60-day window, ",
             "unregulated sites only (n = %d peaks).\n",
             "Marginal criterion  (<%.1f%% of 60d total): ",
             "median multiplier = %.2f (IQR: %.2f–%.2f).\n",
             "Cumulative criterion (>=%.0f%% of 60d total): ",
             "median multiplier = %.2f (IQR: %.2f–%.2f)."),
      SNOW_EM_THRESH, nrow(snow_ems),
      MARGINAL_THRESH_PCT,  med_mult_marg,  q25_marg,  q75_marg,
      CUMUL_TARGET_PCT * 100, med_mult_cumul, q25_cumul, q75_cumul
    )
  )

out_file <- file.path(out_dir, "snowmelt_window.png")
ggsave(out_file, combined, width = 10, height = 14, dpi = 150)
message(sprintf("Saved figure to %s", out_file))