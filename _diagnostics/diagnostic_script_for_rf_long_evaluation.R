#!/usr/bin/env Rscript

# ============================================================================
# Diagnostic: Is rf_long necessary, or does rf_short suffice for snowmelt?
# ============================================================================
# The adaptive classifier uses rf_long (3x search_radius window) exclusively
# to identify snowmelt events (rf_long <= 0.05). This diagnostic asks whether
# rf_short (search_radius window) is already near-zero for those same events,
# or whether rf_long is doing genuine reclassification work.
#
# Key output:
#   1. Joint scatter of rf_short vs rf_long for snowmelt-classified peaks
#   2. Marginal distributions of rf_short for snowmelt peaks — how often is
#      rf_short also <= 0.05, vs. elevated (rf_long is doing work)?
#   3. Counterfactual: how many currently-snowmelt peaks would be lost if
#      the snowmelt criterion switched to rf_short <= 0.05 instead?
#   4. For peaks where rf_long is doing work (rf_short > 0.05 but rf_long
#      <= 0.05), what does the accumulation profile look like?  We summarise
#      with rf_short distribution and flag any seasonal pattern.
# ============================================================================

library(tidyverse)

classified_file <- "_data/pot_peaks_classified_adaptive.rds"

SNOW_EM_THRESH  <- 0.05   # must match 5b config
RAIN_EM_THRESH  <- 0.95

# ============================================================================
# LOAD
# ============================================================================
dat <- readRDS(classified_file)

message("Columns present: ", paste(names(dat), collapse = ", "))
message(sprintf("Total peaks: %d | Sites: %d\n",
                nrow(dat), n_distinct(dat$site_no)))

# Restrict to unregulated, data-sufficient peaks for all analyses
dat_unreg <- dat %>%
  filter(!regulated, flood_type != "unclassified")

message(sprintf("Unregulated data-sufficient peaks: %d\n", nrow(dat_unreg)))

# ============================================================================
# 1. SNOWMELT PEAKS: JOINT rf_short vs rf_long
# ============================================================================
snowmelt <- dat_unreg %>% filter(flood_type == "snowmelt")

message(sprintf("Snowmelt peaks (unregulated): %d\n", nrow(snowmelt)))

# Categorise by whether rf_long is doing unique work
snowmelt <- snowmelt %>%
  mutate(
    rf_long_working = case_when(
      rf_short <= SNOW_EM_THRESH & rf_long <= SNOW_EM_THRESH ~ "Both low (redundant)",
      rf_short >  SNOW_EM_THRESH & rf_long <= SNOW_EM_THRESH ~ "rf_long only (adds value)",
      TRUE ~ "Other"   # shouldn't occur if classification is correct
    )
  )

cat("=== rf_long necessity breakdown for snowmelt peaks ===\n")
snowmelt %>%
  count(rf_long_working) %>%
  mutate(pct = n / sum(n) * 100) %>%
  print()

cat("\n")

# Summary stats for the rf_long-working subset
rf_long_unique <- snowmelt %>% filter(rf_long_working == "rf_long only (adds value)")
message(sprintf(
  "rf_long-unique events: %d (%.1f%% of snowmelt)\n",
  nrow(rf_long_unique),
  100 * nrow(rf_long_unique) / nrow(snowmelt)
))

if (nrow(rf_long_unique) > 0) {
  cat("rf_short distribution for rf_long-unique events:\n")
  print(summary(rf_long_unique$rf_short))
  cat("\n")

  # Seasonal pattern (DOY)
  rf_long_unique <- rf_long_unique %>%
    mutate(doy = as.integer(format(peak_date, "%j")),
           month = as.integer(format(peak_date, "%m")))
  cat("Month distribution for rf_long-unique events:\n")
  print(table(rf_long_unique$month))
  cat("\n")
}

# ============================================================================
# 2. COUNTERFACTUAL — what if snowmelt used rf_short instead of rf_long?
# ============================================================================
# Re-apply the decision tree using rf_short for the snowmelt step
dat_counter <- dat_unreg %>%
  mutate(
    flood_type_counter = case_when(
      rf_short >= RAIN_EM_THRESH  ~ "rainfall",
      rf_short <= SNOW_EM_THRESH  ~ "snowmelt",   # <-- changed criterion
      ros_flux_frac > 0.50        ~ "ROS",
      TRUE                        ~ "mixed"
    )
  )

cat("=== Original classification (rf_long for snowmelt) ===\n")
dat_unreg %>% count(flood_type) %>% mutate(pct = n / sum(n) * 100) %>% print()

cat("\n=== Counterfactual (rf_short for snowmelt) ===\n")
dat_counter %>% count(flood_type_counter) %>% mutate(pct = n / sum(n) * 100) %>% print()

# How many peaks change class?
changes <- dat_counter %>%
  filter(as.character(flood_type) != flood_type_counter) %>%
  count(flood_type, flood_type_counter, name = "n_peaks") %>%
  arrange(desc(n_peaks))

cat("\n=== Classification changes under counterfactual ===\n")
print(changes)
cat("\n")

# ============================================================================
# 3. PLOTS
# ============================================================================
outdir <- "_diagnostics/rf_long_necessity"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# --- Plot 1: Joint scatter rf_short vs rf_long for snowmelt peaks -----------
p1 <- ggplot(snowmelt, aes(x = rf_short, y = rf_long, colour = rf_long_working)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_vline(xintercept = SNOW_EM_THRESH, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = SNOW_EM_THRESH, linetype = "dashed", colour = "grey40") +
  annotate("text", x = SNOW_EM_THRESH + 0.01, y = 0.98,
           label = sprintf("rf_short = %.2f", SNOW_EM_THRESH),
           hjust = 0, size = 3, colour = "grey40") +
  annotate("text", x = 0.98, y = SNOW_EM_THRESH + 0.01,
           label = sprintf("rf_long = %.2f", SNOW_EM_THRESH),
           hjust = 1, vjust = 0, size = 3, colour = "grey40") +
  scale_colour_manual(
    values = c("Both low (redundant)"    = "#377eb8",
               "rf_long only (adds value)" = "#e41a1c",
               "Other"                   = "#999999")
  ) +
  labs(
    title    = "rf_short vs rf_long for snowmelt-classified peaks",
    subtitle = "Red points are events where rf_long is doing unique work (rf_short > threshold)",
    x        = "rf_short (rain fraction over search_radius window)",
    y        = "rf_long (rain fraction over 3x search_radius window)",
    colour   = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "rf_long_scatter.png"), p1,
       width = 7, height = 6, dpi = 150)
message("Saved: _diagnostics/rf_long_necessity/rf_long_scatter.png")

# --- Plot 2: Marginal density of rf_short for snowmelt peaks ----------------
p2 <- ggplot(snowmelt, aes(x = rf_short, fill = rf_long_working)) +
  geom_histogram(bins = 50, position = "stack", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = SNOW_EM_THRESH, linetype = "dashed", colour = "grey30") +
  scale_fill_manual(
    values = c("Both low (redundant)"      = "#377eb8",
               "rf_long only (adds value)" = "#e41a1c",
               "Other"                     = "#999999")
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    title    = "rf_short distribution for snowmelt-classified peaks",
    subtitle = "Red = events where rf_long uniquely enabled snowmelt classification",
    x        = "rf_short",
    y        = "Count",
    fill     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "rf_long_rfshort_hist.png"), p2,
       width = 7, height = 5, dpi = 150)
message("Saved: _diagnostics/rf_long_necessity/rf_long_rfshort_hist.png")

# --- Plot 3: Seasonal pattern for rf_long-unique events ---------------------
if (nrow(rf_long_unique) > 10) {
  p3 <- rf_long_unique %>%
    mutate(month = factor(month, levels = 1:12,
                          labels = month.abb)) %>%
    ggplot(aes(x = month)) +
    geom_bar(fill = "#e41a1c", colour = "white", linewidth = 0.2) +
    labs(
      title    = "Seasonal distribution of rf_long-unique snowmelt events",
      subtitle = "Peaks where rf_short > 0.05 but rf_long <= 0.05",
      x        = "Month",
      y        = "Count"
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(outdir, "rf_long_seasonal.png"), p3,
         width = 7, height = 4, dpi = 150)
  message("Saved: _diagnostics/rf_long_necessity/rf_long_seasonal.png")
}

# --- Plot 4: Counterfactual comparison (stacked bar by class) ---------------
bind_rows(
  dat_unreg   %>% count(flood_type)         %>% rename(type = flood_type)  %>% mutate(scheme = "Original (rf_long)"),
  dat_counter %>% count(flood_type_counter) %>% rename(type = flood_type_counter) %>% mutate(scheme = "Counterfactual (rf_short)")
) %>%
  mutate(type = factor(type, levels = c("rainfall", "snowmelt", "ROS", "mixed"))) %>%
  ggplot(aes(x = scheme, y = n, fill = type)) +
  geom_col(position = "fill", colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(
    rainfall = "#e41a1c", snowmelt = "#377eb8",
    ROS = "#4daf4a", mixed = "#ff7f00"
  )) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Classification composition: original vs. counterfactual",
    x = NULL, y = "Proportion", fill = "Flood type"
  ) +
  theme_bw(base_size = 11) -> p4

ggsave(file.path(outdir, "rf_long_counterfactual.png"), p4,
       width = 6, height = 5, dpi = 150)
message("Saved: _diagnostics/rf_long_necessity/rf_long_counterfactual.png")

message("\nDiagnostic complete.")