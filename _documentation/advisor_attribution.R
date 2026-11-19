#!/usr/bin/env Rscript
# ============================================================================
# What changed the results? Step-by-step attribution, old -> new (INTERNAL)
# ============================================================================
# Re-fits the rain-fraction (3-class) FFA under a sequence of variants, each
# turning on one correction, and tracks the mixture / single-sample ratio.
#
#   Data "old": the initial-results peaks (tag initial-results-prelim-proposal):
#     top-n samples per type (n = max(20, #water years)) and single sample
#     (m = max(25, #water years)); location mu = site candidate threshold.
#   Data "new": M1 pool (fit_eligible peaks in qualifying water years, 215
#     sites, new classification), rate-based sample k = max(20, #qual. years).
#
#   Switches
#     loc_fix  FALSE: pargpa() on excesses estimates a location that is then
#                     discarded (the initial-results defect).
#              TRUE : pargpa(..., xi = 0), lower bound fixed at 0.
#     thresh   "site": u = site candidate threshold (median daily max), far
#                     below the smallest sampled peak.
#              "sample": u just below the smallest sampled peak (D-A rule).
#     annual   FALSE: per-event CDFs, RI = 1 / (1 - F)  (PoT recurrence
#                     interval); mixture S = 1 - prod F_i.
#              TRUE : Poisson–GPD, F = exp(-lambda (1 - H)), annual T.
#
# Variant "S0" must reproduce the old outputs; "S4" must reproduce the new
# outputs. Both are checked below.
#
# Run from the repo root:  Rscript _documentation/advisor_attribution.R
# Output: _documentation/advisor_old_vs_new/10_attribution*.png,
#         _documentation/advisor_old_vs_new/attribution_summary.csv
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
source("functions.R")

TAG     <- "initial-results-prelim-proposal"
fig_dir <- "_documentation/advisor_old_vs_new"
RPS     <- c(2, 5, 10, 25, 50, 100, 200, 500)
TYPES   <- c("rainfall", "snowmelt", "rain_on_snow")

read_tag <- function(path) {
  tmp <- tempfile(fileext = ".rds")
  if (system2("git", c("show", sprintf("%s:%s", TAG, path)), stdout = tmp) != 0)
    stop("Could not read ", path, " at tag ", TAG)
  readRDS(tmp)
}

# ---- Samples ----------------------------------------------------------------
# For each site: mu (candidate threshold), n_years, and for each curve
# (type or single_sample): the full eligible pool and the sample size k.
old_peaks <- read_tag("_data/pot_peaks_classified_frac_thresh_3class.rds")
old_q     <- read_tag("_data/ffa_quantiles_gpd.rds") %>% filter(method == "frac_thresh_3class")
new_q     <- readRDS("_data/ffa_quantiles_pgpd.rds") %>% filter(method == "frac_thresh_3class")

old_sets <- old_peaks %>% filter(is.finite(peak_flow), peak_flow > 0) %>%
  group_split(site_no) %>%
  map(function(d) {
    elig <- d %>% filter(coalesce(data_sufficient & !is.na(rain_frac) & flux_sufficient, FALSE))
    curves <- list()
    for (t in TYPES) {
      smp <- d$peak_flow[as.character(d$flood_type) == t & d$in_type_sample]
      if (length(smp) >= 20)
        curves[[t]] <- list(pool = elig$peak_flow[as.character(elig$flood_type_raw) == t],
                            sample = smp, min_n = 20)
    }
    smp <- d$peak_flow[d$in_single_sample]
    if (length(smp) >= 25) curves$single_sample <- list(pool = elig$peak_flow, sample = smp, min_n = 25)
    list(site_no = d$site_no[1], mu = d$pot_threshold_cfs[1],
         n_years = n_distinct(d$water_year), curves = curves)
  })

new_peaks <- readRDS("_data/pot_peaks_classified_frac_thresh_3class.rds")
qual      <- readRDS("_data/site_qualifying_years.rds") %>% filter(selected, qualifies)
new_sets <- new_peaks %>% filter(fit_eligible, is.finite(peak_flow), peak_flow > 0) %>%
  semi_join(qual, by = c("site_no", "water_year")) %>%
  group_split(site_no) %>%
  map(function(d) {
    ny <- sum(qual$site_no == d$site_no[1])
    k_of <- function(x) min(length(x), max(20L, round(ny)))
    curves <- list()
    for (t in TYPES) {
      x <- d$peak_flow[as.character(d$flood_type_raw) == t]
      if (length(x) >= 20) curves[[t]] <- list(pool = x, sample = sort(x, TRUE)[seq_len(k_of(x))], min_n = 20)
    }
    if (nrow(d) >= 20) curves$single_sample <- list(pool = d$peak_flow,
                                                    sample = sort(d$peak_flow, TRUE)[seq_len(k_of(d$peak_flow))],
                                                    min_n = 20)
    list(site_no = d$site_no[1], mu = d$candidate_threshold_cfs[1], n_years = ny, curves = curves)
  })

# ---- One fit ----------------------------------------------------------------
fit_variant <- function(cv, mu, n_years, loc_fix, thresh) {
  smp <- sort(cv$sample, decreasing = TRUE)
  if (thresh == "site") {
    u <- mu
    y <- smp - u; y <- y[y >= 0]
  } else {
    k <- length(smp)
    lower <- cv$pool[cv$pool < smp[k]]
    nxt <- if (length(lower)) max(lower) else min(mu, smp[k])
    u <- (smp[k] + nxt) / 2
    y <- cv$pool[cv$pool > u] - u   # includes peaks tied with the k-th, as select_threshold_rate()
  }
  if (length(y) < cv$min_n) return(NULL)
  pars <- tryCatch(if (loc_fix) lmomco::pargpa(lmomco::lmoms(y), xi = 0)
                   else lmomco::pargpa(lmomco::lmoms(y)), error = function(e) NULL)
  if (is.null(pars) || anyNA(pars$para)) return(NULL)
  pars$para["xi"] <- 0          # old code discarded the estimated location
  list(u = u, para = pars, lambda = length(y) / n_years, n = length(y))
}

# Per-event (PoT RI) quantiles, as in the initial results
pe_quantile <- function(T, f) f$u + lmomco::quagpa(1 - 1 / T, f$para)
pe_mixture_quantile <- function(T, fits) {
  sf  <- function(x) 1 - Reduce(`*`, lapply(fits, function(f)
    ifelse(x >= f$u, lmomco::cdfgpa(pmax(x - f$u, 0), f$para), 0)))
  lb  <- min(vapply(fits, `[[`, 1, "u"))
  ub  <- max(vapply(fits, function(f) pe_quantile(max(T), f), 1)) * 2
  vapply(T, function(t) {
    if (sign(sf(lb) - 1 / t) == sign(sf(ub) - 1 / t)) return(NA_real_)
    uniroot(function(x) sf(x) - 1 / t, c(lb, ub), tol = 1e-6 * ub)$root
  }, 1)
}

run_variant <- function(sets, loc_fix, thresh, annual) {
  map_dfr(sets, function(s) {
    fits <- compact(map(s$curves, fit_variant, mu = s$mu, n_years = s$n_years,
                        loc_fix = loc_fix, thresh = thresh))
    tf <- fits[intersect(names(fits), TYPES)]
    if (is.null(fits$single_sample) || length(tf) < 2) return(NULL)
    if (annual) {
      ss <- pgpd_quantile(RPS, fits$single_sample)
      mx <- pgpd_mixture_quantile(RPS, tf)
    } else {
      ss <- pe_quantile(RPS, fits$single_sample)
      mx <- pe_mixture_quantile(RPS, tf)
    }
    tibble(site_no = s$site_no, return_period = RPS, mixture = mx, single_sample = ss)
  })
}

steps <- tribble(
  ~step, ~label,                                   ~data, ~loc_fix, ~thresh,  ~annual,
  "S0",  "Old method",                             "old", FALSE,    "site",   FALSE,
  "S1",  "+ fix GPD lower bound",                  "old", TRUE,     "site",   FALSE,
  "S2",  "+ threshold at sampled peaks",           "old", TRUE,     "sample", FALSE,
  "S3",  "+ annual return periods (Poisson)",      "old", TRUE,     "sample", TRUE,
  "S4",  "+ new stations, years, classification",  "new", TRUE,     "sample", TRUE,
  # one change at a time, from the old method
  "A1",  "Old + GPD lower bound only",             "old", TRUE,     "site",   FALSE,
  "A2",  "Old + threshold only",                   "old", FALSE,    "sample", FALSE,
  "A3",  "Old + annual return periods only",       "old", FALSE,    "site",   TRUE
)

message("Fitting variants...")
res <- pmap_dfr(steps, function(step, label, data, loc_fix, thresh, annual) {
  message("  ", step, ": ", label)
  run_variant(if (data == "old") old_sets else new_sets, loc_fix, thresh, annual) %>%
    mutate(step = step, label = label)
})

# ---- Checks: S0 reproduces the old outputs, S4 the new ----------------------
chk <- function(step_id, ref) {
  ref_w <- ref %>% filter(curve_type %in% c("mixture", "single_sample")) %>%
    select(site_no, return_period, curve_type, discharge_cfs) %>%
    pivot_wider(names_from = curve_type, values_from = discharge_cfs)
  j <- res %>% filter(step == step_id) %>% inner_join(ref_w, by = c("site_no", "return_period"), suffix = c("", ".ref"))
  rel <- c(abs(j$single_sample / j$single_sample.ref - 1), abs(j$mixture / j$mixture.ref - 1))
  message(sprintf("  %s vs reference: %d site-RPs, max rel. diff %.2e, %.1f%% within 0.1%%",
                  step_id, nrow(j), max(rel, na.rm = TRUE), 100 * mean(rel < 1e-3, na.rm = TRUE)))
}
message("Checks:")
chk("S0", old_q)
chk("S4", new_q)

# ---- Summary on the sites common to old and new -----------------------------
common <- intersect(map_chr(old_sets, "site_no"), map_chr(new_sets, "site_no"))
summ <- res %>% filter(site_no %in% common, is.finite(mixture), is.finite(single_sample)) %>%
  mutate(ratio = mixture / single_sample) %>%
  group_by(step, label, return_period) %>%
  summarise(n_sites = n(), median = median(ratio), q25 = quantile(ratio, .25), q75 = quantile(ratio, .75),
            share_gt1 = mean(ratio > 1), .groups = "drop")
write_csv(summ, file.path(fig_dir, "attribution_summary.csv"))

summ %>% filter(return_period %in% c(2, 10, 100)) %>%
  transmute(step, label, return_period, txt = sprintf("%.2f (%.0f%%)", median, 100 * share_gt1)) %>%
  pivot_wider(names_from = return_period, values_from = txt, names_prefix = "T") %>%
  print(width = Inf)

# ---- Figure -----------------------------------------------------------------
cum <- summ %>% filter(str_starts(step, "S"), return_period %in% c(2, 10, 100)) %>%
  mutate(label = factor(label, levels = steps$label[1:5]),
         rp = factor(sprintf("%d-year", as.integer(return_period)), levels = c("2-year", "10-year", "100-year")))
p1 <- ggplot(cum, aes(label, median, colour = rp, group = rp)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.15, position = position_dodge(0.4), alpha = 0.6) +
  geom_line(position = position_dodge(0.4), linewidth = 0.9) +
  geom_point(position = position_dodge(0.4), size = 2.4) +
  scale_colour_manual(values = c("2-year" = "#a6cee3", "10-year" = "#1f78b4", "100-year" = "#08306b"), name = NULL) +
  scale_x_discrete(labels = function(x) str_wrap(x, 16)) +
  labs(x = NULL, y = "Median mixture / single sample\n(bars: IQR)") +
  theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank(), legend.position = "top")
p2 <- ggplot(cum, aes(label, share_gt1, colour = rp, group = rp)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
  scale_colour_manual(values = c("2-year" = "#a6cee3", "10-year" = "#1f78b4", "100-year" = "#08306b"), guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_discrete(labels = function(x) str_wrap(x, 16)) +
  labs(x = NULL, y = "Sites where\nmixture > single sample") +
  theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank())
p <- (p1 / p2) + plot_layout(heights = c(3, 2)) +
  plot_annotation(title = "What changed the mixture / single-sample result?",
                  subtitle = sprintf("Corrections applied one at a time, cumulatively; Adaptive (3-class); %d common sites",
                                     length(common)),
                  theme = theme(plot.title = element_text(face = "bold", size = 16),
                                plot.subtitle = element_text(colour = "grey40")))
ggsave(file.path(fig_dir, "10_attribution_steps.png"), p, width = 10, height = 8, dpi = 200)
message("Saved 10_attribution_steps.png and attribution_summary.csv")

# ============================================================================
# Mechanism at one site: the discarded GPD location (initial-results defect)
# ============================================================================
# Same old samples and per-event RI axis in both panels; only the location
# treatment differs. pargpa() on excesses above the site threshold mu
# estimates a location close to (smallest sampled peak - mu); the old code
# dropped it, shifting each curve DOWN by that amount. The single sample
# (top peaks of all types) has the largest gap, so it drops most and the
# mixture / single ratio is inflated, most at short RIs.
MECH_SITE <- "09115500"   # Tomichi Creek at Sargents (poster showcase site)
s_old <- old_sets[[which(map_chr(old_sets, "site_no") == MECH_SITE)]]
RI_D  <- c(seq(1.1, 2, 0.05), seq(2.2, 10, 0.2), seq(11, 500, 1))
lab_of <- c(rainfall = "Rainfall", snowmelt = "Snowmelt", rain_on_snow = "Rain-on-Snow",
            single_sample = "Single Sample", mixture = "Mixture")
col_of <- c("Rainfall" = "#b2182b", "Snowmelt" = "#2166ac", "Rain-on-Snow" = "#9970ab",
            "Single Sample" = "grey50", "Mixture" = "black")

mech <- map_dfr(c(FALSE, TRUE), function(keep_loc) {
  fits <- imap(s_old$curves, function(cv, nm) {
    y <- cv$sample - s_old$mu; y <- y[y >= 0]
    p <- lmomco::pargpa(lmomco::lmoms(y))
    loc <- unname(p$para["xi"])
    if (!keep_loc) p$para["xi"] <- 0
    list(u = s_old$mu, para = p, loc = loc)
  })
  tf <- fits[TYPES]
  curves <- imap_dfr(fits, ~ tibble(curve = .y, ri = RI_D, q = pe_quantile(RI_D, .x))) %>%
    bind_rows(tibble(curve = "mixture", ri = RI_D, q = pe_mixture_quantile(RI_D, tf)))
  curves %>% mutate(panel = if (keep_loc) "Location kept" else "As fitted (location discarded)",
                    loc = map_dbl(curve, ~ if (.x == "mixture") NA_real_ else fits[[.x]]$loc))
})
mech <- mech %>% mutate(lab = factor(lab_of[curve], levels = names(col_of)),
                        panel = factor(panel, levels = c("As fitted (location discarded)", "Location kept")))
pts <- imap_dfr(s_old$curves, function(cv, nm) {
  x <- sort(cv$sample, decreasing = TRUE)
  tibble(curve = nm, ri = (length(x) + 1) / seq_along(x), q = x)
}) %>% mutate(lab = factor(lab_of[curve], levels = names(col_of)))
r_at <- function(pn, T) {
  m <- mech %>% filter(panel == pn)
  approx(m$ri[m$curve == "mixture"], m$q[m$curve == "mixture"], T)$y /
    approx(m$ri[m$curve == "single_sample"], m$q[m$curve == "single_sample"], T)$y
}
loc_txt <- mech %>% filter(curve != "mixture") %>% distinct(curve, loc) %>%
  mutate(txt = sprintf("%s %.0f", lab_of[curve], loc)) %>% pull(txt) %>% paste(collapse = ", ")
strip_lab <- setNames(sprintf("%s\nmixture / single: %.2f at RI 2, %.2f at RI 100", levels(mech$panel),
                              c(r_at(levels(mech$panel)[1], 2), r_at(levels(mech$panel)[2], 2)),
                              c(r_at(levels(mech$panel)[1], 100), r_at(levels(mech$panel)[2], 100))),
                      levels(mech$panel))
pm <- ggplot() +
  geom_point(data = bind_rows(pts %>% mutate(panel = levels(mech$panel)[1]),
                              pts %>% mutate(panel = levels(mech$panel)[2])) %>%
               mutate(panel = factor(panel, levels = levels(mech$panel))),
             aes(ri, q, colour = lab), size = 1.6, alpha = 0.7) +
  geom_line(data = mech %>% filter(is.finite(q), q > 0),
            aes(ri, q, colour = lab, linetype = lab %in% c("Mixture", "Single Sample")), linewidth = 1) +
  facet_wrap(~ panel, labeller = as_labeller(strip_lab)) +
  scale_colour_manual(values = col_of, name = NULL) +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dashed"), guide = "none") +
  scale_x_log10(breaks = c(2, 5, 10, 25, 50, 100, 200, 500), name = "Recurrence interval, initial-results definition (years)") +
  scale_y_log10(labels = scales::label_comma(), name = "Peak discharge (cfs)") +
  labs(title = sprintf("Why the initial results inflated the mixture: Station %s (Tomichi Creek)", MECH_SITE),
       subtitle = sprintf("Old samples and old RI axis in both panels. Site threshold mu = %.0f cfs. Location discarded (cfs): %s",
                          s_old$mu, loc_txt),
       caption = "Points: sampled peaks at within-sample Weibull positions, as plotted in the initial results.") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = "grey40", size = 10.5),
        strip.text = element_text(face = "bold"))
ggsave(file.path(fig_dir, "11_mechanism_discarded_location.png"), pm, width = 13, height = 7, dpi = 200)
message("Saved 11_mechanism_discarded_location.png")
