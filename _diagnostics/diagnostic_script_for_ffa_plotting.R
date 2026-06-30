#!/usr/bin/env Rscript
# =============================================================================
# diagnose-ffa-site.R
#
# Single-site flood frequency analysis diagnostic. Mirrors app.R exactly:
# same fit_gpd(), gpd_sf(), gpd_annual_ep(), mixture_annual_ep(), and
# empirical_ri() logic, including the fixed-mu / pot_threshold_cfs approach.
#
# Usage:
#   source("diagnose-ffa-site.R")                     # uses SITE_ID below
#   Rscript diagnose-ffa-site.R 09085100              # override from CLI
#   Rscript diagnose-ffa-site.R 09085100 plot.png     # also save plot
# =============================================================================

library(dplyr)
library(lmomco)
library(ggplot2)
library(scales)

## ── Config ───────────────────────────────────────────────────────────────────
SITE_ID   <- "09085100"
DATA_FILE <- "_data/pot_peaks_classified.rds"
OUT_FILE  <- NULL   # set to a path e.g. "ffa_09085100.png" to save plot

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) SITE_ID  <- args[1]
if (length(args) >= 2) OUT_FILE <- args[2]

## ── GPD helpers (mirrored exactly from app.R) ─────────────────────────────────
##
## Mixture: F_mix(x) = prod_i F_i(x)  [Yu et al. 2022, Eq. 3]
## RI axis is PoT recurrence interval: RI = 1 / (1 - F(x))
## Empirical plotting positions: within-type Weibull, RI = (n+1)/rank

fit_gpd <- function(x, dates, mu_threshold, min_n = 5) {
  keep  <- is.finite(x) & x > 0
  x     <- x[keep]; dates <- dates[keep]
  if (length(x) < min_n) return(NULL)
  n_years <- as.numeric(diff(range(dates))) / 365.25
  if (n_years <= 0) return(NULL)
  y <- sort(x - mu_threshold); y <- y[y >= 0]
  if (length(y) < min_n) return(NULL)
  pars <- pargpa(lmoms(y))
  list(mu = mu_threshold, sigma = pars$para["alpha"],
       xi = pars$para["kappa"] * -1,
       n = length(x), n_years = n_years)
}

gpd_sf <- function(x, mu, sigma, xi) {
  z <- (x - mu) / sigma; inner <- 1 + xi * z
  if (abs(xi) < 1e-8) pmax(0, exp(-z))
  else pmax(0, ifelse(inner <= 0, 0, inner^(-1 / xi)))
}

gpd_cdf  <- function(x, fit) 1 - gpd_sf(x, fit$mu, fit$sigma, fit$xi)

gpd_quantile <- function(p, fit) {
  if (abs(fit$xi) < 1e-8) fit$mu - fit$sigma * log(p)
  else fit$mu + fit$sigma / fit$xi * (p^(-fit$xi) - 1)
}

gpd_curve <- function(fit, label, ri = c(seq(1.1, 10, 0.2), seq(11, 500, 1))) {
  data.frame(ri = ri, q = pmax(0, sapply(1/ri, gpd_quantile, fit = fit)),
             label = label, stringsAsFactors = FALSE)
}

mixture_sf <- function(x, fits)
  1 - Reduce("*", lapply(fits, function(f) gpd_cdf(x, f)))

mixture_curve <- function(fits, ri = c(seq(1.1, 10, 0.2), seq(11, 500, 1))) {
  min_target <- min(1 / ri)
  ub <- max(sapply(fits, function(f) gpd_quantile(min_target, f)), na.rm = TRUE) * 2
  ub <- max(ub, 1)
  for (i in seq_len(30)) {
    if (isTRUE(mixture_sf(ub, fits) < min_target)) break
    ub <- ub * 2
  }
  q <- sapply(ri, function(r) {
    tgt  <- 1 / r
    lb   <- min(sapply(fits, function(f) f$mu), na.rm = TRUE)
    f_lb <- mixture_sf(lb, fits) - tgt
    f_ub <- mixture_sf(ub, fits) - tgt
    if (is.na(f_lb) || is.na(f_ub) || sign(f_lb) == sign(f_ub)) return(NA_real_)
    tryCatch(uniroot(function(x) mixture_sf(x, fits) - tgt,
                     lower = lb, upper = ub, tol = 0.1)$root,
             error = function(e) NA_real_)
  })
  data.frame(ri = ri, q = pmax(0, q), label = "Mixture", stringsAsFactors = FALSE)
}

empirical_ri <- function(peaks) {
  peaks %>%
    filter(!is.na(peak_flow), peak_flow > 0) %>%
    group_by(flood_type) %>% arrange(desc(peak_flow)) %>%
    mutate(rank = seq_len(n()), n = n(),
           ep = rank / (n + 1), ri = (n + 1) / rank) %>%
    ungroup()
}

## ── Load data ─────────────────────────────────────────────────────────────────
pot   <- readRDS(DATA_FILE)
peaks <- pot %>%
  filter(site_no == SITE_ID, flood_type != "unclassified",
         is.finite(peak_flow), peak_flow > 0)

if (nrow(peaks) == 0) stop(sprintf("No classified peaks found for site %s", SITE_ID))

mu_site <- peaks$pot_threshold_cfs[1]
n_years <- as.numeric(diff(range(peaks$peak_date))) / 365.25

cat(sprintf("\n══ Site: %s ══\n", SITE_ID))
cat(sprintf("Record:    %s to %s (%.1f yr)\n",
            min(peaks$peak_date), max(peaks$peak_date), n_years))
cat(sprintf("Threshold: %.1f cfs  (pot_threshold_cfs used as GPD mu)\n\n", mu_site))

## ── Fit per-type and single-sample GPDs ──────────────────────────────────────
type_names <- c("rainfall", "snowmelt", "rain_on_snow")
fits <- list()

for (typ in type_names) {
  sub <- peaks %>% filter(flood_type == typ)
  f   <- fit_gpd(sub$peak_flow, sub$peak_date, mu_threshold = mu_site)
  if (!is.null(f)) fits[[typ]] <- f
}

f_ss <- fit_gpd(peaks$peak_flow, peaks$peak_date, mu_threshold = mu_site)
if (!is.null(f_ss)) fits[["single_sample"]] <- f_ss

## ── Parameter summary ─────────────────────────────────────────────────────────
cat("── GPD parameters ──\n")
cat(sprintf("%-15s  %4s  %8s  %7s  %7s  %8s\n",
            "Type", "n", "T (yr)", "mu", "sigma", "xi"))
cat(strrep("-", 60), "\n")
for (nm in names(fits)) {
  f <- fits[[nm]]
  cat(sprintf("%-15s  %4d  %7.1f  %7.1f  %7.1f  %+7.3f\n",
              nm, f$n, f$n_years, f$mu, f$sigma, f$xi))
}

## ── Mixture diagnostics ───────────────────────────────────────────────────────
type_fits <- fits[names(fits) != "single_sample"]

cat("\n── Per-type CDF values at selected discharges ──\n")
x_probe <- sort(unique(round(c(
  mu_site,
  sapply(type_fits, function(f) gpd_quantile(0.5, f)),
  sapply(type_fits, function(f) gpd_quantile(0.1, f)),
  sapply(type_fits, function(f) gpd_quantile(0.01, f))
))))
x_probe <- x_probe[is.finite(x_probe) & x_probe > 0]

type_nms <- names(type_fits)
hdr <- paste(sprintf("%-12s", c(paste0("F_", type_nms), "F_mix", "mix_RI")),
             collapse = " ")
cat(sprintf("  %10s  %s\n", "x (cfs)", hdr))
for (x in x_probe) {
  cdfs    <- sapply(type_fits, function(f) gpd_cdf(x, f))
  f_mix   <- prod(cdfs)
  mix_sf  <- 1 - f_mix
  mix_ri  <- if (!is.na(mix_sf) && mix_sf > 0) 1 / mix_sf else Inf
  vals    <- c(cdfs, f_mix, mix_ri)
  fmts    <- c(rep("%-12.4f", length(cdfs) + 1), "%-12.1f")
  cat(sprintf("  %10.0f  %s\n", x,
              paste(mapply(sprintf, fmts, vals), collapse = " ")))
}

cat("\n── Mixture SF at lower bound ──\n")
lb     <- min(sapply(type_fits, function(f) f$mu), na.rm = TRUE)
sf_lb  <- mixture_sf(lb, type_fits)
cat(sprintf("  x = %.1f cfs (common mu):  mix_SF = %.5f  =>  min RI = %.2f yr\n",
            lb, sf_lb, 1 / sf_lb))
for (nm in type_nms) {
  cat(sprintf("    %-15s  F(mu) = %.5f  [sigma=%.1f, xi=%+.3f]\n",
              nm, gpd_cdf(lb, type_fits[[nm]]),
              type_fits[[nm]]$sigma, type_fits[[nm]]$xi))
}

cat("\n── RI grid quantiles ──\n")
for (r in c(2, 5, 10, 25, 50, 100, 200, 500)) {
  tgt <- 1 / r
  q <- tryCatch({
    f_lb <- mixture_sf(lb, type_fits) - tgt
    f_ub <- mixture_sf(lb + 1e8, type_fits) - tgt
    if (sign(f_lb) == sign(f_ub)) NA_real_
    else uniroot(function(x) mixture_sf(x, type_fits) - tgt,
                 lower = lb, upper = lb + 1e8, tol = 0.1)$root
  }, error = function(e) NA_real_)
  cat(sprintf("  RI = %4d yr:  Q = %s cfs\n",
              r, if (is.na(q)) "NA" else formatC(q, format = "f", digits = 0)))
}

## ── Plot ─────────────────────────────────────────────────────────────────────
ri_grid <- c(seq(1.1, 10, 0.2), seq(11, 500, 1))

type_colors <- c(
  rainfall      = "#e41a1c",
  snowmelt      = "#377eb8",
  rain_on_snow  = "#4daf4a",
  single_sample = "#888888",
  Mixture       = "black"
)
type_labels <- c(
  rainfall      = "Rainfall",
  snowmelt      = "Snowmelt",
  rain_on_snow  = "Rain-on-Snow",
  single_sample = "Single Sample",
  Mixture       = "Mixture"
)

curves <- bind_rows(lapply(names(type_fits), function(nm) {
  gpd_curve(type_fits[[nm]], type_labels[nm], ri_grid) %>%
    mutate(color = type_colors[nm], lty = "dashed")
}))

if (length(type_fits) >= 2) {
  mc <- mixture_curve(type_fits, ri_grid) %>%
    mutate(color = "black", lty = "solid")
  curves <- bind_rows(curves, mc)
}

if (!is.null(fits[["single_sample"]])) {
  ss <- gpd_curve(fits[["single_sample"]], "Single Sample", ri_grid) %>%
    mutate(color = type_colors["single_sample"], lty = "solid")
  curves <- bind_rows(curves, ss)
}

emp <- empirical_ri(peaks %>% filter(flood_type %in% names(type_fits))) %>%
  mutate(label = type_labels[flood_type],
         color = type_colors[flood_type])

present <- unique(curves$label)
col_vec <- setNames(unname(curves$color[match(present, curves$label)]), present)
lty_vec <- setNames(unname(curves$lty[match(present, curves$label)]),   present)

## Pretty y-axis: log-spaced breaks formatted with commas
y_vals   <- c(curves$q, emp$peak_flow)
y_range  <- range(y_vals[is.finite(y_vals) & y_vals > 0], na.rm = TRUE)
y_breaks <- scales::breaks_log(n = 6)(y_range)
y_labels <- scales::label_comma(accuracy = 1)(y_breaks)

p <- ggplot() +
  geom_line(data = curves,
            aes(x = ri, y = q, colour = label, linetype = label),
            linewidth = 0.9, na.rm = TRUE) +
  geom_point(data = emp,
             aes(x = ri, y = peak_flow, colour = label),
             size = 2, alpha = 0.75) +
  scale_x_log10(breaks = c(2, 5, 10, 25, 50, 100, 200, 500),
                labels  = c("2","5","10","25","50","100","200","500")) +
  scale_y_log10(breaks = y_breaks, labels = y_labels) +
  coord_cartesian(xlim = c(2, 500)) +
  scale_colour_manual(values = col_vec, name = NULL) +
  scale_linetype_manual(values = lty_vec, name = NULL) +
  labs(
    title   = sprintf("%s — FFA", SITE_ID),
    x       = "Return Period (years)",
    y       = "Peak Discharge (cfs)",
    caption = sprintf(
      "mu = %.1f cfs (pot_threshold_cfs) | GPD via L-moments | Yu et al. 2022",
      mu_site)
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

if (!is.null(OUT_FILE)) {
  ggsave(OUT_FILE, p, width = 8, height = 5.5, dpi = 150)
  cat(sprintf("\nPlot saved to %s\n", OUT_FILE))
} else {
  print(p)
}