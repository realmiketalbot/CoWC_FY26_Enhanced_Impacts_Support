## lmomco is called via `::` throughout (rather than library(lmomco) at the
## top) so that scripts sourcing this file only for the shared utility
## helpers below are not forced to have lmomco installed.

## ── Poisson–GPD helpers (Proposal Eq. 1.3–1.5 and App. A.1; Beguería 2005) ──
##
## For flood type i with fitting threshold u_i:
##   H_i(y)   GPD CDF of excesses y = x - u_i, fit by L-moments with the
##            lower bound FIXED at 0: lmomco::pargpa(lmoms(y), xi = 0).
##   lambda_i exceedances of u_i per qualifying water year (Script 2b).
##   F_i(x) = exp(-lambda_i * (1 - H_i(x - u_i)))   annual-maximum CDF, x >= u_i
##   T(x)   = 1 / (1 - F(x))                        ANNUAL return period
## Mixture (Eq. 1.5): F_mix(x) = prod_i F_i(x). Below u_i, F_i is floored at
## exp(-lambda_i) (decision D4). Quantiles are reported only where valid:
## per type for x >= u_i, mixture for x >= max_i u_i; otherwise NA.
##
## lmomco convention: pargpa() returns para = c(xi = location, alpha = scale,
## kappa), with kappa = -shape in the Coles (2001) convention the proposal
## uses. CDF and quantile evaluation goes through lmomco's cdfgpa()/quagpa()
## so the sign convention is handled in one place.
##
## The initial-results fit_gpd() called pargpa() WITHOUT xi = 0, so lmomco
## estimated a location from the excesses that the code then discarded
## (plan §2). These helpers replace fit_gpd(), the LP3 helpers and the
## family dispatch; the old code is at tag initial-results-prelim-proposal.

## Interim fitting threshold (plan D-A, stand-in for EQD).
## Rate-based: keep the k largest peaks, k = min(n, max(min_count,
## round(rate * n_years))), the closest analogue to the Yu et al. (2022)
## top-n sample. u is placed midway between the k-th largest peak and the
## next lower DISTINCT peak value (or floor_u, the site's candidate
## threshold, when no lower peak exists). Peaks tied with the k-th value are
## therefore all included, so the returned k (the number of peaks above u)
## can exceed the target slightly when there are ties.
## Returns NULL when fewer than min_count peaks are available.
select_threshold_rate <- function(x, n_years, floor_u, min_count = 20L, rate = 1) {
  x <- sort(x[is.finite(x)], decreasing = TRUE)
  n <- length(x)
  k <- min(n, max(min_count, round(rate * n_years)))
  if (k < min_count) return(NULL)
  lower <- x[x < x[k]]
  next_lower <- if (length(lower) > 0) lower[1] else min(floor_u, x[k])
  u <- (x[k] + next_lower) / 2
  list(u = u, k = sum(x > u), rule = sprintf("rate_%g_per_yr", rate))
}

## Fit the Poisson–GPD to the peaks x above u, over n_years qualifying years.
fit_pgpd <- function(x, u, n_years, min_n = 20L) {
  y <- x[is.finite(x) & x > u] - u
  if (length(y) < min_n || !is.finite(n_years) || n_years <= 0) return(NULL)
  pars <- tryCatch(lmomco::pargpa(lmomco::lmoms(y), xi = 0), error = function(e) NULL)
  if (is.null(pars) || anyNA(pars$para) || any(!is.finite(pars$para))) return(NULL)
  list(
    u      = u,
    n_exc  = length(y),
    n_years = n_years,
    lambda = length(y) / n_years,
    sigma  = unname(pars$para["alpha"]),
    xi     = -unname(pars$para["kappa"]),   # Coles (2001) shape
    para   = pars
  )
}

## Annual-maximum CDF of one Poisson–GPD fit, floored below u (D4).
pgpd_cdf <- function(x, fit) {
  h <- ifelse(x >= fit$u, lmomco::cdfgpa(pmax(x - fit$u, 0), fit$para), 0)
  exp(-fit$lambda * (1 - h))
}

## Quantile at annual return period(s) T; NA where the quantile falls below u
## (i.e., where -log(1 - 1/T) > lambda).
pgpd_quantile <- function(T, fit) {
  p_exc <- -log(1 - 1 / T) / fit$lambda       # = 1 - H(y_T)
  q <- rep(NA_real_, length(T))
  ok <- is.finite(p_exc) & p_exc <= 1 & p_exc > 0
  q[ok] <- fit$u + lmomco::quagpa(1 - p_exc[ok], fit$para)
  q
}

## Mixture annual-maximum CDF across a list of fits (Eq. 1.5).
pgpd_mixture_cdf <- function(x, fits) Reduce(`*`, lapply(fits, function(f) pgpd_cdf(x, f)))

## Mixture quantile(s) by root-finding; NA where the quantile would fall below
## max_i u_i (D4) or the root cannot be bracketed.
pgpd_mixture_quantile <- function(T, fits) {
  lb <- max(vapply(fits, `[[`, numeric(1), "u"))
  vapply(T, function(t) {
    target <- 1 - 1 / t
    if (pgpd_mixture_cdf(lb, fits) >= target) return(NA_real_)
    ub <- max(lb * 2, lb + 1)
    for (i in seq_len(60)) {
      if (pgpd_mixture_cdf(ub, fits) >= target) break
      ub <- lb + (ub - lb) * 2
    }
    if (pgpd_mixture_cdf(ub, fits) < target) return(NA_real_)
    tryCatch(uniroot(function(x) pgpd_mixture_cdf(x, fits) - target,
                     lower = lb, upper = ub, tol = 1e-6 * ub)$root,
             error = function(e) NA_real_)
  }, numeric(1))
}

## ── Empirical annual return periods for exceedances ─────────────────────────
##
## For the m-th largest of the peaks above u over N qualifying years, the
## partial-duration recurrence interval is T_pds = (N + 1) / m (Weibull).
## It is converted to an annual return period with Langbein's relation,
## T = 1 / (1 - exp(-1 / T_pds)), the same Poisson link used by F_i(x).
empirical_annual_rp <- function(x, u, n_years) {
  x <- sort(x[is.finite(x) & x > u], decreasing = TRUE)
  m <- seq_along(x)
  t_pds <- (n_years + 1) / m
  data.frame(peak_flow = x, rank = m, ri = 1 / (1 - exp(-1 / t_pds)))
}

## ── Flood-type remap (Script 7, report_figures.R, App) ──────────────────────
##
## Collapses the 3-class rain-fraction labels to 2 classes: snowmelt and
## rain_on_snow become snow_influenced. Works from flood_type_raw, the only
## label column after M1, and overwrites flood_type with the 2-class label.
remap_2class <- function(peaks) {
  snow_influenced_types <- c("rain_on_snow", "snowmelt")
  raw <- as.character(peaks$flood_type_raw)
  peaks$flood_type <- factor(
    ifelse(raw %in% snow_influenced_types, "snow_influenced",
           ifelse(raw == "rainfall", "rainfall", "unclassified")),
    levels = c("rainfall", "snow_influenced", "unclassified")
  )
  peaks
}

## ── Shared utility helpers (Scripts 6a, 6b, 8, 9) ─────────────────────────────

## Water year: WY Y runs Oct 1 of (Y-1) through Sep 30 of Y.
water_year <- function(d) {
  yr <- as.integer(format(d, "%Y"))
  mo <- as.integer(format(d, "%m"))
  ifelse(mo >= 10L, yr + 1L, yr)
}

## Mann-Kendall trend test + Theil-Sen slope, NA-safe (drops NA before
## testing; returns NA fields rather than erroring if a test fails or too
## few points remain). Callers apply their own eligibility threshold to
## n_years before trusting the result.
mk_sens_test <- function(vals) {
  vals_ok <- vals[!is.na(vals)]
  mk_res  <- tryCatch(Kendall::MannKendall(vals_ok), error = function(e) NULL)
  ss_res  <- tryCatch(trend::sens.slope(vals_ok),    error = function(e) NULL)
  list(
    n_years   = length(vals_ok),
    tau       = if (!is.null(mk_res)) unname(mk_res$tau[1])    else NA_real_,
    mk_pvalue = if (!is.null(mk_res)) unname(mk_res$sl[1])     else NA_real_,
    sen_slope = if (!is.null(ss_res)) unname(ss_res$estimates) else NA_real_
  )
}

## Pettitt change-point test, NA-safe. pettitt_index is the 1-based index of
## the most probable change point WITHIN THE NA-STRIPPED SERIES — to map it
## back to a label (e.g. a water year), apply the same `vals[!is.na(vals)]`
## filter to the label vector before indexing, e.g.
##   labels[!is.na(vals)][pettitt_index]
pettitt_shift_test <- function(vals) {
  vals_ok <- vals[!is.na(vals)]
  pt_res  <- tryCatch(trend::pettitt.test(vals_ok), error = function(e) NULL)
  list(
    pettitt_pvalue = if (!is.null(pt_res)) unname(pt_res$p.value)     else NA_real_,
    pettitt_index  = if (!is.null(pt_res)) unname(pt_res$estimate[1]) else NA_integer_
  )
}