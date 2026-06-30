library(lmomco)

## ── GPD helpers (Yu et al. 2022 methodology) ─────────────────────────────────
##
## RI axis is the PoT recurrence interval: RI = 1 / EP where EP = 1 - F(x).
## This matches Yu et al. Fig 2c/2d exactly (confirmed in their D&C section 5).
## Their total PoT rate ≈ 1 event/yr so PoT RI ≈ annual RI numerically, but
## the underlying quantity is PoT EP, not annual exceedance probability.
##
## Mixture: F_mix(x) = prod_i F_i(x)  [Yu et al. Eq. 3, CDF product]
## Per-type RI curve: RI = 1 / (1 - F_i(x)) = 1 / gpd_sf(x)
## Empirical plotting positions: within-type Weibull, RI = (n+1)/rank

fit_gpd <- function(x, dates, mu_threshold, min_n = 5) {
  keep  <- is.finite(x) & x > 0
  x     <- x[keep]
  dates <- dates[keep]
  if (length(x) < min_n) return(NULL)
  n_years <- as.numeric(diff(range(dates))) / 365.25
  if (n_years <= 0) return(NULL)

  ## Fix mu at the site-wide PoT threshold, following Yu et al. 2022.
  ## Shift exceedances to zero-origin; fit sigma and xi only.
  mu <- mu_threshold
  y  <- sort(x - mu); y <- y[y >= 0]
  if (length(y) < min_n) return(NULL)

  pars <- pargpa(lmoms(y))
  if (is.null(pars$para) || anyNA(pars$para) || any(!is.finite(pars$para))) return(NULL)
  list(
    family  = "gpd",
    mu      = mu,
    sigma   = pars$para["alpha"],
    xi      = pars$para["kappa"] * -1,   # lmomco: kappa = -xi
    n       = length(x),
    n_years = n_years
  )
}

## GPD survival function: P(X > x | X > mu)
gpd_sf <- function(x, mu, sigma, xi) {
  z     <- (x - mu) / sigma
  inner <- 1 + xi * z
  if (abs(xi) < 1e-8) pmax(0, exp(-z))
  else pmax(0, ifelse(inner <= 0, 0, inner^(-1 / xi)))
}

## GPD CDF: F(x) = 1 - S(x)
gpd_cdf <- function(x, fit) 1 - gpd_sf(x, fit$mu, fit$sigma, fit$xi)

## Invert GPD SF: given PoT exceedance probability p, return quantile
gpd_quantile <- function(p, fit) {
  if (abs(fit$xi) < 1e-8) fit$mu - fit$sigma * log(p)
  else fit$mu + fit$sigma / fit$xi * (p^(-fit$xi) - 1)
}

## PoT RI curve for a single GPD fit: RI = 1 / gpd_sf(x)
gpd_curve <- function(fit, label, color,
                      ri = c(seq(1.1, 10, 0.2), seq(11, 500, 1))) {
  ep <- 1 / ri
  q  <- pmax(0, sapply(ep, gpd_quantile, fit = fit))
  data.frame(ri = ri, ep = ep, q = q, label = label, color = color,
             stringsAsFactors = FALSE)
}

## Mixture SF: S_mix(x) = 1 - prod_i F_i(x)  [Yu et al. Eq. 3]
mixture_sf <- function(x, fits)
  1 - Reduce("*", lapply(fits, function(f) gpd_cdf(x, f)))

## Invert mixture SF to get PoT RI curve by uniroot.
mixture_curve <- function(fits, ri = c(seq(1.1, 10, 0.2), seq(11, 500, 1))) {
  min_target <- min(1 / ri)
  ub <- max(sapply(fits, function(f) gpd_quantile(min_target, f)),
            na.rm = TRUE) * 2
  ub <- max(ub, 1)
  for (i in seq_len(30)) {
    if (isTRUE(mixture_sf(ub, fits) < min_target)) break
    ub <- ub * 2
  }
  q <- sapply(ri, function(r) {
    target <- 1 / r
    lb     <- min(sapply(fits, function(f) f$mu), na.rm = TRUE)
    f_lb   <- mixture_sf(lb, fits) - target
    f_ub   <- mixture_sf(ub, fits) - target
    if (is.na(f_lb) || is.na(f_ub) || sign(f_lb) == sign(f_ub)) return(NA_real_)
    tryCatch(
      uniroot(function(x) mixture_sf(x, fits) - target,
              lower = lb, upper = ub, tol = 0.1)$root,
      error = function(e) NA_real_
    )
  })
  data.frame(ri = ri, ep = 1 / ri, q = pmax(0, q),
             label = "Mixture", color = "black",
             stringsAsFactors = FALSE)
}

## ── LP3 helpers (Yu et al. 2022, Equations S1–S3) ────────────────────────────
##
## LP3 is fit to log(x) via L-moments, following Yu et al. 2022 and
## Asquith et al. (2017). Parameters are estimated on the log scale:
##   mu_log    = mean of log(x)       [location]
##   sigma_log = scale of log(x)      [spread]
##   gamma     = shape of log(x)      [skewness]
##
## RI axis for LP3 is annual exceedance probability: RI = 1 / AEP.
## This differs conceptually from PoT RI used in GPD, but is numerically
## similar when the PoT rate ≈ 1 event/yr (as in Yu et al.).
##
## Mixture: same CDF product as GPD [Yu et al. Eq. 3], inverted numerically.

fit_lp3 <- function(x, dates, min_n = 5) {
  keep  <- is.finite(x) & x > 0
  x     <- x[keep]
  dates <- dates[keep]
  if (length(x) < min_n) return(NULL)
  n_years <- as.numeric(diff(range(dates))) / 365.25
  if (n_years <= 0) return(NULL)

  lx   <- log(x)
  pars <- parpe3(lmoms(lx))
  if (is.null(pars$para) || anyNA(pars$para) || any(!is.finite(pars$para))) return(NULL)

  list(
    family    = "lp3",
    mu_log    = pars$para["mu"],
    sigma_log = pars$para["sigma"],
    gamma     = pars$para["gamma"],
    n         = length(x),
    n_years   = n_years
  )
}

## LP3 CDF: F(x) = P(X <= x), evaluated via Pearson III CDF on log scale.
lp3_cdf <- function(x, fit) {
  x  <- pmax(x, 1e-6)
  lx <- log(x)
  para <- vec2par(c(fit$mu_log, fit$sigma_log, fit$gamma), type = "pe3")
  sapply(lx, function(v) cdfpe3(v, para))
}

## LP3 quantile: given AEP p, return discharge.
## p = annual exceedance probability, so non-exceedance F = 1 - p.
lp3_quantile <- function(p, fit) {
  para  <- vec2par(c(fit$mu_log, fit$sigma_log, fit$gamma), type = "pe3")
  log_q <- quape3(1 - p, para)
  pmax(0, exp(log_q))
}

## AEP RI curve for a single LP3 fit: RI = 1 / AEP
lp3_curve <- function(fit, label, color,
                      ri = c(seq(1.1, 10, 0.2), seq(11, 500, 1))) {
  ep <- 1 / ri
  q  <- pmax(0, sapply(ep, lp3_quantile, fit = fit))
  data.frame(ri = ri, ep = ep, q = q, label = label, color = color,
             stringsAsFactors = FALSE)
}

## LP3 mixture SF: S_mix(x) = 1 - prod_i F_i(x)  [Yu et al. Eq. 3]
lp3_mixture_sf <- function(x, fits)
  1 - Reduce("*", lapply(fits, function(f) lp3_cdf(x, f)))

## Invert LP3 mixture SF to get RI curve by uniroot.
lp3_mixture_curve <- function(fits, ri = c(seq(1.1, 10, 0.2), seq(11, 500, 1))) {
  min_target <- min(1 / ri)
  ub <- max(sapply(fits, function(f) lp3_quantile(min_target, f)),
            na.rm = TRUE) * 2
  ub <- max(ub, 1)
  for (i in seq_len(30)) {
    if (isTRUE(lp3_mixture_sf(ub, fits) < min_target)) break
    ub <- ub * 2
  }
  q <- sapply(ri, function(r) {
    target <- 1 / r
    lb     <- 1e-3   # LP3 defined for x > 0; GPD uses mu as lower bound
    f_lb   <- lp3_mixture_sf(lb, fits) - target
    f_ub   <- lp3_mixture_sf(ub, fits) - target
    if (is.na(f_lb) || is.na(f_ub) || sign(f_lb) == sign(f_ub)) return(NA_real_)
    tryCatch(
      uniroot(function(x) lp3_mixture_sf(x, fits) - target,
              lower = lb, upper = ub, tol = 0.1)$root,
      error = function(e) NA_real_
    )
  })
  data.frame(ri = ri, ep = 1 / ri, q = pmax(0, q),
             label = "Mixture", color = "black",
             stringsAsFactors = FALSE)
}

## ── Distribution-agnostic dispatch helpers ───────────────────────────────────
##
## These allow app.R to call a single function regardless of which
## distribution family is active.

## Fit one flood-type sample. family is "gpd" or "lp3".
fit_dist <- function(x, dates, family, mu_threshold = NULL, min_n = 5) {
  if (family == "gpd") fit_gpd(x, dates, mu_threshold, min_n)
  else                 fit_lp3(x, dates, min_n)
}

## Theoretical RI curve for a single fit.
dist_curve <- function(fit, label, color,
                       ri = c(seq(1.1, 10, 0.2), seq(11, 500, 1))) {
  if (fit$family == "gpd") gpd_curve(fit, label, color, ri)
  else                     lp3_curve(fit, label, color, ri)
}

## Mixture RI curve across a named list of fits (all same family).
dist_mixture_curve <- function(fits, ri = c(seq(1.1, 10, 0.2), seq(11, 500, 1))) {
  fam <- fits[[1]]$family
  if (fam == "gpd") mixture_curve(fits, ri)
  else              lp3_mixture_curve(fits, ri)
}

## ── Empirical RI (shared by both families) ───────────────────────────────────
##
## Weibull plotting position: RI_k = (n+1)/rank, rank=1 is largest.
## Consistent with both PoT RI (GPD) and annual RI (LP3) at rate ≈ 1/yr.
empirical_ri <- function(peaks) {
  peaks %>%
    filter(!is.na(peak_flow), peak_flow > 0) %>%
    group_by(flood_type) %>%
    arrange(desc(peak_flow)) %>%
    mutate(
      rank = seq_len(n()),
      n    = n(),
      ep   = rank / (n + 1),
      ri   = (n + 1) / rank
    ) %>%
    ungroup()
}