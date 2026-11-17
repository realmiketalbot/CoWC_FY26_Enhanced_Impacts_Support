#!/usr/bin/env Rscript
# Known-answer tests for the Poisson-GPD helpers in functions.R.
# Run from the repo root: Rscript _scripts/checks/test_pgpd_functions.R

source("functions.R")
ok <- function(cond, msg) { cat(if (isTRUE(cond)) "PASS" else "FAIL", "-", msg, "\n"); if (!isTRUE(cond)) quit(status = 1) }
analytic_q <- function(T, u, s, xi, lam) u + s / xi * ((-log(1 - 1/T) / lam)^(-xi) - 1)
mkfit <- function(u, s, xi, lam) list(u = u, sigma = s, xi = xi, lambda = lam,
  para = lmomco::vec2par(c(0, s, -xi), type = "gpa"))
T <- c(2, 5, 10, 25, 50, 100, 200, 500)

for (xi in c(0.2, -0.15)) {
  f <- mkfit(100, 50, xi, 1.3)
  ok(isTRUE(all.equal(pgpd_quantile(T, f), analytic_q(T, 100, 50, xi, 1.3), tolerance = 1e-10)),
     sprintf("quantile matches closed form (xi = %+.2f)", xi))
  ok(isTRUE(all.equal(pgpd_cdf(pgpd_quantile(T, f), f), 1 - 1/T, tolerance = 1e-10)),
     sprintf("CDF inverts quantile (xi = %+.2f)", xi))
}
f <- mkfit(100, 50, 0.2, 0.5)
ok(is.na(pgpd_quantile(1.5, f)) && !is.na(pgpd_quantile(10, f)),
   "quantile is NA below u (T = 1.5 with lambda = 0.5)")
ok(isTRUE(all.equal(pgpd_cdf(50, f), exp(-0.5))), "CDF floored at exp(-lambda) below u (D4)")

f1 <- mkfit(100, 50, 0.1, 0.8)
ok(isTRUE(all.equal(pgpd_mixture_quantile(T, list(f1)), pgpd_quantile(T, f1), tolerance = 1e-5)),
   "mixture of one type equals that type")
fa <- mkfit(100, 50, 0.1, 0.7); fb <- fa; fc <- fa; fc$lambda <- 1.4
ok(isTRUE(all.equal(pgpd_mixture_quantile(T, list(fa, fb)), pgpd_quantile(T, fc), tolerance = 1e-5)),
   "two identical types = one type with doubled lambda")
g1 <- mkfit(80, 40, 0.15, 1.0); g2 <- mkfit(150, 120, 0.05, 0.6)
qm <- pgpd_mixture_quantile(T, list(g1, g2))
ok(isTRUE(all.equal(pgpd_mixture_cdf(qm, list(g1, g2)), 1 - 1/T, tolerance = 1e-6)),
   "mixture root satisfies prod F_i = 1 - 1/T")
ok(all(qm >= pmax(pgpd_quantile(T, g1), pgpd_quantile(T, g2), na.rm = TRUE)),
   "mixture quantile >= each type's quantile")

set.seed(42); y <- lmomco::rlmomco(20000, lmomco::vec2par(c(0, 50, -0.2), type = "gpa"))
fit <- fit_pgpd(y + 100, u = 100, n_years = 20000 / 1.5)
ok(abs(fit$sigma - 50) < 2 && abs(fit$xi - 0.2) < 0.03 && abs(fit$lambda - 1.5) < 1e-9,
   sprintf("fit recovers sigma = 50, xi = 0.2, lambda = 1.5 (got %.1f, %.3f, %.2f)", fit$sigma, fit$xi, fit$lambda))
old <- lmomco::pargpa(lmomco::lmoms(y[1:40]))
ok(old$para["xi"] != 0, sprintf("old call without xi = 0 estimates a location (%.2f): the defect fit_pgpd avoids", old$para["xi"]))

s <- select_threshold_rate(c(10, 9, 8, 7, 6, 5), n_years = 3, floor_u = 1, min_count = 3)
ok(s$k == 3 && s$u == 7.5, "threshold rule: k = n_years when enough peaks (u between 3rd and 4th)")
s <- select_threshold_rate(c(10, 9, 8, 7), n_years = 6, floor_u = 1, min_count = 3)
ok(s$k == 4 && s$u == 4, "threshold rule: all peaks used when fewer than n_years (u midway to floor)")
ok(is.null(select_threshold_rate(c(10, 9), n_years = 6, floor_u = 1, min_count = 3)),
   "threshold rule: NULL below min_count")
s <- select_threshold_rate(c(10, 9, 8, 8, 7), n_years = 3, floor_u = 1, min_count = 3)
ok(s$k == 4 && s$u == 7.5, "threshold rule: ties at the k-th value are all included (u below the tie)")

e <- empirical_annual_rp(c(5, 50, 40, 30), u = 10, n_years = 99)
ok(identical(e$peak_flow, c(50, 40, 30)) && abs(e$ri[1] - 100.5) < 0.01,
   "empirical: exceedances only; T(rank 1, N = 99) ~ 100.5 (Langbein)")
cat("All Poisson-GPD tests passed\n")
