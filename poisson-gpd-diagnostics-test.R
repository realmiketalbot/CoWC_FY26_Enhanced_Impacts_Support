# ---------------------------------------------------------------------------
# Threshold selection by gated weighted ranking  (v3)
#
# Structure:
#   (1) GATE   - retain only candidate thresholds that do not fail the
#                dispersion index chi-square test. DI is a pass/fail
#                admissibility condition on the arrival process.
#   (2) RANK A - among survivors, rank by threshold, lowest ranked best. This
#                encodes the standard POT principle: as low as possible
#                subject to the model holding. DI supplies no natural ordering
#                of its own (a DI nearer 1 is not necessarily preferable), so
#                its contribution is admissibility plus sample size.
#   (3) RANK B - among survivors, rank by R^2 of the mean excess linear fit.
#   (4) RANK C - among survivors, rank by Euclidean distance from the sample
#                (t3, t4) to the theoretical GPA curve.
#   (5) COMBINE - weighted sum of the three normalised ranks; lowest score
#                 wins, ties broken toward the lower threshold.
#
# Note on the mean excess ranking: R^2 is only comparable across thresholds if
# the windows contain the same number of points, so this version fits a FIXED
# window of WINDOW_N points starting at each candidate rather than the
# variable-length window used in v2. Variable-length windows produce R^2
# values that cannot be ranked against one another.
#
# Sign convention: xi (Coles), xi > 0 heavy/unbounded. lmomco kappa = -xi.
# ---------------------------------------------------------------------------

library(tidyverse)
library(lmomco)

# --- Parameters ------------------------------------------------------------

SITE            <- "09081600"
TYPE_SELECT     <- NULL   # NULL = all peaks; or "rainfall", "snowmelt", ...

N_CANDIDATES    <- 40
CAND_MIN_Q      <- 0.00
CAND_MAX_Q      <- 0.95

MIN_EXCEED_ME   <- 10     # points with fewer exceedances excluded from windows
WINDOW_N        <- 10      # fixed mean-excess window length, in candidates
MIN_EXCEED_LMOM <- 10

ALPHA           <- 0.05   # DI chi-square significance level

W_DI            <- 1/3    # weight on the DI rank (by threshold, lowest = best)
W_ME            <- 1/3    # weight on the mean excess rank
W_LMOM          <- 1/3    # weight on the L-moment rank

# --- Load and subset -------------------------------------------------------

peaks_pds_all <- readRDS("_data/pot_peaks_classified_frac_thresh_3class.rds")

site_peaks <- peaks_pds_all |>
  filter(site_no == SITE) |>
  filter(flood_type_raw != "unclassified")

if (!is.null(TYPE_SELECT)) {
  site_peaks <- site_peaks |> filter(flood_type == TYPE_SELECT)
}

site_peaks <- site_peaks |>
  select(site_no, peak_date, water_year, peak_flow, flood_type) |>
  arrange(peak_date)

stopifnot(nrow(site_peaks) > 0)

wy_all  <- seq(min(site_peaks$water_year), max(site_peaks$water_year))
M_years <- length(wy_all)

cat("Site:", SITE, "| type:", ifelse(is.null(TYPE_SELECT), "all peaks", TYPE_SELECT), "\n")
cat("Peaks:", nrow(site_peaks), "| water years:", M_years,
    "(", min(wy_all), "-", max(wy_all), ")\n\n")

# --- Candidate thresholds --------------------------------------------------

candidates <- quantile(
  site_peaks$peak_flow,
  probs = seq(CAND_MIN_Q, CAND_MAX_Q, length.out = N_CANDIDATES),
  names = FALSE, type = 7
) |>
  unique() |>
  sort()

n_cand <- length(candidates)
cat("Candidates:", n_cand, "| range",
    round(min(candidates), 1), "-", round(max(candidates), 1), "cfs\n\n")

# ===========================================================================
# STEP 1 - GATE: dispersion index hypothesis test
# ===========================================================================

di_lower <- qchisq(ALPHA / 2,     df = M_years - 1) / (M_years - 1)
di_upper <- qchisq(1 - ALPHA / 2, df = M_years - 1) / (M_years - 1)

compute_di <- function(u) {

  counts <- site_peaks |>
    filter(peak_flow > u) |>
    count(water_year, name = "n")

  counts_full <- tibble(water_year = wy_all) |>
    left_join(counts, by = "water_year") |>
    mutate(n = replace_na(n, 0))

  lambda_hat <- mean(counts_full$n)
  di <- if (lambda_hat > 0) var(counts_full$n) / lambda_hat else NA_real_

  tibble(
    threshold  = u,
    n_exceed   = sum(counts_full$n),
    lambda_hat = lambda_hat,
    di         = di,
    di_pass    = !is.na(di) & di > di_lower & di < di_upper,
    di_verdict = case_when(
      is.na(di)     ~ "undefined",
      di < di_lower ~ "underdispersed",
      di > di_upper ~ "overdispersed",
      TRUE          ~ "poisson"
    )
  )
}

di_tbl <- map_dfr(candidates, compute_di)

cat("STEP 1 - DI gate\n")
cat("  non-rejection interval: [", round(di_lower, 3), ",", round(di_upper, 3), "]\n")
cat("  verdicts:",
    paste(names(table(di_tbl$di_verdict)), table(di_tbl$di_verdict),
          sep = "=", collapse = ", "), "\n")
cat("  admitted:", sum(di_tbl$di_pass), "of", n_cand, "candidates\n\n")

# ===========================================================================
# STEP 2 - mean excess R^2, fixed-length window
# ===========================================================================

mean_excess <- function(u) {
  ex <- site_peaks$peak_flow[site_peaks$peak_flow > u] - u
  if (length(ex) == 0) return(tibble(threshold = u, n_exceed = 0L, me = NA_real_))
  tibble(threshold = u, n_exceed = length(ex), me = mean(ex))
}

me_points <- map_dfr(candidates, mean_excess)

me_usable <- me_points |>
  filter(n_exceed >= MIN_EXCEED_ME, !is.na(me)) |>
  arrange(threshold)

fixed_window_fit <- function(i) {

  blank <- tibble(threshold = me_usable$threshold[i],
                  me_win_end = NA_real_, me_r2 = NA_real_,
                  me_slope = NA_real_, me_xi = NA_real_)

  j <- i + WINDOW_N - 1
  if (j > nrow(me_usable)) return(blank)

  pts <- me_usable[i:j, ]
  fit <- lm(me ~ threshold, data = pts)
  r2  <- summary(fit)$r.squared
  m   <- unname(coef(fit)[2])

  tibble(
    threshold  = me_usable$threshold[i],
    me_win_end = me_usable$threshold[j],
    me_r2      = r2,
    me_slope   = m,
    me_xi      = m / (1 + m)
  )
}

me_tbl <- map_dfr(seq_len(nrow(me_usable)), fixed_window_fit)

cat("STEP 2 - mean excess\n")
cat("  usable candidates (>=", MIN_EXCEED_ME, "exceedances):", nrow(me_usable), "\n")
cat("  fixed windows of", WINDOW_N, "points fitted:",
    sum(!is.na(me_tbl$me_r2)), "\n")
cat("  R^2 range:", round(min(me_tbl$me_r2, na.rm = TRUE), 3), "-",
    round(max(me_tbl$me_r2, na.rm = TRUE), 3), "\n\n")

# ===========================================================================
# STEP 3 - L-moment distance to the GPA curve
# ===========================================================================

gpa_t4 <- function(t3) t3 * (1 + 5 * t3) / (5 + t3)

curve_grid <- tibble(t3 = seq(-0.99, 0.99, by = 0.0005)) |>
  mutate(t4 = gpa_t4(t3))

distance_to_gpa <- function(t3_hat, t4_hat) {
  if (is.na(t3_hat) || is.na(t4_hat)) return(NA_real_)
  min(sqrt((curve_grid$t3 - t3_hat)^2 + (curve_grid$t4 - t4_hat)^2))
}

compute_lmom <- function(u) {

  ex <- site_peaks$peak_flow[site_peaks$peak_flow > u] - u

  blank <- tibble(threshold = u, lmom_n = length(ex),
                  t3_hat = NA_real_, t4_hat = NA_real_,
                  lmom_dist = NA_real_, lmom_xi = NA_real_)

  if (length(ex) < MIN_EXCEED_LMOM) return(blank)

  lm_est <- lmoms(ex, nmom = 4)
  if (!are.lmom.valid(lm_est)) return(blank)

  gpa_fit <- try(pargpa(lm_est), silent = TRUE)
  xi_hat  <- if (inherits(gpa_fit, "try-error") || is.null(gpa_fit$para)) {
    NA_real_
  } else {
    -unname(gpa_fit$para["kappa"])
  }

  tibble(
    threshold = u,
    lmom_n    = length(ex),
    t3_hat    = lm_est$ratios[3],
    t4_hat    = lm_est$ratios[4],
    lmom_dist = distance_to_gpa(lm_est$ratios[3], lm_est$ratios[4]),
    lmom_xi   = xi_hat
  )
}

lmom_tbl <- map_dfr(candidates, compute_lmom)

cat("STEP 3 - L-moment ratios\n")
cat("  computed for", sum(!is.na(lmom_tbl$lmom_dist)), "of", n_cand, "candidates\n")
cat("  distance range:", round(min(lmom_tbl$lmom_dist, na.rm = TRUE), 4), "-",
    round(max(lmom_tbl$lmom_dist, na.rm = TRUE), 4), "\n\n")

# ===========================================================================
# STEP 4 - intersect, rank, and weight
# ===========================================================================
# Ranks are computed WITHIN the admitted set, so a threshold's rank reflects
# its standing among viable candidates rather than among all candidates.
# Normalised to [0, 1] so the weights mean what they appear to mean when the
# two criteria have different numbers of valid entries.

all_tbl <- di_tbl |>
  left_join(me_tbl   |> select(-any_of("n_exceed")), by = "threshold") |>
  left_join(lmom_tbl |> select(-any_of("n_exceed")), by = "threshold")

admitted <- all_tbl |>
  filter(di_pass, !is.na(me_r2), !is.na(lmom_dist))

used_fallback <- FALSE

if (nrow(admitted) == 0) {
  used_fallback <- TRUE
  cat("!! No candidate passes the DI gate with both metrics available.\n")
  cat("!! Falling back to ranking over all candidates with both metrics.\n\n")
  admitted <- all_tbl |> filter(!is.na(me_r2), !is.na(lmom_dist))
}

stopifnot(nrow(admitted) > 0)

ranked <- admitted |>
  mutate(
    rank_di    = rank(threshold,  ties.method = "min"),   # lowest threshold is best
    rank_me    = rank(-me_r2,     ties.method = "min"),   # high R^2 is good
    rank_lmom  = rank(lmom_dist,  ties.method = "min"),   # low distance is good
    n_ranked   = n(),
    prank_di   = if (n() > 1) (rank_di   - 1) / (n() - 1) else 0,
    prank_me   = if (n() > 1) (rank_me   - 1) / (n() - 1) else 0,
    prank_lmom = if (n() > 1) (rank_lmom - 1) / (n() - 1) else 0,
    score      = W_DI * prank_di + W_ME * prank_me + W_LMOM * prank_lmom,
    xi_diff    = me_xi - lmom_xi
  ) |>
  arrange(score, threshold) |>
  mutate(rank_overall = row_number())

cat("STEP 4 - weighted ranking\n")
cat("  ranked set:", nrow(ranked), "candidates",
    if (used_fallback) "(FALLBACK: DI gate not applied)" else "(DI gate applied)", "\n")
cat("  weights: DI(threshold) =", round(W_DI, 3),
    "| mean excess =", round(W_ME, 3),
    "| L-moments =", round(W_LMOM, 3), "\n\n")

print(
  ranked |>
    select(rank_overall, threshold, n_exceed, lambda_hat,
           di, di_verdict, rank_di,
           me_r2, rank_me, me_xi,
           lmom_dist, rank_lmom, lmom_xi,
           xi_diff, score),
  n = Inf
)

# --- Selected --------------------------------------------------------------

selected <- ranked |> slice_min(order_by = score, n = 1, with_ties = FALSE)

cat("\n=== Selected threshold ===\n")
cat("  threshold:   ", round(selected$threshold, 1), "cfs\n")
cat("  exceedances: ", selected$n_exceed, "| lambda =", round(selected$lambda_hat, 3), "\n")
cat("  DI:          ", round(selected$di, 3), "(", selected$di_verdict, ")",
    "| threshold rank", selected$rank_di, "of", selected$n_ranked, "\n")
cat("  mean excess: R^2 =", round(selected$me_r2, 3),
    "| rank", selected$rank_me, "of", selected$n_ranked, "\n")
cat("  L-moments:   dist =", round(selected$lmom_dist, 4),
    "| rank", selected$rank_lmom, "of", selected$n_ranked, "\n")
cat("  weighted score:", round(selected$score, 4), "\n")
cat("  implied xi:  mean excess =", round(selected$me_xi, 3),
    "| L-moments =", round(selected$lmom_xi, 3),
    "| diff =", round(selected$xi_diff, 3), "\n")
if (used_fallback) cat("  NOTE: selected without a DI gate.\n")

# --- Do the criteria agree? ------------------------------------------------
# Pairwise rank correlation across the admitted set. Strong negative
# correlation means two criteria want opposite ends of the range and the
# weighted score is a compromise neither would endorse. Note that rank_di is
# by construction monotone in threshold, so its correlation with the others
# says which end of the range each of them prefers.

if (nrow(ranked) > 2) {
  rho_me_lmom <- cor(ranked$rank_me,  ranked$rank_lmom, method = "spearman")
  rho_di_me   <- cor(ranked$rank_di,  ranked$rank_me,   method = "spearman")
  rho_di_lmom <- cor(ranked$rank_di,  ranked$rank_lmom, method = "spearman")

  cat("\n=== Pairwise rank correlations (Spearman) ===\n")
  cat("  mean excess vs L-moments:", round(rho_me_lmom, 3),
      if (rho_me_lmom < 0) " (disagree)\n" else " (agree)\n")
  cat("  threshold   vs mean excess:", round(rho_di_me, 3),
      if (rho_di_me < 0) " (mean excess prefers higher thresholds)\n"
      else " (mean excess prefers lower thresholds)\n")
  cat("  threshold   vs L-moments:", round(rho_di_lmom, 3),
      if (rho_di_lmom < 0) " (L-moments prefer higher thresholds)\n"
      else " (L-moments prefer lower thresholds)\n")
}

# --- Plots -----------------------------------------------------------------

p_gate <- ggplot(all_tbl, aes(x = threshold, y = di)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = di_lower, ymax = di_upper,
           fill = "grey80", alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(colour = "grey50") +
  geom_point(aes(colour = di_pass), size = 2) +
  geom_vline(xintercept = selected$threshold, colour = "firebrick") +
  scale_x_log10() +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey70"),
                      name = "admitted") +
  labs(title = paste0("Step 1: DI gate — ", SITE),
       x = "Threshold (cfs, log scale)", y = "DI") +
  theme_minimal()

p_score <- ranked |>
  select(threshold, prank_di, prank_me, prank_lmom, score) |>
  pivot_longer(-threshold, names_to = "component", values_to = "value") |>
  ggplot(aes(x = threshold, y = value, colour = component)) +
  geom_line() +
  geom_point(size = 1.6) +
  geom_vline(xintercept = selected$threshold, colour = "firebrick") +
  scale_x_log10() +
  scale_colour_manual(
    values = c(prank_di = "forestgreen", prank_me = "steelblue",
               prank_lmom = "darkorange", score = "black"),
    labels = c(prank_di = "threshold rank (DI-admitted)",
               prank_lmom = "L-moment rank", prank_me = "mean excess rank",
               score = "weighted score"),
    name = NULL) +
  labs(title = paste0("Step 4: normalised ranks and weighted score — ", SITE),
       subtitle = "Lower is better; ranked within the admitted set",
       x = "Threshold (cfs, log scale)", y = "normalised rank / score") +
  theme_minimal()

win_pts <- me_usable |>
  filter(threshold >= selected$threshold,
         threshold <= replace_na(selected$me_win_end, max(me_usable$threshold)))

p_me <- ggplot(me_usable, aes(x = threshold, y = me)) +
  geom_point(colour = "grey65") +
  geom_point(data = win_pts, colour = "black", size = 2) +
  geom_smooth(data = win_pts, method = "lm", se = FALSE, colour = "firebrick") +
  geom_vline(xintercept = selected$threshold, colour = "firebrick") +
  labs(title = paste0("Mean excess, selected window — ", SITE),
       subtitle = paste0("Fixed window of ", WINDOW_N, " points, R2 = ",
                         round(selected$me_r2, 3)),
       x = "Threshold (cfs)", y = "Mean excess (cfs)") +
  theme_minimal()

p_lmom <- ggplot(lmom_tbl |> filter(!is.na(t3_hat)),
                 aes(x = t3_hat, y = t4_hat)) +
  geom_line(data = curve_grid |> filter(t3 >= -0.2, t3 <= 0.8),
            aes(x = t3, y = t4), inherit.aes = FALSE, linetype = "dashed") +
  geom_path(colour = "grey65") +
  geom_point(aes(colour = threshold), size = 2.2) +
  geom_point(data = lmom_tbl |> filter(threshold == selected$threshold),
             shape = 21, size = 5, stroke = 1.2, colour = "firebrick") +
  scale_colour_viridis_c(trans = "log10", name = "Threshold") +
  labs(title = paste0("L-moment ratio diagram — ", SITE),
       subtitle = "Circled: selected threshold. Dashed: theoretical GPA relationship.",
       x = expression(hat(tau)[3]), y = expression(hat(tau)[4])) +
  theme_minimal()

print(p_gate)
print(p_me)
print(p_lmom)
print(p_score)

# --- Notes -----------------------------------------------------------------
# - The DI gate admits BOTH the low-threshold region (several events per year,
#   genuine year-to-year variability) and the high-threshold region (events
#   rare enough for the Poisson limit). These pass for different reasons. If
#   the ranking keeps landing in the low region, consider gating on the upper
#   passing region only, or adding a minimum-threshold constraint.
# - Ranks are computed within the admitted set. Comparing scores across sites
#   is therefore not meaningful; compare the selected thresholds and their
#   underlying R^2 / distance values instead.
# - WINDOW_N controls what "linear" means for the mean excess ranking. Short
#   windows fit almost anything; long windows may span curvature. Vary it in
#   the sensitivity analysis alongside W_ME / W_LMOM.
# - Equal weights are a placeholder. If one criterion proves more reliable
#   across the network, reweight and report both.