#!/usr/bin/env Rscript

# ============================================================================
# Canonical GMM + SVM boundary derivation (v4)
# ============================================================================
# Fixes from v3:
#   - B1/B2 lines now clipped to full [0,1] x [sw_min, sw_max] range,
#     not to the interior subset range. The interior subset determines
#     which events are used for fitting; the line is then extended across
#     the full feature space for plotting and classification.
#   - Per-boundary SVM cost: B1 and B2 can use different cost values.
#     Default remains 0.01 for both but can be tuned independently.
#
# End-members (rain_frac <= 0.01, >= 0.99) excluded from all SVM fitting.
# Both boundaries fit on the same interior subset.
# Panels B and C show only the two relevant classes for each boundary.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(mclust)
  library(e1071)
  library(patchwork)
})

select <- dplyr::select

# ── Config ───────────────────────────────────────────────────────────────────
proj_root <- "/Users/mtalbot/Desktop/GitHub/CoWC_FY26_Enhanced_Impacts/"
data_dir  <- file.path(proj_root, "_data")
out_dir   <- file.path(proj_root, "_diagnostics", "threshold_evaluation",
                       "gmm_rain_frac_swe", "svm_canonical")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

peaks_file <- file.path(data_dir, "pot_peaks_classified_adaptive.rds")
labels_out <- file.path(data_dir, "gmm_canonical_labels.rds")

END_MEMBER_LO <- 0.01
END_MEMBER_HI <- 0.99

# Per-boundary cost — tune these independently if needed
B1_COST <- 1
B2_COST <- 1

# Current 6b values — orientation panel only
CUR_B1_INT <- 0.4973;  CUR_B1_SLP <- 5.3433
CUR_B2_INT <- 0.1078;  CUR_B2_SLP <- 0.3503

# ── Load data ────────────────────────────────────────────────────────────────
cat("Loading data...\n")
peaks_adap <- readRDS(peaks_file)

gmm_data <- peaks_adap %>%
  dplyr::filter(data_sufficient, !is.na(rain_frac), !is.na(log10_swe)) %>%
  dplyr::select(site_no, peak_id, rain_frac, log10_swe)

cat(sprintf("  Peaks for GMM: %d\n", nrow(gmm_data)))

# ── Canonical GMM fit ────────────────────────────────────────────────────────
cat("\nFitting GMM (G=4, EEV, unscaled, seed=42)...\n")

set.seed(42L)
mc <- mclust::Mclust(
  data       = cbind(gmm_data$rain_frac, gmm_data$log10_swe),
  G          = 4,
  modelNames = "EEV"
)

gmm_data <- gmm_data %>%
  dplyr::mutate(component = as.integer(mc$classification))

comp_summary <- gmm_data %>%
  dplyr::group_by(component) %>%
  dplyr::summarise(
    n             = dplyr::n(),
    rain_frac_med = median(rain_frac),
    log10_swe_med = median(log10_swe),
    .groups       = "drop"
  ) %>%
  dplyr::arrange(rain_frac_med)

cat("\nComponent summary:\n")
print(comp_summary)

comp_snowmelt  <- comp_summary$component[1]
comp_ros       <- comp_summary$component[2]
comp_rain_low  <- comp_summary$component[3]
comp_rain_high <- comp_summary$component[4]

gmm_data <- gmm_data %>%
  dplyr::mutate(
    gmm_class = dplyr::case_when(
      component %in% c(comp_rain_low, comp_rain_high) ~ "rainfall",
      component == comp_ros                           ~ "ROS",
      component == comp_snowmelt                      ~ "snowmelt"
    ),
    gmm_class = factor(gmm_class, levels = c("snowmelt", "ROS", "rainfall"))
  )

saveRDS(gmm_data, labels_out)
cat(sprintf("\nSaved canonical GMM labels → %s\n", labels_out))
cat("\nFull dataset class counts:\n")
print(gmm_data %>% dplyr::count(gmm_class))

# ── Interior subset ──────────────────────────────────────────────────────────
interior <- gmm_data %>%
  dplyr::filter(rain_frac > END_MEMBER_LO,
                rain_frac < END_MEMBER_HI)

cat(sprintf(
  "\nInterior subset (%.2f < rain_frac < %.2f): %d of %d peaks\n",
  END_MEMBER_LO, END_MEMBER_HI, nrow(interior), nrow(gmm_data)))
cat("Interior class counts:\n")
print(interior %>% dplyr::count(gmm_class))

# ── SVM fitting ──────────────────────────────────────────────────────────────
# Rather than manually extracting w and b from the SVM internals (which has
# sign ambiguity depending on e1071 internal class encoding), we use the
# model's own decision.values to score points and derive the boundary.
#
# The decision function for a linear SVM is:
#   f(x) = w^T x - b
# with the boundary at f(x) = 0. We recover w and b by fitting a linear
# regression of f(x) ~ rain_frac + log10_swe on a grid, which is equivalent
# to reading off the coefficients from the hyperplane directly.
#
# Specifically: e1071 stores the bias as model$rho with the convention that
# the decision function is w^T x - rho. The sign of w relative to the class
# labels depends on the internal encoding. We correct for this by checking
# the sign of the decision value at the class_a centroid and flipping if
# needed so that class_a is always on the positive side.

fit_svm <- function(df, class_a, class_b, cost) {
  d <- df %>%
    dplyr::filter(gmm_class %in% c(class_a, class_b)) %>%
    dplyr::mutate(cls = factor(gmm_class, levels = c(class_a, class_b)))

  model <- e1071::svm(
    cls ~ rain_frac + log10_swe,
    data = d, kernel = "linear", cost = cost, scale = FALSE,
    decision.values = TRUE
  )

  # Extract decision values from the model's own predictions
  pred_obj <- predict(model, d, decision.values = TRUE)
  dv       <- as.numeric(attr(pred_obj, "decision.values"))

  # Fit a linear model of decision values on the features to recover w and b
  # in a sign-consistent way: lm(dv ~ rain_frac + log10_swe)
  # gives intercept = -b, coef[rain_frac] = w_rf, coef[log10_swe] = w_sw
  lm_fit <- lm(dv ~ rain_frac + log10_swe, data = d)
  w       <- coef(lm_fit)[c("rain_frac", "log10_swe")]
  b       <- -coef(lm_fit)["(Intercept)"]
  names(w) <- c("rain_frac", "log10_swe")

  # Ensure class_a is on the positive side
  mean_a <- colMeans(d[d$cls == class_a, c("rain_frac", "log10_swe")])
  if (sum(w * mean_a) - b < 0) {
    w <- -w
    b <- -b
  }

  pred     <- as.character(pred_obj)
  acc      <- mean(pred == as.character(d$cls))
  recall_a <- mean(pred[d$cls == class_a] == class_a)
  recall_b <- mean(pred[d$cls == class_b] == class_b)

  list(
    n_a = sum(d$cls == class_a), n_b = sum(d$cls == class_b),
    n_sv = nrow(model$SV),
    w = w, b = unname(b), cost = cost,
    accuracy = acc, recall_a = recall_a, recall_b = recall_b,
    class_a = class_a, class_b = class_b,
    model = model
  )
}

cat(sprintf("\nFitting SVMs on interior subset...\n"))

b1 <- fit_svm(interior, "snowmelt", "ROS",     cost = B1_COST)
b2 <- fit_svm(interior, "ROS",      "rainfall", cost = B2_COST)

# ── Boundary validation ───────────────────────────────────────────────────────
# Compare manual w/b evaluation against SVM's own predictions to verify that
# the extracted weights and bias correctly reproduce the fitted decision boundary.
cat("
", strrep("=", 60), "
", sep="")
cat("BOUNDARY VALIDATION
")
cat(strrep("=", 60), "
", sep="")

for (bnd in list(b1, b2)) {
  d_val <- interior %>%
    dplyr::filter(gmm_class %in% c(bnd$class_a, bnd$class_b)) %>%
    dplyr::mutate(cls = factor(gmm_class,
                               levels = c(bnd$class_a, bnd$class_b)))

  # Manual evaluation using extracted w and b
  decision    <- bnd$w["rain_frac"] * d_val$rain_frac +
                 bnd$w["log10_swe"] * d_val$log10_swe - bnd$b
  pred_manual <- ifelse(decision > 0, bnd$class_a, bnd$class_b)
  misclass_manual <- mean(pred_manual != d_val$gmm_class)

  # SVM model own predictions
  pred_svm     <- as.character(predict(bnd$model, d_val))
  misclass_svm <- mean(pred_svm != d_val$gmm_class)

  cat(sprintf("
%s vs. %s:
", bnd$class_a, bnd$class_b))
  cat(sprintf("  SVM model predictions misclassification : %.4f
", misclass_svm))
  cat(sprintf("  Manual w/b evaluation misclassification : %.4f
", misclass_manual))
  cat(sprintf("  Match: %s
",
              ifelse(abs(misclass_manual - misclass_svm) < 0.001, "YES", "NO — extraction bug")))
  cat(sprintf("  model$rho = %.6f,  b = -rho = %.6f
", bnd$model$rho, -bnd$model$rho))
  cat(sprintf("  Decision value range: [%.4f, %.4f]
", min(decision), max(decision)))
  cat("  Decision value by class:
")
  print(
    d_val %>%
      dplyr::mutate(decision = decision) %>%
      dplyr::group_by(gmm_class) %>%
      dplyr::summarise(min_dec  = min(decision),
                       mean_dec = mean(decision),
                       max_dec  = max(decision), .groups = "drop")
  )
}
cat(strrep("=", 60), "
", sep="")


for (bnd in list(b1, b2)) {
  cat(sprintf(
    "\n%s vs. %s (cost=%.3f): n=%d+%d, n_sv=%d\n",
    bnd$class_a, bnd$class_b, bnd$cost, bnd$n_a, bnd$n_b, bnd$n_sv))
  cat(sprintf("  w=(%.4f, %.4f), b=%.4f\n",
              bnd$w["rain_frac"], bnd$w["log10_swe"], bnd$b))
  cat(sprintf("  accuracy=%.4f | %s recall=%.4f | %s recall=%.4f\n",
              bnd$accuracy, bnd$class_a, bnd$recall_a,
              bnd$class_b, bnd$recall_b))
  # Report line in slope/intercept form for both axes
  cat(sprintf("  log10_swe = %.4f + %.4f * rain_frac\n",
              bnd$b / bnd$w["log10_swe"],
              -bnd$w["rain_frac"] / bnd$w["log10_swe"]))
  cat(sprintf("  rain_frac = %.4f + %.4f * log10_swe\n",
              bnd$b / bnd$w["rain_frac"],
              -bnd$w["log10_swe"] / bnd$w["rain_frac"]))
}

# ── Overall agreement (sequential tree, full dataset) ───────────────────────
classified <- gmm_data %>%
  dplyr::mutate(
    pred_class = dplyr::case_when(
      rain_frac <= END_MEMBER_LO                                          ~ "snowmelt",
      rain_frac >= END_MEMBER_HI                                          ~ "rainfall",
      b1$w["rain_frac"] * rain_frac + b1$w["log10_swe"] * log10_swe > b1$b ~ "snowmelt",
      b2$w["rain_frac"] * rain_frac + b2$w["log10_swe"] * log10_swe <= b2$b ~ "ROS",
      TRUE                                                                ~ "rainfall"
    ),
    pred_class = factor(pred_class, levels = c("snowmelt", "ROS", "rainfall"))
  )

overall_agree <- mean(classified$pred_class == classified$gmm_class)
cat(sprintf("\nOverall agreement (full dataset, sequential tree): %.4f\n",
            overall_agree))

agree_tbl <- purrr::map_dfr(c("snowmelt", "ROS", "rainfall"), function(cls) {
  tp <- sum(classified$pred_class == cls & classified$gmm_class == cls)
  fp <- sum(classified$pred_class == cls & classified$gmm_class != cls)
  fn <- sum(classified$pred_class != cls & classified$gmm_class == cls)
  tibble::tibble(
    class     = cls,
    n_true    = sum(classified$gmm_class == cls),
    precision = if (tp + fp > 0) tp / (tp + fp) else NA_real_,
    recall    = if (tp + fn > 0) tp / (tp + fn) else NA_real_
  )
})
print(agree_tbl)

# ── Figure ───────────────────────────────────────────────────────────────────
cat("\nBuilding figure...\n")

class_colors <- c("snowmelt" = "#5B8DB8",
                  "ROS"      = "#6A9E6F",
                  "rainfall" = "#E07B39")

# 95% confidence ellipses
chi2_95     <- qchisq(0.95, df = 2)
theta_seq   <- seq(0, 2 * pi, length.out = 300)
unit_circle <- rbind(cos(theta_seq), sin(theta_seq))

ellipse_df <- purrr::map_dfr(seq_len(4L), function(k) {
  mu    <- mc$parameters$mean[, k]
  sigma <- mc$parameters$variance$sigma[, , k]
  eig   <- eigen(sigma)
  axes  <- eig$vectors %*% diag(sqrt(eig$values * chi2_95))
  pts   <- t(axes %*% unit_circle)
  cls   <- dplyr::case_when(
    k %in% c(comp_rain_low, comp_rain_high) ~ "rainfall",
    k == comp_ros                           ~ "ROS",
    k == comp_snowmelt                      ~ "snowmelt"
  )
  tibble::tibble(x = pts[, 1] + mu[1], y = pts[, 2] + mu[2],
                 component = factor(k), gmm_class = cls)
})

# Line builder: w_rf * rain_frac + w_sw * log10_swe = b
# Clip to FULL data range (not interior range) so lines render across
# the whole plot, not just where interior events happen to exist.
sw_full_range <- range(gmm_data$log10_swe)

build_line <- function(w, b, rf_lo = 0, rf_hi = 1) {
  rf_seq <- seq(rf_lo, rf_hi, length.out = 400)
  tibble::tibble(
    rain_frac = rf_seq,
    log10_swe = (b - w["rain_frac"] * rf_seq) / w["log10_swe"]
  ) %>%
    dplyr::filter(log10_swe >= sw_full_range[1],
                  log10_swe <= sw_full_range[2])
}

# Current 6b lines expressed as w/b for consistent rendering
# B1: log10_swe = CUR_B1_INT + CUR_B1_SLP * rain_frac
#     => -CUR_B1_SLP * rain_frac + 1 * log10_swe = CUR_B1_INT
cur_b1_w <- c(rain_frac = -CUR_B1_SLP, log10_swe = 1)
cur_b1_b <- CUR_B1_INT
# B2: rain_frac = CUR_B2_INT + CUR_B2_SLP * log10_swe
#     => 1 * rain_frac - CUR_B2_SLP * log10_swe = CUR_B2_INT
cur_b2_w <- c(rain_frac = 1, log10_swe = -CUR_B2_SLP)
cur_b2_b <- CUR_B2_INT

base_scatter <- function(data, highlight, title_str, subtitle_str) {
  ell <- ellipse_df %>% dplyr::filter(gmm_class %in% highlight)
  ggplot(data %>% dplyr::filter(gmm_class %in% highlight),
         aes(x = rain_frac, y = log10_swe)) +
    geom_point(aes(color = gmm_class), alpha = 0.20, size = 0.5) +
    geom_path(data = ell,
              aes(x = x, y = y, group = component, color = gmm_class),
              linewidth = 0.7, inherit.aes = FALSE) +
    scale_color_manual(values = class_colors, name = "GMM class",
                       limits = highlight) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       name = "Rain fraction") +
    scale_y_continuous(name = "log10(SWE + 1)",
                       limits = sw_full_range) +
    labs(title = title_str, subtitle = subtitle_str) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

p_cur <- base_scatter(
  gmm_data, c("snowmelt", "ROS", "rainfall"),
  "Current 6b (orientation)",
  sprintf("B1: log10_swe = %.4f + %.4f*rf  |  B2: rf = %.4f + %.4f*log10_swe",
          CUR_B1_INT, CUR_B1_SLP, CUR_B2_INT, CUR_B2_SLP)
) +
  geom_line(data = build_line(cur_b1_w, cur_b1_b),
            aes(x = rain_frac, y = log10_swe), color = "black",
            linewidth = 0.7, linetype = "dashed", inherit.aes = FALSE) +
  geom_line(data = build_line(cur_b2_w, cur_b2_b),
            aes(x = rain_frac, y = log10_swe), color = "black",
            linewidth = 0.7, linetype = "dashed", inherit.aes = FALSE)

p_b1 <- base_scatter(
  gmm_data, c("snowmelt", "ROS"),
  sprintf("SVM B1 — Snowmelt / ROS  (cost=%.3f)", B1_COST),
  sprintf("w=(%.3f, %.3f), b=%.3f | accuracy=%.3f | n_sv=%d",
          b1$w["rain_frac"], b1$w["log10_swe"], b1$b, b1$accuracy, b1$n_sv)
) +
  geom_line(data = build_line(b1$w, b1$b),
            aes(x = rain_frac, y = log10_swe), color = "black",
            linewidth = 0.8, inherit.aes = FALSE)

p_b2 <- base_scatter(
  gmm_data, c("ROS", "rainfall"),
  sprintf("SVM B2 — ROS / rainfall  (cost=%.3f)", B2_COST),
  sprintf("w=(%.3f, %.3f), b=%.3f | accuracy=%.3f | n_sv=%d",
          b2$w["rain_frac"], b2$w["log10_swe"], b2$b, b2$accuracy, b2$n_sv)
) +
  geom_line(data = build_line(b2$w, b2$b),
            aes(x = rain_frac, y = log10_swe), color = "black",
            linewidth = 0.8, inherit.aes = FALSE)

fig <- (p_cur | p_b1 | p_b2) +
  plot_annotation(
    title    = "Canonical GMM + linear SVM (end-members excluded a priori)",
    subtitle = sprintf(
      "Interior: rain_frac in (%.2f, %.2f), n=%d | B1 cost=%.3f | B2 cost=%.3f | Overall agreement: %.4f",
      END_MEMBER_LO, END_MEMBER_HI, nrow(interior),
      B1_COST, B2_COST, overall_agree),
    theme = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

fig_path <- file.path(out_dir, "canonical_boundary_comparison.png")
ggsave(fig_path, plot = fig, width = 16, height = 5.5, dpi = 150)
cat(sprintf("Wrote: %s\n", fig_path))

# ── Text summary ─────────────────────────────────────────────────────────────
sink(file.path(out_dir, "canonical_boundary_summary.txt"))
cat("Canonical GMM + SVM (v4 — end-members excluded a priori)\n")
cat(format(Sys.time()), "\n\n")
cat(sprintf("GMM: G=4, EEV, seed=42, unscaled, n=%d\n", nrow(gmm_data)))
cat(sprintf("Interior: rain_frac in (%.2f, %.2f), n=%d\n\n",
            END_MEMBER_LO, END_MEMBER_HI, nrow(interior)))
cat("Component summary:\n"); print(comp_summary)
cat("\nMerged class counts (full):\n")
print(gmm_data %>% dplyr::count(gmm_class))
cat("\nInterior class counts:\n")
print(interior %>% dplyr::count(gmm_class))

for (bnd in list(b1, b2)) {
  cat(sprintf("\n%s vs. %s (cost=%.3f):\n", bnd$class_a, bnd$class_b, bnd$cost))
  cat(sprintf("  n_%s=%d, n_%s=%d, n_sv=%d\n",
              bnd$class_a, bnd$n_a, bnd$class_b, bnd$n_b, bnd$n_sv))
  cat(sprintf("  w=(%.4f, %.4f), b=%.4f\n",
              bnd$w["rain_frac"], bnd$w["log10_swe"], bnd$b))
  cat(sprintf("  accuracy=%.4f | %s recall=%.4f | %s recall=%.4f\n",
              bnd$accuracy, bnd$class_a, bnd$recall_a,
              bnd$class_b, bnd$recall_b))
  cat(sprintf("  log10_swe = %.4f + %.4f * rain_frac\n",
              bnd$b / bnd$w["log10_swe"],
              -bnd$w["rain_frac"] / bnd$w["log10_swe"]))
  cat(sprintf("  rain_frac = %.4f + %.4f * log10_swe\n",
              bnd$b / bnd$w["rain_frac"],
              -bnd$w["log10_swe"] / bnd$w["rain_frac"]))
}

cat(sprintf("\nOverall agreement (full dataset, sequential tree): %.4f\n",
            overall_agree))
cat("\nPer-class agreement:\n"); print(agree_tbl, n = Inf)
sink()

cat(sprintf("Wrote: %s\n", file.path(out_dir, "canonical_boundary_summary.txt")))
cat("\nDone.\n")