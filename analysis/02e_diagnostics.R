# Phase 2.d: cheap pre-tests on the Phase 2.2 level-based MS-AR(1) regime model,
# reusing Phase 2.2/2.3 outputs with NO refitting.
#
# Two questions, before any Phase-3 forecasting/simulation work:
#
#   1. Is filt_prob_high (the real-time expanding-window filtered P(high-fear))
#      a trustworthy probability?
#        1a. Rigorous, observable-based: PIT + interval coverage of the model's
#            one-step-ahead predictive density for log(VIX) (a regime mixture).
#            Needs no regime label.
#        1b. Internal consistency: does filt_prob_high track the SMOOTHED
#            (full-sample, retrospective) high-fear classification? The regime
#            is a LATENT state with no observed ground truth, so smoothed is a
#            proxy, not truth.
#
#   2. Does the real-time regime carry forward information about VIX at all?
#      Purely empirical (model-free) forward log(VIX) behavior conditional on
#      the real-time regime state. This is the benchmark a Phase-3 simulator
#      would need to reproduce, and a pre-test of whether the forecasting
#      application has legs.
#
# CAVEAT (documented, not hidden): Diagnostic 1a needs regime-conditional means/
# AR/sigma and the transition matrix. The Phase 2.2 expanding-window OOS loop
# *refits* these every ~21 trading days, so the parameters used to produce
# filt_prob_high[t] vary slightly over OOS history. Only the FINAL full-sample
# fit1 parameters were persisted (02c_regime_summary.csv). Diagnostic 1a applies
# those fixed final parameters throughout OOS as the only no-refit option.
# filt_prob_high itself remains the honest real-time series; only the
# regime-conditional densities used to *score* it carry this approximation.

library(tidyverse)

source("R/hysteresis.R") # schmitt_trigger()
source("R/oos_forecast.R") # versioned_path()

# ---- Parameters --------------------------------------------------------------

OUTPUT_DIR <- "outputs"
HYST_UPPER <- 0.6 # 0.6/0.4 hysteresis def, reused from Phase 2.3
HYST_LOWER <- 0.4
POST_TRANSITION_WINDOW <- 10 # trading days after a High->Low flip classified "post_transition_low"
HORIZONS <- c(1, 5, 10, 21) # forward horizons (trading days) for Diagnostic 2
PIT_BINS <- 20
RELIABILITY_BINS <- 10 # decile bins of filt_prob_high

# ---- Load Phase 2.2 outputs (no refit) -----------------------------------------

oos <- readr::read_csv("outputs/02c_oos_filtered.csv", show_col_types = FALSE) |>
  dplyr::arrange(date) |>
  dplyr::mutate(log_vix = log(vix_close))

probs1 <- readr::read_csv("outputs/02c_full_sample_probs.csv", show_col_types = FALSE)

regime_summary <- readr::read_csv("outputs/02c_regime_summary.csv", show_col_types = FALSE) |>
  dplyr::filter(model == "fit1_common_AR")

cat(sprintf(
  "Loaded real-time filtered series (oos): n = %d, %s to %s\n",
  nrow(oos), min(oos$date), max(oos$date)
))
cat(sprintf(
  "Loaded full-sample probs (probs1): n = %d, %s to %s\n",
  nrow(probs1), min(probs1$date), max(probs1$date)
))

# ---- Regime parameters + transition matrix (reconstructed, no refit) -----------

low_row <- regime_summary |> dplyr::filter(regime == "low_fear")
high_row <- regime_summary |> dplyr::filter(regime == "high_fear")

ar_coef <- low_row$ar_coef # common across regimes by construction (fit1)
mean_low <- low_row$mean_log_vix
mean_high <- high_row$mean_log_vix
sigma_low <- low_row$sigma
sigma_high <- high_row$sigma
duration_low <- low_row$duration_days
duration_high <- high_row$duration_days

# duration_days = 1 / (1 - P_ii); P is column-stochastic (R/level_regime.R::transition_stats)
P_low_low <- 1 - 1 / duration_low
P_high_high <- 1 - 1 / duration_high
P_low_high <- 1 - P_high_high # column "high" sums to 1
P_high_low <- 1 - P_low_low # column "low" sums to 1
P <- matrix(c(P_low_low, P_high_low, P_low_high, P_high_high), nrow = 2) # cols = (low, high)

cat("\n--- Reconstructed transition matrix (cols = (low, high), column-stochastic) ---\n")
print(P)
cat(sprintf("Column sums (should be ~1): %s\n", paste(signif(colSums(P), 6), collapse = ", ")))
cat(sprintf(
  "Implied durations: low = %.2f (input %.2f), high = %.2f (input %.2f)\n",
  1 / (1 - P[1, 1]), duration_low, 1 / (1 - P[2, 2]), duration_high
))

# =================================================================================
# Diagnostic 1a: PIT + interval coverage of the one-step predictive density
# =================================================================================

n <- nrow(oos)

pi_filt_low <- 1 - oos$filt_prob_high[1:(n - 1)]
pi_filt_high <- oos$filt_prob_high[1:(n - 1)]

pi_pred_low <- P[1, 1] * pi_filt_low + P[1, 2] * pi_filt_high
pi_pred_high <- P[2, 1] * pi_filt_low + P[2, 2] * pi_filt_high

cond_mean_low <- mean_low * (1 - ar_coef) + ar_coef * oos$log_vix[1:(n - 1)]
cond_mean_high <- mean_high * (1 - ar_coef) + ar_coef * oos$log_vix[1:(n - 1)]

y_realized <- oos$log_vix[2:n]

pit <- pi_pred_low * pnorm(y_realized, cond_mean_low, sigma_low) +
  pi_pred_high * pnorm(y_realized, cond_mean_high, sigma_high)

pit_df <- tibble::tibble(date = oos$date[2:n], pit = pit)

cat(sprintf(
  "\n--- Diagnostic 1a: PIT summary (n = %d) ---\nMean = %.4f (uniform ref = 0.5), SD = %.4f (uniform ref = %.4f)\n",
  nrow(pit_df), mean(pit), sd(pit), sqrt(1 / 12)
))

p_pit <- ggplot(pit_df, aes(x = pit)) +
  geom_histogram(bins = PIT_BINS, fill = "#0072B2", color = "white", boundary = 0) +
  geom_hline(yintercept = nrow(pit_df) / PIT_BINS, linetype = "dashed", color = "#D55E00") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "PIT of Realized log(VIX) under One-Step Regime-Mixture Predictive Density",
    subtitle = sprintf(
      "n = %d, OOS %s to %s. Approximately uniform on [0,1] if calibrated; dashed line = uniform reference count",
      nrow(pit_df), min(pit_df$date), max(pit_df$date)
    ),
    x = "Probability Integral Transform (PIT)", y = "Count"
  ) +
  theme_minimal()

pit_plot_path <- versioned_path(OUTPUT_DIR, "02e_pit_histogram", "png")
ggsave(pit_plot_path, p_pit, width = 10, height = 6, dpi = 150)
cat(sprintf("Saved: %s\n", pit_plot_path))

# Coverage: F is monotonic, so the central alpha-interval contains y iff
# pit falls in [(1-alpha)/2, (1+alpha)/2] -- no uniroot needed.
coverage <- tibble::tibble(nominal = c(0.5, 0.8, 0.95)) |>
  dplyr::mutate(
    lower = (1 - nominal) / 2,
    upper = (1 + nominal) / 2,
    empirical_coverage = purrr::map2_dbl(lower, upper, ~ mean(pit >= .x & pit <= .y))
  )

cat("\n--- Diagnostic 1a: interval coverage (nominal vs empirical) ---\n")
print(as.data.frame(coverage))

coverage_path <- versioned_path(OUTPUT_DIR, "02e_pit_coverage", "csv")
readr::write_csv(coverage, coverage_path)
cat(sprintf("Saved: %s\n", coverage_path))

# =================================================================================
# Diagnostic 1b: reliability curve (filt_prob_high deciles vs smoothed-high freq)
# =================================================================================

rel <- oos |>
  dplyr::select(date, filt_prob_high) |>
  dplyr::inner_join(probs1 |> dplyr::select(date, smooth_prob_high), by = "date") |>
  dplyr::mutate(
    smoothed_high = smooth_prob_high > 0.5,
    decile = dplyr::ntile(filt_prob_high, RELIABILITY_BINS)
  )

reliability <- rel |>
  dplyr::group_by(decile) |>
  dplyr::summarize(
    mean_filt_prob_high = mean(filt_prob_high),
    freq_smoothed_high = mean(smoothed_high),
    n = dplyr::n(),
    .groups = "drop"
  )

cat(sprintf("\n--- Diagnostic 1b: reliability table (n = %d matched dates) ---\n", nrow(rel)))
print(as.data.frame(reliability))

reliability_table_path <- versioned_path(OUTPUT_DIR, "02e_reliability_table", "csv")
readr::write_csv(reliability, reliability_table_path)
cat(sprintf("Saved: %s\n", reliability_table_path))

p_rel <- ggplot(reliability, aes(x = mean_filt_prob_high, y = freq_smoothed_high)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(color = "#0072B2") +
  geom_point(aes(size = n), color = "#0072B2") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Reliability Curve: Real-Time filt_prob_high vs Smoothed High-Fear Frequency",
    subtitle = "Decile bins of real-time filt_prob_high; y = fraction of days in bin with SMOOTHED P(high-fear) > 0.5.\nSmoothed is a retrospective proxy, NOT ground truth -- the regime is latent. Dashed = 45-degree line.",
    x = "Mean real-time P(high-fear), filtered (decile bin)",
    y = "Frequency smoothed P(high-fear) > 0.5",
    size = "n days"
  ) +
  theme_minimal()

reliability_plot_path <- versioned_path(OUTPUT_DIR, "02e_reliability_curve", "png")
ggsave(reliability_plot_path, p_rel, width = 10, height = 6, dpi = 150)
cat(sprintf("Saved: %s\n", reliability_plot_path))

# =================================================================================
# Diagnostic 2: regime-conditional forward log(VIX) behavior (empirical, model-free)
# =================================================================================

is_high <- schmitt_trigger(oos$filt_prob_high, HYST_UPPER, HYST_LOWER)

group <- character(n)
last_drop_idx <- -Inf
for (t in seq_len(n)) {
  if (t > 1 && is_high[t - 1] && !is_high[t]) {
    last_drop_idx <- t
  }
  if (is_high[t]) {
    group[t] <- "stable_high"
  } else if ((t - last_drop_idx) < POST_TRANSITION_WINDOW) {
    group[t] <- "post_transition_low"
  } else {
    group[t] <- "stable_low"
  }
}

cat("\n--- Diagnostic 2: real-time regime-state classification (0.6/0.4 hysteresis) ---\n")
print(table(group))

forward_rows <- vector("list", length(HORIZONS))
for (i in seq_along(HORIZONS)) {
  k <- HORIZONS[i]
  valid_t <- seq_len(n - k)

  d_k <- oos$log_vix[valid_t + k] - oos$log_vix[valid_t]
  # RMS of forward daily log-changes; defined for k = 1 (= |d|) as well.
  rv_k <- vapply(valid_t, function(t) sqrt(mean(diff(oos$log_vix[t:(t + k)])^2)), numeric(1))

  forward_rows[[i]] <- tibble::tibble(
    group = group[valid_t],
    horizon = k,
    d = d_k,
    rv = rv_k
  )
}
forward_data <- dplyr::bind_rows(forward_rows)

forward_summary <- forward_data |>
  dplyr::group_by(group, horizon) |>
  dplyr::summarize(
    n = dplyr::n(),
    mean = mean(d),
    median = median(d),
    q10 = quantile(d, 0.10),
    q25 = quantile(d, 0.25),
    q75 = quantile(d, 0.75),
    q90 = quantile(d, 0.90),
    p_rise = mean(d > 0),
    fwd_realized_vol_mean = mean(rv),
    .groups = "drop"
  ) |>
  dplyr::mutate(group = factor(group, levels = c("stable_low", "post_transition_low", "stable_high"))) |>
  dplyr::arrange(horizon, group)

cat("\n--- Diagnostic 2: forward log(VIX) change by regime state and horizon ---\n")
print(as.data.frame(forward_summary))

forward_summary_path <- versioned_path(OUTPUT_DIR, "02e_forward_vix_by_regime", "csv")
readr::write_csv(forward_summary, forward_summary_path)
cat(sprintf("Saved: %s\n", forward_summary_path))

forward_plot_data <- forward_data |>
  dplyr::mutate(
    group = factor(group, levels = c("stable_low", "post_transition_low", "stable_high")),
    horizon_label = factor(
      horizon,
      levels = HORIZONS,
      labels = paste0("k = ", HORIZONS, " days")
    )
  )

p_fwd <- ggplot(forward_plot_data, aes(x = group, y = d)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(fill = "#0072B2", alpha = 0.4, outlier.size = 0.5) +
  facet_wrap(~horizon_label, ncol = 2) +
  labs(
    title = "Forward log(VIX) Change by Real-Time Regime State (Empirical, Model-Free)",
    subtitle = "d = log(VIX[t+k]) - log(VIX[t]), grouped by real-time filt_prob_high state at t (0.6/0.4 hysteresis,\npost-transition window = 10 days). Not a model-based simulation.",
    x = "Real-time regime state at t", y = "Forward change in log(VIX), d"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

forward_plot_path <- versioned_path(OUTPUT_DIR, "02e_forward_vix_by_regime", "png")
ggsave(forward_plot_path, p_fwd, width = 12, height = 8, dpi = 150)
cat(sprintf("Saved: %s\n", forward_plot_path))

# =================================================================================
# Verdict
# =================================================================================

cat("\n========================== VERDICT ==========================\n")

# (i) PIT / coverage calibration
max_cov_dev <- max(abs(coverage$empirical_coverage - coverage$nominal))
pit_calibrated <- max_cov_dev < 0.07

cat(sprintf(
  "(i) One-step predictive density: PIT mean = %.3f (ref 0.5), SD = %.3f (ref %.3f).\n",
  mean(pit), sd(pit), sqrt(1 / 12)
))
cat(sprintf(
  "    Coverage (50/80/95%%): %.1f%% / %.1f%% / %.1f%% (nominal 50/80/95%%, max dev %.1f pts) -> %s\n",
  100 * coverage$empirical_coverage[1], 100 * coverage$empirical_coverage[2], 100 * coverage$empirical_coverage[3],
  100 * max_cov_dev, if (pit_calibrated) "approximately calibrated" else "MISCALIBRATED"
))

# (ii) Reliability curve vs 45-degree line
reliability_dev <- mean(abs(reliability$freq_smoothed_high - reliability$mean_filt_prob_high))
reliability_tracks <- reliability_dev < 0.10

cat(sprintf(
  "(ii) Real-time filt_prob_high vs smoothed-high frequency: mean |deviation from 45-degree| = %.3f -> %s\n",
  reliability_dev, if (reliability_tracks) "tracks the smoothed (retrospective) state reasonably well" else "does NOT track the smoothed state well"
)
)
cat("     (Smoothed is a retrospective proxy, not ground truth -- the regime is latent.)\n")

# (iii) Regime-conditional forward VIX information
hi21 <- forward_summary |> dplyr::filter(group == "stable_high", horizon == 21)
lo21 <- forward_summary |> dplyr::filter(group == "stable_low", horizon == 21)
pt21 <- forward_summary |> dplyr::filter(group == "post_transition_low", horizon == 21)

sep_21 <- hi21$mean - lo21$mean
reversion_visible <- hi21$mean < -0.05 && abs(sep_21) > 0.05

cat(sprintf(
  "(iii) Forward log(VIX) change at k = 21: stable_high mean = %.4f, stable_low mean = %.4f, post_transition_low mean = %.4f.\n",
  hi21$mean, lo21$mean, pt21$mean
))
cat(sprintf(
  "      Separation (stable_high - stable_low) = %.4f -> %s\n",
  sep_21, if (reversion_visible) "meaningful mean-reversion asymmetry visible" else "distributions look similar -- little/no forward information"
))

cat("\n--------------------------------------------------------------\n")
if (pit_calibrated && reliability_tracks && reversion_visible) {
  cat("OVERALL: filt_prob_high is reasonably calibrated, tracks the retrospective\n")
  cat("regime estimate, and the regime carries real forward VIX information\n")
  cat("(mean-reversion asymmetry across states). The forecasting application is\n")
  cat("worth pursuing (Phase 3).\n")
} else {
  cat("OVERALL: at least one diagnostic failed -- see (i)/(ii)/(iii) above.\n")
  if (!pit_calibrated) cat("  - Predictive density is miscalibrated; treat filt_prob_high probabilities with caution.\n")
  if (!reliability_tracks) cat("  - Real-time filt_prob_high diverges from the smoothed retrospective estimate.\n")
  if (!reversion_visible) cat("  - Regime-conditional forward log(VIX) distributions do not differ meaningfully:\n    the regime carries little/no forward information, and a Phase-3 forecasting\n    simulator built on this signal is NOT supported by this evidence.\n")
}
cat("================================================================\n")
