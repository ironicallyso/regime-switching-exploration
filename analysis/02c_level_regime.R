# Phase 2.2: 2-state Markov-switching AR(1) model on log(VIX) *levels*.
#
# Pivots the target from Phase 2 (MS-GARCH on VIX log-returns, which captures
# vol-of-vol regimes) to the level of log(VIX) itself, to test whether a
# level-based Markov-switching model produces persistent, real-time-usable
# high-fear / low-fear regimes. Not part of the SPEC.md phase sequence (same
# status as 02b_oos_forecast_compare.R for Phase 2.1) — a standalone
# diagnostic on whether level-based regimes have legs.

library(tidyverse)
library(MSwM)
library(tseries)
library(caret)

source("R/data.R")
source("R/regimes.R")
source("R/level_regime.R")
source("R/oos_forecast.R") # for versioned_path()

# ---- Parameters -----------------------------------------------------------

VIX_THRESHOLD <- 20
TRAIN_FRAC <- 0.6
REFIT_EVERY <- 21 # ~1 month, trading days
SEED <- 2024
MIN_EPISODE_DAYS <- 5 # flag high-fear runs shorter than this as implausible
OUTPUT_DIR <- "outputs"

STRESS_WINDOWS <- tibble::tribble(
  ~label, ~start, ~end,
  "2008 GFC", "2008-09-01", "2009-03-31",
  "Feb 2018 Volmageddon", "2018-02-01", "2018-02-28",
  "2020 COVID", "2020-02-15", "2020-04-30",
  "2022", "2022-01-01", "2022-12-31"
) |>
  dplyr::mutate(start = as.Date(start), end = as.Date(end))

# ---- Step 0: data & stationarity check -------------------------------------

vix_data <- fetch_vix()
log_vix <- log(vix_data$vix_close)
n <- length(log_vix)

cat("Series: log(vix_close) — LEVELS, not differenced\n")
cat(sprintf(
  "Date range: %s to %s (n = %d)\n",
  min(vix_data$date), max(vix_data$date), n
))
cat(sprintf("Mean log(VIX): %.4f (exp = %.2f)\n", mean(log_vix), exp(mean(log_vix))))

adf_result <- tseries::adf.test(log_vix)
kpss_result <- tseries::kpss.test(log_vix)

cat("\n--- ADF test (H0: unit root / non-stationary) ---\n")
print(adf_result)
cat("\n--- KPSS test (H0: stationary) ---\n")
print(kpss_result)

adf_stationary <- adf_result$p.value < 0.05
kpss_stationary <- kpss_result$p.value >= 0.05

if (!adf_stationary) {
  stop(
    "STOP: ADF test does not reject a unit root (p = ", signif(adf_result$p.value, 4),
    "). log(VIX) does not appear stationary. Per project rules, do NOT difference ",
    "the series automatically — flag this result and revisit before fitting a ",
    "level-based model."
  )
}

if (!kpss_stationary) {
  cat("\nNOTE: ADF rejects a unit root (log(VIX) is stationary), but KPSS rejects\n")
  cat("stationarity around a single constant mean (KPSS p = ", signif(kpss_result$p.value, 4), ").\n", sep = "")
  cat("This is the expected signature of a REGIME-DEPENDENT MEAN, not a unit root:\n")
  cat("KPSS's null assumes one fixed mean over the full 1990-2026 sample, which a\n")
  cat("2-state mean-switching process will violate even though each regime is itself\n")
  cat("stationary. This is consistent with (not contrary to) the premise of this\n")
  cat("phase and is not treated as a stop condition. Proceeding with the MS-AR fit\n")
  cat("without differencing.\n")
} else {
  cat("\nStationarity confirmed (ADF rejects unit root, KPSS does not reject stationarity).\n")
}

# ---- Step 1: primary fit — mean-switching AR(1), common AR coefficient ----

SW_COMMON_AR <- c(TRUE, FALSE, TRUE) # intercept switches, AR(1) common, variance switches

fit1 <- fit_ms_ar(log_vix, sw = SW_COMMON_AR, p = 1, seed = SEED)

cat("\n=== Fit 1: mean-switching AR(1), common AR coefficient ===\n")
print(summary(fit1))

high1 <- label_high_fear_state(fit1)
low1 <- setdiff(1:2, high1)

cat(sprintf("\nHigh-fear state index: %d, low-fear state index: %d\n", high1, low1))

probs1 <- extract_probs(fit1, vix_data$date, high1)
cat(sprintf(
  "extract_probs: filt rows = %d, smooth rows (after dropping init row) = %d, dates aligned = %d\n",
  nrow(fit1@Fit@filtProb), n - 1, nrow(probs1)
))

# ---- Step 2: regime characterization (in-sample, full fit) ----------------

ts1 <- transition_stats(fit1@transMat)

cat("\n--- Fit 1: transition matrix (column-stochastic: P[i,j] = P(state_t=i | state_t-1=j)) ---\n")
print(fit1@transMat)
cat(sprintf("Column sums (should be ~1): %s\n", paste(signif(colSums(fit1@transMat), 6), collapse = ", ")))

cat("\n--- Fit 1: implied expected durations (trading days) ---\n")
cat(sprintf("State %d (low-fear):  %.1f days (%.1f weeks)\n", low1, ts1$durations[low1], ts1$durations[low1] / 5))
cat(sprintf("State %d (high-fear): %.1f days (%.1f weeks)\n", high1, ts1$durations[high1], ts1$durations[high1] / 5))

cat("\n--- Fit 1: stable (stationary) probabilities ---\n")
cat(sprintf("State %d (low-fear):  %.4f\n", low1, ts1$stable_probs[low1]))
cat(sprintf("State %d (high-fear): %.4f\n", high1, ts1$stable_probs[high1]))

coefs1 <- fit1@Coef
ar1_coef <- coefs1[["y_1"]][1] # common across states by construction (sw = FALSE)
means1 <- coefs1[["(Intercept)"]] / (1 - coefs1[["y_1"]])
sigmas1 <- fit1@std

cat("\n--- Fit 1: regime means (unconditional mean = intercept / (1 - AR)) ---\n")
cat(sprintf(
  "Low-fear:  mean log(VIX) = %.4f (VIX ~ %.2f), sigma = %.4f\n",
  means1[low1], exp(means1[low1]), sigmas1[low1]
))
cat(sprintf(
  "High-fear: mean log(VIX) = %.4f (VIX ~ %.2f), sigma = %.4f\n",
  means1[high1], exp(means1[high1]), sigmas1[high1]
))
cat(sprintf(
  "Separation: %.4f log-points (VIX ratio ~ %.2fx), common AR(1) coef = %.4f\n",
  means1[high1] - means1[low1], exp(means1[high1] - means1[low1]), ar1_coef
))

# High-fear episodes from smoothed probabilities (full-sample, in-sample only)
episodes1 <- regime_episodes(probs1$smooth_prob_high > 0.5, probs1$date, min_days = MIN_EPISODE_DAYS)

cat(sprintf("\n--- Fit 1: high-fear episodes (smoothed P > 0.5), n = %d ---\n", nrow(episodes1)))
cat(sprintf(
  "Duration (days) summary: min = %d, median = %.0f, mean = %.1f, max = %d\n",
  min(episodes1$duration_days), median(episodes1$duration_days),
  mean(episodes1$duration_days), max(episodes1$duration_days)
))
n_short <- sum(episodes1$short)
cat(sprintf(
  "Episodes shorter than %d days (implausibly short): %d / %d (%.1f%%)\n",
  MIN_EPISODE_DAYS, n_short, nrow(episodes1), 100 * n_short / nrow(episodes1)
))

cat("\n--- Coverage of known stress windows by high-fear episodes (smoothed) ---\n")
for (i in seq_len(nrow(STRESS_WINDOWS))) {
  w <- STRESS_WINDOWS[i, ]
  overlapping <- episodes1 |>
    dplyr::filter(start_date <= w$end, end_date >= w$start)

  if (nrow(overlapping) == 0) {
    cat(sprintf("%-22s [%s, %s]: NOT covered by any high-fear episode\n", w$label, w$start, w$end))
  } else {
    for (j in seq_len(nrow(overlapping))) {
      cat(sprintf(
        "%-22s [%s, %s]: covered by episode %s -> %s (%d days)%s\n",
        w$label, w$start, w$end,
        overlapping$start_date[j], overlapping$end_date[j], overlapping$duration_days[j],
        if (overlapping$short[j]) " [SHORT]" else ""
      ))
    }
  }
}

# Money plot: log(VIX) with high-fear regions shaded, log(20) reference line
money_data <- tibble::tibble(date = probs1$date, log_vix = log_vix[2:n])

p_money <- ggplot(money_data, aes(x = date, y = log_vix)) +
  geom_rect(
    data = episodes1, inherit.aes = FALSE,
    aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
    fill = "#D55E00", alpha = 0.18
  ) +
  geom_line(linewidth = 0.3) +
  geom_hline(yintercept = log(VIX_THRESHOLD), linetype = "dashed", color = "#0072B2", linewidth = 0.6) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  labs(
    title = "log(VIX) with High-Fear Regimes Shaded (Smoothed P(high-fear) > 0.5)",
    subtitle = sprintf("Dashed line = log(%d); y-axis is on a log scale", VIX_THRESHOLD),
    x = NULL, y = "log(VIX close)"
  ) +
  theme_minimal()

money_plot_path <- versioned_path(OUTPUT_DIR, "02c_money_plot", "png")
ggsave(money_plot_path, p_money, width = 12, height = 6, dpi = 150)
cat(sprintf("\nSaved: %s\n", money_plot_path))

# ---- Step 3: secondary fit — AR coefficient also switches -----------------

SW_SWITCHING_AR <- c(TRUE, TRUE, TRUE)

fit2 <- fit_ms_ar(log_vix, sw = SW_SWITCHING_AR, p = 1, seed = SEED)

cat("\n=== Fit 2: mean- and AR-switching AR(1) ===\n")
print(summary(fit2))

high2 <- label_high_fear_state(fit2)
low2 <- setdiff(1:2, high2)
ts2 <- transition_stats(fit2@transMat)

coefs2 <- fit2@Coef
ar2_coefs <- coefs2[["y_1"]]
means2 <- coefs2[["(Intercept)"]] / (1 - ar2_coefs)
sigmas2 <- fit2@std

cat("\n--- Fit 2: transition matrix ---\n")
print(fit2@transMat)
cat(sprintf(
  "Durations: low-fear (state %d) = %.1f days, high-fear (state %d) = %.1f days\n",
  low2, ts2$durations[low2], high2, ts2$durations[high2]
))
cat(sprintf(
  "AR coefficients: low-fear = %.4f, high-fear = %.4f (fit 1 common AR = %.4f)\n",
  ar2_coefs[low2], ar2_coefs[high2], ar1_coef
))

n_eff <- nrow(fit1@Fit@filtProb) # n - 1, same for both fits
logLik1 <- fit1@Fit@logLikel
logLik2 <- fit2@Fit@logLikel
k1 <- 7 # 2 intercepts + 1 common AR + 2 sigmas + 2 free transition probs
k2 <- 8 # + 1 extra AR coef
aic1 <- -2 * logLik1 + 2 * k1
aic2 <- -2 * logLik2 + 2 * k2
bic1 <- -2 * logLik1 + k1 * log(n_eff)
bic2 <- -2 * logLik2 + k2 * log(n_eff)

cat("\n--- Model comparison ---\n")
cat(sprintf("Fit 1 (common AR):    logLik = %.2f, AIC = %.2f, BIC = %.2f\n", logLik1, aic1, bic1))
cat(sprintf("Fit 2 (switching AR): logLik = %.2f, AIC = %.2f, BIC = %.2f\n", logLik2, aic2, bic2))
cat(sprintf(
  "Duration comparison (days): fit1 low/high = %.1f / %.1f, fit2 low/high = %.1f / %.1f\n",
  ts1$durations[low1], ts1$durations[high1], ts2$durations[low2], ts2$durations[high2]
))

# ---- Step 4: benchmark validation vs VIX > 20 (in-sample, descriptive) ----

vix_close_aligned <- vix_data$vix_close[2:n]
threshold_regime_insample <- label_regime(vix_close_aligned, VIX_THRESHOLD)
ms_regime_insample <- ifelse(probs1$smooth_prob_high > 0.5, "high", "low")

cm_insample <- caret::confusionMatrix(
  data = factor(ms_regime_insample, levels = c("low", "high")),
  reference = factor(threshold_regime_insample, levels = c("low", "high")),
  positive = "high"
)

cat("\n=== Step 4: in-sample MS regime (smoothed) vs VIX > 20 ===\n")
print(cm_insample$table)
cat(sprintf("Cohen's kappa:      %.3f\n", cm_insample$overall["Kappa"]))
cat(sprintf("Balanced accuracy:  %.3f\n", cm_insample$byClass["Balanced Accuracy"]))
cat(sprintf("Sensitivity (high): %.3f\n", cm_insample$byClass["Sensitivity"]))
cat(sprintf("Specificity (low):  %.3f\n", cm_insample$byClass["Specificity"]))
cat(sprintf("(Plain accuracy:    %.3f — reported for reference only)\n", cm_insample$overall["Accuracy"]))

# ---- Step 5: real-time robustness — expanding-window filtered classification ----

# Sanity check: hand-rolled Hamilton filter (hamilton_filter_update) should
# reproduce fit1's own filtProb when run forward with fit1's parameters.
# This validates the transMat / state-ordering convention before it's used
# for the OOS expanding-window filter below.
sigmas1_vec <- fit1@std
ar1_vec <- coefs1[["y_1"]] # length-2, identical values (common AR)
means1_int <- coefs1[["(Intercept)"]]
P1 <- fit1@transMat

n_filt1 <- nrow(fit1@Fit@filtProb)
hand_filt <- matrix(NA_real_, nrow = n_filt1, ncol = 2)
hand_filt[1, ] <- fit1@Fit@filtProb[1, ]
pi_cur <- hand_filt[1, ]
for (i in 2:n_filt1) {
  t_idx <- i + 1
  upd <- hamilton_filter_update(pi_cur, log_vix[t_idx], log_vix[t_idx - 1], means1_int, ar1_vec, sigmas1_vec, P1)
  pi_cur <- upd$pi_filtered
  hand_filt[i, ] <- pi_cur
}
max_diff <- max(abs(hand_filt - fit1@Fit@filtProb))
cat(sprintf("\n=== Step 5: Hamilton filter sanity check ===\n"))
cat(sprintf("Max abs diff vs fit1@Fit@filtProb: %.2e\n", max_diff))
if (max_diff > 1e-4) {
  cat("WARNING: hand-rolled filter does not reproduce MSwM's filtProb — check\n")
  cat("the transMat row/column convention before trusting OOS kappa below.\n")
} else {
  cat("OK: hand-rolled filter matches MSwM's filtProb.\n")
}

# Expanding-window backtest: refit every REFIT_EVERY days; roll the filtered
# probability forward one observation at a time between refits using
# hamilton_filter_update() with the fixed (last-refit) parameters — no
# look-ahead, O(1) per new day.
n_train <- floor(TRAIN_FRAC * n)
refit_points <- seq(n_train, n - 1, by = REFIT_EVERY)

oos_list <- vector("list", length(refit_points))
prev_fit <- NULL

for (k in seq_along(refit_points)) {
  rp <- refit_points[k]
  seg_end <- min(rp + REFIT_EVERY, n)
  if (seg_end <= rp) next

  fit_k <- tryCatch(
    fit_ms_ar(log_vix[1:rp], sw = SW_COMMON_AR, p = 1, seed = SEED),
    error = function(e) NULL
  )
  if (is.null(fit_k)) {
    message(sprintf("msmFit failed at refit point %d (obs 1:%d); reusing previous fit", k, rp))
    fit_k <- prev_fit
  } else {
    prev_fit <- fit_k
  }

  high_k <- label_high_fear_state(fit_k)
  means_k <- fit_k@Coef[["(Intercept)"]]
  ar_k <- fit_k@Coef[["y_1"]]
  sigmas_k <- fit_k@std
  P_k <- fit_k@transMat

  pi_cur <- as.numeric(fit_k@Fit@filtProb[rp - 1, ])

  lm_ar_k <- lm(y2 ~ y1, data = data.frame(y2 = log_vix[2:rp], y1 = log_vix[1:(rp - 1)]))
  ar1_coefs_k <- coef(lm_ar_k)

  seg_n <- seg_end - rp
  seg_filt <- numeric(seg_n)
  seg_ms_fc <- numeric(seg_n)
  seg_ar_fc <- numeric(seg_n)

  for (j in seq_len(seg_n)) {
    t_idx <- rp + j
    upd <- hamilton_filter_update(pi_cur, log_vix[t_idx], log_vix[t_idx - 1], means_k, ar_k, sigmas_k, P_k)
    seg_ms_fc[j] <- sum(upd$pi_predicted * (means_k + ar_k * log_vix[t_idx - 1]))
    seg_ar_fc[j] <- ar1_coefs_k[1] + ar1_coefs_k[2] * log_vix[t_idx - 1]
    pi_cur <- upd$pi_filtered
    seg_filt[j] <- pi_cur[high_k]
  }

  oos_list[[k]] <- tibble::tibble(
    date = vix_data$date[(rp + 1):seg_end],
    vix_close = vix_data$vix_close[(rp + 1):seg_end],
    log_vix = log_vix[(rp + 1):seg_end],
    filt_prob_high = seg_filt,
    ms_forecast = seg_ms_fc,
    ar1_forecast = seg_ar_fc
  )

  if (k %% 20 == 0 || k == length(refit_points)) {
    message(sprintf("Expanding window refit %d / %d (training obs 1:%d)", k, length(refit_points), rp))
  }
}

oos <- dplyr::bind_rows(oos_list) |>
  dplyr::mutate(
    ms_regime = ifelse(filt_prob_high > 0.5, "high", "low"),
    threshold_regime = label_regime(vix_close, VIX_THRESHOLD)
  )

cm_oos <- caret::confusionMatrix(
  data = factor(oos$ms_regime, levels = c("low", "high")),
  reference = factor(oos$threshold_regime, levels = c("low", "high")),
  positive = "high"
)

cat(sprintf("\n=== Step 5: real-time filtered regime (expanding window, n_oos = %d) vs VIX > 20 ===\n", nrow(oos)))
print(cm_oos$table)
cat(sprintf("Cohen's kappa:      %.3f\n", cm_oos$overall["Kappa"]))
cat(sprintf("Balanced accuracy:  %.3f\n", cm_oos$byClass["Balanced Accuracy"]))
cat(sprintf("Sensitivity (high): %.3f\n", cm_oos$byClass["Sensitivity"]))
cat(sprintf("Specificity (low):  %.3f\n", cm_oos$byClass["Specificity"]))
cat(sprintf("(Plain accuracy:    %.3f — reported for reference only)\n", cm_oos$overall["Accuracy"]))

# Filtered (real-time) vs smoothed (in-sample) over 2020, to show lag
fs_2020 <- dplyr::inner_join(
  probs1 |> dplyr::filter(date >= "2020-01-01", date <= "2020-12-31") |> dplyr::select(date, smooth_prob_high),
  oos |> dplyr::filter(date >= "2020-01-01", date <= "2020-12-31") |> dplyr::select(date, filt_prob_high),
  by = "date"
) |>
  tidyr::pivot_longer(cols = c(smooth_prob_high, filt_prob_high), names_to = "type", values_to = "prob")

p_fs <- ggplot(fs_2020, aes(x = date, y = prob, color = type)) +
  geom_line(linewidth = 0.6) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey40") +
  scale_color_manual(
    values = c(smooth_prob_high = "#0072B2", filt_prob_high = "#E69F00"),
    labels = c(smooth_prob_high = "Smoothed (in-sample)", filt_prob_high = "Filtered (real-time)")
  ) +
  labs(
    title = "P(High-Fear Regime): Filtered (Real-Time) vs Smoothed (In-Sample), 2020",
    x = NULL, y = "P(high-fear)", color = NULL
  ) +
  theme_minimal()

fs_plot_path <- versioned_path(OUTPUT_DIR, "02c_filtered_vs_smoothed_2020", "png")
ggsave(fs_plot_path, p_fs, width = 12, height = 6, dpi = 150)
cat(sprintf("\nSaved: %s\n", fs_plot_path))

# Lag at known regime shifts: first date filtered crosses 0.5 vs smoothed
first_cross <- function(df, prob_col, start_date, end_date) {
  sub <- df |>
    dplyr::filter(date >= as.Date(start_date), date <= as.Date(end_date), .data[[prob_col]] > 0.5)
  if (nrow(sub) == 0) return(NA)
  min(sub$date)
}

lag_windows <- tibble::tribble(
  ~label, ~start, ~end,
  "Feb 2018 Volmageddon", "2018-01-15", "2018-03-15",
  "2020 COVID onset", "2020-02-01", "2020-04-01"
)

cat("\n--- Filtered (real-time) vs smoothed (in-sample) lag at known regime shifts ---\n")
for (i in seq_len(nrow(lag_windows))) {
  w <- lag_windows[i, ]
  smooth_date <- first_cross(probs1, "smooth_prob_high", w$start, w$end)
  filt_date <- first_cross(oos, "filt_prob_high", w$start, w$end)
  if (is.na(smooth_date) || is.na(filt_date)) {
    cat(sprintf("%-22s smoothed = %s, filtered = %s (one or both never crossed 0.5 in window)\n",
      w$label, smooth_date, filt_date))
  } else {
    lag_days <- as.numeric(filt_date - smooth_date)
    cat(sprintf(
      "%-22s smoothed crosses 0.5 on %s, filtered crosses on %s (lag = %d days)\n",
      w$label, smooth_date, filt_date, lag_days
    ))
  }
}

# ---- Step 6 (secondary): one-step log(VIX) forecast RMSE ------------------

rmse_ms <- sqrt(mean((oos$log_vix - oos$ms_forecast)^2))
rmse_ar <- sqrt(mean((oos$log_vix - oos$ar1_forecast)^2))

cat("\n=== Step 6 (secondary): one-step log(VIX) forecast RMSE, OOS ===\n")
cat(sprintf("MS-AR (regime-weighted): RMSE = %.5f\n", rmse_ms))
cat(sprintf("Plain AR(1):             RMSE = %.5f\n", rmse_ar))
cat("(Tiebreaker only — point forecasts of a persistent level are AR-dominated.)\n")

# ---- Save summary tables ----------------------------------------------------

summarize_fit <- function(fit, model_name) {
  high <- label_high_fear_state(fit)
  low <- setdiff(1:2, high)
  order_idx <- c(low, high)

  coefs <- fit@Coef
  intercepts <- coefs[["(Intercept)"]]
  ar_coefs <- coefs[["y_1"]]
  means <- intercepts / (1 - ar_coefs)
  sigmas <- fit@std
  ts <- transition_stats(fit@transMat)

  logLik <- fit@Fit@logLikel
  n_eff_local <- nrow(fit@Fit@filtProb)
  k_params <- if (length(unique(round(ar_coefs, 10))) == 1) 7 else 8
  aic <- -2 * logLik + 2 * k_params
  bic <- -2 * logLik + k_params * log(n_eff_local)

  tibble::tibble(
    model = model_name,
    regime = c("low_fear", "high_fear"),
    mean_log_vix = means[order_idx],
    mean_vix = exp(means[order_idx]),
    ar_coef = ar_coefs[order_idx],
    sigma = sigmas[order_idx],
    duration_days = ts$durations[order_idx],
    stable_prob = ts$stable_probs[order_idx],
    logLik = logLik,
    AIC = aic,
    BIC = bic
  )
}

regime_summary <- dplyr::bind_rows(
  summarize_fit(fit1, "fit1_common_AR"),
  summarize_fit(fit2, "fit2_switching_AR")
)

regime_summary_path <- versioned_path(OUTPUT_DIR, "02c_regime_summary", "csv")
readr::write_csv(regime_summary, regime_summary_path)
cat(sprintf("\nSaved: %s\n", regime_summary_path))

benchmark_summary <- tibble::tibble(
  validation = c("in_sample_smoothed", "real_time_filtered_oos"),
  n = c(nrow(probs1), nrow(oos)),
  kappa = c(cm_insample$overall["Kappa"], cm_oos$overall["Kappa"]),
  balanced_accuracy = c(cm_insample$byClass["Balanced Accuracy"], cm_oos$byClass["Balanced Accuracy"]),
  sensitivity = c(cm_insample$byClass["Sensitivity"], cm_oos$byClass["Sensitivity"]),
  specificity = c(cm_insample$byClass["Specificity"], cm_oos$byClass["Specificity"]),
  accuracy = c(cm_insample$overall["Accuracy"], cm_oos$overall["Accuracy"])
)

benchmark_summary_path <- versioned_path(OUTPUT_DIR, "02c_benchmark_metrics", "csv")
readr::write_csv(benchmark_summary, benchmark_summary_path)
cat(sprintf("Saved: %s\n", benchmark_summary_path))

forecast_rmse <- tibble::tibble(model = c("MS-AR", "AR(1)"), rmse = c(rmse_ms, rmse_ar))
forecast_rmse_path <- versioned_path(OUTPUT_DIR, "02c_forecast_rmse", "csv")
readr::write_csv(forecast_rmse, forecast_rmse_path)
cat(sprintf("Saved: %s\n", forecast_rmse_path))

# ---- Verdict ----------------------------------------------------------------

cat("\n========================== VERDICT ==========================\n")

separated <- abs(means1[high1] - means1[low1]) > 0.3 # ~35% VIX ratio
persistent <- ts1$durations[low1] >= 5 & ts1$durations[high1] >= 5 # >= 1 week
kappa_oos <- cm_oos$overall["Kappa"]
value_add <- kappa_oos > 0.1

cat(sprintf(
  "Regime separation: %.4f log-points (VIX ratio %.2fx) -> %s\n",
  means1[high1] - means1[low1], exp(means1[high1] - means1[low1]),
  if (separated) "clearly separated" else "weakly separated"
))
cat(sprintf(
  "Persistence: low-fear %.1f days, high-fear %.1f days -> %s\n",
  ts1$durations[low1], ts1$durations[high1],
  if (persistent) "persistent (>= ~1 week)" else "NOT persistent (sub-week)"
))
cat(sprintf(
  "Real-time filtered kappa vs VIX>20: %.3f -> %s\n",
  kappa_oos,
  if (value_add) "meaningfully > 0" else "near zero / weak"
))

if (separated && persistent && value_add) {
  cat("\nOVERALL: Level-based MS-AR regimes appear to have legs — separated,\n")
  cat("persistent, and the real-time filtered classification adds signal\n")
  cat("beyond the VIX > 20 threshold. Worth pursuing further.\n")
} else {
  cat("\nOVERALL: SOFT-NEGATIVE. One or more criteria (separation, persistence,\n")
  cat("real-time value-add over VIX>20) were not clearly met — see details above\n")
  cat("before treating this as a usable real-time fear-regime signal.\n")
}
cat("================================================================\n")
