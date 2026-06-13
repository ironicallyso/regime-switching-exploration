# Phase 2.1: Out-of-sample one-step-ahead variance forecast comparison.
#
# Go/no-go check: does the 2-regime MS-GARCH model produce better OOS
# variance forecasts than single-regime baselines, fit/forecast with the
# same MSGARCH interface (apples-to-apples)? Not part of the SPEC.md phase
# sequence (which covers regime-probability forecasting) — this is a
# standalone diagnostic on forecast accuracy.

library(tidyverse)
library(MSGARCH)

source("R/data.R")
source("R/oos_forecast.R")

# ---- Parameters -------------------------------------------------------

TRAIN_FRAC <- 0.6 # initial training window as a fraction of the full series
REFIT_EVERY <- 21 # refit cadence in trading days (~1 month)
SEED <- 2024
OUTPUT_DIR <- "outputs"

MODEL_LABELS <- c(
  sgarch_norm = "sGARCH (normal)",
  sgarch_std = "sGARCH (Student-t)",
  msgarch_2s_norm = "2-regime MS-GARCH (normal)"
)

# Okabe-Ito colorblind-safe palette
MODEL_COLORS <- c(
  sgarch_norm = "#0072B2",
  sgarch_std = "#E69F00",
  msgarch_2s_norm = "#009E73"
)

# ---- Data ---------------------------------------------------------------

vix_data <- fetch_vix()

cat("Series: vix_log_return (decimal units, e.g. 0.01 = 1%)\n")
cat(sprintf(
  "Date range: %s to %s (n = %d)\n",
  min(vix_data$date), max(vix_data$date), nrow(vix_data)
))
cat(sprintf("Mean log-return: %.6f\n", mean(vix_data$vix_log_return)))

# No realized-variance series exists in this project, so the OOS variance
# proxy is squared demeaned returns: (r_t - mean(r))^2.

# ---- Specs & expanding-window backtest -----------------------------------

specs <- build_specs()
model_names <- names(specs)

set.seed(SEED)
oos <- run_expanding_window(
  y = vix_data$vix_log_return,
  specs = specs,
  train_frac = TRAIN_FRAC,
  refit_every = REFIT_EVERY,
  seed = SEED
)

n <- nrow(vix_data)
n_train <- floor(TRAIN_FRAC * n)
oos_dates <- vix_data$date[(n_train + 1):n]
oos <- oos |> mutate(date = oos_dates, .before = 1)

cat(sprintf(
  "OOS period: %s to %s (n = %d, %.1f%% of series), refit every %d days\n",
  min(oos$date), max(oos$date), nrow(oos),
  100 * nrow(oos) / n, REFIT_EVERY
))

# ---- Losses & results table -----------------------------------------------

losses <- compute_losses(oos, model_names)
results <- summarize_results(losses)

dm_rows <- model_names[model_names != "msgarch_2s_norm"] |>
  purrr::map_dfr(function(baseline) {
    target_qlike <- losses |> filter(model == "msgarch_2s_norm") |> pull(qlike)
    baseline_qlike <- losses |> filter(model == baseline) |> pull(qlike)
    dm <- dm_test(target_qlike, baseline_qlike)
    tibble(
      model = baseline,
      dm_stat_vs_msgarch = dm$statistic,
      dm_pvalue_vs_msgarch = dm$p.value
    )
  })

results <- results |>
  left_join(dm_rows, by = "model") |>
  mutate(model_label = MODEL_LABELS[model], .after = model)

print(results)

results_path <- versioned_path(OUTPUT_DIR, "02b_oos_results", "csv")
write_csv(results, results_path)
cat(sprintf("Results table written to %s\n", results_path))

# ---- Plot 1: OOS forecast volatility vs proxy -----------------------------

vol_long <- oos |>
  pivot_longer(
    cols = all_of(model_names),
    names_to = "model",
    values_to = "forecast_var"
  ) |>
  mutate(
    forecast_vol = sqrt(forecast_var),
    proxy_vol = sqrt(proxy),
    model_label = factor(MODEL_LABELS[model], levels = MODEL_LABELS)
  )

# Squared-return proxy spikes are far larger than any model's forecast
# (it's a single-day realization, not a smoothed estimate). Cap the y-axis
# to the forecast range so the model lines stay legible; the proxy area is
# clipped at the top rather than stretching the whole plot.
vol_ylim <- quantile(vol_long$forecast_vol, 0.995, na.rm = TRUE)

p_vol <- ggplot(vol_long, aes(x = date)) +
  geom_area(aes(y = proxy_vol), fill = "grey70", alpha = 0.3) +
  geom_line(aes(y = forecast_vol, color = model_label), linewidth = 0.4, alpha = 0.9) +
  scale_color_manual(values = setNames(MODEL_COLORS[model_names], MODEL_LABELS[model_names])) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  coord_cartesian(ylim = c(0, vol_ylim)) +
  labs(
    title = "One-step-ahead volatility forecasts vs. realized-return proxy",
    subtitle = sprintf(
      "Models: %s\nOOS window: %s to %s | grey area = squared-return proxy",
      paste(MODEL_LABELS[model_names], collapse = " / "),
      min(oos$date), max(oos$date)
    ),
    x = NULL,
    y = "Volatility (decimal, sqrt of variance)",
    color = NULL
  ) +
  theme_minimal()

vol_plot_path <- versioned_path(OUTPUT_DIR, "02b_oos_vol_forecasts", "png")
ggsave(vol_plot_path, p_vol, width = 11, height = 6, dpi = 150)
cat(sprintf("Volatility forecast plot written to %s\n", vol_plot_path))

# ---- Plot 2: cumulative QLIKE loss differential ----------------------------

qlike_wide <- losses |>
  select(date, model, qlike) |>
  pivot_wider(names_from = model, values_from = qlike)

cum_diff <- qlike_wide |>
  arrange(date) |>
  mutate(
    across(
      all_of(model_names[model_names != "msgarch_2s_norm"]),
      ~ cumsum(msgarch_2s_norm - .x),
      .names = "cum_{.col}"
    )
  ) |>
  select(date, starts_with("cum_")) |>
  pivot_longer(
    cols = starts_with("cum_"),
    names_to = "baseline",
    names_prefix = "cum_",
    values_to = "cum_qlike_diff"
  ) |>
  mutate(baseline_label = factor(MODEL_LABELS[baseline], levels = MODEL_LABELS[model_names != "msgarch_2s_norm"]))

BASELINE_COLORS <- setNames(
  MODEL_COLORS[model_names[model_names != "msgarch_2s_norm"]],
  MODEL_LABELS[model_names[model_names != "msgarch_2s_norm"]]
)

p_cum_full <- ggplot(cum_diff, aes(x = date, y = cum_qlike_diff, color = baseline_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.5, alpha = 0.9) +
  scale_color_manual(values = BASELINE_COLORS) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  labs(
    title = "Cumulative QLIKE loss differential: MS-GARCH minus baseline",
    subtitle = "Downward slope = MS-GARCH winning over that period",
    x = NULL,
    y = "Cumulative QLIKE diff",
    color = "Baseline"
  ) +
  theme_minimal()

p_cum_2020 <- p_cum_full %+%
  filter(cum_diff, date >= as.Date("2020-01-01"), date <= as.Date("2020-12-31")) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(title = "2020 (COVID) zoom", subtitle = NULL) +
  theme(legend.position = "none")

p_cum <- cowplot::plot_grid(p_cum_full, p_cum_2020, ncol = 1, rel_heights = c(2, 1))

cum_plot_path <- versioned_path(OUTPUT_DIR, "02b_cumulative_qlike_diff", "png")
ggsave(cum_plot_path, p_cum, width = 11, height = 6, dpi = 150)
cat(sprintf("Cumulative QLIKE plot written to %s\n", cum_plot_path))

# ---- Verdict ----------------------------------------------------------------

msgarch_qlike <- results |> filter(model == "msgarch_2s_norm") |> pull(mean_qlike)
baselines <- results |> filter(model != "msgarch_2s_norm")
beats <- setNames(msgarch_qlike < baselines$mean_qlike, baselines$model)
stronger_baseline <- baselines$model[which.min(baselines$mean_qlike)]
stronger_dm_p <- baselines |> filter(model == stronger_baseline) |> pull(dm_pvalue_vs_msgarch)

cat("\n---- VERDICT ----\n")
cat(sprintf(
  "MS-GARCH (2-regime, normal) mean QLIKE: %.5f\n", msgarch_qlike
))
for (b in baselines$model) {
  cat(sprintf(
    "%s mean QLIKE: %.5f (MS-GARCH %s, DM p = %.3f)\n",
    MODEL_LABELS[b],
    baselines |> filter(model == b) |> pull(mean_qlike),
    if (beats[b]) "lower (better)" else "not lower",
    baselines |> filter(model == b) |> pull(dm_pvalue_vs_msgarch)
  ))
}
if (all(beats)) {
  cat(sprintf(
    "MS-GARCH beats BOTH single-regime baselines on mean QLIKE. Against the\nstronger baseline (%s), DM p = %.3f — %s.\n",
    MODEL_LABELS[stronger_baseline], stronger_dm_p,
    if (stronger_dm_p < 0.10) "significant at the 10% level, two regimes have legs" else "not significant at 10%, the edge is suggestive but not yet conclusive"
  ))
} else {
  cat("MS-GARCH does NOT beat both single-regime baselines on mean QLIKE.\n")
  cat(sprintf(
    "It ties or loses to %s. Regime-switching is not justified by this OOS\ncomparison alone — do not spin this as a win for two regimes.\n",
    paste(MODEL_LABELS[baselines$model[!beats]], collapse = " and ")
  ))
}
