# Phase 2.3: hysteresis (Schmitt-trigger) reliability pass on the Phase 2.2
# level-based high-fear regime signal.
#
# Phase 2.2 (02c_level_regime.R) produced filt_prob_high, the real-time
# expanding-window filtered P(high-fear), and showed that a hard
# filt_prob_high > 0.5 classifier is whipsaw-prone (high-fear expected
# duration ~6.3 days; 49% of high-fear episodes < 5 days). This phase keeps
# filt_prob_high as the underlying object and tests whether a two-threshold
# Schmitt trigger (enter high-fear above U, exit below L, hold in between)
# produces a materially more stable discrete signal, and at what cost in
# detection lag. Model-quality / signal-stability only — no refitting, no
# underlying/ARIMAX/trading logic.
#
# Inputs (reused from Phase 2.2, not recomputed):
#   outputs/02c_oos_filtered.csv      real-time expanding-window filtered
#                                      P(high-fear), OOS ~2011-2026 (the
#                                      "honest" signal; stability + lag use
#                                      this)
#   outputs/02c_full_sample_probs.csv full-sample filtered + smoothed
#                                      P(high-fear), 1990-2026 (carries
#                                      look-ahead in parameters; used only
#                                      for descriptive coverage / as a lag
#                                      reference)

library(tidyverse)

source("R/level_regime.R") # regime_episodes()
source("R/hysteresis.R") # schmitt_trigger(), regime_stability_metrics(), find_whipsaw_window()
source("R/oos_forecast.R") # versioned_path()

# ---- Parameters --------------------------------------------------------------

OUTPUT_DIR <- "outputs"
MIN_EPISODE_DAYS <- 5 # matches Phase 2.2's "implausibly short" threshold
ONSET_LAG_THRESHOLD_DAYS <- 10 # VERDICT cutoff for "a few days" of onset lag

# Pre-registered hysteresis grid: (label, upper, lower). "hard_0.5" is the
# Phase 2.2 baseline; classified via prob > 0.5 directly (see
# classify_is_high() below and R/hysteresis.R for why schmitt_trigger() is
# not used for this row).
HYSTERESIS_GRID <- tibble::tribble(
  ~label,     ~upper, ~lower,
  "hard_0.5",   0.5,    0.5,
  "0.6/0.4",    0.6,    0.4,
  "0.7/0.3",    0.7,    0.3,
  "0.8/0.2",    0.8,    0.2
)
MAIN_UL_LABEL <- "0.7/0.3" # headline U/L for plots 1 and 2

STRESS_WINDOWS <- tibble::tribble(
  ~label, ~start, ~end,
  "2008 GFC", "2008-09-01", "2009-03-31",
  "Feb 2018 Volmageddon", "2018-02-01", "2018-02-28",
  "2020 COVID", "2020-02-15", "2020-04-30",
  "2022", "2022-01-01", "2022-12-31"
) |>
  dplyr::mutate(start = as.Date(start), end = as.Date(end))

LAG_WINDOWS <- tibble::tribble(
  ~label, ~start, ~end,
  "Feb 2018 Volmageddon", "2018-01-15", "2018-03-15",
  "2020 COVID onset", "2020-02-01", "2020-04-01"
) |>
  dplyr::mutate(start = as.Date(start), end = as.Date(end))

# ---- Load Phase 2.2 outputs ----------------------------------------------------

oos <- readr::read_csv("outputs/02c_oos_filtered.csv", show_col_types = FALSE)
probs1 <- readr::read_csv("outputs/02c_full_sample_probs.csv", show_col_types = FALSE)

cat(sprintf(
  "Loaded real-time series (oos): n = %d, %s to %s\n",
  nrow(oos), min(oos$date), max(oos$date)
))
cat(sprintf(
  "Loaded full-sample series (probs1): n = %d, %s to %s\n",
  nrow(probs1), min(probs1$date), max(probs1$date)
))
cat("NOTE: probs1$filt_prob_high and probs1$smooth_prob_high carry\n")
cat("look-ahead-in-parameters (full-sample fit); used here only as a\n")
cat("coverage/lag REFERENCE, not as a real-time-usable signal.\n")

# ---- Helpers -------------------------------------------------------------------

# Classify a probability vector under a HYSTERESIS_GRID row. hard_0.5 uses
# prob > 0.5 directly to avoid the boundary ambiguity schmitt_trigger() would
# have at exactly upper == lower == 0.5.
classify_is_high <- function(prob, label, upper, lower) {
  if (label == "hard_0.5") prob > 0.5 else schmitt_trigger(prob, upper, lower)
}

# Slugify a label for use in column-name suffixes, e.g.
# "Feb 2018 Volmageddon" -> "feb_2018_volmageddon".
slug <- function(x) {
  s <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("^_|_$", "", s)
}

# First date within [w_start, w_end] where is_high is TRUE (onset), and the
# first date after that episode ends where is_high returns to FALSE (exit).
# Uses min_days = 0 in regime_episodes() so no episode is dropped as "short".
onset_exit_in_window <- function(is_high, dates, w_start, w_end) {
  in_window <- dates >= w_start & dates <= w_end
  if (!any(is_high[in_window])) {
    return(list(onset = as.Date(NA), exit = as.Date(NA)))
  }

  onset_idx <- which(in_window & is_high)[1]
  onset <- dates[onset_idx]

  episodes <- regime_episodes(is_high, dates, min_days = 0)
  ep <- episodes |> dplyr::filter(start_date <= onset, end_date >= onset)
  if (nrow(ep) == 0) return(list(onset = onset, exit = as.Date(NA)))

  exit_idx <- which(dates == ep$end_date[1]) + 1
  exit_date <- if (exit_idx <= length(dates)) dates[exit_idx] else as.Date(NA)

  list(onset = onset, exit = exit_date)
}

# ---- A) Stability metrics (real-time series) ------------------------------------

stability_rows <- vector("list", nrow(HYSTERESIS_GRID))
for (g in seq_len(nrow(HYSTERESIS_GRID))) {
  grid_label <- HYSTERESIS_GRID$label[g]
  is_high <- classify_is_high(oos$filt_prob_high, grid_label, HYSTERESIS_GRID$upper[g], HYSTERESIS_GRID$lower[g])

  stab <- regime_stability_metrics(is_high, oos$date, min_days = MIN_EPISODE_DAYS)
  names(stab) <- paste0("rt_", names(stab))

  stability_rows[[g]] <- dplyr::bind_cols(tibble::tibble(label = grid_label), stab)
}
stability_results <- dplyr::bind_rows(stability_rows)

cat("\n=== A) Stability metrics (real-time, expanding-window filtered) ===\n")
print(as.data.frame(stability_results))

# ---- B) Coverage of known stress windows (full-sample series) -------------------

stress_slugs <- slug(STRESS_WINDOWS$label)

coverage_rows <- vector("list", nrow(HYSTERESIS_GRID))
for (g in seq_len(nrow(HYSTERESIS_GRID))) {
  grid_label <- HYSTERESIS_GRID$label[g]
  is_high_fs <- classify_is_high(probs1$filt_prob_high, grid_label, HYSTERESIS_GRID$upper[g], HYSTERESIS_GRID$lower[g])
  episodes_fs <- regime_episodes(is_high_fs, probs1$date, min_days = MIN_EPISODE_DAYS)

  row <- tibble::tibble(label = grid_label)
  for (i in seq_len(nrow(STRESS_WINDOWS))) {
    w <- STRESS_WINDOWS[i, ]
    overlapping <- episodes_fs |>
      dplyr::filter(start_date <= w$end, end_date >= w$start)

    flagged <- nrow(overlapping) > 0
    dur <- if (flagged) max(overlapping$duration_days) else NA_real_

    row[[paste0("cov_flagged_", stress_slugs[i])]] <- flagged
    row[[paste0("cov_duration_days_", stress_slugs[i])]] <- dur
  }
  coverage_rows[[g]] <- row
}
coverage_results <- dplyr::bind_rows(coverage_rows)

cat("\n=== B) Coverage of known stress windows (full-sample filtered, look-ahead in params) ===\n")
print(as.data.frame(coverage_results))

# ---- C) Detection lag (real-time series, vs hard-0.5 and vs smoothed) -----------

smoothed_is_high <- probs1$smooth_prob_high > 0.5
hard05_is_high_rt <- oos$filt_prob_high > 0.5

lag_window_slugs <- slug(LAG_WINDOWS$label)

lag_rows <- list()
idx <- 1
for (g in seq_len(nrow(HYSTERESIS_GRID))) {
  grid_label <- HYSTERESIS_GRID$label[g]
  is_high_rt <- classify_is_high(oos$filt_prob_high, grid_label, HYSTERESIS_GRID$upper[g], HYSTERESIS_GRID$lower[g])

  for (w in seq_len(nrow(LAG_WINDOWS))) {
    w_start <- LAG_WINDOWS$start[w]
    w_end <- LAG_WINDOWS$end[w]

    oe_rt <- onset_exit_in_window(is_high_rt, oos$date, w_start, w_end)
    oe_hard05 <- onset_exit_in_window(hard05_is_high_rt, oos$date, w_start, w_end)
    oe_smooth <- onset_exit_in_window(smoothed_is_high, probs1$date, w_start, w_end)

    lag_rows[[idx]] <- tibble::tibble(
      label = grid_label,
      window_slug = lag_window_slugs[w],
      onset_date = oe_rt$onset,
      exit_date = oe_rt$exit,
      lag_onset_vs_hard05_days = as.numeric(oe_rt$onset - oe_hard05$onset),
      lag_exit_vs_hard05_days = as.numeric(oe_rt$exit - oe_hard05$exit),
      lag_onset_vs_smoothed_days = as.numeric(oe_rt$onset - oe_smooth$onset),
      lag_exit_vs_smoothed_days = as.numeric(oe_rt$exit - oe_smooth$exit)
    )
    idx <- idx + 1
  }
}
lag_results <- dplyr::bind_rows(lag_rows)

lag_wide <- lag_results |>
  tidyr::pivot_wider(
    id_cols = label,
    names_from = window_slug,
    values_from = c(
      onset_date, exit_date,
      lag_onset_vs_hard05_days, lag_exit_vs_hard05_days,
      lag_onset_vs_smoothed_days, lag_exit_vs_smoothed_days
    ),
    names_glue = "{.value}_{window_slug}"
  )

cat("\n=== C) Detection lag at known regime-shift onsets (real-time vs hard-0.5 / smoothed) ===\n")
print(as.data.frame(lag_results))

# ---- Assemble metrics table ------------------------------------------------------

metrics <- HYSTERESIS_GRID |>
  dplyr::left_join(stability_results, by = "label") |>
  dplyr::left_join(coverage_results, by = "label") |>
  dplyr::left_join(lag_wide, by = "label")

cat("\n=== Phase 2.3 hysteresis metrics table (rows = classifier) ===\n")
print(as.data.frame(metrics))

metrics_path <- versioned_path(OUTPUT_DIR, "02d_hysteresis_metrics", "csv")
readr::write_csv(metrics, metrics_path)
cat(sprintf("\nSaved: %s\n", metrics_path))

# ---- Plot 1: filt_prob_high with hysteresis band and shaded episodes -----------

main_row <- HYSTERESIS_GRID |> dplyr::filter(label == MAIN_UL_LABEL)
is_high_main <- schmitt_trigger(oos$filt_prob_high, main_row$upper, main_row$lower)
episodes_main <- regime_episodes(is_high_main, oos$date, min_days = MIN_EPISODE_DAYS)

p_thresh <- ggplot(oos, aes(x = date, y = filt_prob_high)) +
  geom_rect(
    data = episodes_main, inherit.aes = FALSE,
    aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
    fill = "#D55E00", alpha = 0.18
  ) +
  geom_line(linewidth = 0.3, color = "#0072B2") +
  geom_hline(yintercept = main_row$upper, linetype = "dashed", color = "#E69F00", linewidth = 0.6) +
  geom_hline(yintercept = main_row$lower, linetype = "dashed", color = "#009E73", linewidth = 0.6) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title = sprintf("Real-Time Filtered P(High-Fear) with Hysteresis Band (U = %.1f, L = %.1f)", main_row$upper, main_row$lower),
    subtitle = "Orange dashed = upper threshold (enter high-fear); green dashed = lower threshold (exit high-fear); red shading = resulting high-fear episodes",
    x = NULL, y = "P(high-fear), filtered (real-time, expanding window)"
  ) +
  theme_minimal()

p1_path <- versioned_path(OUTPUT_DIR, "02d_filt_prob_with_thresholds", "png")
ggsave(p1_path, p_thresh, width = 12, height = 6, dpi = 150)
cat(sprintf("\nSaved: %s\n", p1_path))

# ---- Plot 2/3: hard-0.5 vs hysteresis regime shading, full series + zoom --------

shading_data <- dplyr::bind_rows(
  oos |> dplyr::mutate(classifier = "hard_0.5", is_high = hard05_is_high_rt),
  oos |> dplyr::mutate(classifier = MAIN_UL_LABEL, is_high = is_high_main)
) |>
  dplyr::mutate(classifier = factor(classifier, levels = c("hard_0.5", MAIN_UL_LABEL)))

episodes_hard05 <- regime_episodes(hard05_is_high_rt, oos$date, min_days = MIN_EPISODE_DAYS) |>
  dplyr::mutate(classifier = "hard_0.5")
episodes_main_labeled <- episodes_main |> dplyr::mutate(classifier = MAIN_UL_LABEL)
episodes_all <- dplyr::bind_rows(episodes_hard05, episodes_main_labeled) |>
  dplyr::mutate(classifier = factor(classifier, levels = c("hard_0.5", MAIN_UL_LABEL)))

p_full <- ggplot(shading_data, aes(x = date, y = filt_prob_high)) +
  geom_rect(
    data = episodes_all, inherit.aes = FALSE,
    aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
    fill = "#D55E00", alpha = 0.18
  ) +
  geom_line(linewidth = 0.3, color = "#0072B2") +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey40") +
  facet_wrap(~classifier, ncol = 1) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title = sprintf("High-Fear Regime Shading: Hard-0.5 vs Hysteresis (U/L = %.1f/%.1f)", main_row$upper, main_row$lower),
    subtitle = "Real-time expanding-window filtered P(high-fear); dotted line = 0.5; red shading = high-fear episodes per classifier",
    x = NULL, y = "P(high-fear), filtered (real-time)"
  ) +
  theme_minimal()

p2_path <- versioned_path(OUTPUT_DIR, "02d_regime_shading_comparison", "png")
ggsave(p2_path, p_full, width = 12, height = 8, dpi = 150)
cat(sprintf("\nSaved: %s\n", p2_path))

whipsaw <- find_whipsaw_window(hard05_is_high_rt, oos$date, window_days = 252)
cat(sprintf(
  "\nWhipsaw zoom window (most hard-0.5 transitions in ~1 trading year): %s to %s\n",
  whipsaw$start, whipsaw$end
))

shading_zoom <- shading_data |> dplyr::filter(date >= whipsaw$start, date <= whipsaw$end)
episodes_zoom <- episodes_all |> dplyr::filter(start_date <= whipsaw$end, end_date >= whipsaw$start)

p_zoom <- ggplot(shading_zoom, aes(x = date, y = filt_prob_high)) +
  geom_rect(
    data = episodes_zoom, inherit.aes = FALSE,
    aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
    fill = "#D55E00", alpha = 0.18
  ) +
  geom_line(linewidth = 0.4, color = "#0072B2") +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey40") +
  facet_wrap(~classifier, ncol = 1) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(
    title = "High-Fear Regime Shading: Zoom on Whipsaw Window",
    subtitle = sprintf(
      "%s to %s — the ~1-trading-year window with the most hard-0.5 state transitions",
      whipsaw$start, whipsaw$end
    ),
    x = NULL, y = "P(high-fear), filtered (real-time)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3_path <- versioned_path(OUTPUT_DIR, "02d_regime_shading_comparison_zoom", "png")
ggsave(p3_path, p_zoom, width = 12, height = 8, dpi = 150)
cat(sprintf("Saved: %s\n", p3_path))

# ---- Plot 4: episode duration histograms, hard-0.5 vs hysteresis grid ----------

duration_rows <- list()
for (g in seq_len(nrow(HYSTERESIS_GRID))) {
  grid_label <- HYSTERESIS_GRID$label[g]
  is_high <- classify_is_high(oos$filt_prob_high, grid_label, HYSTERESIS_GRID$upper[g], HYSTERESIS_GRID$lower[g])
  episodes <- regime_episodes(is_high, oos$date, min_days = MIN_EPISODE_DAYS)
  if (nrow(episodes) > 0) {
    duration_rows[[g]] <- tibble::tibble(label = grid_label, duration_days = episodes$duration_days)
  }
}
duration_data <- dplyr::bind_rows(duration_rows) |>
  dplyr::mutate(label = factor(label, levels = HYSTERESIS_GRID$label))

p_hist <- ggplot(duration_data, aes(x = duration_days)) +
  geom_histogram(bins = 30, fill = "#0072B2", color = "white") +
  geom_vline(xintercept = MIN_EPISODE_DAYS, linetype = "dashed", color = "#D55E00") +
  facet_wrap(~label, ncol = 2) +
  labs(
    title = "High-Fear Episode Duration: Hard-0.5 vs Hysteresis Grid (Real-Time Series)",
    subtitle = sprintf("Dashed line = %d-day 'short episode' threshold used in pct_episodes_short", MIN_EPISODE_DAYS),
    x = "Episode duration (trading days)", y = "Count"
  ) +
  theme_minimal()

p4_path <- versioned_path(OUTPUT_DIR, "02d_episode_duration_histograms", "png")
ggsave(p4_path, p_hist, width = 12, height = 8, dpi = 150)
cat(sprintf("Saved: %s\n", p4_path))

# ---- Verdict -----------------------------------------------------------------

cat("\n========================== VERDICT ==========================\n")

hard05 <- metrics |> dplyr::filter(label == "hard_0.5")

cat(sprintf(
  "Hard-0.5 baseline (real-time): %d episodes (%.1f/yr), median duration %.1f days, %.1f%% short (<%d days), %.1f%% time high-fear\n",
  hard05$rt_n_episodes, hard05$rt_episodes_per_year, hard05$rt_median_duration_days,
  hard05$rt_pct_episodes_short, MIN_EPISODE_DAYS, hard05$rt_pct_time_high_fear
))

stress_flag_cols <- paste0("cov_flagged_", stress_slugs)
onset_lag_cols <- paste0("lag_onset_vs_hard05_days_", lag_window_slugs)

candidates <- character(0)

for (lbl in setdiff(HYSTERESIS_GRID$label, "hard_0.5")) {
  row <- metrics |> dplyr::filter(label == lbl)

  whipsaw_cut <- row$rt_pct_episodes_short < hard05$rt_pct_episodes_short &&
    row$rt_episodes_per_year < hard05$rt_episodes_per_year

  coverage_ok <- all(unlist(row[stress_flag_cols]))

  onset_lags <- unlist(row[onset_lag_cols])
  max_onset_lag <- max(abs(onset_lags), na.rm = TRUE)
  lag_ok <- max_onset_lag <= ONSET_LAG_THRESHOLD_DAYS

  cat(sprintf(
    "\n%-10s: %d episodes (%.1f/yr), median duration %.1f days, %.1f%% short, %.1f%% time high-fear\n",
    lbl, row$rt_n_episodes, row$rt_episodes_per_year, row$rt_median_duration_days,
    row$rt_pct_episodes_short, row$rt_pct_time_high_fear
  ))
  cat(sprintf(
    "  Whipsaw reduction vs hard-0.5: %s (episodes/yr %.1f -> %.1f, %%short %.1f%% -> %.1f%%)\n",
    if (whipsaw_cut) "YES" else "no", hard05$rt_episodes_per_year, row$rt_episodes_per_year,
    hard05$rt_pct_episodes_short, row$rt_pct_episodes_short
  ))
  cat(sprintf("  Coverage of all 4 stress windows: %s\n", if (coverage_ok) "YES (all flagged)" else "NO (gap)"))
  cat(sprintf(
    "  Max onset lag vs hard-0.5 (Feb 2018 / COVID 2020): %.0f days -> %s (threshold <= %d)\n",
    max_onset_lag, if (lag_ok) "OK" else "TOO SLOW", ONSET_LAG_THRESHOLD_DAYS
  ))

  if (whipsaw_cut && coverage_ok && lag_ok) {
    candidates <- c(candidates, lbl)
  }
}

cat("\n--------------------------------------------------------------\n")
if (length(candidates) > 0) {
  rec_label <- candidates[1]
  rec <- metrics |> dplyr::filter(label == rec_label)
  rec_max_onset_lag <- max(abs(unlist(rec[onset_lag_cols])), na.rm = TRUE)

  cat(sprintf(
    "RECOMMENDATION: %s — cuts whipsaw (episodes/yr %.1f -> %.1f, %%short %.1f%% -> %.1f%%),\n",
    rec_label, hard05$rt_episodes_per_year, rec$rt_episodes_per_year,
    hard05$rt_pct_episodes_short, rec$rt_pct_episodes_short
  ))
  cat(sprintf(
    "covers all 4 stress windows, at an onset-lag cost of %.0f days vs hard-0.5 (<= %d-day threshold).\n",
    rec_max_onset_lag, ONSET_LAG_THRESHOLD_DAYS
  ))
  cat(sprintf("Default hysteresis band: upper = %.1f, lower = %.1f.\n", rec$upper, rec$lower))
} else {
  cat(sprintf(
    "RECOMMENDATION: none of the hysteresis grid options pass all three criteria\n(whipsaw reduction, full stress-window coverage, onset lag <= %d days).\n",
    ONSET_LAG_THRESHOLD_DAYS
  ))
  cat("Hard-0.5 remains the default; revisit grid spacing or thresholds.\n")
}
cat("================================================================\n")
