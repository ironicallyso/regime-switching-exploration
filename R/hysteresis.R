# Helper functions for Phase 2.3: hysteresis (Schmitt-trigger) reliability
# pass on the Phase 2.2 high-fear regime signal (filt_prob_high). Adds a
# dead-band between entry and exit thresholds to reduce whipsaw from a
# single hard 0.5 cutoff, and quantifies the stability/coverage/lag
# trade-offs versus the hard-0.5 baseline.
#
# No library() calls here per project convention — callers load packages.

# Schmitt-trigger (two-threshold hysteresis) classifier on a probability
# series.
#
# prob:  numeric vector of P(high-fear) in [0, 1], e.g. filt_prob_high
# upper: enter high-fear (TRUE) when prob crosses ABOVE this threshold
# lower: exit high-fear (FALSE) when prob crosses BELOW this threshold
#        (lower <= upper)
#
# The initial state is prob[1] > 0.5 (the midpoint), independent of
# upper/lower, so the first observation is classified on its own merits
# rather than requiring a crossing from an arbitrary bootstrap state.
#
# The degenerate case upper == lower == 0.5 (the hard-0.5 baseline) is NOT
# handled by this function — callers should classify that case directly as
# `prob > 0.5` to avoid ambiguity when prob is exactly 0.5 (see
# analysis/02d_hysteresis.R).
#
# Returns a logical vector, same length as prob: TRUE = high-fear.
schmitt_trigger <- function(prob, upper, lower) {
  stopifnot(lower <= upper)

  n <- length(prob)
  state <- logical(n)
  state[1] <- prob[1] > 0.5

  for (i in seq_len(n)[-1]) {
    if (state[i - 1]) {
      state[i] <- !(prob[i] < lower)
    } else {
      state[i] <- prob[i] > upper
    }
  }

  state
}

# Stability metrics for a high-fear classification, built on
# regime_episodes() (R/level_regime.R).
#
# is_high:  logical vector, TRUE = high-fear
# dates:    Date vector, same length as is_high
# min_days: passed through to regime_episodes() for the `short` flag and
#           pct_episodes_short
#
# Returns a one-row tibble:
#   n_episodes, episodes_per_year, median_duration_days, min_duration_days,
#   max_duration_days, pct_episodes_short, pct_time_high_fear
regime_stability_metrics <- function(is_high, dates, min_days = 5) {
  episodes <- regime_episodes(is_high, dates, min_days = min_days)
  span_years <- as.numeric(max(dates) - min(dates)) / 365.25

  if (nrow(episodes) == 0) {
    return(tibble::tibble(
      n_episodes = 0L,
      episodes_per_year = 0,
      median_duration_days = NA_real_,
      min_duration_days = NA_real_,
      max_duration_days = NA_real_,
      pct_episodes_short = NA_real_,
      pct_time_high_fear = 100 * mean(is_high)
    ))
  }

  tibble::tibble(
    n_episodes = nrow(episodes),
    episodes_per_year = nrow(episodes) / span_years,
    median_duration_days = median(episodes$duration_days),
    min_duration_days = min(episodes$duration_days),
    max_duration_days = max(episodes$duration_days),
    pct_episodes_short = 100 * mean(episodes$duration_days < min_days),
    pct_time_high_fear = 100 * mean(is_high)
  )
}

# Find the ~window_days-long stretch with the most state transitions
# (flips) in a logical vector — used to pick a "whippy" zoom window for
# plotting. Transitions are counted via diff(as.integer(is_high)) != 0.
#
# is_high:     logical vector (typically the hard-0.5 classification)
# dates:       Date vector, same length as is_high
# window_days: approximate window width in calendar days (default ~1
#              trading year)
# pad_days:    calendar-day padding added to the start/end of the chosen
#              window, for plotting context
#
# Returns list(start = Date, end = Date).
find_whipsaw_window <- function(is_high, dates, window_days = 252, pad_days = 10) {
  flips <- c(0, diff(as.integer(is_high)) != 0)
  n <- length(dates)

  best_start_idx <- 1L
  best_end_idx <- n
  best_count <- -1L

  for (i in seq_len(n)) {
    end_date <- dates[i] + window_days
    if (end_date > dates[n]) end_date <- dates[n]
    j <- max(which(dates <= end_date))
    if (j <= i) next

    cnt <- sum(flips[i:j])
    if (cnt > best_count) {
      best_count <- cnt
      best_start_idx <- i
      best_end_idx <- j
    }
  }

  list(
    start = dates[best_start_idx] - pad_days,
    end = dates[best_end_idx] + pad_days
  )
}
