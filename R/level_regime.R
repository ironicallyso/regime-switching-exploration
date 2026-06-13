# Helper functions for Phase 2.2: 2-state Markov-switching AR(1) model on
# log(VIX) *levels* (mean-switching regimes), as opposed to the Phase 2
# MS-GARCH model on VIX log-returns (volatility-of-volatility regimes).
#
# No library() calls here per project convention — callers load packages.

# Fit a k=2 Markov-switching AR(p) model on y via MSwM::msmFit.
#
# sw: logical vector passed to msmFit's `sw` argument, in the order
#     (intercept, AR coefficient(s), residual variance). E.g.
#     c(TRUE, FALSE, TRUE) = switching mean, common AR(1), switching variance.
# p:  AR order (lags of y added automatically by msmFit).
fit_ms_ar <- function(y, sw, p = 1, seed) {
  set.seed(seed)
  base_model <- lm(y ~ 1, data = data.frame(y = y))
  MSwM::msmFit(base_model, k = 2, sw = sw, p = p)
}

# Identify which of the 2 states is "high-fear" by unconditional mean
# (mean = intercept / (1 - sum(AR coefs))), not by raw state index — the
# optimizer's state ordering is arbitrary.
label_high_fear_state <- function(fit) {
  coefs <- fit@Coef
  intercepts <- coefs[["(Intercept)"]]
  ar_cols <- grep("^y_", names(coefs))

  if (length(ar_cols) > 0) {
    ar_sums <- rowSums(as.matrix(coefs[, ar_cols, drop = FALSE]))
    means <- intercepts / (1 - ar_sums)
  } else {
    means <- intercepts
  }

  which.max(means)
}

# Align filtProb/smoProb (which lose the first observation to the AR lag)
# with the original date vector.
#
# With p = 1: filtProb has n-1 rows, corresponding to dates[2:n]. smoProb has
# n rows; row 1 is an initial-state row (not associated with a date), and
# rows 2:n correspond to dates[2:n] — same alignment as filtProb. Both are
# stopifnot-checked here so a future p != 1 (or an MSwM version change)
# fails loudly rather than silently misaligning dates.
extract_probs <- function(fit, dates, high_state) {
  n <- length(dates)
  filt <- fit@Fit@filtProb
  smoo <- fit@Fit@smoProb

  stopifnot(nrow(filt) == n - 1)
  stopifnot(nrow(smoo) == n)

  tibble::tibble(
    date = dates[2:n],
    filt_prob_high = filt[, high_state],
    smooth_prob_high = smoo[2:n, high_state]
  )
}

# Expected regime durations (1 / (1 - P_ii)) and the stable (stationary)
# distribution of a 2-state transition matrix P. MSwM's transMat is
# column-stochastic: P[i, j] = P(state_t = i | state_{t-1} = j), so columns
# (not rows) sum to 1 — verified empirically against fit@Fit@filtProb via
# the Hamilton-filter sanity check in the main script.
transition_stats <- function(P) {
  durations <- 1 / (1 - diag(P))

  # Stationary column-vector pi satisfies pi = P %*% pi, i.e. pi is a right
  # eigenvector of P with eigenvalue 1.
  eig <- eigen(P)
  idx <- which.min(abs(eig$values - 1))
  vec <- Re(eig$vectors[, idx])
  stable_probs <- vec / sum(vec)

  list(durations = durations, stable_probs = stable_probs)
}

# Table of contiguous high-fear episodes (rle on a logical vector), with
# start/end dates, duration in trading days, and a flag for episodes shorter
# than min_days (implausibly short for a "regime").
regime_episodes <- function(is_high, dates, min_days = 5) {
  runs <- rle(is_high)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1

  tibble::tibble(
    is_high = runs$values,
    start_date = dates[starts],
    end_date = dates[ends],
    duration_days = runs$lengths
  ) |>
    dplyr::filter(is_high) |>
    dplyr::mutate(short = duration_days < min_days) |>
    dplyr::select(-is_high)
}

# One Hamilton-filter step for a 2-state mean/AR/variance-switching model.
#
# pi_prev:  filtered state probabilities at t-1, P(S_{t-1} = i | y_1..y_{t-1})
# y_t, y_lag: y_t and y_{t-1}
# means, ar_coefs, sigmas: length-2 vectors of per-state intercepts, AR(1)
#   coefficients, and residual std devs (state order matches P's rows/cols)
# P: 2x2 column-stochastic transition matrix, P[i, j] = P(state_t = i | state_{t-1} = j)
#
# Returns the predicted (pre-observation) and filtered (post-observation)
# state probability vectors.
hamilton_filter_update <- function(pi_prev, y_t, y_lag, means, ar_coefs, sigmas, P) {
  pi_pred <- as.numeric(P %*% pi_prev)

  cond_means <- means + ar_coefs * y_lag
  dens <- dnorm(y_t, mean = cond_means, sd = sigmas)

  pi_filtered <- pi_pred * dens
  pi_filtered <- pi_filtered / sum(pi_filtered)

  list(pi_filtered = pi_filtered, pi_predicted = pi_pred)
}
