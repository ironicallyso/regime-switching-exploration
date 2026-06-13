# Helper functions for Phase 2.1: expanding-window OOS comparison of
# MS-GARCH vs single-regime baselines on one-step-ahead variance forecasts.
#
# No library() calls here per project convention — callers load packages.

# Build the three specs to compare. Single-regime models are expressed as
# K = 1 MSGARCH specs so estimation/forecasting stays apples-to-apples with
# the 2-regime model. Add a 4th entry here (e.g. a 2-regime Student-t spec)
# to extend the comparison.
build_specs <- function() {
  list(
    sgarch_norm = CreateSpec(
      variance.spec = list(model = "sGARCH"),
      distribution.spec = list(distribution = "norm"),
      switch.spec = list(K = 1)
    ),
    sgarch_std = CreateSpec(
      variance.spec = list(model = "sGARCH"),
      distribution.spec = list(distribution = "std"),
      switch.spec = list(K = 1)
    ),
    msgarch_2s_norm = CreateSpec(
      variance.spec = list(model = "sGARCH"),
      distribution.spec = list(distribution = "norm"),
      switch.spec = list(K = 2)
    )
  )
}

# Expanding-window, one-step-ahead variance forecast backtest.
#
# y:           numeric vector of log-returns (full series)
# specs:       named list of MSGARCH_SPEC objects (from build_specs())
# train_frac:  fraction of the series used for the initial training window
# refit_every: refit cadence in trading days; parameters held fixed between
#              refits while the variance filter is rolled forward via
#              predict(fit, newdata = ..., nahead = 1)
# seed:        passed to set.seed() before each FitML() call, per the
#              reproducibility request (empirically FitML() on this data is
#              deterministic given the data, but the seed is set anyway so
#              behavior is reproducible if that ever changes)
#
# Returns a tibble with one row per OOS day: `proxy` (squared demeaned
# return for that day) and one variance-forecast column per model name.
run_expanding_window <- function(y, specs, train_frac, refit_every, seed) {
  n <- length(y)
  n_train <- floor(train_frac * n)
  oos_idx <- (n_train + 1):n
  n_oos <- length(oos_idx)

  y_mean <- mean(y)

  model_names <- names(specs)
  forecasts <- matrix(NA_real_, nrow = n_oos, ncol = length(specs))
  colnames(forecasts) <- model_names

  fits <- vector("list", length(specs))
  names(fits) <- model_names

  for (j in seq_len(n_oos)) {
    i <- oos_idx[j]
    train_data <- y[1:(i - 1)]

    needs_refit <- (j - 1) %% refit_every == 0

    for (m in model_names) {
      if (needs_refit) {
        set.seed(seed)
        fits[[m]] <- FitML(specs[[m]], train_data)
      }
      pred <- predict(fits[[m]], newdata = train_data, nahead = 1)
      forecasts[j, m] <- pred$vol^2
    }

    if (j %% 100 == 0 || j == n_oos) {
      message(sprintf("OOS day %d / %d (%.1f%%)", j, n_oos, 100 * j / n_oos))
    }
  }

  tibble::as_tibble(forecasts) |>
    dplyr::mutate(
      proxy = (y[oos_idx] - y_mean)^2,
      .before = 1
    )
}

# Long-format losses: one row per OOS day x model, with QLIKE and MSE.
compute_losses <- function(forecast_tbl, model_names) {
  forecast_tbl |>
    dplyr::select(dplyr::all_of(c("date", "proxy", model_names))) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(model_names),
      names_to = "model",
      values_to = "forecast"
    ) |>
    dplyr::mutate(
      qlike = log(forecast) + proxy / forecast,
      mse = (forecast - proxy)^2
    )
}

# Wide summary table: model x mean QLIKE, mean MSE.
summarize_results <- function(loss_tbl) {
  loss_tbl |>
    dplyr::group_by(model) |>
    dplyr::summarise(
      mean_qlike = mean(qlike),
      mean_mse = mean(mse),
      .groups = "drop"
    )
}

# Diebold-Mariano test on a loss differential d = loss_target - loss_baseline.
# Negative mean_diff with a small p-value means `target` has lower average
# loss than `baseline` (target wins).
#
# Uses a Newey-West/Bartlett HAC variance estimate of the mean of d, since
# d can be autocorrelated (e.g. due to the monthly refit cadence). No extra
# package dependency (sandwich/multDM not available in this project).
dm_test <- function(loss_target, loss_baseline) {
  d <- loss_target - loss_baseline
  d <- d[!is.na(d)]
  n <- length(d)
  d_bar <- mean(d)

  max_lag <- max(1, floor(4 * (n / 100)^(2 / 9)))

  gamma0 <- stats::var(d)
  v <- gamma0
  for (k in seq_len(max_lag)) {
    cov_k <- stats::cov(d[1:(n - k)], d[(1 + k):n])
    weight <- 1 - k / (max_lag + 1)
    v <- v + 2 * weight * cov_k
  }

  dm_stat <- d_bar / sqrt(v / n)
  p_value <- 2 * stats::pt(-abs(dm_stat), df = n - 1)

  list(statistic = dm_stat, p.value = p_value, mean_diff = d_bar)
}

# Return dir/basename.ext, or dir/basename_v2.ext, _v3, ... if the file
# already exists. Never overwrites an existing output.
versioned_path <- function(dir, basename, ext) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  candidate <- file.path(dir, paste0(basename, ".", ext))
  if (!file.exists(candidate)) {
    return(candidate)
  }

  v <- 2
  repeat {
    candidate <- file.path(dir, paste0(basename, "_v", v, ".", ext))
    if (!file.exists(candidate)) {
      return(candidate)
    }
    v <- v + 1
  }
}
