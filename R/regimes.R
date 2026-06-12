label_regime <- function(vix_close, threshold = 20) {
  dplyr::if_else(vix_close > threshold, "high", "low")
}

regime_persistence <- function(regime) {
  runs <- rle(regime)

  tibble::tibble(
    regime = runs$values,
    length = runs$lengths
  ) |>
    dplyr::group_by(regime) |>
    dplyr::summarise(avg_consecutive_days = mean(length), .groups = "drop")
}

regime_transition_matrix <- function(regime) {
  n <- length(regime)

  tibble::tibble(
    from = regime[-n],
    to = regime[-1]
  ) |>
    dplyr::count(from, to, name = "n") |>
    dplyr::group_by(from) |>
    dplyr::mutate(prob = n / sum(n)) |>
    dplyr::ungroup()
}
