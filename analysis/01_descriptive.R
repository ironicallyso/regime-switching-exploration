library(tidyverse)
library(quantmod)
library(lubridate)

source("R/data.R")
source("R/regimes.R")

VIX_THRESHOLD <- 20

vix_data <- fetch_vix()

# Phase 1, Requirement 2: Plot VIX close prices over full history with a
# horizontal reference line at VIX = 20
print(
  ggplot(vix_data, aes(x = date, y = vix_close)) +
    geom_line() +
    geom_hline(yintercept = VIX_THRESHOLD, linetype = "dashed", color = "red") +
    labs(title = "VIX Close Price", x = NULL, y = "VIX Close")
)

vix_labeled <- vix_data |>
  mutate(
    regime = label_regime(vix_close, VIX_THRESHOLD),
    year = year(date)
  )

# Phase 1, Requirement 3: Compute % of trading days above and below 20,
# for full history and broken out by calendar year
overall_pct <- vix_labeled |>
  count(regime) |>
  mutate(pct = 100 * n / sum(n))

by_year_pct <- vix_labeled |>
  count(year, regime) |>
  group_by(year) |>
  mutate(pct = 100 * n / sum(n)) |>
  ungroup()

print(overall_pct)
print(by_year_pct, n = Inf)

# Phase 1, Requirement 4: Compute regime persistence — average number of
# consecutive days spent in each threshold-based regime (above/below 20)
persistence <- regime_persistence(vix_labeled$regime)

print(persistence)

# Phase 1, Requirement 5: Compute empirical transition matrix from
# threshold-based labels: P(low->low), P(low->high), P(high->low), P(high->high)
transition <- regime_transition_matrix(vix_labeled$regime)

print(transition)
