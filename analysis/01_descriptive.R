library(tidyverse)
library(quantmod)
library(lubridate)

source("R/data.R")

vix_data <- fetch_vix()

# Phase 1, Requirement 2: Plot VIX close prices over full history with a
# horizontal reference line at VIX = 20

# Phase 1, Requirement 3: Compute % of trading days above and below 20,
# for full history and broken out by calendar year

# Phase 1, Requirement 4: Compute regime persistence — average number of
# consecutive days spent in each threshold-based regime (above/below 20)

# Phase 1, Requirement 5: Compute empirical transition matrix from
# threshold-based labels: P(low->low), P(low->high), P(high->low), P(high->high)
