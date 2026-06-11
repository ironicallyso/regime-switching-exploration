library(tidyverse)
library(MSGARCH)

source("R/data.R")

vix_data <- fetch_vix()

# Phase 3, Requirement 10: Produce one-step-ahead filtered probability
# forecast: P(state 1 tomorrow), P(state 2 tomorrow)

# Phase 3, Requirement 11: Label states interpretively (low-vol / high-vol)
# based on which state has higher conditional variance

# Phase 3, Requirement 12: Backtest — rolling out-of-sample evaluation; for
# each day, predict next-day regime and compare to next day's threshold-based
# label; report accuracy

# Phase 3, Requirement 13: Report backtest accuracy overall and separately
# for low-vol and high-vol days (class imbalance matters)
