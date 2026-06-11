library(tidyverse)
library(MSGARCH)

source("R/data.R")

vix_data <- fetch_vix()

# Phase 2, Requirement 6: Fit a 2-state Markov-Switching GARCH model on
# VIX log-returns using the MSGARCH package

# Phase 2, Requirement 7: Extract smoothed state probabilities and filtered
# state probabilities for the full history

# Phase 2, Requirement 8: Plot smoothed state probabilities over time,
# overlaid with the VIX = 20 threshold from Phase 1

# Phase 2, Requirement 9: Diagnostic — compare model-assigned latent states
# to threshold-based labels from Phase 1; document alignment and divergence
