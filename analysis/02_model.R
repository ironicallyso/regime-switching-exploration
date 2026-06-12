library(tidyverse)
library(MSGARCH)

source("R/data.R")
source("R/regimes.R")

VIX_THRESHOLD <- 20

vix_data <- fetch_vix()

# Phase 2, Requirement 6: Fit a 2-state Markov-Switching GARCH model on
# VIX log-returns using the MSGARCH package
spec <- CreateSpec()
fit <- FitML(spec = spec, data = vix_data$vix_log_return)

print(summary(fit))

# Phase 2, Requirement 7: Extract smoothed state probabilities and filtered
# state probabilities for the full history
state <- State(fit)

n_obs <- nrow(vix_data)

state_probs <- vix_data |>
  mutate(
    regime = label_regime(vix_close, VIX_THRESHOLD),
    filt_prob_1 = state$FiltProb[1:n_obs, 1, 1],
    filt_prob_2 = state$FiltProb[1:n_obs, 1, 2],
    smooth_prob_1 = state$SmoothProb[1:n_obs, 1, 1],
    smooth_prob_2 = state$SmoothProb[1:n_obs, 1, 2],
    viterbi_state = state$Viterbi[1:n_obs, 1]
  )

# Label latent states by inspecting conditional (unconditional) variances
# post-fit, not by state index — the higher-volatility state is "high"
state_fits <- ExtractStateFit(fit)
unc_vols <- sapply(state_fits, function(f) UncVol(f))
high_vol_state <- which.max(unc_vols)
low_vol_state <- which.min(unc_vols)

state_probs <- state_probs |>
  mutate(
    smooth_prob_high = if (high_vol_state == 1) smooth_prob_1 else smooth_prob_2,
    viterbi_regime = if_else(viterbi_state == high_vol_state, "high", "low")
  )

# Phase 2, Requirement 8: Plot smoothed state probabilities over time,
# overlaid with the VIX = 20 threshold from Phase 1
plot_data <- bind_rows(
  state_probs |>
    transmute(date, series = "Smoothed P(high-vol)", value = smooth_prob_high),
  state_probs |>
    transmute(date, series = "Threshold regime (VIX > 20)", value = as.numeric(regime == "high"))
)

print(
  ggplot(plot_data, aes(x = date, y = value)) +
    geom_line() +
    facet_wrap(~series, ncol = 1, scales = "free_y") +
    labs(title = "MS-GARCH Smoothed State Probabilities vs. VIX = 20 Threshold",
         x = NULL, y = NULL)
)

# Phase 2, Requirement 9: Diagnostic — compare model-assigned latent states
# to threshold-based labels from Phase 1; document alignment and divergence
diagnostic <- state_probs |>
  count(regime, viterbi_regime) |>
  group_by(regime) |>
  mutate(pct = 100 * n / sum(n)) |>
  ungroup()

agreement_rate <- mean(state_probs$regime == state_probs$viterbi_regime)

print(diagnostic)
print(agreement_rate)

# Diagnostic note: the MS-GARCH latent states characterize regimes in the
# VIX *return* volatility process, while the threshold label is based on the
# VIX *level*. Some divergence between `viterbi_regime` and `regime` is
# therefore expected and is not a model bug — see `diagnostic` /
# `agreement_rate` above for the observed overlap.
