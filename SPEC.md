# Regime Switching Exploration — SPEC.md

> Last updated: 2026-06-11

## Goal

Exploratory R project to fit a 2-state Markov-Switching GARCH model on VIX log-returns,
characterize volatility regimes, and produce next-day regime probability forecasts.
Standalone research repo; findings will inform whether and how to integrate regime signals
into a separate options modeler.

---

## Requirements

### Phase 1: Descriptive Analysis
1. Load VIX close prices from `data/vix.csv`; if file does not exist, download via `quantmod` and write it
2. Plot VIX close prices over full history with a horizontal reference line at VIX = 20
3. Compute % of trading days above and below 20, for full history and broken out by calendar year
4. Compute regime persistence: average number of consecutive days spent in each threshold-based regime (above/below 20)
5. Compute empirical transition matrix from threshold-based labels: P(low→low), P(low→high), P(high→low), P(high→high)

### Phase 2: Fit 2-State MS-GARCH Model
6. Fit a 2-state Markov-Switching GARCH model on VIX log-returns using the `MSGARCH` package
7. Extract smoothed state probabilities and filtered state probabilities for the full history
8. Plot smoothed state probabilities over time, overlaid with the VIX = 20 threshold from Phase 1
9. Diagnostic: compare model-assigned latent states to threshold-based labels from Phase 1 — document alignment and divergence; do not expect them to be identical

### Phase 3: Prediction Interface
10. Produce one-step-ahead filtered probability forecast: P(state 1 tomorrow), P(state 2 tomorrow)
11. Label states interpretively (low-vol / high-vol) based on which state has higher conditional variance
12. Backtest: rolling out-of-sample evaluation — for each day, predict next-day regime and compare to next day's threshold-based label; report accuracy
13. Report backtest accuracy overall and separately for low-vol and high-vol days (class imbalance matters)

---

## Out of Scope
- 3-state models
- Integration into the options modeler (this repo is standalone; integration is a future decision)
- Intraday data or anything other than daily VIX close
- Automated trading signals

---

## Acceptance Criteria
- Phase 1: plots render without error; transition matrix rows sum to 1.0; yearly breakdown covers full VIX history
- Phase 2: model converges; smoothed probability plot produced; diagnostic alignment documented in a comment or markdown cell
- Phase 3: one-step-ahead forecast runs on the most recent observation and returns two probabilities summing to 1.0; backtest produces a confusion matrix and per-class accuracy

---

## Data
- Source: `data/vix.csv` (gitignored) — VIX daily close, columns: `date`, `vix_close`
- If `data/vix.csv` is absent, `R/data.R` downloads full history via `quantmod` (`^VIX`) and writes the CSV
- No re-download if file exists (incremental update only needed if stale)

---

## Known Pitfalls
- Do not model VIX levels — model VIX log-returns (`log(vix_close / lag(vix_close))`); levels are non-stationary
- Do not expect MS-GARCH latent states to map 1:1 to the VIX > 20 threshold — the model finds volatility regimes in the return process, not price level regimes; document the divergence rather than treating it as a bug
- State labeling (low-vol / high-vol) must be done post-fit by inspecting conditional variances, not assumed from state index order
- `MSGARCH` uses its own simulation/forecasting interface — do not use `forecast::forecast()` on MS-GARCH objects
