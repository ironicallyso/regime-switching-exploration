# Regime Switching Exploration

## Project Overview
Exploratory R project fitting a 2-state Markov-Switching GARCH model on VIX log-returns.
Goal: characterize volatility regimes and produce next-day regime probability forecasts.
Standalone research repo; not connected to the options modeler.

## Project Spec
At the start of each session, read `SPEC.md` to understand current requirements and phase
status before planning or writing code.

## Tech Stack
- R (via RStudio)
- Key packages: `tidyverse`, `quantmod`, `MSGARCH`, `lubridate`
- Data: `data/vix.csv` (gitignored daily VIX close); auto-downloaded if absent

## Project Structure
```
R/
  data.R          # fetch_vix(): download VIX via quantmod if data/vix.csv absent; return tibble
analysis/
  01_descriptive.R   # Phase 1: plots, threshold stats, transition matrix
  02_model.R         # Phase 2: fit MS-GARCH, extract probabilities, diagnostic
  03_predict.R       # Phase 3: one-step-ahead forecast + backtest
data/               # gitignored; holds vix.csv
```

## Dev Workflow
```r
# Install deps (first time)
install.packages(c("tidyverse", "quantmod", "MSGARCH", "lubridate"))

# Run phases in order
source("analysis/01_descriptive.R")
source("analysis/02_model.R")
source("analysis/03_predict.R")
```

## Key Rules
- Always model VIX log-returns, never raw levels — levels are non-stationary
- `data/vix.csv` is gitignored — never commit it; `R/data.R` handles download on first run
- MS-GARCH latent states ≠ VIX > 20 threshold labels — document divergence, don't fix it
- Label states (low-vol / high-vol) by inspecting conditional variances post-fit, not by state index
- Use `MSGARCH` forecasting interface only — do not call `forecast::forecast()` on MS-GARCH objects
- No `library()` calls inside `R/` functions — callers (analysis scripts) load packages
