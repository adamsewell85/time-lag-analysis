# Time Lag Analysis

A tool for exploring inter-game variability in sport performance data.
Available in both **R** and **Python**.

| Language | Core file | Example script |
|----------|-----------|----------------|
| R | `r/TimeLagAnalysis.R` | `r/example_analysis.R` |
| Python | `python/time_lag_analysis.py` | `python/example_analysis.py` |

---

## What is Time Lag Analysis?

Time Lag Analysis quantifies how performance diverges or stabilises across increasing temporal distances between matches. For each pair of games within a team's season, it computes the Euclidean distance between standardised performance indicators and regresses that distance against the square-root of the lag (number of games between the pair). The resulting slope and significance are used to classify each team–indicator combination into one of four performance profiles.

---

## Data Format

Your data should be structured as one row per match per team, with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `team` | character | Group identifier (e.g. team name) |
| `date` | date | Match date — used to order games within each group |
| `result` | character | Match outcome, e.g. `"Win"` / `"Loss"` — **optional**, enables win/loss colour coding in plots |
| KPI columns | numeric | Any number of performance indicators (e.g. shots, passes, tackles) |

An example dataset (`Euro2024M.xlsx`, sheet `per90`) is included. It contains per-90-minute match statistics for the four UEFA Euro 2024 (Men's) semifinalist teams.

---

## Dependencies

**R**
```r
install.packages(c("dplyr", "tidyr", "purrr", "vegan", "ggplot2", "readxl"))
```

**Python**
```bash
pip install -r python/requirements.txt
# then set your working directory to python/ before running
```

---

## Quick Start

**R** *(run from the `r/` folder)*
```r
source("TimeLagAnalysis.R")
library(readxl)
library(dplyr)

df <- read_excel("your_data.xlsx", sheet = "your_sheet")

teams <- c("Team A", "Team B")
df_filtered <- df %>% filter(team %in% teams)

results <- timeLagAnalysis(df_filtered, resultsIncluded = TRUE, poly = TRUE)

plot_df     <- results$lagData
summary_df  <- results$profileSummary
descriptive <- results$descriptiveStatistics

timeLagLinearPlot(plot_df)
timeLagPolyPlot(plot_df)
```

**Python** *(run from the `python/` folder)*
```python
import pandas as pd
from time_lag_analysis import time_lag_analysis, time_lag_linear_plot, time_lag_poly_plot

df = pd.read_excel("your_data.xlsx", sheet_name="your_sheet")
df_filtered = df[df["team"].isin(["Team A", "Team B"])]

results = time_lag_analysis(df_filtered, results_included=True, poly=True)

time_lag_linear_plot(results["lag_data"])
time_lag_poly_plot(results["lag_data"])
```

See `r/example_analysis.R` / `python/example_analysis.py` for full worked examples.

---

## Interpreting Results

### Linear Profiles (`profileSummary`)

Each team–indicator combination is assigned one of four profiles based on the linear regression slope and p-value:

| Profile | Condition | Interpretation |
|---------|-----------|----------------|
| Stochastic | p > 0.05 | No significant relationship between lag and performance distance |
| Convergence | slope < -0.02 | Performance becomes more similar as lag increases |
| Stable | -0.02 ≤ slope ≤ 0.02 | Consistent level of variability regardless of lag |
| Directional Change | slope > 0.02 | Performance diverges as lag increases |

### Descriptive Statistics (`descriptiveStatistics`)

| Column | Description |
|--------|-------------|
| `Mean` | Mean raw value across all matches |
| `SD` | Standard deviation |
| `Lag1_MAD` | Mean absolute difference between consecutive z-scored values |
| `CoV` | Coefficient of variation (%) |

---

## Citation

If you use this tool in your research, please cite:

> [Author(s), Year. Title. *Journal*, Volume(Issue), pp. DOI]

*(Citation will be updated upon publication.)*
