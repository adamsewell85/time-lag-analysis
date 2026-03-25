"""
example_analysis.py
===================
Worked example of Time Lag Analysis using UEFA Euro 2024 (Men's) data.
The four semifinalist teams are used to demonstrate the tool.

To apply this to your own data:
    1. Replace the read_excel() call with your own file and sheet name.
    2. Update the 'teams' list to match your group/team identifiers.
    3. Ensure your data contains: team, date, and numeric KPI columns.
       A 'result' column (e.g. "Win" / "Loss") is optional but enables
       win/loss colour coding in plots.
"""

import pandas as pd
import matplotlib.pyplot as plt
from time_lag_analysis import time_lag_analysis, time_lag_linear_plot, time_lag_poly_plot

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
df = pd.read_excel("../Euro2024M.xlsx", sheet_name="per90")

# -----------------------------------------------------------------------------
# 2. Filter to teams of interest
# -----------------------------------------------------------------------------
teams = ["Spain", "England", "France", "Netherlands"]
df_filtered = df[df["team"].isin(teams)].copy()

# -----------------------------------------------------------------------------
# 3. Run the analysis
#    results_included = True   if your data has a 'result' column
#    poly             = True   to also fit and compare polynomial models
# -----------------------------------------------------------------------------
results = time_lag_analysis(df_filtered, results_included=True, poly=True)

plot_df     = results["lag_data"]               # pairwise distances and lag values
summary_df  = results["profile_summary"]        # model coefficients and profiles
descriptive = results["descriptive_statistics"] # mean, SD, Lag-1 MAD, CoV

# -----------------------------------------------------------------------------
# 4. Generate plots
# -----------------------------------------------------------------------------

# Linear plot — regression line coloured by significance
fig_linear = time_lag_linear_plot(plot_df)
plt.show()

# Polynomial plot — compares linear (dashed) vs polynomial (solid) fit
# with ANOVA result and delta-R2 annotated per panel
fig_poly = time_lag_poly_plot(plot_df)
plt.show()
