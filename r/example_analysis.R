# =============================================================================
# example_analysis.R
# =============================================================================
# Worked example of Time Lag Analysis using UEFA Euro 2024 (Men's) data.
# The four semifinalist teams are used to demonstrate the tool.
#
# To apply this to your own data:
#   1. Replace the read_excel() call with your own file and sheet name.
#   2. Update the 'teams' vector to match your group/team identifiers.
#   3. Ensure your data contains: team, date, and numeric KPI columns.
#      A 'result' column (e.g. "Win" / "Loss") is optional but enables
#      win/loss colour coding in plots.
# =============================================================================

source("TimeLagAnalysis.R")
library(readxl)
library(dplyr)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
df <- read_excel("../Euro2024M.xlsx", sheet = "per90")

# -----------------------------------------------------------------------------
# 2. Filter to teams of interest
# -----------------------------------------------------------------------------
teams <- c("Spain", "England", "France", "Netherlands")
df_filtered <- df %>% filter(team %in% teams)

# -----------------------------------------------------------------------------
# 3. Run the analysis
#    resultsIncluded = TRUE   if your data has a 'result' column
#    poly            = TRUE   to also fit and compare polynomial models
# -----------------------------------------------------------------------------
results <- timeLagAnalysis(df_filtered, resultsIncluded = TRUE, poly = TRUE)

plot_df     <- results$lagData              # pairwise distances and lag values
summary_df  <- results$profileSummary       # model coefficients and profiles
descriptive <- results$descriptiveStatistics # mean, SD, Lag-1 MAD, CoV

# -----------------------------------------------------------------------------
# 4. Generate plots
# -----------------------------------------------------------------------------

# Linear plot — regression line coloured by significance
timeLagLinearPlot(plot_df)

# Polynomial plot — compares linear (dashed) vs polynomial (solid) fit
# with ANOVA result and delta-R2 annotated per panel
timeLagPolyPlot(plot_df)
