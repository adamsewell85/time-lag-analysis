# =============================================================================
# generate_basketball_figures.R
# =============================================================================
# Reproduces the two basketball figures for the 2023 FIBA Basketball World Cup:
#   Figure 2 - Efficiency (OFFRTG, DEFRTG)
#   Figure 3 - Technical  (eFG, TOV%, OREB%, FTR)
# from the corrected data (FIBA_WC_2023.xlsx), with permutation p-values shown
# per panel. Run from the r/ folder.
#
# Output directory can be overridden with the FIG_OUT_DIR environment variable;
# defaults to ./figures.
# =============================================================================

source("TimeLagAnalysis.R")
suppressMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
})

out_dir <- Sys.getenv("FIG_OUT_DIR", unset = "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load corrected data; label USA and fix team display order
# -----------------------------------------------------------------------------
team_order <- c("Germany", "Serbia", "Canada", "USA")

df <- read_excel("../FIBA_WC_2023.xlsx") %>%
  mutate(team = ifelse(team == "United States", "USA", team))

# -----------------------------------------------------------------------------
# 2. Run time-lag analysis (uniform grey points -> resultsIncluded = FALSE)
# -----------------------------------------------------------------------------
results <- timeLagAnalysis(df, resultsIncluded = FALSE, poly = FALSE)
lag_data <- results$lagData %>%
  mutate(team = factor(team, levels = team_order))

# -----------------------------------------------------------------------------
# 3. Helper: build one figure for a set of KPIs with display labels
# -----------------------------------------------------------------------------
# Physical size is kept small (with high dpi) so the theme's fixed point-size
# fonts render large relative to the panels, matching the published figures.
# pixels = width_in * dpi ; height 9 in * 600 dpi = 5400 px (as in the originals).
make_figure <- function(kpi_keys, kpi_labels, file, width_in) {
  plot_df <- lag_data %>%
    filter(performanceIndicator %in% kpi_keys) %>%
    mutate(performanceIndicator = factor(performanceIndicator,
                                         levels = kpi_keys,
                                         labels = kpi_labels))
  p <- timeLagLinearPlot(plot_df)
  ggsave(
    filename = file.path(out_dir, file),
    plot     = p,
    width    = width_in,
    height   = 9,
    units    = "in",
    dpi      = 600,
    device   = "jpeg"
  )
  message("Wrote ", file.path(out_dir, file))
}

# -----------------------------------------------------------------------------
# 4. Figure 2 (Efficiency) and Figure 3 (Technical)
# -----------------------------------------------------------------------------
make_figure(
  kpi_keys   = c("offrtg", "defrtg"),
  kpi_labels = c("OFFRTG", "DEFRTG"),
  file       = "Figure 2 - Efficiency (OFFRTG DEFRTG).jpg",
  width_in   = 5.5
)

make_figure(
  kpi_keys   = c("efg_pct", "tov_pct", "oreb_pct", "ftr"),
  kpi_labels = c("eFG", "TOV%", "OREB%", "FTR"),
  file       = "Figure 3 - Technical (eFG TOV OREB FTR).jpg",
  width_in   = 10
)
