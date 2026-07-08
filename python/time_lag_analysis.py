"""
time_lag_analysis.py
====================
Time Lag Analysis: A tool for exploring inter-game variability in sport.

Quantifies how performance diverges or stabilises across increasing temporal
distances between matches, using pairwise Euclidean distance on standardised
KPIs regressed against the square-root of the lag period.

REQUIRED INPUT COLUMNS
    team   : str  — group identifier (e.g. team name)
    date   : date — match date, used to order games within each group
    <KPIs> : float — one or more numeric performance indicator columns
    result : str  (optional) — match outcome; enables win/loss colouring

OUTPUTS  (returned as a dict)
    lag_data              : pairwise distances, lag values, and model labels
    profile_summary       : regression coefficients, R2, p-values, and profiles
    descriptive_statistics: mean, SD, Lag-1 MAD, and CoV per team x KPI

LINEAR PROFILES
    Stochastic         : p > 0.05  — no significant relationship with lag
    Convergence        : slope < -0.02 — performance stabilising over time
    Stable             : -0.02 <= slope <= 0.02 — consistent performance
    Directional Change : slope > 0.02  — performance diverging over time

SIGNIFICANCE
    The pairwise distances are not mutually independent (each game contributes
    to multiple comparisons), so parametric p-values are inflated. The linear
    slope significance that drives profile classification is therefore assessed
    by a permutation test: game order is reshuffled n_perm times to build an
    empirical null distribution of the slope, and p is the proportion of
    |permuted slopes| >= |observed slope|. The linear-vs-polynomial comparison
    keeps the parametric ANOVA F-test with visual inspection (a plain order
    shuffle is not a valid null for an isolated nested term).

DEPENDENCIES
    pip install pandas numpy scipy matplotlib
"""

import numpy as np
import pandas as pd
from itertools import combinations
from scipy.stats import f as f_dist
from scipy.spatial.distance import euclidean
from scipy import stats


# =============================================================================
# Core analysis function
# =============================================================================

def time_lag_analysis(df: pd.DataFrame, results_included: bool = False, poly: bool = False,
                      n_perm: int = 5000, seed: int = 42) -> dict:
    """
    Run time lag analysis on a tidy match-level dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain 'team', 'date', and at least one numeric KPI column.
        An optional 'result' column enables win/loss colouring.
    results_included : bool
        Set True if df contains a 'result' column.
    poly : bool
        Set True to also fit and compare a degree-2 polynomial model.
    n_perm : int
        Number of permutations for the linear-slope significance test.
    seed : int
        Random seed for the permutation test (reproducibility).

    Returns
    -------
    dict with keys: 'lag_data', 'profile_summary', 'descriptive_statistics'
    """
    non_kpi = {"team", "date", "result", "game"}
    kpi_cols = [c for c in df.select_dtypes(include="number").columns if c not in non_kpi]

    # Sort and number games per team
    df = df.copy()
    df["date"] = pd.to_datetime(df["date"])
    df = df.sort_values(["team", "date"]).reset_index(drop=True)
    df["game"] = df.groupby("team").cumcount() + 1

    # Z-score standardise KPIs (per whole dataset, consistent with R version)
    df_scaled = df.copy()
    for col in kpi_cols:
        mu = df[col].mean()
        sd = df[col].std(ddof=1)
        df_scaled[col] = (df[col] - mu) / sd if sd > 0 else 0.0

    all_lag_rows = []
    all_summary_rows = []

    for team, team_df in df_scaled.groupby("team"):
        team_df = team_df.reset_index(drop=True)
        n = len(team_df)
        if n < 2:
            continue

        for kpi in kpi_cols:
            values = team_df[kpi].values
            pairs = list(combinations(range(n), 2))

            rows = []
            for i, j in pairs:
                dist = abs(values[i] - values[j])   # Euclidean on single column
                lag = abs(i - j)
                sqrt_lag = np.sqrt(lag)

                row = {
                    "team": team,
                    "performance_indicator": kpi,
                    "game_i": i + 1,
                    "game_j": j + 1,
                    "lag": lag,
                    "sqrt_lag": sqrt_lag,
                    "distance": dist,
                }

                if results_included and "result" in team_df.columns:
                    row["win_pair"] = (
                        f"{team_df['result'].iloc[i]} and {team_df['result'].iloc[j]}"
                    )

                rows.append(row)

            kpi_df = pd.DataFrame(rows)

            # --- Linear model ---
            x = kpi_df["sqrt_lag"].values
            y = kpi_df["distance"].values
            slope, intercept, r, _p_parametric, _ = stats.linregress(x, y)
            r2_linear = r ** 2

            # Permutation test of the slope (game-order shuffle) -> empirical p
            p_linear = _perm_slope_pvalue(values, slope, n_perm, seed)

            linear_profile = _classify_profile(slope, p_linear)
            linear_formula = _format_formula(slope, intercept)
            p_linear_fmt = _format_p(p_linear)
            r2_linear_fmt = f"{r2_linear:.2f}"

            # Significance label per row
            kpi_df["linear_p_value"] = p_linear_fmt
            kpi_df["linear_profile"] = linear_profile
            kpi_df["linear_formula"] = linear_formula
            kpi_df["r2"] = r2_linear_fmt
            kpi_df["significance"] = "Significant" if p_linear <= 0.05 else "Not Significant"

            summary_row = {
                "team": team,
                "performance_indicator": kpi,
                "r2": r2_linear_fmt,
                "linear_p_value": p_linear_fmt,
                "linear_profile": linear_profile,
                "linear_formula": linear_formula,
            }

            # --- Optional polynomial model ---
            if poly:
                poly_coeffs = np.polyfit(x, y, 2)   # [a, b, c] for ax^2 + bx + c
                a, b, c = poly_coeffs
                y_poly_pred = np.polyval(poly_coeffs, x)
                ss_res_poly = np.sum((y - y_poly_pred) ** 2)
                ss_tot = np.sum((y - y.mean()) ** 2)
                r2_poly = 1 - ss_res_poly / ss_tot if ss_tot > 0 else np.nan

                poly_formula = _format_poly_formula(a, b, c)
                delta_r2 = r2_poly - r2_linear

                # ANOVA F-test: linear vs polynomial
                anova_report, model_recommendation = _anova_comparison(
                    x, y, slope, intercept, poly_coeffs, r2_linear, r2_poly, delta_r2
                )

                kpi_df["poly_r2"] = f"{r2_poly:.2f}" if not np.isnan(r2_poly) else None
                kpi_df["poly_formula"] = poly_formula
                kpi_df["anova_report"] = anova_report
                kpi_df["delta_r2"] = f"{delta_r2:.2f}" if not np.isnan(delta_r2) else None
                kpi_df["model_recommendation"] = model_recommendation

                summary_row.update({
                    "poly_r2": f"{r2_poly:.2f}" if not np.isnan(r2_poly) else None,
                    "poly_formula": poly_formula,
                    "anova_report": anova_report,
                    "delta_r2": f"{delta_r2:.2f}" if not np.isnan(delta_r2) else None,
                    "model_recommendation": model_recommendation,
                })

            all_lag_rows.append(kpi_df)
            all_summary_rows.append(summary_row)

    lag_data = pd.concat(all_lag_rows, ignore_index=True) if all_lag_rows else pd.DataFrame()
    profile_summary = pd.DataFrame(all_summary_rows)

    # --- Descriptive statistics (on original unstandardised values) ---
    descriptive_rows = []
    for team, team_df in df.groupby("team"):
        for kpi in kpi_cols:
            vals = team_df[kpi].dropna().values
            if len(vals) == 0:
                continue
            mean_val = vals.mean()
            sd_val = vals.std(ddof=1)
            cov = (sd_val / mean_val * 100) if mean_val != 0 else np.nan
            z = (vals - mean_val) / sd_val if sd_val > 0 else np.zeros_like(vals)
            lag1_mad = np.mean(np.abs(np.diff(z))) if len(z) > 1 else np.nan
            descriptive_rows.append({
                "team": team,
                "performance_indicator": kpi,
                "mean": round(mean_val, 1),
                "sd": round(sd_val, 2),
                "lag1_mad": round(lag1_mad, 3),
                "cov": round(cov, 1),
            })

    descriptive_statistics = pd.DataFrame(descriptive_rows)

    return {
        "lag_data": lag_data,
        "profile_summary": profile_summary,
        "descriptive_statistics": descriptive_statistics,
    }


# =============================================================================
# Plotting functions
# =============================================================================

def time_lag_linear_plot(df: pd.DataFrame):
    """
    Faceted scatter plot with linear regression lines.
    Regression lines are coloured by significance (black = significant,
    grey = not significant). Points are optionally coloured by win_pair.

    Parameters
    ----------
    df : pd.DataFrame — the 'lag_data' output from time_lag_analysis()

    Returns
    -------
    matplotlib Figure
    """
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    teams = df["team"].unique()
    indicators = df["performance_indicator"].unique()
    has_winpair = "win_pair" in df.columns

    fig, axes = plt.subplots(
        nrows=len(teams),
        ncols=len(indicators),
        figsize=(3.5 * len(indicators), 3 * len(teams)),
        squeeze=False,
    )

    win_colours = {
        "Win and Win": "#2C7BB6",
        "Loss and Win": "#2C7BB6",
        "Win and Loss": "#D7191C",
        "Loss and Loss": "#D7191C",
    }

    for row_idx, team in enumerate(teams):
        for col_idx, kpi in enumerate(indicators):
            ax = axes[row_idx][col_idx]
            subset = df[(df["team"] == team) & (df["performance_indicator"] == kpi)]

            if subset.empty:
                ax.set_visible(False)
                continue

            x = subset["sqrt_lag"].values
            y = subset["distance"].values

            # Scatter
            if has_winpair:
                for wp, colour in win_colours.items():
                    mask = subset["win_pair"] == wp
                    ax.scatter(x[mask.values], y[mask.values], color=colour, s=10, alpha=0.7)
            else:
                ax.scatter(x, y, color="darkgrey", s=10, alpha=0.7)

            # Regression line
            significant = subset["significance"].iloc[0] == "Significant"
            line_colour = "black" if significant else "grey"
            if len(x) >= 2:
                m, b, *_ = stats.linregress(x, y)
                x_line = np.linspace(x.min(), x.max(), 100)
                ax.plot(x_line, m * x_line + b, color=line_colour, linewidth=1.2)

            ax.set_ylim(0, 6)
            ax.set_xlim(left=0)
            ax.tick_params(labelsize=8)
            ax.spines[["top", "right"]].set_color("grey")
            ax.spines[["bottom", "left"]].set_color("grey")

            if row_idx == 0:
                ax.set_title(kpi, fontsize=9, pad=4)
            if col_idx == len(indicators) - 1:
                ax.set_ylabel(team, fontsize=9, rotation=270, labelpad=14)
                ax.yaxis.set_label_position("right")
            if row_idx == len(teams) - 1:
                ax.set_xlabel("Time Lag (sqrt)", fontsize=9)

    fig.tight_layout()
    return fig


def time_lag_poly_plot(df: pd.DataFrame):
    """
    Faceted plot comparing linear (dashed) and polynomial (solid) regression.
    ANOVA result and delta-R2 are annotated in each panel.

    Parameters
    ----------
    df : pd.DataFrame — the 'lag_data' output from time_lag_analysis()
                        with poly=True (must contain 'anova_report', 'delta_r2')

    Returns
    -------
    matplotlib Figure
    """
    import matplotlib.pyplot as plt

    required = {"anova_report", "delta_r2"}
    if not required.issubset(df.columns):
        raise ValueError("poly plot requires 'anova_report' and 'delta_r2' columns. "
                         "Run time_lag_analysis() with poly=True.")

    teams = df["team"].unique()
    indicators = df["performance_indicator"].unique()

    fig, axes = plt.subplots(
        nrows=len(teams),
        ncols=len(indicators),
        figsize=(3.5 * len(indicators), 3 * len(teams)),
        squeeze=False,
    )

    for row_idx, team in enumerate(teams):
        for col_idx, kpi in enumerate(indicators):
            ax = axes[row_idx][col_idx]
            subset = df[(df["team"] == team) & (df["performance_indicator"] == kpi)]

            if subset.empty:
                ax.set_visible(False)
                continue

            x = subset["sqrt_lag"].values
            y = subset["distance"].values

            ax.scatter(x, y, color="grey", s=10, alpha=0.7)

            if len(x) >= 2:
                x_line = np.linspace(x.min(), x.max(), 100)
                # Linear (dashed)
                m, b, *_ = stats.linregress(x, y)
                ax.plot(x_line, m * x_line + b, color="grey", linewidth=1, linestyle="--")
                # Polynomial (solid)
                coeffs = np.polyfit(x, y, 2)
                ax.plot(x_line, np.polyval(coeffs, x_line), color="grey", linewidth=1, linestyle="-")

            # Annotation
            anova_txt = subset["anova_report"].iloc[0]
            delta_txt = subset["delta_r2"].iloc[0]
            if anova_txt and delta_txt:
                ax.text(
                    0.98, 0.98,
                    f"{anova_txt}\nΔR² = {delta_txt}",
                    transform=ax.transAxes,
                    fontsize=7,
                    ha="right", va="top",
                )

            ax.set_ylim(0, 6)
            ax.set_xlim(left=0)
            ax.tick_params(labelsize=8)
            ax.spines[["top", "right"]].set_color("grey")
            ax.spines[["bottom", "left"]].set_color("grey")

            if row_idx == 0:
                ax.set_title(kpi, fontsize=9, pad=4)
            if col_idx == len(indicators) - 1:
                ax.set_ylabel(team, fontsize=9, rotation=270, labelpad=14)
                ax.yaxis.set_label_position("right")
            if row_idx == len(teams) - 1:
                ax.set_xlabel("Time Lag (sqrt)", fontsize=9)

    fig.tight_layout()
    return fig


# =============================================================================
# Internal helpers
# =============================================================================

def _perm_slope_pvalue(values: np.ndarray, obs_slope: float,
                       n_perm: int = 5000, seed: int = 42) -> float:
    """
    Empirical p-value for the time-lag regression slope via a game-order
    permutation. The square-root lag values are fixed across permutations, so
    only the pairwise distances are recomputed; the slope is obtained in closed
    form (slope = sum(w * distance)) for speed.
    """
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]
    n = len(values)
    if n < 3:
        return np.nan

    idx = np.array(list(combinations(range(n), 2)))
    i, j = idx[:, 0], idx[:, 1]
    lagv = np.sqrt(np.abs(i - j))
    lc = lagv - lagv.mean()
    w = lc / np.sum(lc ** 2)              # slope = sum(w * distance)

    rng = np.random.default_rng(seed)
    abs_obs = abs(obs_slope)
    count = 0
    for _ in range(n_perm):
        vp = rng.permutation(values)
        dp = np.abs(vp[i] - vp[j])
        if abs(np.sum(w * dp)) >= abs_obs:
            count += 1
    return (count + 1) / (n_perm + 1)


def _classify_profile(slope: float, p_value: float) -> str:
    if np.isnan(p_value):
        return "Unknown"
    if p_value > 0.05:
        return "Stochastic"
    if slope < -0.02:
        return "Convergence"
    if slope <= 0.02:
        return "Stable"
    return "Directional Change"


def _format_p(p: float) -> str:
    if np.isnan(p):
        return "NA"
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


def _format_formula(slope: float, intercept: float) -> str:
    sign = "+" if intercept >= 0 else "-"
    return f"y = {slope:.2f}x {sign} {abs(intercept):.2f}"


def _format_poly_formula(a: float, b: float, c: float) -> str:
    b_sign = "+" if b >= 0 else "-"
    c_sign = "+" if c >= 0 else "-"
    return f"y = {a:.2f}x² {b_sign} {abs(b):.2f}x {c_sign} {abs(c):.2f}"


def _anova_comparison(x, y, slope, intercept, poly_coeffs, r2_lin, r2_poly, delta_r2):
    n = len(x)
    if n < 4:
        return None, "Linear"

    y_lin_pred = slope * x + intercept
    y_poly_pred = np.polyval(poly_coeffs, x)

    ss_res_lin = np.sum((y - y_lin_pred) ** 2)
    ss_res_poly = np.sum((y - y_poly_pred) ** 2)

    df_diff = 1          # polynomial adds 1 parameter
    df_res = n - 3       # polynomial uses 3 parameters

    if df_res <= 0 or ss_res_poly <= 0:
        return None, "Linear"

    f_val = ((ss_res_lin - ss_res_poly) / df_diff) / (ss_res_poly / df_res)
    p_val = 1 - f_dist.cdf(f_val, df_diff, df_res)

    report = f"F({df_diff}, {df_res}) = {f_val:.2f}, p = {_format_p(p_val)}"
    recommendation = "Poly" if (p_val < 0.05 and delta_r2 > 0.01) else "Linear"

    return report, recommendation
