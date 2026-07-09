# =============================================================================
# TimeLagAnalysis.R
# =============================================================================
# Time Lag Analysis: A tool for exploring inter-game variability in sport.
#
# Quantifies how performance diverges or stabilises across increasing temporal
# distances between matches, using pairwise Euclidean distance on standardised
# KPIs regressed against the square-root of the lag period.
#
# REQUIRED INPUT COLUMNS
#   team   : character — group identifier (e.g. team name)
#   date   : date     — match date, used to order games within each group
#   <KPIs> : numeric  — one or more performance indicator columns
#   result : character (optional) — match outcome; enables win/loss colouring
#
# OUTPUTS  (returned as a named list)
#   lagData              : pairwise distances, lag values, and model labels
#   profileSummary       : regression coefficients, R², p-values, and profiles
#   descriptiveStatistics: mean, SD, Lag-1 MAD, and CoV per team × KPI
#
# LINEAR PROFILES
#   Stochastic        : p > 0.05  — no significant relationship with lag
#   Convergence       : slope < -0.02 — performance stabilising over time
#   Stable            : -0.02 <= slope <= 0.02 — consistent performance
#   Directional Change: slope > 0.02  — performance diverging over time
#
# SIGNIFICANCE
#   Because the pairwise distances are not mutually independent (each game
#   contributes to multiple comparisons), parametric p-values are inflated.
#   The linear-slope significance that drives profile classification is therefore
#   assessed by a permutation test on game order, with the empirical p-value the
#   proportion of |permuted slopes| >= |observed slope|. When the number of
#   orderings (n!) is small (n <= 8, the usual tournament case) the null is
#   enumerated exactly, so the p-value is exact and seed-independent; longer
#   series fall back to n_perm random reshuffles.
#   The linear-vs-polynomial comparison uses the parametric ANOVA F-test (a plain
#   order shuffle is not a valid null for an isolated nested term) together with
#   visual inspection, as recommended for guarding against overfitting.
#
# DEPENDENCIES
#   install.packages(c("dplyr", "tidyr", "purrr", "vegan", "ggplot2", "readxl"))
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(vegan)
library(ggplot2)
library(readxl)

# -----------------------------------------------------------------------------
# All orderings of 1:n as a matrix (one permutation per row). Used to enumerate
# the exact permutation null for short series.
# -----------------------------------------------------------------------------
.allPerms <- function(n) {
  if (n == 1L) return(matrix(1L, 1L, 1L))
  prev <- .allPerms(n - 1L)
  m <- nrow(prev)
  res <- matrix(0L, m * n, n)
  r <- 1L
  for (a in seq_len(m)) {
    row <- prev[a, ]
    for (pos in 0:(n - 1L)) {
      res[r, ] <- append(row, n, after = pos)
      r <- r + 1L
    }
  }
  res
}

# -----------------------------------------------------------------------------
# Permutation test of significance for a single ordered (z-scored) KPI series.
# Reshuffles game order to build a null distribution of the slope. When the
# number of orderings (n!) is small enough the null is enumerated exactly, so
# the p-value is exact and seed-independent — the usual case for tournament
# series (n <= 8). Longer series fall back to n_perm random permutations. The
# square-root lag values are fixed, so only pairwise distances are recomputed
# and the slope is obtained in closed form (slope = sum(w * distance)).
# -----------------------------------------------------------------------------
timeLagPermTest <- function(values, n_perm = 5000, seed = 42, max_exact = 50000) {
  values <- values[!is.na(values)]
  n <- length(values)
  if (n < 3) return(list(slope = NA_real_, p_linear = NA_real_))

  cb   <- utils::combn(n, 2)
  i    <- cb[1, ]; j <- cb[2, ]
  lagv <- sqrt(abs(i - j))
  lc     <- lagv - mean(lagv)
  wslope <- lc / sum(lc^2)

  obs_slope <- sum(wslope * abs(values[i] - values[j]))

  if (factorial(n) <= max_exact) {                 # exact: enumerate all orderings
    P <- .allPerms(n)
    exact <- TRUE
  } else {                                         # approximate: random permutations
    set.seed(seed)
    P <- matrix(0L, n_perm, n)
    for (b in seq_len(n_perm)) P[b, ] <- sample(n)
    exact <- FALSE
  }

  VP     <- matrix(values[as.vector(t(P))], ncol = n, byrow = TRUE)
  slopes <- as.vector(abs(VP[, i] - VP[, j]) %*% wslope)
  hits   <- sum(abs(slopes) >= abs(obs_slope) - 1e-12)

  list(
    slope    = obs_slope,
    p_linear = if (exact) hits / nrow(P) else (hits + 1) / (n_perm + 1)
  )
}

timeLagAnalysis <- function(df, resultsIncluded = FALSE, poly = FALSE,
                            n_perm = 5000, seed = 42) {
  # --- Identify numeric columns (KPIs) ---
  numeric_cols <- names(df %>% select(where(is.numeric)))
  numeric_cols_kpi <- setdiff(numeric_cols, c("game", "date", "result"))      
  
  # --- Keep mandatory columns + numeric KPIs ---
  keep_cols <- unique(c("date", "team", "result", numeric_cols))
  df <- df %>% select(any_of(keep_cols))      
  
  # --- Arrange by team/date and add game numbers ---
  df <- df %>%
    arrange(team, date) %>%
    group_by(team) %>%
    mutate(game = row_number()) %>%
    ungroup()      
  
  # --- Standardize numeric KPIs ---
  df_scaled <- df %>%
    mutate(across(all_of(numeric_cols_kpi),
                  ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)))      
  
  # Storage
  all_results <- tibble()
  summary_results <- tibble()      
  
  # ============================================================
  # LOOP OVER EACH TEAM
  # ============================================================
  for (t in unique(df$team)) {          
    team_data <- df_scaled %>% filter(team == t)
    numeric_cols_team <- intersect(names(team_data %>% select(where(is.numeric))), numeric_cols_kpi)
    if(length(numeric_cols_team) == 0) next
    
    # ---------- Compute Distance Matrices ----------
    team_results <- map_df(numeric_cols_team, function(colname) {              
      numeric_df <- team_data %>% select(all_of(colname))
      dist_matrix <- vegan::vegdist(numeric_df, method = "euclidean") %>% as.matrix()
      idx <- which(lower.tri(dist_matrix), arr.ind = TRUE)              
      
      lag_df <- tibble(
        row = idx[,1],
        col = idx[,2],
        abs_diff = abs(idx[,1] - idx[,2]),
        distance = dist_matrix[idx]
      ) %>%
        mutate(
          lag = factor(paste0("Lag ", abs_diff)),
          sqrtLag = sqrt(abs_diff),
          performanceIndicator = colname,
          team = t
        )              
      
      # --- Optional winPair feature ---
      if (resultsIncluded && "result" %in% names(team_data)) {
        lag_df <- lag_df %>%
          mutate(winPair = paste(team_data$result[idx[,2]],
                                 team_data$result[idx[,1]],
                                 sep = " and "))
      }              
      lag_df
    })          
    
    # ------------ LINEAR & POLY MODEL INFORMATION PER KPI ------------
    compute_model_info <- function(data, poly = FALSE) {              
      linear_fit <- lm(distance ~ sqrtLag, data = data)
      linear_sum <- summary(linear_fit)
      
      linear_slope <- coef(linear_fit)["sqrtLag"]
      linear_p <- linear_sum$coefficients["sqrtLag", "Pr(>|t|)"]
      linear_r2 <- linear_sum$r.squared
      
      linear_profile <- if (is.na(linear_p)) {
        NA
      } else if (linear_p > 0.05) {
        "Stochastic"
      } else if (linear_slope < -0.02) {
        "Convergence"
      } else if (linear_slope >= -0.02 & linear_slope <= 0.02) {
        "Stable"
      } else {
        "Directional Change"
      }
      
      linear_formula_text <- paste0(
        "y = ", round(linear_slope, 2), "x ",
        if (coef(linear_fit)[1] < 0) {
          paste0("- ", abs(round(coef(linear_fit)[1], 2)))
        } else {
          paste0("+ ", round(coef(linear_fit)[1], 2))
        }
      )
      
      p_value_formatted <- if (!is.na(linear_p)) {
        if (linear_p < 0.001) "<0.001" else sprintf("%.3f", linear_p)
      } else NA
      
      R2_formatted <- if (!is.na(linear_r2)) sprintf("%.2f", linear_r2) else NA
      
      out <- tibble(
        R2 = R2_formatted,
        Linear_p_value = p_value_formatted,
        Linear_Profile = linear_profile,
        Linear_Formula = linear_formula_text
      )
      
      if (poly) {                  
        poly_fit <- lm(distance ~ poly(sqrtLag, 2, raw = TRUE), data = data)
        poly_sum <- summary(poly_fit)
        
        a <- poly_sum$coefficients[3, 1]
        b <- poly_sum$coefficients[2, 1]
        c <- poly_sum$coefficients[1, 1]
        
        poly_formula_text <- paste0(
          "y = ", round(a, 2), "x² ",
          if (b < 0) paste0("- ", abs(round(b, 2)), "x") else paste0("+ ", round(b, 2), "x"),
          " ",
          if (c < 0) paste0("- ", abs(round(c, 2))) else paste0("+ ", round(c, 2))
        )
        
        poly_r2 <- poly_sum$r.squared
        poly_p <- poly_sum$coefficients[3, "Pr(>|t|)"]
        
        poly_R2_formatted <- if (!is.na(poly_r2)) sprintf("%.2f", poly_r2) else NA
        poly_p_formatted <- if (!is.na(poly_p)) {
          if (poly_p < 0.001) "<0.001" else sprintf("%.3f", poly_p)
        } else NA
        
        ANOVA_Report <- NA
        Model_Recommendation <- "Linear"
        Delta_R2 <- NA
        
        try({
          anova_result <- anova(linear_fit, poly_fit)
          if (!is.null(anova_result) && nrow(anova_result) >= 2) {
            df_diff <- anova_result$Df[2]
            df_res  <- anova_result$Res.Df[2]
            F_val   <- anova_result$F[2]
            p_val   <- anova_result$`Pr(>F)`[2]
            
            if (!is.na(F_val) && is.finite(F_val) && df_res > 0) {
              ANOVA_Report <- sprintf("F(%d, %d) = %.2f, p = %s",
                                      df_diff, df_res, F_val,
                                      ifelse(!is.na(p_val) && p_val < 0.001, "<0.001",
                                             sprintf("%.3f", p_val)))
            }
            
            Delta_R2 <- poly_r2 - linear_r2
            if (!is.na(p_val) && p_val < 0.05 && Delta_R2 > 0.01) {
              Model_Recommendation <- "Poly"
            }
          }
        }, silent = TRUE)
        
        out <- out %>% mutate(
          poly_R2 = poly_R2_formatted,
          poly_p_value = poly_p_formatted,
          poly_Formula = poly_formula_text,
          ANOVA_Report = ANOVA_Report,
          Model_Recommendation = Model_Recommendation,
          Delta_R2 = if (!is.na(Delta_R2)) sprintf("%.2f", Delta_R2) else NA
        )
      }
      
      return(out)
    }
    
    model_info_df <- team_results %>%
      group_by(performanceIndicator, team) %>%
      group_modify(~ compute_model_info(.x, poly = poly)) %>%
      ungroup()

    # ---- Permutation test of significance for the linear slope ----
    # A game-order shuffle is the correct null for the directional slope. The
    # linear-vs-polynomial comparison is left on the parametric ANOVA F-test
    # (see compute_model_info), since a plain order shuffle does not isolate the
    # quadratic term in a nested model. Profiles below use the empirical p.
    perm_df <- map_df(numeric_cols_team, function(colname) {
      pr <- timeLagPermTest(team_data[[colname]], n_perm = n_perm, seed = seed)
      tibble(
        performanceIndicator = colname,
        team                 = t,
        obs_slope            = pr$slope,
        Linear_p_perm        = pr$p_linear
      )
    })

    model_info_df <- model_info_df %>%
      left_join(perm_df, by = c("performanceIndicator", "team")) %>%
      mutate(
        Linear_Profile = dplyr::case_when(
          is.na(Linear_p_perm) ~ NA_character_,
          Linear_p_perm > 0.05 ~ "Stochastic",
          obs_slope < -0.02    ~ "Convergence",
          obs_slope <=  0.02   ~ "Stable",
          TRUE                 ~ "Directional Change"
        ),
        Linear_p_value = ifelse(is.na(Linear_p_perm), NA_character_,
                                ifelse(Linear_p_perm < 0.001, "<0.001",
                                       sprintf("%.3f", Linear_p_perm)))
      ) %>%
      select(-any_of(c("obs_slope", "Linear_p_perm")))

    team_results <- team_results %>%
      left_join(
        model_info_df %>% select(any_of(c(
          "performanceIndicator",
          "team",
          "Model_Recommendation",
          "Linear_p_value",
          "Linear_Profile",
          "Linear_Formula",
          "poly_Formula",
          "poly_p_value",
          "ANOVA_Report",
          "Delta_R2"
        ))),
        by = c("performanceIndicator", "team")
      ) %>%
      mutate(
        p_value_num = suppressWarnings(as.numeric(Linear_p_value)),
        Significance = ifelse(
          Linear_p_value == "<0.001" | (!is.na(p_value_num) & p_value_num < 0.05),
          "Significant", "Not Significant"
        )
      ) %>%
      select(-p_value_num) %>%
      arrange(performanceIndicator, lag)
    
    summary_results <- bind_rows(summary_results, model_info_df)
    team_results <- team_results %>% select(-row, -col, -abs_diff)
    all_results <- bind_rows(all_results, team_results)
  }
  
  # DESCRIPTIVE STATISTICS TABLE
  # DESCRIPTIVE STATISTICS TABLE
  descriptive_table <- df %>%
    select(team, all_of(numeric_cols_kpi)) %>%
    pivot_longer(cols = -team,
                 names_to = "performanceIndicator",
                 values_to = "value") %>%
    group_by(team, performanceIndicator) %>%
    summarise(
      Mean = mean(value, na.rm = TRUE),
      SD   = sd(value, na.rm = TRUE),
      Lag1_MAD = {
        # Lag-1 difference on z-scored values
        z <- (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)
        lag_diff <- abs(z - dplyr::lag(z))
        mean(lag_diff[-1], na.rm = TRUE)  # remove first NA from lag
      },
      CoV = ifelse(mean(value, na.rm = TRUE) != 0,
                   sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE) * 100,
                   NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      Mean = round(Mean, 1),
      SD   = round(SD, 2),
      Lag1_MAD  = round(Lag1_MAD, 3),
      CoV  = round(CoV, 1)
    )
  
  
  summary_table <- summary_results %>%
    mutate(
      p_value_num = suppressWarnings(as.numeric(Linear_p_value)),
      Linear_p_value = ifelse(!is.na(p_value_num),
                              ifelse(p_value_num < 0.001, "<0.001",
                                     sprintf("%.3f", p_value_num)),
                              Linear_p_value),
      R2_num = suppressWarnings(as.numeric(R2)),
      R2 = ifelse(!is.na(R2_num), sprintf("%.2f", R2_num), R2)
    ) %>%
    select(-p_value_num, -R2_num)
  
  list(
    lagData = all_results,
    profileSummary = summary_table,
    descriptiveStatistics = descriptive_table
  )
}

theme_time_lag <- function() {
  theme(
    text = element_text(family = "Times New Roman"),
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
    axis.title.y.right = element_text(size = 12, margin = margin(l = 10)),
    
    # Panel border (soft grey)
    panel.border = element_rect(
      color = "grey70",
      fill = NA,
      linewidth = 0.4
    ),
    panel.background = element_blank(),
    panel.spacing = unit(1, "lines"),
    
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(size = 12),
    strip.text.y.left = element_text(angle = 90),
    
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    
    panel.grid = element_blank(),
    
    # Axis lines — match panel border
    axis.line = element_line(colour = "grey70", linewidth = 0.4),

    # Axis ticks — match panel border
    axis.ticks = element_line(colour = "grey70", linewidth = 0.4),

    # Tidy tick length
    axis.ticks.length = unit(3, "pt"),
    
    strip.placement = "outside"
  )
}

timeLagLinearPlot <- function(df) {
  
  
  # Validation check
  #--------------------------------------------------
  required_cols <- c(
    "sqrtLag",
    "distance",
    "Significance",
    "team",
    "performanceIndicator"
  )
  
  if (!all(required_cols %in% names(df))) {
    stop(
      "Input dataframe must contain: ",
      paste(required_cols, collapse = ", ")
    )
  }
  
  #--------------------------------------------------
  # winPair colour mapping (only if winPair exists)
  #--------------------------------------------------
  has_winpair <- "winPair" %in% names(df)
  
  if (has_winpair) {
    winpair_colours <- c(
      "Win and Win"   = "#2C7BB6",
      "Loss and Win"  = "#2C7BB6",
      "Win and Loss"  = "#D7191C",
      "Loss and Loss" = "#D7191C"
    )
  }
  
  #--------------------------------------------------
  # Begin plot
  #--------------------------------------------------
  p <- ggplot(df, aes(x = sqrtLag, y = distance))
  
  # Points
  if (has_winpair) {
    p <- p + geom_point(aes(color = winPair), size = 1.2, alpha = 0.8)
  } else {
    p <- p + geom_point(color = "darkgrey", size = 1.2, alpha = 0.8)
  }
  
  # Regression line (Significance)
  p <- p +
    geom_smooth(
      aes(color = Significance),
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.9
    )
  
  # Facets
  p <- p + facet_grid(
    team ~ performanceIndicator,
    scales = "free_x",
    switch = "y"
  )
  
  # Axes & limits
  p <- p +
    coord_cartesian(ylim = c(0, 6)) +
    scale_y_continuous(position = "right") +
    labs(
      x = "Time Lag (sqrt)",
      y = "Distance (z-score)"
    )

  # Per-panel p-value labels (uses the permutation p in Linear_p_value, if present)
  if ("Linear_p_value" %in% names(df)) {
    pval_df <- df %>%
      group_by(team, performanceIndicator) %>%
      summarise(Linear_p_value = dplyr::first(Linear_p_value), .groups = "drop") %>%
      mutate(label = ifelse(grepl("^<", Linear_p_value),
                            paste0("p", Linear_p_value),
                            paste0("p=", Linear_p_value)))
    p <- p + geom_text(
      data = pval_df,
      aes(x = Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = 1.15, vjust = 1.5,
      size = 4.5,
      family = "Times New Roman"
    )
  }

  # Colour scales
  if (has_winpair) {
    p <- p + scale_color_manual(
      values = c(
        winpair_colours,
        "Significant" = "black",
        "Not Significant" = "grey60"
      ),
      breaks = c(
        "Significant", "Not Significant",
        names(winpair_colours)
      ),
      name = NULL
    )
  } else {
    p <- p + scale_color_manual(
      values = c(
        "Significant" = "black",
        "Not Significant" = "grey60"
      ),
      breaks = c("Significant", "Not Significant"),
      name = NULL
    )
  }
  
  #--------------------------------------------------
  # Apply shared theme
  #--------------------------------------------------
  p <- p +
    theme_time_lag() +
    theme(
      legend.position = "top"
    )
  
  return(p)
}

timeLagPolyPlot <- function(df) {
  #---------------------Validation-----------------------------
  required_cols <- c(
    "sqrtLag",
    "distance",
    "team",
    "performanceIndicator",
    "ANOVA_Report",
    "Delta_R2"
  )
  
  if (!all(required_cols %in% names(df))) {
    stop(
      "Input dataframe must contain: ",
      paste(required_cols, collapse = ", ")
    )
  }
  
  #---------------------Annotation labels----------------------
  annotation_df <- df %>%
    group_by(team, performanceIndicator) %>%
    summarise(
      ANOVA_Report = unique(ANOVA_Report),
      Delta_R2 = unique(Delta_R2),
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0(ANOVA_Report, "\nΔR² = ", Delta_R2),
      y_pos = Inf - 10  # small offset below top
    )
  
  #----------------------Begin plot-----------------------------
  p <- ggplot(df, aes(x = sqrtLag, y = distance)) +
    
    # Points
    geom_point(
      color = "grey50",
      size = 1,
      alpha = 0.8
    ) +
    
    # Linear regression
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 0.8,
      color = "grey20",
      linetype = "dashed"
    ) +
    
    # Polynomial regression (degree 2)
    geom_smooth(
      method = "lm",
      formula = y ~ poly(x, 2),
      se = FALSE,
      linewidth = 0.8,
      color = "grey70",
      linetype = "solid"
    ) +
    
    # Facets
    facet_grid(
      team ~ performanceIndicator,
      scales = "free_x",
      switch = "y"
    ) +
    
    # Annotation in top-right, slightly below top
    geom_text(
      data = annotation_df,
      aes(x = Inf, y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 3.6,
      hjust = 1,
      vjust = 1
    ) +
    
    # Axes & limits
    coord_cartesian(ylim = c(0, 6)) +
    scale_y_continuous(position = "right") +
    labs(
      x = "Time Lag (sqrt)",
      y = NULL
    ) +
    
    # Shared theme
    theme_time_lag()
  
  return(p)
}

