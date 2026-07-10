# =============================================================================
# Method validation for time-lag analysis (simulation study)
# -----------------------------------------------------------------------------
# Purpose: test whether the time-lag classification recovers KNOWN temporal
# profiles at the small sample sizes typical of tournaments (n = 6 and n = 8).
# We generate data with a deliberately imposed structure, run it through the
# identical pipeline used on the real data, and check the verdict matches.
#
# Reproducible: set.seed() fixed below. Uses only base R (no packages).
# Companion write-up: docs/method_validation.md
# =============================================================================

set.seed(2024)

# ---- Core time-lag routine (identical logic to the main analysis) -----------
# Precompute the fixed design (lags and sqrt-lags) for a series of length n.
make_design <- function(n) {
  idx <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)  # all game pairs i<j
  sx  <- sqrt(abs(idx[, "col"] - idx[, "row"]))             # sqrt(lag) per pair
  list(row = idx[, "row"], col = idx[, "col"], sx = sx,
       sxc = sx - mean(sx), Sxx = sum((sx - mean(sx))^2))
}

# Closed-form slope of dissimilarity ~ sqrt(lag) (fast; avoids refitting lm).
slope_from_z <- function(z, D) sum(D$sxc * abs(z[D$col] - z[D$row])) / D$Sxx

# Permutation test on the slope: shuffle game order, as in Collins et al. (2000).
perm_p <- function(x, D, nperm = 5000) {
  z <- as.numeric(scale(x)); n <- length(z)
  obs  <- slope_from_z(z, D)
  null <- replicate(nperm, slope_from_z(z[sample.int(n)], D))
  list(slope = obs, p = mean(abs(null) >= abs(obs)))
}

# Four-way classification following the manuscript's decision tree.
classify <- function(slope, p, alpha = 0.05, thr = 0.02) {
  if (p > alpha)          return("stochastic")
  if (abs(slope) <= thr)  return("stable")
  if (slope > thr)        return("directional")
  "convergence"
}

# Linear-vs-quadratic model preference (Case Study 2 non-linear step).
poly_pref <- function(x, D) {
  z <- as.numeric(scale(x)); diss <- abs(z[D$col] - z[D$row])
  anova(lm(diss ~ D$sx), lm(diss ~ D$sx + I(D$sx^2)))$`Pr(>F)`[2] < 0.05
}

# ---- Generators: series with a KNOWN temporal profile -----------------------
gen <- list(
  # No temporal structure: independent noise.
  stochastic  = function(n) rnorm(n),
  # Systematic trend across games (rate b): produces rising dissimilarity.
  directional = function(n, b = 0.6) b * (1:n) + rnorm(n),
  # Inverted-U: flat -> shift up -> return, reproducing the football pattern.
  invU        = function(n) {
    f <- floor(n / 3)
    base <- c(rep(0, f), rep(2.5, n - 2 * f), rep(0, f))[1:n]
    base + rnorm(n, sd = 0.5)
  }
)

# ---- Monte Carlo evaluation -------------------------------------------------
LEVELS <- c("stochastic", "stable", "directional", "convergence")
REPS   <- 2000    # replicates per profile (rates, not single outcomes)
NPERM  <- 5000    # permutations per test

run_cell <- function(n) {
  D <- make_design(n)
  cat(sprintf("\n==== n = %d  (REPS = %d, NPERM = %d) ====\n", n, REPS, NPERM))
  for (nm in names(gen)) {
    g <- gen[[nm]]
    res <- t(replicate(REPS, {
      x  <- g(n)
      pr <- perm_p(x, D, NPERM)
      c(classify(pr$slope, pr$p), poly_pref(x, D))
    }))
    tb   <- table(factor(res[, 1], levels = LEVELS))
    poly <- mean(res[, 2] == "TRUE")
    cat(sprintf("truth = %-11s | stoch=%4d stable=%3d direct=%4d conv=%3d | poly-preferred = %4.1f%%\n",
                nm, tb["stochastic"], tb["stable"], tb["directional"], tb["convergence"], 100 * poly))
  }
}

for (n in c(6, 8)) run_cell(n)

# ---- Slope bias vs series length (true rate b = 0.6) ------------------------
cat("\n==== slope estimate vs series length (true b = 0.6) ====\n")
bias <- do.call(rbind, lapply(c(6, 8, 12, 20, 40), function(n) {
  D <- make_design(n)
  s <- replicate(REPS, slope_from_z(as.numeric(scale(gen$directional(n))), D))
  data.frame(n = n, mean_slope = round(mean(s), 3), sd_slope = round(sd(s), 3))
}))
print(bias, row.names = FALSE)
cat("\nDONE\n")
