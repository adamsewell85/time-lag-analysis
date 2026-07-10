# Method validation by simulation

Companion documentation for `r/validation_simulation.R`. This is the detailed
version kept for the thesis (appendix / supplementary results); the journal
manuscript reports only a condensed summary.

## Rationale
The temporal structure of real performance data is never known with certainty,
so the classification cannot be validated on the tournament data itself. Instead
the method is tested on simulated series in which the underlying profile is fixed
in advance. Because the generating pattern is known, agreement between the
classification and that pattern gives a direct measure of accuracy. Critically,
the simulated series pass through the *identical* pipeline used on the real data
(standardise → pairwise dissimilarity → regress on √lag → permutation test →
classify), so the pipeline is blind to how each series was produced.

## Generators (known-truth profiles)
| Profile | Construction | Represents |
|---|---|---|
| Stochastic | `rnorm(n)` | No temporal structure |
| Directional | `0.6 * (1:n) + rnorm(n)` | Systematic game-to-game trend |
| Inverted-U | flat → +2.5 shift → return, + noise | The football attacking pattern (Case Study 2) |

Convergence (a negative dissimilarity–lag slope) is an uncommon anti-persistent
structure that did not occur in the data and is not simulated; this is stated
rather than contrived.

## Evaluation
Each profile was generated 2,000 times at n = 6 and n = 8, run through the
pipeline (5,000 permutations per test), and the classification tallied. Four
properties were assessed: false-positive rate, statistical power, non-linear
detection (polynomial preference), and slope stability.

## Results (seed = 2024; REPS = 2000; NPERM = 5000)

### Classification outcomes
| Truth | n | stochastic | stable | directional | convergence | polynomial preferred |
|---|---|---|---|---|---|---|
| stochastic  | 6 | 1899 | 0 | 95   | 6 | 2.2% |
| directional | 6 | 1095 | 0 | 905  | 0 | 1.3% |
| inverted-U  | 6 | 2000 | 0 | 0    | 0 | 60.9% |
| stochastic  | 8 | 1882 | 0 | 111  | 7 | 3.8% |
| directional | 8 | 287  | 0 | 1713 | 0 | 11.6% |
| inverted-U  | 8 | 1999 | 0 | 1    | 0 | 87.6% |

Derived headline figures:
- **False-positive rate:** 5.05% (n = 6), 5.90% (n = 8) — near the nominal 5%.
- **Power (directional):** 45.3% (n = 6), 85.7% (n = 8).
- **Inverted-U detection:** 60.9% (n = 6), 87.6% (n = 8); false polynomial
  preference under noise 2.2–3.8%.
- The inverted-U is (correctly) classified stochastic by the *linear* step —
  it is recovered only by the polynomial ANOVA, which is why Case Study 2 needs
  both steps.

### Slope estimate vs series length (true b = 0.6)
| n | mean slope | sd |
|---|---|---|
| 6  | 0.834 | 0.519 |
| 8  | 0.893 | 0.282 |
| 12 | 0.888 | 0.104 |
| 20 | 0.769 | 0.027 |
| 40 | 0.570 | 0.004 |

Slope estimates are upwardly biased and highly variable at n = 6–8, supporting
the manuscript's emphasis on profile classification rather than slope magnitude.

## Reproduce
```
Rscript r/validation_simulation.R
```
Runtime ~2 min on a laptop. Numbers are deterministic given the seed.
