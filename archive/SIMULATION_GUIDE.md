# Simulation Guide

## Table 2 Replication

### What the simulation does

The script `simulation_table2.R` replicates Table 2 from the paper by running 500 Monte Carlo simulations under two data-generating processes (DGPs):

1. **Homogeneous effects**: All cohorts have the same constant treatment effect (`beta_g = -5`, `r_g = 0`). Under this DGP, all estimators should be unbiased and the GMM estimator is expected to have the lowest variance.

2. **Heterogeneous effects**: Cohorts have different baseline effects and growth rates (`beta_g in {-16,-12,-10,-9,-2}`, `r_g in {0.01,0.04,0.08,0.10,0.07}`). Under this DGP, TWFE is biased due to forbidden comparisons, while the heterogeneity-robust estimators (CS, SA, DCdH, Gardner) and GMM remain unbiased. The GMM estimator is expected to have the lowest variance among the unbiased estimators.

### Simulation parameters

| Parameter | Value |
|---|---|
| N (units) | 51 |
| T (periods) | 33 |
| Treated cohorts | 5 (8 states each), first treated at {10, 13, 16, 19, 22} |
| Never-treated | 11 states |
| Simulations | 500 |
| Fixed effects | alpha_i, lambda_t ~ N(0,1) |
| Error | epsilon_it ~ N(0,1) |
| Aggregation | Unit-weighted |

### Estimators

| Label | Method | Implementation |
|---|---|---|
| TWFE | Two-way fixed effects | `fixest::feols()` |
| CS | Callaway & Sant'Anna (2021) | `did::att_gt()` + `aggte()` |
| SA | Sun & Abraham (2021) | `fixest::feols()` with `sunab()` |
| DCdH | De Chaisemartin & d'Haultfoeuille (2020) | `DIDmultiplegt::did_multiplegt()` or manual switchers |
| Gardner | Gardner et al. (2024) | `did2s::did2s()` |
| GMM | Proposed (this paper) | Custom: builds Q_H incidence matrix, iterative GMM |

### How the GMM estimator works

1. **Enumerate all 2x2 DiD estimates**: For each CATT `beta_{g,g+k}`, we compute all possible 2x2 DiD estimates using never-treated, not-yet-treated, and already-treated (forbidden) comparisons across all valid pre-treatment reference periods.

2. **Build the incidence matrix Q_H**: Under homogeneous effects, Q_H = Q (simple incidence). Under heterogeneous effects, the rows corresponding to forbidden comparisons are corrected: entries of +1 and -1 are added to absorb the contaminating CATTs from the already-treated control group.

3. **GMM estimation**: Minimize `phi(beta)' A phi(beta)` where `phi(beta) = Delta - Q_H * beta`. Start with A = I (OLS), then iterate with residual-based diagonal weighting.

4. **Aggregate**: Compute ATT as a unit-weighted average of estimated CATTs.

### Expected results (from paper)

| Estimator | Hom. Bias | Hom. Var | Het. Bias | Het. Var |
|---|---|---|---|---|
| TWFE | 0.0020 | 0.009 | 3.5600 | 0.009 |
| CS | 0.0170 | 0.045 | -0.018 | 0.540 |
| SA | 0.0080 | 0.049 | -0.006 | 0.406 |
| DCdH | -0.0110 | 0.045 | 0.0025 | 0.431 |
| Gardner | -0.0070 | 0.0104 | -0.007 | 0.0104 |
| GMM | 0.0072 | 0.0078 | 0.0100 | 0.0098 |

### Running

```r
# Full simulation (takes several hours)
source("simulation_table2.R")

# Reduced simulation (for testing)
n_sims <- 20
source("simulation_table2.R")
```

Results are saved to `simulation_results.RData`.
