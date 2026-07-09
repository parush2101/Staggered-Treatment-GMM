# Estimating Treatment Effects under Staggered Timing and Non-Spherical Errors

Replication code for Arora and Bijani (2026), *Estimating Treatment Effects under Staggered Timing and Non-Spherical Errors.*

The paper develops a GMM framework for the ATT under staggered adoption with heterogeneous treatment effects and **non-spherical** (serially and/or cross-sectionally correlated) errors. Forbidden comparisons are shown to be redundant, the clean comparisons collapse to a design-sized minimal set, and optimally weighting that set delivers efficiency under arbitrary within-unit serial correlation.

## Repository structure

```
simulations/            Main Monte Carlo studies (Section 5)
appendix_simulations/   Extension / robustness studies (Appendices B and C)
archive/                Earlier exploratory scripts and superseded drafts
README.md
```

## Main simulations (`simulations/`)

All three scripts share the same design: a balanced panel over **T = 33** periods with **five treated cohorts** first treated at t ∈ {10, 13, 16, 19, 22} plus a **never-treated** cohort, and

```
Y_it = alpha_i + lambda_t + tau_it + eps_it,   tau_it = beta_g (1 + r_g)^(t - g) D_it
```

with heterogeneous `beta_g ∈ {-16,-12,-10,-9,-2}`, `r_g ∈ {0.01,0.04,0.08,0.10,0.07}`, `alpha_i, lambda_t ~ N(0,1)`, and **500 replications**. They differ only in the error covariance:

| Script | Error structure | Paper |
|---|---|---|
| `hetcov_dgp2_n300_r09_het_01_06_copy.R` | AR(1), ρ ∈ {0, .1, .3, .5, .7, .9} | Table 1 |
| `hetcov_dgp2_arbitrary_serial.R` | Arbitrary serial (positive mixture of AR kernels), swept low / medium / high | Table 3 |
| `hetcov_dgp1_cross_sectional_sweep.R` | Cross-sectional common-factor, σ_Φ × ρ (15 cells) | Tables 4–5 |

Each script computes the three GMM weightings and the parallel iterated-GLS estimators, alongside the existing heterogeneity-robust estimators (CS, SA, Gardner / flexible TWFE) and TWFE. **Code labels map to the paper as:**

| Code label | Paper | Weighting |
|---|---|---|
| GMM I / GLS I | GMM<sub>T</sub> / GLS<sub>T</sub> | homoskedastic Toeplitz |
| GMM II / GLS II | GMM<sub>HT</sub> / GLS<sub>HT</sub> | heteroskedastic Toeplitz |
| GMM III / GLS III | GMM<sub>U</sub> / GLS<sub>U</sub> | unrestricted within-cohort |

Sample size is set by `cohort_size` (50 → N = 300 large-sample tables; set `cohort_size = 10, n_never = 10` for the small-sample tables).

## Appendix simulations (`appendix_simulations/`)

Self-contained **base-R** scripts (no external packages).

| Script | Study | Paper |
|---|---|---|
| `jtest_pretrends_simulation.R` | Size and power of the J specification test for parallel trends / no anticipation; full placebo set (df = 70) vs. local pre-window (df = 14) | Appendix C |
| `switchers_arbitrary_serial.R` | Reversible-treatment (switchers) design; DID<sub>M</sub> vs. over-identified GMM<sub>eq</sub> / GMM<sub>sph</sub> / GMM<sub>T</sub>, with the variance-gain decomposition | Appendix B |

Run a quick pass with `QUICK=1 Rscript <file>`. Environment overrides: `N_SIMS`, `STRENGTH` ∈ {low, medium, high}; `jtest_pretrends_simulation.R` also takes `PRE_WINDOW` (placebo-window length).

## Estimators

| Abbrev. | Method | R package |
|---|---|---|
| TWFE | Two-way fixed effects | `fixest` |
| CS | Callaway and Sant'Anna (2021) | `did` |
| SA | Sun and Abraham (2021) | `fixest` (`sunab`) |
| DID<sub>M</sub> | de Chaisemartin and D'Haultfœuille (2020) | base R (appendix) |
| Gardner / Flex TWFE | Gardner (2022); Wooldridge (2025); Borusyak et al. (2024) | `did2s` |
| GMM<sub>T/HT/U</sub> | Proposed estimator (this paper) | custom |
| GLS<sub>T/HT/U</sub> | Iterated feasible-GLS benchmark | custom |

## Requirements

R (≥ 4.0). The **main** scripts use:

```r
install.packages(c("data.table", "fixest", "did", "did2s", "MASS"))
```

The **appendix** scripts use base R only.

## Usage

```r
# main simulations (parallel; set N_CORES for cores used)
N_CORES=4 Rscript simulations/hetcov_dgp2_arbitrary_serial.R

# appendix simulations (quick check)
QUICK=1 Rscript appendix_simulations/jtest_pretrends_simulation.R
Rscript appendix_simulations/switchers_arbitrary_serial.R
```

The main scripts run 500 replications per configuration and print the paper tables; they are computationally heavy (the GLS benchmarks are the slowest component).

## `archive/`

Earlier exploratory scripts and superseded LaTeX section drafts, retained for provenance. Not needed to reproduce the paper's results.

## Reference

Arora, P. and Bijani, R. (2026). *Estimating Treatment Effects under Staggered Timing and Non-Spherical Errors.* Working Paper, Ashoka University.
