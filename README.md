# Estimating Treatment Effects under Staggered Timing using GMM

Replication code for the simulation study in Arora and Bijani (2026), "Estimating Treatment Effects under Staggered Timing using GMM."

## Overview

This repository replicates **Table 2** from the paper, which compares the bias and variance of the proposed GMM estimator against existing DiD estimators under both homogeneous and heterogeneous treatment effects with staggered adoption.

## Simulation Design

A balanced panel of **N = 51 units** (inspired by US states) observed over **T = 33 time periods**. The outcome is generated as:

```
Y_it = alpha_i + lambda_t + tau_it + epsilon_it
```

where `alpha_i ~ N(0,1)`, `lambda_t ~ N(0,1)`, `epsilon_it ~ N(0,1)`, and the treatment effect is:

```
tau_it = beta_g * (1 + r_g)^(t - g_i) * D_it
```

**Cohort structure:**
- 5 treated cohorts (8 states each), first treated at T_g in {10, 13, 16, 19, 22}
- 1 never-treated cohort (11 states)
- 500 Monte Carlo simulations

**Homogeneous effects:** `beta_g = -5`, `r_g = 0` for all cohorts.

**Heterogeneous effects:**
- `beta_g in {-16, -12, -10, -9, -2}`
- `r_g in {0.01, 0.04, 0.08, 0.10, 0.07}`

ATT is aggregated using unit-weighted aggregation.

## Estimators Compared

| Abbreviation | Method | R Package |
|---|---|---|
| TWFE | Two-Way Fixed Effects | `fixest` |
| CS | Callaway and Sant'Anna (2021) | `did` |
| SA | Sun and Abraham (2021) | `fixest` (sunab) |
| DCdH | De Chaisemartin and d'Haultfoeuille (2020) | `DIDmultiplegt` |
| Gardner | Gardner et al. (2024) | `did2s` |
| GMM | Proposed estimator (this paper) | Custom implementation |

## Files

- `simulation_table2.R` — Main simulation script that replicates Table 2
- `README.md` — This file

## Requirements

R (>= 4.0) with the following packages:

```r
install.packages(c("data.table", "fixest", "did", "did2s",
                    "DIDmultiplegt", "MASS", "Matrix"))
```

## Usage

```r
source("simulation_table2.R")
```

The script will run 500 simulations for each setting (homogeneous and heterogeneous) and print the replicated Table 2 alongside the paper's original values. Results are saved to `simulation_results.RData`.

**Note:** The full simulation takes considerable time due to 500 iterations x 6 estimators x 2 settings. The DCdH estimator is the slowest component.

## Reference

Arora, P. and Bijani, R. (2026). Estimating Treatment Effects under Staggered Timing using GMM. Working Paper, Ashoka University.
