###############################################################################
# Pooled TWFE vs Goodman-Bacon Decomposition — Homogeneous Effects DGP
#
# DGP: Y_it = alpha_i + lambda_t + tau * D_it + epsilon_it
#   where tau = 1 (homogeneous, constant treatment effect)
#         D_it = 1 if unit i is treated and t >= cohort_g_i
#         epsilon_it ~ iid N(0, 1)
#
# Purpose:
#   Under homogeneous (constant) treatment effects, pooled TWFE is consistent
#   and recovers tau. The Goodman-Bacon decomposition expresses the TWFE
#   estimate as a weighted average of ALL pairwise 2x2 DiD estimates
#   (including "forbidden" comparisons that use already-treated units as
#   controls). The two should be numerically identical since the Goodman-Bacon
#   decomposition is an exact algebraic identity for TWFE, not an approximation.
#
# Estimators:
#   1. TWFE:         pooled two-way FE (feols with unit + time FE)
#   2. GB:           Goodman-Bacon weighted average (via bacon() from bacondecomp)
#
# Simulation: 100 iterations, track point estimates and check:
#   - Bias of each estimator relative to true tau = 1
#   - Variance across replications
#   - Whether TWFE == GB numerically in every replication
###############################################################################

library(data.table)
library(fixest)
library(bacondecomp)   # install.packages("bacondecomp") if needed

set.seed(42)

# ===========================================================================
# 1. Parameters
# ===========================================================================

tau_true        <- 1.0      # homogeneous treatment effect
cohort_size     <- 20       # units per treatment cohort
n_never         <- 20       # never-treated units
treatment_times <- c(5, 8, 11, 14)   # cohort treatment start periods
n_cohorts       <- length(treatment_times)
T_total         <- 20       # total time periods
N_total         <- n_cohorts * cohort_size + n_never
n_sims          <- 100

# Unit-to-cohort mapping: 0 means never-treated
unit_cohort <- c(rep(treatment_times, each = cohort_size), rep(0L, n_never))
unit_ids    <- seq_len(N_total)
time_ids    <- seq_len(T_total)

# ===========================================================================
# 2. Data generation
# ===========================================================================

generate_data <- function() {
  # Unit and time fixed effects
  alpha_i <- rnorm(N_total)                   # unit FE
  lambda_t <- rnorm(T_total)                  # time FE
  eps      <- rnorm(N_total * T_total)        # iid N(0,1) errors

  # Build panel (long format)
  dt <- data.table(
    id   = rep(unit_ids,  each  = T_total),
    time = rep(time_ids,  times = N_total),
    g    = rep(unit_cohort, each = T_total)   # cohort (0 = never)
  )

  # Treatment indicator
  dt[, D := as.integer(g > 0 & time >= g)]

  # Outcome
  dt[, Y := alpha_i[id] + lambda_t[time] + tau_true * D + eps]

  return(dt)
}

# ===========================================================================
# 3. Simulation loop
# ===========================================================================

results <- data.table(
  sim    = seq_len(n_sims),
  TWFE   = NA_real_,
  GB     = NA_real_
)

for (s in seq_len(n_sims)) {
  dt <- generate_data()

  # --- TWFE ------------------------------------------------------------------
  fit_twfe  <- feols(Y ~ D | id + time, data = dt, warn = FALSE, notes = FALSE)
  twfe_est  <- coef(fit_twfe)["D"]

  # --- Goodman-Bacon decomposition -------------------------------------------
  # bacon() requires: outcome, treatment dummy, unit id, time id
  # It returns a data frame of pairwise 2x2 DiD estimates and their weights.
  # The weighted average of those equals the TWFE coefficient exactly.
  gb_out <- bacon(Y ~ D, data = dt, id_var = "id", time_var = "time")
  gb_est <- sum(gb_out$estimate * gb_out$weight)   # weighted average = TWFE

  results[s, `:=`(TWFE = twfe_est, GB = gb_est)]
}

# ===========================================================================
# 4. Summary statistics
# ===========================================================================

cat("\n=== Simulation Results (n_sims =", n_sims, ", tau_true =", tau_true, ") ===\n\n")

cat(sprintf("%-25s %10s %10s %12s\n", "Estimator", "Mean", "Bias", "Std Dev"))
cat(strrep("-", 60), "\n")

for (col in c("TWFE", "GB")) {
  vals <- results[[col]]
  cat(sprintf("%-25s %10.6f %10.6f %12.6f\n",
              col,
              mean(vals),
              mean(vals) - tau_true,
              sd(vals)))
}

cat("\n=== Numerical Equivalence Check ===\n")
diffs    <- abs(results$TWFE - results$GB)
cat(sprintf("Max |TWFE - GB| across %d replications: %.2e\n", n_sims, max(diffs)))
cat(sprintf("Mean |TWFE - GB|:                        %.2e\n", mean(diffs)))
cat(sprintf("All differences < 1e-8: %s\n", ifelse(all(diffs < 1e-8), "YES", "NO")))

cat("\n=== Per-Replication Detail (first 10 sims) ===\n")
cat(sprintf("%-6s %12s %12s %14s\n", "Sim", "TWFE", "GB", "|diff|"))
cat(strrep("-", 48), "\n")
for (s in seq_len(min(10, n_sims))) {
  cat(sprintf("%-6d %12.6f %12.6f %14.2e\n",
              s,
              results$TWFE[s],
              results$GB[s],
              abs(results$TWFE[s] - results$GB[s])))
}
