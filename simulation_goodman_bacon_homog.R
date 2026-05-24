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
#   The third estimator (GB_Clean) uses only the "clean" subset of 2x2 DiDs —
#   those whose control group is either never-treated or not-yet-treated. The
#   "forbidden" comparison (type = "Later vs Earlier Treated") is dropped and
#   the remaining weights are renormalized to sum to 1. Under homogeneous
#   effects, every clean 2x2 DiD also estimates tau, so GB_Clean is unbiased
#   with zero bias and comparable variance, demonstrating that forbidden
#   comparisons carry no additional information (they are algebraically
#   redundant linear combinations of clean DiDs).
#
# Estimators:
#   1. TWFE:         pooled two-way FE (feols with unit + time FE)
#   2. GB:           Goodman-Bacon weighted average — ALL 2x2 DiDs (exact = TWFE)
#   3. GB_Clean:     Goodman-Bacon weighted average — clean DiDs only (renormalized)
#                    Forbidden type "Later vs Earlier Treated" is dropped.
#
# Simulation: 100 iterations, track point estimates and check:
#   - Bias of each estimator relative to true tau = 1
#   - Variance across replications
#   - Whether TWFE == GB numerically in every replication
#   - Whether GB_Clean recovers tau despite discarding forbidden comparisons
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
  sim      = seq_len(n_sims),
  TWFE     = NA_real_,
  GB       = NA_real_,
  GB_Clean = NA_real_
)

for (s in seq_len(n_sims)) {
  dt <- generate_data()

  # --- TWFE ------------------------------------------------------------------
  fit_twfe  <- feols(Y ~ D | id + time, data = dt, warn = FALSE, notes = FALSE)
  twfe_est  <- coef(fit_twfe)["D"]

  # --- Goodman-Bacon decomposition -------------------------------------------
  # bacon() returns one row per pairwise 2x2 DiD with its weight and type.
  # Types produced by bacondecomp:
  #   "Earlier Treated vs Untreated"  — early cohort vs never-treated (CLEAN)
  #   "Later Treated vs Untreated"    — late cohort vs never-treated  (CLEAN)
  #   "Earlier vs Later Treated"      — early cohort vs not-yet-treated late (CLEAN)
  #   "Later vs Earlier Treated"      — late cohort vs already-treated early (FORBIDDEN)
  gb_out  <- bacon(Y ~ D, data = dt, id_var = "id", time_var = "time")
  gb_est  <- sum(gb_out$estimate * gb_out$weight)  # weighted average = TWFE exactly

  # --- GB_Clean: drop forbidden comparisons, renormalize weights -------------
  # "Later vs Earlier Treated" uses an already-treated cohort as control —
  # the one "forbidden" type. Dropping it and renormalizing tests whether
  # forbidden comparisons contain any information beyond the clean ones.
  gb_clean     <- gb_out[gb_out$type != "Later vs Earlier Treated", ]
  w_clean_norm <- gb_clean$weight / sum(gb_clean$weight)  # renormalize to sum = 1
  gb_clean_est <- sum(gb_clean$estimate * w_clean_norm)

  results[s, `:=`(TWFE = twfe_est, GB = gb_est, GB_Clean = gb_clean_est)]
}

# ===========================================================================
# 4. Summary statistics
# ===========================================================================

cat("\n=== Simulation Results (n_sims =", n_sims, ", tau_true =", tau_true, ") ===\n\n")

cat(sprintf("%-25s %10s %10s %12s\n", "Estimator", "Mean", "Bias", "Std Dev"))
cat(strrep("-", 60), "\n")

for (col in c("TWFE", "GB", "GB_Clean")) {
  vals <- results[[col]]
  cat(sprintf("%-25s %10.6f %10.6f %12.6f\n",
              col,
              mean(vals),
              mean(vals) - tau_true,
              sd(vals)))
}

cat("\n=== Numerical Equivalence: TWFE vs GB (all DiDs) ===\n")
diffs_gb <- abs(results$TWFE - results$GB)
cat(sprintf("Max |TWFE - GB| across %d replications: %.2e\n", n_sims, max(diffs_gb)))
cat(sprintf("Mean |TWFE - GB|:                        %.2e\n", mean(diffs_gb)))
cat(sprintf("All differences < 1e-8: %s\n", ifelse(all(diffs_gb < 1e-8), "YES", "NO")))
cat("(Exact to machine precision — GB is an algebraic identity for TWFE)\n")

cat("\n=== Bias comparison: TWFE vs GB_Clean (clean DiDs only) ===\n")
diffs_clean <- results$TWFE - results$GB_Clean
cat(sprintf("Mean (TWFE - GB_Clean):   %8.6f  (should be ~0 under homogeneous effects)\n",
            mean(diffs_clean)))
cat(sprintf("Std (TWFE - GB_Clean):    %8.6f\n", sd(diffs_clean)))
cat(sprintf("Max |TWFE - GB_Clean|:    %8.6f\n", max(abs(diffs_clean))))
cat("(GB_Clean is unbiased but not algebraically identical; weights differ)\n")

cat("\n=== Per-Replication Detail (first 10 sims) ===\n")
cat(sprintf("%-6s %10s %10s %10s %12s %12s\n",
            "Sim", "TWFE", "GB", "GB_Clean", "|TWFE-GB|", "|TWFE-Cln|"))
cat(strrep("-", 65), "\n")
for (s in seq_len(min(10, n_sims))) {
  cat(sprintf("%-6d %10.6f %10.6f %10.6f %12.2e %12.6f\n",
              s,
              results$TWFE[s],
              results$GB[s],
              results$GB_Clean[s],
              abs(results$TWFE[s] - results$GB[s]),
              abs(results$TWFE[s] - results$GB_Clean[s])))
}
