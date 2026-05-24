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
#   The third estimator (Gardner_Clean) is Gardner's two-stage DiD (did2s).
#   It uses only clean comparisons by construction: Stage 1 estimates unit and
#   time fixed effects from untreated/not-yet-treated observations only, so no
#   already-treated unit ever serves as a control. Stage 2 regresses residuals
#   on D using all observations.
#
#   Why not simple GB weight redistribution?
#   The intuition "forbidden DiDs are linear combinations of clean DiDs, so
#   redistribute their weights" works in the fine-grained GMM moment-conditions
#   framework (one moment per cohort × time cell), where the identity holds
#   algebraically. In the Goodman-Bacon 2x2 framework each DiD aggregates across
#   periods with DIFFERENT time windows per pair, so the forbidden 2x2 DiD
#   (g_l vs already-treated g_e) equals (g_l vs never-treated) minus a residual
#   term (g_e vs never-treated in g_l's window) that is NOT in the bacon() output.
#   There are no exact c_{ji} coefficients within the GB 2x2 set.
#
#   Gardner operates at the fine-grained moment level (within-period variation),
#   where redundancy holds exactly. Under homogeneous effects + spherical errors,
#   Gardner is asymptotically BLUE = TWFE (Gauss-Markov). In finite samples the
#   two-stage estimation introduces an O(1/N) deviation from TWFE's one-stage
#   estimate, so |TWFE - Gardner| is small but nonzero in each replication.
#
# Estimators:
#   1. TWFE:          pooled two-way FE (feols with unit + time FE)
#   2. GB:            Goodman-Bacon weighted average — ALL 2x2 DiDs (exact = TWFE)
#   3. Gardner_Clean: two-stage DiD (did2s) — clean comparisons only by construction
#
# Simulation: 100 iterations, track point estimates and check:
#   - Bias of each estimator relative to true tau = 1
#   - Variance across replications
#   - Whether TWFE == GB numerically in every replication
#   - Whether Gardner_Clean recovers tau despite using only clean comparisons
#   - Size of finite-sample gap between TWFE and Gardner_Clean (vanishes as N grows)
###############################################################################

library(data.table)
library(fixest)
library(bacondecomp)   # install.packages("bacondecomp") if needed
library(did2s)         # install.packages("did2s") if needed

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
  sim           = seq_len(n_sims),
  TWFE          = NA_real_,
  GB            = NA_real_,
  Gardner_Clean = NA_real_
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

  # --- Gardner_Clean: two-stage DiD using only clean comparisons -------------
  # Stage 1: regress Y on unit FE + time FE from D=0 (untreated) obs only.
  # Stage 2: regress residuals on D using all observations.
  # No already-treated unit ever enters as control — all comparisons are clean.
  # Converges to TWFE as N→∞ (both are BLUE under homogeneous + spherical errors).
  # Finite-sample gap with TWFE is O(1/N) from two-stage vs one-stage partialling.
  fit_gardner   <- suppressMessages(
    did2s(dt, yname = "Y", first_stage = ~ 0 | id + time,
          second_stage = ~ D, treatment = "D", cluster_var = "id")
  )
  gardner_est <- coef(fit_gardner)["D"]

  results[s, `:=`(TWFE = twfe_est, GB = gb_est, Gardner_Clean = gardner_est)]
}

# ===========================================================================
# 4. Summary statistics
# ===========================================================================

cat("\n=== Simulation Results (n_sims =", n_sims, ", tau_true =", tau_true, ") ===\n\n")

cat(sprintf("%-25s %10s %10s %12s\n", "Estimator", "Mean", "Bias", "Std Dev"))
cat(strrep("-", 60), "\n")

for (col in c("TWFE", "GB", "Gardner_Clean")) {
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

cat("\n=== Bias comparison: TWFE vs Gardner_Clean (clean comparisons only) ===\n")
diffs_gardner <- results$TWFE - results$Gardner_Clean
cat(sprintf("Mean (TWFE - Gardner_Clean):  %8.6f  (should be ~0; both unbiased)\n",
            mean(diffs_gardner)))
cat(sprintf("Std (TWFE - Gardner_Clean):   %8.6f  (finite-sample two-stage gap)\n",
            sd(diffs_gardner)))
cat(sprintf("Max |TWFE - Gardner_Clean|:   %8.6f  (shrinks as N grows)\n",
            max(abs(diffs_gardner))))
cat("(Gardner_Clean is unbiased and converges to TWFE as N->Inf,\n")
cat(" but differs by O(1/N) in finite samples due to two-stage vs one-stage FE)\n")

cat("\n=== Per-Replication Detail (first 10 sims) ===\n")
cat(sprintf("%-6s %10s %10s %12s %12s %12s\n",
            "Sim", "TWFE", "GB", "Gard_Clean", "|TWFE-GB|", "|TWFE-Gard|"))
cat(strrep("-", 70), "\n")
for (s in seq_len(min(10, n_sims))) {
  cat(sprintf("%-6d %10.6f %10.6f %12.6f %12.2e %12.6f\n",
              s,
              results$TWFE[s],
              results$GB[s],
              results$Gardner_Clean[s],
              abs(results$TWFE[s] - results$GB[s]),
              abs(results$TWFE[s] - results$Gardner_Clean[s])))
}
