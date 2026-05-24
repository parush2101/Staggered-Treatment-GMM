###############################################################################
# Pooled TWFE vs Goodman-Bacon Decomposition — Homogeneous Effects DGP
#
# DGP: Y_it = alpha_i + lambda_t + tau * D_it + epsilon_it
#   where tau = 1 (homogeneous, constant treatment effect)
#         D_it = 1 if unit i is treated and t >= cohort_g_i
#         epsilon_it ~ iid N(0, 1)
#
# Purpose:
#   Under homogeneous effects, TWFE is consistent and equals the Goodman-Bacon
#   weighted average of ALL 2x2 DiDs (algebraic identity, exact to machine
#   precision in every replication).
#
#   The third estimator (GB_Clean) uses only the clean 2x2 DiDs from the
#   Goodman-Bacon decomposition — those whose control group is never-treated
#   ("Treated vs Untreated") or not-yet-treated ("Earlier vs Later Treated").
#   Forbidden DiDs ("Later vs Earlier Treated", where an already-treated cohort
#   serves as control) are dropped.
#
#   Weight redistribution (corrected):
#   Simply renormalizing the remaining clean weights gives the wrong estimator
#   because it misallocates the weight that was on forbidden DiDs. The correct
#   redistribution uses the identity:
#
#     E[Phi_j] = E[Delta_{l,nt}]
#
#   where Phi_j is the forbidden DiD (g_l vs already-treated g_e) and
#   Delta_{l,nt} is the clean "Treated vs Untreated" DiD for the same focal
#   cohort g_l. Under homogeneous effects, g_e's treatment effect cancels
#   completely in Phi_j (g_e is treated in BOTH the pre and post periods of
#   the comparison window), so the forbidden DiD is informationally equivalent
#   to Delta_{l,nt} in expectation. Redistributing weight v_j to Delta_{l,nt}
#   preserves the total weight = 1 and gives E[GB_Clean] = E[TWFE] = tau.
#
#   In finite samples, GB_Clean differs from TWFE by the noise term
#   Sigma_j v_j * (Delta_{l,nt} - Phi_j), which has mean zero. Both
#   estimators are unbiased with comparable variance.
#
# Estimators:
#   1. TWFE:     pooled two-way FE (feols with unit + time FE)
#   2. GB:       Goodman-Bacon weighted average — ALL 2x2 DiDs (exact = TWFE)
#   3. GB_Clean: Goodman-Bacon weighted average — clean DiDs only, with
#                forbidden weights redistributed to the focal cohort's
#                "Treated vs Untreated" DiD (corrected weights, sum = 1)
#
# Simulation: 100 iterations, check:
#   - Bias of each estimator relative to true tau = 1
#   - Variance across replications
#   - TWFE == GB to machine precision in every replication
#   - GB_Clean recovers tau with ~zero bias despite discarding forbidden DiDs
###############################################################################

library(data.table)
library(fixest)
library(bacondecomp)   # install.packages("bacondecomp") if needed

set.seed(42)

# ===========================================================================
# 1. Parameters
# ===========================================================================

tau_true        <- 1.0
cohort_size     <- 20
n_never         <- 20
treatment_times <- c(5, 8, 11, 14)
n_cohorts       <- length(treatment_times)
T_total         <- 20
N_total         <- n_cohorts * cohort_size + n_never
n_sims          <- 100

unit_cohort <- c(rep(treatment_times, each = cohort_size), rep(0L, n_never))
unit_ids    <- seq_len(N_total)
time_ids    <- seq_len(T_total)

# ===========================================================================
# 2. Data generation
# ===========================================================================

generate_data <- function() {
  alpha_i <- rnorm(N_total)
  lambda_t <- rnorm(T_total)
  eps      <- rnorm(N_total * T_total)

  dt <- data.table(
    id   = rep(unit_ids,  each  = T_total),
    time = rep(time_ids,  times = N_total),
    g    = rep(unit_cohort, each = T_total)
  )
  dt[, D := as.integer(g > 0 & time >= g)]
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
  fit_twfe <- feols(Y ~ D | id + time, data = dt, warn = FALSE, notes = FALSE)
  twfe_est <- coef(fit_twfe)["D"]

  # --- Goodman-Bacon decomposition -------------------------------------------
  # bacon() columns: treated, untreated, estimate, weight, type
  # Types:
  #   "Treated vs Untreated"    — any cohort vs never-treated          (CLEAN)
  #   "Earlier vs Later Treated"— early cohort vs not-yet-treated late (CLEAN)
  #   "Later vs Earlier Treated"— late cohort vs already-treated early (FORBIDDEN)
  gb_out <- suppressMessages(
    bacon(Y ~ D, data = dt, id_var = "id", time_var = "time")
  )
  gb_est <- sum(gb_out$estimate * gb_out$weight)   # = TWFE exactly

  # --- GB_Clean: corrected weight redistribution -----------------------------
  # For each forbidden DiD j (g_l vs g_e, weight v_j):
  #   E[Phi_j] = E[Delta_{l,nt}] because g_e is already treated in both
  #   the pre- and post-period of g_l's comparison window, so g_e's treatment
  #   effect cancels and the forbidden DiD reduces to (g_l vs never-treated)
  #   plus mean-zero noise. Redistribute v_j onto Delta_{l,nt} where
  #   l = gb_out$treated for that forbidden row.
  gb_corr <- gb_out
  forbidden_idx <- which(gb_corr$type == "Later vs Earlier Treated")
  for (r in forbidden_idx) {
    g_l   <- gb_corr$treated[r]    # focal later cohort
    v_j   <- gb_corr$weight[r]
    nt_r  <- which(gb_corr$type == "Treated vs Untreated" &
                   gb_corr$treated == g_l)
    gb_corr$weight[nt_r] <- gb_corr$weight[nt_r] + v_j   # redistribute
  }
  # Drop forbidden rows (their weight is now zero by redistribution)
  gb_clean_rows <- gb_corr[gb_corr$type != "Later vs Earlier Treated", ]
  gb_clean_est  <- sum(gb_clean_rows$estimate * gb_clean_rows$weight)
  # Weights already sum to 1 — no renormalization needed

  results[s, `:=`(TWFE = twfe_est, GB = gb_est, GB_Clean = gb_clean_est)]
}

# ===========================================================================
# 4. Summary statistics
# ===========================================================================

cat("\n=== Simulation Results (n_sims =", n_sims, ", tau_true =", tau_true, ") ===\n\n")

cat(sprintf("%-20s %10s %10s %12s\n", "Estimator", "Mean", "Bias", "Std Dev"))
cat(strrep("-", 55), "\n")
for (col in c("TWFE", "GB", "GB_Clean")) {
  vals <- results[[col]]
  cat(sprintf("%-20s %10.6f %10.6f %12.6f\n",
              col, mean(vals), mean(vals) - tau_true, sd(vals)))
}

cat("\n=== Numerical Equivalence: TWFE vs GB (all DiDs) ===\n")
diffs_gb <- abs(results$TWFE - results$GB)
cat(sprintf("Max |TWFE - GB|:   %.2e  (algebraic identity — exact to machine precision)\n",
            max(diffs_gb)))
cat(sprintf("All < 1e-8: %s\n", ifelse(all(diffs_gb < 1e-8), "YES", "NO")))

cat("\n=== Bias comparison: TWFE vs GB_Clean (clean DiDs, corrected weights) ===\n")
diffs_clean <- results$TWFE - results$GB_Clean
cat(sprintf("Mean (TWFE - GB_Clean):  %8.6f  (should be ~0; both unbiased for tau)\n",
            mean(diffs_clean)))
cat(sprintf("Std  (TWFE - GB_Clean):  %8.6f\n", sd(diffs_clean)))
cat(sprintf("Max |TWFE - GB_Clean|:   %8.6f\n", max(abs(diffs_clean))))
cat("Note: GB_Clean differs from TWFE by Sigma_j v_j*(Delta_l_nt - Phi_j),\n")
cat("      a mean-zero noise term; the two are NOT algebraically identical.\n")

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
