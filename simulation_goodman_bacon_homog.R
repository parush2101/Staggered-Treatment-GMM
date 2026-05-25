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
#   Exact algebraic identity (key insight):
#   For each forbidden DiD Phi_j = DiD(g_l, g_e, post, mid) with weight v_j,
#   the following holds exactly in every sample (ȳ_U terms cancel):
#
#     DiD(g_l, g_e, post, mid) = DiD(g_l, U, post, mid) - DiD(g_e, U, post, mid)
#
#   where U = never-treated units, mid = [g_e, g_l-1], post = [g_l, T].
#   Both RHS terms use never-treated control (clean DiDs). Substituting for
#   every forbidden DiD and collecting terms gives:
#
#     TWFE = Sigma_{clean} w_i * Delta_i
#          + Sigma_{forbidden} v_j * [DiD(g_l,U,post,mid) - DiD(g_e,U,post,mid)]
#
#   This decomposition uses ONLY clean (never-treated or not-yet-treated) DiDs
#   and recovers TWFE exactly, to machine precision, in every replication.
#   Note: the DiD(g_e, U, post, mid) terms enter with NEGATIVE weight (-v_j),
#   because in [mid, post], g_e is already treated, so its clean DiD captures
#   g_e's own treatment effect (≈ 0 under homogeneous constant effects, but
#   nonzero in finite samples), introducing noise that reduces efficiency.
#
# Estimators:
#   1. TWFE:     pooled two-way FE (feols with unit + time FE)
#   2. GB:       Goodman-Bacon weighted average — ALL 2x2 DiDs (exact = TWFE)
#   3. GB_Clean: uses the algebraic identity DiD(g_l,g_e) = DiD(g_l,U) - DiD(g_e,U)
#                to substitute every forbidden DiD with two clean DiDs computed
#                over the same time window. Recovers TWFE exactly (machine precision).
#
# Simulation: 100 iterations, check:
#   - Bias of each estimator relative to true tau = 1
#   - Variance across replications
#   - TWFE == GB to machine precision in every replication
#   - GB_Clean == TWFE to machine precision in every replication (exact identity)
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

  # --- GB_Clean: exact recovery via algebraic identity -----------------------
  # Start with sum over all clean rows.
  gb_clean_est <- sum(gb_out$estimate[gb_out$type != "Later vs Earlier Treated"] *
                      gb_out$weight[gb_out$type != "Later vs Earlier Treated"])

  # For each forbidden DiD (g_l vs already-treated g_e, weight v_j):
  #   DiD(g_l, g_e, post, mid) = DiD(g_l, U, post, mid) - DiD(g_e, U, post, mid)
  # where mid = [g_e, g_l-1], post = [g_l, T_total].
  # The ȳ_U terms cancel algebraically, so the identity holds to machine precision.
  # Add v_j * DiD(g_l, U) - v_j * DiD(g_e, U) in place of each forbidden DiD.
  for (r in which(gb_out$type == "Later vs Earlier Treated")) {
    g_l_r <- gb_out$treated[r]    # focal later cohort
    g_e_r <- gb_out$untreated[r]  # control earlier cohort (already treated)
    v_j   <- gb_out$weight[r]

    mid_periods  <- g_e_r:(g_l_r - 1)  # g_e treated, g_l not yet treated
    post_periods <- g_l_r:T_total       # both cohorts treated

    mean_gl_post <- mean(dt[g == g_l_r & time %in% post_periods, Y])
    mean_gl_mid  <- mean(dt[g == g_l_r & time %in% mid_periods,  Y])
    mean_ge_post <- mean(dt[g == g_e_r & time %in% post_periods, Y])
    mean_ge_mid  <- mean(dt[g == g_e_r & time %in% mid_periods,  Y])
    mean_U_post  <- mean(dt[g == 0     & time %in% post_periods, Y])
    mean_U_mid   <- mean(dt[g == 0     & time %in% mid_periods,  Y])

    did_gl_U <- (mean_gl_post - mean_gl_mid) - (mean_U_post - mean_U_mid)
    did_ge_U <- (mean_ge_post - mean_ge_mid) - (mean_U_post - mean_U_mid)

    gb_clean_est <- gb_clean_est + v_j * did_gl_U - v_j * did_ge_U
  }
  # Result equals TWFE to machine precision (ȳ_U cancels; same linear combination)

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

cat("\n=== Numerical Equivalence: TWFE vs GB_Clean (clean DiDs, exact identity) ===\n")
diffs_clean <- abs(results$TWFE - results$GB_Clean)
cat(sprintf("Max |TWFE - GB_Clean|:   %.2e  (algebraic identity — exact to machine precision)\n",
            max(diffs_clean)))
cat(sprintf("All < 1e-8: %s\n", ifelse(all(diffs_clean < 1e-8), "YES", "NO")))
cat("Note: GB_Clean substitutes each forbidden DiD via DiD(g_l,g_e) = DiD(g_l,U) - DiD(g_e,U).\n")
cat("      The ȳ_U terms cancel algebraically, so GB_Clean = TWFE in every sample.\n")

cat("\n=== Per-Replication Detail (first 10 sims) ===\n")
cat(sprintf("%-6s %10s %10s %10s %12s %12s\n",
            "Sim", "TWFE", "GB", "GB_Clean", "|TWFE-GB|", "|TWFE-Cln|"))
cat(strrep("-", 65), "\n")
for (s in seq_len(min(10, n_sims))) {
  cat(sprintf("%-6d %10.6f %10.6f %10.6f %12.2e %12.2e\n",
              s,
              results$TWFE[s],
              results$GB[s],
              results$GB_Clean[s],
              abs(results$TWFE[s] - results$GB[s]),
              abs(results$TWFE[s] - results$GB_Clean[s])))
}
