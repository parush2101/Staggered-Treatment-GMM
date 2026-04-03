###############################################################################
# Simulation with AR(1) Errors — Always-Treated Cohort Extension
#
# Based on simulation_ar1.R. Adds one always-treated cohort (g=1, treated
# from the first sample period) with 50 observations. Runs TWFE, CS, SA, and
# Gardner only (no GMM). 100 simulation iterations.
#
# True ATT is defined over within-sample cohorts only (g ∈ {10,13,16,19,22}).
# The always-treated cohort is present in the data and affects estimators via:
#   - TWFE  : contaminated time fixed effects (D=1 always → zero within-unit
#             variation in D, but outcomes inflate time FEs for other cohorts)
#   - CS    : always-treated never appear as not-yet-treated controls; CS skips
#             cohort g=1 (no pre-treatment period) and uses never-treated only
#   - SA    : cohort g=1 has only post-treatment event times; sunab may drop it
#             or produce distorted interaction weights
#   - Gardner: always-treated never appear in the first-stage (D=1 always);
#              excluded from FE estimation → Gardner unaffected
###############################################################################

library(data.table)
library(fixest)
library(did)
library(did2s)

set.seed(42)

# =============================================================================
# 1. Parameters
# =============================================================================

T_total         <- 33
n_sims          <- 100
cohort_size     <- 8
n_never         <- 11
treatment_times <- c(10, 13, 16, 19, 22)
n_cohorts       <- length(treatment_times)

g_at            <- 1L    # always-treated: treated from period 1 (all sample obs)
n_at            <- 50    # size of always-treated cohort
beta_at_val     <- -5    # treatment effect for always-treated cohort

N_within        <- cohort_size * n_cohorts + n_never   # 51
N_total         <- N_within + n_at                     # 101

rho             <- 0.5
sigma_u         <- 1.0

# Unit cohort vector: units 1-40 within-sample, 41-51 never-treated, 52-101 always-treated
unit_cohort <- c(rep(treatment_times, each = cohort_size),
                 rep(0L,   n_never),
                 rep(g_at, n_at))

# =============================================================================
# 2. DGP with AR(1) errors
# =============================================================================

generate_data_ar1 <- function(beta_g_vec, r_g_vec, beta_at, r_at_val) {

  unit_id <- rep(1:N_total, each = T_total)
  time_id <- rep(1:T_total, times = N_total)

  alpha  <- rnorm(N_total)
  lambda <- rnorm(T_total)

  eps <- numeric(N_total * T_total)
  for (i in 1:N_total) {
    idx <- (i - 1L) * T_total + 1L
    eps[idx] <- rnorm(1, 0, sigma_u / sqrt(1 - rho^2))
    for (t in 2:T_total)
      eps[idx + t - 1L] <- rho * eps[idx + t - 2L] + rnorm(1, 0, sigma_u)
  }

  g_vec <- unit_cohort[unit_id]
  # g_at = 1 > 0 and time_id >= 1 always → D = 1 for all always-treated obs
  D_vec <- as.integer(g_vec > 0L & time_id >= g_vec)

  tau_vec <- numeric(N_total * T_total)

  # Within-sample cohorts (g ∈ {10, 13, 16, 19, 22})
  for (c_idx in 1:n_cohorts) {
    g_c  <- treatment_times[c_idx]
    mask <- (g_vec == g_c) & (D_vec == 1L)
    tau_vec[mask] <- beta_g_vec[c_idx] *
                     (1 + r_g_vec[c_idx])^(time_id[mask] - g_c)
  }

  # Always-treated cohort (g = 1)
  mask_at <- (g_vec == g_at) & (D_vec == 1L)
  tau_vec[mask_at] <- beta_at * (1 + r_at_val)^(time_id[mask_at] - g_at)

  Y_vec <- alpha[unit_id] + lambda[time_id] + tau_vec + eps

  dt <- data.table(unit = unit_id, time = time_id, Y = Y_vec,
                   D = D_vec, g = g_vec)
  # g_inf = Inf for never-treated (g=0); actual g for all others incl. g=1
  dt[, g_inf := fifelse(g == 0L, Inf, as.numeric(g))]
  # g_cs  = 0  for never-treated; actual g for all others incl. g=1
  dt[, g_cs  := fifelse(g == 0L, 0L, as.integer(g))]
  return(dt)
}

# =============================================================================
# 3. True ATT (within-sample cohorts only)
# =============================================================================

compute_true_att <- function(beta_g_vec, r_g_vec) {
  total_te <- 0; total_obs <- 0
  for (c_idx in 1:n_cohorts) {
    g_c <- treatment_times[c_idx]
    for (t in g_c:T_total) {
      te        <- beta_g_vec[c_idx] * (1 + r_g_vec[c_idx])^(t - g_c)
      total_te  <- total_te  + cohort_size * te
      total_obs <- total_obs + cohort_size
    }
  }
  return(total_te / total_obs)
}

# =============================================================================
# 4. Estimators
# =============================================================================

estimate_twfe <- function(dt) {
  tryCatch(
    as.numeric(coef(feols(Y ~ D | unit + time, data = dt))["D"]),
    error = function(e) NA_real_
  )
}

estimate_cs <- function(dt) {
  tryCatch({
    out <- att_gt(yname = "Y", tname = "time", idname = "unit", gname = "g_cs",
                  data = as.data.frame(dt), control_group = "nevertreated",
                  print_details = FALSE, bstrap = FALSE, cband = FALSE)
    aggte(out, type = "simple", na.rm = TRUE)$overall.att
  }, error = function(e) NA_real_)
}

estimate_sa <- function(dt) {
  tryCatch({
    mod <- feols(Y ~ sunab(g_inf, time) | unit + time, data = dt)
    summary(mod, agg = "ATT")$coeftable[1, 1]
  }, error = function(e) NA_real_)
}

estimate_gardner <- function(dt) {
  tryCatch({
    dt_g <- copy(dt)
    # did2s first stage uses D=0 obs only → always-treated (D=1 always) excluded
    dt_g[, first_treat := fifelse(g == 0L, Inf, as.numeric(g))]
    as.numeric(coef(
      did2s(as.data.frame(dt_g), yname = "Y",
            first_stage  = ~ 0 | unit + time,
            second_stage = ~ i(D, ref = 0),
            treatment    = "D", cluster_var = "unit", verbose = FALSE)
    )["D::1"])
  }, error = function(e) NA_real_)
}

# =============================================================================
# 5. Run simulation
# =============================================================================

run_simulation <- function(beta_g_vec, r_g_vec, beta_at, r_at_val, label) {

  true_att <- compute_true_att(beta_g_vec, r_g_vec)
  cat(sprintf("\n=== %s (rho=%.2f) ===\n", label, rho))
  cat(sprintf("True ATT (within-sample cohorts): %.4f\n", true_att))

  results <- data.table(sim     = integer(),
                        TWFE    = numeric(),
                        CS      = numeric(),
                        SA      = numeric(),
                        Gardner = numeric())

  for (s in 1:n_sims) {
    if (s %% 25 == 0) cat(sprintf("  sim %d / %d\n", s, n_sims))
    dt <- generate_data_ar1(beta_g_vec, r_g_vec, beta_at, r_at_val)
    results <- rbindlist(list(results, data.table(
      sim     = s,
      TWFE    = estimate_twfe(dt),
      CS      = estimate_cs(dt),
      SA      = estimate_sa(dt),
      Gardner = estimate_gardner(dt)
    )))
  }

  ests <- c("TWFE", "CS", "SA", "Gardner")
  summary_dt <- data.table(
    Estimator = ests,
    Bias = sapply(ests, function(e) {
      v <- results[[e]][!is.na(results[[e]])]
      round(mean(v) - true_att, 4)
    }),
    Variance = sapply(ests, function(e) {
      v <- results[[e]][!is.na(results[[e]])]
      round(var(v), 4)
    }),
    RMSE = sapply(ests, function(e) {
      v <- results[[e]][!is.na(results[[e]])]
      round(sqrt(mean((v - true_att)^2)), 4)
    }),
    N_valid = sapply(ests, function(e) sum(!is.na(results[[e]])))
  )

  cat(sprintf("\nResults (%s):\n", label))
  print(summary_dt)
  return(list(results = results, summary = summary_dt, true_att = true_att))
}

# =============================================================================
# 6. Execute
# =============================================================================

cat("=================================================================\n")
cat(sprintf("  AR(1) SIMULATION — ALWAYS-TREATED COHORT (n_at=%d, g=1)\n", n_at))
cat(sprintf("  rho=%.2f  |  n_sims=%d  |  beta_at=%.1f\n",
            rho, n_sims, beta_at_val))
cat("  Estimators: TWFE, CS, SA, Gardner\n")
cat("  True ATT: within-sample cohorts (g ∈ {10,13,16,19,22}) only\n")
cat("=================================================================\n")

# --- Homogeneous ---
res_hom <- run_simulation(
  beta_g_vec = rep(-5, n_cohorts),
  r_g_vec    = rep(0,  n_cohorts),
  beta_at    = beta_at_val,
  r_at_val   = 0,
  label      = "Homogeneous"
)

# --- Heterogeneous ---
res_het <- run_simulation(
  beta_g_vec = c(-16, -12, -10, -9, -2),
  r_g_vec    = c(0.01, 0.04, 0.08, 0.10, 0.07),
  beta_at    = beta_at_val,
  r_at_val   = 0,
  label      = "Heterogeneous"
)

# =============================================================================
# 7. Combined table
# =============================================================================

cat("\n\n")
cat("=================================================================\n")
cat(sprintf("  TABLE: Bias, Variance, RMSE — Always-Treated (n_at=%d)\n", n_at))
cat(sprintf("  AR(1) errors (rho=%.2f), %d simulations\n", rho, n_sims))
cat("=================================================================\n\n")

table_out <- merge(
  res_hom$summary[, .(Estimator,
                      Hom_Bias = Bias, Hom_Var = Variance, Hom_RMSE = RMSE)],
  res_het$summary[, .(Estimator,
                      Het_Bias = Bias, Het_Var = Variance, Het_RMSE = RMSE)],
  by = "Estimator"
)
table_out <- table_out[match(c("TWFE", "CS", "SA", "Gardner"),
                              table_out$Estimator)]

cat(sprintf("%-10s  %9s %9s %9s    %9s %9s %9s\n",
            "Estimator",
            "Hom.Bias", "Hom.Var", "Hom.RMSE",
            "Het.Bias", "Het.Var", "Het.RMSE"))
cat(paste(rep("-", 72), collapse = ""), "\n")
for (i in 1:nrow(table_out)) {
  cat(sprintf("%-10s  %9.4f %9.4f %9.4f    %9.4f %9.4f %9.4f\n",
              table_out$Estimator[i],
              table_out$Hom_Bias[i],  table_out$Hom_Var[i],  table_out$Hom_RMSE[i],
              table_out$Het_Bias[i],  table_out$Het_Var[i],  table_out$Het_RMSE[i]))
}

save(res_hom, res_het, table_out,
     file = "simulation_ar1_always_treated_results.RData")
cat("\nResults saved to simulation_ar1_always_treated_results.RData\n")
