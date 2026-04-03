###############################################################################
# Simulation with AR(1) Errors — Always-Treated Cohort Extension
#
# Based on simulation_ar1.R. Adds one always-treated cohort (g=1, treated
# from the first sample period) with 50 observations. Runs TWFE, CS, SA,
# Gardner, and GMM-Eff with 100 simulation iterations.
#
# True ATT is defined over within-sample cohorts only (g ∈ {10,13,16,19,22}).
# The always-treated cohort is present in the data and affects estimators via:
#   - TWFE    : contaminated time fixed effects (D=1 always → zero within-unit
#               variation in D, but outcomes inflate time FEs for other cohorts)
#   - CS      : always-treated skipped (no pre-period); uses never-treated only
#   - SA      : cohort g=1 has only post-treatment event times; sunab may drop
#               it or distort interaction weights
#   - Gardner : always-treated excluded from first stage (D=1 always) → robust
#   - GMM-Eff : g=1 not in treatment_times → excluded from all moment conditions
#               by construction; GMM targets within-sample cohorts only
###############################################################################

library(data.table)
library(fixest)
library(did)
library(did2s)
library(Matrix)
library(MASS)

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

g_at        <- 1L   # always-treated: treated from period 1 (all sample obs)
n_at        <- 50   # size of always-treated cohort
beta_at_val <- -5   # treatment effect for always-treated cohort

N_within <- cohort_size * n_cohorts + n_never   # 51
N_total  <- N_within + n_at                     # 101

rho     <- 0.5
sigma_u <- 1.0

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
  # g_at=1 > 0 and time_id >= 1 always → D=1 for all always-treated obs
  D_vec <- as.integer(g_vec > 0L & time_id >= g_vec)

  tau_vec <- numeric(N_total * T_total)

  # Within-sample cohorts (g ∈ {10, 13, 16, 19, 22})
  for (c_idx in 1:n_cohorts) {
    g_c  <- treatment_times[c_idx]
    mask <- (g_vec == g_c) & (D_vec == 1L)
    tau_vec[mask] <- beta_g_vec[c_idx] *
                     (1 + r_g_vec[c_idx])^(time_id[mask] - g_c)
  }

  # Always-treated cohort (g=1)
  mask_at <- (g_vec == g_at) & (D_vec == 1L)
  tau_vec[mask_at] <- beta_at * (1 + r_at_val)^(time_id[mask_at] - g_at)

  Y_vec <- alpha[unit_id] + lambda[time_id] + tau_vec + eps

  dt <- data.table(unit = unit_id, time = time_id, Y = Y_vec,
                   D = D_vec, g = g_vec)
  dt[, g_inf := fifelse(g == 0L, Inf,  as.numeric(g))]
  dt[, g_cs  := fifelse(g == 0L, 0L,   as.integer(g))]
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
# 4. GMM system: Q_H matrix and Delta vector
#
# Moments used (g=1 is NOT in treatment_times → excluded by construction):
#   (a) never-treated (g=0) as control
#   (b) not-yet-treated cohorts as control
#   (c) already-treated cohorts from {10,13,16,19,22} with bias correction
# =============================================================================

build_gmm_system <- function(dt) {
  cohort_means <- dt[, .(Y_mean = mean(Y)), by = .(g, time)]
  setkey(cohort_means, g, time)
  get_mean <- function(g_val, t_val) {
    res <- cohort_means[.(g_val, t_val), Y_mean]
    if (length(res) == 0 || is.na(res)) return(NA_real_)
    res
  }

  treated_g <- sort(treatment_times)   # {10,13,16,19,22}; g=1 excluded
  catt_list <- list()
  for (g_c in treated_g)
    for (k in 0:(T_total - g_c))
      catt_list[[length(catt_list) + 1]] <- c(g_c, g_c + k)
  n_catt <- length(catt_list)

  did_estimates <- list()

  for (catt_idx in 1:n_catt) {
    g_c    <- catt_list[[catt_idx]][1]
    t_post <- catt_list[[catt_idx]][2]
    k      <- t_post - g_c

    for (m in 1:(g_c - 1)) {
      t_pre <- g_c - m

      # (a) Never-treated (g=0) as control
      vals <- c(get_mean(g_c, t_post), get_mean(g_c, t_pre),
                get_mean(0,   t_post), get_mean(0,   t_pre))
      if (!anyNA(vals))
        did_estimates[[length(did_estimates) + 1]] <- list(
          delta    = (vals[1] - vals[2]) - (vals[3] - vals[4]),
          catt_idx = catt_idx, type = "never",
          focal_g  = g_c, ctrl_g = 0, t_post = t_post, t_pre = t_pre
        )

      # (b) Not-yet-treated cohorts as control
      for (g_l in treated_g) {
        if (g_l - g_c <= k) next
        if (g_l <= t_post)  next
        vals <- c(get_mean(g_c, t_post), get_mean(g_c, t_pre),
                  get_mean(g_l, t_post), get_mean(g_l, t_pre))
        if (!anyNA(vals))
          did_estimates[[length(did_estimates) + 1]] <- list(
            delta    = (vals[1] - vals[2]) - (vals[3] - vals[4]),
            catt_idx = catt_idx, type = "notyet",
            focal_g  = g_c, ctrl_g = g_l, t_post = t_post, t_pre = t_pre
          )
      }

      # (c) Already-treated cohorts from {10,13,16,19,22} with bias correction
      #     (g=1 never enters here as it is not in treated_g)
      for (g_j in treated_g) {
        j <- g_c - g_j
        if (j <= m)    next
        if (g_j >= g_c) next
        vals <- c(get_mean(g_c, t_post), get_mean(g_c, t_pre),
                  get_mean(g_j, t_post), get_mean(g_j, t_pre))
        if (!anyNA(vals)) {
          bias_neg <- NULL; bias_pos <- NULL
          for (ci in 1:n_catt) {
            if (catt_list[[ci]][1] == g_j && catt_list[[ci]][2] == t_post) bias_neg <- ci
            if (catt_list[[ci]][1] == g_j && catt_list[[ci]][2] == t_pre)  bias_pos <- ci
          }
          did_estimates[[length(did_estimates) + 1]] <- list(
            delta    = (vals[1] - vals[2]) - (vals[3] - vals[4]),
            catt_idx = catt_idx, type = "already",
            bias_neg = bias_neg, bias_pos = bias_pos,
            focal_g  = g_c, ctrl_g = g_j, t_post = t_post, t_pre = t_pre
          )
        }
      }
    }
  }

  n_did <- length(did_estimates)
  if (n_did == 0) return(NULL)

  Delta <- numeric(n_did)
  Q_H   <- matrix(0, nrow = n_did, ncol = n_catt)

  for (s in 1:n_did) {
    est        <- did_estimates[[s]]
    Delta[s]   <- est$delta
    Q_H[s, est$catt_idx] <- 1
    if (est$type == "already") {
      if (!is.null(est$bias_neg)) Q_H[s, est$bias_neg] <- Q_H[s, est$bias_neg] - 1
      if (!is.null(est$bias_pos)) Q_H[s, est$bias_pos] <- Q_H[s, est$bias_pos] + 1
    }
  }

  list(Delta = Delta, Q_H = Q_H, n_did = n_did, n_catt = n_catt,
       catt_list = catt_list, did_estimates = did_estimates)
}

# =============================================================================
# 5. R matrix (maps unit-level observations to DiD estimates)
#    Always-treated units (unit_cohort == g_at) never appear as focal or
#    control cohort, so their columns in R are identically zero.
# =============================================================================

build_R_matrix <- function(gmm_sys) {
  n_did <- gmm_sys$n_did
  NT    <- N_total * T_total   # 101 * 33 = 3333

  R_i <- c(); R_j <- c(); R_x <- c()

  for (s in 1:n_did) {
    est      <- gmm_sys$did_estimates[[s]]
    focal_g  <- est$focal_g
    ctrl_g   <- est$ctrl_g
    t_post   <- est$t_post
    t_pre    <- est$t_pre

    focal_units <- which(unit_cohort == focal_g)
    ctrl_units  <- which(unit_cohort == ctrl_g)
    N_focal <- length(focal_units)
    N_ctrl  <- length(ctrl_units)

    for (u in focal_units) {
      R_i <- c(R_i, s, s)
      R_j <- c(R_j, (u - 1L) * T_total + t_post,
                     (u - 1L) * T_total + t_pre)
      R_x <- c(R_x,  1/N_focal, -1/N_focal)
    }
    for (u in ctrl_units) {
      R_i <- c(R_i, s, s)
      R_j <- c(R_j, (u - 1L) * T_total + t_post,
                     (u - 1L) * T_total + t_pre)
      R_x <- c(R_x, -1/N_ctrl,  1/N_ctrl)
    }
  }

  sparseMatrix(i = R_i, j = R_j, x = R_x, dims = c(n_did, NT))
}

# =============================================================================
# 6. Unit-level Omega_v (block-diagonal; blocks for always-treated units are
#    computed but never accessed by R since those units are not in moments)
# =============================================================================

estimate_omega_v <- function(dt) {
  dt[, resid := residuals(feols(Y ~ D | unit + time, data = dt))]
  lapply(1:N_total, function(i) {
    v <- dt[unit == i, resid]
    v %o% v
  })
}

# =============================================================================
# 7. Var[Delta] = R * Omega_v * R'
# =============================================================================

compute_var_delta <- function(R_mat, omega_blocks) {
  Omega_v   <- bdiag(lapply(omega_blocks, function(b) Matrix(b, sparse = TRUE)))
  R_Om      <- R_mat %*% Omega_v
  as.matrix(R_Om %*% t(R_mat))
}

# =============================================================================
# 8. Estimators
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
    dt_g[, first_treat := fifelse(g == 0L, Inf, as.numeric(g))]
    as.numeric(coef(
      did2s(as.data.frame(dt_g), yname = "Y",
            first_stage  = ~ 0 | unit + time,
            second_stage = ~ i(D, ref = 0),
            treatment    = "D", cluster_var = "unit", verbose = FALSE)
    )["D::1"])
  }, error = function(e) NA_real_)
}

# Efficient GMM: iterates Var[Delta] = R Omega_v R', up to 10 iterations.
# ATT = cohort-size-weighted average of CATTs for cohorts {10,13,16,19,22}.
# g=1 excluded from moment conditions → always-treated cohort does not
# contaminate the GMM estimates.
estimate_gmm_eff <- function(dt, gmm_sys, R_mat) {
  tryCatch({
    Q_H   <- gmm_sys$Q_H
    Delta <- gmm_sys$Delta
    n_catt <- gmm_sys$n_catt
    n_did  <- gmm_sys$n_did

    # Initial OLS estimate
    QtQ     <- crossprod(Q_H)
    beta    <- tryCatch(solve(QtQ, crossprod(Q_H, Delta)),
                        error = function(e) ginv(QtQ) %*% crossprod(Q_H, Delta))

    for (iter in 1:10) {
      beta_old     <- beta
      omega_blocks <- estimate_omega_v(dt)
      Var_Delta    <- compute_var_delta(R_mat, omega_blocks)

      V_inv <- tryCatch(solve(Var_Delta + diag(1e-6, n_did)),
                        error = function(e) NULL)
      if (is.null(V_inv)) break

      QtVQ <- crossprod(Q_H, V_inv %*% Q_H)
      QtVD <- crossprod(Q_H, V_inv %*% Delta)
      if (rcond(QtVQ) < 1e-15) break

      beta <- tryCatch(solve(QtVQ, QtVD), error = function(e) beta_old)
      if (max(abs(beta - beta_old)) < 1e-8) break
    }

    # Cohort-size-weighted ATT
    weights <- rep(cohort_size, n_catt) / (cohort_size * n_catt)
    as.numeric(sum(weights * beta))
  }, error = function(e) NA_real_)
}

# =============================================================================
# 9. Run simulation
# =============================================================================

run_simulation <- function(beta_g_vec, r_g_vec, beta_at, r_at_val, label) {

  true_att <- compute_true_att(beta_g_vec, r_g_vec)
  cat(sprintf("\n=== %s (rho=%.2f) ===\n", label, rho))
  cat(sprintf("True ATT (within-sample cohorts): %.4f\n", true_att))

  results <- data.table(sim     = integer(),
                        TWFE    = numeric(),
                        CS      = numeric(),
                        SA      = numeric(),
                        Gardner = numeric(),
                        GMM_Eff = numeric())

  for (s in 1:n_sims) {
    if (s %% 25 == 0) cat(sprintf("  sim %d / %d\n", s, n_sims))

    dt      <- generate_data_ar1(beta_g_vec, r_g_vec, beta_at, r_at_val)
    gmm_sys <- build_gmm_system(dt)
    R_mat   <- if (!is.null(gmm_sys)) build_R_matrix(gmm_sys) else NULL

    results <- rbindlist(list(results, data.table(
      sim     = s,
      TWFE    = estimate_twfe(dt),
      CS      = estimate_cs(dt),
      SA      = estimate_sa(dt),
      Gardner = estimate_gardner(dt),
      GMM_Eff = if (!is.null(gmm_sys)) estimate_gmm_eff(dt, gmm_sys, R_mat)
                else NA_real_
    )))
  }

  ests <- c("TWFE", "CS", "SA", "Gardner", "GMM_Eff")
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
# 10. Execute
# =============================================================================

cat("=================================================================\n")
cat(sprintf("  AR(1) SIMULATION — ALWAYS-TREATED COHORT (n_at=%d, g=1)\n", n_at))
cat(sprintf("  rho=%.2f  |  n_sims=%d  |  beta_at=%.1f\n",
            rho, n_sims, beta_at_val))
cat("  Estimators: TWFE, CS, SA, Gardner, GMM-Eff\n")
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
# 11. Combined table
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
table_out <- table_out[match(c("TWFE", "CS", "SA", "Gardner", "GMM_Eff"),
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
