###############################################################################
# Simulation with Dependent Errors (AR(1) within units)
# Compares three GMM weighting approaches:
#   1. GMM with A = I (identity)
#   2. GMM with iterative diagonal-block A (non-zero covariances within unit)
#   3. GMM with iterative GLS using structured Var[Delta] = R * Omega_v * R'
#
# DGP: Y_it = alpha_i + lambda_t + tau_it + epsilon_it
#   where epsilon_it = rho * epsilon_{i,t-1} + u_it, u_it ~ N(0, sigma_u^2)
#   so errors are AR(1) within units, independent across units.
###############################################################################

library(data.table)
library(fixest)
library(did)
library(did2s)
library(MASS)
library(Matrix)

set.seed(42)

# ===========================================================================
# 1. Parameters
# ===========================================================================

N_total         <- 51
T_total         <- 33
n_sims          <- 500
cohort_size     <- 8
n_never         <- 11
treatment_times <- c(10, 13, 16, 19, 22)
n_cohorts       <- length(treatment_times)

# AR(1) error parameters
rho     <- 0.5     # autocorrelation
sigma_u <- 1.0     # innovation sd

unit_cohort <- c(rep(treatment_times, each = cohort_size), rep(0, n_never))

# ===========================================================================
# 2. DGP with AR(1) errors
# ===========================================================================

generate_data_ar1 <- function(beta_g_vec, r_g_vec) {
  unit_id <- rep(1:N_total, each = T_total)
  time_id <- rep(1:T_total, times = N_total)

  alpha <- rnorm(N_total)
  lambda <- rnorm(T_total)

  # Generate AR(1) errors for each unit
  eps <- numeric(N_total * T_total)
  for (i in 1:N_total) {
    idx_start <- (i - 1) * T_total + 1
    eps[idx_start] <- rnorm(1, 0, sigma_u / sqrt(1 - rho^2))  # stationary init
    for (t in 2:T_total) {
      eps[idx_start + t - 1] <- rho * eps[idx_start + t - 2] + rnorm(1, 0, sigma_u)
    }
  }

  g_vec <- unit_cohort[unit_id]
  D_vec <- as.integer(g_vec > 0 & time_id >= g_vec)

  # Treatment effect
  tau_vec <- numeric(N_total * T_total)
  for (c_idx in 1:n_cohorts) {
    g_c <- treatment_times[c_idx]
    mask <- (g_vec == g_c) & (D_vec == 1)
    tau_vec[mask] <- beta_g_vec[c_idx] * (1 + r_g_vec[c_idx])^(time_id[mask] - g_c)
  }

  Y_vec <- alpha[unit_id] + lambda[time_id] + tau_vec + eps

  dt <- data.table(
    unit = unit_id, time = time_id, Y = Y_vec, D = D_vec,
    g = g_vec, eps = eps
  )
  dt[, g_inf := fifelse(g == 0, Inf, as.numeric(g))]
  dt[, g_cs  := fifelse(g == 0, 0L, as.integer(g))]
  return(dt)
}

# ===========================================================================
# 3. True ATT
# ===========================================================================

compute_true_att <- function(beta_g_vec, r_g_vec) {
  total_te <- 0; total_obs <- 0
  for (c_idx in 1:n_cohorts) {
    g_c <- treatment_times[c_idx]
    for (t in g_c:T_total) {
      te <- beta_g_vec[c_idx] * (1 + r_g_vec[c_idx])^(t - g_c)
      total_te <- total_te + cohort_size * te
      total_obs <- total_obs + cohort_size
    }
  }
  return(total_te / total_obs)
}

# ===========================================================================
# 4. Build Q_H and Delta (shared across GMM variants)
# ===========================================================================

build_gmm_system <- function(dt) {
  cohort_means <- dt[, .(Y_mean = mean(Y)), by = .(g, time)]
  setkey(cohort_means, g, time)
  get_mean <- function(g_val, t_val) {
    res <- cohort_means[.(g_val, t_val), Y_mean]
    if (length(res) == 0 || is.na(res)) return(NA_real_)
    return(res)
  }

  treated_g <- sort(treatment_times)
  catt_list <- list()
  for (g_c in treated_g) {
    for (k in 0:(T_total - g_c)) {
      catt_list[[length(catt_list) + 1]] <- c(g_c, g_c + k)
    }
  }
  n_catt <- length(catt_list)

  # Store metadata for R matrix construction
  did_estimates <- list()

  for (catt_idx in 1:n_catt) {
    g_c <- catt_list[[catt_idx]][1]
    t_post <- catt_list[[catt_idx]][2]
    k <- t_post - g_c

    for (m in 1:(g_c - 1)) {
      t_pre <- g_c - m

      # Never-treated
      vals <- c(get_mean(g_c, t_post), get_mean(g_c, t_pre),
                get_mean(0, t_post), get_mean(0, t_pre))
      if (!any(is.na(vals))) {
        did_estimates[[length(did_estimates) + 1]] <- list(
          delta = (vals[1] - vals[2]) - (vals[3] - vals[4]),
          catt_idx = catt_idx, type = "never",
          focal_g = g_c, ctrl_g = 0, t_post = t_post, t_pre = t_pre
        )
      }

      # Not-yet-treated
      for (g_l in treated_g) {
        if (g_l - g_c <= k) next
        if (g_l <= t_post) next
        vals <- c(get_mean(g_c, t_post), get_mean(g_c, t_pre),
                  get_mean(g_l, t_post), get_mean(g_l, t_pre))
        if (!any(is.na(vals))) {
          did_estimates[[length(did_estimates) + 1]] <- list(
            delta = (vals[1] - vals[2]) - (vals[3] - vals[4]),
            catt_idx = catt_idx, type = "notyet",
            focal_g = g_c, ctrl_g = g_l, t_post = t_post, t_pre = t_pre
          )
        }
      }

      # Already-treated (forbidden)
      for (g_j in treated_g) {
        j <- g_c - g_j
        if (j <= m) next
        if (g_j >= g_c) next
        vals <- c(get_mean(g_c, t_post), get_mean(g_c, t_pre),
                  get_mean(g_j, t_post), get_mean(g_j, t_pre))
        if (!any(is.na(vals))) {
          bias_neg <- NULL; bias_pos <- NULL
          for (ci in 1:n_catt) {
            if (catt_list[[ci]][1] == g_j && catt_list[[ci]][2] == t_post) bias_neg <- ci
            if (catt_list[[ci]][1] == g_j && catt_list[[ci]][2] == t_pre)  bias_pos <- ci
          }
          did_estimates[[length(did_estimates) + 1]] <- list(
            delta = (vals[1] - vals[2]) - (vals[3] - vals[4]),
            catt_idx = catt_idx, type = "already",
            bias_neg = bias_neg, bias_pos = bias_pos,
            focal_g = g_c, ctrl_g = g_j, t_post = t_post, t_pre = t_pre
          )
        }
      }
    }
  }

  n_did <- length(did_estimates)
  if (n_did == 0) return(NULL)

  Delta <- numeric(n_did)
  Q_H <- matrix(0, nrow = n_did, ncol = n_catt)

  for (s in 1:n_did) {
    est <- did_estimates[[s]]
    Delta[s] <- est$delta
    Q_H[s, est$catt_idx] <- 1
    if (est$type == "already") {
      if (!is.null(est$bias_neg) && !is.na(est$bias_neg))
        Q_H[s, est$bias_neg] <- Q_H[s, est$bias_neg] - 1
      if (!is.null(est$bias_pos) && !is.na(est$bias_pos))
        Q_H[s, est$bias_pos] <- Q_H[s, est$bias_pos] + 1
    }
  }

  return(list(Delta = Delta, Q_H = Q_H, n_did = n_did, n_catt = n_catt,
              catt_list = catt_list, did_estimates = did_estimates))
}

# ===========================================================================
# 5. Build R matrix (maps unit-level observations to 2x2 DiD estimates)
# ===========================================================================

build_R_matrix <- function(gmm_sys, dt) {
  # R is n_did x (N*T) matrix
  # Each row corresponds to a 2x2 DiD estimate
  # Delta_s = sum over units of (+/- 1/N_g) * Y_{i,t}

  n_did <- gmm_sys$n_did
  NT <- N_total * T_total

  # Sparse construction
  R_i <- c(); R_j <- c(); R_x <- c()

  for (s in 1:n_did) {
    est <- gmm_sys$did_estimates[[s]]
    focal_g <- est$focal_g
    ctrl_g  <- est$ctrl_g
    t_post  <- est$t_post
    t_pre   <- est$t_pre

    # Get units in focal cohort
    focal_units <- which(unit_cohort == focal_g)
    N_focal <- length(focal_units)

    # Get units in control cohort
    ctrl_units <- which(unit_cohort == ctrl_g)
    N_ctrl <- length(ctrl_units)

    # DiD = (1/N_f) sum Y_{i,t_post} - (1/N_f) sum Y_{i,t_pre}
    #      -(1/N_c) sum Y_{j,t_post} + (1/N_c) sum Y_{j,t_pre}

    for (u in focal_units) {
      col_post <- (u - 1) * T_total + t_post
      col_pre  <- (u - 1) * T_total + t_pre
      R_i <- c(R_i, s, s)
      R_j <- c(R_j, col_post, col_pre)
      R_x <- c(R_x, 1/N_focal, -1/N_focal)
    }
    for (u in ctrl_units) {
      col_post <- (u - 1) * T_total + t_post
      col_pre  <- (u - 1) * T_total + t_pre
      R_i <- c(R_i, s, s)
      R_j <- c(R_j, col_post, col_pre)
      R_x <- c(R_x, -1/N_ctrl, 1/N_ctrl)
    }
  }

  R <- sparseMatrix(i = R_i, j = R_j, x = R_x,
                    dims = c(n_did, NT))
  return(R)
}

# ===========================================================================
# 6. Estimate unit-level Omega_v (block diagonal, within-unit covariance)
# ===========================================================================

estimate_omega_v <- function(dt, beta_hat, gmm_sys) {
  # Step 1: Get residuals from TWFE
  mod <- feols(Y ~ D | unit + time, data = dt)
  dt[, resid := residuals(mod)]

  # Step 2: Build block-diagonal Omega_v
  # Each block is T x T covariance matrix for unit i
  # estimated from outer product of residuals

  # Return as list of T x T blocks (for efficiency)
  omega_blocks <- list()
  for (i in 1:N_total) {
    v_i <- dt[unit == i, resid]
    omega_blocks[[i]] <- v_i %o% v_i  # T x T
  }

  return(omega_blocks)
}

# ===========================================================================
# 7. Compute Var[Delta] = R * Omega_v * R' using block structure
# ===========================================================================

compute_var_delta <- function(R_mat, omega_blocks) {
  # Omega_v is block diagonal: block i is T x T for unit i
  # R * Omega_v * R' can be computed block-wise

  n_did <- nrow(R_mat)
  NT <- ncol(R_mat)

  # Build full Omega_v as sparse block diagonal
  block_list <- lapply(omega_blocks, function(b) Matrix(b, sparse = TRUE))
  Omega_v <- bdiag(block_list)

  # Var[Delta] = R %*% Omega_v %*% t(R)
  R_Omega <- R_mat %*% Omega_v
  Var_Delta <- as.matrix(R_Omega %*% t(R_mat))

  return(Var_Delta)
}

# ===========================================================================
# 8. Three GMM estimators
# ===========================================================================

compute_att <- function(beta_hat, n_catt) {
  weights <- rep(cohort_size, n_catt)
  weights <- weights / sum(weights)
  return(as.numeric(sum(weights * beta_hat)))
}

# Method 1: GMM with A = I
gmm_identity <- function(gmm_sys) {
  Q_H <- gmm_sys$Q_H; Delta <- gmm_sys$Delta; n_catt <- gmm_sys$n_catt
  QtQ <- crossprod(Q_H)
  if (rcond(QtQ) < 1e-15) return(NA_real_)
  beta_hat <- solve(QtQ, crossprod(Q_H, Delta))
  return(compute_att(beta_hat, n_catt))
}

# Method 2: Iterative GMM with structured A
# A is based on Var[Delta] with block-diagonal Omega_v
gmm_iterative_structured <- function(gmm_sys, dt) {
  Q_H <- gmm_sys$Q_H; Delta <- gmm_sys$Delta; n_catt <- gmm_sys$n_catt
  n_did <- gmm_sys$n_did

  # Step 1: Initial estimate with A = I
  QtQ <- crossprod(Q_H)
  if (rcond(QtQ) < 1e-15) return(NA_real_)
  beta_hat <- solve(QtQ, crossprod(Q_H, Delta))

  # Build R matrix once
  R_mat <- build_R_matrix(gmm_sys, dt)

  # Step 2: Iterate
  for (iter in 1:5) {
    # Estimate Omega_v from residuals
    omega_blocks <- estimate_omega_v(dt, beta_hat, gmm_sys)

    # Compute Var[Delta] = R Omega_v R'
    Var_Delta <- compute_var_delta(R_mat, omega_blocks)

    # A = Var[Delta]^{-1}
    A <- tryCatch(solve(Var_Delta + diag(1e-6, n_did)),
                  error = function(e) NULL)
    if (is.null(A)) break

    QtAQ <- crossprod(Q_H, A %*% Q_H)
    QtAD <- crossprod(Q_H, A %*% Delta)
    if (rcond(QtAQ) < 1e-15) break
    beta_hat <- solve(QtAQ, QtAD)
  }

  return(compute_att(beta_hat, n_catt))
}

# Method 3: Iterative GLS on moment conditions
# Directly run GLS: beta = (Q_H' Var[Delta]^{-1} Q_H)^{-1} Q_H' Var[Delta]^{-1} Delta
# with Var[Delta] = R * Omega_v * R', iterating Omega_v estimation
gmm_gls <- function(gmm_sys, dt) {
  Q_H <- gmm_sys$Q_H; Delta <- gmm_sys$Delta; n_catt <- gmm_sys$n_catt
  n_did <- gmm_sys$n_did

  # Build R matrix once
  R_mat <- build_R_matrix(gmm_sys, dt)

  # Initial: OLS
  QtQ <- crossprod(Q_H)
  if (rcond(QtQ) < 1e-15) return(NA_real_)
  beta_hat <- solve(QtQ, crossprod(Q_H, Delta))

  for (iter in 1:10) {
    beta_old <- beta_hat

    # Estimate TWFE with current beta to get residuals
    # Reconstruct Y_hat and get residuals
    mod <- feols(Y ~ D | unit + time, data = dt)
    dt[, resid := residuals(mod)]

    # Build Omega_v blocks from residuals
    omega_blocks <- list()
    for (i in 1:N_total) {
      v_i <- dt[unit == i, resid]
      omega_blocks[[i]] <- v_i %o% v_i
    }

    # Compute structured Var[Delta]
    Var_Delta <- compute_var_delta(R_mat, omega_blocks)

    # GLS: invert Var[Delta]
    Var_Delta_inv <- tryCatch(
      solve(Var_Delta + diag(1e-6, n_did)),
      error = function(e) NULL
    )
    if (is.null(Var_Delta_inv)) break

    # GLS estimate
    QtVQ <- crossprod(Q_H, Var_Delta_inv %*% Q_H)
    QtVD <- crossprod(Q_H, Var_Delta_inv %*% Delta)

    if (rcond(QtVQ) < 1e-15) break
    beta_hat <- solve(QtVQ, QtVD)

    # Check convergence
    if (max(abs(beta_hat - beta_old)) < 1e-8) break
  }

  return(compute_att(beta_hat, n_catt))
}

# ===========================================================================
# 9. Other estimators (TWFE, CS, SA, Gardner)
# ===========================================================================

estimate_twfe <- function(dt) {
  tryCatch({
    mod <- feols(Y ~ D | unit + time, data = dt)
    return(coef(mod)["D"])
  }, error = function(e) NA_real_)
}

estimate_cs <- function(dt) {
  tryCatch({
    out <- att_gt(yname = "Y", tname = "time", idname = "unit", gname = "g_cs",
                  data = as.data.frame(dt), control_group = "nevertreated",
                  print_details = FALSE, bstrap = FALSE, cband = FALSE)
    agg <- aggte(out, type = "simple")
    return(agg$overall.att)
  }, error = function(e) NA_real_)
}

estimate_sa <- function(dt) {
  tryCatch({
    mod <- feols(Y ~ sunab(g_inf, time) | unit + time, data = dt)
    agg <- summary(mod, agg = "ATT")
    return(agg$coeftable[1, 1])
  }, error = function(e) NA_real_)
}

estimate_gardner <- function(dt) {
  tryCatch({
    dt_g <- copy(dt)
    dt_g[, first_treat := fifelse(g == 0, Inf, as.numeric(g))]
    mod <- did2s(data = as.data.frame(dt_g), yname = "Y",
                 first_stage = ~ 0 | unit + time, second_stage = ~ i(D, ref = 0),
                 treatment = "D", cluster_var = "unit", verbose = FALSE)
    return(coef(mod)["D::1"])
  }, error = function(e) NA_real_)
}

# ===========================================================================
# 10. Run Simulation
# ===========================================================================

run_simulation <- function(beta_g_vec, r_g_vec, label) {
  true_att <- compute_true_att(beta_g_vec, r_g_vec)
  cat(sprintf("\n=== %s (rho = %.2f) ===\n", label, rho))
  cat(sprintf("True ATT: %.4f\n", true_att))

  results <- data.table(
    sim = integer(), TWFE = numeric(), CS = numeric(),
    SA = numeric(), Gardner = numeric(),
    GMM_I = numeric(), GMM_iter = numeric(), GMM_GLS = numeric()
  )

  for (s in 1:n_sims) {
    if (s %% 50 == 0) cat(sprintf("  Simulation %d/%d\n", s, n_sims))

    dt <- generate_data_ar1(beta_g_vec, r_g_vec)

    att_twfe    <- estimate_twfe(dt)
    att_cs      <- estimate_cs(dt)
    att_sa      <- estimate_sa(dt)
    att_gardner <- estimate_gardner(dt)

    # Build GMM system once, use for all three methods
    gmm_sys <- build_gmm_system(dt)
    if (is.null(gmm_sys)) {
      att_gmm_i <- att_gmm_iter <- att_gmm_gls <- NA_real_
    } else {
      att_gmm_i    <- tryCatch(gmm_identity(gmm_sys), error = function(e) NA_real_)
      att_gmm_iter <- tryCatch(gmm_iterative_structured(gmm_sys, dt), error = function(e) NA_real_)
      att_gmm_gls  <- tryCatch(gmm_gls(gmm_sys, dt), error = function(e) NA_real_)
    }

    results <- rbindlist(list(results, data.table(
      sim = s, TWFE = att_twfe, CS = att_cs, SA = att_sa,
      Gardner = att_gardner,
      GMM_I = att_gmm_i, GMM_iter = att_gmm_iter, GMM_GLS = att_gmm_gls
    )))
  }

  # Summary
  est_names <- c("TWFE", "CS", "SA", "Gardner", "GMM_I", "GMM_iter", "GMM_GLS")
  summary_dt <- data.table(
    Estimator = est_names,
    Bias      = numeric(length(est_names)),
    Variance  = numeric(length(est_names))
  )
  for (i in seq_along(est_names)) {
    vals <- results[[est_names[i]]]
    vals <- vals[!is.na(vals)]
    summary_dt$Bias[i]     <- round(mean(vals) - true_att, 4)
    summary_dt$Variance[i] <- round(var(vals), 4)
  }

  cat(sprintf("\nResults (%s):\n", label))
  print(summary_dt)
  return(list(results = results, summary = summary_dt, true_att = true_att))
}

# ===========================================================================
# 11. Execute
# ===========================================================================

cat("================================================================\n")
cat("  SIMULATION WITH AR(1) ERRORS (rho =", rho, ")\n")
cat("================================================================\n")

# --- Homogeneous ---
beta_hom <- rep(-5, n_cohorts)
r_hom    <- rep(0, n_cohorts)
res_hom  <- run_simulation(beta_hom, r_hom, "Homogeneous")

# --- Heterogeneous ---
beta_het <- c(-16, -12, -10, -9, -2)
r_het    <- c(0.01, 0.04, 0.08, 0.10, 0.07)
res_het  <- run_simulation(beta_het, r_het, "Heterogeneous")

# ===========================================================================
# 12. Print Table
# ===========================================================================

cat("\n\n")
cat("====================================================================\n")
cat("  TABLE: Bias and Variance with AR(1) Errors (rho =", rho, ")\n")
cat("====================================================================\n\n")

table_out <- merge(
  res_hom$summary[, .(Estimator, Hom_Bias = Bias, Hom_Var = Variance)],
  res_het$summary[, .(Estimator, Het_Bias = Bias, Het_Var = Variance)],
  by = "Estimator"
)
table_out <- table_out[match(c("TWFE", "CS", "SA", "Gardner",
                                "GMM_I", "GMM_iter", "GMM_GLS"),
                              table_out$Estimator)]

cat(sprintf("%-12s  %12s  %12s  %12s  %12s\n",
            "Estimator", "Hom. Bias", "Hom. Var", "Het. Bias", "Het. Var"))
cat(paste(rep("-", 66), collapse = ""), "\n")
for (i in 1:nrow(table_out)) {
  cat(sprintf("%-12s  %12.4f  %12.4f  %12.4f  %12.4f\n",
              table_out$Estimator[i],
              table_out$Hom_Bias[i], table_out$Hom_Var[i],
              table_out$Het_Bias[i], table_out$Het_Var[i]))
}

save(res_hom, res_het, table_out, file = "simulation_ar1_results.RData")
cat("\nResults saved to simulation_ar1_results.RData\n")
