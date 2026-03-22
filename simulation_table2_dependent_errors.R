###############################################################################
# Table 2 Replication with AR(1) Dependent Errors
#
# DGP: Y_it = alpha_i + lambda_t + tau_it + epsilon_it
#   where epsilon_i ~ MVN(0, Sigma) with Sigma[t,s] = rho^|t-s|
#   Errors are independent across units i but dependent within unit across t.
#
# Estimators:
#   1. TWFE
#   2. Callaway-Sant'Anna (CS)
#   3. Sun-Abraham (SA)
#   4. Gardner (did2s)
#   5. GMM with A = I
#   6. Flexible TWFE (Wooldridge 2025) - cohort-time specific treatment effects
#   7. Flexible TWFE + Iterated GLS - with block-diagonal covariance
#
# 100 simulations
###############################################################################

library(data.table)
library(fixest)
library(did)
library(did2s)
library(MASS)

set.seed(42)

# ===========================================================================
# 1. Simulation Parameters
# ===========================================================================

N_total     <- 51
T_total     <- 33
n_sims      <- 100
cohort_size <- 8
n_never     <- 11
rho         <- 0.5  # AR(1) autocorrelation parameter

treatment_times <- c(10, 13, 16, 19, 22)
n_cohorts       <- length(treatment_times)

unit_cohort <- c(rep(treatment_times, each = cohort_size), rep(0, n_never))

# AR(1) covariance matrix: Sigma[t,s] = rho^|t-s|
Sigma_ar1 <- toeplitz(rho^(0:(T_total - 1)))

# ===========================================================================
# 2. Data Generating Process
# ===========================================================================

generate_data <- function(beta_g_vec, r_g_vec) {
  unit_id <- rep(1:N_total, each = T_total)
  time_id <- rep(1:T_total, times = N_total)

  alpha <- rnorm(N_total)
  lambda <- rnorm(T_total)

  # AR(1) errors: independent across i, dependent within i across t
  eps_mat <- mvrnorm(N_total, mu = rep(0, T_total), Sigma = Sigma_ar1)
  eps <- as.vector(t(eps_mat))  # row-major flatten

  g_vec <- unit_cohort[unit_id]
  D_vec <- as.integer(g_vec > 0 & time_id >= g_vec)

  tau_vec <- numeric(N_total * T_total)
  for (c_idx in 1:n_cohorts) {
    g_c  <- treatment_times[c_idx]
    mask <- (g_vec == g_c) & (D_vec == 1)
    tau_vec[mask] <- beta_g_vec[c_idx] * (1 + r_g_vec[c_idx])^(time_id[mask] - g_c)
  }

  Y_vec <- alpha[unit_id] + lambda[time_id] + tau_vec + eps

  dt <- data.table(
    unit = unit_id, time = time_id, Y = Y_vec, D = D_vec, g = g_vec
  )
  dt[, g_inf := fifelse(g == 0, Inf, as.numeric(g))]
  dt[, g_cs  := fifelse(g == 0, 0L, as.integer(g))]
  return(dt)
}

# ===========================================================================
# 3. True ATT (unit-weighted)
# ===========================================================================

compute_true_att <- function(beta_g_vec, r_g_vec) {
  total_te <- 0; total_obs <- 0
  for (c_idx in 1:n_cohorts) {
    g_c <- treatment_times[c_idx]
    for (t in g_c:T_total) {
      te <- beta_g_vec[c_idx] * (1 + r_g_vec[c_idx])^(t - g_c)
      total_te  <- total_te + cohort_size * te
      total_obs <- total_obs + cohort_size
    }
  }
  return(total_te / total_obs)
}

# ===========================================================================
# 4. GMM Estimator with A = I
# ===========================================================================

gmm_estimator <- function(dt) {
  cohort_means <- dt[, .(Y_mean = mean(Y)), by = .(g, time)]
  setkey(cohort_means, g, time)

  get_mean <- function(g_val, t_val) {
    res <- cohort_means[.(g_val, t_val), Y_mean]
    if (length(res) == 0 || is.na(res)) return(NA_real_)
    return(res)
  }

  treated_g <- sort(treatment_times)

  # Enumerate CATTs
  catt_list <- list()
  for (g_c in treated_g) {
    for (k in 0:(T_total - g_c)) {
      catt_list[[length(catt_list) + 1]] <- c(g_c, g_c + k)
    }
  }
  n_catt <- length(catt_list)

  # Enumerate 2x2 DiD estimates
  did_estimates <- list()

  for (catt_idx in 1:n_catt) {
    g_c    <- catt_list[[catt_idx]][1]
    t_post <- catt_list[[catt_idx]][2]
    k      <- t_post - g_c

    for (m in 1:(g_c - 1)) {
      t_pre <- g_c - m

      Y_g_post <- get_mean(g_c, t_post)
      Y_g_pre  <- get_mean(g_c, t_pre)

      # Never-treated control
      Y_0_post <- get_mean(0, t_post)
      Y_0_pre  <- get_mean(0, t_pre)
      if (!any(is.na(c(Y_g_post, Y_g_pre, Y_0_post, Y_0_pre)))) {
        did_estimates[[length(did_estimates) + 1]] <- list(
          delta = (Y_g_post - Y_g_pre) - (Y_0_post - Y_0_pre),
          catt_idx = catt_idx, type = "never"
        )
      }

      # Not-yet-treated control
      for (g_l in treated_g) {
        if (g_l <= t_post) next
        Y_l_post <- get_mean(g_l, t_post)
        Y_l_pre  <- get_mean(g_l, t_pre)
        if (!any(is.na(c(Y_g_post, Y_g_pre, Y_l_post, Y_l_pre)))) {
          did_estimates[[length(did_estimates) + 1]] <- list(
            delta = (Y_g_post - Y_g_pre) - (Y_l_post - Y_l_pre),
            catt_idx = catt_idx, type = "notyet"
          )
        }
      }

      # Already-treated control (forbidden comparison, bias-corrected via Q_H)
      for (g_j in treated_g) {
        j <- g_c - g_j
        if (j <= m) next
        if (g_j >= g_c) next
        Y_j_post <- get_mean(g_j, t_post)
        Y_j_pre  <- get_mean(g_j, t_pre)
        if (!any(is.na(c(Y_g_post, Y_g_pre, Y_j_post, Y_j_pre)))) {
          bias_neg_idx <- NULL; bias_pos_idx <- NULL
          for (ci in 1:n_catt) {
            if (catt_list[[ci]][1] == g_j && catt_list[[ci]][2] == t_post) bias_neg_idx <- ci
            if (catt_list[[ci]][1] == g_j && catt_list[[ci]][2] == t_pre)  bias_pos_idx <- ci
          }
          did_estimates[[length(did_estimates) + 1]] <- list(
            delta = (Y_g_post - Y_g_pre) - (Y_j_post - Y_j_pre),
            catt_idx = catt_idx, type = "already",
            bias_neg = bias_neg_idx, bias_pos = bias_pos_idx
          )
        }
      }
    }
  }

  n_did <- length(did_estimates)
  if (n_did == 0) return(list(att = NA))

  Delta <- numeric(n_did)
  Q_H   <- matrix(0, nrow = n_did, ncol = n_catt)

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

  # GMM with A = I
  beta_hat <- tryCatch(
    solve(crossprod(Q_H), crossprod(Q_H, Delta)),
    error = function(e) ginv(crossprod(Q_H)) %*% crossprod(Q_H, Delta)
  )

  weights <- rep(cohort_size, n_catt)
  weights <- weights / sum(weights)
  att <- sum(weights * beta_hat)

  return(list(att = as.numeric(att)))
}

# ===========================================================================
# 5. Standard Estimators
# ===========================================================================

estimate_twfe <- function(dt) {
  tryCatch({
    mod <- feols(Y ~ D | unit + time, data = dt)
    return(coef(mod)["D"])
  }, error = function(e) NA_real_)
}

estimate_cs <- function(dt) {
  tryCatch({
    out <- att_gt(
      yname = "Y", tname = "time", idname = "unit", gname = "g_cs",
      data = as.data.frame(dt), control_group = "nevertreated",
      print_details = FALSE, bstrap = FALSE, cband = FALSE
    )
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
    mod <- did2s(
      data = as.data.frame(dt_g), yname = "Y",
      first_stage = ~ 0 | unit + time, second_stage = ~ i(D, ref = 0),
      treatment = "D", cluster_var = "unit", verbose = FALSE
    )
    return(coef(mod)["D::1"])
  }, error = function(e) NA_real_)
}

# ===========================================================================
# 6. Flexible TWFE (Wooldridge 2025)
#
#    Y_it = alpha_i + lambda_t + sum_{g,t: t>=g} delta_{g,t} * 1{g_i=g, time=t}
#    Each (cohort, time) cell gets its own treatment coefficient.
#    ATT = unit-weighted average of all delta_{g,t}.
# ===========================================================================

estimate_flex_twfe <- function(dt) {
  tryCatch({
    dt_f <- copy(dt)
    # Unique numeric ID for each treated (g, time) cell; 0 = untreated reference
    dt_f[, treat_gt := fifelse(D == 1, g * 100L + as.integer(time), 0L)]

    mod <- feols(Y ~ i(treat_gt, ref = 0) | unit + time, data = dt_f)
    coefs <- coef(mod)

    # Unit-weighted ATT: simple average (all cohorts same size)
    att <- mean(coefs)
    return(att)
  }, error = function(e) NA_real_)
}

# ===========================================================================
# 7. Flexible TWFE + Iterated GLS
#
#    Same model as Flex TWFE, but estimated via iterated Feasible GLS.
#    The variance-covariance matrix Omega is block-diagonal:
#      - Across units i: zero covariance
#      - Within unit i across time t: non-zero (estimated T x T matrix Sigma)
#    Sigma is estimated from OLS/GLS residuals and iterated until convergence.
# ===========================================================================

estimate_flex_twfe_gls <- function(dt, max_iter = 50, tol = 1e-6) {
  tryCatch({
    dt_f <- copy(dt)
    setorder(dt_f, unit, time)

    Y     <- dt_f$Y
    N_obs <- nrow(dt_f)
    n_u   <- N_total
    n_t   <- T_total

    # --- Build explicit design matrix ---
    # Treatment cells: unique (g, time) where D=1
    treat_cells <- unique(dt_f[D == 1, .(g, time)])
    setorder(treat_cells, g, time)
    n_tr <- nrow(treat_cells)

    # Columns: intercept (1) + unit FE (n_u-1) + time FE (n_t-1) + treatment (n_tr)
    n_col  <- 1 + (n_u - 1) + (n_t - 1) + n_tr
    col_u  <- 2                         # unit FE start at col 2
    col_t  <- col_u + (n_u - 1)         # time FE start
    col_tr <- col_t + (n_t - 1)         # treatment dummies start

    X <- matrix(0, nrow = N_obs, ncol = n_col)
    X[, 1] <- 1  # intercept

    # Unit FE (drop last unit as reference)
    for (j in 1:(n_u - 1)) {
      X[dt_f$unit == j, col_u + j - 1] <- 1
    }

    # Time FE (drop last time as reference)
    for (j in 1:(n_t - 1)) {
      X[dt_f$time == j, col_t + j - 1] <- 1
    }

    # Treatment cell dummies
    for (j in 1:n_tr) {
      X[dt_f$g == treat_cells$g[j] & dt_f$time == treat_cells$time[j],
        col_tr + j - 1] <- 1
    }

    # --- Initial OLS ---
    beta_hat <- solve(crossprod(X), crossprod(X, Y))

    # --- Iterated GLS ---
    for (iter in 1:max_iter) {
      beta_old <- beta_hat

      # Residuals in original space
      resid_vec <- as.numeric(Y - X %*% beta_hat)

      # Reshape: T x N matrix (each column = one unit's residual vector)
      resid_mat <- matrix(resid_vec, nrow = n_t, ncol = n_u)

      # Estimate within-unit covariance: Sigma_hat = (1/N) sum_i e_i e_i'
      Sigma_hat <- tcrossprod(resid_mat) / n_u

      # Ensure symmetry and positive definiteness
      Sigma_hat <- (Sigma_hat + t(Sigma_hat)) / 2
      eig <- eigen(Sigma_hat, symmetric = TRUE)
      eig$values <- pmax(eig$values, 1e-8)
      Sigma_hat <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)

      # Cholesky: R'R = Sigma_hat, then P = (R')^{-1} so P Sigma P' = I
      R     <- chol(Sigma_hat)
      R_inv <- backsolve(R, diag(n_t))
      P     <- t(R_inv)  # lower triangular transformation matrix

      # Transform data unit by unit
      Y_star <- numeric(N_obs)
      X_star <- matrix(0, nrow = N_obs, ncol = n_col)
      for (i in 1:n_u) {
        idx <- ((i - 1) * n_t + 1):(i * n_t)
        Y_star[idx]    <- P %*% Y[idx]
        X_star[idx, ]  <- P %*% X[idx, ]
      }

      # GLS estimate = OLS on transformed data
      beta_hat <- solve(crossprod(X_star), crossprod(X_star, Y_star))

      # Convergence check
      if (max(abs(beta_hat - beta_old)) < tol) break
    }

    # Extract treatment coefficients (last n_tr columns)
    treat_coefs <- beta_hat[col_tr:(col_tr + n_tr - 1)]

    # Unit-weighted ATT: simple average (all cohorts have equal size)
    att <- mean(treat_coefs)
    return(as.numeric(att))

  }, error = function(e) NA_real_)
}

# ===========================================================================
# 8. Run Simulations
# ===========================================================================

run_simulation <- function(beta_g_vec, r_g_vec, label) {
  true_att <- compute_true_att(beta_g_vec, r_g_vec)
  cat(sprintf("\n=== %s (rho = %.2f) ===\n", label, rho))
  cat(sprintf("True ATT: %.4f\n", true_att))

  est_names <- c("TWFE", "CS", "SA", "Gardner", "GMM",
                 "Flex_TWFE", "Flex_TWFE_GLS")
  n_est <- length(est_names)

  results <- data.table(sim = integer())
  for (nm in est_names) results[, (nm) := numeric()]

  for (s in 1:n_sims) {
    if (s %% 10 == 0) cat(sprintf("  Simulation %d/%d\n", s, n_sims))

    dt <- generate_data(beta_g_vec, r_g_vec)

    att_twfe          <- estimate_twfe(dt)
    att_cs            <- estimate_cs(dt)
    att_sa            <- estimate_sa(dt)
    att_gardner       <- estimate_gardner(dt)
    att_gmm           <- gmm_estimator(dt)$att
    att_flex_twfe     <- estimate_flex_twfe(dt)
    att_flex_twfe_gls <- estimate_flex_twfe_gls(dt)

    results <- rbindlist(list(results, data.table(
      sim = s, TWFE = att_twfe, CS = att_cs, SA = att_sa,
      Gardner = att_gardner, GMM = att_gmm,
      Flex_TWFE = att_flex_twfe, Flex_TWFE_GLS = att_flex_twfe_gls
    )))
  }

  # Compute bias and variance
  summary_dt <- data.table(
    Estimator = est_names,
    Bias      = numeric(n_est),
    Variance  = numeric(n_est)
  )
  for (i in 1:n_est) {
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
# 9. Execute
# ===========================================================================

cat("================================================================\n")
cat("  SIMULATION: AR(1) DEPENDENT ERRORS (rho =", rho, ")\n")
cat("  100 simulations, 7 estimators\n")
cat("================================================================\n")

# --- Homogeneous Effects ---
# Y_it = alpha_i + lambda_t - 5*D_it + eps_it
beta_hom <- rep(-5, n_cohorts)
r_hom    <- rep(0, n_cohorts)
res_hom  <- run_simulation(beta_hom, r_hom, "Homogeneous")

# --- Heterogeneous Effects ---
# Y_it = alpha_i + lambda_t + [beta_g*(1+r_g)^(t-g_i)] * D_it + eps_it
cat("\n\n")
beta_het <- c(-16, -12, -10, -9, -2)
r_het    <- c(0.01, 0.04, 0.08, 0.10, 0.07)
res_het  <- run_simulation(beta_het, r_het, "Heterogeneous")

# ===========================================================================
# 10. Combine and Display Table 2
# ===========================================================================

cat("\n\n")
cat("========================================================================\n")
cat("  TABLE 2: Bias and Variance of Estimators — AR(1) Errors (rho =", rho, ")\n")
cat("  Errors: independent across i, AR(1) within i (Sigma[t,s] = rho^|t-s|)\n")
cat("  GMM uses A = I. Flex TWFE GLS uses iterated block-diagonal FGLS.\n")
cat("========================================================================\n\n")

table2 <- merge(
  res_hom$summary[, .(Estimator, Hom_Bias = Bias, Hom_Var = Variance)],
  res_het$summary[, .(Estimator, Het_Bias = Bias, Het_Var = Variance)],
  by = "Estimator"
)

est_order <- c("TWFE", "CS", "SA", "Gardner", "GMM", "Flex_TWFE", "Flex_TWFE_GLS")
table2 <- table2[match(est_order, table2$Estimator)]

cat(sprintf("%-16s  %12s  %12s  %12s  %12s\n",
            "Estimator", "Hom. Bias", "Hom. Var", "Het. Bias", "Het. Var"))
cat(paste(rep("-", 70), collapse = ""), "\n")
for (i in 1:nrow(table2)) {
  cat(sprintf("%-16s  %12.4f  %12.4f  %12.4f  %12.4f\n",
              table2$Estimator[i],
              table2$Hom_Bias[i], table2$Hom_Var[i],
              table2$Het_Bias[i], table2$Het_Var[i]))
}

save(res_hom, res_het, table2,
     file = "simulation_dependent_errors_results.RData")
cat("\nResults saved to simulation_dependent_errors_results.RData\n")
