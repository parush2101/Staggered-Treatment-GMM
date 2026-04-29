###############################################################################
# Compare J-test calibration: A = I vs A = (Omega + eps*I)^{-1}
# AR1 errors, paired comparison across 1000 independent simulations
#
# For each simulation we:
#   1. Generate one fresh AR(1) panel
#   2. Estimate sigma^2_eps and rho once (from A=I residuals, two-step style)
#   3. Draw ONE stratified subset of K=15 moments (shared across both variants)
#   4. Run the Satterthwaite-Box J-test twice on that subset:
#        VARIANT I — beta = (Q'Q)^{-1} Q' Delta              (current code)
#        VARIANT O — beta = (Q'AQ)^{-1} Q'A Delta,
#                    A = (Omega + eps*I)^{-1}, Omega = full m x m covariance
#   5. Record naive and Satterthwaite-corrected p-values for both
#
# Pairing on the same data + subset removes Monte Carlo noise from the
# comparison: every difference in rejection rate is attributable to A.
###############################################################################

library(data.table)
library(MASS)

set.seed(123)

# ===========================================================================
# 1. Parameters
# ===========================================================================

cohort_size <- 50
n_never     <- 60
N_total     <- 5 * cohort_size + n_never
T_total     <- 33
K_draw      <- 15L
n_per_type  <- 5L
n_sims      <- 1000L
rho_true    <- 0.5
eps_reg     <- 1e-6   # ridge factor for Omega regularization (relative to max eigenvalue)

treatment_times <- c(10, 13, 16, 19, 22)
n_cohorts       <- length(treatment_times)
unit_cohort     <- c(rep(treatment_times, each = cohort_size), rep(0, n_never))
get_cohort_size <- function(g_val) ifelse(g_val == 0, n_never, cohort_size)

all_groups  <- c(0, treatment_times)
n_groups    <- length(all_groups)
group_sizes <- sapply(all_groups, get_cohort_size)
group_masks <- lapply(all_groups, function(g) which(unit_cohort == g))

# ===========================================================================
# 2. CATT and DiD-moment structure
# ===========================================================================

treated_g <- sort(treatment_times)
catt_list <- list()
for (g_c in treated_g) {
  for (k in 0:(T_total - g_c)) {
    catt_list[[length(catt_list) + 1]] <- c(g_c, g_c + k)
  }
}
n_catt <- length(catt_list)
catt_g <- sapply(catt_list, `[`, 1)
catt_t <- sapply(catt_list, `[`, 2)

did_meta <- list()
for (catt_idx in 1:n_catt) {
  g_c    <- catt_list[[catt_idx]][1]
  t_post <- catt_list[[catt_idx]][2]
  for (m in 1:(g_c - 1)) {
    t_pre <- g_c - m
    did_meta[[length(did_meta) + 1]] <- list(
      catt_idx = catt_idx, type = "never",
      focal_g = g_c, ctrl_g = 0, t_post = t_post, t_pre = t_pre, lag_m = m)
    for (g_l in treated_g) {
      if (g_l <= t_post) next
      did_meta[[length(did_meta) + 1]] <- list(
        catt_idx = catt_idx, type = "notyet",
        focal_g = g_c, ctrl_g = g_l, t_post = t_post, t_pre = t_pre, lag_m = m)
    }
    for (g_j in treated_g) {
      j <- g_c - g_j
      if (j <= m || g_j >= g_c) next
      bias_neg <- NULL; bias_pos <- NULL
      for (ci in 1:n_catt) {
        if (catt_list[[ci]][1] == g_j && catt_list[[ci]][2] == t_post) bias_neg <- ci
        if (catt_list[[ci]][1] == g_j && catt_list[[ci]][2] == t_pre)  bias_pos <- ci
      }
      did_meta[[length(did_meta) + 1]] <- list(
        catt_idx = catt_idx, type = "already",
        bias_neg = bias_neg, bias_pos = bias_pos,
        focal_g = g_c, ctrl_g = g_j, t_post = t_post, t_pre = t_pre, lag_m = m)
    }
  }
}
n_did <- length(did_meta)

Q_H <- matrix(0, nrow = n_did, ncol = n_catt)
for (s in 1:n_did) {
  est <- did_meta[[s]]
  Q_H[s, est$catt_idx] <- 1
  if (est$type == "already") {
    if (!is.null(est$bias_neg) && !is.na(est$bias_neg))
      Q_H[s, est$bias_neg] <- Q_H[s, est$bias_neg] - 1
    if (!is.null(est$bias_pos) && !is.na(est$bias_pos))
      Q_H[s, est$bias_pos] <- Q_H[s, est$bias_pos] + 1
  }
}

# Identity-weighted shortcuts (used by VARIANT I)
QtQ     <- crossprod(Q_H)
QtQ_inv <- solve(QtQ)

meta_type   <- sapply(did_meta, `[[`, "type")
idx_never   <- which(meta_type == "never")
idx_notyet  <- which(meta_type == "notyet")
idx_already <- which(meta_type == "already")

cat(sprintf("n_catt = %d, n_did = %d\n", n_catt, n_did))
cat(sprintf("  never: %d, notyet: %d, already: %d\n",
            length(idx_never), length(idx_notyet), length(idx_already)))

# ===========================================================================
# 3. Vectorized Delta indices and Z_g matrices
# ===========================================================================

meta_focal <- sapply(did_meta, `[[`, "focal_g")
meta_ctrl  <- sapply(did_meta, `[[`, "ctrl_g")
meta_tp    <- sapply(did_meta, `[[`, "t_post")
meta_tr    <- sapply(did_meta, `[[`, "t_pre")
focal_gi <- match(meta_focal, all_groups)
ctrl_gi  <- match(meta_ctrl, all_groups)
idx_fp <- cbind(focal_gi, meta_tp)
idx_fr <- cbind(focal_gi, meta_tr)
idx_cp <- cbind(ctrl_gi,  meta_tp)
idx_cr <- cbind(ctrl_gi,  meta_tr)

Z_group <- list()
for (gi in 1:n_groups) {
  g_val <- all_groups[gi]
  Z_g <- matrix(0, nrow = n_did, ncol = T_total)
  for (s in 1:n_did) {
    e <- did_meta[[s]]
    if (e$focal_g == g_val) {
      Z_g[s, e$t_post] <- Z_g[s, e$t_post] + 1
      Z_g[s, e$t_pre]  <- Z_g[s, e$t_pre]  - 1
    }
    if (e$ctrl_g == g_val) {
      Z_g[s, e$t_post] <- Z_g[s, e$t_post] - 1
      Z_g[s, e$t_pre]  <- Z_g[s, e$t_pre]  + 1
    }
  }
  Z_group[[gi]] <- Z_g
}

# Pre-compute P_g for VARIANT I: (Q'Q)^{-1} Q' Z_g
P_group_I <- lapply(Z_group, function(Z) QtQ_inv %*% crossprod(Q_H, Z))

cat("Precomputation done.\n")

# ===========================================================================
# 4. Treatment effects
# ===========================================================================

beta_het <- c(-16, -12, -10, -9, -2)
r_het    <- c(0.01, 0.04, 0.08, 0.10, 0.07)
tau_lookup <- matrix(0, n_cohorts, T_total)
for (c_idx in 1:n_cohorts) {
  g_c <- treatment_times[c_idx]
  for (t in g_c:T_total) {
    tau_lookup[c_idx, t] <- beta_het[c_idx] * (1 + r_het[c_idx])^(t - g_c)
  }
}
unit_cohort_idx <- match(unit_cohort, treatment_times)

# ===========================================================================
# 5. Helpers for VARIANT O (optimal-weighted beta)
# ===========================================================================

build_Omega_full <- function(sigma2_eps, R_ar1) {
  Omega <- matrix(0, n_did, n_did)
  for (gi in 1:n_groups) {
    Z_g <- Z_group[[gi]]
    Omega <- Omega + (sigma2_eps / group_sizes[gi]) * (Z_g %*% R_ar1) %*% t(Z_g)
  }
  Omega
}

# Regularized inverse via eigendecomposition: A = V diag(1 / (lambda + eps)) V'
reg_inverse <- function(M, eps) {
  e <- eigen(M, symmetric = TRUE)
  e$vectors %*% (1 / (e$values + eps) * t(e$vectors))
}

# ===========================================================================
# 6. Containers
# ===========================================================================

# VARIANT I: A = I
T_raw_I    <- numeric(n_sims)
p_naive_I  <- numeric(n_sims)
p_satter_I <- numeric(n_sims)
nu_eff_I   <- numeric(n_sims)
c_scale_I  <- numeric(n_sims)
valid_I    <- logical(n_sims)

# VARIANT O: A = (Omega + eps I)^{-1}
T_raw_O    <- numeric(n_sims)
p_naive_O  <- numeric(n_sims)
p_satter_O <- numeric(n_sims)
nu_eff_O   <- numeric(n_sims)
c_scale_O  <- numeric(n_sims)
valid_O    <- logical(n_sims)

rho_vec     <- numeric(n_sims)
sigma2e_vec <- numeric(n_sims)

# ===========================================================================
# 7. Main loop
# ===========================================================================

cat(sprintf("\nRunning %d simulations (rho=%.1f, ridge factor=%.0e)...\n",
            n_sims, rho_true, eps_reg))
t_start <- Sys.time()

for (sim in 1:n_sims) {

  # --- Generate AR(1) data -------------------------------------------------
  alpha_i  <- rnorm(N_total)
  lambda_t <- rnorm(T_total)
  eps_mat <- matrix(0, N_total, T_total)
  eps_mat[, 1] <- rnorm(N_total, 0, sd = 1 / sqrt(1 - rho_true^2))
  for (t in 2:T_total) {
    eps_mat[, t] <- rho_true * eps_mat[, t - 1] + rnorm(N_total)
  }
  Y_mat <- outer(alpha_i, rep(1, T_total)) + outer(rep(1, N_total), lambda_t) + eps_mat
  for (c_idx in 1:n_cohorts) {
    units_c <- which(unit_cohort_idx == c_idx)
    for (t in treatment_times[c_idx]:T_total) {
      Y_mat[units_c, t] <- Y_mat[units_c, t] + tau_lookup[c_idx, t]
    }
  }

  # --- Group-time means and Delta ------------------------------------------
  means_mat <- matrix(0, n_groups, T_total)
  for (gi in 1:n_groups) {
    means_mat[gi, ] <- colMeans(Y_mat[group_masks[[gi]], , drop = FALSE])
  }
  Delta <- means_mat[idx_fp] - means_mat[idx_fr] - means_mat[idx_cp] + means_mat[idx_cr]

  # --- First-step beta with A = I (used by VARIANT I and to seed rho/sigma)-
  beta_I <- as.numeric(QtQ_inv %*% crossprod(Q_H, Delta))

  # --- Two-way demeaned residuals -> sigma^2_eps and rho -------------------
  Y_adj_mat <- Y_mat
  for (c_idx in 1:n_cohorts) {
    units_c <- which(unit_cohort_idx == c_idx)
    for (ci in 1:n_catt) {
      if (catt_g[ci] == treatment_times[c_idx]) {
        Y_adj_mat[units_c, catt_t[ci]] <- Y_adj_mat[units_c, catt_t[ci]] - beta_I[ci]
      }
    }
  }
  row_m <- rowMeans(Y_adj_mat)
  col_m <- colMeans(Y_adj_mat)
  grd_m <- mean(Y_adj_mat)
  resid_mat <- Y_adj_mat - outer(row_m, rep(1, T_total)) -
               outer(rep(1, N_total), col_m) + grd_m

  resid_lag  <- resid_mat[, 1:(T_total - 1)]
  resid_cur  <- resid_mat[, 2:T_total]
  rho_hat    <- sum(resid_cur * resid_lag) / sum(resid_lag^2)
  eta_hat    <- resid_cur - rho_hat * resid_lag
  sigma2_eta <- mean(eta_hat^2)
  sigma2_eps <- sigma2_eta / max(1 - rho_hat^2, 0.01)

  rho_vec[sim]     <- rho_hat
  sigma2e_vec[sim] <- sigma2_eps

  R_ar1 <- toeplitz(rho_hat^(0:(T_total - 1)))

  # --- Common K=15 stratified subset (shared across variants) --------------
  idx_K <- c(sample(idx_never,   n_per_type),
             sample(idx_notyet,  n_per_type),
             sample(idx_already, n_per_type))
  K   <- length(idx_K)
  Q_S <- Q_H[idx_K, , drop = FALSE]

  # --- Subset Omega_S (depends only on sigma^2_eps, rho, idx_K) ------------
  Omega_S <- matrix(0, K, K)
  for (gi in 1:n_groups) {
    w_g  <- sigma2_eps / group_sizes[gi]
    Z_gS <- Z_group[[gi]][idx_K, , drop = FALSE]
    Omega_S <- Omega_S + w_g * (Z_gS %*% R_ar1) %*% t(Z_gS)
  }
  eig_O <- eigen(Omega_S, symmetric = TRUE)

  if (min(eig_O$values) < 1e-12) {
    valid_I[sim] <- FALSE; valid_O[sim] <- FALSE
    T_raw_I[sim] <- NA; p_naive_I[sim] <- NA; p_satter_I[sim] <- NA
    T_raw_O[sim] <- NA; p_naive_O[sim] <- NA; p_satter_O[sim] <- NA
    next
  }
  Omega_S_inv <- solve(Omega_S)
  sqrt_inv    <- eig_O$vectors %*%
                 diag(1 / sqrt(eig_O$values)) %*% t(eig_O$vectors)

  # =========================================================================
  # VARIANT I: A = I
  # =========================================================================
  e_K_I <- (Delta - Q_H %*% beta_I)[idx_K]
  T_raw_I[sim]   <- as.numeric(t(e_K_I) %*% Omega_S_inv %*% e_K_I)
  p_naive_I[sim] <- 1 - pchisq(T_raw_I[sim], df = K)

  V_S_I <- matrix(0, K, K)
  for (gi in 1:n_groups) {
    w_g  <- sigma2_eps / group_sizes[gi]
    Z_gS <- Z_group[[gi]][idx_K, , drop = FALSE]
    R_gS <- Z_gS - Q_S %*% P_group_I[[gi]]
    V_S_I <- V_S_I + w_g * (R_gS %*% R_ar1) %*% t(R_gS)
  }
  M_I  <- sqrt_inv %*% V_S_I %*% sqrt_inv
  lam  <- pmax(eigen(M_I, symmetric = TRUE)$values, 0)
  s1   <- sum(lam); s2 <- sum(lam^2)
  if (s1 < 1e-10) {
    p_satter_I[sim] <- 1; nu_eff_I[sim] <- 0; c_scale_I[sim] <- NA
  } else {
    c_scale_I[sim]  <- s2 / s1
    nu_eff_I[sim]   <- s1^2 / s2
    p_satter_I[sim] <- 1 - pchisq(T_raw_I[sim] / c_scale_I[sim], df = nu_eff_I[sim])
  }
  valid_I[sim] <- TRUE

  # =========================================================================
  # VARIANT O: A = (Omega_full + eps_use * I)^{-1}
  # =========================================================================
  Omega_full <- build_Omega_full(sigma2_eps, R_ar1)
  max_eig    <- max(abs(eigen(Omega_full, symmetric = TRUE,
                              only.values = TRUE)$values))
  eps_use    <- eps_reg * max_eig
  A_opt      <- reg_inverse(Omega_full, eps_use)

  QtA      <- crossprod(Q_H, A_opt)        # P x n_did
  QtAQ     <- QtA %*% Q_H                  # P x P
  QtAQ_inv <- tryCatch(solve(QtAQ), error = function(e) NULL)

  if (is.null(QtAQ_inv)) {
    valid_O[sim] <- FALSE
    T_raw_O[sim] <- NA; p_naive_O[sim] <- NA; p_satter_O[sim] <- NA
  } else {
    beta_O <- as.numeric(QtAQ_inv %*% (QtA %*% Delta))
    e_K_O  <- (Delta - Q_H %*% beta_O)[idx_K]
    T_raw_O[sim]   <- as.numeric(t(e_K_O) %*% Omega_S_inv %*% e_K_O)
    p_naive_O[sim] <- 1 - pchisq(T_raw_O[sim], df = K)

    V_S_O <- matrix(0, K, K)
    for (gi in 1:n_groups) {
      w_g  <- sigma2_eps / group_sizes[gi]
      Z_gS <- Z_group[[gi]][idx_K, , drop = FALSE]
      P_gO <- QtAQ_inv %*% (QtA %*% Z_group[[gi]])   # P x T  (depends on A)
      R_gS <- Z_gS - Q_S %*% P_gO
      V_S_O <- V_S_O + w_g * (R_gS %*% R_ar1) %*% t(R_gS)
    }
    M_O  <- sqrt_inv %*% V_S_O %*% sqrt_inv
    lam  <- pmax(eigen(M_O, symmetric = TRUE)$values, 0)
    s1   <- sum(lam); s2 <- sum(lam^2)
    if (s1 < 1e-10) {
      p_satter_O[sim] <- 1; nu_eff_O[sim] <- 0; c_scale_O[sim] <- NA
    } else {
      c_scale_O[sim]  <- s2 / s1
      nu_eff_O[sim]   <- s1^2 / s2
      p_satter_O[sim] <- 1 - pchisq(T_raw_O[sim] / c_scale_O[sim], df = nu_eff_O[sim])
    }
    valid_O[sim] <- TRUE
  }

  if (sim %% 100 == 0) cat(sprintf("  sim %d/%d\n", sim, n_sims))
}

elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
cat(sprintf("Done. %.1f seconds.\n", elapsed))

# ===========================================================================
# 8. Comparison table
# ===========================================================================

alphas <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90)

reject_naive_I  <- sapply(alphas, function(a) mean(p_naive_I[valid_I]  < a, na.rm = TRUE))
reject_satter_I <- sapply(alphas, function(a) mean(p_satter_I[valid_I] < a, na.rm = TRUE))
reject_naive_O  <- sapply(alphas, function(a) mean(p_naive_O[valid_O]  < a, na.rm = TRUE))
reject_satter_O <- sapply(alphas, function(a) mean(p_satter_O[valid_O] < a, na.rm = TRUE))

cat("\n")
cat("==========================================================================\n")
cat("  J-TEST CALIBRATION COMPARISON: A = I  vs  A = (Omega + eps*I)^{-1}\n")
cat(sprintf("  N=%d, T=%d, rho=%.1f, K=%d stratified moments, %d sims\n",
            N_total, T_total, rho_true, K_draw, n_sims))
cat(sprintf("  Ridge factor: eps = %.0e * max eigenvalue of Omega\n", eps_reg))
cat("==========================================================================\n\n")

cat(sprintf("Valid sims:  A=I: %d/%d   A=Omega^-1: %d/%d\n",
            sum(valid_I), n_sims, sum(valid_O), n_sims))
cat(sprintf("rho   estimates:  mean=%.4f  sd=%.4f  (true=%.1f)\n",
            mean(rho_vec), sd(rho_vec), rho_true))
cat(sprintf("sigma^2_eps est:  mean=%.4f  sd=%.4f  (true=%.4f)\n\n",
            mean(sigma2e_vec), sd(sigma2e_vec), 1 / (1 - rho_true^2)))

cat("--- SATTERTHWAITE PARAMETERS ---\n\n")
cat(sprintf("  %-25s  %12s  %12s\n", "", "A = I", "A = Omega^-1"))
cat(paste(rep("-", 55), collapse = ""), "\n")
cat(sprintf("  %-25s  %12.2f  %12.2f\n", "mean nu (eff df)",
            mean(nu_eff_I[valid_I],  na.rm = TRUE),
            mean(nu_eff_O[valid_O],  na.rm = TRUE)))
cat(sprintf("  %-25s  %12.4f  %12.4f\n", "mean c (scale)",
            mean(c_scale_I[valid_I], na.rm = TRUE),
            mean(c_scale_O[valid_O], na.rm = TRUE)))
cat(sprintf("  %-25s  %12.2f  %12.2f\n", "mean T raw",
            mean(T_raw_I[valid_I],   na.rm = TRUE),
            mean(T_raw_O[valid_O],   na.rm = TRUE)))
cat(sprintf("  %-25s  %12.2f  %12.2f\n", "sd   T raw",
            sd(T_raw_I[valid_I],     na.rm = TRUE),
            sd(T_raw_O[valid_O],     na.rm = TRUE)))

cat("\n--- REJECTION RATES (Satterthwaite-corrected) ---\n\n")
cat(sprintf("  %-12s  %10s  %10s  %14s  %12s\n",
            "Nominal a", "Expected", "A = I", "A = Omega^-1", "Diff (O-I)"))
cat(paste(rep("-", 65), collapse = ""), "\n")
for (i in seq_along(alphas)) {
  cat(sprintf("  %-12.2f  %9.1f%%  %9.1f%%  %13.1f%%  %+10.1f pp\n",
              alphas[i], 100 * alphas[i],
              100 * reject_satter_I[i],
              100 * reject_satter_O[i],
              100 * (reject_satter_O[i] - reject_satter_I[i])))
}

cat("\n--- REJECTION RATES (Naive chi^2(K), for reference) ---\n\n")
cat(sprintf("  %-12s  %10s  %10s  %14s\n",
            "Nominal a", "Expected", "A = I", "A = Omega^-1"))
cat(paste(rep("-", 50), collapse = ""), "\n")
for (i in seq_along(alphas)) {
  cat(sprintf("  %-12.2f  %9.1f%%  %9.1f%%  %13.1f%%\n",
              alphas[i], 100 * alphas[i],
              100 * reject_naive_I[i],
              100 * reject_naive_O[i]))
}

cat("\n--- INTERPRETATION ---\n\n")
cat("  - 'A = I' is the current procedure (identity weighting for beta).\n")
cat("  - 'A = Omega^-1' uses the regularized full inverse covariance for beta.\n")
cat("  - Rejection rates closer to the Expected column = better size calibration.\n")
cat("  - Diff (O-I) > 0  means optimal-weighting variant rejects MORE often.\n")
cat("  - Diff (O-I) < 0  means optimal-weighting variant rejects LESS often.\n")
cat("  - Comparison is paired: same data and same moment subset per sim,\n")
cat("    so differences are attributable to A, not Monte Carlo noise.\n")

# Save
dir.create(path.expand("~/econ_sims"), showWarnings = FALSE, recursive = TRUE)
save(T_raw_I, p_naive_I, p_satter_I, nu_eff_I, c_scale_I, valid_I,
     T_raw_O, p_naive_O, p_satter_O, nu_eff_O, c_scale_O, valid_O,
     rho_vec, sigma2e_vec,
     file = path.expand("~/econ_sims/j_test_compare_A.RData"))
cat("\nResults saved to ~/econ_sims/j_test_compare_A.RData\n")
