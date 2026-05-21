###############################################################################
# Table 2 Replication — Spherical Errors DGP
#
# DGP: Y_it = alpha_i + lambda_t + tau_it + epsilon_it
#   where epsilon_it ~ iid N(0, sigma_eps^2)  (spherical / homoskedastic errors)
#
# Purpose: Under spherical errors, Gauss-Markov implies Gardner = Flex TWFE = BLUE.
#   The two new spherical-optimal GMM estimators (GMM_Sph_Clean, GMM_Sph_All) use
#   Omega_phi constructed under the spherical assumption (only sigma_0 non-zero).
#   Key predictions to verify:
#     (i)  GMM_Sph_Clean == GMM_Sph_All  (forbidden comparisons redundant)
#     (ii) Both == Gardner == Flex_TWFE  (BLUE uniqueness under spherical errors)
#     (iii) GMM_Eff (which estimates full AR(1) Omega from residuals) has slightly
#           higher finite-sample variance due to noise from over-parameterising Omega.
#
# Estimators:
#   1. TWFE
#   2. Callaway-Sant'Anna
#   3. Sun-Abraham
#   4. Gardner (did2s)
#   5. Flex TWFE (Wooldridge 2025)
#   6. Flex TWFE FGLS (iterated feasible GLS)
#   7. GMM_Clean:      efficient GMM, clean DiDs, full AR(1) Omega estimated from residuals
#   8. GMM_Eff:        efficient GMM, all DiDs,   full AR(1) Omega estimated from residuals
#   9. GMM_Eff_HetVar: GMM, all DiDs, cohort-specific sigma, cross-cov = 0
#  10. GMM_Eff_HetCov: GMM, all DiDs, cohort-specific sigma + cross-cohort gamma
#  11. GLS_Eq15:       GLS on eq (15); algebraically == GMM_Eff
#  12. GMM_Sph_Clean:  GMM, clean DiDs only, Omega built under spherical assumption
#                      (sigma_0 estimated; sigma_d = 0 for d > 0)
#  13. GMM_Sph_All:    GMM, all DiDs,        Omega built under spherical assumption
#
# CATTs saved for: Flex_TWFE, GMM_Clean, GMM_Eff, GLS_Eq15,
#                  GMM_Sph_Clean, GMM_Sph_All
#
# Key optimization: panel structure (Q_H, Q_clean, C_mat, lag indices) is
# pre-computed once before the simulation loop.
###############################################################################

library(data.table)
library(fixest)
library(did)
library(did2s)
library(MASS)

set.seed(42)

# ===========================================================================
# 1. Parameters
# ===========================================================================

cohort_size <- 10
n_never     <- 10
N_total     <- 5 * cohort_size + n_never
T_total     <- 33
n_sims      <- 100
sigma_eps   <- 1.0   # iid error standard deviation (spherical errors)

treatment_times <- c(10, 13, 16, 19, 22)
n_cohorts       <- length(treatment_times)
unit_cohort     <- c(rep(treatment_times, each = cohort_size), rep(0, n_never))

get_cohort_size <- function(g_val) ifelse(g_val == 0, n_never, cohort_size)

# ===========================================================================
# 2. Pre-compute GMM Structure (ONCE — invariant across simulations)
# ===========================================================================

cat("Pre-computing GMM structure (Q_H, Q_clean, C_mat, lag indices)...\n")

treated_g <- sort(treatment_times)

# Enumerate CATTs
catt_list <- list()
for (g_c in treated_g) {
  for (k in 0:(T_total - g_c)) {
    catt_list[[length(catt_list) + 1]] <- c(g_c, g_c + k)
  }
}
n_catt <- length(catt_list)

# CATT to Flex TWFE coefficient mapping: treat_gt = g*100 + t
catt_gt_key <- sapply(catt_list, function(x) x[1] * 100L + as.integer(x[2]))

# Enumerate all 2x2 DiD estimate metadata
did_meta <- list()
for (catt_idx in 1:n_catt) {
  g_c    <- catt_list[[catt_idx]][1]
  t_post <- catt_list[[catt_idx]][2]
  k      <- t_post - g_c

  for (m in 1:(g_c - 1)) {
    t_pre <- g_c - m

    # Never-treated
    did_meta[[length(did_meta) + 1]] <- list(
      catt_idx = catt_idx, type = "never",
      focal_g = g_c, ctrl_g = 0, t_post = t_post, t_pre = t_pre
    )

    # Not-yet-treated
    for (g_l in treated_g) {
      if (g_l <= t_post) next
      did_meta[[length(did_meta) + 1]] <- list(
        catt_idx = catt_idx, type = "notyet",
        focal_g = g_c, ctrl_g = g_l, t_post = t_post, t_pre = t_pre
      )
    }

    # Already-treated (forbidden, bias-corrected)
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
        focal_g = g_c, ctrl_g = g_j, t_post = t_post, t_pre = t_pre
      )
    }
  }
}

n_did <- length(did_meta)

# Identify clean DiDs (never-treated or not-yet-treated control)
clean_idx <- which(sapply(did_meta, function(e) e$type != "already"))
n_clean   <- length(clean_idx)

cat(sprintf("  n_catt = %d, n_did = %d  (%d clean, %d forbidden)\n",
            n_catt, n_did, n_clean, n_did - n_clean))

# Build Q_H (incidence matrix, invariant)
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

# Q_clean: restrict Q_H to clean DiD rows
Q_clean        <- Q_H[clean_idx, , drop = FALSE]
QtQ_clean      <- crossprod(Q_clean)
QtQ_clean_inv  <- tryCatch(solve(QtQ_clean), error = function(e) ginv(QtQ_clean))

# (Keep QtQ for gmm_efficient initialization — uses all DiDs)
QtQ     <- crossprod(Q_H)
QtQ_inv <- tryCatch(solve(QtQ), error = function(e) ginv(QtQ))

# Extract metadata vectors
meta_focal <- integer(n_did)
meta_ctrl  <- integer(n_did)
meta_tp    <- integer(n_did)
meta_tr    <- integer(n_did)
for (s in 1:n_did) {
  meta_focal[s] <- did_meta[[s]]$focal_g
  meta_ctrl[s]  <- did_meta[[s]]$ctrl_g
  meta_tp[s]    <- did_meta[[s]]$t_post
  meta_tr[s]    <- did_meta[[s]]$t_pre
}

# Build C_mat (structural cohort factor, invariant)
N_f <- sapply(meta_focal, get_cohort_size)
N_c <- sapply(meta_ctrl, get_cohort_size)

cat("  Building C_mat (cohort factor matrix)...\n")
gg   <- outer(meta_focal, meta_focal, "==")
gc_m <- outer(meta_focal, meta_ctrl,  "==")
cg   <- outer(meta_ctrl,  meta_focal, "==")
cc   <- outer(meta_ctrl,  meta_ctrl,  "==")
C_mat <- sweep(gg - gc_m, 1, 1/N_f, "*") + sweep(cc - cg, 1, 1/N_c, "*")
rm(gg, gc_m, cg, cc); invisible(gc())

# C_mat_clean: submatrix restricted to clean DiDs (for GMM_Clean)
C_mat_clean <- C_mat[clean_idx, clean_idx]

# Pre-compute time lag index vectors for all DiDs (for GMM_Eff)
pp_v <- as.vector(abs(outer(meta_tp, meta_tp, "-")))
pr_v <- as.vector(abs(outer(meta_tp, meta_tr, "-")))
rp_v <- as.vector(abs(outer(meta_tr, meta_tp, "-")))
rr_v <- as.vector(abs(outer(meta_tr, meta_tr, "-")))

# Pre-compute time lag index vectors for clean DiDs (for GMM_Clean)
meta_tp_clean <- meta_tp[clean_idx]
meta_tr_clean <- meta_tr[clean_idx]
pp_v_clean <- as.vector(abs(outer(meta_tp_clean, meta_tp_clean, "-")))
pr_v_clean <- as.vector(abs(outer(meta_tp_clean, meta_tr_clean, "-")))
rp_v_clean <- as.vector(abs(outer(meta_tr_clean, meta_tp_clean, "-")))
rr_v_clean <- as.vector(abs(outer(meta_tr_clean, meta_tr_clean, "-")))

# Pre-compute cohort-level bookkeeping for GMM_Eff_HetVar / GMM_Eff_HetCov
all_cohorts_hv      <- c(0L, sort(as.integer(treatment_times)))
n_cohort_grp_hv     <- length(all_cohorts_hv)
N_by_cohort_hv      <- sapply(all_cohorts_hv, get_cohort_size)
units_by_cohort_hv  <- lapply(all_cohorts_hv, function(g) which(unit_cohort == g))
active_by_cohort_hv <- lapply(all_cohorts_hv,
  function(g) which(meta_focal == g | meta_ctrl == g))
sign_by_cohort_hv   <- lapply(seq_along(all_cohorts_hv), function(gi) {
  g   <- all_cohorts_hv[gi]
  idx <- active_by_cohort_hv[[gi]]
  ifelse(meta_focal[idx] == g, 1L, -1L)
})

cat("  Pre-computation done.\n\n")

# ===========================================================================
# 3. DGP + True ATT / True CATTs
# ===========================================================================

generate_data <- function(beta_g_vec, r_g_vec) {
  unit_id <- rep(1:N_total, each = T_total)
  time_id <- rep(1:T_total, times = N_total)
  alpha <- rnorm(N_total); lambda <- rnorm(T_total)
  eps <- rnorm(N_total * T_total, sd = sigma_eps)   # iid spherical errors
  g_vec <- unit_cohort[unit_id]
  D_vec <- as.integer(g_vec > 0 & time_id >= g_vec)
  tau_vec <- numeric(N_total * T_total)
  for (c_idx in 1:n_cohorts) {
    g_c <- treatment_times[c_idx]; mask <- (g_vec == g_c) & (D_vec == 1)
    tau_vec[mask] <- beta_g_vec[c_idx] * (1 + r_g_vec[c_idx])^(time_id[mask] - g_c)
  }
  dt <- data.table(unit = unit_id, time = time_id,
                   Y = alpha[unit_id] + lambda[time_id] + tau_vec + eps,
                   D = D_vec, g = g_vec)
  dt[, g_inf := fifelse(g == 0, Inf, as.numeric(g))]
  dt[, g_cs  := fifelse(g == 0, 0L, as.integer(g))]
  return(dt)
}

compute_true_att <- function(beta_g_vec, r_g_vec) {
  total_te <- 0; total_obs <- 0
  for (c_idx in 1:n_cohorts) {
    g_c <- treatment_times[c_idx]
    for (t in g_c:T_total) {
      total_te  <- total_te + cohort_size * beta_g_vec[c_idx] * (1 + r_g_vec[c_idx])^(t - g_c)
      total_obs <- total_obs + cohort_size
    }
  }
  return(total_te / total_obs)
}

compute_true_catt <- function(beta_g_vec, r_g_vec) {
  true_catt <- numeric(n_catt)
  for (ci in 1:n_catt) {
    g_c <- catt_list[[ci]][1]; t_c <- catt_list[[ci]][2]
    c_idx <- which(treatment_times == g_c)
    true_catt[ci] <- beta_g_vec[c_idx] * (1 + r_g_vec[c_idx])^(t_c - g_c)
  }
  return(true_catt)
}

# ===========================================================================
# 4. Compute Delta from data (per simulation)
# ===========================================================================

compute_delta <- function(dt) {
  cohort_means <- dt[, .(Y_mean = mean(Y)), by = .(g, time)]
  setkey(cohort_means, g, time)
  get_mean <- function(g_val, t_val) cohort_means[.(g_val, t_val), Y_mean]

  Delta <- numeric(n_did)
  for (s in 1:n_did) {
    e <- did_meta[[s]]
    Yf_post <- get_mean(e$focal_g, e$t_post)
    Yf_pre  <- get_mean(e$focal_g, e$t_pre)
    Yc_post <- get_mean(e$ctrl_g,  e$t_post)
    Yc_pre  <- get_mean(e$ctrl_g,  e$t_pre)
    Delta[s] <- (Yf_post - Yf_pre) - (Yc_post - Yc_pre)
  }
  return(Delta)
}

# ===========================================================================
# 5. GMM_Clean: Iterative efficient GMM on clean DiDs only
#    A = Omega_clean^{-1}, where Omega_clean is estimated from residuals
#    and restricted to the clean DiD subset.
#    Omega_clean[s,s'] = C_mat_clean[s,s'] * S_mat_clean[s,s']
#    where S_mat_clean uses only the time lags of clean DiDs.
# ===========================================================================

gmm_efficient_clean <- function(Delta, dt, max_iter = 3, tol = 1e-6) {
  Delta_clean <- Delta[clean_idx]
  beta_hat    <- as.numeric(QtQ_clean_inv %*% crossprod(Q_clean, Delta_clean))

  dt_r <- copy(dt)
  setorder(dt_r, unit, time)

  for (iter in 1:max_iter) {
    beta_old <- beta_hat

    dt_r[, tau_hat := 0]
    for (ci in 1:n_catt) {
      dt_r[g == catt_list[[ci]][1] & time == catt_list[[ci]][2],
           tau_hat := beta_hat[ci]]
    }
    dt_r[, Y_adj := Y - tau_hat]

    resid_mat <- matrix(residuals(feols(Y_adj ~ 1 | unit + time, data = dt_r)),
                        nrow = T_total, ncol = N_total)

    sigma_d <- numeric(T_total)
    for (d in 0:(T_total - 1)) {
      r1 <- 1:(T_total - d); r2 <- (1 + d):T_total
      sigma_d[d + 1] <- sum(resid_mat[r1, ] * resid_mat[r2, ]) / (N_total * (T_total - d))
    }

    # Omega for clean DiDs only
    S_vec_clean <- sigma_d[pp_v_clean + 1] - sigma_d[pr_v_clean + 1] -
                   sigma_d[rp_v_clean + 1] + sigma_d[rr_v_clean + 1]
    Omega_clean <- C_mat_clean * matrix(S_vec_clean, nrow = n_clean)
    Omega_clean <- (Omega_clean + t(Omega_clean)) / 2
    diag(Omega_clean) <- diag(Omega_clean) + 1e-6

    OQ <- tryCatch(solve(Omega_clean, Q_clean), error = function(e) NULL)
    if (is.null(OQ)) break
    OD <- solve(Omega_clean, Delta_clean)

    QtAQ <- crossprod(Q_clean, OQ)
    QtAD <- crossprod(Q_clean, OD)
    beta_hat <- as.numeric(tryCatch(solve(QtAQ, QtAD), error = function(e) beta_old))

    if (max(abs(beta_hat - beta_old)) < tol) break
  }

  list(att = mean(beta_hat), catt = beta_hat)
}

# ===========================================================================
# 6. GMM_Eff: Iterative Efficient GMM on all DiDs (Paper Eq. 29-31)
#    A = Omega_phi^{-1}, iterated from residuals
# ===========================================================================

gmm_efficient <- function(Delta, dt, max_iter = 3, tol = 1e-6) {
  beta_hat <- as.numeric(QtQ_inv %*% crossprod(Q_H, Delta))

  dt_r <- copy(dt)
  setorder(dt_r, unit, time)

  for (iter in 1:max_iter) {
    beta_old <- beta_hat

    dt_r[, tau_hat := 0]
    for (ci in 1:n_catt) {
      dt_r[g == catt_list[[ci]][1] & time == catt_list[[ci]][2],
           tau_hat := beta_hat[ci]]
    }
    dt_r[, Y_adj := Y - tau_hat]

    resid_mat <- matrix(residuals(feols(Y_adj ~ 1 | unit + time, data = dt_r)),
                        nrow = T_total, ncol = N_total)

    sigma_d <- numeric(T_total)
    for (d in 0:(T_total - 1)) {
      r1 <- 1:(T_total - d); r2 <- (1 + d):T_total
      sigma_d[d + 1] <- sum(resid_mat[r1, ] * resid_mat[r2, ]) / (N_total * (T_total - d))
    }

    S_vec <- sigma_d[pp_v + 1] - sigma_d[pr_v + 1] - sigma_d[rp_v + 1] + sigma_d[rr_v + 1]
    Omega_phi <- C_mat * matrix(S_vec, nrow = n_did)
    Omega_phi <- (Omega_phi + t(Omega_phi)) / 2
    diag(Omega_phi) <- diag(Omega_phi) + 1e-6

    OQ <- tryCatch(solve(Omega_phi, Q_H), error = function(e) NULL)
    if (is.null(OQ)) break
    OD <- solve(Omega_phi, Delta)

    QtAQ <- crossprod(Q_H, OQ)
    QtAD <- crossprod(Q_H, OD)
    beta_hat <- as.numeric(tryCatch(solve(QtAQ, QtAD), error = function(e) beta_old))

    if (max(abs(beta_hat - beta_old)) < tol) break
  }

  list(att = mean(beta_hat), catt = beta_hat)
}

# ===========================================================================
# 6b. GMM_Eff_HetVar: Efficient GMM, all DiDs, cohort-specific sigma (unit-level)
#     Cross-cohort covariances are set to zero (Assumption 1 preserved).
#     Omega_phi = sum_g  outer(sgn_g, sgn_g) * S_g / N_g
#     where S_g uses sigma_d estimated from units in cohort g only.
# ===========================================================================

gmm_efficient_hetvar <- function(Delta, dt, max_iter = 3, tol = 1e-6) {
  beta_hat <- as.numeric(QtQ_inv %*% crossprod(Q_H, Delta))

  dt_r <- copy(dt)
  setorder(dt_r, unit, time)

  for (iter in 1:max_iter) {
    beta_old <- beta_hat

    dt_r[, tau_hat := 0]
    for (ci in 1:n_catt) {
      dt_r[g == catt_list[[ci]][1] & time == catt_list[[ci]][2],
           tau_hat := beta_hat[ci]]
    }
    dt_r[, Y_adj := Y - tau_hat]

    resid_mat <- matrix(residuals(feols(Y_adj ~ 1 | unit + time, data = dt_r)),
                        nrow = T_total, ncol = N_total)

    Omega_phi <- matrix(0.0, nrow = n_did, ncol = n_did)

    for (gi in seq_along(all_cohorts_hv)) {
      N_g     <- N_by_cohort_hv[gi]
      units_g <- units_by_cohort_hv[[gi]]
      active  <- active_by_cohort_hv[[gi]]
      if (length(active) == 0L) next
      sgn <- sign_by_cohort_hv[[gi]]

      rm_g <- resid_mat[, units_g, drop = FALSE]
      sigma_d_g <- numeric(T_total)
      for (d in 0L:(T_total - 1L)) {
        r1 <- 1L:(T_total - d); r2 <- (1L + d):T_total
        sigma_d_g[d + 1L] <- sum(rm_g[r1, ] * rm_g[r2, ]) / (N_g * (T_total - d))
      }

      tp_a <- meta_tp[active]; tr_a <- meta_tr[active]
      pp <- abs(outer(tp_a, tp_a, "-")); pr <- abs(outer(tp_a, tr_a, "-"))
      rp <- abs(outer(tr_a, tp_a, "-")); rr <- abs(outer(tr_a, tr_a, "-"))
      S_g <- sigma_d_g[pp + 1L] - sigma_d_g[pr + 1L] -
             sigma_d_g[rp + 1L] + sigma_d_g[rr + 1L]

      Omega_phi[active, active] <- Omega_phi[active, active] + outer(sgn, sgn) * S_g / N_g
    }

    Omega_phi <- (Omega_phi + t(Omega_phi)) / 2
    diag(Omega_phi) <- diag(Omega_phi) + 1e-6

    OQ <- tryCatch(solve(Omega_phi, Q_H), error = function(e) NULL)
    if (is.null(OQ)) break
    OD <- solve(Omega_phi, Delta)

    QtAQ <- crossprod(Q_H, OQ)
    QtAD <- crossprod(Q_H, OD)
    beta_hat <- as.numeric(tryCatch(solve(QtAQ, QtAD), error = function(e) beta_old))

    if (max(abs(beta_hat - beta_old)) < tol) break
  }

  list(att = mean(beta_hat), catt = beta_hat)
}

# ===========================================================================
# 6c. GMM_Eff_HetCov: Efficient GMM, all DiDs, cohort-specific sigma (unit-level)
#     plus non-zero symmetric cross-cohort gamma (no regularization).
#     gamma_d^{g,g'} estimated as cohort-mean-level lag-d cross-covariance,
#     symmetrised: (1/(2(T-d))) * [sum_t bar_e_{g,t}*bar_e_{g',t+d} + reverse].
# ===========================================================================

gmm_efficient_hetcov <- function(Delta, dt, max_iter = 3, tol = 1e-6) {
  beta_hat <- as.numeric(QtQ_inv %*% crossprod(Q_H, Delta))

  dt_r <- copy(dt)
  setorder(dt_r, unit, time)

  for (iter in 1:max_iter) {
    beta_old <- beta_hat

    dt_r[, tau_hat := 0]
    for (ci in 1:n_catt) {
      dt_r[g == catt_list[[ci]][1] & time == catt_list[[ci]][2],
           tau_hat := beta_hat[ci]]
    }
    dt_r[, Y_adj := Y - tau_hat]

    resid_mat <- matrix(residuals(feols(Y_adj ~ 1 | unit + time, data = dt_r)),
                        nrow = T_total, ncol = N_total)

    # Cohort-mean residuals: T_total x n_cohort_grp_hv
    bar_e <- matrix(0.0, nrow = T_total, ncol = n_cohort_grp_hv)
    for (gi in seq_along(all_cohorts_hv)) {
      units_g <- units_by_cohort_hv[[gi]]
      bar_e[, gi] <- rowMeans(resid_mat[, units_g, drop = FALSE])
    }

    Omega_phi <- matrix(0.0, nrow = n_did, ncol = n_did)

    # Within-cohort contributions (identical to HetVar)
    for (gi in seq_along(all_cohorts_hv)) {
      N_g     <- N_by_cohort_hv[gi]
      units_g <- units_by_cohort_hv[[gi]]
      active  <- active_by_cohort_hv[[gi]]
      if (length(active) == 0L) next
      sgn <- sign_by_cohort_hv[[gi]]

      rm_g <- resid_mat[, units_g, drop = FALSE]
      sigma_d_g <- numeric(T_total)
      for (d in 0L:(T_total - 1L)) {
        r1 <- 1L:(T_total - d); r2 <- (1L + d):T_total
        sigma_d_g[d + 1L] <- sum(rm_g[r1, ] * rm_g[r2, ]) / (N_g * (T_total - d))
      }

      tp_a <- meta_tp[active]; tr_a <- meta_tr[active]
      pp <- abs(outer(tp_a, tp_a, "-")); pr <- abs(outer(tp_a, tr_a, "-"))
      rp <- abs(outer(tr_a, tp_a, "-")); rr <- abs(outer(tr_a, tr_a, "-"))
      S_g <- sigma_d_g[pp + 1L] - sigma_d_g[pr + 1L] -
             sigma_d_g[rp + 1L] + sigma_d_g[rr + 1L]

      Omega_phi[active, active] <- Omega_phi[active, active] + outer(sgn, sgn) * S_g / N_g
    }

    # Cross-cohort contributions: symmetric gamma_{g,g'} for each unordered pair
    if (n_cohort_grp_hv > 1L) {
      for (gi in 1L:(n_cohort_grp_hv - 1L)) {
        active_g <- active_by_cohort_hv[[gi]]
        sgn_g    <- sign_by_cohort_hv[[gi]]
        if (length(active_g) == 0L) next

        for (gpi in (gi + 1L):n_cohort_grp_hv) {
          active_gp <- active_by_cohort_hv[[gpi]]
          sgn_gp    <- sign_by_cohort_hv[[gpi]]
          if (length(active_gp) == 0L) next

          gamma_d <- numeric(T_total)
          for (d in 0L:(T_total - 1L)) {
            r1 <- 1L:(T_total - d); r2 <- (1L + d):T_total
            gamma_d[d + 1L] <- (sum(bar_e[r1, gi] * bar_e[r2, gpi]) +
                                 sum(bar_e[r1, gpi] * bar_e[r2, gi])) /
                                (2.0 * (T_total - d))
          }

          tp_g  <- meta_tp[active_g];  tr_g  <- meta_tr[active_g]
          tp_gp <- meta_tp[active_gp]; tr_gp <- meta_tr[active_gp]
          pp_c <- abs(outer(tp_g, tp_gp, "-")); pr_c <- abs(outer(tp_g, tr_gp, "-"))
          rp_c <- abs(outer(tr_g, tp_gp, "-")); rr_c <- abs(outer(tr_g, tr_gp, "-"))
          S_cross <- gamma_d[pp_c + 1L] - gamma_d[pr_c + 1L] -
                     gamma_d[rp_c + 1L] + gamma_d[rr_c + 1L]

          cross_block <- outer(sgn_g, sgn_gp) * S_cross
          Omega_phi[active_g,  active_gp] <- Omega_phi[active_g,  active_gp] + cross_block
          Omega_phi[active_gp, active_g]  <- Omega_phi[active_gp, active_g]  + t(cross_block)
        }
      }
    }

    Omega_phi <- (Omega_phi + t(Omega_phi)) / 2
    diag(Omega_phi) <- diag(Omega_phi) + 1e-6

    OQ <- tryCatch(solve(Omega_phi, Q_H), error = function(e) NULL)
    if (is.null(OQ)) break
    OD <- solve(Omega_phi, Delta)

    QtAQ <- crossprod(Q_H, OQ)
    QtAD <- crossprod(Q_H, OD)
    beta_hat <- as.numeric(tryCatch(solve(QtAQ, QtAD), error = function(e) beta_old))

    if (max(abs(beta_hat - beta_old)) < tol) break
  }

  list(att = mean(beta_hat), catt = beta_hat)
}

# ===========================================================================
# 6d. GLS_Eq15: Direct GLS on the moment system (paper eq. 15)
#     The system is Delta = Q_H * beta + xi,  Var[xi] = Omega_phi.
#     GLS minimises (Delta - Q_H*beta)' Omega_phi^{-1} (Delta - Q_H*beta),
#     giving beta_hat = (Q_H' Omega_phi^{-1} Q_H)^{-1} Q_H' Omega_phi^{-1} Delta.
#     This is algebraically identical to GMM_Eff with A = Omega_phi^{-1}
#     (paper Section 3.3.2, p.11). Included to verify numerical equivalence.
# ===========================================================================

gmm_gls_eq15 <- function(Delta, dt, max_iter = 3, tol = 1e-6) {
  beta_hat <- as.numeric(QtQ_inv %*% crossprod(Q_H, Delta))

  dt_r <- copy(dt)
  setorder(dt_r, unit, time)

  for (iter in 1:max_iter) {
    beta_old <- beta_hat

    dt_r[, tau_hat := 0]
    for (ci in 1:n_catt) {
      dt_r[g == catt_list[[ci]][1] & time == catt_list[[ci]][2],
           tau_hat := beta_hat[ci]]
    }
    dt_r[, Y_adj := Y - tau_hat]

    resid_mat <- matrix(residuals(feols(Y_adj ~ 1 | unit + time, data = dt_r)),
                        nrow = T_total, ncol = N_total)

    sigma_d <- numeric(T_total)
    for (d in 0:(T_total - 1)) {
      r1 <- 1:(T_total - d); r2 <- (1 + d):T_total
      sigma_d[d + 1] <- sum(resid_mat[r1, ] * resid_mat[r2, ]) / (N_total * (T_total - d))
    }

    S_vec     <- sigma_d[pp_v + 1] - sigma_d[pr_v + 1] - sigma_d[rp_v + 1] + sigma_d[rr_v + 1]
    Omega_phi <- C_mat * matrix(S_vec, nrow = n_did)
    Omega_phi <- (Omega_phi + t(Omega_phi)) / 2
    diag(Omega_phi) <- diag(Omega_phi) + 1e-6

    # GLS closed form: (Q_H' Omega^{-1} Q_H)^{-1} Q_H' Omega^{-1} Delta
    OQ <- tryCatch(solve(Omega_phi, Q_H), error = function(e) NULL)
    if (is.null(OQ)) break
    OD <- solve(Omega_phi, Delta)

    QtOQ <- crossprod(Q_H, OQ)   # Q_H' Omega^{-1} Q_H
    QtOD <- crossprod(Q_H, OD)   # Q_H' Omega^{-1} Delta
    beta_hat <- as.numeric(tryCatch(solve(QtOQ, QtOD), error = function(e) beta_old))

    if (max(abs(beta_hat - beta_old)) < tol) break
  }

  list(att = mean(beta_hat), catt = beta_hat)
}

# ===========================================================================
# 6e. GMM_Sph_Clean: Efficient GMM on clean DiDs, Omega built under spherical
#     assumption — only sigma_0 = Var[eps_it] is non-zero; sigma_d = 0 for d > 0.
#     Under spherical errors this is the true optimal A; should coincide with
#     Gardner / Flex TWFE asymptotically (unique BLUE).
# ===========================================================================

gmm_sph_clean <- function(Delta, dt, max_iter = 3, tol = 1e-6) {
  Delta_clean <- Delta[clean_idx]
  beta_hat    <- as.numeric(QtQ_clean_inv %*% crossprod(Q_clean, Delta_clean))

  dt_r <- copy(dt)
  setorder(dt_r, unit, time)

  for (iter in 1:max_iter) {
    beta_old <- beta_hat

    dt_r[, tau_hat := 0]
    for (ci in 1:n_catt) {
      dt_r[g == catt_list[[ci]][1] & time == catt_list[[ci]][2],
           tau_hat := beta_hat[ci]]
    }
    dt_r[, Y_adj := Y - tau_hat]

    resid_mat <- matrix(residuals(feols(Y_adj ~ 1 | unit + time, data = dt_r)),
                        nrow = T_total, ncol = N_total)

    # Spherical: estimate sigma_0 only; set all higher-lag covariances to zero
    sigma2_hat <- sum(resid_mat^2) / (N_total * T_total)
    sigma_d    <- c(sigma2_hat, rep(0.0, T_total - 1L))

    S_vec_clean <- sigma_d[pp_v_clean + 1L] - sigma_d[pr_v_clean + 1L] -
                   sigma_d[rp_v_clean + 1L] + sigma_d[rr_v_clean + 1L]
    Omega_clean <- C_mat_clean * matrix(S_vec_clean, nrow = n_clean)
    Omega_clean <- (Omega_clean + t(Omega_clean)) / 2
    diag(Omega_clean) <- diag(Omega_clean) + 1e-6

    OQ <- tryCatch(solve(Omega_clean, Q_clean), error = function(e) NULL)
    if (is.null(OQ)) break
    OD <- solve(Omega_clean, Delta_clean)

    QtAQ <- crossprod(Q_clean, OQ)
    QtAD <- crossprod(Q_clean, OD)
    beta_hat <- as.numeric(tryCatch(solve(QtAQ, QtAD), error = function(e) beta_old))

    if (max(abs(beta_hat - beta_old)) < tol) break
  }

  list(att = mean(beta_hat), catt = beta_hat)
}

# ===========================================================================
# 6f. GMM_Sph_All: Efficient GMM on all DiDs, Omega built under spherical
#     assumption — same sigma_d structure as GMM_Sph_Clean but uses Q_H (all DiDs).
#     Under Assumption 1, forbidden comparisons are redundant, so this should be
#     numerically identical to GMM_Sph_Clean.
# ===========================================================================

gmm_sph_all <- function(Delta, dt, max_iter = 3, tol = 1e-6) {
  beta_hat <- as.numeric(QtQ_inv %*% crossprod(Q_H, Delta))

  dt_r <- copy(dt)
  setorder(dt_r, unit, time)

  for (iter in 1:max_iter) {
    beta_old <- beta_hat

    dt_r[, tau_hat := 0]
    for (ci in 1:n_catt) {
      dt_r[g == catt_list[[ci]][1] & time == catt_list[[ci]][2],
           tau_hat := beta_hat[ci]]
    }
    dt_r[, Y_adj := Y - tau_hat]

    resid_mat <- matrix(residuals(feols(Y_adj ~ 1 | unit + time, data = dt_r)),
                        nrow = T_total, ncol = N_total)

    # Spherical: sigma_0 only
    sigma2_hat <- sum(resid_mat^2) / (N_total * T_total)
    sigma_d    <- c(sigma2_hat, rep(0.0, T_total - 1L))

    S_vec     <- sigma_d[pp_v + 1L] - sigma_d[pr_v + 1L] -
                 sigma_d[rp_v + 1L] + sigma_d[rr_v + 1L]
    Omega_phi <- C_mat * matrix(S_vec, nrow = n_did)
    Omega_phi <- (Omega_phi + t(Omega_phi)) / 2
    diag(Omega_phi) <- diag(Omega_phi) + 1e-6

    OQ <- tryCatch(solve(Omega_phi, Q_H), error = function(e) NULL)
    if (is.null(OQ)) break
    OD <- solve(Omega_phi, Delta)

    QtAQ <- crossprod(Q_H, OQ)
    QtAD <- crossprod(Q_H, OD)
    beta_hat <- as.numeric(tryCatch(solve(QtAQ, QtAD), error = function(e) beta_old))

    if (max(abs(beta_hat - beta_old)) < tol) break
  }

  list(att = mean(beta_hat), catt = beta_hat)
}

# ===========================================================================
# 7. Package-based Estimators (scalar ATT only)
# ===========================================================================

estimate_twfe <- function(dt) {
  tryCatch({ coef(feols(Y ~ D | unit + time, data = dt))["D"] },
           error = function(e) NA_real_)
}

estimate_cs <- function(dt) {
  tryCatch({
    out <- att_gt(yname = "Y", tname = "time", idname = "unit", gname = "g_cs",
                  data = as.data.frame(dt), control_group = "nevertreated",
                  print_details = FALSE, bstrap = FALSE, cband = FALSE)
    aggte(out, type = "simple")$overall.att
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
    dt_g <- copy(dt); dt_g[, first_treat := fifelse(g == 0, Inf, as.numeric(g))]
    coef(did2s(data = as.data.frame(dt_g), yname = "Y",
               first_stage = ~ 0 | unit + time, second_stage = ~ i(D, ref = 0),
               treatment = "D", cluster_var = "unit", verbose = FALSE))["D::1"]
  }, error = function(e) NA_real_)
}

# ===========================================================================
# 8. Flexible TWFE (Wooldridge 2025) — returns list: att + catt vector
# ===========================================================================

estimate_flex_twfe <- function(dt) {
  tryCatch({
    dt_f <- copy(dt)
    dt_f[, treat_gt := fifelse(D == 1, g * 100L + as.integer(time), 0L)]
    mod <- feols(Y ~ i(treat_gt, ref = 0) | unit + time, data = dt_f)
    coef_vals <- coef(mod)
    coef_keys <- as.integer(gsub("treat_gt::", "", names(coef_vals)))
    idx <- match(catt_gt_key, coef_keys)
    catt_vec <- ifelse(is.na(idx), NA_real_, coef_vals[idx])
    list(att = mean(catt_vec, na.rm = TRUE), catt = catt_vec)
  }, error = function(e) list(att = NA_real_, catt = rep(NA_real_, n_catt)))
}

# ===========================================================================
# Flex TWFE with Iterated Feasible GLS (exploiting estimated error structure)
# ===========================================================================

estimate_flex_twfe_fgls <- function(dt, max_iter = 5, tol = 1e-4) {
  tryCatch({
    dt_f <- copy(dt)
    setorder(dt_f, unit, time)

    var_per_unit_old <- NULL

    for (iter in 1:max_iter) {
      dt_f[, treat_gt := fifelse(D == 1, g * 100L + as.integer(time), 0L)]

      if (iter == 1) {
        dt_f[, obs_wt := 1.0]
      }

      mod <- tryCatch(
        feols(Y ~ i(treat_gt, ref = 0) | unit + time, data = dt_f, weights = ~obs_wt),
        error = function(e) NULL
      )
      if (is.null(mod)) break

      resid_vec <- residuals(mod)
      resid_mat <- matrix(resid_vec, nrow = T_total, ncol = N_total)

      sigma_d <- numeric(T_total)
      for (d in 0:(T_total - 1)) {
        r1 <- 1:(T_total - d)
        r2 <- (1 + d):T_total
        sigma_d[d + 1] <- sum(resid_mat[r1, ] * resid_mat[r2, ]) / (N_total * (T_total - d))
      }

      var_per_unit <- numeric(N_total)
      for (i in 1:N_total) {
        var_per_unit[i] <- mean(resid_mat[, i]^2)
      }

      if (!is.null(var_per_unit_old)) {
        if (max(abs(var_per_unit - var_per_unit_old)) / mean(var_per_unit) < tol) break
      }

      dt_f[, obs_wt := 1 / (var_per_unit[unit] + 1e-6)]
      dt_f[, obs_wt := obs_wt / mean(obs_wt)]

      var_per_unit_old <- var_per_unit
    }

    coef_vals <- coef(mod)
    coef_keys <- as.integer(gsub("treat_gt::", "", names(coef_vals)))
    idx <- match(catt_gt_key, coef_keys)
    catt_vec <- ifelse(is.na(idx), NA_real_, coef_vals[idx])
    list(att = mean(catt_vec, na.rm = TRUE), catt = catt_vec)
  }, error = function(e) list(att = NA_real_, catt = rep(NA_real_, n_catt)))
}

# ===========================================================================
# 9. Run Simulation — saves CATTs for Flex_TWFE, GMM_Clean, GMM_Eff
# ===========================================================================

run_simulation <- function(beta_g_vec, r_g_vec, label) {
  true_att  <- compute_true_att(beta_g_vec, r_g_vec)
  true_catt <- compute_true_catt(beta_g_vec, r_g_vec)
  cat(sprintf("\n=== %s (spherical errors, sigma=%.2f, N=%d) ===\nTrue ATT: %.4f\n",
              label, sigma_eps, N_total, true_att))

  est_names <- c("TWFE", "CS", "SA", "Gardner",
                 "Flex_TWFE", "Flex_TWFE_FGLS",
                 "GMM_Clean", "GMM_Eff", "GMM_Eff_HetVar", "GMM_Eff_HetCov",
                 "GLS_Eq15", "GMM_Sph_Clean", "GMM_Sph_All")
  n_est <- length(est_names)
  results <- data.table(sim = integer())
  for (nm in est_names) results[, (nm) := numeric()]

  # CATT storage matrices (n_sims × n_catt)
  catt_est  <- c("Flex_TWFE", "GMM_Clean", "GMM_Eff", "GLS_Eq15",
                 "GMM_Sph_Clean", "GMM_Sph_All")
  catt_mats <- setNames(
    lapply(catt_est, function(x) matrix(NA_real_, nrow = n_sims, ncol = n_catt)),
    catt_est
  )

  na_catt <- rep(NA_real_, n_catt)

  for (s in 1:n_sims) {
    t0_sim <- proc.time()[3]
    dt    <- generate_data(beta_g_vec, r_g_vec)
    Delta <- compute_delta(dt)

    att_twfe    <- estimate_twfe(dt)
    att_cs      <- estimate_cs(dt)
    att_sa      <- estimate_sa(dt)
    att_gardner <- estimate_gardner(dt)

    gmm_cl_res   <- tryCatch(gmm_efficient_clean(Delta, dt),
                              error = function(e) list(att = NA_real_, catt = na_catt))
    gmm_eff_res  <- tryCatch(gmm_efficient(Delta, dt),
                              error = function(e) list(att = NA_real_, catt = na_catt))
    gmm_hv_res   <- tryCatch(gmm_efficient_hetvar(Delta, dt),
                              error = function(e) list(att = NA_real_, catt = na_catt))
    gmm_hc_res   <- tryCatch(gmm_efficient_hetcov(Delta, dt),
                              error = function(e) list(att = NA_real_, catt = na_catt))
    gls_eq15_res  <- tryCatch(gmm_gls_eq15(Delta, dt),
                               error = function(e) list(att = NA_real_, catt = na_catt))
    gmm_sph_cl_res <- tryCatch(gmm_sph_clean(Delta, dt),
                               error = function(e) list(att = NA_real_, catt = na_catt))
    gmm_sph_al_res <- tryCatch(gmm_sph_all(Delta, dt),
                               error = function(e) list(att = NA_real_, catt = na_catt))
    flex_res      <- estimate_flex_twfe(dt)
    flex_fgls_res <- estimate_flex_twfe_fgls(dt)

    # Store CATTs
    catt_mats[["Flex_TWFE"]][s, ]      <- flex_res$catt
    catt_mats[["GMM_Clean"]][s, ]      <- gmm_cl_res$catt
    catt_mats[["GMM_Eff"]][s, ]        <- gmm_eff_res$catt
    catt_mats[["GLS_Eq15"]][s, ]       <- gls_eq15_res$catt
    catt_mats[["GMM_Sph_Clean"]][s, ]  <- gmm_sph_cl_res$catt
    catt_mats[["GMM_Sph_All"]][s, ]    <- gmm_sph_al_res$catt

    elapsed <- round(proc.time()[3] - t0_sim, 1)
    cat(sprintf("  Sim %d/%d  (%.1fs)\n", s, n_sims, elapsed))
    flush(stdout())

    results <- rbindlist(list(results, data.table(
      sim = s, TWFE = att_twfe, CS = att_cs, SA = att_sa,
      Gardner = att_gardner,
      Flex_TWFE = flex_res$att, Flex_TWFE_FGLS = flex_fgls_res$att,
      GMM_Clean = gmm_cl_res$att, GMM_Eff = gmm_eff_res$att,
      GMM_Eff_HetVar = gmm_hv_res$att, GMM_Eff_HetCov = gmm_hc_res$att,
      GLS_Eq15 = gls_eq15_res$att,
      GMM_Sph_Clean = gmm_sph_cl_res$att, GMM_Sph_All = gmm_sph_al_res$att
    )))
  }

  # Summarize ATT bias and variance
  summary_dt <- data.table(Estimator = est_names, Bias = numeric(n_est), Variance = numeric(n_est))
  for (i in 1:n_est) {
    vals <- results[[est_names[i]]]; vals <- vals[!is.na(vals)]
    summary_dt$Bias[i]     <- round(mean(vals) - true_att, 4)
    summary_dt$Variance[i] <- round(var(vals), 4)
  }

  # Summarize CATT bias and variance (per-CATT, averaged)
  catt_summary <- list()
  for (nm in catt_est) {
    mat <- catt_mats[[nm]]
    catt_bias <- colMeans(mat, na.rm = TRUE) - true_catt
    catt_var  <- apply(mat, 2, var, na.rm = TRUE)
    catt_summary[[nm]] <- data.table(
      catt_idx = 1:n_catt,
      g = sapply(catt_list, `[`, 1),
      t = sapply(catt_list, `[`, 2),
      true_catt = true_catt,
      bias = catt_bias,
      variance = catt_var
    )
  }

  cat(sprintf("\nResults (%s):\n", label)); print(summary_dt)
  cat(sprintf("\nCATT summary (mean |bias| / mean variance across %d CATTs):\n", n_catt))
  for (nm in catt_est) {
    cs <- catt_summary[[nm]]
    cat(sprintf("  %-12s  mean|bias|=%.4f  meanVar=%.4f\n",
                nm, mean(abs(cs$bias)), mean(cs$variance)))
  }

  return(list(results = results, summary = summary_dt, true_att = true_att,
              true_catt = true_catt, catt_mats = catt_mats, catt_summary = catt_summary))
}

# ===========================================================================
# 10. Execute
# ===========================================================================

cat("================================================================\n")
cat(sprintf("  N=%d (%d/cohort), T=%d, sigma=%.1f, %d sims, 13 estimators\n",
            N_total, cohort_size, T_total, sigma_eps, n_sims))
cat("  DGP: SPHERICAL errors (iid N(0, sigma^2))\n")
cat("  GMM_Sph_Clean: optimal GMM for spherical errors, clean DiDs only\n")
cat("  GMM_Sph_All:   optimal GMM for spherical errors, all DiDs\n")
cat("  Predictions: GMM_Sph_Clean == GMM_Sph_All (redundancy)\n")
cat("               Both == Gardner == Flex_TWFE (unique BLUE)\n")
cat("================================================================\n")

beta_hom <- rep(-5, n_cohorts); r_hom <- rep(0, n_cohorts)
res_hom  <- run_simulation(beta_hom, r_hom, "Homogeneous")

cat("\n\n")
beta_het <- c(-16, -12, -10, -9, -2); r_het <- c(0.01, 0.04, 0.08, 0.10, 0.07)
res_het  <- run_simulation(beta_het, r_het, "Heterogeneous")

# ===========================================================================
# 11. Display Table
# ===========================================================================

cat("\n\n")
cat("==========================================================================\n")
cat(sprintf("  TABLE: Bias and Variance — Spherical Errors (sigma=%.1f), N=%d (%d/cohort)\n",
            sigma_eps, N_total, cohort_size))
cat("  GMM_Sph_Clean:  A=Omega_sph^{-1}, clean DiDs, sigma_0 only\n")
cat("  GMM_Sph_All:    A=Omega_sph^{-1}, all DiDs,   sigma_0 only\n")
cat("  GMM_Clean/Eff:  A=Omega_ar1^{-1} estimated from residuals (over-parameterised)\n")
cat("  Expected: GMM_Sph_Clean==GMM_Sph_All==Gardner==Flex_TWFE (BLUE uniqueness)\n")
cat("==========================================================================\n\n")

table2 <- merge(
  res_hom$summary[, .(Estimator, Hom_Bias = Bias, Hom_Var = Variance)],
  res_het$summary[, .(Estimator, Het_Bias = Bias, Het_Var = Variance)],
  by = "Estimator"
)
est_order <- c("TWFE", "CS", "SA", "Gardner",
               "Flex_TWFE", "Flex_TWFE_FGLS",
               "GMM_Clean", "GMM_Eff", "GMM_Eff_HetVar", "GMM_Eff_HetCov",
               "GLS_Eq15", "GMM_Sph_Clean", "GMM_Sph_All")
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
     catt_list, catt_gt_key,
     file = "simulation_spherical_errors_results.RData")
cat("\nResults (incl. CATTs) saved to simulation_spherical_errors_results.RData\n")

################## Plots ###############
library(ggplot2)
data_plot <- data.frame(
  rbind(cbind(apply(res_het$catt_mats$Flex_TWFE,     2, var), 1),
        cbind(apply(res_het$catt_mats$GMM_Eff,        2, var), 2),
        cbind(apply(res_het$catt_mats$GMM_Sph_Clean,  2, var), 3),
        cbind(apply(res_het$catt_mats$GMM_Sph_All,    2, var), 4)),
  c(1:n_catt, 1:n_catt, 1:n_catt, 1:n_catt)
)
colnames(data_plot) <- c("Variance", "Estimator", "Index")
data_plot[data_plot$Estimator == 1, 2] <- "Flex_TWFE"
data_plot[data_plot$Estimator == 2, 2] <- "GMM_Eff"
data_plot[data_plot$Estimator == 3, 2] <- "GMM_Sph_Clean"
data_plot[data_plot$Estimator == 4, 2] <- "GMM_Sph_All"

ggplot(data = data_plot, aes(x = Index, y = Variance, group = Estimator)) +
  geom_line(aes(color = Estimator)) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  labs(title = "Per-CATT variance under spherical errors",
       subtitle = "GMM_Sph_Clean, GMM_Sph_All, Flex_TWFE should overlap")

## =========================================================
## Verification checks (key predictions under spherical errors)
## =========================================================

cat("\n=== Verification: Redundancy — GMM_Sph_Clean vs GMM_Sph_All ===\n")
max_diff_sph_hom <- max(abs(res_hom$results$GMM_Sph_Clean - res_hom$results$GMM_Sph_All),
                        na.rm = TRUE)
max_diff_sph_het <- max(abs(res_het$results$GMM_Sph_Clean - res_het$results$GMM_Sph_All),
                        na.rm = TRUE)
cat(sprintf("  Max |GMM_Sph_Clean - GMM_Sph_All| (Hom): %.2e\n", max_diff_sph_hom))
cat(sprintf("  Max |GMM_Sph_Clean - GMM_Sph_All| (Het): %.2e\n", max_diff_sph_het))
cat("  (Should be ~machine precision: forbidden comparisons algebraically redundant)\n")

cat("\n=== Verification: BLUE uniqueness — GMM_Sph_Clean vs Gardner ===\n")
var_sph  <- mean(apply(res_het$catt_mats$GMM_Sph_Clean, 2, var))
var_gard <- var(res_het$results$Gardner, na.rm = TRUE)
var_flex <- mean(apply(res_het$catt_mats$Flex_TWFE, 2, var))
var_gmm  <- mean(apply(res_het$catt_mats$GMM_Eff, 2, var))
cat(sprintf("  Mean CATT var — GMM_Sph_Clean: %.6f\n", var_sph))
cat(sprintf("  Mean CATT var — Flex_TWFE:     %.6f\n", var_flex))
cat(sprintf("  ATT   var    — Gardner:         %.6f\n", var_gard))
cat(sprintf("  Mean CATT var — GMM_Eff (AR1):  %.6f\n", var_gmm))
cat("  (GMM_Sph_Clean should match Flex_TWFE; GMM_Eff may be slightly higher\n")
cat("   in finite samples due to noise from estimating the full AR(1) Omega)\n")

cat("\n=== Verification: GLS_Eq15 == GMM_Eff ===\n")
max_gls_hom <- max(abs(res_hom$results$GLS_Eq15 - res_hom$results$GMM_Eff), na.rm = TRUE)
max_gls_het <- max(abs(res_het$results$GLS_Eq15 - res_het$results$GMM_Eff), na.rm = TRUE)
cat(sprintf("  Max |GLS_Eq15 - GMM_Eff| (Hom): %.2e\n", max_gls_hom))
cat(sprintf("  Max |GLS_Eq15 - GMM_Eff| (Het): %.2e\n", max_gls_het))
cat("  (Should be ~machine precision)\n")
