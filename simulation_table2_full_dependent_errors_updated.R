###############################################################################
# Table 2 Replication with Full (NT×NT) Dependent Errors — Large N
# UPDATED VERSION: Removes GMM_I, GMM_Diag, GMM_Eff_HomoX;
#                  Varies rho_fac in {0.3, 0.5, 0.7};
#                  Varies sd(Phi) in {0.5, 0.9, 2.0}
#
# DGP: Y_it = alpha_i + lambda_t + tau_it + epsilon_it
#   where epsilon_{i,t} = Phi[i,] %*% F_time[t,] + sqrt(D_unit[i]) * N(0,1)
#   Phi    : N x k_fac unit-specific factor loadings (fixed per phi_sd setting)
#   F_time : T x k_fac common AR(1) factors (rho_fac), drawn each simulation
#   D_unit : N unit-specific idiosyncratic variances (fixed)
#
# Estimators retained:
#   1. TWFE
#   2. Callaway-Sant'Anna
#   3. Sun-Abraham
#   4. Gardner (did2s)
#   5. Flex TWFE (Wooldridge 2025)
#   6. Flex TWFE FGLS
#   7. Iterative GMM (efficient A)
#
# CATTs saved for: Flex_TWFE, Flex_TWFE_FGLS, GMM_Eff
#
# Key optimization: panel structure (Q_H, C_mat, lag indices) pre-computed
# once before the simulation loop; Phi regenerated per phi_sd value.
###############################################################################

library(data.table)
library(fixest)
library(did)
library(did2s)
library(Matrix)

set.seed(42)

# ===========================================================================
# 1. Parameters
# ===========================================================================

cohort_size <- 10
n_never     <- 10
N_total     <- 5 * cohort_size + n_never
T_total     <- 33
n_sims      <- 100

treatment_times <- c(10, 13, 16, 19, 22)
n_cohorts       <- length(treatment_times)
unit_cohort     <- c(rep(treatment_times, each = cohort_size), rep(0, n_never))

# ---------------------------------------------------------------------------
# Varying parameters
rho_fac_vals <- c(0.3, 0.5, 0.7)   # AR(1) persistence of common factors
phi_sd_vals  <- c(0.5, 0.9, 2.0)   # SD of unit-specific factor loadings Phi
# ---------------------------------------------------------------------------

k_fac  <- 5L
NT_all <- N_total * T_total   # 1980

# Unit-specific idiosyncratic variances (fixed across all runs)
set.seed(999)
D_unit <- runif(N_total, min = 0.5, max = 2.0)

get_cohort_size <- function(g_val) ifelse(g_val == 0, n_never, cohort_size)

# ===========================================================================
# 2. Pre-compute GMM Structure (ONCE — invariant across all parameter settings)
# ===========================================================================

cat("Pre-computing GMM structure (Q_H, C_mat, lag indices)...\n")

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
cat(sprintf("  n_catt = %d, n_did = %d\n", n_catt, n_did))

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

# Pre-compute for GMM_Eff initialisation
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

# Pre-compute time lag index vectors (for S_mat construction in GMM_Eff)
pp_v <- as.vector(abs(outer(meta_tp, meta_tp, "-")))
pr_v <- as.vector(abs(outer(meta_tp, meta_tr, "-")))
rp_v <- as.vector(abs(outer(meta_tr, meta_tp, "-")))
rr_v <- as.vector(abs(outer(meta_tr, meta_tr, "-")))

cat("  Pre-computation done.\n\n")

# ===========================================================================
# 3. DGP + True ATT / True CATTs
# ===========================================================================

generate_data <- function(beta_g_vec, r_g_vec, Phi, rho_fac) {
  unit_id <- rep(1:N_total, each = T_total)
  time_id <- rep(1:T_total, times = N_total)
  alpha <- rnorm(N_total); lambda <- rnorm(T_total)
  # AR(1) common factors (T x k_fac): each column AR(1) with persistence rho_fac
  F_time <- matrix(0, T_total, k_fac)
  F_time[1, ] <- rnorm(k_fac) / sqrt(1 - rho_fac^2)   # stationary initialisation
  for (tt in 2:T_total)
    F_time[tt, ] <- rho_fac * F_time[tt - 1, ] + rnorm(k_fac)
  # e_{i,t} = Phi[i,] . F_time[t,]  +  sqrt(D_unit[i]) * N(0,1)
  E_cross <- Phi %*% t(F_time)                                     # N x T
  E_idio  <- matrix(rnorm(NT_all), N_total, T_total) * sqrt(D_unit) # N x T
  eps     <- as.vector(t(E_cross + E_idio))                        # NT unit-major
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
# 5. Iterative GMM with Efficient A (Paper Eq. 29-31)
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
# 6. Package-based Estimators (scalar ATT only)
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
# 7. Flexible TWFE (Wooldridge 2025) — returns list: att + catt vector
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
# 8. Flex TWFE with Iterated Feasible GLS
# ===========================================================================

estimate_flex_twfe_fgls <- function(dt, max_iter = 5, tol = 1e-4) {
  tryCatch({
    dt_f <- copy(dt)
    setorder(dt_f, unit, time)

    var_per_unit_old <- NULL
    mod <- NULL

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

    if (is.null(mod)) return(list(att = NA_real_, catt = rep(NA_real_, n_catt)))

    coef_vals <- coef(mod)
    coef_keys <- as.integer(gsub("treat_gt::", "", names(coef_vals)))
    idx <- match(catt_gt_key, coef_keys)
    catt_vec <- ifelse(is.na(idx), NA_real_, coef_vals[idx])
    list(att = mean(catt_vec, na.rm = TRUE), catt = catt_vec)
  }, error = function(e) list(att = NA_real_, catt = rep(NA_real_, n_catt)))
}

# ===========================================================================
# 9. Run Simulation — saves CATTs for Flex_TWFE, Flex_TWFE_FGLS, GMM_Eff
# ===========================================================================

run_simulation <- function(beta_g_vec, r_g_vec, label, Phi, rho_fac, phi_sd) {
  true_att  <- compute_true_att(beta_g_vec, r_g_vec)
  true_catt <- compute_true_catt(beta_g_vec, r_g_vec)
  cat(sprintf("\n=== %s (rho=%.1f, phi_sd=%.1f, N=%d) ===\nTrue ATT: %.4f\n",
              label, rho_fac, phi_sd, N_total, true_att))

  est_names <- c("TWFE", "CS", "SA", "Gardner", "Flex_TWFE", "Flex_TWFE_FGLS", "GMM_Eff")
  n_est <- length(est_names)
  results <- data.table(sim = integer())
  for (nm in est_names) results[, (nm) := numeric()]

  catt_est  <- c("Flex_TWFE", "Flex_TWFE_FGLS", "GMM_Eff")
  catt_mats <- setNames(
    lapply(catt_est, function(x) matrix(NA_real_, nrow = n_sims, ncol = n_catt)),
    catt_est
  )

  na_catt <- rep(NA_real_, n_catt)

  for (s in 1:n_sims) {
    t0_sim <- proc.time()[3]
    dt    <- generate_data(beta_g_vec, r_g_vec, Phi, rho_fac)
    Delta <- compute_delta(dt)

    att_twfe    <- estimate_twfe(dt)
    att_cs      <- estimate_cs(dt)
    att_sa      <- estimate_sa(dt)
    att_gardner <- estimate_gardner(dt)

    gmm_e_res     <- tryCatch(gmm_efficient(Delta, dt),
                              error = function(e) list(att = NA_real_, catt = na_catt))
    flex_res      <- estimate_flex_twfe(dt)
    flex_fgls_res <- estimate_flex_twfe_fgls(dt)

    catt_mats[["GMM_Eff"]][s, ]        <- gmm_e_res$catt
    catt_mats[["Flex_TWFE"]][s, ]      <- flex_res$catt
    catt_mats[["Flex_TWFE_FGLS"]][s, ] <- flex_fgls_res$catt

    elapsed <- round(proc.time()[3] - t0_sim, 1)
    cat(sprintf("  Sim %d/%d  (%.1fs)\n", s, n_sims, elapsed))
    flush(stdout())

    results <- rbindlist(list(results, data.table(
      sim = s, TWFE = att_twfe, CS = att_cs, SA = att_sa,
      Gardner = att_gardner, Flex_TWFE = flex_res$att,
      Flex_TWFE_FGLS = flex_fgls_res$att, GMM_Eff = gmm_e_res$att
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
      catt_idx  = 1:n_catt,
      g         = sapply(catt_list, `[`, 1),
      t         = sapply(catt_list, `[`, 2),
      true_catt = true_catt,
      bias      = catt_bias,
      variance  = catt_var
    )
  }

  cat(sprintf("\nResults (%s):\n", label)); print(summary_dt)
  cat(sprintf("\nCATT summary (mean |bias| / mean variance across %d CATTs):\n", n_catt))
  for (nm in catt_est) {
    cs <- catt_summary[[nm]]
    cat(sprintf("  %-18s  mean|bias|=%.4f  meanVar=%.4f\n",
                nm, mean(abs(cs$bias)), mean(cs$variance)))
  }

  return(list(results = results, summary = summary_dt, true_att = true_att,
              true_catt = true_catt, catt_mats = catt_mats, catt_summary = catt_summary))
}

# ===========================================================================
# 10. Execute over full parameter grid: phi_sd x rho_fac x DGP
# ===========================================================================

cat("================================================================\n")
cat(sprintf("  N=%d (%d/cohort), T=%d, %d sims, 7 estimators, full NT×NT Sigma\n",
            N_total, cohort_size, T_total, n_sims))
cat(sprintf("  rho_fac in {%s}\n", paste(rho_fac_vals, collapse = ", ")))
cat(sprintf("  phi_sd  in {%s}\n", paste(phi_sd_vals,  collapse = ", ")))
cat("================================================================\n")

beta_hom <- rep(-5, n_cohorts);          r_hom <- rep(0, n_cohorts)
beta_het <- c(-16, -12, -10, -9, -2);   r_het <- c(0.01, 0.04, 0.08, 0.10, 0.07)

all_results <- list()

for (phi_sd in phi_sd_vals) {
  # Regenerate Phi for this phi_sd value (fixed seed for reproducibility)
  set.seed(123)
  Phi <- matrix(rnorm(N_total * k_fac, sd = phi_sd), nrow = N_total, ncol = k_fac)

  for (rho_fac in rho_fac_vals) {
    key <- sprintf("phi%.1f_rho%.1f", phi_sd, rho_fac)
    cat(sprintf("\n\n>>> Running: phi_sd=%.1f | rho=%.1f <<<\n", phi_sd, rho_fac))

    res_hom <- run_simulation(
      beta_hom, r_hom,
      sprintf("Homogeneous | rho=%.1f | phi_sd=%.1f", rho_fac, phi_sd),
      Phi, rho_fac, phi_sd
    )
    res_het <- run_simulation(
      beta_het, r_het,
      sprintf("Heterogeneous | rho=%.1f | phi_sd=%.1f", rho_fac, phi_sd),
      Phi, rho_fac, phi_sd
    )

    all_results[[key]] <- list(
      hom     = res_hom,
      het     = res_het,
      rho_fac = rho_fac,
      phi_sd  = phi_sd
    )
  }
}

# ===========================================================================
# 11. Display Summary Tables
# ===========================================================================

est_order <- c("TWFE", "CS", "SA", "Gardner", "Flex_TWFE", "Flex_TWFE_FGLS", "GMM_Eff")

cat("\n\n")
cat("============================================================\n")
cat("  TABLE 2: Bias and Variance — Full NT×NT Sigma\n")
cat(sprintf("  N=%d (%d/cohort), T=%d, %d simulations\n",
            N_total, cohort_size, T_total, n_sims))
cat("  Estimators: TWFE | CS | SA | Gardner | Flex_TWFE | Flex_TWFE_FGLS | GMM_Eff\n")
cat("============================================================\n")

for (phi_sd in phi_sd_vals) {
  cat(sprintf("\n\n========== phi_sd (SD of Phi) = %.1f ==========\n", phi_sd))

  for (rho_fac in rho_fac_vals) {
    key <- sprintf("phi%.1f_rho%.1f", phi_sd, rho_fac)
    res <- all_results[[key]]

    table_k <- merge(
      res$hom$summary[, .(Estimator, Hom_Bias = Bias, Hom_Var = Variance)],
      res$het$summary[, .(Estimator, Het_Bias = Bias, Het_Var = Variance)],
      by = "Estimator"
    )
    table_k <- table_k[match(est_order, table_k$Estimator)]

    cat(sprintf("\n  --- rho = %.1f ---\n", rho_fac))
    cat(sprintf("  %-18s  %12s  %12s  %12s  %12s\n",
                "Estimator", "Hom. Bias", "Hom. Var", "Het. Bias", "Het. Var"))
    cat("  ", paste(rep("-", 72), collapse = ""), "\n", sep = "")
    for (i in 1:nrow(table_k)) {
      cat(sprintf("  %-18s  %12.4f  %12.4f  %12.4f  %12.4f\n",
                  table_k$Estimator[i],
                  table_k$Hom_Bias[i], table_k$Hom_Var[i],
                  table_k$Het_Bias[i], table_k$Het_Var[i]))
    }
  }
}

# Aggregate long-format table (for programmatic use)
agg_rows <- list()
for (key in names(all_results)) {
  res <- all_results[[key]]
  for (dgp in c("hom", "het")) {
    sm <- res[[dgp]]$summary
    sm[, rho_fac := res$rho_fac]
    sm[, phi_sd  := res$phi_sd]
    sm[, DGP     := ifelse(dgp == "hom", "Homogeneous", "Heterogeneous")]
    agg_rows[[length(agg_rows) + 1]] <- sm
  }
}
table2_long <- rbindlist(agg_rows)
table2_long <- table2_long[, .(rho_fac, phi_sd, DGP, Estimator, Bias, Variance)]
setorder(table2_long, phi_sd, rho_fac, DGP, Estimator)

cat("\n\n--- Long-format aggregate table (table2_long) ---\n")
print(table2_long)

# Save all results
save(all_results, table2_long, est_order,
     catt_list, catt_gt_key,
     rho_fac_vals, phi_sd_vals,
     file = "simulation_full_dependent_errors_updated_results.RData")
cat("\nResults saved to simulation_full_dependent_errors_updated_results.RData\n")

# ===========================================================================
# 12. Variance Gain Plots (GMM_Eff vs Flex_TWFE)
#     Shown for each parameter combination in all_results
# ===========================================================================

library(ggplot2)

k_list <- data.frame(idx = 1:n_catt,
                     k   = sapply(catt_list, function(x) x[2] - x[1]))

for (key in names(all_results)) {
  res_het_k <- all_results[[key]]$het
  rho_k     <- all_results[[key]]$rho_fac
  phi_k     <- all_results[[key]]$phi_sd

  # Per-CATT variance comparison plot
  data_plot <- data.frame(
    rbind(
      cbind(apply(res_het_k$catt_mats$GMM_Eff,    2, var), 1),
      cbind(apply(res_het_k$catt_mats$Flex_TWFE,  2, var), 2)
    ),
    c(1:n_catt, 1:n_catt)
  )
  colnames(data_plot) <- c("Variance", "Estimator", "Index")
  data_plot[data_plot$Estimator == 1, "Estimator"] <- "GMM Efficient"
  data_plot[data_plot$Estimator == 2, "Estimator"] <- "Flexible TWFE"

  p <- ggplot(data = data_plot, aes(x = Index, y = Variance, group = Estimator)) +
    geom_line(aes(color = Estimator)) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    labs(title = sprintf("CATT Variance: rho=%.1f, phi_sd=%.1f (Heterogeneous DGP)",
                         rho_k, phi_k),
         x = "CATT index", y = "Variance")
  print(p)

  # % variance gain by event-time k
  var_k <- data.frame(k = integer(0), gain_pct = numeric(0),
                      var_gmm = numeric(0), var_flex = numeric(0), n_catts = integer(0))
  for (kk in 0:20) {
    kl <- k_list$idx[k_list$k == kk]
    if (length(kl) == 0) next
    v_gmm  <- mean(apply(res_het_k$catt_mats$GMM_Eff[, kl, drop = FALSE],  2, var))
    v_flex <- mean(apply(res_het_k$catt_mats$Flex_TWFE[, kl, drop = FALSE], 2, var))
    var_k  <- rbind(var_k, data.frame(
      k = kk, gain_pct = -100 * (v_gmm - v_flex) / v_flex,
      var_gmm = v_gmm, var_flex = v_flex, n_catts = length(kl)
    ))
  }

  plot(var_k$k, var_k$gain_pct, type = "l",
       xlab = "k (event time)",
       ylab = "% Variance Gain of GMM_Eff over Flex_TWFE",
       main = sprintf("Variance Gains: rho=%.1f, phi_sd=%.1f", rho_k, phi_k))
  abline(h = 0, lty = 2)

  diff_v <- apply(res_het_k$catt_mats$GMM_Eff,   2, var) -
            apply(res_het_k$catt_mats$Flex_TWFE,  2, var)
  den_v  <- apply(res_het_k$catt_mats$Flex_TWFE,  2, var)
  gain   <- -100 * diff_v / den_v
  cat(sprintf("\nMean %% variance gain GMM_Eff vs Flex_TWFE [rho=%.1f, phi_sd=%.1f]: %.2f%%\n",
              rho_k, phi_k, mean(gain)))
}
