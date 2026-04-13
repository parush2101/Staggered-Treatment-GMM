###############################################################################
# Table 2 Replication with AR(1) Dependent Errors — Large N
#
# DGP: Y_it = alpha_i + lambda_t + tau_it + epsilon_it
#   where epsilon_i ~ MVN(0, Sigma) with Sigma[t,s] = rho^|t-s|
#
# Estimators:
#   1. TWFE                    5. GMM (A = I)
#   2. Callaway-Sant'Anna      6. Flex TWFE (Wooldridge 2025)
#   3. Sun-Abraham             7. GMM (diagonal A, frequency-weighted)
#   4. Gardner (did2s)         8. Iterative GMM (efficient A)
#
# GMM_Diag: diagonal A where w_s ∝ 1/(sum of frequencies of the 4 cell means
#   used by DiD s). DiDs relying on rare cohort×time means get higher weight.
#
# CATTs saved for: GMM_I, Flex_TWFE, GMM_Diag, GMM_Eff (n_sims × n_catt)
#
# Key optimization: panel structure (Q_H, C_mat, lag indices, diagonal weights)
# is pre-computed once before the simulation loop.
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
rho         <- 0.7

treatment_times <- c(10, 13, 16, 19, 22)
n_cohorts       <- length(treatment_times)
unit_cohort     <- c(rep(treatment_times, each = cohort_size), rep(0, n_never))

Sigma_ar1 <- toeplitz(rho^(0:(T_total - 1)))

get_cohort_size <- function(g_val) ifelse(g_val == 0, n_never, cohort_size)

# ===========================================================================
# 2. Pre-compute GMM Structure (ONCE — invariant across simulations)
# ===========================================================================

cat("Pre-computing GMM structure (Q_H, C_mat, lag indices, diagonal weights)...\n")

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

# Pre-compute for GMM (A=I)
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
gg <- outer(meta_focal, meta_focal, "==")
gc_m <- outer(meta_focal, meta_ctrl,  "==")
cg <- outer(meta_ctrl,  meta_focal, "==")
cc <- outer(meta_ctrl,  meta_ctrl,  "==")
C_mat <- sweep(gg - gc_m, 1, 1/N_f, "*") + sweep(cc - cg, 1, 1/N_c, "*")
rm(gg, gc_m, cg, cc); invisible(gc())

# Pre-compute time lag index vectors (for S_mat construction)
pp_v <- as.vector(abs(outer(meta_tp, meta_tp, "-")))
pr_v <- as.vector(abs(outer(meta_tp, meta_tr, "-")))
rp_v <- as.vector(abs(outer(meta_tr, meta_tp, "-")))
rr_v <- as.vector(abs(outer(meta_tr, meta_tr, "-")))

# ---------------------------------------------------------------------------
# Diagonal A: frequency-based weights
#   Each 2x2 DiD uses 4 cell means: (focal,t_post), (focal,t_pre),
#   (ctrl,t_post), (ctrl,t_pre). Count how many DiDs each cell appears in.
#   Weight w_s = 1 / (sum of frequencies of s's 4 cells).
#   DiDs using rare cell means get upweighted; those using heavily-reused
#   means (e.g., never-treated at common times) get downweighted.
# ---------------------------------------------------------------------------
cat("  Computing diagonal frequency weights...\n")

cell_keys_fp <- paste(meta_focal, meta_tp, sep = "_")  # focal, t_post
cell_keys_fr <- paste(meta_focal, meta_tr, sep = "_")  # focal, t_pre
cell_keys_cp <- paste(meta_ctrl,  meta_tp, sep = "_")  # ctrl,  t_post
cell_keys_cr <- paste(meta_ctrl,  meta_tr, sep = "_")  # ctrl,  t_pre

freq_table <- table(c(cell_keys_fp, cell_keys_fr, cell_keys_cp, cell_keys_cr))

total_freq <- as.numeric(freq_table[cell_keys_fp]) +
              as.numeric(freq_table[cell_keys_fr]) +
              as.numeric(freq_table[cell_keys_cp]) +
              as.numeric(freq_table[cell_keys_cr])

diag_weights <- 1 / total_freq
diag_weights <- diag_weights / sum(diag_weights)  # normalize to sum to 1

# Pre-compute for GMM (diagonal A): (Q_H' diag(w) Q_H)^{-1}
QtAQ_diag     <- crossprod(Q_H, diag_weights * Q_H)
QtAQ_diag_inv <- solve(QtAQ_diag)

cat(sprintf("  Diagonal weights: min=%.6f, max=%.6f, ratio=%.1f\n",
            min(diag_weights), max(diag_weights),
            max(diag_weights) / min(diag_weights)))

cat("  Pre-computation done.\n\n")

# ===========================================================================
# 3. DGP + True ATT / True CATTs
# ===========================================================================

generate_data <- function(beta_g_vec, r_g_vec) {
  unit_id <- rep(1:N_total, each = T_total)
  time_id <- rep(1:T_total, times = N_total)
  alpha <- rnorm(N_total); lambda <- rnorm(T_total)
  eps <- as.vector(t(mvrnorm(N_total, rep(0, T_total), Sigma_ar1)))
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
# 5. GMM with A = I  (returns list: att + catt vector)
# ===========================================================================

gmm_identity <- function(Delta) {
  beta_hat <- as.numeric(QtQ_inv %*% crossprod(Q_H, Delta))
  list(att = mean(beta_hat), catt = beta_hat)
}

# ===========================================================================
# 6. GMM with Diagonal A (frequency-weighted)
#    A = diag(w) where w_s ∝ 1/(sum of cell-mean frequencies for DiD s)
#    Pre-computed: diag_weights, QtAQ_diag_inv
# ===========================================================================

gmm_diagonal <- function(Delta) {
  beta_hat <- as.numeric(QtAQ_diag_inv %*% crossprod(Q_H, diag_weights * Delta))
  list(att = mean(beta_hat), catt = beta_hat)
}

# ===========================================================================
# 7. Iterative GMM with Efficient A (Paper Eq. 29-31)
#    A = Omega_phi^{-1}, iterated from residuals
# ===========================================================================

gmm_efficient <- function(Delta, dt, max_iter = 3, tol = 1e-6) {
  beta_hat <- as.numeric(QtQ_inv %*% crossprod(Q_H, Delta))

  dt_r <- copy(dt)
  setorder(dt_r, unit, time)

  for (iter in 1:max_iter) {
    beta_old <- beta_hat

    # Reconstruct treatment effects from CATTs
    dt_r[, tau_hat := 0]
    for (ci in 1:n_catt) {
      dt_r[g == catt_list[[ci]][1] & time == catt_list[[ci]][2],
           tau_hat := beta_hat[ci]]
    }
    dt_r[, Y_adj := Y - tau_hat]

    # Residuals from FE regression on adjusted Y
    resid_mat <- matrix(residuals(feols(Y_adj ~ 1 | unit + time, data = dt_r)),
                        nrow = T_total, ncol = N_total)

    # Vectorized autocovariance estimation
    sigma_d <- numeric(T_total)
    for (d in 0:(T_total - 1)) {
      r1 <- 1:(T_total - d); r2 <- (1 + d):T_total
      sigma_d[d + 1] <- sum(resid_mat[r1, ] * resid_mat[r2, ]) / (N_total * (T_total - d))
    }

    # S_mat via vectorized sigma lookup, then Omega_phi = C_mat * S_mat
    S_vec <- sigma_d[pp_v + 1] - sigma_d[pr_v + 1] - sigma_d[rp_v + 1] + sigma_d[rr_v + 1]
    Omega_phi <- C_mat * matrix(S_vec, nrow = n_did)
    Omega_phi <- (Omega_phi + t(Omega_phi)) / 2
    diag(Omega_phi) <- diag(Omega_phi) + 1e-6

    # Solve Omega_phi X = Q_H and Omega_phi y = Delta
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
# 8. Package-based Estimators (scalar ATT only)
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
# 9. Flexible TWFE (Wooldridge 2025) — returns list: att + catt vector
# ===========================================================================

estimate_flex_twfe <- function(dt) {
  tryCatch({
    dt_f <- copy(dt)
    dt_f[, treat_gt := fifelse(D == 1, g * 100L + as.integer(time), 0L)]
    mod <- feols(Y ~ i(treat_gt, ref = 0) | unit + time, data = dt_f)
    coef_vals <- coef(mod)
    coef_keys <- as.integer(gsub("treat_gt::", "", names(coef_vals)))
    # Map to CATT ordering
    idx <- match(catt_gt_key, coef_keys)
    catt_vec <- ifelse(is.na(idx), NA_real_, coef_vals[idx])
    list(att = mean(catt_vec, na.rm = TRUE), catt = catt_vec)
  }, error = function(e) list(att = NA_real_, catt = rep(NA_real_, n_catt)))
}

# ===========================================================================
# 10. Run Simulation — saves CATTs for GMM_I, Flex_TWFE, GMM_Diag, GMM_Eff
# ===========================================================================

run_simulation <- function(beta_g_vec, r_g_vec, label) {
  true_att  <- compute_true_att(beta_g_vec, r_g_vec)
  true_catt <- compute_true_catt(beta_g_vec, r_g_vec)
  cat(sprintf("\n=== %s (rho=%.2f, N=%d) ===\nTrue ATT: %.4f\n",
              label, rho, N_total, true_att))

  est_names <- c("TWFE", "CS", "SA", "Gardner", "GMM_I",
                 "Flex_TWFE", "GMM_Diag", "GMM_Eff")
  n_est <- length(est_names)
  results <- data.table(sim = integer())
  for (nm in est_names) results[, (nm) := numeric()]

  # CATT storage matrices (n_sims × n_catt) for 4 estimators
  catt_est  <- c("GMM_I", "Flex_TWFE", "GMM_Diag", "GMM_Eff")
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

    gmm_i_res <- tryCatch(gmm_identity(Delta),
                           error = function(e) list(att = NA_real_, catt = na_catt))
    gmm_d_res <- tryCatch(gmm_diagonal(Delta),
                           error = function(e) list(att = NA_real_, catt = na_catt))
    gmm_e_res <- tryCatch(gmm_efficient(Delta, dt),
                           error = function(e) list(att = NA_real_, catt = na_catt))
    flex_res  <- estimate_flex_twfe(dt)

    # Store CATTs
    catt_mats[["GMM_I"]][s, ]     <- gmm_i_res$catt
    catt_mats[["GMM_Diag"]][s, ]  <- gmm_d_res$catt
    catt_mats[["GMM_Eff"]][s, ]   <- gmm_e_res$catt
    catt_mats[["Flex_TWFE"]][s, ] <- flex_res$catt

    elapsed <- round(proc.time()[3] - t0_sim, 1)
    cat(sprintf("  Sim %d/%d  (%.1fs)\n", s, n_sims, elapsed))
    flush(stdout())

    results <- rbindlist(list(results, data.table(
      sim = s, TWFE = att_twfe, CS = att_cs, SA = att_sa,
      Gardner = att_gardner, GMM_I = gmm_i_res$att,
      Flex_TWFE = flex_res$att, GMM_Diag = gmm_d_res$att,
      GMM_Eff = gmm_e_res$att
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
    cat(sprintf("  %-10s  mean|bias|=%.4f  meanVar=%.4f\n",
                nm, mean(abs(cs$bias)), mean(cs$variance)))
  }

  return(list(results = results, summary = summary_dt, true_att = true_att,
              true_catt = true_catt, catt_mats = catt_mats, catt_summary = catt_summary))
}

# ===========================================================================
# 11. Execute
# ===========================================================================

cat("================================================================\n")
cat(sprintf("  N=%d (%d/cohort), T=%d, rho=%.1f, %d sims, 8 estimators\n",
            N_total, cohort_size, T_total, rho, n_sims))
cat("================================================================\n")

beta_hom <- rep(-5, n_cohorts); r_hom <- rep(0, n_cohorts)
res_hom  <- run_simulation(beta_hom, r_hom, "Homogeneous")

cat("\n\n")
beta_het <- c(-16, -12, -10, -9, -2); r_het <- c(0.01, 0.04, 0.08, 0.10, 0.07)
res_het  <- run_simulation(beta_het, r_het, "Heterogeneous")
 
# ===========================================================================
# 12. Display Table
# ===========================================================================

plot(res_het$catt_summary$GMM_Eff$variance-res_het$catt_summary$Flex_TWFE$variance)

cat("\n\n")
cat("==========================================================================\n")
cat(sprintf("  TABLE 2: Bias and Variance — AR(1) (rho=%.1f), N=%d (%d/cohort)\n",
            rho, N_total, cohort_size))
cat("  GMM_I: A=I | GMM_Diag: freq-weighted diagonal | GMM_Eff: Omega_phi^{-1}\n")
cat("==========================================================================\n\n")

table2 <- merge(
  res_hom$summary[, .(Estimator, Hom_Bias = Bias, Hom_Var = Variance)],
  res_het$summary[, .(Estimator, Het_Bias = Bias, Het_Var = Variance)],
  by = "Estimator"
)
est_order <- c("TWFE", "CS", "SA", "Gardner", "GMM_I",
               "Flex_TWFE", "GMM_Diag", "GMM_Eff")
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
     catt_list, catt_gt_key, diag_weights,
     file = "simulation_dependent_errors_results.RData")
cat("\nResults (incl. CATTs) saved to simulation_dependent_errors_results.RData\n")

################## Plots ###############
library(ggplot2)
data_plot = data.frame(rbind(cbind(apply(res_het$catt_mats$GMM_Eff, 2, var),1),
                        cbind(apply(res_het$catt_mats$Flex_TWFE, 2, var),2)),
                       c(1:90,1:90))
colnames(data_plot) = c("Variance", "Estimator", "Index")
data_plot[data_plot$Estimator == 1,2] = "GMM Efficient"
data_plot[data_plot$Estimator == 2,2] = "Flexible TWFE"


ggplot(data=data_plot, aes(x = Index, y=Variance, group = Estimator)) + 
  geom_line(aes(color=Estimator))  +
  scale_color_brewer(palette="Dark2") +
  theme_minimal() + theme(legend.position="bottom",legend.title=element_blank())

k_list = data.frame(0,0)
for (i in 1:90){
  k_list[i,] = c(i,catt_list[[i]][2] - catt_list[[i]][1] )
}


## Gains in GMM vs Flex TWFE

var_k = data.frame(0,0,0,0,0)
for (i in 0:20){
kl = k_list[k_list[,2] == i,1]   ## CATTs based on k

  var_k[i+1,] = c(i,mean(apply(res_het$catt_mats$GMM_Eff[,kl],2,var)-
apply(res_het$catt_mats$Flex_TWFE[,kl],2,var)),mean(apply(res_het$catt_mats$GMM_Eff[,kl],2,var)),
          mean(apply(res_het$catt_mats$Flex_TWFE[,kl],2,var)),ncol(res_het$catt_mats$GMM_Eff[,kl])  )
}

plot(var_k[,1], -100*var_k[,2]/var_k[,4], type = "l")

var_k[,6] = -100*var_k[,2]/var_k[,4]

diff = apply(res_het$catt_mats$GMM_Eff,2,var)-apply(res_het$catt_mats$Flex_TWFE,2,var)
den = apply(res_het$catt_mats$Flex_TWFE,2,var)

gain = c()
for (i in 1:90){
 gain[i] = -100*diff[i]/den[i]
}

mean(gain)
