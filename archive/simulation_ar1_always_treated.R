###############################################################################
# Simulation with AR(1) Errors — Always-Treated Cohort Extension
#
# Based on simulation_ar1.R. Adds one always-treated cohort (g=1, treated
# from the first sample period) with 50 observations. Runs TWFE, CS, SA,
# Gardner, GMM-Eff, and GMM-Eff-AT with 100 simulation iterations.
#
# True ATT is defined over within-sample cohorts only (g ∈ {10,13,16,19,22}).
# The always-treated cohort affects estimators via:
#   - TWFE       : contaminated time FEs (D=1 always → zero within-unit
#                  variation in D, but outcomes inflate time FEs)
#   - CS         : always-treated skipped (no pre-period); uses never-treated
#   - SA         : cohort g=1 has only post-treatment event times; sunab may
#                  drop it or distort interaction weights
#   - Gardner    : always-treated excluded from first stage (D=1 always)
#   - GMM-Eff    : g=1 not in treatment_times → excluded from all moment
#                  conditions by construction
#   - GMM-Eff-AT : augments the base GMM system with bias-corrected moments
#                  using g=1 as control; β̃₂ = CATT(1,t) nuisance parameters
#                  (normalised to CATT(1,1)=0) are jointly estimated but only
#                  β₁ (within-sample CATTs) enters the reported ATT
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
# 4. Base GMM system: Q_H matrix and Delta vector
#
# Moments: never-treated (g=0), not-yet-treated, already-treated with bias
# correction — all from treatment_times = {10,13,16,19,22}.
# g=1 (always-treated) is NOT in treatment_times → excluded by construction.
# =============================================================================

build_gmm_system <- function(dt) {
  cohort_means <- dt[, .(Y_mean = mean(Y)), by = .(g, time)]
  setkey(cohort_means, g, time)
  get_mean <- function(g_val, t_val) {
    res <- cohort_means[.(g_val, t_val), Y_mean]
    if (length(res) == 0 || is.na(res)) return(NA_real_)
    res
  }

  treated_g <- sort(treatment_times)
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

      # (a) Never-treated (g=0)
      vals <- c(get_mean(g_c, t_post), get_mean(g_c, t_pre),
                get_mean(0,   t_post), get_mean(0,   t_pre))
      if (!anyNA(vals))
        did_estimates[[length(did_estimates) + 1]] <- list(
          delta = (vals[1]-vals[2])-(vals[3]-vals[4]),
          catt_idx=catt_idx, type="never",
          focal_g=g_c, ctrl_g=0, t_post=t_post, t_pre=t_pre)

      # (b) Not-yet-treated
      for (g_l in treated_g) {
        if (g_l - g_c <= k || g_l <= t_post) next
        vals <- c(get_mean(g_c, t_post), get_mean(g_c, t_pre),
                  get_mean(g_l, t_post), get_mean(g_l, t_pre))
        if (!anyNA(vals))
          did_estimates[[length(did_estimates) + 1]] <- list(
            delta = (vals[1]-vals[2])-(vals[3]-vals[4]),
            catt_idx=catt_idx, type="notyet",
            focal_g=g_c, ctrl_g=g_l, t_post=t_post, t_pre=t_pre)
      }

      # (c) Already-treated with bias correction (g=1 excluded as not in treated_g)
      for (g_j in treated_g) {
        if (g_j - g_c >= -m || g_j >= g_c) next
        vals <- c(get_mean(g_c, t_post), get_mean(g_c, t_pre),
                  get_mean(g_j, t_post), get_mean(g_j, t_pre))
        if (!anyNA(vals)) {
          bias_neg <- NULL; bias_pos <- NULL
          for (ci in 1:n_catt) {
            if (catt_list[[ci]][1]==g_j && catt_list[[ci]][2]==t_post) bias_neg <- ci
            if (catt_list[[ci]][1]==g_j && catt_list[[ci]][2]==t_pre)  bias_pos <- ci
          }
          did_estimates[[length(did_estimates) + 1]] <- list(
            delta = (vals[1]-vals[2])-(vals[3]-vals[4]),
            catt_idx=catt_idx, type="already",
            bias_neg=bias_neg, bias_pos=bias_pos,
            focal_g=g_c, ctrl_g=g_j, t_post=t_post, t_pre=t_pre)
        }
      }
    }
  }

  n_did <- length(did_estimates)
  if (n_did == 0) return(NULL)

  Delta <- numeric(n_did)
  Q_H   <- matrix(0, nrow=n_did, ncol=n_catt)
  for (s in 1:n_did) {
    est <- did_estimates[[s]]
    Delta[s] <- est$delta
    Q_H[s, est$catt_idx] <- 1
    if (est$type == "already") {
      if (!is.null(est$bias_neg)) Q_H[s, est$bias_neg] <- Q_H[s, est$bias_neg] - 1
      if (!is.null(est$bias_pos)) Q_H[s, est$bias_pos] <- Q_H[s, est$bias_pos] + 1
    }
  }

  list(Delta=Delta, Q_H=Q_H, n_did=n_did, n_catt=n_catt,
       catt_list=catt_list, did_estimates=did_estimates)
}

# =============================================================================
# 5. Augmented GMM system: adds bias-corrected moments using g=1 as control
#
# For each focal cohort g ∈ {10,13,16,19,22}, pre-period s ∈ {1,...,g-1},
# post-period t ∈ {g,...,T_total}, the bias-corrected moment is:
#
#   DiD_{g,1,s,t} = CATT(g,t) - [CATT(1,t) - CATT(1,s)] + ε
#
# Normalisation: CATT(1,1) = 0 (treatment period for g=1; pins the level).
# β̃₂ = { CATT(1,t) : t = 2,...,T_total } → 32 nuisance parameters.
#
# Q_H column layout: [ β₁ (n_catt cols) | β̃₂ (T_total−1 cols) ]
# Q_H row for new moment:
#   +1 on β₁[catt_idx]
#   −1 on β̃₂[t_post − 1]          (CATT(1,t_post); t_post ≥ g ≥ 10 > 1 always)
#   +1 on β̃₂[t_pre  − 1]  if t_pre ≥ 2  (else CATT(1,1)=0 → no entry)
#
# build_R_matrix works unchanged: which(unit_cohort == g_at) = units 52-101.
# =============================================================================

build_gmm_system_augmented <- function(dt) {
  base  <- build_gmm_system(dt)
  if (is.null(base)) return(NULL)

  n_catt    <- base$n_catt
  catt_list <- base$catt_list

  # β̃₂ indexing: CATT(1, t) for t = 2,...,T_total → column n_catt + (t-1)
  n_at_par  <- T_total - 1L          # 32 nuisance parameters
  at_col    <- function(t) n_catt + (t - 1L)   # t=2 → col n_catt+1, etc.

  # Always-treated cohort means
  cm_at <- dt[g == g_at, .(Y_at = mean(Y)), by = time]
  setkey(cm_at, time)
  get_at <- function(t) {
    r <- cm_at[.(t), Y_at]
    if (length(r) == 0 || is.na(r)) NA_real_ else r
  }

  # Focal cohort means
  cm_fc <- dt[, .(Y_mean = mean(Y)), by = .(g, time)]
  setkey(cm_fc, g, time)
  get_fc <- function(g_val, t_val) {
    r <- cm_fc[.(g_val, t_val), Y_mean]
    if (length(r) == 0 || is.na(r)) NA_real_ else r
  }

  new_did <- list()

  for (catt_idx in seq_len(n_catt)) {
    g_c    <- catt_list[[catt_idx]][1]
    t_post <- catt_list[[catt_idx]][2]

    for (s in 1L:(g_c - 1L)) {
      vals <- c(get_fc(g_c, t_post), get_fc(g_c, s),
                get_at(t_post),      get_at(s))
      if (!anyNA(vals))
        new_did[[length(new_did) + 1L]] <- list(
          delta    = (vals[1] - vals[2]) - (vals[3] - vals[4]),
          catt_idx = catt_idx,
          focal_g  = g_c, ctrl_g = g_at,
          t_post   = t_post, t_pre = s)
    }
  }

  n_new     <- length(new_did)
  n_did_tot <- base$n_did + n_new
  n_params  <- n_catt + n_at_par

  # Augmented Q_H: (n_did_tot) x (n_params)
  Q_aug <- matrix(0, nrow = n_did_tot, ncol = n_params)
  Q_aug[1:base$n_did, 1:n_catt] <- base$Q_H   # copy base block

  for (k in seq_len(n_new)) {
    row    <- base$n_did + k
    nd     <- new_did[[k]]
    t_post <- nd$t_post
    t_pre  <- nd$t_pre

    Q_aug[row, nd$catt_idx]    <-  1   # +1 on CATT(g, t_post)
    Q_aug[row, at_col(t_post)] <- -1   # -1 on CATT(1, t_post); t_post >= 10 > 1 always
    if (t_pre >= 2L)
      Q_aug[row, at_col(t_pre)] <- 1  # +1 on CATT(1, t_pre); if t_pre=1 → 0 by normalisation
  }

  # Augmented Delta and did_estimates (build_R_matrix uses did_estimates generically)
  list(
    Delta        = c(base$Delta, vapply(new_did, `[[`, numeric(1), "delta")),
    Q_H          = Q_aug,
    n_did        = n_did_tot,
    n_catt       = n_catt,
    n_at_par     = n_at_par,
    n_params     = n_params,
    catt_list    = catt_list,
    did_estimates = c(base$did_estimates, new_did)
  )
}

# =============================================================================
# 6. R matrix — maps unit observations to DiD estimates (shared by both GMMs)
#    For new moments ctrl_g = g_at = 1: which(unit_cohort == 1) = units 52-101
# =============================================================================

build_R_matrix <- function(gmm_sys) {
  n_did <- gmm_sys$n_did
  NT    <- N_total * T_total

  R_i <- integer(0); R_j <- integer(0); R_x <- numeric(0)

  for (s in seq_len(n_did)) {
    est    <- gmm_sys$did_estimates[[s]]
    fu     <- which(unit_cohort == est$focal_g)
    cu     <- which(unit_cohort == est$ctrl_g)
    Nf     <- length(fu); Nc <- length(cu)
    tp     <- est$t_post;  tr <- est$t_pre

    R_i <- c(R_i, rep(s, 2L * Nf), rep(s, 2L * Nc))
    R_j <- c(R_j,
             (fu - 1L) * T_total + tp, (fu - 1L) * T_total + tr,
             (cu - 1L) * T_total + tp, (cu - 1L) * T_total + tr)
    R_x <- c(R_x,
             rep( 1/Nf, Nf), rep(-1/Nf, Nf),
             rep(-1/Nc, Nc), rep( 1/Nc, Nc))
  }

  sparseMatrix(i=R_i, j=R_j, x=R_x, dims=c(n_did, NT))
}

# =============================================================================
# 7. Unit-level Omega_v and Var[Delta]
# =============================================================================

estimate_omega_v <- function(dt) {
  dt[, resid := residuals(feols(Y ~ D | unit + time, data = dt))]
  lapply(seq_len(N_total), function(i) { v <- dt[unit == i, resid]; v %o% v })
}

compute_var_delta <- function(R_mat, omega_blocks) {
  Omv  <- bdiag(lapply(omega_blocks, function(b) Matrix(b, sparse = TRUE)))
  as.matrix(R_mat %*% Omv %*% t(R_mat))
}

# =============================================================================
# 8. Estimators
# =============================================================================

estimate_twfe <- function(dt) {
  tryCatch(
    as.numeric(coef(feols(Y ~ D | unit + time, data = dt))["D"]),
    error = function(e) NA_real_)
}

estimate_cs <- function(dt) {
  tryCatch({
    out <- att_gt(yname="Y", tname="time", idname="unit", gname="g_cs",
                  data=as.data.frame(dt), control_group="nevertreated",
                  print_details=FALSE, bstrap=FALSE, cband=FALSE)
    aggte(out, type="simple", na.rm=TRUE)$overall.att
  }, error = function(e) NA_real_)
}

estimate_sa <- function(dt) {
  tryCatch({
    mod <- feols(Y ~ sunab(g_inf, time) | unit + time, data = dt)
    summary(mod, agg="ATT")$coeftable[1, 1]
  }, error = function(e) NA_real_)
}

estimate_gardner <- function(dt) {
  tryCatch({
    dt_g <- copy(dt)
    dt_g[, first_treat := fifelse(g == 0L, Inf, as.numeric(g))]
    as.numeric(coef(
      did2s(as.data.frame(dt_g), yname="Y",
            first_stage=~0|unit+time, second_stage=~i(D, ref=0),
            treatment="D", cluster_var="unit", verbose=FALSE))["D::1"])
  }, error = function(e) NA_real_)
}

# Efficient GMM — base moments only (g=1 excluded)
# ATT = cohort-size-weighted average over β₁
gmm_eff_core <- function(Q_H, Delta, n_catt, n_did, R_mat, dt) {
  QtQ  <- crossprod(Q_H)
  beta <- tryCatch(solve(QtQ, crossprod(Q_H, Delta)),
                   error = function(e) as.numeric(ginv(QtQ) %*% crossprod(Q_H, Delta)))

  for (iter in seq_len(10L)) {
    beta_old  <- beta
    Var_Delta <- compute_var_delta(R_mat, estimate_omega_v(dt))
    V_inv     <- tryCatch(solve(Var_Delta + diag(1e-6, n_did)), error=function(e) NULL)
    if (is.null(V_inv)) break
    QtVQ <- crossprod(Q_H, V_inv %*% Q_H)
    QtVD <- crossprod(Q_H, V_inv %*% Delta)
    if (rcond(QtVQ) < 1e-15) break
    beta <- tryCatch(solve(QtVQ, QtVD), error=function(e) beta_old)
    if (max(abs(beta - beta_old)) < 1e-8) break
  }
  beta
}

estimate_gmm_eff <- function(dt, gmm_sys, R_mat) {
  tryCatch({
    beta   <- gmm_eff_core(gmm_sys$Q_H, gmm_sys$Delta,
                           gmm_sys$n_catt, gmm_sys$n_did, R_mat, dt)
    w      <- rep(cohort_size, gmm_sys$n_catt) / (cohort_size * gmm_sys$n_catt)
    as.numeric(sum(w * beta))
  }, error = function(e) NA_real_)
}

# Efficient GMM — augmented moments (g=1 used as bias-corrected control)
# β = (β₁, β̃₂); only β₁ (first n_catt entries) enters the reported ATT.
# Normalisation CATT(1,1)=0 is embedded in Q_aug column structure.
estimate_gmm_eff_at <- function(dt, gmm_aug, R_aug) {
  tryCatch({
    n_catt <- gmm_aug$n_catt
    beta   <- gmm_eff_core(gmm_aug$Q_H, gmm_aug$Delta,
                           gmm_aug$n_params, gmm_aug$n_did, R_aug, dt)
    # Report ATT from β₁ only (within-sample CATTs; discard β̃₂ nuisance)
    beta1  <- beta[seq_len(n_catt)]
    w      <- rep(cohort_size, n_catt) / (cohort_size * n_catt)
    as.numeric(sum(w * beta1))
  }, error = function(e) NA_real_)
}

# =============================================================================
# 9. Run simulation
# =============================================================================

run_simulation <- function(beta_g_vec, r_g_vec, beta_at, r_at_val, label) {

  true_att <- compute_true_att(beta_g_vec, r_g_vec)
  cat(sprintf("\n=== %s (rho=%.2f) ===\n", label, rho))
  cat(sprintf("True ATT (within-sample cohorts): %.4f\n", true_att))

  results <- data.table(sim        = integer(),
                        TWFE       = numeric(),
                        CS         = numeric(),
                        SA         = numeric(),
                        Gardner    = numeric(),
                        GMM_Eff    = numeric(),
                        GMM_Eff_AT = numeric())

  for (s in seq_len(n_sims)) {
    if (s %% 25 == 0) cat(sprintf("  sim %d / %d\n", s, n_sims))

    dt      <- generate_data_ar1(beta_g_vec, r_g_vec, beta_at, r_at_val)

    # Base GMM system (g=1 excluded)
    gmm_sys <- build_gmm_system(dt)
    R_base  <- if (!is.null(gmm_sys)) build_R_matrix(gmm_sys)  else NULL

    # Augmented GMM system (g=1 used as bias-corrected control)
    gmm_aug <- if (!is.null(gmm_sys)) build_gmm_system_augmented(dt) else NULL
    R_aug   <- if (!is.null(gmm_aug)) build_R_matrix(gmm_aug)        else NULL

    results <- rbindlist(list(results, data.table(
      sim        = s,
      TWFE       = estimate_twfe(dt),
      CS         = estimate_cs(dt),
      SA         = estimate_sa(dt),
      Gardner    = estimate_gardner(dt),
      GMM_Eff    = if (!is.null(gmm_sys)) estimate_gmm_eff(dt, gmm_sys, R_base)
                   else NA_real_,
      GMM_Eff_AT = if (!is.null(gmm_aug)) estimate_gmm_eff_at(dt, gmm_aug, R_aug)
                   else NA_real_
    )))
  }

  ests <- c("TWFE","CS","SA","Gardner","GMM_Eff","GMM_Eff_AT")
  summary_dt <- data.table(
    Estimator = ests,
    Bias = sapply(ests, function(e) {
      v <- results[[e]][!is.na(results[[e]])]
      round(mean(v) - true_att, 4) }),
    Variance = sapply(ests, function(e) {
      v <- results[[e]][!is.na(results[[e]])]
      round(var(v), 4) }),
    RMSE = sapply(ests, function(e) {
      v <- results[[e]][!is.na(results[[e]])]
      round(sqrt(mean((v - true_att)^2)), 4) }),
    N_valid = sapply(ests, function(e) sum(!is.na(results[[e]])))
  )

  cat(sprintf("\nResults (%s):\n", label))
  print(summary_dt)
  return(list(results=results, summary=summary_dt, true_att=true_att))
}

# =============================================================================
# 10. Execute
# =============================================================================

cat("=================================================================\n")
cat(sprintf("  AR(1) SIMULATION — ALWAYS-TREATED COHORT (n_at=%d, g=1)\n", n_at))
cat(sprintf("  rho=%.2f  |  n_sims=%d  |  beta_at=%.1f\n",
            rho, n_sims, beta_at_val))
cat("  Estimators: TWFE, CS, SA, Gardner, GMM-Eff, GMM-Eff-AT\n")
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
                      Hom_Bias=Bias, Hom_Var=Variance, Hom_RMSE=RMSE)],
  res_het$summary[, .(Estimator,
                      Het_Bias=Bias, Het_Var=Variance, Het_RMSE=RMSE)],
  by = "Estimator"
)
table_out <- table_out[match(c("TWFE","CS","SA","Gardner","GMM_Eff","GMM_Eff_AT"),
                              table_out$Estimator)]

cat(sprintf("%-12s  %9s %9s %9s    %9s %9s %9s\n",
            "Estimator","Hom.Bias","Hom.Var","Hom.RMSE",
            "Het.Bias","Het.Var","Het.RMSE"))
cat(paste(rep("-", 74), collapse=""), "\n")
for (i in seq_len(nrow(table_out)))
  cat(sprintf("%-12s  %9.4f %9.4f %9.4f    %9.4f %9.4f %9.4f\n",
              table_out$Estimator[i],
              table_out$Hom_Bias[i], table_out$Hom_Var[i], table_out$Hom_RMSE[i],
              table_out$Het_Bias[i], table_out$Het_Var[i], table_out$Het_RMSE[i]))

save(res_hom, res_het, table_out,
     file="simulation_ar1_always_treated_results.RData")
cat("\nResults saved to simulation_ar1_always_treated_results.RData\n")

# =============================================================================
# 12. Comprehensive per-scenario summary
#
# summarize_simulation() accepts the list returned by run_simulation() and
# prints a full statistics table plus a focused efficiency comparison between
# GMM-Eff and GMM-Eff-AT.  Call it for any scenario after results exist.
#
# Statistics reported per estimator:
#   Mean      — average estimate across valid simulations
#   Bias      — Mean − True ATT
#   Std Dev   — standard deviation of estimates
#   Variance  — Var of estimates (= Bias^2 + Variance decomposition target)
#   RMSE      — root mean squared error = sqrt(Bias^2 + Variance)
#   MAE       — mean absolute error (robust to outliers)
#   N_valid   — simulations where estimator returned a non-NA value
#
# Efficiency block:
#   Var ratio  — Var(GMM-Eff) / Var(GMM-Eff-AT); >1 means AT variant is tighter
#   RMSE reduc — percentage RMSE reduction of AT variant over base GMM-Eff
# =============================================================================

summarize_simulation <- function(res, label) {
  true_att <- res$true_att
  results  <- res$results
  ests     <- c("TWFE", "CS", "SA", "Gardner", "GMM_Eff", "GMM_Eff_AT")

  cat(sprintf("\n%s\n", paste(rep("=", 74), collapse = "")))
  cat(sprintf("  DETAILED SUMMARY: %-30s  True ATT = %.4f\n", label, true_att))
  cat(sprintf("%s\n\n", paste(rep("=", 74), collapse = "")))

  stats <- data.table(
    Estimator = ests,
    Mean    = sapply(ests, function(e) round(mean(na.omit(results[[e]])), 4)),
    Bias    = sapply(ests, function(e) round(mean(na.omit(results[[e]])) - true_att, 4)),
    Std_Dev = sapply(ests, function(e) round(sd(na.omit(results[[e]])), 4)),
    Variance= sapply(ests, function(e) round(var(na.omit(results[[e]])), 4)),
    RMSE    = sapply(ests, function(e) {
                v <- na.omit(results[[e]])
                round(sqrt(mean((v - true_att)^2)), 4) }),
    MAE     = sapply(ests, function(e) {
                v <- na.omit(results[[e]])
                round(mean(abs(v - true_att)), 4) }),
    N_valid = sapply(ests, function(e) sum(!is.na(results[[e]])))
  )

  # Main table
  cat(sprintf("%-12s  %8s  %8s  %8s  %8s  %8s  %8s  %6s\n",
              "Estimator", "Mean", "Bias", "Std.Dev",
              "Variance", "RMSE", "MAE", "N"))
  cat(paste(rep("-", 80), collapse = ""), "\n")
  for (i in seq_len(nrow(stats)))
    cat(sprintf("%-12s  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %6d\n",
                stats$Estimator[i], stats$Mean[i],    stats$Bias[i],
                stats$Std_Dev[i],   stats$Variance[i], stats$RMSE[i],
                stats$MAE[i],       stats$N_valid[i]))

  # Bias-variance decomposition check: MSE = Bias^2 + Variance
  cat(sprintf("\n  MSE decomposition  (Bias^2 + Variance = RMSE^2):\n"))
  cat(sprintf("  %-12s  %8s  %8s  %8s\n", "Estimator", "Bias^2", "Variance", "RMSE^2"))
  cat(paste(rep("-", 46), collapse = ""), "\n")
  for (i in seq_len(nrow(stats)))
    cat(sprintf("  %-12s  %8.4f  %8.4f  %8.4f\n",
                stats$Estimator[i],
                round(stats$Bias[i]^2, 4),
                stats$Variance[i],
                round(stats$RMSE[i]^2, 4)))

  # Efficiency comparison: GMM-Eff-AT vs GMM-Eff
  e_base <- stats[Estimator == "GMM_Eff"]
  e_at   <- stats[Estimator == "GMM_Eff_AT"]
  if (nrow(e_base) == 1 && nrow(e_at) == 1) {
    var_ratio  <- e_base$Variance / e_at$Variance
    rmse_reduc <- 100 * (e_base$RMSE - e_at$RMSE) / e_base$RMSE
    cat(sprintf("\n  Efficiency gain — GMM-Eff-AT vs GMM-Eff:\n"))
    cat(sprintf("    Variance(GMM-Eff)  = %.4f\n", e_base$Variance))
    cat(sprintf("    Variance(GMM-Eff-AT) = %.4f\n", e_at$Variance))
    cat(sprintf("    Variance ratio (base / AT) = %.3f  ",  var_ratio))
    cat(sprintf("[>1 means AT is more efficient]\n"))
    cat(sprintf("    RMSE reduction             = %+.1f%%\n", rmse_reduc))
  }

  invisible(stats)
}

# --- Call summary for both scenarios ---
cat("\n")
sum_hom <- summarize_simulation(res_hom, "Homogeneous")
sum_het <- summarize_simulation(res_het, "Heterogeneous")
