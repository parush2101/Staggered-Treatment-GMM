###############################################################################
# Section 3.1 Example: Four Cohorts, Three Time Periods
# Arora & Bijani (2026), Figure 1 and Equation (7)
#
# Cohort-time mean outcomes:
#   Never-treated   (g = 0): {Y_1, Y_2, Y_3} = {100, 100, 100}
#   Always-treated  (g = 1): {Y_1, Y_2, Y_3} = {150, 155,  65}  [unidentifiable CATTs]
#   Early-treated   (g = 2): {Y_1, Y_2, Y_3} = {110, 130, 125}
#   Late-treated    (g = 3): {Y_1, Y_2, Y_3} = {140, 140, 165}
#
# True CATTs under parallel trends (verified on page 7):
#   beta_22 = 20  (g=2, t=2: impact period)
#   beta_23 = 15  (g=2, t=3: one period after impact)
#   beta_33 = 25  (g=3, t=3: impact period)
#
# DGP: Y_igt = mu_g + alpha_i + tau_gt + epsilon_it
#   mu_g      : cohort baseline intercept (100 / 110 / 140)
#   alpha_i   : unit fixed effect ~ N(0, 1)
#   tau_gt    : true CATT (0 in pre-treatment periods)
#   epsilon_it: idiosyncratic error ~ MVN(0, Sigma_nt)
#              Sigma_nt: full 1200x1200, heteroskedastic + non-zero cross-unit/time cov
#              (factor model: F_load %*% t(F_load) + diag(D_var), guaranteed PD)
#
# Cohort size: 100 units each (400 total).
# Always-treated (g=1): D_it = 1 for all t; CATTs unidentifiable (no pre-period).
#   CS and Gardner exclude g=1. Flex TWFE and GMM include g=1 in data but
#   ATT is aggregated only over identifiable cohorts g=2 and g=3.
# Estimation: Flexible TWFE (Wooldridge 2025).
# Simulation: 10 iterations; reports bias and variance of the aggregated ATT.
###############################################################################

library(data.table)
library(fixest)
library(did)
library(did2s)

set.seed(42)

# ===========================================================================
# 1. Parameters
# ===========================================================================

n_sims      <- 10
cohort_size <- 100
T_total     <- 3
N_total     <- 4L * cohort_size   # 400 (4 cohorts)

# Cohort membership:
#   g = 0 : never-treated
#   g = 1 : always-treated (D_it = 1 for all t; CATTs unidentifiable)
#   g = 2 : treated from t = 2
#   g = 3 : treated from t = 3
g_cohort <- c(rep(0L, cohort_size), rep(1L, cohort_size),
              rep(2L, cohort_size), rep(3L, cohort_size))

# Cohort baseline intercepts (match pre-treatment means / t=1 levels)
mu_g <- c("0" = 100, "1" = 150, "2" = 110, "3" = 140)

# True CATTs from the paper
true_catt <- c(beta_22 = 20, beta_23 = 15, beta_33 = 25)
true_att  <- mean(true_catt)   # equal-weighted ATT = 20

catt_keys  <- c("g2_t2", "g2_t3", "g3_t3")
catt_names <- c("beta_22", "beta_23", "beta_33")

# Fixed observation-level indices (invariant across simulations)
unit_id <- rep(seq_len(N_total), each  = T_total)
time_id <- rep(seq_len(T_total), times = N_total)
g_id    <- g_cohort[unit_id]
baseline <- mu_g[as.character(g_id)]

D_it <- as.integer( (g_id == 1L)                 |   # always treated (all periods)
                    (g_id == 2L & time_id >= 2L) |
                    (g_id == 3L & time_id >= 3L))

tau_it <- numeric(N_total * T_total)
tau_it[g_id == 2L & time_id == 2L] <- true_catt["beta_22"]
tau_it[g_id == 2L & time_id == 3L] <- true_catt["beta_23"]
tau_it[g_id == 3L & time_id == 3L] <- true_catt["beta_33"]
# Always-treated (g=1): tau encodes observed mean deviations from mu_g["1"]=150
# Observed means {150, 155, 65} => increments: t=1: 0, t=2: +5, t=3: -85
tau_it[g_id == 1L & time_id == 2L] <- 5
tau_it[g_id == 1L & time_id == 3L] <- -85

# Full NT x NT error covariance: heteroskedastic diagonal + non-zero covariance
# across ALL (i,t) and (j,s) pairs (cross-unit AND cross-time).
# Factor model: Sigma_nt = F_load %*% t(F_load) + diag(D_var)  [guaranteed PD]
# Pre-computed once; fixed across all simulation draws.
{
  NT      <- N_total * T_total   # 1200
  k_fac   <- 5L
  set.seed(123)                  # separate seed so Sigma is stable
  F_load  <- matrix(rnorm(NT * k_fac, sd = 0.4), nrow = NT, ncol = k_fac)
  D_var   <- runif(NT, min = 0.5, max = 2.0)
  Sigma_nt <- tcrossprod(F_load) + diag(D_var)
  L_chol   <- t(chol(Sigma_nt))   # lower-triangular Cholesky; epsilon = L_chol %*% z
  rm(F_load, D_var)               # free memory; only L_chol needed
}

# ===========================================================================
# 2. Data-Generating Function
# ===========================================================================

generate_data <- function() {
  alpha_i <- rnorm(N_total, mean = 0, sd = 1)   # unit FEs ~ N(0,1)

  # Full NT x NT covariance: draw epsilon ~ MVN(0, Sigma_nt) via Cholesky
  # L_chol pre-computed above; ordering is (unit1,t1),(unit1,t2),...,(unitN,tT)
  epsilon <- as.numeric(L_chol %*% rnorm(N_total * T_total))

  Y_it <- baseline + alpha_i[unit_id] + tau_it + epsilon

  dt <- data.table(unit = unit_id, time = time_id,
                   g = g_id, D = D_it, Y = Y_it)
  dt[, cohort_time := fifelse(D == 1L,
                              paste0("g", g, "_t", time),
                              "ref")]
  dt
}

# ===========================================================================
# 3. Verify DGP: cohort-time means should approximate Eq. (7)
# ===========================================================================

dt_check <- generate_data()
means_check <- dcast(
  dt_check[, .(mean_Y = round(mean(Y), 3)), by = .(g, time)],
  g ~ time, value.var = "mean_Y"
)
setnames(means_check, as.character(1:T_total), paste0("t=", 1:T_total))

cat("=======================================================\n")
cat("  DGP Verification: Cohort-Time Means (one draw)\n")
cat("=======================================================\n")
cat("  g=0: {100,100,100}  g=1: {150,155,65}\n")
cat("  g=2: {110,130,125}  g=3: {140,140,165}\n\n")
print(means_check)
cat("\n")

# ===========================================================================
# 4. Estimator Functions
# ===========================================================================

# --- Flexible TWFE (Wooldridge 2025) ---------------------------------------
estimate_flex_twfe <- function(dt) {
  tryCatch({
    mod <- feols(Y ~ i(cohort_time, ref = "ref") | unit + time,
                 data = dt, warn = FALSE, notes = FALSE)
    coef_raw        <- coef(mod)
    names(coef_raw) <- gsub("cohort_time::", "", names(coef_raw))
    mean(coef_raw[catt_keys])
  }, error = function(e) NA_real_)
}

# --- Callaway & Sant'Anna (2021) -------------------------------------------
# Uses never-treated (g = 0) as control; aggregates to simple ATT.
# Always-treated (g = 1) excluded: no pre-treatment period available, so
# att_gt cannot form a pre/post DiD for that cohort.
estimate_cs <- function(dt) {
  tryCatch({
    out <- att_gt(yname    = "Y",
                  tname    = "time",
                  idname   = "unit",
                  gname    = "g",
                  data     = as.data.frame(dt[g != 1L]),
                  control_group = "nevertreated",
                  print_details = FALSE,
                  bstrap   = FALSE,
                  cband    = FALSE)
    aggte(out, type = "simple")$overall.att
  }, error = function(e) NA_real_)
}

# --- Gardner (did2s, 2024) -------------------------------------------------
# Two-stage DiD: first stage removes unit + time FEs from untreated obs;
# second stage regresses on treatment indicator.
# Always-treated (g = 1) excluded: their unit FEs cannot be estimated in
# stage 1 (no D = 0 observations), biasing stage-2 residuals.
estimate_gardner <- function(dt) {
  tryCatch({
    dt_g <- dt[g != 1L]   # exclude always-treated before copying
    dt_g[, first_treat := fifelse(g == 0L, Inf, as.numeric(g))]
    mod <- did2s(data         = as.data.frame(dt_g),
                 yname        = "Y",
                 first_stage  = ~ 0 | unit + time,
                 second_stage = ~ i(D, ref = 0),
                 treatment    = "D",
                 cluster_var  = "unit",
                 verbose      = FALSE)
    as.numeric(coef(mod)["D::1"])
  }, error = function(e) NA_real_)
}

# --- GMM with Efficient Weighting (A = Omega_phi^{-1}) ---------------------
# Builds the 6 x 3 incidence matrix Q_H for the 3-cohort 3-period example
# (Table 1 of the paper). The 6 DiD estimates are:
#   D1: beta_22, never-treated control          -> Q_H[1,] = [ 1,  0,  0]
#   D2: beta_23, never-treated control          -> Q_H[2,] = [ 0,  1,  0]
#   D3: beta_33, never-treated control (m=2)    -> Q_H[3,] = [ 0,  0,  1]
#   D4: beta_33, never-treated control (m=1)    -> Q_H[4,] = [ 0,  0,  1]
#   D5: beta_22, not-yet-treated (g=3) control  -> Q_H[5,] = [ 1,  0,  0]
#   D6: forbidden (g=2 already treated control) -> Q_H[6,] = [ 1, -1,  1]
#        E[D6] = beta_33 - (beta_23 - beta_22); bias-corrected in Q_H.
#
# Iterated GMM (Hansen et al. 1996):
#   Step 0 : initial beta from A = I
#   While not converged:
#     (a) subtract current beta to get treatment-free Y_adj
#     (b) estimate autocovariances sigma_d from Y_adj
#     (c) build Omega_phi via structural covariance formulas (paper Eq. 31)
#     (d) update beta = (Q_H' Omega^{-1} Q_H)^{-1} Q_H' Omega^{-1} Delta
#     (e) check convergence on BOTH beta and Omega_phi; stop when both
#         change by less than tol across successive iterations.
estimate_gmm_efficient <- function(dt, max_iter = 10, tol = 1e-8) {
  tryCatch({
    # ---- Fixed metadata for the 6 DiDs in this DGP ----
    meta_fg <- c(2L, 2L, 3L, 3L, 2L, 3L)   # focal cohort
    meta_cg <- c(0L, 0L, 0L, 0L, 3L, 2L)   # control cohort
    meta_tp <- c(2L, 3L, 3L, 3L, 2L, 3L)   # t_post
    meta_tr <- c(1L, 1L, 1L, 2L, 1L, 2L)   # t_pre
    n_did   <- 6L

    # ---- Q_H: 6 x 3 incidence matrix ----
    Q_H <- matrix(c(
      1,  0,  0,
      0,  1,  0,
      0,  0,  1,
      0,  0,  1,
      1,  0,  0,
      1, -1,  1
    ), nrow = n_did, byrow = TRUE)

    # ---- Delta: 6 DiD estimates from cohort-time means ----
    cms <- dt[, .(Ybar = mean(Y)), by = .(g, time)]
    m   <- function(g_val, t_val) cms[g == g_val & time == t_val, Ybar]

    Delta <- c(
      (m(2,2) - m(2,1)) - (m(0,2) - m(0,1)),   # D1
      (m(2,3) - m(2,1)) - (m(0,3) - m(0,1)),   # D2
      (m(3,3) - m(3,1)) - (m(0,3) - m(0,1)),   # D3
      (m(3,3) - m(3,2)) - (m(0,3) - m(0,2)),   # D4
      (m(2,2) - m(2,1)) - (m(3,2) - m(3,1)),   # D5
      (m(3,3) - m(3,2)) - (m(2,3) - m(2,2))    # D6
    )

    # ---- Step 0: initial estimate with A = I ----
    beta_hat  <- as.numeric(solve(crossprod(Q_H), crossprod(Q_H, Delta)))
    Omega_old <- diag(n_did)   # placeholder for convergence check

    # ---- Lag-distance matrices (fixed across all iterations) ----
    pp_v <- abs(outer(meta_tp, meta_tp, "-"))
    pr_v <- abs(outer(meta_tp, meta_tr, "-"))
    rp_v <- abs(outer(meta_tr, meta_tp, "-"))
    rr_v <- abs(outer(meta_tr, meta_tr, "-"))

    dt_r <- copy(dt)
    setorder(dt_r, unit, time)

    iter <- 0L
    repeat {
      iter     <- iter + 1L
      beta_old <- beta_hat

      # (a) Subtract treatment effects to recover treatment-free outcomes
      dt_r[, Y_adj := Y]
      dt_r[g == 2L & time == 2L, Y_adj := Y - beta_hat[1L]]
      dt_r[g == 2L & time == 3L, Y_adj := Y - beta_hat[2L]]
      dt_r[g == 3L & time == 3L, Y_adj := Y - beta_hat[3L]]

      # (b) Estimate autocovariances sigma_d, d = 0, ..., T_total
      sigma_d <- numeric(T_total + 1L)
      for (d in 0:T_total) {
        dt_r[, Y_lag := shift(Y_adj, d, type = "lag"), by = unit]
        sigma_d[d + 1L] <- dt_r[
          !is.na(Y_lag),
          mean((Y_adj - mean(Y_adj)) * (Y_lag - mean(Y_lag)))
        ]
      }

      # (c) Build Omega_phi via structural covariance formulas (paper Eq. 31)
      #   Two DiDs covary when they share cohort observations:
      #     +cov_term / N  if focal cohorts match   (fg == fg2)
      #     +cov_term / N  if control cohorts match  (cg == cg2)
      #     -cov_term / N  if focal of s1 = control of s2 (fg == cg2)
      #     -cov_term / N  if control of s1 = focal of s2 (cg == fg2)
      Omega_phi <- matrix(0, nrow = n_did, ncol = n_did)
      for (s1 in seq_len(n_did)) {
        for (s2 in seq_len(n_did)) {
          cov_term <- (sigma_d[pp_v[s1, s2] + 1L] - sigma_d[pr_v[s1, s2] + 1L]
                     - sigma_d[rp_v[s1, s2] + 1L] + sigma_d[rr_v[s1, s2] + 1L])
          val <- 0
          if (meta_fg[s1] == meta_fg[s2]) val <- val + cov_term / cohort_size
          if (meta_cg[s1] == meta_cg[s2]) val <- val + cov_term / cohort_size
          if (meta_fg[s1] == meta_cg[s2]) val <- val - cov_term / cohort_size
          if (meta_cg[s1] == meta_fg[s2]) val <- val - cov_term / cohort_size
          Omega_phi[s1, s2] <- val
        }
      }
      Omega_phi <- (Omega_phi + t(Omega_phi)) / 2          # symmetrize
      diag(Omega_phi) <- diag(Omega_phi) + 1e-8            # ridge

      # (d) Update beta = (Q_H' Omega^{-1} Q_H)^{-1} Q_H' Omega^{-1} Delta
      OQ <- tryCatch(solve(Omega_phi, Q_H), error = function(e) NULL)
      if (is.null(OQ)) break
      OD <- solve(Omega_phi, Delta)

      beta_hat <- as.numeric(tryCatch(
        solve(crossprod(Q_H, OQ), crossprod(Q_H, OD)),
        error = function(e) beta_old
      ))

      # (e) Convergence: both beta and Omega_phi must stabilise
      beta_conv  <- max(abs(beta_hat - beta_old))
      omega_conv <- max(abs(Omega_phi - Omega_old))
      Omega_old  <- Omega_phi

      if ((beta_conv < tol && omega_conv < tol) || iter >= max_iter) break
    }

    mean(beta_hat)   # equal-weighted ATT
  }, error = function(e) NA_real_)
}

# ===========================================================================
# 5. Simulation Loop (10 iterations)
# ===========================================================================

est_names  <- c("Flex_TWFE", "CS", "Gardner", "GMM_Efficient")
att_matrix <- matrix(NA_real_, nrow = n_sims, ncol = length(est_names),
                     dimnames = list(NULL, est_names))

cat(sprintf("Running %d simulation iterations...\n", n_sims))

for (s in seq_len(n_sims)) {
  dt_s <- generate_data()

  att_matrix[s, "Flex_TWFE"]     <- estimate_flex_twfe(dt_s)
  att_matrix[s, "CS"]            <- estimate_cs(dt_s)
  att_matrix[s, "Gardner"]       <- estimate_gardner(dt_s)
  att_matrix[s, "GMM_Efficient"] <- estimate_gmm_efficient(dt_s[g != 1L])

  cat(sprintf(
    "  Iter %2d:  Flex=%.3f  CS=%.3f  Gardner=%.3f  GMM_Eff=%.3f\n",
    s,
    att_matrix[s, "Flex_TWFE"],
    att_matrix[s, "CS"],
    att_matrix[s, "Gardner"],
    att_matrix[s, "GMM_Efficient"]
  ))
}

# ===========================================================================
# 6. Summary: Bias and Variance of Aggregated ATT across Estimators
# ===========================================================================

cat("\n")
cat("==============================================================\n")
cat(sprintf("  Summary: ATT Bias & Variance (%d iters | True ATT = %.0f)\n",
            n_sims, true_att))
cat("==============================================================\n")
cat(sprintf("  %-14s  %9s  %9s  %9s  %9s\n",
            "Estimator", "Mean ATT", "Bias", "Variance", "RMSE"))
cat(paste(rep("-", 62), collapse = ""), "\n")

for (nm in est_names) {
  draws <- att_matrix[, nm]
  draws <- draws[!is.na(draws)]
  bias  <- mean(draws) - true_att
  vr    <- var(draws)
  rmse  <- sqrt(bias^2 + vr)
  cat(sprintf("  %-14s  %9.4f  %9.4f  %9.4f  %9.4f\n",
              nm, mean(draws), bias, vr, rmse))
}

cat("==============================================================\n")
