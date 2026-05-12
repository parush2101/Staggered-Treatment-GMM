###############################################################################
# Beck et al. (2010) Replication — Panel A (Full Sample, 49 States)
#
# Source: Beck, T., Levine, R., & Levkov, A. (2010). Big Bad Banks? The Winners
#   and Losers from Bank Deregulation in the United States. Journal of Finance.
#
# Specification (Arora & Bijani eq. 27):
#   ln(Gini)_st = alpha_s + lambda_t + beta * D_branch_st + eps_st
#
# Panel A: 49 U.S. states, 1976-2006, 1,519 state-year obs
#   Includes 12 always-treated states (deregulated before 1976)
#   Always-treated handling:
#     CS        -> g_cs  = 0    (coded as never-treated control)
#     SA / Flex -> g_inf = Inf  (coded as never-treated reference)
#     Gardner   -> first_treat = Inf (contributes to first-stage FE as control)
#     GMM       -> g_idx = 0   (plays role of never-treated control cohort)
#
# Estimators (Table 6, Panel A targets):
#   1. Pooled TWFE           -0.0213  (0.0076)
#   2. Callaway-Sant'Anna    -0.0101  (0.0078)
#   3. Sun-Abraham           -0.0354  (0.0114)
#   4. Gardner / did2s       +0.0195  (0.0067)
#   5. Flexible TWFE         -0.0165  (0.0141)
#   6. GMM Efficient (Full)  +0.0002  (0.0141)
###############################################################################

library(data.table)
library(fixest)
library(did)
library(did2s)
library(MASS)
library(Matrix)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
dt <- fread("panel_A_beck_replication (1).csv")

cat(sprintf("Panel A: %d obs, %d states, %d years\n",
            nrow(dt), dt[, uniqueN(state)], dt[, uniqueN(wrkyr)]))

first_yr <- min(dt$wrkyr)   # 1976

# Cohort variables — always-treated states have branch_reform < 1976
dt[, g_cs        := fifelse(branch_reform < first_yr, 0L,  as.integer(branch_reform))]
dt[, g_inf       := fifelse(branch_reform < first_yr, Inf, as.numeric(branch_reform))]
dt[, first_treat := fifelse(branch_reform < first_yr, Inf, as.numeric(branch_reform))]

# ===========================================================================
# 1. Pooled TWFE
# ===========================================================================
twfe <- feols(ln_gini ~ D_branch | state + wrkyr,
              data = dt, cluster = ~state)

# ===========================================================================
# 2. Callaway-Sant'Anna (2021)
#    control_group = "notyettreated"; always-treated (g_cs=0) act as controls
# ===========================================================================
cs_out <- att_gt(
  yname         = "ln_gini",
  tname         = "wrkyr",
  idname        = "statefip",
  gname         = "g_cs",
  data          = as.data.frame(dt),
  control_group = "notyettreated",
  est_method    = "dr",
  print_details = FALSE,
  bstrap        = FALSE,
  cband         = FALSE
)
cs_agg <- aggte(cs_out, type = "simple")

# ===========================================================================
# 3. Sun-Abraham (2021)
#    Always-treated coded as g_inf = Inf (never-treated reference in sunab)
# ===========================================================================
sa_mod <- feols(ln_gini ~ sunab(g_inf, wrkyr) | state + wrkyr,
                data = dt, cluster = ~state)
sa_agg <- summary(sa_mod, agg = "ATT")$coeftable[1, ]

# ===========================================================================
# 4. Gardner / did2s (2024)
#    first_stage uses only untreated obs; always-treated (first_treat=Inf)
#    contribute all periods to the first-stage FE as if never-treated
# ===========================================================================
gardner_mod <- did2s(
  data         = as.data.frame(dt),
  yname        = "ln_gini",
  first_stage  = ~ 0 | state + wrkyr,
  second_stage = ~ i(D_branch, ref = 0),
  treatment    = "D_branch",
  cluster_var  = "state",
  verbose      = FALSE
)

# ===========================================================================
# 5. Flexible TWFE (Wooldridge 2025)
#    Fully interacts treatment with cohort x time for within-sample cohorts.
#    treat_gt = g*100 + t when unit is actively treated by a within-sample
#    cohort, 0 otherwise (reference: untreated + always-treated observations).
#    ATT aggregated as equal-weighted mean; SE via delta method.
# ===========================================================================
dt[, treat_gt := fifelse(
  D_branch == 1L & !is.infinite(g_inf),
  as.integer(g_inf) * 100L + as.integer(wrkyr),
  0L
)]

flex_mod  <- feols(ln_gini ~ i(treat_gt, ref = 0) | state + wrkyr,
                   data = dt, cluster = ~state)

flex_coefs <- coef(flex_mod)
flex_vcov  <- vcov(flex_mod)
n_fc       <- length(flex_coefs)
w_flex     <- rep(1 / n_fc, n_fc)
flex_att   <- as.numeric(w_flex %*% flex_coefs)
flex_se    <- as.numeric(sqrt(t(w_flex) %*% flex_vcov %*% w_flex))

# ===========================================================================
# 6. Iterative Efficient GMM  (Arora & Bijani, Paper Eq. 29-31)
#
#   Adapted from simulation_table2_dependent_errors.R.
#   Key differences vs. simulation:
#     - Time re-indexed to 1..T_gmm (abstract), cohorts likewise
#     - g_idx = 0 for always-treated states (role of never-treated control)
#     - Variable cohort sizes (not fixed at 10)
#     - SE computed via efficient GMM formula after convergence
# ===========================================================================

cat("\n--- GMM Pre-computation ---\n")

# --- 6a. Re-index calendar time and cohorts to abstract integers 1..T_gmm ---
T_gmm    <- dt[, uniqueN(wrkyr)]      # 31
N_gmm    <- dt[, uniqueN(statefip)]   # 49
year_min <- first_yr                   # 1976

dt[, time_idx := as.integer(wrkyr - year_min + 1L)]
dt[, g_idx    := fifelse(branch_reform < year_min, 0L,
                          as.integer(branch_reform - year_min + 1L))]

# Standardise column names; sort unit-major (required for residual matrix)
dt_gmm <- dt[order(statefip, time_idx),
             .(unit = statefip, time = time_idx, Y = ln_gini, g = g_idx)]

# --- 6b. Cohort sizes ---
cohort_sz  <- dt_gmm[, .(N_g = uniqueN(unit)), by = g]
setkey(cohort_sz, g)
get_csize  <- function(g_val) cohort_sz[.(g_val), N_g]

treated_g_gmm <- sort(unique(dt_gmm$g[dt_gmm$g > 0L]))
cat(sprintf("  %d treated cohorts  |  T=%d  |  N=%d\n",
            length(treated_g_gmm), T_gmm, N_gmm))

# --- 6c. Enumerate CATTs ---
catt_list_gmm <- list()
for (g_c in treated_g_gmm)
  for (k in 0L:(T_gmm - g_c))
    catt_list_gmm[[length(catt_list_gmm) + 1L]] <- c(g_c, g_c + k)
n_catt_gmm <- length(catt_list_gmm)

# --- 6d. Enumerate DiD moments (never / not-yet / already-treated) ---
did_meta_gmm <- list()
for (catt_idx in seq_len(n_catt_gmm)) {
  g_c    <- catt_list_gmm[[catt_idx]][1]
  t_post <- catt_list_gmm[[catt_idx]][2]
  if (g_c <= 1L) next          # no pre-periods for cohort treated at t=1

  for (m in seq_len(g_c - 1L)) {
    t_pre <- g_c - m

    # Never / always-treated as control (g=0)
    did_meta_gmm[[length(did_meta_gmm) + 1L]] <- list(
      catt_idx = catt_idx, type = "never",
      focal_g = g_c, ctrl_g = 0L, t_post = t_post, t_pre = t_pre
    )
    # Not-yet-treated cohorts
    for (g_l in treated_g_gmm) {
      if (g_l <= t_post) next
      did_meta_gmm[[length(did_meta_gmm) + 1L]] <- list(
        catt_idx = catt_idx, type = "notyet",
        focal_g = g_c, ctrl_g = g_l, t_post = t_post, t_pre = t_pre
      )
    }
    # Already-treated cohorts (bias-corrected)
    for (g_j in treated_g_gmm) {
      j <- g_c - g_j
      if (j <= m || g_j >= g_c) next
      bias_neg <- NULL; bias_pos <- NULL
      for (ci in seq_len(n_catt_gmm)) {
        if (catt_list_gmm[[ci]][1] == g_j && catt_list_gmm[[ci]][2] == t_post) bias_neg <- ci
        if (catt_list_gmm[[ci]][1] == g_j && catt_list_gmm[[ci]][2] == t_pre)  bias_pos <- ci
      }
      did_meta_gmm[[length(did_meta_gmm) + 1L]] <- list(
        catt_idx = catt_idx, type = "already",
        bias_neg = bias_neg, bias_pos = bias_pos,
        focal_g = g_c, ctrl_g = g_j, t_post = t_post, t_pre = t_pre
      )
    }
  }
}
n_did_gmm <- length(did_meta_gmm)
cat(sprintf("  n_catt = %d  |  n_did = %d\n", n_catt_gmm, n_did_gmm))

# --- 6e. Build Q_H (modified incidence matrix) ---
Q_H_gmm <- matrix(0, nrow = n_did_gmm, ncol = n_catt_gmm)
for (s in seq_len(n_did_gmm)) {
  est <- did_meta_gmm[[s]]
  Q_H_gmm[s, est$catt_idx] <- 1L
  if (est$type == "already") {
    if (!is.null(est$bias_neg) && !is.na(est$bias_neg))
      Q_H_gmm[s, est$bias_neg] <- Q_H_gmm[s, est$bias_neg] - 1L
    if (!is.null(est$bias_pos) && !is.na(est$bias_pos))
      Q_H_gmm[s, est$bias_pos] <- Q_H_gmm[s, est$bias_pos] + 1L
  }
}
QtQ_gmm     <- crossprod(Q_H_gmm)
QtQ_inv_gmm <- tryCatch(solve(QtQ_gmm), error = function(e) ginv(QtQ_gmm))

# --- 6f. Metadata vectors and SPARSE C_mat (cohort-size structural factor) ---
# C_mat[s,s'] is non-zero only when moments s and s' share a cohort.
# Build by iterating over cohorts and accumulating triplets — avoids
# materialising the full n_did x n_did dense matrix (~20GB with all pre-periods).
meta_focal_gmm <- as.integer(sapply(did_meta_gmm, `[[`, "focal_g"))
meta_ctrl_gmm  <- as.integer(sapply(did_meta_gmm, `[[`, "ctrl_g"))
meta_tp_gmm    <- as.integer(sapply(did_meta_gmm, `[[`, "t_post"))
meta_tr_gmm    <- as.integer(sapply(did_meta_gmm, `[[`, "t_pre"))

moments_by_focal_gmm <- split(seq_len(n_did_gmm), meta_focal_gmm)
moments_by_ctrl_gmm  <- split(seq_len(n_did_gmm), meta_ctrl_gmm)
all_cohorts_C        <- unique(c(meta_focal_gmm, meta_ctrl_gmm))

# Pass 1: count total triplets for pre-allocation
n_sp_total <- 0L
for (g in all_cohorts_C) {
  n_ff <- length(moments_by_focal_gmm[[as.character(g)]])
  n_cc <- length(moments_by_ctrl_gmm[[as.character(g)]])
  n_sp_total <- n_sp_total + n_ff * n_ff + n_cc * n_cc + 2L * n_ff * n_cc
}
cat(sprintf("  Allocating %d sparse triplets (~%.0f MB)\n",
            n_sp_total, n_sp_total * 16 / 1e6))

sp_i <- integer(n_sp_total)
sp_j <- integer(n_sp_total)
sp_c <- numeric(n_sp_total)
pos  <- 1L

# Pass 2: fill triplets
# C_mat[s,s'] contributions per cohort g (with size N_g):
#   focal(s)==g & focal(s')==g  ->  +1/N_g
#   ctrl(s)==g  & ctrl(s')==g   ->  +1/N_g
#   focal(s)==g & ctrl(s')==g   ->  -1/N_g
#   ctrl(s)==g  & focal(s')==g  ->  -1/N_g
for (g in all_cohorts_C) {
  inv_Ng <- 1 / get_csize(g)
  ff  <- moments_by_focal_gmm[[as.character(g)]]
  cc  <- moments_by_ctrl_gmm[[as.character(g)]]
  n_ff <- length(ff); n_cc <- length(cc)

  if (n_ff > 0L) {
    n_blk <- n_ff * n_ff
    idx <- pos:(pos + n_blk - 1L)
    sp_i[idx] <- rep(ff, each = n_ff)
    sp_j[idx] <- rep(ff, times = n_ff)
    sp_c[idx] <- inv_Ng
    pos <- pos + n_blk
  }
  if (n_cc > 0L) {
    n_blk <- n_cc * n_cc
    idx <- pos:(pos + n_blk - 1L)
    sp_i[idx] <- rep(cc, each = n_cc)
    sp_j[idx] <- rep(cc, times = n_cc)
    sp_c[idx] <- inv_Ng
    pos <- pos + n_blk
  }
  if (n_ff > 0L && n_cc > 0L) {
    n_blk <- n_ff * n_cc
    idx <- pos:(pos + n_blk - 1L)          # focal(s)=g, ctrl(s')=g
    sp_i[idx] <- rep(ff, each = n_cc)
    sp_j[idx] <- rep(cc, times = n_ff)
    sp_c[idx] <- -inv_Ng
    pos <- pos + n_blk
    idx <- pos:(pos + n_blk - 1L)          # ctrl(s)=g, focal(s')=g
    sp_i[idx] <- rep(cc, each = n_ff)
    sp_j[idx] <- rep(ff, times = n_cc)
    sp_c[idx] <- -inv_Ng
    pos <- pos + n_blk
  }
}

# sparseMatrix sums duplicate (i,j) entries automatically
C_mat_sp <- sparseMatrix(i = sp_i, j = sp_j, x = sp_c,
                         dims = c(n_did_gmm, n_did_gmm))
rm(sp_i, sp_j, sp_c); invisible(gc())

# Extract consolidated non-zero positions for per-iteration S_mat lookup
C_nz    <- summary(C_mat_sp)
sp_i_nz <- C_nz$i
sp_j_nz <- C_nz$j
sp_c_nz <- C_nz$x
cat(sprintf("  C_mat: %d non-zeros (%.2f%% fill in %dx%d)\n",
            length(sp_i_nz), 100 * length(sp_i_nz) / n_did_gmm^2,
            n_did_gmm, n_did_gmm))

# --- 6g. Lag indices at non-zero C_mat positions (for S_mat lookup each iter) ---
# Only compute lags where C_mat is non-zero; avoids n_did^2-length vectors.
sp_pp_nz <- abs(meta_tp_gmm[sp_i_nz] - meta_tp_gmm[sp_j_nz]) + 1L
sp_pr_nz <- abs(meta_tp_gmm[sp_i_nz] - meta_tr_gmm[sp_j_nz]) + 1L
sp_rp_nz <- abs(meta_tr_gmm[sp_i_nz] - meta_tp_gmm[sp_j_nz]) + 1L
sp_rr_nz <- abs(meta_tr_gmm[sp_i_nz] - meta_tr_gmm[sp_j_nz]) + 1L

# --- 6h. Unit aggregation weights (in line with Callaway & Sant'Anna 2021) ---
# w[ci] = N_{g_ci} / sum_{all CATTs ci'} N_{g_ci'}
# Cohorts with more states receive proportionally more weight.
# (In the simulation all cohorts have size 10, so this equals 1/n_catt there.)
catt_cohort_n <- sapply(seq_len(n_catt_gmm),
                        function(ci) get_csize(catt_list_gmm[[ci]][1]))
w_unit_gmm    <- catt_cohort_n / sum(catt_cohort_n)

cat("  Pre-computation done.\n")

# --- 6h. Compute Delta (vector of 2x2 DiD estimates) ---
compute_delta_gmm <- function(dt_g) {
  cmeans <- dt_g[, .(Y_mean = mean(Y)), by = .(g, time)]
  setkey(cmeans, g, time)
  Delta <- numeric(n_did_gmm)
  for (s in seq_len(n_did_gmm)) {
    e <- did_meta_gmm[[s]]
    Delta[s] <- (cmeans[.(e$focal_g, e$t_post), Y_mean] -
                 cmeans[.(e$focal_g, e$t_pre),  Y_mean]) -
                (cmeans[.(e$ctrl_g,  e$t_post), Y_mean] -
                 cmeans[.(e$ctrl_g,  e$t_pre),  Y_mean])
  }
  Delta
}

# --- 6i. Iterative efficient GMM ---
gmm_efficient_beck <- function(Delta, dt_g, max_iter = 3, tol = 1e-6) {
  beta_hat  <- as.numeric(QtQ_inv_gmm %*% crossprod(Q_H_gmm, Delta))
  dt_r      <- copy(dt_g)
  setorder(dt_r, unit, time)
  Omega_phi <- NULL

  for (iter in seq_len(max_iter)) {
    beta_old <- beta_hat

    # Subtract estimated treatment effects from outcome
    dt_r[, tau_hat := 0]
    for (ci in seq_len(n_catt_gmm))
      dt_r[g == catt_list_gmm[[ci]][1] & time == catt_list_gmm[[ci]][2],
           tau_hat := beta_hat[ci]]
    dt_r[, Y_adj := Y - tau_hat]

    # Two-way FE residuals, arranged as T x N matrix
    resid_mat <- matrix(
      residuals(feols(Y_adj ~ 1 | unit + time, data = dt_r)),
      nrow = T_gmm, ncol = N_gmm
    )

    # Pooled autocovariance at each lag
    sigma_d <- numeric(T_gmm)
    for (d in 0L:(T_gmm - 1L)) {
      r1 <- 1L:(T_gmm - d); r2 <- (1L + d):T_gmm
      sigma_d[d + 1L] <- sum(resid_mat[r1, ] * resid_mat[r2, ]) /
                         (N_gmm * (T_gmm - d))
    }

    # Sparse Omega_phi = C_mat ⊙ S_mat evaluated only at non-zero positions
    s_vals    <- sigma_d[sp_pp_nz] - sigma_d[sp_pr_nz] -
                 sigma_d[sp_rp_nz] + sigma_d[sp_rr_nz]
    Omega_phi <- sparseMatrix(i = sp_i_nz, j = sp_j_nz,
                               x = sp_c_nz * s_vals,
                               dims = c(n_did_gmm, n_did_gmm))
    Omega_phi <- forceSymmetric((Omega_phi + t(Omega_phi)) / 2)
    Omega_phi <- Omega_phi + 1e-6 * Diagonal(n_did_gmm)

    OQ <- tryCatch(Matrix::solve(Omega_phi, Q_H_gmm), error = function(e) NULL)
    if (is.null(OQ)) break
    OD <- as.numeric(Matrix::solve(Omega_phi, Delta))

    beta_hat <- as.numeric(tryCatch(
      solve(as.matrix(crossprod(Q_H_gmm, OQ)),
            as.numeric(crossprod(Q_H_gmm, OD))),
      error = function(e) beta_old
    ))
    if (max(abs(beta_hat - beta_old)) < tol) break
  }

  # ATT: unit-weighted mean of CATTs (w[ci] = N_{g_ci} / total treated obs)
  att <- sum(w_unit_gmm * beta_hat)

  # SE via efficient GMM formula (Table 6 note: "(Q'Omega^{-1}Q)^{-1}")
  #   Var[theta_hat] = w' (Q'_H Omega^{-1} Q_H)^{-1} w,  w = unit weights
  se <- tryCatch({
    OmInvQ       <- Matrix::solve(Omega_phi, Q_H_gmm)
    QtOmInvQ     <- as.matrix(crossprod(Q_H_gmm, OmInvQ))
    QtOmInvQ_inv <- solve(QtOmInvQ)
    sqrt(as.numeric(t(w_unit_gmm) %*% QtOmInvQ_inv %*% w_unit_gmm))
  }, error = function(e) NA_real_)

  list(att = att, se = se, beta_hat = beta_hat)
}

# --- 6j. Run ---
cat("Computing Delta...\n")
Delta_beck <- compute_delta_gmm(dt_gmm)

cat("Running iterative efficient GMM (max 3 iterations)...\n")
gmm_res <- tryCatch(
  gmm_efficient_beck(Delta_beck, dt_gmm),
  error = function(e) {
    cat("GMM error:", conditionMessage(e), "\n")
    list(att = NA_real_, se = NA_real_, beta_hat = NULL)
  }
)
cat(sprintf("GMM ATT = %.4f  SE = %.4f\n", gmm_res$att, gmm_res$se))

# ===========================================================================
# Summary table
# ===========================================================================
paper_att <- c(-0.0213, -0.0101, -0.0354,  0.0195, -0.0165,  0.0002)
paper_se  <- c( 0.0076,  0.0078,  0.0114,  0.0067,  0.0141,  0.0141)

est_names <- c("Pooled TWFE", "CS", "Sun-Abraham", "Gardner",
               "Flex TWFE", "GMM Efficient")
att_vals  <- c(
  coef(twfe)["D_branch"],
  cs_agg$overall.att,
  sa_agg[1],
  coef(gardner_mod)["D_branch::1"],
  flex_att,
  gmm_res$att
)
se_vals <- c(
  se(twfe)["D_branch"],
  cs_agg$overall.se,
  sa_agg[2],
  se(gardner_mod)["D_branch::1"],
  flex_se,
  gmm_res$se
)

cat("\n")
cat("=================================================================\n")
cat("  Table 6 Replication — Panel A (49 states, 1,519 obs)\n")
cat("  Outcome: ln(Gini)  |  Treatment: D_branch  |  SE by state\n")
cat("=================================================================\n")
cat(sprintf("  %-14s  %9s  %9s  %9s  %9s\n",
            "Estimator", "ATT", "SE", "Paper ATT", "Paper SE"))
cat("  ", paste(rep("-", 58), collapse = ""), "\n", sep = "")
for (i in seq_along(est_names)) {
  cat(sprintf("  %-14s  %9.4f  %9.4f  %9.4f  %9.4f\n",
              est_names[i], att_vals[i], se_vals[i],
              paper_att[i], paper_se[i]))
}
cat("\n")
