###############################################################################
# Beck (2010) Application: Effect of Branch Banking Deregulation on Inequality
#
# Data:   panel_A_beck_replication.csv
#         49 US states x 31 years (1976-2006)
#         Outcome:   ln_gini (log Gini coefficient)
#         Treatment: D_branch (branch deregulation indicator)
#         Cohort:    branch_reform (year of deregulation)
#
# Estimators:
#   1. TWFE (vanilla)              — single scalar on D_branch
#   2. TWFE Flexible               — cohort x year cell dummies (interacted TWFE)
#   3. Gardner (did2s)             — two-stage DiD, scalar ATT
#   4. Callaway-Sant'Anna (2021)   — att_gt with not-yet-treated control
#   5. Sun-Abraham (2021)          — interaction-weighted estimator via sunab()
#   6. Efficient GMM               — iterated optimal weighting (Arora-Bijani 2025)
#
# Reference: Beck, T., Levine, R., & Levkov, A. (2010). Big Bad Banks?
#   The Winners and Losers from Bank Deregulation in the United States.
#   Journal of Finance, 65(5), 1637-1667.
###############################################################################

library(data.table)
library(fixest)
library(did2s)
library(did)
library(MASS)

# ===========================================================================
# 1. Load Data
# ===========================================================================

dt <- fread("panel_A_beck_replication.csv")

# 49 states x 31 years (1976-2006); no never-treated states
# branch_reform: year of deregulation (some pre-sample: 1960, 1970, 1975, 1976)
# Within-sample treated cohorts: branch_reform in {1977, ..., 1999}

cat("States:", uniqueN(dt$state), "\n")
cat("Years: ", paste(range(dt$wrkyr), collapse = "-"), "\n")
cat("Within-sample treated cohorts:",
    paste(sort(unique(dt[branch_reform >= 1977, branch_reform])), collapse = " "), "\n\n")

# ===========================================================================
# 2. TWFE (Vanilla)
# ===========================================================================
# Single scalar coefficient on D_branch; unit + time FE.
# Aggregates CATTs using FWL-residualized weights — can be negative when
# already-treated units act as implicit controls.

m_twfe <- feols(ln_gini ~ D_branch | state + wrkyr, data = dt, cluster = ~state)

att_twfe <- as.numeric(coef(m_twfe)["D_branch"])
se_twfe  <- as.numeric(se(m_twfe)["D_branch"])

# ===========================================================================
# 3. TWFE Flexible (Interacted TWFE)
# ===========================================================================
# Separate CATT per (cohort, post-period) cell via cohort x year interactions.
# ATT reported as simple (cell-weighted) mean over all identified cells.
# Restricted to within-sample treated cohorts (branch_reform >= 1977).

dt_tw <- dt[branch_reform >= 1977]
dt_tw[, cell_id := fifelse(wrkyr >= branch_reform,
                            paste0(branch_reform, "_", wrkyr), "pre")]

m_flex <- feols(ln_gini ~ i(cell_id, ref = "pre") | state + wrkyr,
                data = dt_tw, cluster = ~state)

fc      <- coef(m_flex)
fn      <- names(fc)
keep    <- grepl("cell_id::", fn, fixed = TRUE)
cell_str <- sub("cell_id::", "", fn[keep])
fc_g    <- as.integer(sub("_.*", "", cell_str))
fc_t    <- as.integer(sub(".*_", "", cell_str))
post    <- fc_t >= fc_g

post_names    <- fn[keep][post]
cell_catts    <- fc[post_names]
n_cells       <- sum(post)
att_flex_cell <- mean(cell_catts)

# SE via delta method: ATT = w'beta, w = 1/n => Var(ATT) = w' V w
# V is the full clustered VCV of post-treatment cell CATTs
V_post   <- vcov(m_flex)[post_names, post_names]
w        <- rep(1 / n_cells, n_cells)
se_flex  <- as.numeric(sqrt(t(w) %*% V_post %*% w))

# ===========================================================================
# 4. Gardner (2021) Two-Stage DiD
# ===========================================================================
# Stage 1: estimate unit + time FE on never/not-yet-treated observations only.
# Stage 2: regress residualised outcome on D_branch to recover ATT.
# Corrects for contamination bias from already-treated implicit controls.

m_gard <- did2s(data         = dt,
                yname        = "ln_gini",
                first_stage  = ~0 | state + wrkyr,
                second_stage = ~D_branch,
                treatment    = "D_branch",
                cluster_var  = "state")

att_gard <- as.numeric(coef(m_gard)["D_branch"])
se_gard  <- as.numeric(se(m_gard)["D_branch"])

# ===========================================================================
# 5. Callaway-Sant'Anna (2021)
# ===========================================================================
# att_gt estimates CATT(g,t) for each cohort-period pair using not-yet-treated
# units as the comparison group. Pre-sample-treated states (branch_reform < 1977)
# are recoded as 0 (never-treated) so CS treats them as always-controls.
# ATT aggregated via simple average over all (g,t) post-treatment cells.

dt_cs <- copy(dt)
dt_cs[branch_reform < 1977, branch_reform := 0L]

cs_out <- att_gt(yname         = "ln_gini",
                 tname         = "wrkyr",
                 idname        = "statefip",
                 gname         = "branch_reform",
                 control_group = "notyettreated",
                 data          = as.data.frame(dt_cs),
                 est_method    = "reg",
                 print_details = FALSE)

cs_agg   <- aggte(cs_out, type = "simple")
att_cs   <- cs_agg$overall.att
se_cs    <- cs_agg$overall.se

# ===========================================================================
# 6. Sun-Abraham (2021)
# ===========================================================================
# Interaction-weighted estimator via fixest::sunab(). Pre-sample-treated
# states recoded to 0 (used as never-treated base group).
# ATT aggregated as the standard interaction-weighted average across cohorts.

dt[, g_sa := fifelse(branch_reform >= 1977 & branch_reform <= 2006,
                      branch_reform, 0L)]

m_sa    <- feols(ln_gini ~ sunab(g_sa, wrkyr) | state + wrkyr,
                 data = dt, cluster = ~state)
att_sa  <- as.numeric(aggregate(m_sa, agg = "ATT")[1])
se_sa   <- as.numeric(aggregate(m_sa, agg = "ATT")[2])

# ===========================================================================
# 7. Efficient GMM (Arora-Bijani 2025)
# ===========================================================================
# Adapts the simulation estimator from sim4_all_estimators.R to real data.
# Identification: not-yet-treated comparisons only (no never-treated exists).
# Identified CATTs: t_post <= 1998 (last year with a not-yet-treated control).
# Omega estimated iteratively from autocovariances of two-way-demeaned residuals.
# ATT = simple cell mean of all identified CATTs (cell-weighted).

T_yr   <- 1976:2006          # sample years
T_len  <- length(T_yr)       # 31
yr_idx <- function(y) y - 1975L   # maps year -> 1-based column index

# Focal cohorts and identifiable post-periods
treated_g  <- sort(unique(dt[branch_reform >= 1977, branch_reform]))
T_end_id   <- max(treated_g) - 1L   # 1998

# State-level cohort mapping (49 states in order of statefip)
states_ord  <- dt[, .(statefip = first(statefip),
                       cohort   = first(branch_reform)), by = state]
setorder(states_ord, statefip)
N_st        <- nrow(states_ord)    # 49
unit_cohort <- states_ord$cohort   # cohort year for each state (index 1:49)
cohort_n    <- dt[branch_reform %in% treated_g,
                   .(n = uniqueN(state)), by = branch_reform]
setkey(cohort_n, branch_reform)
get_N <- function(g) cohort_n[.(g), n]

# CATT list: (g_c, t_post) with t_post in [g_c, T_end_id]
catt_list_g <- list()
for (g_c in treated_g) {
  if (g_c > T_end_id) next
  for (tp in g_c:T_end_id)
    catt_list_g[[length(catt_list_g) + 1]] <- c(g_c, tp)
}
n_catt_g   <- length(catt_list_g)   # 236
catt_g_vec <- sapply(catt_list_g, `[`, 1)
catt_t_vec <- sapply(catt_list_g, `[`, 2)

# CATT index lookup
catt_key_g <- setNames(seq_len(n_catt_g),
  paste(catt_g_vec, catt_t_vec, sep = "-"))

# Build DiD meta (not-yet-treated only, all pre-periods)
did_meta_g <- vector("list", 8000L); cnt_g <- 0L
for (ci in seq_len(n_catt_g)) {
  g_c  <- catt_list_g[[ci]][1]; t_post <- catt_list_g[[ci]][2]
  notyet <- treated_g[treated_g > t_post]
  if (!length(notyet)) next
  for (t_pre in 1976L:(g_c - 1L)) {
    for (g_l in notyet) {
      cnt_g <- cnt_g + 1L
      did_meta_g[[cnt_g]] <- list(catt_idx = ci, type = "notyet",
                                   focal_g = g_c, ctrl_g = g_l,
                                   t_post = t_post, t_pre = t_pre)
    }
  }
}
did_meta_g <- did_meta_g[seq_len(cnt_g)]
n_did_g    <- cnt_g   # 5829

meta_fg <- sapply(did_meta_g, `[[`, "focal_g")
meta_cg <- sapply(did_meta_g, `[[`, "ctrl_g")
meta_tp <- sapply(did_meta_g, `[[`, "t_post")
meta_tr <- sapply(did_meta_g, `[[`, "t_pre")

# Q_H matrix (n_did x n_catt), one 1 per row (notyet: no bias corrections)
Q_H_g   <- matrix(0L, nrow = n_did_g, ncol = n_catt_g)
for (s in seq_len(n_did_g)) Q_H_g[s, did_meta_g[[s]]$catt_idx] <- 1L
Q_H_g   <- matrix(as.numeric(Q_H_g), nrow = n_did_g)
QtQ_g   <- crossprod(Q_H_g)
QtQ_inv_g <- tryCatch(solve(QtQ_g), error = function(e) ginv(QtQ_g))

# C_mat: cohort-pair covariance structure
N_f_g <- sapply(meta_fg, get_N)
N_c_g <- sapply(meta_cg, get_N)
gg_m  <- outer(meta_fg, meta_fg, "==")
gc_m  <- outer(meta_fg, meta_cg, "==")
cg_m  <- outer(meta_cg, meta_fg, "==")
cc_m  <- outer(meta_cg, meta_cg, "==")
C_mat_g <- sweep(gg_m - gc_m, 1L, 1 / N_f_g, "*") +
           sweep(cc_m - cg_m, 1L, 1 / N_c_g, "*")
rm(gg_m, gc_m, cg_m, cc_m); invisible(gc())

# Lag-index vectors for autocovariance lookup
pp_v_g <- as.vector(abs(outer(meta_tp, meta_tp, "-")))
pr_v_g <- as.vector(abs(outer(meta_tp, meta_tr, "-")))
rp_v_g <- as.vector(abs(outer(meta_tr, meta_tp, "-")))
rr_v_g <- as.vector(abs(outer(meta_tr, meta_tr, "-")))

# Build Y matrix (N_st x T_len) from ln_gini
Y_mat_g <- matrix(NA_real_, nrow = N_st, ncol = T_len)
for (i in seq_len(N_st)) {
  st <- states_ord$state[i]
  rows_i <- dt[state == st][order(wrkyr)]
  Y_mat_g[i, ] <- rows_i$ln_gini
}

# Cohort means Ybar (n_cohorts x T_len)
all_cohorts_g <- treated_g
Ybar_g <- matrix(NA_real_, nrow = length(all_cohorts_g), ncol = T_len)
for (ci in seq_along(all_cohorts_g)) {
  units <- which(unit_cohort == all_cohorts_g[ci])
  Ybar_g[ci, ] <- colMeans(Y_mat_g[units, , drop = FALSE])
}
coh_idx_g <- function(g_val) which(all_cohorts_g == g_val)

# Delta vector
Delta_g <- numeric(n_did_g)
for (s in seq_len(n_did_g)) {
  e  <- did_meta_g[[s]]
  gi <- coh_idx_g(e$focal_g); ci <- coh_idx_g(e$ctrl_g)
  Delta_g[s] <- (Ybar_g[gi, yr_idx(e$t_post)] - Ybar_g[gi, yr_idx(e$t_pre)]) -
                (Ybar_g[ci, yr_idx(e$t_post)] - Ybar_g[ci, yr_idx(e$t_pre)])
}

# Iterative efficient GMM (3 iterations)
cat("Running Efficient GMM...\n"); flush(stdout())

# Identity-weighted OLS GMM: beta = (Q'Q)^{-1} Q' Delta
beta_ols_gmm <- as.numeric(QtQ_inv_g %*% crossprod(Q_H_g, Delta_g))
att_ols_gmm  <- mean(beta_ols_gmm)
beta_gmm     <- beta_ols_gmm
Omega        <- NULL

for (iter in 1:3) {
  beta_old <- beta_gmm

  # Subtract estimated treatment effects from Y
  Y_adj <- Y_mat_g
  for (k in seq_len(n_catt_g)) {
    units_k <- which(unit_cohort == catt_g_vec[k])
    t_col   <- yr_idx(catt_t_vec[k])
    Y_adj[units_k, t_col] <- Y_adj[units_k, t_col] - beta_gmm[k]
  }

  # Two-way demean (remove unit + time FE)
  row_m <- rowMeans(Y_adj); col_m <- colMeans(Y_adj); grand_m <- mean(Y_adj)
  resid_mat <- Y_adj - outer(row_m, rep(1, T_len)) -
               outer(rep(1, N_st), col_m) + grand_m

  # Autocovariances sigma(d) for d = 0, ..., T-1
  sigma_d <- numeric(T_len)
  for (d in 0:(T_len - 1)) {
    r1 <- seq_len(T_len - d); r2 <- r1 + d
    sigma_d[d + 1] <- sum(resid_mat[, r1] * resid_mat[, r2]) / (N_st * (T_len - d))
  }

  # Omega = C_mat * S_vec (element-wise)
  S_vec <- sigma_d[pp_v_g + 1] - sigma_d[pr_v_g + 1] -
           sigma_d[rp_v_g + 1] + sigma_d[rr_v_g + 1]
  Omega <- C_mat_g * matrix(S_vec, nrow = n_did_g)
  Omega <- (Omega + t(Omega)) / 2
  diag(Omega) <- diag(Omega) + 1e-6

  OQ <- tryCatch(solve(Omega, Q_H_g), error = function(e) NULL)
  if (is.null(OQ)) { cat("  Omega solve failed at iter", iter, "\n"); break }
  OD       <- solve(Omega, Delta_g)
  beta_gmm <- as.numeric(tryCatch(
    solve(crossprod(Q_H_g, OQ), crossprod(Q_H_g, OD)),
    error = function(e) beta_old))

  conv <- max(abs(beta_gmm - beta_old))
  cat(sprintf("  Iter %d: max|delta_beta| = %.2e\n", iter, conv))
  if (conv < 1e-6) break
}

att_gmm <- mean(beta_gmm)
cat(sprintf("Efficient GMM: ATT = %.4f\n\n", att_gmm))

# ---- Standard errors (delta method: Var(ATT) = w' V_beta w, w = 1/n) ----
w_g <- rep(1 / n_catt_g, n_catt_g)

# OLS GMM sandwich: V_ols = (Q'Q)^{-1} Q' Omega Q (Q'Q)^{-1}
V_ols_gmm <- QtQ_inv_g %*% crossprod(Q_H_g, Omega %*% Q_H_g) %*% QtQ_inv_g
se_ols_gmm <- as.numeric(sqrt(t(w_g) %*% V_ols_gmm %*% w_g))

# Efficient GMM: V_eff = (Q' Omega^{-1} Q)^{-1}
OQ_fin   <- solve(Omega, Q_H_g)
V_eff    <- solve(crossprod(Q_H_g, OQ_fin))
se_gmm   <- as.numeric(sqrt(t(w_g) %*% V_eff %*% w_g))

# ===========================================================================
# 8. Results Table
# ===========================================================================

sig <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", "")))
pval <- function(est, se) 2 * pnorm(-abs(est / se))

p_twfe    <- pval(att_twfe,      se_twfe)
p_flex    <- pval(att_flex_cell, se_flex)
p_gard    <- pval(att_gard,      se_gard)
p_cs      <- pval(att_cs,        se_cs)
p_sa      <- pval(att_sa,        se_sa)
p_ols_gmm <- pval(att_ols_gmm,   se_ols_gmm)
p_gmm     <- pval(att_gmm,       se_gmm)

fmt <- "  %-30s  %8.4f  %8.4f  %7.4f  %s\n"

cat(strrep("=", 76), "\n")
cat("  Beck (2010): Effect of Branch Banking Deregulation on ln(Gini)\n")
cat("  N = 49 states x 31 years (1976-2006)\n")
cat(strrep("=", 76), "\n")
cat(sprintf("  %-30s  %8s  %8s  %7s\n", "Estimator", "ATT", "SE", "p"))
cat(strrep("-", 60), "\n")
cat(sprintf(fmt, "TWFE (vanilla)",            att_twfe,      se_twfe,    p_twfe,    sig(p_twfe)))
cat(sprintf(fmt, "Sun-Abraham (2021)",        att_sa,        se_sa,      p_sa,      sig(p_sa)))
cat(sprintf(fmt, "Callaway-Sant'Anna (2021)", att_cs,        se_cs,      p_cs,      sig(p_cs)))
cat(sprintf(fmt, "Gardner (2021)",            att_gard,      se_gard,    p_gard,    sig(p_gard)))
cat(sprintf(fmt, sprintf("TWFE Flexible (%d cells)", n_cells),
            att_flex_cell, se_flex,    p_flex,    sig(p_flex)))
cat(strrep("-", 60), "\n")
cat(sprintf(fmt, sprintf("GMM-OLS (%d cells)", n_catt_g),
            att_ols_gmm,   se_ols_gmm, p_ols_gmm, sig(p_ols_gmm)))
cat(sprintf(fmt, sprintf("GMM-Efficient (%d cells)", n_catt_g),
            att_gmm,       se_gmm,     p_gmm,     sig(p_gmm)))
cat(strrep("=", 76), "\n")
cat("  SE method — TWFE Flex: delta, V clustered by state\n")
cat("             GMM-OLS:   delta, sandwich (Q'Q)^{-1} Q'OmegaQ (Q'Q)^{-1}\n")
cat("             GMM-Eff:   delta, V = (Q'Omega^{-1}Q)^{-1}\n")
cat("  CS: not-yet-treated control; SA: interaction-weighted\n")
cat("  GMM: notyet comparisons, all pre-periods, Omega from autocovariances\n\n")

cat("Cell CATT distribution (TWFE Flexible):\n")
cat(sprintf("  Min=%.4f  p25=%.4f  Median=%.4f  p75=%.4f  Max=%.4f\n\n",
            min(cell_catts), quantile(cell_catts, 0.25),
            median(cell_catts), quantile(cell_catts, 0.75),
            max(cell_catts)))
