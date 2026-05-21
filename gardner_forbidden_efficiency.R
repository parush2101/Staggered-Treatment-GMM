###############################################################################
# gardner_forbidden_efficiency.R
#
# Purpose:
#   Show that "extended Gardner" — Gardner's GMM form but with forbidden
#   comparisons included (debiased via Q_H) — is more efficient than the
#   original Gardner estimator, even under spherical (iid) errors.
#
# Three estimators compared for each CATT:
#   1. Gardner_pkg  : Gardner (2024) via did2s package
#   2. Gardner_W    : Gardner re-expressed as GMM with A = W'(WW')^{-1}W
#                     (A has zero rows for forbidden comparisons — identical
#                     to Gardner_pkg up to numerical precision)
#   3. Ext_Gardner  : Same GMM form but with A = Q_H'Q_H (under spherical
#                     errors, Omega_phi^{-1} ∝ I, so optimal A ∝ I; here we
#                     use the OLS normal equations over ALL DiDs including
#                     forbidden ones debiased by Q_H)
#
# DGP:
#   Y_it = alpha_i + lambda_t + tau_it + eps_it
#   eps_it ~ iid N(0, sigma^2)   [spherical errors]
#   Heterogeneous treatment effects: tau_it = beta_g * (1 + r_g)^{t - g}
#
# Design:
#   - 3 treated cohorts (g = 4, 7, 10) + 1 never-treated cohort
#   - 12 time periods
#   - 10 units per cohort, 10 never-treated units
#   - Variance across simulation draws reported for each CATT
###############################################################################

library(data.table)
library(fixest)
library(did2s)   # for Gardner_pkg

set.seed(2025)

# ===========================================================================
# 1. DGP Parameters
# ===========================================================================

T_total       <- 12
cohort_size   <- 10
n_never       <- 10
treat_times   <- c(4L, 7L, 10L)
n_cohorts     <- length(treat_times)
N_total       <- n_cohorts * cohort_size + n_never
sigma_eps     <- 1.0          # spherical errors

# Heterogeneous treatment effects
beta_g <- c(g4 = 3.0, g7 = -2.0, g10 = 5.0)
r_g    <- c(g4 = 0.10, g7 = 0.05, g10 = 0.00)

unit_cohort <- c(rep(treat_times, each = cohort_size), rep(0L, n_never))
unit_ids    <- seq_len(N_total)

# True CATTs: beta_g * (1 + r_g)^{t - g}  for t >= g, else 0
true_catt <- function(g, t) {
  if (g == 0L || t < g) return(0.0)
  gchar <- paste0("g", g)
  beta_g[gchar] * (1 + r_g[gchar])^(t - g)
}

# ===========================================================================
# 2. Pre-compute GMM structure: Q_H and DiD meta-data
# ===========================================================================

# --- 2a. CATT list ---
catt_list <- list()
for (g in treat_times) {
  for (t in g:T_total) {
    catt_list[[length(catt_list) + 1]] <- c(g, t)
  }
}
n_catt <- length(catt_list)

# Helper: CATT index given (g, t)
catt_index <- function(g, t) {
  for (k in seq_len(n_catt)) {
    if (catt_list[[k]][1] == g && catt_list[[k]][2] == t) return(k)
  }
  return(NA_integer_)
}

# --- 2b. Enumerate all 2x2 DiDs (clean + forbidden) ---
# Each DiD: (focal_g, ctrl_g, t_post, t_pre)
# Clean: ctrl_g == 0 (never-treated) or ctrl_g > t_post (not-yet-treated)
# Forbidden: ctrl_g < focal_g and ctrl_g <= t_post (already treated)

did_list <- list()

for (focal_g in treat_times) {
  for (t_post in focal_g:T_total) {
    catt_idx <- catt_index(focal_g, t_post)

    for (m in 1:(focal_g - 1L)) {
      t_pre <- focal_g - m

      # Never-treated control (always clean)
      did_list[[length(did_list) + 1]] <- list(
        catt_idx = catt_idx,
        focal_g  = focal_g,
        ctrl_g   = 0L,
        t_post   = t_post,
        t_pre    = t_pre,
        type     = "clean_never"
      )

      # Not-yet-treated controls (clean)
      for (ctrl_g in treat_times) {
        if (ctrl_g <= t_post) next   # already treated at t_post → skip clean
        did_list[[length(did_list) + 1]] <- list(
          catt_idx = catt_idx,
          focal_g  = focal_g,
          ctrl_g   = ctrl_g,
          t_post   = t_post,
          t_pre    = t_pre,
          type     = "clean_notyet"
        )
      }

      # Already-treated controls (forbidden)
      for (ctrl_g in treat_times) {
        if (ctrl_g >= focal_g) next       # same or later cohort
        if (ctrl_g > t_post)  next        # not yet treated → handled above
        # ctrl_g is treated before focal_g and already treated at t_post
        did_list[[length(did_list) + 1]] <- list(
          catt_idx = catt_idx,
          focal_g  = focal_g,
          ctrl_g   = ctrl_g,
          t_post   = t_post,
          t_pre    = t_pre,
          type     = "forbidden"
        )
      }
    }
  }
}

n_did <- length(did_list)
cat(sprintf("Total DiDs: %d  (inspecting types...)\n", n_did))
types <- sapply(did_list, `[[`, "type")
cat(sprintf("  clean_never=%d  clean_notyet=%d  forbidden=%d\n",
            sum(types == "clean_never"),
            sum(types == "clean_notyet"),
            sum(types == "forbidden")))

# --- 2c. Build Q_H (n_did x n_catt) ---
# Clean row s:   Q_H[s, catt_idx(s)] = 1
# Forbidden row: Q_H[s, catt_idx(s)] = 1
#                Q_H[s, catt_idx(ctrl_g, t_post)] = -1   (bias term +)
#                Q_H[s, catt_idx(ctrl_g, t_pre)]  = +1   (bias term -)
# i.e. E[Delta_s] = beta_{focal_g,t_post}
#                  - beta_{ctrl_g,t_post} + beta_{ctrl_g,t_pre}
# so Q_H row maps Delta_s -> beta_{focal_g,t_post} after bias correction

Q_H <- matrix(0.0, nrow = n_did, ncol = n_catt)

for (s in seq_len(n_did)) {
  d <- did_list[[s]]
  Q_H[s, d$catt_idx] <- 1.0
  if (d$type == "forbidden") {
    neg_idx <- catt_index(d$ctrl_g, d$t_post)
    pos_idx <- catt_index(d$ctrl_g, d$t_pre)
    if (!is.na(neg_idx)) Q_H[s, neg_idx] <- Q_H[s, neg_idx] - 1.0
    if (!is.na(pos_idx)) Q_H[s, pos_idx] <- Q_H[s, pos_idx] + 1.0
  }
}

# ===========================================================================
# 3. Build Gardner's W matrix (clean DiDs only)
# ===========================================================================
# Gardner's first-stage joint FE regression on untreated cells yields
# hat_beta^Gard = W * Delta.
# W is derived from the normal equations of the stage-1 OLS.
# For the simple DGP here we build W analytically:
#   - Identify which time periods t appear in the untreated pool for each
#     cohort-period cell (focal_g, t_post).
#   - Solve for the implicit weights using the FE normal equations.
#
# Shortcut used here: since Gardner = OLS of Delta on Q_H restricted to
# clean rows, we compute W as the OLS estimator on clean rows only.

clean_idx <- which(types != "forbidden")
Q_clean   <- Q_H[clean_idx, , drop = FALSE]   # n_clean x n_catt

# W = (Q_clean' Q_clean)^{-1} Q_clean'   [n_catt x n_clean]
# so hat_beta^Gard = W * Delta_clean
QcQc      <- t(Q_clean) %*% Q_clean
W_mat     <- solve(QcQc) %*% t(Q_clean)     # n_catt x n_clean

# Sanity check: W Q_clean' should be I
stopifnot(max(abs(W_mat %*% Q_clean - diag(n_catt))) < 1e-8)

# A matrix for Gardner_W:  A_Gard = W' (WW')^{-1} W  (n_clean x n_clean)
WW_inv    <- solve(W_mat %*% t(W_mat))
A_gard    <- t(W_mat) %*% WW_inv %*% W_mat   # n_clean x n_clean

# ===========================================================================
# 4. Extended Gardner: A = I over ALL DiDs (optimal under spherical errors)
# ===========================================================================
# Ext_Gardner uses GMM with A = I_{n_did}:
#   hat_beta = (Q_H' Q_H)^{-1} Q_H' Delta
# This is OLS of Delta (all DiDs, clean + forbidden) on Q_H.

QhQh      <- t(Q_H) %*% Q_H
QhQh_inv  <- solve(QhQh)   # n_catt x n_catt

# ===========================================================================
# 5. Helper: compute cohort-period mean outcomes from panel data
# ===========================================================================

# Given data.table with columns (id, time, y, cohort),
# return Delta vector (length n_did).
compute_Delta <- function(dt) {
  gp_means <- dt[, .(ybar = mean(y)), by = .(cohort, time)]
  setkey(gp_means, cohort, time)

  Delta <- numeric(n_did)
  for (s in seq_len(n_did)) {
    d  <- did_list[[s]]
    fg <- d$focal_g; cg <- d$ctrl_g
    tp <- d$t_post;  tr <- d$t_pre

    yf_post <- gp_means[.(fg, tp), ybar]
    yf_pre  <- gp_means[.(fg, tr), ybar]
    yc_post <- gp_means[.(cg, tp), ybar]
    yc_pre  <- gp_means[.(cg, tr), ybar]

    if (any(is.na(c(yf_post, yf_pre, yc_post, yc_pre)))) {
      Delta[s] <- NA_real_
    } else {
      Delta[s] <- (yf_post - yf_pre) - (yc_post - yc_pre)
    }
  }
  Delta
}

# ===========================================================================
# 6. Simulation
# ===========================================================================

n_sims <- 1000

# Storage: for each CATT, store estimates from all three estimators
results <- list(
  Gardner_pkg = matrix(NA_real_, nrow = n_sims, ncol = n_catt),
  Gardner_W   = matrix(NA_real_, nrow = n_sims, ncol = n_catt),
  Ext_Gardner = matrix(NA_real_, nrow = n_sims, ncol = n_catt)
)

# Fixed unit and time FEs
alpha_i <- rnorm(N_total, mean = 0, sd = 2)
lambda_t <- rnorm(T_total, mean = 0, sd = 1)

for (sim in seq_len(n_sims)) {

  # --- 6a. Generate panel ---
  dt <- data.table(
    id     = rep(unit_ids, each = T_total),
    time   = rep(seq_len(T_total), times = N_total),
    cohort = rep(unit_cohort, each = T_total)
  )
  dt[, alpha := alpha_i[id]]
  dt[, lambda := lambda_t[time]]
  dt[, tau := mapply(true_catt, cohort, time)]
  dt[, eps := rnorm(.N, mean = 0, sd = sigma_eps)]
  dt[, y := alpha + lambda + tau + eps]

  # Treatment indicator (for did2s)
  dt[, treated  := (cohort > 0L) & (time >= cohort)]
  dt[, first_treat := ifelse(cohort == 0L, 0L, cohort)]

  # --- 6b. Gardner_pkg via did2s ---
  tryCatch({
    fit_g <- did2s(
      data          = dt,
      yname         = "y",
      first_stage   = ~ 0 | id + time,
      second_stage  = ~ i(cohort, time, ref = FALSE),
      treatment     = "treated",
      cluster_var   = "id"
    )
    coef_g <- coef(fit_g)

    for (k in seq_len(n_catt)) {
      g_val <- catt_list[[k]][1]; t_val <- catt_list[[k]][2]
      nm <- paste0("cohort::", g_val, ":time::", t_val)
      if (nm %in% names(coef_g)) {
        results$Gardner_pkg[sim, k] <- coef_g[nm]
      }
    }
  }, error = function(e) NULL)

  # --- 6c. Compute Delta vector ---
  Delta_all   <- compute_Delta(dt)
  Delta_clean <- Delta_all[clean_idx]

  # Handle any NA (should not occur in balanced panel)
  if (anyNA(Delta_all) || anyNA(Delta_clean)) next

  # --- 6d. Gardner_W: hat_beta = W * Delta_clean ---
  beta_gard_W <- as.vector(W_mat %*% Delta_clean)
  results$Gardner_W[sim, ] <- beta_gard_W

  # --- 6e. Ext_Gardner: hat_beta = (Q_H' Q_H)^{-1} Q_H' Delta_all ---
  beta_ext <- as.vector(QhQh_inv %*% t(Q_H) %*% Delta_all)
  results$Ext_Gardner[sim, ] <- beta_ext
}

# ===========================================================================
# 7. Results: Bias and Variance for each CATT
# ===========================================================================

cat("\n=== True CATTs ===\n")
true_vals <- sapply(seq_len(n_catt), function(k) {
  true_catt(catt_list[[k]][1], catt_list[[k]][2])
})
for (k in seq_len(n_catt)) {
  cat(sprintf("  ATT(%d,%d) = %.4f\n", catt_list[[k]][1], catt_list[[k]][2], true_vals[k]))
}

cat("\n=== Bias (should be ~0 for all estimators under parallel trends) ===\n")
cat(sprintf("%-20s %-20s %-20s %-20s\n", "CATT", "Gardner_pkg", "Gardner_W", "Ext_Gardner"))
for (k in seq_len(n_catt)) {
  b_pkg <- mean(results$Gardner_pkg[, k], na.rm = TRUE) - true_vals[k]
  b_W   <- mean(results$Gardner_W[, k],   na.rm = TRUE) - true_vals[k]
  b_ext <- mean(results$Ext_Gardner[, k], na.rm = TRUE) - true_vals[k]
  cat(sprintf("ATT(%d,%d)%-12s %+.6f             %+.6f             %+.6f\n",
              catt_list[[k]][1], catt_list[[k]][2], "", b_pkg, b_W, b_ext))
}

cat("\n=== Variance (lower = more efficient) ===\n")
cat(sprintf("%-20s %-20s %-20s %-20s\n", "CATT", "Gardner_pkg", "Gardner_W", "Ext_Gardner"))
for (k in seq_len(n_catt)) {
  v_pkg <- var(results$Gardner_pkg[, k], na.rm = TRUE)
  v_W   <- var(results$Gardner_W[, k],   na.rm = TRUE)
  v_ext <- var(results$Ext_Gardner[, k], na.rm = TRUE)
  ratio <- v_W / v_ext   # > 1 means Ext_Gardner more efficient
  cat(sprintf("ATT(%d,%d)%-12s %.6f             %.6f             %.6f  [ratio Gard/Ext=%.3f]\n",
              catt_list[[k]][1], catt_list[[k]][2], "", v_pkg, v_W, v_ext, ratio))
}

cat("\n=== Numerical equivalence: Gardner_pkg vs Gardner_W ===\n")
for (k in seq_len(n_catt)) {
  max_diff <- max(abs(results$Gardner_pkg[, k] - results$Gardner_W[, k]), na.rm = TRUE)
  cat(sprintf("  ATT(%d,%d): max|Gardner_pkg - Gardner_W| = %.2e\n",
              catt_list[[k]][1], catt_list[[k]][2], max_diff))
}
