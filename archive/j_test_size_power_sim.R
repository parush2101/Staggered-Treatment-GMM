###############################################################################
# Size and Power Simulation for the Iterated Satterthwaite-Box J-Test
#
# DGP: AR(1) errors, N=310, T=33, 5 treated cohorts + 1 never-treated,
# heterogeneous treatment effects (same as j_test_satterthwaite_ar1.R).
#
# Three violations of the parallel-trends assumption (set magnitude=0 for null):
#   1. trend         — cohort-specific linear pre/post trends (DIFFUSE)
#   2. shock         — one-time level shift in cohort 13 at t=11 (LOCALIZED)
#   3. anticipation  — fraction of treatment effect in 3 pre-periods (REALISTIC)
#
# For each (violation, magnitude) cell, n_sims independent datasets. For each
# dataset, draw B_max stratified subsets at each K in K_values and compute the
# Satterthwaite-corrected p-value per subset.
#
# Decision rules for B > 1 (paper's "iterated" procedure):
#   single   : reject if first p < α                   (B=1 baseline)
#   fraction : reject if mean(p_b < α)  > α            (paper's stated rule)
#   min_bonf : reject if min(p_b)       < α/B          (Bonferroni — strict)
#   median   : reject if median(p_b)    < α            (robust middle-ground)
#
# Outputs (printed at end):
#   Table 1  Size calibration (null DGP) by K, decision rule, α
#   Table 2  Power vs. magnitude — diffuse trend violation
#   Table 3  Power vs. magnitude — localized shock violation
#   Table 4  Power vs. magnitude — anticipation violation
#   Table 5  K-sensitivity at one fixed violation magnitude
#   Table 6  B-sensitivity at one fixed violation magnitude (key headline)
#
# DOES NOT execute on its own — source the file in R and run main().
###############################################################################

library(data.table)

set.seed(123)

# ===========================================================================
# 1. CONFIGURATION (tune these knobs)
# ===========================================================================

cfg <- list(
  # DGP
  cohort_size      = 50,
  n_never          = 60,
  T_total          = 33,
  treatment_times  = c(10, 13, 16, 19, 22),
  rho_true         = 0.5,
  beta_het         = c(-16, -12, -10, -9, -2),
  r_het            = c(0.01, 0.04, 0.08, 0.10, 0.07),

  # Test parameters
  K_values         = c(15L, 30L, 45L),
  B_values         = c(1L, 20L, 100L),
  alpha_levels     = c(0.01, 0.05, 0.10),

  # Simulation budget (drop n_sims to e.g. 100 for a quick smoke test)
  n_sims           = 300L,

  # Violation magnitude grids
  trend_mags       = c(0.02, 0.04, 0.06, 0.08, 0.10),
  shock_mags       = c(0.5, 1.0, 2.0, 4.0),
  anticip_mags     = c(0.10, 0.25, 0.50),

  # Sensitivity-table headline cells (must lie inside the grids above)
  k_sens_violation = list(type = "anticipation", mag = 0.25),
  b_sens_violation = list(type = "shock",        mag = 2.0)
)

# ===========================================================================
# 2. DERIVED CONSTANTS AND PANEL STRUCTURE (same as AR(1) script)
# ===========================================================================

N_total     <- 5 * cfg$cohort_size + cfg$n_never
n_cohorts   <- length(cfg$treatment_times)
unit_cohort <- c(rep(cfg$treatment_times, each = cfg$cohort_size),
                 rep(0, cfg$n_never))
get_cohort_size <- function(g_val) ifelse(g_val == 0, cfg$n_never, cfg$cohort_size)

all_groups  <- c(0, cfg$treatment_times)
n_groups    <- length(all_groups)
group_sizes <- sapply(all_groups, get_cohort_size)
group_masks <- lapply(all_groups, function(g) which(unit_cohort == g))

T_total <- cfg$T_total

# CATTs
treated_g <- sort(cfg$treatment_times)
catt_list <- list()
for (g_c in treated_g) {
  for (k in 0:(T_total - g_c)) {
    catt_list[[length(catt_list) + 1]] <- c(g_c, g_c + k)
  }
}
n_catt <- length(catt_list)
catt_g <- sapply(catt_list, `[`, 1)
catt_t <- sapply(catt_list, `[`, 2)

# DiD-moment metadata
did_meta <- list()
for (catt_idx in 1:n_catt) {
  g_c    <- catt_list[[catt_idx]][1]
  t_post <- catt_list[[catt_idx]][2]
  for (m in 1:(g_c - 1)) {
    t_pre <- g_c - m
    did_meta[[length(did_meta) + 1]] <- list(
      catt_idx = catt_idx, type = "never",
      focal_g = g_c, ctrl_g = 0, t_post = t_post, t_pre = t_pre)
    for (g_l in treated_g) {
      if (g_l <= t_post) next
      did_meta[[length(did_meta) + 1]] <- list(
        catt_idx = catt_idx, type = "notyet",
        focal_g = g_c, ctrl_g = g_l, t_post = t_post, t_pre = t_pre)
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
        focal_g = g_c, ctrl_g = g_j, t_post = t_post, t_pre = t_pre)
    }
  }
}
n_did <- length(did_meta)

# Incidence matrix Q_H
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
QtQ     <- crossprod(Q_H)
QtQ_inv <- solve(QtQ)

meta_type   <- sapply(did_meta, `[[`, "type")
idx_never   <- which(meta_type == "never")
idx_notyet  <- which(meta_type == "notyet")
idx_already <- which(meta_type == "already")

cat(sprintf("n_catt = %d, n_did = %d  (never: %d, notyet: %d, already: %d)\n",
            n_catt, n_did, length(idx_never), length(idx_notyet), length(idx_already)))

# Vectorized Delta indices
meta_focal <- sapply(did_meta, `[[`, "focal_g")
meta_ctrl  <- sapply(did_meta, `[[`, "ctrl_g")
meta_tp    <- sapply(did_meta, `[[`, "t_post")
meta_tr    <- sapply(did_meta, `[[`, "t_pre")
focal_gi   <- match(meta_focal, all_groups)
ctrl_gi    <- match(meta_ctrl,  all_groups)
idx_fp <- cbind(focal_gi, meta_tp)
idx_fr <- cbind(focal_gi, meta_tr)
idx_cp <- cbind(ctrl_gi,  meta_tp)
idx_cr <- cbind(ctrl_gi,  meta_tr)

# Z_g matrices and projector P_g for A=I weighting
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
P_group <- lapply(Z_group, function(Z) QtQ_inv %*% crossprod(Q_H, Z))

# Treatment-effect lookup
tau_lookup <- matrix(0, n_cohorts, T_total)
for (c_idx in 1:n_cohorts) {
  g_c <- cfg$treatment_times[c_idx]
  for (t in g_c:T_total) {
    tau_lookup[c_idx, t] <- cfg$beta_het[c_idx] * (1 + cfg$r_het[c_idx])^(t - g_c)
  }
}
unit_cohort_idx <- match(unit_cohort, cfg$treatment_times)

cat("Precomputation done.\n\n")

# ===========================================================================
# 3. DGP WITH OPTIONAL VIOLATION
# ===========================================================================

generate_data <- function(violation = "none", magnitude = 0) {
  alpha_i  <- rnorm(N_total)
  lambda_t <- rnorm(T_total)
  eps_mat  <- matrix(0, N_total, T_total)
  eps_mat[, 1] <- rnorm(N_total, 0, sd = 1 / sqrt(1 - cfg$rho_true^2))
  for (t in 2:T_total) {
    eps_mat[, t] <- cfg$rho_true * eps_mat[, t - 1] + rnorm(N_total)
  }
  Y_mat <- outer(alpha_i, rep(1, T_total)) +
           outer(rep(1, N_total), lambda_t) + eps_mat

  # Baseline treatment effects (heterogeneous, growing)
  for (c_idx in 1:n_cohorts) {
    units_c <- which(unit_cohort_idx == c_idx)
    for (t in cfg$treatment_times[c_idx]:T_total) {
      Y_mat[units_c, t] <- Y_mat[units_c, t] + tau_lookup[c_idx, t]
    }
  }

  if (magnitude == 0 || violation == "none") return(Y_mat)

  if (violation == "trend") {
    # Cohort-specific linear trends in Y (varies across g => parallel trends fail)
    # gamma_g for treated cohorts: spread around zero so no single cohort dominates
    gamma_treated <- magnitude * c(-2, -1, 0, 1, 2)
    for (c_idx in 1:n_cohorts) {
      units_c <- which(unit_cohort_idx == c_idx)
      Y_mat[units_c, ] <- Y_mat[units_c, ] +
        gamma_treated[c_idx] *
        matrix(1:T_total, nrow = length(units_c), ncol = T_total, byrow = TRUE)
    }
  } else if (violation == "shock") {
    # Localized: one cohort × one period (cohort treated at t=13, pre-period t=11)
    units_g13 <- which(unit_cohort == 13)
    Y_mat[units_g13, 11] <- Y_mat[units_g13, 11] + magnitude
  } else if (violation == "anticipation") {
    # Treated units respond in the 3 periods just before treatment date
    for (c_idx in 1:n_cohorts) {
      units_c <- which(unit_cohort_idx == c_idx)
      g_c     <- cfg$treatment_times[c_idx]
      for (lead in 1:3) {
        t_anticip <- g_c - lead
        if (t_anticip >= 1) {
          Y_mat[units_c, t_anticip] <- Y_mat[units_c, t_anticip] +
            magnitude * cfg$beta_het[c_idx]
        }
      }
    }
  } else {
    stop("Unknown violation: ", violation)
  }

  Y_mat
}

# ===========================================================================
# 4. SINGLE-DATASET SUMMARIES (beta, sigma2_eps, rho, R_ar1)
# ===========================================================================

compute_dataset_summaries <- function(Y_mat) {
  # Group-time means and Delta
  means_mat <- matrix(0, n_groups, T_total)
  for (gi in 1:n_groups) {
    means_mat[gi, ] <- colMeans(Y_mat[group_masks[[gi]], , drop = FALSE])
  }
  Delta <- means_mat[idx_fp] - means_mat[idx_fr] -
           means_mat[idx_cp] + means_mat[idx_cr]

  beta_hat <- as.numeric(QtQ_inv %*% crossprod(Q_H, Delta))

  # Two-way demeaned residuals -> rho_hat, sigma2_eps
  Y_adj_mat <- Y_mat
  for (c_idx in 1:n_cohorts) {
    units_c <- which(unit_cohort_idx == c_idx)
    for (ci in 1:n_catt) {
      if (catt_g[ci] == cfg$treatment_times[c_idx]) {
        Y_adj_mat[units_c, catt_t[ci]] <-
          Y_adj_mat[units_c, catt_t[ci]] - beta_hat[ci]
      }
    }
  }
  row_m <- rowMeans(Y_adj_mat)
  col_m <- colMeans(Y_adj_mat)
  grd_m <- mean(Y_adj_mat)
  resid_mat <- Y_adj_mat -
               outer(row_m, rep(1, T_total)) -
               outer(rep(1, N_total), col_m) + grd_m

  resid_lag  <- resid_mat[, 1:(T_total - 1)]
  resid_cur  <- resid_mat[, 2:T_total]
  rho_hat    <- sum(resid_cur * resid_lag) / sum(resid_lag^2)
  eta_hat    <- resid_cur - rho_hat * resid_lag
  sigma2_eta <- mean(eta_hat^2)
  sigma2_eps <- sigma2_eta / max(1 - rho_hat^2, 0.01)

  R_ar1 <- toeplitz(rho_hat^(0:(T_total - 1)))

  list(Delta = Delta, beta_hat = beta_hat,
       rho_hat = rho_hat, sigma2_eps = sigma2_eps, R_ar1 = R_ar1)
}

# ===========================================================================
# 5. SATTERTHWAITE p-VALUE FOR ONE SUBSET
# ===========================================================================

draw_stratified_subset <- function(K) {
  n_per_type <- K %/% 3
  c(sample(idx_never,   n_per_type),
    sample(idx_notyet,  n_per_type),
    sample(idx_already, n_per_type))
}

p_value_for_subset <- function(idx_K, ds) {
  Delta      <- ds$Delta
  beta_hat   <- ds$beta_hat
  R_ar1      <- ds$R_ar1
  sigma2_eps <- ds$sigma2_eps

  K   <- length(idx_K)
  Q_S <- Q_H[idx_K, , drop = FALSE]

  Omega_S <- matrix(0, K, K)
  V_S     <- matrix(0, K, K)
  for (gi in 1:n_groups) {
    w_g  <- sigma2_eps / group_sizes[gi]
    Z_gS <- Z_group[[gi]][idx_K, , drop = FALSE]
    R_gS <- Z_gS - Q_S %*% P_group[[gi]]
    Z_gSR <- Z_gS %*% R_ar1
    R_gSR <- R_gS %*% R_ar1
    Omega_S <- Omega_S + w_g * tcrossprod(Z_gSR, Z_gS)
    V_S     <- V_S     + w_g * tcrossprod(R_gSR, R_gS)
  }

  eig_O <- eigen(Omega_S, symmetric = TRUE)
  if (min(eig_O$values) < 1e-12) return(NA_real_)

  Omega_S_inv <- solve(Omega_S)
  e_K   <- (Delta - Q_H %*% beta_hat)[idx_K]
  T_raw <- as.numeric(t(e_K) %*% Omega_S_inv %*% e_K)

  sqrt_inv <- eig_O$vectors %*%
              diag(1 / sqrt(eig_O$values)) %*% t(eig_O$vectors)
  M   <- sqrt_inv %*% V_S %*% sqrt_inv
  lam <- pmax(eigen(M, symmetric = TRUE)$values, 0)
  s1  <- sum(lam); s2 <- sum(lam^2)
  if (s1 < 1e-10) return(1)

  c_d  <- s2 / s1
  nu_d <- s1^2 / s2
  1 - pchisq(T_raw / c_d, df = nu_d)
}

# ===========================================================================
# 6. RUN ONE (violation × magnitude) CELL
#
# Returns a 4-D array of p-values: [sim, K_index, B_index] where B_index runs
# over the SUBSET draws (1..B_max), so we can later derive any B in B_values.
# ===========================================================================

run_cell <- function(violation, magnitude, n_sims, K_values, B_max,
                     verbose_prefix = "") {
  pmat <- array(NA_real_,
                dim = c(n_sims, length(K_values), B_max),
                dimnames = list(NULL, paste0("K", K_values), NULL))

  for (sim in 1:n_sims) {
    Y_mat <- generate_data(violation, magnitude)
    ds    <- compute_dataset_summaries(Y_mat)
    for (ki in seq_along(K_values)) {
      K <- K_values[ki]
      for (b in 1:B_max) {
        idx_K <- draw_stratified_subset(K)
        pmat[sim, ki, b] <- p_value_for_subset(idx_K, ds)
      }
    }
    if (sim %% 50 == 0) {
      cat(sprintf("    %s sim %d/%d\n", verbose_prefix, sim, n_sims))
    }
  }
  pmat
}

# ===========================================================================
# 7. DECISION RULES — given a vector of p-values from B subsets
# ===========================================================================

apply_rule <- function(p_vec, B, alpha, rule) {
  pv <- p_vec[1:B]
  pv <- pv[!is.na(pv)]
  if (length(pv) == 0) return(NA)
  switch(rule,
         "single"   = as.numeric(pv[1] < alpha),
         "fraction" = as.numeric(mean(pv < alpha) > alpha),
         "min_bonf" = as.numeric(min(pv) < alpha / B),
         "median"   = as.numeric(median(pv) < alpha))
}

rejection_rate <- function(pmat_K, B, alpha, rule) {
  rejects <- apply(pmat_K, 1, function(p_vec) apply_rule(p_vec, B, alpha, rule))
  mean(rejects, na.rm = TRUE)
}

# ===========================================================================
# 8. MAIN: run all cells, build summary tables
# ===========================================================================

main <- function() {
  B_max <- max(cfg$B_values)
  scenarios <- list(
    null         = list(type = "trend",        mags = 0),
    trend        = list(type = "trend",        mags = cfg$trend_mags),
    shock        = list(type = "shock",        mags = cfg$shock_mags),
    anticipation = list(type = "anticipation", mags = cfg$anticip_mags)
  )

  results <- list()  # results[[scenario_name]][[mag_str]] = pmat
  t_start <- Sys.time()
  for (scn in names(scenarios)) {
    s <- scenarios[[scn]]
    results[[scn]] <- list()
    for (mag in s$mags) {
      cat(sprintf("[%s] violation=%s, magnitude=%g\n",
                  format(Sys.time(), "%H:%M:%S"), s$type, mag))
      results[[scn]][[as.character(mag)]] <-
        run_cell(s$type, mag, cfg$n_sims, cfg$K_values, B_max,
                 verbose_prefix = sprintf("[%s mag=%g]", s$type, mag))
    }
  }
  cat(sprintf("\nAll cells done in %.1f minutes.\n",
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

  # =========================================================================
  # Helpers for table construction
  # =========================================================================
  rules_to_show <- c("single", "fraction", "min_bonf", "median")

  # Returns rejection rate for given (pmat, B, alpha, rule), with B clamped
  rate <- function(pmat_full, ki, B, alpha, rule) {
    rejection_rate(pmat_full[, ki, , drop = FALSE][, 1, ], B, alpha, rule)
  }

  # =========================================================================
  # TABLE 1 — Size calibration (null DGP)
  # =========================================================================
  pmat_null <- results$null[["0"]]
  cat("\n\n")
  cat("==========================================================================\n")
  cat("TABLE 1.  SIZE CALIBRATION (parallel trends hold; magnitude = 0)\n")
  cat(sprintf("           DGP: AR(1) rho=%.1f, N=%d, T=%d, %d sims\n",
              cfg$rho_true, N_total, T_total, cfg$n_sims))
  cat("==========================================================================\n")
  cat(sprintf("%-6s %-10s %-8s %8s %8s %8s\n",
              "K", "Rule", "B",
              sprintf("a=%.2f", cfg$alpha_levels[1]),
              sprintf("a=%.2f", cfg$alpha_levels[2]),
              sprintf("a=%.2f", cfg$alpha_levels[3])))
  cat(paste(rep("-", 60), collapse = ""), "\n")
  for (ki in seq_along(cfg$K_values)) {
    K <- cfg$K_values[ki]
    for (rule in rules_to_show) {
      Bset <- if (rule == "single") 1L else cfg$B_values
      for (B in Bset) {
        if (rule == "single" && B != 1L) next
        row <- sapply(cfg$alpha_levels, function(a) rate(pmat_null, ki, B, a, rule))
        cat(sprintf("%-6d %-10s %-8d %7.1f%% %7.1f%% %7.1f%%\n",
                    K, rule, B, 100*row[1], 100*row[2], 100*row[3]))
      }
    }
    cat("\n")
  }
  cat("Reference: rejection rate should be close to nominal alpha.\n")

  # =========================================================================
  # TABLES 2-4 — Power vs. magnitude, one per violation type
  # =========================================================================
  power_table <- function(scn_name, table_no, title) {
    pmats <- results[[scn_name]]
    mags  <- as.numeric(names(pmats))
    cat("\n\n")
    cat("==========================================================================\n")
    cat(sprintf("TABLE %d.  POWER vs. MAGNITUDE — %s\n", table_no, title))
    cat(sprintf("           Reject at alpha=0.05; %d sims per cell\n", cfg$n_sims))
    cat("==========================================================================\n")
    header_mag <- paste(sprintf("%8s", sprintf("m=%.2f", mags)), collapse = " ")
    cat(sprintf("%-6s %-10s %-6s %s\n", "K", "Rule", "B", header_mag))
    cat(paste(rep("-", 30 + 9 * length(mags)), collapse = ""), "\n")
    for (ki in seq_along(cfg$K_values)) {
      K <- cfg$K_values[ki]
      for (rule in rules_to_show) {
        Bset <- if (rule == "single") 1L else cfg$B_values
        for (B in Bset) {
          if (rule == "single" && B != 1L) next
          row <- sapply(mags, function(m)
            rate(pmats[[as.character(m)]], ki, B, 0.05, rule))
          row_str <- paste(sprintf("%7.1f%%", 100 * row), collapse = " ")
          cat(sprintf("%-6d %-10s %-6d %s\n", K, rule, B, row_str))
        }
      }
      cat("\n")
    }
  }
  power_table("trend",        2, "diffuse trend (cohort-specific linear trends)")
  power_table("shock",        3, "localized shock (cohort 13, t=11)")
  power_table("anticipation", 4, "anticipation (3 pre-periods, fraction kappa)")

  # =========================================================================
  # TABLE 5 — K-sensitivity at one fixed (violation, magnitude)
  # =========================================================================
  ksv <- cfg$k_sens_violation
  pmats_k <- results[[ksv$type]]
  pmat_k  <- pmats_k[[as.character(ksv$mag)]]
  cat("\n\n")
  cat("==========================================================================\n")
  cat(sprintf("TABLE 5.  K-SENSITIVITY  (violation=%s, magnitude=%g, alpha=0.05)\n",
              ksv$type, ksv$mag))
  cat("==========================================================================\n")
  cat(sprintf("%-10s %-6s %s\n", "Rule", "B",
              paste(sprintf("%8s", paste0("K=", cfg$K_values)), collapse = " ")))
  cat(paste(rep("-", 22 + 9 * length(cfg$K_values)), collapse = ""), "\n")
  for (rule in rules_to_show) {
    Bset <- if (rule == "single") 1L else cfg$B_values
    for (B in Bset) {
      if (rule == "single" && B != 1L) next
      row <- sapply(seq_along(cfg$K_values),
                    function(ki) rate(pmat_k, ki, B, 0.05, rule))
      row_str <- paste(sprintf("%7.1f%%", 100 * row), collapse = " ")
      cat(sprintf("%-10s %-6d %s\n", rule, B, row_str))
    }
  }
  cat("\nLarger K should give more power, but raises chance of singular Omega_S.\n")

  # =========================================================================
  # TABLE 6 — B-sensitivity at one fixed (violation, magnitude)  [HEADLINE]
  # =========================================================================
  bsv <- cfg$b_sens_violation
  pmats_b <- results[[bsv$type]]
  pmat_b  <- pmats_b[[as.character(bsv$mag)]]
  cat("\n\n")
  cat("==========================================================================\n")
  cat(sprintf("TABLE 6.  B-SENSITIVITY  (violation=%s, magnitude=%g, alpha=0.05)\n",
              bsv$type, bsv$mag))
  cat("           [Headline: does iterating subsets recover power for localized\n")
  cat("            violations that single subset draws miss?]\n")
  cat("==========================================================================\n")
  cat(sprintf("%-6s %-10s %s\n", "K", "Rule",
              paste(sprintf("%8s", paste0("B=", cfg$B_values)), collapse = " ")))
  cat(paste(rep("-", 22 + 9 * length(cfg$B_values)), collapse = ""), "\n")
  for (ki in seq_along(cfg$K_values)) {
    K <- cfg$K_values[ki]
    for (rule in rules_to_show) {
      if (rule == "single") {
        # Only report B=1 column for "single"; show "—" for larger B
        row_strs <- character(length(cfg$B_values))
        for (bi in seq_along(cfg$B_values)) {
          B <- cfg$B_values[bi]
          row_strs[bi] <- if (B == 1L) {
            sprintf("%7.1f%%", 100 * rate(pmat_b, ki, 1L, 0.05, rule))
          } else "       —"
        }
      } else {
        row_strs <- sapply(cfg$B_values, function(B) {
          if (B == 1L && rule != "single") return("       —")  # single==first p
          sprintf("%7.1f%%", 100 * rate(pmat_b, ki, B, 0.05, rule))
        })
      }
      cat(sprintf("%-6d %-10s %s\n", K, rule, paste(row_strs, collapse = " ")))
    }
    cat("\n")
  }
  cat("Comparing 'single' (B=1) row to 'fraction'/'median' rows at B=100\n")
  cat("isolates the value of iterating the subset draw.\n")

  # =========================================================================
  # SAVE
  # =========================================================================
  dir.create(path.expand("~/econ_sims"), showWarnings = FALSE, recursive = TRUE)
  save(results, cfg,
       file = path.expand("~/econ_sims/j_test_size_power.RData"))
  cat("\nResults saved to ~/econ_sims/j_test_size_power.RData\n")

  invisible(results)
}

# ===========================================================================
# 9. ENTRY POINT
# ===========================================================================
# To run:    source("j_test_size_power_sim.R"); main()
#
# Quick smoke test (~1-2 min):  set cfg$n_sims <- 30; cfg$B_values <- c(1L, 20L)
# Full run     (~15-30 min):    defaults
#
# Don't auto-execute on source.
if (FALSE) {
  main()
}
