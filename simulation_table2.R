###############################################################################
# Replication of Table 2: Bias and Variances of Estimators
# Paper: "Estimating Treatment Effects under Staggered Timing using GMM"
# Authors: Parush Arora and Rishabh Bijani
#
# This script replicates Table 2 from the paper, comparing the proposed GMM
# estimator against TWFE, Callaway-Sant'Anna (CS), Sun-Abraham (SA),
# De Chaisemartin-d'Haultfoeuille (DCdH), and Gardner (did2s) estimators
# under both homogeneous and heterogeneous treatment effects.
###############################################################################

# --- Load required packages ---
library(data.table)
library(fixest)       # TWFE + Sun-Abraham
library(did)          # Callaway-Sant'Anna
library(did2s)        # Gardner (two-stage DiD)
library(MASS)
library(Matrix)

# DIDmultiplegt may fail to install on some systems; handle gracefully
has_dcdh <- requireNamespace("DIDmultiplegt", quietly = TRUE)
if (has_dcdh) {
  library(DIDmultiplegt)
} else {
  cat("NOTE: DIDmultiplegt package not available. DCdH estimator will be\n")
  cat("      implemented manually using the switchers-based approach.\n\n")
}

set.seed(42)

# ===========================================================================
# 1. Simulation Parameters
# ===========================================================================

N_total     <- 51          # Total units (states)
T_total     <- 33          # Total time periods
n_sims      <- 500         # Number of simulations
cohort_size <- 8           # States per treated cohort
n_never     <- 11          # Never-treated states

# Treatment timing for 5 treated cohorts
treatment_times <- c(10, 13, 16, 19, 22)
n_cohorts       <- length(treatment_times)

# Assign units to cohorts
# Cohorts: g=10 (units 1-8), g=13 (9-16), g=16 (17-24), g=19 (25-32), g=22 (33-40)
# Never-treated: units 41-51 (g=0 in our notation, g=Inf for packages)
unit_cohort <- c(rep(treatment_times, each = cohort_size),
                 rep(0, n_never))

# ===========================================================================
# 2. Data Generating Process
# ===========================================================================

generate_data <- function(beta_g_vec, r_g_vec) {
  # beta_g_vec: baseline treatment effects for each of 5 cohorts
  # r_g_vec:    growth rates for each of 5 cohorts

  dt_list <- vector("list", N_total)

  for (i in 1:N_total) {
    g_i     <- unit_cohort[i]
    alpha_i <- rnorm(1)           # unit fixed effect

    for (t in 1:T_total) {
      lambda_t <- rnorm(1)        # time fixed effect (drawn fresh per unit-time for idiosyncratic component)
      eps_it   <- rnorm(1)        # error

      # Treatment indicator
      if (g_i == 0) {
        D_it  <- 0
        tau_it <- 0
      } else {
        D_it <- as.integer(t >= g_i)
        # Find which cohort index this unit belongs to
        cohort_idx <- which(treatment_times == g_i)
        beta_g <- beta_g_vec[cohort_idx]
        r_g    <- r_g_vec[cohort_idx]
        tau_it <- beta_g * (1 + r_g)^(t - g_i) * D_it
      }

      dt_list[[length(unlist(dt_list))/5 + 1]] <- NULL
    }
  }

  # More efficient: build directly
  unit_id <- rep(1:N_total, each = T_total)
  time_id <- rep(1:T_total, times = N_total)

  alpha <- rnorm(N_total)
  lambda <- rnorm(T_total)
  eps <- rnorm(N_total * T_total)

  g_vec  <- unit_cohort[unit_id]
  D_vec  <- as.integer(g_vec > 0 & time_id >= g_vec)

  # Treatment effect
  tau_vec <- numeric(N_total * T_total)
  for (c_idx in 1:n_cohorts) {
    g_c <- treatment_times[c_idx]
    mask <- (g_vec == g_c) & (D_vec == 1)
    tau_vec[mask] <- beta_g_vec[c_idx] * (1 + r_g_vec[c_idx])^(time_id[mask] - g_c)
  }

  Y_vec <- alpha[unit_id] + lambda[time_id] + tau_vec + eps

  dt <- data.table(
    unit = unit_id,
    time = time_id,
    Y    = Y_vec,
    D    = D_vec,
    g    = g_vec  # cohort (first treatment period; 0 = never treated)
  )

  # For packages that need g=0 as Inf or NA
  dt[, g_inf := fifelse(g == 0, Inf, as.numeric(g))]
  dt[, g_cs  := fifelse(g == 0, 0L, as.integer(g))]  # did package: 0 = never-treated

  return(dt)
}

# ===========================================================================
# 3. Compute True ATT (unit-weighted)
# ===========================================================================

compute_true_att <- function(beta_g_vec, r_g_vec) {
  # Unit-weighted ATT: average of all unit-time-level treatment effects
  total_te  <- 0
  total_obs <- 0

  for (c_idx in 1:n_cohorts) {
    g_c    <- treatment_times[c_idx]
    beta_g <- beta_g_vec[c_idx]
    r_g    <- r_g_vec[c_idx]
    n_g    <- cohort_size

    for (t in g_c:T_total) {
      te <- beta_g * (1 + r_g)^(t - g_c)
      total_te  <- total_te + n_g * te
      total_obs <- total_obs + n_g
    }
  }
  return(total_te / total_obs)
}

# ===========================================================================
# 4. GMM Estimator (from the paper)
# ===========================================================================

gmm_estimator <- function(dt) {
  # Enumerate all 2x2 DiD estimates and build Q_H, Delta

  # Compute cohort-time means
  cohort_means <- dt[, .(Y_mean = mean(Y)), by = .(g, time)]
  setkey(cohort_means, g, time)

  get_mean <- function(g_val, t_val) {
    res <- cohort_means[.(g_val, t_val), Y_mean]
    if (length(res) == 0 || is.na(res)) return(NA_real_)
    return(res)
  }

  # Treated cohorts
  treated_g <- sort(treatment_times)
  # All CATTs: beta_{g, g+k} for each treated cohort g, k = 0, ..., T-g
  catt_list <- list()
  for (g_c in treated_g) {
    for (k in 0:(T_total - g_c)) {
      catt_list[[length(catt_list) + 1]] <- c(g_c, g_c + k)
    }
  }
  n_catt <- length(catt_list)
  catt_labels <- sapply(catt_list, function(x) paste0("b_", x[1], "_", x[2]))

  # Enumerate all 2x2 DiD estimates
  did_estimates <- list()  # each: list(delta, g, k, type, control_g, m, catt_idx, bias_catts)

  for (catt_idx in 1:n_catt) {
    g_c <- catt_list[[catt_idx]][1]
    t_post <- catt_list[[catt_idx]][2]
    k <- t_post - g_c

    # Pre-treatment periods: m = 1, ..., g_c - 1
    for (m in 1:(g_c - 1)) {
      t_pre <- g_c - m

      # (a) Never-treated control (c = 0)
      Y_g_post <- get_mean(g_c, t_post)
      Y_g_pre  <- get_mean(g_c, t_pre)
      Y_0_post <- get_mean(0, t_post)
      Y_0_pre  <- get_mean(0, t_pre)

      if (!any(is.na(c(Y_g_post, Y_g_pre, Y_0_post, Y_0_pre)))) {
        delta <- (Y_g_post - Y_g_pre) - (Y_0_post - Y_0_pre)
        did_estimates[[length(did_estimates) + 1]] <- list(
          delta = delta, catt_idx = catt_idx, type = "never",
          bias_catts = NULL  # unbiased
        )
      }

      # (b) Not-yet-treated control
      for (g_l in treated_g) {
        l <- g_l - g_c
        if (l <= k) next  # need l > k
        if (g_l <= t_post) next  # control must not be treated at t_post

        Y_l_post <- get_mean(g_l, t_post)
        Y_l_pre  <- get_mean(g_l, t_pre)

        if (!any(is.na(c(Y_g_post, Y_g_pre, Y_l_post, Y_l_pre)))) {
          delta <- (Y_g_post - Y_g_pre) - (Y_l_post - Y_l_pre)
          did_estimates[[length(did_estimates) + 1]] <- list(
            delta = delta, catt_idx = catt_idx, type = "notyet",
            bias_catts = NULL  # unbiased under heterogeneous effects
          )
        }
      }

      # (c) Already-treated control (forbidden comparison)
      for (g_j in treated_g) {
        j <- g_c - g_j
        if (j <= m) next    # need j > m (i.e., g_j < g_c - m = t_pre)
        if (g_j >= g_c) next  # must be already treated

        Y_j_post <- get_mean(g_j, t_post)
        Y_j_pre  <- get_mean(g_j, t_pre)

        if (!any(is.na(c(Y_g_post, Y_g_pre, Y_j_post, Y_j_pre)))) {
          delta <- (Y_g_post - Y_g_pre) - (Y_j_post - Y_j_pre)

          # Under heterogeneous effects, E[delta] = beta_{g,g+k} - beta_{g-j,g+k} + beta_{g-j,g-m}
          # So bias catts: subtract beta_{g-j,g+k} and add beta_{g-j,g-m}
          # g-j = g_j, post period = t_post = g+k, pre period = t_pre = g-m
          bias_neg_idx <- NULL  # index of beta_{g_j, t_post}
          bias_pos_idx <- NULL  # index of beta_{g_j, t_pre}

          for (ci in 1:n_catt) {
            if (catt_list[[ci]][1] == g_j && catt_list[[ci]][2] == t_post) bias_neg_idx <- ci
            if (catt_list[[ci]][1] == g_j && catt_list[[ci]][2] == t_pre)  bias_pos_idx <- ci
          }

          did_estimates[[length(did_estimates) + 1]] <- list(
            delta = delta, catt_idx = catt_idx, type = "already",
            bias_neg = bias_neg_idx,  # subtract this CATT
            bias_pos = bias_pos_idx   # add this CATT
          )
        }
      }
    }
  }

  n_did <- length(did_estimates)
  if (n_did == 0) return(list(att = NA, catts = rep(NA, n_catt)))

  # Build Delta vector and Q_H matrix
  Delta <- numeric(n_did)
  Q_H   <- matrix(0, nrow = n_did, ncol = n_catt)

  for (s in 1:n_did) {
    est <- did_estimates[[s]]
    Delta[s] <- est$delta
    Q_H[s, est$catt_idx] <- 1  # maps to its own CATT

    if (est$type == "already") {
      # Bias correction
      if (!is.null(est$bias_neg) && !is.na(est$bias_neg)) {
        Q_H[s, est$bias_neg] <- Q_H[s, est$bias_neg] - 1
      }
      if (!is.null(est$bias_pos) && !is.na(est$bias_pos)) {
        Q_H[s, est$bias_pos] <- Q_H[s, est$bias_pos] + 1
      }
    }
  }

  # GMM estimation: beta_hat = (Q_H' A Q_H)^{-1} Q_H' A Delta
  # Start with A = I (OLS), then iterate

  # Step 1: OLS estimate
  QtQ <- crossprod(Q_H)
  QtD <- crossprod(Q_H, Delta)

  # Check if system is solvable
  if (rcond(QtQ) < 1e-15) {
    # Use pseudoinverse
    beta_hat <- tryCatch(solve(QtQ, QtD), error = function(e) rep(NA, n_catt))
  } else {
    beta_hat <- solve(QtQ, QtD)
  }

  if (any(is.na(beta_hat))) {
    return(list(att = NA, catts = beta_hat))
  }

  # Step 2: Iterative GMM - estimate variance and re-weight
  for (iter in 1:5) {
    residuals <- Delta - Q_H %*% beta_hat
    # Estimate Omega_hat (variance of moment conditions)
    Omega_hat <- crossprod(residuals) / n_did * diag(n_did)

    # Use residual-based diagonal weighting
    resid_sq <- as.numeric(residuals)^2
    resid_sq[resid_sq < 1e-10] <- 1e-10  # floor
    A <- diag(1 / resid_sq)

    QtAQ <- crossprod(Q_H, A %*% Q_H)
    QtAD <- crossprod(Q_H, A %*% Delta)

    if (rcond(QtAQ) < 1e-15) break
    beta_hat <- solve(QtAQ, QtAD)
  }

  # Compute ATT using unit-weighted aggregation
  # Weight each CATT by cohort_size (all same here = 8)
  # ATT = sum over all treated (g,t) of beta_{g,t} / count
  weights <- numeric(n_catt)
  for (ci in 1:n_catt) {
    g_c <- catt_list[[ci]][1]
    weights[ci] <- cohort_size
  }
  weights <- weights / sum(weights)

  att <- sum(weights * beta_hat)

  return(list(att = as.numeric(att), catts = as.numeric(beta_hat)))
}

# ===========================================================================
# 5. Wrapper functions for each estimator
# ===========================================================================

estimate_twfe <- function(dt) {
  tryCatch({
    mod <- feols(Y ~ D | unit + time, data = dt)
    return(coef(mod)["D"])
  }, error = function(e) NA_real_)
}

estimate_cs <- function(dt) {
  tryCatch({
    out <- att_gt(
      yname  = "Y",
      tname  = "time",
      idname = "unit",
      gname  = "g_cs",
      data   = as.data.frame(dt),
      control_group = "nevertreated",
      print_details = FALSE,
      bstrap = FALSE,
      cband  = FALSE
    )
    agg <- aggte(out, type = "simple")
    return(agg$overall.att)
  }, error = function(e) NA_real_)
}

estimate_sa <- function(dt) {
  tryCatch({
    mod <- feols(Y ~ sunab(g_inf, time) | unit + time, data = dt)
    agg <- summary(mod, agg = "ATT")
    return(agg$coeftable[1, 1])
  }, error = function(e) NA_real_)
}

estimate_dcdh <- function(dt) {
  # De Chaisemartin and d'Haultfoeuille (2020) estimator
  # Uses switchers: units whose treatment status changes between t-1 and t
  # ATT_DM = weighted average of DiD for switchers vs non-switchers
  if (has_dcdh) {
    tryCatch({
      df <- as.data.frame(dt[, .(unit, time, Y, D)])
      out <- did_multiplegt(
        df    = df,
        Y     = "Y",
        G     = "unit",
        T     = "time",
        D     = "D",
        brep  = 0,
        parallel = FALSE
      )
      if (is.list(out)) {
        att_val <- out$effect
        if (is.null(att_val)) att_val <- out[1]
      } else {
        att_val <- out[1]
      }
      return(as.numeric(att_val))
    }, error = function(e) NA_real_)
  } else {
    # Manual implementation of the DCdH instantaneous treatment effect
    # For each period t, compare switchers (D changes from 0 to 1) against
    # non-switchers (D stays at 0) using the change in Y
    tryCatch({
      dt_wide <- dcast(dt, unit ~ time, value.var = c("Y", "D"))
      effects <- c()
      weights <- c()

      for (t in 2:T_total) {
        Y_t   <- dt[time == t, .(unit, Y_t = Y, D_t = D)]
        Y_tm1 <- dt[time == t - 1, .(unit, Y_tm1 = Y, D_tm1 = D)]
        merged <- merge(Y_t, Y_tm1, by = "unit")

        # Switchers: D goes from 0 to 1
        switchers     <- merged[D_tm1 == 0 & D_t == 1]
        # Non-switchers who stay at D=0
        non_switchers <- merged[D_tm1 == 0 & D_t == 0]

        if (nrow(switchers) > 0 & nrow(non_switchers) > 0) {
          dy_switch    <- mean(switchers$Y_t - switchers$Y_tm1)
          dy_nonswitch <- mean(non_switchers$Y_t - non_switchers$Y_tm1)
          effects <- c(effects, dy_switch - dy_nonswitch)
          weights <- c(weights, nrow(switchers))
        }
      }

      if (length(effects) == 0) return(NA_real_)
      return(sum(effects * weights) / sum(weights))
    }, error = function(e) NA_real_)
  }
}

estimate_gardner <- function(dt) {
  tryCatch({
    # did2s expects a first_treat variable
    dt_g <- copy(dt)
    dt_g[, first_treat := fifelse(g == 0, Inf, as.numeric(g))]

    mod <- did2s(
      data          = as.data.frame(dt_g),
      yname         = "Y",
      first_stage   = ~ 0 | unit + time,
      second_stage  = ~ i(D, ref = 0),
      treatment     = "D",
      cluster_var   = "unit",
      verbose       = FALSE
    )
    return(coef(mod)["D::1"])
  }, error = function(e) NA_real_)
}

# ===========================================================================
# 6. Run Simulations
# ===========================================================================

run_simulation <- function(beta_g_vec, r_g_vec, label) {
  true_att <- compute_true_att(beta_g_vec, r_g_vec)
  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("True ATT: %.4f\n", true_att))

  results <- data.table(
    sim     = integer(),
    TWFE    = numeric(),
    CS      = numeric(),
    SA      = numeric(),
    DCdH    = numeric(),
    Gardner = numeric(),
    GMM     = numeric()
  )

  for (s in 1:n_sims) {
    if (s %% 50 == 0) cat(sprintf("  Simulation %d/%d\n", s, n_sims))

    dt <- generate_data(beta_g_vec, r_g_vec)

    att_twfe    <- estimate_twfe(dt)
    att_cs      <- estimate_cs(dt)
    att_sa      <- estimate_sa(dt)
    att_dcdh    <- estimate_dcdh(dt)
    att_gardner <- estimate_gardner(dt)
    att_gmm     <- gmm_estimator(dt)$att

    results <- rbindlist(list(results, data.table(
      sim = s, TWFE = att_twfe, CS = att_cs, SA = att_sa,
      DCdH = att_dcdh, Gardner = att_gardner, GMM = att_gmm
    )))
  }

  # Compute bias and variance
  summary_dt <- data.table(
    Estimator = c("TWFE", "CS", "SA", "DCdH", "Gardner", "GMM"),
    Bias      = numeric(6),
    Variance  = numeric(6)
  )

  for (i in 1:6) {
    est_name <- summary_dt$Estimator[i]
    vals     <- results[[est_name]]
    vals     <- vals[!is.na(vals)]
    summary_dt$Bias[i]     <- round(mean(vals) - true_att, 4)
    summary_dt$Variance[i] <- round(var(vals), 4)
  }

  cat(sprintf("\nResults (%s):\n", label))
  print(summary_dt)

  return(list(results = results, summary = summary_dt, true_att = true_att))
}

# --- Homogeneous Effects ---
cat("========================================\n")
cat("  HOMOGENEOUS EFFECTS SIMULATION\n")
cat("========================================\n")

beta_hom <- rep(-5, n_cohorts)
r_hom    <- rep(0, n_cohorts)

res_hom <- run_simulation(beta_hom, r_hom, "Homogeneous")

# --- Heterogeneous Effects ---
cat("\n\n========================================\n")
cat("  HETEROGENEOUS EFFECTS SIMULATION\n")
cat("========================================\n")

beta_het <- c(-16, -12, -10, -9, -2)
r_het    <- c(0.01, 0.04, 0.08, 0.10, 0.07)

res_het <- run_simulation(beta_het, r_het, "Heterogeneous")

# ===========================================================================
# 7. Combine and Display Table 2
# ===========================================================================

cat("\n\n")
cat("============================================================\n")
cat("  TABLE 2: Bias and Variances of Estimators (Replication)\n")
cat("============================================================\n\n")

table2 <- merge(
  res_hom$summary[, .(Estimator, Hom_Bias = Bias, Hom_Var = Variance)],
  res_het$summary[, .(Estimator, Het_Bias = Bias, Het_Var = Variance)],
  by = "Estimator"
)

# Reorder to match paper
table2 <- table2[match(c("TWFE", "CS", "SA", "DCdH", "Gardner", "GMM"),
                       table2$Estimator)]

cat(sprintf("%-10s  %12s  %12s  %12s  %12s\n",
            "Estimator", "Hom. Bias", "Hom. Var", "Het. Bias", "Het. Var"))
cat(paste(rep("-", 62), collapse = ""), "\n")
for (i in 1:nrow(table2)) {
  cat(sprintf("%-10s  %12.4f  %12.4f  %12.4f  %12.4f\n",
              table2$Estimator[i],
              table2$Hom_Bias[i], table2$Hom_Var[i],
              table2$Het_Bias[i], table2$Het_Var[i]))
}

cat("\n\nPaper's Table 2 (for reference):\n")
cat(sprintf("%-10s  %12s  %12s  %12s  %12s\n",
            "Estimator", "Hom. Bias", "Hom. Var", "Het. Bias", "Het. Var"))
cat(paste(rep("-", 62), collapse = ""), "\n")
ref <- data.frame(
  est = c("TWFE", "CS", "SA", "DCdH", "Gardner", "GMM"),
  hb  = c(0.0020, 0.0170, 0.0080, -0.0110, -0.0070, 0.0072),
  hv  = c(0.009, 0.045, 0.049, 0.045, 0.0104, 0.0078),
  htb = c(3.5600, -0.018, -0.0060, 0.0025, -0.0070, 0.0100),
  htv = c(0.009, 0.540, 0.406, 0.431, 0.0104, 0.0098)
)
for (i in 1:6) {
  cat(sprintf("%-10s  %12.4f  %12.4f  %12.4f  %12.4f\n",
              ref$est[i], ref$hb[i], ref$hv[i], ref$htb[i], ref$htv[i]))
}

# Save results
save(res_hom, res_het, table2, file = "simulation_results.RData")
cat("\nResults saved to simulation_results.RData\n")
