###############################################################################
# gardner_forbidden_efficiency_ar1.R
#
# Extends gardner_forbidden_efficiency.R to an AR(1) error DGP.
#
# Primary question: Does AR1_Clean = AR1_All under serial correlation?
#   Under iid (spherical) errors, the Redundancy Proposition proves
#   Sph_Clean = Sph_All: forbidden DiDs add no information.
#   The proof exploits cross-cohort independence (Assumption 1), NOT iid errors.
#   If Assumption 1 is the key condition, redundancy should persist under AR(1).
#
# DGP: Y_it = alpha_i + lambda_t + tau_it + eps_it
#   eps_it = rho * eps_{i,t-1} + sqrt(1-rho^2) * eta_it,  eta_it ~ iid N(0,sigma^2)
#   => Cov(eps_it, eps_it') = sigma^2 * rho^|t-t'|  (stationary AR(1))
#
# Analytical Omega_phi under AR(1) (paper's formula, sigma_d = sigma^2 * rho^d):
#   Cov(Delta_s, Delta_s') = sigma^2 * sum_g sum_{t,t'} z_{s,g,t} * z_{s',g,t'} * rho^|t-t'|
#   where z_{s,g,t} = signed weight of cell (g,t) in DiD s (sign/N_g).
#   Under iid (rho=0): rho^|t-t'| = 1{t=t'}, recovering the spherical formula.
#
# Five estimators:
#   1. Gardner     : Two-stage feols (OLS on residuals — optimal under iid,
#                    generally suboptimal under AR(1))
#   2. EqWt_Clean  : A=I on clean DiDs only (unweighted)
#   3. AR1_Clean   : Optimal GMM on clean DiDs, A = pinv(Omega_cl^AR1)
#   4. EqWt_All    : A=I on all DiDs (clean + forbidden)
#   5. AR1_All     : Optimal GMM on all DiDs, A = pinv(Omega_phi^AR1)
#
# Key diagnostic: per-sim max|AR1_Clean - AR1_All| to test redundancy under AR(1).
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(MASS)   # for ginv
})

set.seed(2025)

# ===========================================================================
# 1. DGP Parameters
# ===========================================================================

rho_ar1     <- 0.5          # AR(1) serial correlation parameter
T_total     <- 12L
cohort_size <- 10L
n_never     <- 10L
treat_times <- c(4L, 7L, 10L)
N_total     <- length(treat_times) * cohort_size + n_never   # 40
sigma_eps   <- 1.0

beta_g_val  <- c("4" = 3.0, "7" = -2.0, "10" = 5.0)
r_g_val     <- c("4" = 0.10, "7" = 0.05, "10" = 0.00)

unit_cohort <- c(rep(treat_times, each = cohort_size), rep(0L, n_never))
unit_ids    <- seq_len(N_total)

true_catt_val <- function(g, t) {
  if (g == 0L || t < g) return(0.0)
  beta_g_val[as.character(g)] * (1 + r_g_val[as.character(g)])^(t - g)
}

# Fixed unit and time FEs (held constant across sims — only eps varies)
set.seed(2025)
alpha_i  <- rnorm(N_total, sd = 2.0)
lambda_t <- rnorm(T_total, sd = 1.0)

cat(sprintf("AR(1) DGP: rho = %.2f, sigma = %.2f, T = %d, N = %d\n",
            rho_ar1, sigma_eps, T_total, N_total))

# ===========================================================================
# 2. CATT list and fast lookup
# ===========================================================================

catt_df <- rbindlist(lapply(treat_times, function(g) {
  data.table(g = g, t = g:T_total,
             true_val = sapply(g:T_total, function(t) true_catt_val(g, t)))
}))
catt_df[, catt_idx := .I]
n_catt <- nrow(catt_df)

catt_key  <- catt_df$g * 100L + catt_df$t
catt_lkup <- integer(max(catt_key) + 1L)
catt_lkup[catt_key + 1L] <- catt_df$catt_idx

get_cidx <- function(g, t) {
  k <- as.integer(g) * 100L + as.integer(t)
  if (k < 1L || (k + 1L) > length(catt_lkup)) return(NA_integer_)
  v <- catt_lkup[k + 1L]
  if (v == 0L) NA_integer_ else v
}

# ===========================================================================
# 3. Enumerate all 2x2 DiDs and build Q_H
# ===========================================================================

did_rows <- list()
for (focal_g in treat_times) {
  for (t_post in focal_g:T_total) {
    ci <- get_cidx(focal_g, t_post)
    for (m in seq_len(focal_g - 1L)) {
      t_pre <- focal_g - m
      # Never-treated (always clean)
      did_rows[[length(did_rows) + 1]] <- list(
        focal_g = focal_g, ctrl_g = 0L, t_post = t_post, t_pre = t_pre,
        catt_idx = ci, type = "clean")
      for (ctrl_g in treat_times) {
        if (ctrl_g == focal_g) next
        if (ctrl_g > t_post) {               # not-yet-treated: clean
          did_rows[[length(did_rows) + 1]] <- list(
            focal_g = focal_g, ctrl_g = ctrl_g, t_post = t_post, t_pre = t_pre,
            catt_idx = ci, type = "clean")
        } else if (ctrl_g < focal_g) {        # already-treated: forbidden
          did_rows[[length(did_rows) + 1]] <- list(
            focal_g = focal_g, ctrl_g = ctrl_g, t_post = t_post, t_pre = t_pre,
            catt_idx = ci, type = "forbidden")
        }
      }
    }
  }
}

did_dt    <- rbindlist(did_rows)
n_did     <- nrow(did_dt)
is_clean  <- did_dt$type == "clean"
clean_idx <- which(is_clean)
n_clean   <- length(clean_idx)

cat(sprintf("DiDs: %d total  (%d clean,  %d forbidden)\n",
            n_did, n_clean, n_did - n_clean))

# Build Q_H [n_did x n_catt]
Q_H <- matrix(0.0, nrow = n_did, ncol = n_catt)
for (s in seq_len(n_did)) {
  d <- did_dt[s]
  Q_H[s, d$catt_idx] <- 1.0
  if (d$type == "forbidden") {
    ni <- get_cidx(d$ctrl_g, d$t_post)
    pi <- get_cidx(d$ctrl_g, d$t_pre)
    if (!is.na(ni)) Q_H[s, ni] <- Q_H[s, ni] - 1.0
    if (!is.na(pi)) Q_H[s, pi] <- Q_H[s, pi] + 1.0
  }
}

# ===========================================================================
# 4. Analytical Omega_phi under AR(1) panel errors
#
#    Each DiD s = Ybar_{fg,tp} - Ybar_{fg,tr} - Ybar_{cg,tp} + Ybar_{cg,tr}
#    Under AR(1): Cov(Ybar_{g,t}, Ybar_{g',t'}) = (sigma^2 * rho^|t-t'| / N_g) * 1{g=g'}
#
#    => Cov(Delta_s, Delta_s') = sigma^2 *
#       sum_{(g,t) in s, (g,t') in s', g shared}
#         (sign_s[g,t]/N_g) * (sign_s'[g,t']/1) * rho^|t-t'|
#
#    which simplifies to:
#       sigma^2 * sum_g sum_{ta in terms_s, tb in terms_s': ta$g==tb$g}
#         ta$w * tb$w * N_g * rho^|ta$t - tb$t|
#
#    Setting rho=0 recovers the iid formula (rho^0 = 1 => only same-time terms).
#    Setting rho=0 with same-time only recovers build_Omega_phi exactly.
# ===========================================================================

build_Omega_phi_ar1 <- function(rho = 0.0, sigma2 = 1.0) {
  get_n <- function(g) if (g == 0L) n_never else cohort_size

  Omega <- matrix(0.0, nrow = n_did, ncol = n_did)
  for (s in seq_len(n_did)) {
    ds <- did_dt[s]
    terms_s <- list(
      list(g = ds$focal_g, t = ds$t_post,  w = +1.0 / get_n(ds$focal_g)),
      list(g = ds$focal_g, t = ds$t_pre,   w = -1.0 / get_n(ds$focal_g)),
      list(g = ds$ctrl_g,  t = ds$t_post,  w = -1.0 / get_n(ds$ctrl_g)),
      list(g = ds$ctrl_g,  t = ds$t_pre,   w = +1.0 / get_n(ds$ctrl_g))
    )
    for (sp in s:n_did) {
      dsp <- did_dt[sp]
      terms_sp <- list(
        list(g = dsp$focal_g, t = dsp$t_post,  w = +1.0 / get_n(dsp$focal_g)),
        list(g = dsp$focal_g, t = dsp$t_pre,   w = -1.0 / get_n(dsp$focal_g)),
        list(g = dsp$ctrl_g,  t = dsp$t_post,  w = -1.0 / get_n(dsp$ctrl_g)),
        list(g = dsp$ctrl_g,  t = dsp$t_pre,   w = +1.0 / get_n(dsp$ctrl_g))
      )
      cov_ss <- 0.0
      for (ta in terms_s) {
        for (tb in terms_sp) {
          if (ta$g == tb$g) {
            # Cov(Ybar_{g,t}, Ybar_{g,t'}) = sigma2 * rho^|t-t'| / N_g
            # ta$w and tb$w already include 1/N_g, so multiply back by N_g once
            cov_ss <- cov_ss + ta$w * tb$w * get_n(ta$g) * rho^abs(ta$t - tb$t)
          }
        }
      }
      Omega[s, sp] <- Omega[sp, s] <- sigma2 * cov_ss
    }
  }
  Omega
}

cat(sprintf("Computing analytical Omega_phi under AR(1) errors (rho = %.2f)...\n", rho_ar1))
Omega_phi_ar1 <- build_Omega_phi_ar1(rho = rho_ar1, sigma2 = sigma_eps^2)

# Verify: rho=0 should match the iid formula from the original file.
# (Uncomment to check during development)
# Omega_phi_iid <- build_Omega_phi_ar1(rho = 0.0, sigma2 = sigma_eps^2)
# stopifnot(max(abs(Omega_phi_iid - build_Omega_phi(sigma2 = sigma_eps^2))) < 1e-12)
# cat("  Verification: rho=0 matches iid formula. OK.\n")

# ===========================================================================
# 5. Pre-compute GMM matrices
# ===========================================================================

# EqWt_Clean: A = I restricted to clean rows (same as original)
Q_clean   <- Q_H[clean_idx, , drop = FALSE]
QcQc_inv  <- solve(t(Q_clean) %*% Q_clean)
L_clean   <- QcQc_inv %*% t(Q_clean)

# EqWt_All: A = I over all rows (same as original)
QhQh_inv  <- solve(t(Q_H) %*% Q_H)
L_all     <- QhQh_inv %*% t(Q_H)

stopifnot(max(abs(L_clean %*% Q_clean - diag(n_catt))) < 1e-8)
stopifnot(max(abs(L_all   %*% Q_H     - diag(n_catt))) < 1e-8)

# EqWt_All analytical sandwich variance uses the AR(1) Omega_phi
AVar_clean <- L_clean %*% Omega_phi_ar1[clean_idx, clean_idx] %*% t(L_clean)
AVar_all   <- L_all   %*% Omega_phi_ar1                       %*% t(L_all)

# AR1_Clean: A = pinv(Omega_cl^AR1) — optimal GMM on clean DiDs under AR(1)
Omega_clean_ar1   <- Omega_phi_ar1[clean_idx, clean_idx]
Omega_cl_ar1_pinv <- ginv(Omega_clean_ar1)
QcOQc_ar1         <- t(Q_clean) %*% Omega_cl_ar1_pinv %*% Q_clean
L_ar1_clean       <- solve(QcOQc_ar1) %*% t(Q_clean) %*% Omega_cl_ar1_pinv
AVar_ar1_clean    <- L_ar1_clean %*% Omega_clean_ar1 %*% t(L_ar1_clean)
stopifnot(max(abs(L_ar1_clean %*% Q_clean - diag(n_catt))) < 1e-8)

# AR1_All: A = pinv(Omega_phi^AR1) — optimal GMM on all DiDs under AR(1)
Omega_phi_ar1_pinv <- ginv(Omega_phi_ar1)
QhOQh_ar1          <- t(Q_H) %*% Omega_phi_ar1_pinv %*% Q_H
QhOQh_ar1_inv      <- solve(QhOQh_ar1)
L_ar1_all          <- QhOQh_ar1_inv %*% t(Q_H) %*% Omega_phi_ar1_pinv
AVar_ar1_all       <- L_ar1_all %*% Omega_phi_ar1 %*% t(L_ar1_all)

cat("GMM matrices computed.\n")

# ===========================================================================
# 6. Helper: compute cohort-time mean outcomes -> Delta vector
# ===========================================================================

compute_Delta <- function(dt_sim) {
  gp <- dt_sim[, .(ybar = mean(y)), keyby = .(cohort, time)]
  Delta <- numeric(n_did)
  for (s in seq_len(n_did)) {
    d <- did_dt[s]
    yf_post <- gp[.(d$focal_g, d$t_post), ybar]
    yf_pre  <- gp[.(d$focal_g, d$t_pre),  ybar]
    yc_post <- gp[.(d$ctrl_g,  d$t_post), ybar]
    yc_pre  <- gp[.(d$ctrl_g,  d$t_pre),  ybar]
    v <- c(yf_post, yf_pre, yc_post, yc_pre)
    Delta[s] <- if (anyNA(v)) NA_real_ else (yf_post - yf_pre) - (yc_post - yc_pre)
  }
  Delta
}

# ===========================================================================
# 7. Gardner two-stage via feols
# ===========================================================================

estimate_gardner <- function(dt_sim) {
  untreated <- dt_sim[cohort == 0L | time < cohort]
  m1  <- feols(y ~ 0 | id + time, data = untreated, warn = FALSE, notes = FALSE)
  fe  <- fixef(m1)
  alpha_hat  <- fe$id[as.character(dt_sim$id)]
  lambda_hat <- fe$time[as.character(dt_sim$time)]
  dt_sim[, resid := y - alpha_hat - lambda_hat]
  dt_sim[cohort > 0L & time >= cohort, .(catt = mean(resid)), keyby = .(cohort, time)]
}

# ===========================================================================
# 8. Simulation (100 draws)
# ===========================================================================

n_sims   <- 100L
innov_sd <- sigma_eps * sqrt(1 - rho_ar1^2)   # innovation SD for stationary AR(1)

res_gard      <- matrix(NA_real_, nrow = n_sims, ncol = n_catt)
res_clean     <- matrix(NA_real_, nrow = n_sims, ncol = n_catt)
res_ar1_clean <- matrix(NA_real_, nrow = n_sims, ncol = n_catt)
res_all       <- matrix(NA_real_, nrow = n_sims, ncol = n_catt)
res_ar1_all   <- matrix(NA_real_, nrow = n_sims, ncol = n_catt)

# Per-sim diagnostic: max|AR1_Clean - AR1_All| to test redundancy under AR(1)
diff_ar1_per_sim <- numeric(n_sims)

cat(sprintf("Running %d simulations with rho = %.2f...\n", n_sims, rho_ar1))
set.seed(2025)

for (sim in seq_len(n_sims)) {

  # Generate stationary AR(1) errors: eps_{i,1} ~ N(0, sigma^2) marginal start
  eps_mat <- matrix(0.0, nrow = N_total, ncol = T_total)
  eps_mat[, 1L] <- rnorm(N_total, sd = sigma_eps)
  for (tt in 2L:T_total) {
    eps_mat[, tt] <- rho_ar1 * eps_mat[, tt - 1L] + rnorm(N_total, sd = innov_sd)
  }

  dt_sim <- data.table(
    id     = rep(unit_ids, each = T_total),
    time   = rep(seq_len(T_total), times = N_total),
    cohort = rep(unit_cohort, each = T_total)
  )
  dt_sim[, y := alpha_i[id] + lambda_t[time] +
           mapply(true_catt_val, cohort, time) +
           eps_mat[cbind(id, time)]]

  # Gardner (feols 2-stage)
  catts_g <- estimate_gardner(dt_sim)
  for (k in seq_len(n_catt)) {
    v <- catts_g[.(catt_df$g[k], catt_df$t[k]), catt]
    if (length(v) == 1L && !is.na(v)) res_gard[sim, k] <- v
  }

  # Delta vectors
  Delta_all <- compute_Delta(dt_sim)
  if (anyNA(Delta_all)) next
  Delta_clean <- Delta_all[clean_idx]

  # GMM estimators
  res_clean[sim, ]     <- as.vector(L_clean     %*% Delta_clean)
  res_ar1_clean[sim, ] <- as.vector(L_ar1_clean %*% Delta_clean)
  res_all[sim, ]       <- as.vector(L_all       %*% Delta_all)
  res_ar1_all[sim, ]   <- as.vector(L_ar1_all   %*% Delta_all)

  # Sim-by-sim redundancy diagnostic
  diff_ar1_per_sim[sim] <- max(abs(res_ar1_clean[sim, ] - res_ar1_all[sim, ]))

  if (sim %% 20L == 0L) cat(sprintf("  sim %d / %d\n", sim, n_sims))
}

# ===========================================================================
# 9. Aggregate to cohort-level ATT
# ===========================================================================

agg_cohort <- function(mat) {
  sapply(seq_along(treat_times), function(ci) {
    rowMeans(mat[, catt_df$g == treat_times[ci], drop = FALSE], na.rm = TRUE)
  })
}

true_cohort  <- sapply(treat_times, function(g) mean(catt_df$true_val[catt_df$g == g]))
n_post       <- sapply(treat_times, function(g) sum(catt_df$g == g))
wts          <- n_post / sum(n_post)
true_overall <- sum(true_cohort * wts)

ca_gard      <- agg_cohort(res_gard)
ca_clean     <- agg_cohort(res_clean)
ca_ar1_clean <- agg_cohort(res_ar1_clean)
ca_all       <- agg_cohort(res_all)
ca_ar1_all   <- agg_cohort(res_ar1_all)

oa <- function(ca) as.vector(ca %*% wts)
oa_gard      <- oa(ca_gard)
oa_clean     <- oa(ca_clean)
oa_ar1_clean <- oa(ca_ar1_clean)
oa_all       <- oa(ca_all)
oa_ar1_all   <- oa(ca_ar1_all)

# Analytical cohort-level variances (diagonal of AVar summed/averaged per cohort)
agg_avar <- function(AV) {
  sapply(seq_along(treat_times), function(ci) {
    idx <- which(catt_df$g == treat_times[ci])
    mean(diag(AV)[idx])
  })
}

av_clean_c     <- agg_avar(AVar_clean)
av_ar1_clean_c <- agg_avar(AVar_ar1_clean)
av_all_c       <- agg_avar(AVar_all)
av_ar1_all_c   <- agg_avar(AVar_ar1_all)

c_oa <- rep(1.0 / n_catt, n_catt)
av_clean_oa     <- as.numeric(t(c_oa) %*% AVar_clean     %*% c_oa)
av_ar1_clean_oa <- as.numeric(t(c_oa) %*% AVar_ar1_clean %*% c_oa)
av_all_oa       <- as.numeric(t(c_oa) %*% AVar_all       %*% c_oa)
av_ar1_all_oa   <- as.numeric(t(c_oa) %*% AVar_ar1_all   %*% c_oa)

# ===========================================================================
# 10. Print results
# ===========================================================================

SEP  <- strrep("-", 100)
SEP2 <- strrep("=", 100)

cat("\n\n")
cat(SEP2, "\n", sep = "")
cat(sprintf("AR(1) SIMULATION RESULTS  (rho = %.2f, sigma = %.2f, T = %d, N = %d)\n",
            rho_ar1, sigma_eps, T_total, N_total))
cat(SEP2, "\n", sep = "")

cat("\n", SEP2, "\n", sep = "")
cat("TABLE 1: BIAS — all five estimators (should all be ~0 under parallel trends)\n")
cat(SEP2, "\n", sep = "")
cat(sprintf("%-16s  %8s  %10s  %10s  %12s  %10s  %12s\n",
            "CATT", "True", "Gardner", "EqWt_Cln", "AR1_Clean", "EqWt_All", "AR1_All"))
cat(SEP, "\n", sep = "")
for (ci in seq_along(treat_times)) {
  g  <- treat_times[ci]; tv <- true_cohort[ci]
  cat(sprintf("ATT(g=%2d, avg)  %8.4f  %+10.6f  %+10.6f  %+12.6f  %+10.6f  %+12.6f\n",
              g, tv,
              mean(ca_gard[, ci],      na.rm = TRUE) - tv,
              mean(ca_clean[, ci],     na.rm = TRUE) - tv,
              mean(ca_ar1_clean[, ci], na.rm = TRUE) - tv,
              mean(ca_all[, ci],       na.rm = TRUE) - tv,
              mean(ca_ar1_all[, ci],   na.rm = TRUE) - tv))
}
cat(SEP, "\n", sep = "")
cat(sprintf("%-16s  %8.4f  %+10.6f  %+10.6f  %+12.6f  %+10.6f  %+12.6f\n",
            "Overall ATT", true_overall,
            mean(oa_gard,      na.rm = TRUE) - true_overall,
            mean(oa_clean,     na.rm = TRUE) - true_overall,
            mean(oa_ar1_clean, na.rm = TRUE) - true_overall,
            mean(oa_all,       na.rm = TRUE) - true_overall,
            mean(oa_ar1_all,   na.rm = TRUE) - true_overall))

cat("\n", SEP2, "\n", sep = "")
cat("TABLE 2: VARIANCE — empirical (100 sims) and analytical (sandwich formula)\n")
cat("  Under AR(1), Gardner is no longer equivalent to the spherical-optimal estimator.\n")
cat(SEP2, "\n", sep = "")
cat(sprintf("%-16s  %-9s  %-9s  %-12s  %-9s  %-12s  %-12s\n",
            "CATT",
            "Grd(emp)", "Cln(emp)", "AR1Cln(emp)", "All(emp)",
            "AR1Cln(th)", "AR1All(th)"))
cat(SEP, "\n", sep = "")
for (ci in seq_along(treat_times)) {
  g      <- treat_times[ci]
  vg     <- var(ca_gard[, ci],      na.rm = TRUE)
  vc     <- var(ca_clean[, ci],     na.rm = TRUE)
  vac    <- var(ca_ar1_clean[, ci], na.rm = TRUE)
  va     <- var(ca_all[, ci],       na.rm = TRUE)
  th_ac  <- av_ar1_clean_c[ci]
  th_ao  <- av_ar1_all_c[ci]
  cat(sprintf(
    "ATT(g=%2d, avg)  %9.5f  %9.5f  %12.5f  %9.5f  %10.5f  %10.5f\n",
    g, vg, vc, vac, va, th_ac, th_ao))
}
cat(SEP, "\n", sep = "")
vg_o   <- var(oa_gard,      na.rm = TRUE)
vc_o   <- var(oa_clean,     na.rm = TRUE)
vac_o  <- var(oa_ar1_clean, na.rm = TRUE)
va_o   <- var(oa_all,       na.rm = TRUE)
cat(sprintf(
  "%-16s  %9.5f  %9.5f  %12.5f  %9.5f  %10.5f  %10.5f\n",
  "Overall ATT", vg_o, vc_o, vac_o, va_o, av_ar1_clean_oa, av_ar1_all_oa))

cat("\n", SEP2, "\n", sep = "")
cat("TABLE 3: ANALYTICAL VARIANCE RATIOS (under AR(1) errors, rho =",
    rho_ar1, ")\n")
cat("  Ratio > 1 means numerator has HIGHER variance (denominator more efficient).\n")
cat("  AR1Cln/AR1All = 1.000 iff redundancy holds under AR(1).\n")
cat(SEP2, "\n", sep = "")
cat(sprintf("%-24s  %14s  %14s  %14s  %14s\n",
            "CATT", "Cln/AR1Cln", "Cln/All", "AR1Cln/AR1All", "Cln/AR1All"))
cat(SEP, "\n", sep = "")

print_ratios <- function(k, label) {
  c_ac  <- AVar_clean[k, k]     / AVar_ar1_clean[k, k]
  c_a   <- AVar_clean[k, k]     / AVar_all[k, k]
  ac_ao <- AVar_ar1_clean[k, k] / AVar_ar1_all[k, k]
  c_ao  <- AVar_clean[k, k]     / AVar_ar1_all[k, k]
  cat(sprintf("%-24s  %14.4f  %14.4f  %14.4f  %14.4f\n",
              label, c_ac, c_a, ac_ao, c_ao))
}

for (ci in seq_along(treat_times)) {
  g  <- treat_times[ci]
  k0 <- which(catt_df$g == g & catt_df$t == g)
  print_ratios(k0, sprintf("ATT(g=%2d, t=%2d) k=0", g, g))
}
for (ci in seq_along(treat_times)) {
  g      <- treat_times[ci]
  k_last <- which(catt_df$g == g & catt_df$t == T_total)
  print_ratios(k_last, sprintf("ATT(g=%2d, t=%2d) k=%2d", g, T_total, T_total - g))
}
cat(SEP, "\n", sep = "")

c_ac_oa  <- av_clean_oa     / av_ar1_clean_oa
c_a_oa   <- av_clean_oa     / av_all_oa
ac_ao_oa <- av_ar1_clean_oa / av_ar1_all_oa
c_ao_oa  <- av_clean_oa     / av_ar1_all_oa
cat(sprintf("%-24s  %14.4f  %14.4f  %14.4f  %14.4f\n",
            "Overall ATT", c_ac_oa, c_a_oa, ac_ao_oa, c_ao_oa))
cat(SEP2, "\n", sep = "")

# ===========================================================================
# 11. Key diagnostic: Does AR1_Clean = AR1_All under AR(1) errors?
# ===========================================================================

cat("\n", SEP2, "\n", sep = "")
cat("DIAGNOSTIC: Redundancy under AR(1) — Are AR1_Clean and AR1_All identical?\n")
cat("  Under spherical errors: Sph_Clean = Sph_All exactly (Redundancy Proposition).\n")
cat("  Proof relies on cross-cohort independence (Assumption 1), not iid errors.\n")
cat("  AR(1) preserves within-unit serial correlation but keeps units in\n")
cat("  different cohorts independent => Assumption 1 still holds.\n")
cat("  Prediction: AR1_Clean = AR1_All under AR(1) as well.\n")
cat(SEP2, "\n", sep = "")

cat(sprintf("\n  Max |AR1_Clean - AR1_All| per simulation:\n"))
cat(sprintf("    Mean across %d sims: %.2e\n", n_sims, mean(diff_ar1_per_sim)))
cat(sprintf("    Median:              %.2e\n",          median(diff_ar1_per_sim)))
cat(sprintf("    Max:                 %.2e\n",          max(diff_ar1_per_sim)))
cat(sprintf("    Min:                 %.2e\n",          min(diff_ar1_per_sim)))
cat(sprintf("\n  CONCLUSION: AR1_Clean %s AR1_All under AR(1) (rho = %.2f)\n",
            if (max(diff_ar1_per_sim) < 1e-6) "IS IDENTICAL TO" else "DIFFERS FROM",
            rho_ar1))

cat(sprintf("\n  Analytical check — Var ratio AR1Cln/AR1All (should = 1 if redundant):\n"))
cat(sprintf("    Overall ATT: %.6f\n", ac_ao_oa))
ratio_range <- range(AVar_ar1_clean[cbind(seq_len(n_catt), seq_len(n_catt))] /
                     AVar_ar1_all[cbind(seq_len(n_catt), seq_len(n_catt))])
cat(sprintf("    CATT range:  [%.6f, %.6f]\n", ratio_range[1], ratio_range[2]))

cat(SEP2, "\n", sep = "")
cat("\nSummary:\n")
cat("  Cln/AR1Cln:    gain from AR(1)-optimal weighting within clean DiDs.\n")
cat("  Cln/All:       effect of adding forbidden DiDs with A=I.\n")
cat("  AR1Cln/AR1All: marginal value of forbidden DiDs with AR(1)-optimal A.\n")
cat("                 = 1 confirms redundancy holds; > 1 means they add information.\n")
cat("  Cln/AR1All:    total gain from optimal weighting + forbidden DiDs.\n")
