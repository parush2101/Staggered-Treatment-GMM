###############################################################################
# gardner_forbidden_efficiency.R
#
# Question: Is extended Gardner — Gardner with forbidden comparisons included
# (debiased via Q_H, A = I) — more efficient than original Gardner?
# And: can a spherical-errors-optimal A on clean DiDs beat equal weighting?
#
# Five estimators compared:
#   1. Gardner     : Two-stage feols (stage 1 on untreated, stage 2 averages
#                    residuals per cohort×time). Standard Gardner estimator.
#   2. EqWt_Clean  : GMM with A = I restricted to clean DiDs only.
#                    hat_beta = (Q_clean'Q_clean)^{-1} Q_clean' Delta_clean.
#                    NOTE: differs from Gardner due to cross-CATT coupling in
#                    Gardner's joint FE regression (paper Section 3c).
#   3. Sph_Clean   : Optimal GMM on clean DiDs under spherical (iid) panel
#                    errors.  A = pinv(Omega_clean) where Omega_clean is the
#                    analytically-computed covariance of clean DiDs under iid
#                    errors. Efficient among all clean-DiD-only estimators.
#                    Accounts for DiD-level correlations (shared cohort-time
#                    cells) without requiring any data-based error estimation.
#   4. EqWt_All    : GMM with A = I over ALL DiDs (clean + forbidden).
#                    hat_beta = (Q_H'Q_H)^{-1} Q_H' Delta_all.
#                    Forbidden moments are unbiased via Q_H off-diagonals.
#                    This is "Extended Gardner": adds forbidden comparisons.
#   5. Opt_All     : Optimal GMM: A = pinv(Omega_phi) over ALL DiDs.
#                    Efficient among all linear unbiased estimators under iid.
#
# 2x2 table of (DiD set) x (A matrix):
#   Clean / A=I         => EqWt_Clean
#   Clean / A=Omega^+   => Sph_Clean   [NEW]
#   All   / A=I         => EqWt_All
#   All   / A=Omega^+   => Opt_All
#
# Key efficiency comparisons (variance ratios, ratio > 1 means denominator
# estimator is more efficient):
#   EqWt_Clean / Sph_Clean — gain from optimal weighting within clean DiDs
#   EqWt_Clean / EqWt_All  — gain from adding forbidden DiDs (A=I fixed)
#   Sph_Clean  / Opt_All   — marginal cost of excluding forbidden DiDs
#                            when both use the spherical-optimal A
#   EqWt_All   / Opt_All   — gain from optimal vs equal weighting (all DiDs)
#
# DGP: Y_it = alpha_i + lambda_t + tau_it + eps_it,  eps_it ~ iid N(0,1)
# Design: 3 cohorts (g = 4, 7, 10), T = 12, 10 units/cohort + 10 never-treated
# Heterogeneous treatment effects: tau_it = beta_g * (1 + r_g)^{t - g}
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
# 4. Pre-compute GMM matrices
# ===========================================================================

# EqWt_Clean: A = I restricted to clean rows
Q_clean   <- Q_H[clean_idx, , drop = FALSE]
QcQc_inv  <- solve(t(Q_clean) %*% Q_clean)   # n_catt x n_catt
L_clean   <- QcQc_inv %*% t(Q_clean)          # n_catt x n_clean

# EqWt_All: A = I over all rows
QhQh_inv  <- solve(t(Q_H) %*% Q_H)            # n_catt x n_catt
L_all     <- QhQh_inv %*% t(Q_H)              # n_catt x n_did

# Verify unbiasedness: L * Q = I
stopifnot(max(abs(L_clean %*% Q_clean - diag(n_catt))) < 1e-8)
stopifnot(max(abs(L_all   %*% Q_H     - diag(n_catt))) < 1e-8)

# ===========================================================================
# 5. Analytical Omega_phi under iid panel errors
#    Cov(Delta_s, Delta_s') = sigma^2 * phi(s, s')
#    phi(s,s') = sum over {focal/ctrl} x {post/pre} of shared-cell indicators
#    (with signs: +focal_post, -focal_pre, -ctrl_post, +ctrl_pre)
# ===========================================================================

build_Omega_phi <- function(sigma2 = 1.0) {
  # Each DiD s = Y_{fg,tp} - Y_{fg,tr} - Y_{cg,tp} + Y_{cg,tr}   (cohort means)
  # sign_s[cohort][time] contribution:
  #   focal_g at t_post: +1/N_{fg}
  #   focal_g at t_pre:  -1/N_{fg}
  #   ctrl_g  at t_post: -1/N_{cg}
  #   ctrl_g  at t_pre:  +1/N_{cg}
  # Cov(Y_{g,t}, Y_{g',t'}) = sigma^2/N_g * 1{g=g', t=t'} (iid across units & time)
  # => Cov(Delta_s, Delta_{s'}) = sigma^2 * sum of (sign_s[g,t] * sign_{s'}[g,t] / N_g)
  #    over all (g, t) cells that appear in both

  get_n <- function(g) if (g == 0L) n_never else cohort_size

  Omega <- matrix(0.0, nrow = n_did, ncol = n_did)
  for (s in seq_len(n_did)) {
    ds <- did_dt[s]
    # Signed contributions for DiD s: list of (cohort, time, sign/N)
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
          if (ta$g == tb$g && ta$t == tb$t) {
            cov_ss <- cov_ss + ta$w * tb$w * get_n(ta$g)
          }
        }
      }
      Omega[s, sp] <- Omega[sp, s] <- sigma2 * cov_ss
    }
  }
  Omega
}

cat("Computing analytical Omega_phi under iid errors...\n")
Omega_phi <- build_Omega_phi(sigma2 = sigma_eps^2)

# Omega_phi is rank-deficient: 225 DiDs span only 48 cohort-time cells.
# Matrix products L %*% Omega_phi %*% t(L) are well-defined without inversion.
AVar_clean <- L_clean %*% Omega_phi[clean_idx, clean_idx] %*% t(L_clean)
AVar_all   <- L_all   %*% Omega_phi                       %*% t(L_all)

# Sph_Clean: A = pinv(Omega_clean) — optimal GMM on clean DiDs under iid errors.
# Omega_clean is the submatrix of Omega_phi for clean DiDs only.
# pinv(Omega_clean) != submatrix of pinv(Omega_phi): genuinely different from Opt_All.
# This A is design-based: determined entirely by panel structure, no data estimation.
Omega_clean    <- Omega_phi[clean_idx, clean_idx]
Omega_cl_pinv  <- ginv(Omega_clean)
QcOQc          <- t(Q_clean) %*% Omega_cl_pinv %*% Q_clean   # n_catt x n_catt (full rank)
L_sph_clean    <- solve(QcOQc) %*% t(Q_clean) %*% Omega_cl_pinv
AVar_sph_clean <- L_sph_clean %*% Omega_clean %*% t(L_sph_clean)
stopifnot(max(abs(L_sph_clean %*% Q_clean - diag(n_catt))) < 1e-8)

# Opt_All: A = Omega_phi^+ (Moore-Penrose pseudoinverse).
# The pseudoinverse is well-defined on the column space of Q_H.
Omega_pinv <- ginv(Omega_phi)
QhOQh      <- t(Q_H) %*% Omega_pinv %*% Q_H
QhOQh_inv  <- solve(QhOQh)
L_opt      <- QhOQh_inv %*% t(Q_H) %*% Omega_pinv
AVar_opt   <- L_opt %*% Omega_phi %*% t(L_opt)

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

n_sims <- 100L

res_gard      <- matrix(NA_real_, nrow = n_sims, ncol = n_catt)
res_clean     <- matrix(NA_real_, nrow = n_sims, ncol = n_catt)
res_sph_clean <- matrix(NA_real_, nrow = n_sims, ncol = n_catt)
res_all       <- matrix(NA_real_, nrow = n_sims, ncol = n_catt)
res_opt       <- matrix(NA_real_, nrow = n_sims, ncol = n_catt)

cat(sprintf("Running %d simulations...\n", n_sims))
set.seed(2025)

for (sim in seq_len(n_sims)) {
  dt_sim <- data.table(
    id     = rep(unit_ids, each = T_total),
    time   = rep(seq_len(T_total), times = N_total),
    cohort = rep(unit_cohort, each = T_total)
  )
  dt_sim[, y := alpha_i[id] + lambda_t[time] +
           mapply(true_catt_val, cohort, time) +
           rnorm(.N, sd = sigma_eps)]

  # Gardner (feols 2-stage)
  catts_g <- estimate_gardner(dt_sim)
  for (k in seq_len(n_catt)) {
    v <- catts_g[.(catt_df$g[k], catt_df$t[k]), catt]
    if (length(v) == 1L && !is.na(v)) res_gard[sim, k] <- v
  }

  # Delta vectors
  Delta_all   <- compute_Delta(dt_sim)
  if (anyNA(Delta_all)) next
  Delta_clean <- Delta_all[clean_idx]

  # EqWt_Clean, Sph_Clean, EqWt_All, Opt_All
  res_clean[sim, ]     <- as.vector(L_clean     %*% Delta_clean)
  res_sph_clean[sim, ] <- as.vector(L_sph_clean %*% Delta_clean)
  res_all[sim, ]       <- as.vector(L_all       %*% Delta_all)
  res_opt[sim, ]       <- as.vector(L_opt       %*% Delta_all)

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

true_cohort <- sapply(treat_times, function(g) mean(catt_df$true_val[catt_df$g == g]))
n_post      <- sapply(treat_times, function(g) sum(catt_df$g == g))
wts         <- n_post / sum(n_post)
true_overall <- sum(true_cohort * wts)

ca_gard      <- agg_cohort(res_gard)
ca_clean     <- agg_cohort(res_clean)
ca_sph_clean <- agg_cohort(res_sph_clean)
ca_all       <- agg_cohort(res_all)
ca_opt       <- agg_cohort(res_opt)

oa <- function(ca) as.vector(ca %*% wts)
oa_gard      <- oa(ca_gard)
oa_clean     <- oa(ca_clean)
oa_sph_clean <- oa(ca_sph_clean)
oa_all       <- oa(ca_all)
oa_opt       <- oa(ca_opt)

# Analytical cohort-level variances (diagonal of AVar summed/averaged per cohort)
agg_avar <- function(AV) {
  sapply(seq_along(treat_times), function(ci) {
    idx <- which(catt_df$g == treat_times[ci])
    mean(diag(AV)[idx])
  })
}

av_clean_c     <- agg_avar(AVar_clean)
av_sph_clean_c <- agg_avar(AVar_sph_clean)
av_all_c       <- agg_avar(AVar_all)
av_opt_c       <- agg_avar(AVar_opt)

# Overall ATT aggregation weights at CATT level: equal weight 1/n_catt
c_oa <- rep(1.0 / n_catt, n_catt)
av_clean_oa     <- as.numeric(t(c_oa) %*% AVar_clean     %*% c_oa)
av_sph_clean_oa <- as.numeric(t(c_oa) %*% AVar_sph_clean %*% c_oa)
av_all_oa       <- as.numeric(t(c_oa) %*% AVar_all       %*% c_oa)
av_opt_oa       <- as.numeric(t(c_oa) %*% AVar_opt       %*% c_oa)

# ===========================================================================
# 10. Print results
# ===========================================================================

SEP  <- strrep("-", 100)
SEP2 <- strrep("=", 100)

cat("\n\n")
cat(SEP2, "\n", sep = "")
cat("TABLE 1: CROSS-CATT COUPLING — Gardner (feols) vs EqWt_Clean\n")
cat("  Gardner's joint FE regression creates cross-cohort coupling (paper Section 3c);\n")
cat("  EqWt_Clean is OLS on clean DiDs only (no coupling).  They differ when\n")
cat("  not-yet-treated cohorts appear as clean controls (g=4 and g=7 here);\n")
cat("  they agree when only never-treated controls exist (g=10 here).\n")
cat(SEP2, "\n", sep = "")
cat(sprintf("%-18s  %12s  %12s  %16s\n",
            "Cohort ATT", "Gardner", "EqWt_Clean", "Max|diff|/sim"))
cat(SEP, "\n", sep = "")
for (ci in seq_along(treat_times)) {
  g    <- treat_times[ci]
  bg   <- mean(ca_gard[, ci],  na.rm = TRUE) - true_cohort[ci]
  bc   <- mean(ca_clean[, ci], na.rm = TRUE) - true_cohort[ci]
  dif  <- max(abs(ca_gard[, ci] - ca_clean[, ci]), na.rm = TRUE)
  coupling <- if (g == 10L) "(only never-treated ctrl: no coupling)" else
              "(not-yet-treated ctrl present: coupling exists)"
  cat(sprintf("ATT(g=%2d)        bias=%+.4f  bias=%+.4f  max|diff|=%.2e  %s\n",
              g, bg, bc, dif, coupling))
}
cat(SEP, "\n", sep = "")

cat("\n", SEP2, "\n", sep = "")
cat("TABLE 2: BIAS — all five estimators (should all be ~0 under parallel trends)\n")
cat(SEP2, "\n", sep = "")
cat(sprintf("%-16s  %8s  %10s  %12s  %12s  %12s  %12s\n",
            "CATT", "True", "Gardner", "EqWt_Clean", "Sph_Clean", "EqWt_All", "Opt_All"))
cat(SEP, "\n", sep = "")
for (ci in seq_along(treat_times)) {
  g  <- treat_times[ci]; tv <- true_cohort[ci]
  cat(sprintf("ATT(g=%2d, avg)  %8.4f  %+10.6f  %+12.6f  %+12.6f  %+12.6f  %+12.6f\n",
              g, tv,
              mean(ca_gard[,ci],      na.rm=TRUE) - tv,
              mean(ca_clean[,ci],     na.rm=TRUE) - tv,
              mean(ca_sph_clean[,ci], na.rm=TRUE) - tv,
              mean(ca_all[,ci],       na.rm=TRUE) - tv,
              mean(ca_opt[,ci],       na.rm=TRUE) - tv))
}
cat(SEP, "\n", sep = "")
cat(sprintf("%-16s  %8.4f  %+10.6f  %+12.6f  %+12.6f  %+12.6f  %+12.6f\n",
            "Overall ATT", true_overall,
            mean(oa_gard,      na.rm=TRUE) - true_overall,
            mean(oa_clean,     na.rm=TRUE) - true_overall,
            mean(oa_sph_clean, na.rm=TRUE) - true_overall,
            mean(oa_all,       na.rm=TRUE) - true_overall,
            mean(oa_opt,       na.rm=TRUE) - true_overall))

cat("\n", SEP2, "\n", sep = "")
cat("TABLE 3: VARIANCE — empirical (100 sims) and analytical (sandwich formula)\n")
cat("  Sph_Clean uses A = pinv(Omega_clean): optimal weighting within clean DiDs.\n")
cat("  Ratio > 1 means denominator estimator is more efficient.\n")
cat("  Columns: [ratio] = Clean(th) / that column's analytical variance.\n")
cat(SEP2, "\n", sep = "")
cat(sprintf("%-16s  %-9s  %-9s  %-12s  %-12s  %-9s  %-12s  %-12s\n",
            "CATT",
            "Grd(emp)", "Cln(emp)", "Sph_Cln(emp)", "All(emp)",
            "Cln(th)", "SphCln(th)", "Opt(th)"))
cat(SEP, "\n", sep = "")
for (ci in seq_along(treat_times)) {
  g      <- treat_times[ci]
  vg     <- var(ca_gard[, ci],      na.rm = TRUE)
  vc     <- var(ca_clean[, ci],     na.rm = TRUE)
  vsc    <- var(ca_sph_clean[, ci], na.rm = TRUE)
  va     <- var(ca_all[, ci],       na.rm = TRUE)
  th_c   <- av_clean_c[ci]
  th_sc  <- av_sph_clean_c[ci]
  th_o   <- av_opt_c[ci]
  cat(sprintf(
    "ATT(g=%2d, avg)  %9.5f  %9.5f  %12.5f  %12.5f  %9.5f  %10.5f[%5.3f]  %10.5f[%5.3f]\n",
    g, vg, vc, vsc, va, th_c, th_sc, th_c/th_sc, th_o, th_c/th_o))
}
cat(SEP, "\n", sep = "")
vg_o  <- var(oa_gard,      na.rm = TRUE)
vc_o  <- var(oa_clean,     na.rm = TRUE)
vsc_o <- var(oa_sph_clean, na.rm = TRUE)
va_o  <- var(oa_all,       na.rm = TRUE)
cat(sprintf(
  "%-16s  %9.5f  %9.5f  %12.5f  %12.5f  %9.5f  %10.5f[%5.3f]  %10.5f[%5.3f]\n",
  "Overall ATT", vg_o, vc_o, vsc_o, va_o,
  av_clean_oa, av_sph_clean_oa, av_clean_oa/av_sph_clean_oa,
  av_opt_oa, av_clean_oa/av_opt_oa))

cat("\n", SEP2, "\n", sep = "")
cat("TABLE 4: ANALYTICAL VARIANCE RATIOS (exact under iid errors)\n")
cat("  Ratio > 1 means numerator estimator has HIGHER variance (denominator more efficient).\n")
cat("  Clean/SphCln = gain from optimal A within clean DiDs (Sph_Clean over EqWt_Clean)\n")
cat("  Clean/All    = gain from adding forbidden DiDs (A=I fixed)\n")
cat("  SphCln/Opt   = marginal gain from adding forbidden DiDs, both using Omega^+ weighting\n")
cat("  Clean/Opt    = total gain from both optimal weighting AND forbidden DiDs\n")
cat(SEP2, "\n", sep = "")
cat(sprintf("%-24s  %14s  %14s  %14s  %14s\n",
            "CATT", "Clean/SphCln", "Clean/All", "SphCln/Opt", "Clean/Opt"))
cat(SEP, "\n", sep = "")

print_ratios <- function(k, g, label) {
  c_sc <- AVar_clean[k, k]     / AVar_sph_clean[k, k]
  c_a  <- AVar_clean[k, k]     / AVar_all[k, k]
  sc_o <- AVar_sph_clean[k, k] / AVar_opt[k, k]
  c_o  <- AVar_clean[k, k]     / AVar_opt[k, k]
  cat(sprintf("%-24s  %14.4f  %14.4f  %14.4f  %14.4f\n",
              label, c_sc, c_a, sc_o, c_o))
}

for (ci in seq_along(treat_times)) {
  g  <- treat_times[ci]
  k0 <- which(catt_df$g == g & catt_df$t == g)
  print_ratios(k0, g, sprintf("ATT(g=%2d, t=%2d) k=0", g, g))
}
for (ci in seq_along(treat_times)) {
  g      <- treat_times[ci]
  k_last <- which(catt_df$g == g & catt_df$t == T_total)
  print_ratios(k_last, g, sprintf("ATT(g=%2d, t=%2d) k=%2d", g, T_total, T_total - g))
}
cat(SEP, "\n", sep = "")

# Overall ATT ratios
c_sc_oa <- av_clean_oa     / av_sph_clean_oa
c_a_oa  <- av_clean_oa     / av_all_oa
sc_o_oa <- av_sph_clean_oa / av_opt_oa
c_o_oa  <- av_clean_oa     / av_opt_oa
cat(sprintf("%-24s  %14.4f  %14.4f  %14.4f  %14.4f\n",
            "Overall ATT", c_sc_oa, c_a_oa, sc_o_oa, c_o_oa))
cat(SEP2, "\n", sep = "")

cat("\nSummary:\n")
cat("  Clean/SphCln: gain from optimal weighting within clean DiDs.\n")
cat("    At k=0, exchangeable covariance => EqWt_Clean is already optimal => ratio = 1.\n")
cat("    At later horizons, Sph_Clean improves over EqWt_Clean.\n")
cat("  Clean/All: adding forbidden DiDs with A=I can hurt at k=0 (ratio < 1)\n")
cat("    because positive cross-DiD correlation at k=0 makes equal-weighting suboptimal.\n")
cat("  SphCln/Opt: marginal value of forbidden DiDs given optimal A.\n")
cat("    Shows whether forbidden DiDs add independent information beyond clean DiDs.\n")
