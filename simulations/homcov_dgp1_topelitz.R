# GMM I (Toeplitz homosk.) + GMM II (Toeplitz heterosk.) + GMM III (HetCov) | DGP2 (AR(1) only)
# rho in {0,0.1,0.3,0.5,0.7,0.9}  cohort_size=50  n_never=50  n_sims=500 | HETEROGENEOUS
# GMM I and GMM II use CLEAN comparisons only (never-treated + not-yet-treated).

library(data.table)
library(fixest)
library(did)
library(did2s)
library(MASS)
library(parallel)

options(warn = 1)

# Parallelize estimation across simulations. Data are generated sequentially
# (below) so the RNG stream is identical to the serial run; estimation is
# deterministic, so results are unchanged regardless of n_cores.
#
# n_cores controls how many cores run flat-out at once -- the main driver of CPU
# heat/fan on a laptop. Set it explicitly via the N_CORES environment variable:
#   N_CORES=2 Rscript thisfile.R   # cool, slow         (laptop)
#   N_CORES=48 Rscript thisfile.R  # full speed         (cluster)
# If unset, we leave 2 cores free rather than maxing the machine.
n_cores <- {
  v <- suppressWarnings(as.integer(Sys.getenv("N_CORES")))
  if (is.na(v) || v < 1L) max(1L, detectCores() - 2L) else v
}
setFixest_nthreads(1)   # avoid thread oversubscription under forked workers

# Keep BLAS single-threaded so it does not oversubscribe cores under the
# mclapply fork (cores x BLAS-threads). This matters once R is linked against an
# optimized multithreaded BLAS (OpenBLAS/MKL/Accelerate) -- which is the main
# lever for GLS speed: the GLS bottleneck is a dense crossprod that runs ~20-50x
# faster on an optimized BLAS than on the reference libRblas, with no change to
# results (~1e-14). Switch BLAS at the R-installation level, e.g. on macOS:
#   sudo ln -sf libRblas.vecLib.dylib \
#     /Library/Frameworks/R.framework/Resources/lib/libRblas.dylib
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE))
  RhpcBLASctl::blas_set_num_threads(1)

log_file <- "/Users/parusharora/Downloads/hetcov_dgp2_n300_r09_het_01_06.log"
con <- file(log_file, open = "wt")
sink(con, split = TRUE)
sink(con, type = "message")

cat(format(Sys.time()), "\n\n")

set.seed(312844)

cohort_size     <- 10L
n_never         <- 10L
N_total         <- 5L * cohort_size + n_never
T_total         <- 33L
n_sims          <- 500L
treatment_times <- c(10, 13, 16, 19, 22)
n_cohorts       <- length(treatment_times)
unit_cohort     <- c(rep(treatment_times, each = cohort_size), rep(0L, n_never))

rho_grid  <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9)   # AR(1) persistence values to sweep
Sigma_ar1 <- NULL                            # set per rho inside the loop below


build_cells <- function(start_yrs, T_max) {
  rows <- vector("list", sum(T_max - start_yrs + 1L))
  id_c <- 0L
  for (g in start_yrs)
    for (tt in g:T_max) {
      id_c <- id_c + 1L
      rows[[id_c]] <- data.frame(g = as.integer(g), t = as.integer(tt),
                                 col_id = id_c)
    }
  do.call(rbind, rows)
}


# Enumerates CLEAN comparisons only: never-treated (NT) and not-yet-treated
# (NYT). Already-treated (forbidden) comparisons are never formed, so every
# estimator in this script uses clean comparisons exclusively.
enumerate_contrasts_clean <- function(cells, all_cohorts, start_yrs, T_min,
                                      has_nt) {
  buf <- vector("list", 1000000L)
  n   <- 0L
  
  for (idx in seq_len(nrow(cells))) {
    g  <- cells$g[idx]
    tt <- cells$t[idx]
    fc <- cells$col_id[idx]
    
    if (g > T_min) {
      if (has_nt) {
        for (s in T_min:(g - 1L)) {
          n <- n + 1L
          buf[[n]] <- list(g = g, t_post = tt, c = 0L, t_pre = s,
                           type = "NT", focal_col = fc,
                           bias_neg_col = NA_integer_,
                           bias_pos_col = NA_integer_)
        }
      }
      for (h in all_cohorts[all_cohorts > tt & all_cohorts != 0L]) {
        for (s in T_min:(g - 1L)) {
          n <- n + 1L
          buf[[n]] <- list(g = g, t_post = tt, c = as.integer(h), t_pre = s,
                           type = "NYT", focal_col = fc,
                           bias_neg_col = NA_integer_,
                           bias_pos_col = NA_integer_)
        }
      }
    }
  }
  
  buf <- buf[seq_len(n)]
  data.frame(
    g            = vapply(buf, `[[`, integer(1L),   "g"),
    t_post       = vapply(buf, `[[`, integer(1L),   "t_post"),
    c            = vapply(buf, `[[`, integer(1L),   "c"),
    t_pre        = vapply(buf, `[[`, integer(1L),   "t_pre"),
    type         = vapply(buf, `[[`, character(1L), "type"),
    focal_col    = vapply(buf, `[[`, integer(1L),   "focal_col"),
    bias_neg_col = vapply(buf, `[[`, integer(1L),   "bias_neg_col"),
    bias_pos_col = vapply(buf, `[[`, integer(1L),   "bias_pos_col"),
    stringsAsFactors = FALSE
  )
}


build_Q_H <- function(contrasts_df, N_2x2, N_beta) {
  Q_H <- matrix(0, nrow = N_2x2, ncol = N_beta)
  for (r in seq_len(N_2x2)) {
    Q_H[r, contrasts_df$focal_col[r]] <- 1
    if (contrasts_df$type[r] == "AT") {
      bn <- contrasts_df$bias_neg_col[r]
      bp <- contrasts_df$bias_pos_col[r]
      if (!is.na(bn)) Q_H[r, bn] <- Q_H[r, bn] - 1
      if (!is.na(bp)) Q_H[r, bp] <- Q_H[r, bp] + 1
    }
  }
  Q_H
}


# Reduced-space optimal-GMM solve. The N_2x2 moments are linear combinations of
# only the G*T cohort-time means, so Omega_phi = U W U' is rank-deficient with a
# fixed (data-independent) column space spanned by the SVD basis u of U. Because
# u, Q_H and Delta do not change across GMM iterations, the projections
#   uQ = u' Q_H   (m x N_beta),   uDe = u' Delta   (m)
# and the SVD of U are precomputed ONCE by gmm_reduced_factor() below; each
# iteration then works entirely in the m = G*T dimensional space (m << N_2x2),
# never forming an N_2x2-sized object.
#
# With Omega_phi^+ = u M^+ u', M = diag(d) (v' W v) diag(d), the optimal-GMM
# normal-equation pieces are
#   Q_H' Omega_phi^+ Q_H = uQ' M^+ uQ ,   Q_H' Omega_phi^+ Delta = uQ' M^+ uDe ,
# where M^+ is the truncated (generalized) inverse keeping eigen-directions above
# tol_rel * lambda_max. This is algebraically identical to weighting every clean
# comparison by Omega_phi^+ (equivalently, to optimal GMM on a maximal linearly
# independent subset of the comparisons).

# One-off factorization of the fixed incidence factor U.
gmm_reduced_factor <- function(U, Q_H, Delta) {
  sv <- svd(U)
  list(d = sv$d, v = sv$v,
       uQ = crossprod(sv$u, Q_H),               # m x N_beta
       uDe = as.numeric(crossprod(sv$u, Delta)) # m
  )
}

# Per-iteration solve given the precomputed factor `f` and the current W.
# Returns the N_beta-sized normal-equation blocks (QtAQ, QtAD); no N_2x2 work.
gmm_reduced_solve <- function(f, W, tol_rel = 1e-8) {
  M  <- outer(f$d, f$d) * crossprod(f$v, W %*% f$v)   # diag(d) (v'Wv) diag(d)
  M  <- (M + t(M)) / 2
  eg <- eigen(M, symmetric = TRUE)
  lam <- eg$values
  keep <- lam > tol_rel * max(lam)
  if (!any(keep)) return(NULL)
  Vk <- eg$vectors[, keep, drop = FALSE]
  dk <- 1 / lam[keep]
  Mp_uQ  <- Vk %*% (dk * crossprod(Vk, f$uQ))         # M^+ uQ   (m x N_beta)
  Mp_uDe <- Vk %*% (dk * crossprod(Vk, f$uDe))        # M^+ uDe  (m)
  list(QtAQ = crossprod(f$uQ, Mp_uQ),
       QtAD = as.numeric(crossprod(f$uQ, Mp_uDe)))
}


gmm_staggered_hetcov <- function(data, yname, tname, idname, gname,
                                 has_nt   = FALSE,
                                 eig_tol  = 1e-8,
                                 max_iter = 10L,
                                 tol      = 1e-6) {
  
  Y_raw  <- as.numeric(data[[yname]])
  t_raw  <- as.integer(data[[tname]])
  id_raw <- data[[idname]]
  g_raw  <- as.integer(data[[gname]])
  
  unit_ids  <- sort(unique(id_raw))
  N_units   <- length(unit_ids)
  T_min     <- min(t_raw)
  T_max     <- max(t_raw)
  TT        <- T_max - T_min + 1L
  unit_row  <- match(id_raw, unit_ids)
  time_idx  <- t_raw - T_min + 1L
  
  cohort_of   <- as.integer(tapply(g_raw, unit_row, function(x) x[1L]))
  start_yrs   <- sort(unique(cohort_of[cohort_of > 0L & cohort_of <= T_max]))
  K           <- length(start_yrs)
  all_cohorts <- if (has_nt) c(0L, start_yrs) else start_yrs
  
  N_g_lookup <- setNames(
    vapply(start_yrs, function(g) sum(cohort_of == g), integer(1L)),
    as.character(start_yrs)
  )
  if (has_nt)
    N_g_lookup["0"] <- sum(cohort_of == 0L)
  
  Y_mat <- matrix(NA_real_, nrow = N_units, ncol = TT)
  Y_mat[cbind(unit_row, time_idx)] <- Y_raw
  
  dt_long <- data.table(
    unit     = id_raw,
    unit_row = unit_row,
    time     = t_raw,
    time_idx = time_idx,
    Y        = Y_raw,
    g        = g_raw
  )
  setorder(dt_long, unit, time)
  
  cells        <- build_cells(start_yrs, T_max)
  N_beta       <- nrow(cells)
  contrasts_df <- enumerate_contrasts_clean(cells, all_cohorts, start_yrs,
                                            T_min, has_nt)
  N_2x2      <- nrow(contrasts_df)
  meta_focal <- contrasts_df$g
  meta_ctrl  <- contrasts_df$c
  meta_tp    <- contrasts_df$t_post
  meta_tr    <- contrasts_df$t_pre
  
  coh_row_map    <- setNames(seq_along(all_cohorts), as.character(all_cohorts))
  coh_ymeans_mat <- matrix(NA_real_, nrow = length(all_cohorts), ncol = TT)
  for (g in all_cohorts) {
    idx <- which(cohort_of == g)
    if (length(idx) > 0L)
      coh_ymeans_mat[coh_row_map[as.character(g)], ] <-
        colMeans(Y_mat[idx, , drop = FALSE], na.rm = TRUE)
  }
  
  tp_idx_v <- meta_tp - T_min + 1L
  tr_idx_v <- meta_tr - T_min + 1L
  g_row_v  <- as.integer(coh_row_map[as.character(meta_focal)])
  c_row_v  <- as.integer(coh_row_map[as.character(meta_ctrl)])
  
  Delta <- (coh_ymeans_mat[cbind(g_row_v, tp_idx_v)] -
              coh_ymeans_mat[cbind(g_row_v, tr_idx_v)]) -
    (coh_ymeans_mat[cbind(c_row_v, tp_idx_v)] -
       coh_ymeans_mat[cbind(c_row_v, tr_idx_v)])
  
  Q_H      <- build_Q_H(contrasts_df, N_2x2, N_beta)
  QtQ      <- crossprod(Q_H)
  QtQ_inv  <- tryCatch(solve(QtQ), error = function(e) ginv(QtQ))
  Q_H_rank <- qr(Q_H)$rank
  beta_hat <- as.numeric(QtQ_inv %*% crossprod(Q_H, Delta))
  
  n_grp         <- length(all_cohorts)
  units_by_grp  <- lapply(all_cohorts, function(g) which(cohort_of == g))
  active_by_grp <- lapply(all_cohorts,
                          function(g) which(meta_focal == g | meta_ctrl == g))
  sign_by_grp   <- lapply(seq_len(n_grp), function(gi) {
    idx <- active_by_grp[[gi]]
    ifelse(meta_focal[idx] == all_cohorts[gi], 1L, -1L)
  })
  
  T_g_count <- setNames(
    vapply(start_yrs, function(g) sum(cells$g == g), integer(1L)),
    as.character(start_yrs)
  )
  N_total_post <- sum(vapply(start_yrs, function(g)
    N_g_lookup[as.character(g)] * T_g_count[as.character(g)], numeric(1L)))
  
  w_EW <- numeric(N_beta)
  w_CW <- numeric(N_beta)
  for (ci in seq_len(N_beta)) {
    g_c      <- as.character(cells$g[ci])
    w_EW[ci] <- (1 / K) / T_g_count[g_c]
    w_CW[ci] <- N_g_lookup[g_c] / N_total_post
  }
  w_EW <- w_EW / sum(w_EW)
  w_CW <- w_CW / sum(w_CW)
  
  eff_ok           <- FALSE
  cohort_cov_out   <- NULL
  QtAQ_out         <- NULL
  n_iter           <- 0L
  
  # Fixed incidence factor of Omega_phi = U W U' (U is data-independent: it
  # depends only on the design). Its SVD basis and the projections u'Q_H, u'Delta
  # are precomputed ONCE; each GMM iteration then works only in the
  # m = n_grp*TT dimensional space, never forming an N_2x2-sized object.
  m_cols <- n_grp * TT
  U_fac  <- matrix(0.0, nrow = N_2x2, ncol = m_cols)
  for (gi in seq_len(n_grp)) {
    active <- active_by_grp[[gi]]
    if (length(active) == 0L) next
    sgn  <- sign_by_grp[[gi]]
    cols <- ((gi - 1L) * TT + 1L):(gi * TT)
    tp_i <- meta_tp[active] - T_min + 1L
    tr_i <- meta_tr[active] - T_min + 1L
    U_fac[cbind(active, cols[tp_i])] <- U_fac[cbind(active, cols[tp_i])] + sgn
    U_fac[cbind(active, cols[tr_i])] <- U_fac[cbind(active, cols[tr_i])] - sgn
  }
  red_fac <- gmm_reduced_factor(U_fac, Q_H, Delta)
  
  for (iter in seq_len(max_iter)) {
    n_iter   <- iter
    beta_old <- beta_hat
    
    idx <- match(paste(dt_long$g, dt_long$time),
                 paste(cells$g,   cells$t))
    dt_long[, tau_hat := ifelse(is.na(idx), 0, beta_hat[idx])]
    dt_long[, Y_adj   := Y - tau_hat]
    
    fe_mod <- feols(Y_adj ~ 1 | unit + time, data = dt_long)
    R_mat  <- matrix(NA_real_, nrow = N_units, ncol = TT)
    R_mat[cbind(dt_long$unit_row, dt_long$time_idx)] <- residuals(fe_mod)
    
    # Block-diagonal W (changes with beta through the residuals). Each cohort
    # uses an UNRESTRICTED within-unit covariance (full T x T): arbitrary serial
    # correlation is allowed, cross-sectional dependence is ignored.
    W               <- matrix(0.0, nrow = m_cols, ncol = m_cols)
    cohort_cov_list <- vector("list", n_grp)   # full T x T cohort covariance
    for (gi in seq_len(n_grp)) {
      g_val   <- all_cohorts[gi]
      N_g     <- N_g_lookup[as.character(g_val)]
      units_g <- units_by_grp[[gi]]
      active  <- active_by_grp[[gi]]
      cols    <- ((gi - 1L) * TT + 1L):(gi * TT)
      if (length(active) == 0L || length(units_g) == 0L) {
        cohort_cov_list[[gi]] <- matrix(NA_real_, TT, TT)
        next
      }
      rm_g  <- R_mat[units_g, , drop = FALSE]          # N_g x TT
      Sig_g <- crossprod(rm_g) / N_g                   # TT x TT unrestricted cov
      cohort_cov_list[[gi]] <- Sig_g
      W[cols, cols] <- Sig_g / N_g
    }
    cohort_cov_out <- cohort_cov_list
    
    # Optimal-GMM normal equations via the generalized inverse, in reduced space.
    red <- gmm_reduced_solve(red_fac, W, tol_rel = eig_tol)
    if (is.null(red)) break
    QtAQ     <- red$QtAQ
    QtAD     <- red$QtAD
    beta_new <- tryCatch(solve(QtAQ, QtAD), error = function(e) NULL)
    if (is.null(beta_new)) break
    
    QtAQ_out <- QtAQ
    eff_ok   <- TRUE
    beta_hat <- beta_new
    
    if (max(abs(beta_hat - beta_old)) < tol) break
  }
  
  if (eff_ok) {
    V_beta <- tryCatch(
      solve(QtAQ_out),
      error = function(e) ginv(as.matrix(QtAQ_out))
    )
  } else {
    V_beta <- ginv(crossprod(Q_H))   # degenerate fallback (not reached in practice)
  }
  
  SE_catt <- sqrt(pmax(0, diag(V_beta)))
  ATT_EW  <- as.numeric(crossprod(w_EW, beta_hat))
  ATT_CW  <- as.numeric(crossprod(w_CW, beta_hat))
  SE_EW   <- sqrt(max(0, as.numeric(crossprod(w_EW, V_beta %*% w_EW))))
  SE_CW   <- sqrt(max(0, as.numeric(crossprod(w_CW, V_beta %*% w_CW))))
  
  names(cohort_cov_out) <- as.character(all_cohorts)
  
  list(
    beta_hat    = beta_hat,
    SE_catt     = SE_catt,
    catt_out    = data.frame(g = cells$g, t = cells$t,
                             beta_hat = beta_hat, SE = SE_catt),
    ATT_EW      = ATT_EW,
    ATT_CW      = ATT_CW,
    SE_EW       = SE_EW,
    SE_CW       = SE_CW,
    w_EW        = w_EW,
    w_CW        = w_CW,
    cells       = cells,
    eff_ok      = eff_ok,
    n_iter      = n_iter,
    N_beta      = N_beta,
    N_2x2       = N_2x2,
    Q_H_rank    = Q_H_rank,
    cohort_cov  = cohort_cov_out
  )
}


# =====================================================================
# GMM I  (variance = "homo")  : Toeplitz with a single pooled set of
#                               autocovariances (homoskedastic).
# GMM II (variance = "hetero"): Toeplitz with cohort-specific
#                               autocovariances (heteroskedastic).
# Both use CLEAN comparisons only (NT + NYT) and assume cross-sectional
# independence -> the moment covariance is built from within-cohort
# blocks only (no cross-cohort / factor blocks, unlike GMM III).
# =====================================================================
gmm_staggered_toeplitz <- function(data, yname, tname, idname, gname,
                                   variance = c("homo", "hetero"),
                                   has_nt   = FALSE,
                                   eig_tol  = 1e-8,
                                   max_iter = 10L,
                                   tol      = 1e-6) {
  
  variance <- match.arg(variance)
  
  Y_raw  <- as.numeric(data[[yname]])
  t_raw  <- as.integer(data[[tname]])
  id_raw <- data[[idname]]
  g_raw  <- as.integer(data[[gname]])
  
  unit_ids  <- sort(unique(id_raw))
  N_units   <- length(unit_ids)
  T_min     <- min(t_raw)
  T_max     <- max(t_raw)
  TT        <- T_max - T_min + 1L
  unit_row  <- match(id_raw, unit_ids)
  time_idx  <- t_raw - T_min + 1L
  
  cohort_of   <- as.integer(tapply(g_raw, unit_row, function(x) x[1L]))
  start_yrs   <- sort(unique(cohort_of[cohort_of > 0L & cohort_of <= T_max]))
  K           <- length(start_yrs)
  all_cohorts <- if (has_nt) c(0L, start_yrs) else start_yrs
  
  N_g_lookup <- setNames(
    vapply(start_yrs, function(g) sum(cohort_of == g), integer(1L)),
    as.character(start_yrs)
  )
  if (has_nt)
    N_g_lookup["0"] <- sum(cohort_of == 0L)
  
  Y_mat <- matrix(NA_real_, nrow = N_units, ncol = TT)
  Y_mat[cbind(unit_row, time_idx)] <- Y_raw
  
  dt_long <- data.table(
    unit     = id_raw,
    unit_row = unit_row,
    time     = t_raw,
    time_idx = time_idx,
    Y        = Y_raw,
    g        = g_raw
  )
  setorder(dt_long, unit, time)
  
  cells        <- build_cells(start_yrs, T_max)
  N_beta       <- nrow(cells)
  contrasts_df <- enumerate_contrasts_clean(cells, all_cohorts, start_yrs,
                                            T_min, has_nt)
  N_2x2      <- nrow(contrasts_df)
  meta_focal <- contrasts_df$g
  meta_ctrl  <- contrasts_df$c
  meta_tp    <- contrasts_df$t_post
  meta_tr    <- contrasts_df$t_pre
  
  coh_row_map    <- setNames(seq_along(all_cohorts), as.character(all_cohorts))
  coh_ymeans_mat <- matrix(NA_real_, nrow = length(all_cohorts), ncol = TT)
  for (g in all_cohorts) {
    idx <- which(cohort_of == g)
    if (length(idx) > 0L)
      coh_ymeans_mat[coh_row_map[as.character(g)], ] <-
        colMeans(Y_mat[idx, , drop = FALSE], na.rm = TRUE)
  }
  
  tp_idx_v <- meta_tp - T_min + 1L
  tr_idx_v <- meta_tr - T_min + 1L
  g_row_v  <- as.integer(coh_row_map[as.character(meta_focal)])
  c_row_v  <- as.integer(coh_row_map[as.character(meta_ctrl)])
  
  Delta <- (coh_ymeans_mat[cbind(g_row_v, tp_idx_v)] -
              coh_ymeans_mat[cbind(g_row_v, tr_idx_v)]) -
    (coh_ymeans_mat[cbind(c_row_v, tp_idx_v)] -
       coh_ymeans_mat[cbind(c_row_v, tr_idx_v)])
  
  # With clean comparisons only, Q_H is a pure selection matrix (+1 on focal).
  Q_H      <- build_Q_H(contrasts_df, N_2x2, N_beta)
  QtQ      <- crossprod(Q_H)
  QtQ_inv  <- tryCatch(solve(QtQ), error = function(e) ginv(QtQ))
  Q_H_rank <- qr(Q_H)$rank
  beta_hat <- as.numeric(QtQ_inv %*% crossprod(Q_H, Delta))
  
  n_grp         <- length(all_cohorts)
  units_by_grp  <- lapply(all_cohorts, function(g) which(cohort_of == g))
  active_by_grp <- lapply(all_cohorts,
                          function(g) which(meta_focal == g | meta_ctrl == g))
  sign_by_grp   <- lapply(seq_len(n_grp), function(gi) {
    idx <- active_by_grp[[gi]]
    ifelse(meta_focal[idx] == all_cohorts[gi], 1L, -1L)
  })
  
  T_g_count <- setNames(
    vapply(start_yrs, function(g) sum(cells$g == g), integer(1L)),
    as.character(start_yrs)
  )
  N_total_post <- sum(vapply(start_yrs, function(g)
    N_g_lookup[as.character(g)] * T_g_count[as.character(g)], numeric(1L)))
  
  w_EW <- numeric(N_beta)
  w_CW <- numeric(N_beta)
  for (ci in seq_len(N_beta)) {
    g_c      <- as.character(cells$g[ci])
    w_EW[ci] <- (1 / K) / T_g_count[g_c]
    w_CW[ci] <- N_g_lookup[g_c] / N_total_post
  }
  w_EW <- w_EW / sum(w_EW)
  w_CW <- w_CW / sum(w_CW)
  
  eff_ok           <- FALSE
  sigma_d_list_out <- NULL
  QtAQ_out         <- NULL
  n_iter           <- 0L
  
  # Fixed incidence factor of Omega_phi = U W U' (U depends only on the design);
  # SVD basis and projections u'Q_H, u'Delta are precomputed ONCE. Each iteration
  # then works only in the m = n_grp*TT dimensional space.
  m_cols <- n_grp * TT
  U_fac  <- matrix(0.0, nrow = N_2x2, ncol = m_cols)
  for (gi in seq_len(n_grp)) {
    active <- active_by_grp[[gi]]
    if (length(active) == 0L) next
    sgn  <- sign_by_grp[[gi]]
    cols <- ((gi - 1L) * TT + 1L):(gi * TT)
    tp_i <- meta_tp[active] - T_min + 1L
    tr_i <- meta_tr[active] - T_min + 1L
    U_fac[cbind(active, cols[tp_i])] <- U_fac[cbind(active, cols[tp_i])] + sgn
    U_fac[cbind(active, cols[tr_i])] <- U_fac[cbind(active, cols[tr_i])] - sgn
  }
  red_fac <- gmm_reduced_factor(U_fac, Q_H, Delta)
  
  for (iter in seq_len(max_iter)) {
    n_iter   <- iter
    beta_old <- beta_hat
    
    idx <- match(paste(dt_long$g, dt_long$time),
                 paste(cells$g,   cells$t))
    dt_long[, tau_hat := ifelse(is.na(idx), 0, beta_hat[idx])]
    dt_long[, Y_adj   := Y - tau_hat]
    
    fe_mod <- feols(Y_adj ~ 1 | unit + time, data = dt_long)
    R_mat  <- matrix(NA_real_, nrow = N_units, ncol = TT)
    R_mat[cbind(dt_long$unit_row, dt_long$time_idx)] <- residuals(fe_mod)
    
    # Homoskedastic (GMM I): single pooled autocovariance sequence over all units.
    sigma_d_pool <- NULL
    if (variance == "homo") {
      sigma_d_pool <- numeric(TT)
      for (d in 0L:(TT - 1L)) {
        r1 <- seq_len(TT - d)
        r2 <- (1L + d):TT
        sigma_d_pool[d + 1L] <- sum(R_mat[, r1] * R_mat[, r2], na.rm = TRUE) /
          (N_units * (TT - d))
      }
    }
    
    # Block-diagonal W with W block = toeplitz(sigma_d_g)/N_g (homoskedastic ->
    # common sigma_d across cohorts; heteroskedastic -> cohort-specific).
    W            <- matrix(0.0, nrow = m_cols, ncol = m_cols)
    sigma_d_list <- vector("list", n_grp)
    for (gi in seq_len(n_grp)) {
      g_val   <- all_cohorts[gi]
      N_g     <- N_g_lookup[as.character(g_val)]
      units_g <- units_by_grp[[gi]]
      active  <- active_by_grp[[gi]]
      cols    <- ((gi - 1L) * TT + 1L):(gi * TT)
      if (length(active) == 0L || length(units_g) == 0L) {
        sigma_d_list[[gi]] <- rep(NA_real_, TT)
        next
      }
      if (variance == "hetero") {
        rm_g      <- R_mat[units_g, , drop = FALSE]
        sigma_d_g <- numeric(TT)
        for (d in 0L:(TT - 1L)) {
          r1 <- seq_len(TT - d)
          r2 <- (1L + d):TT
          sigma_d_g[d + 1L] <- sum(rm_g[, r1] * rm_g[, r2], na.rm = TRUE) /
            (N_g * (TT - d))
        }
      } else {
        sigma_d_g <- sigma_d_pool          # homoskedastic: common across cohorts
      }
      sigma_d_list[[gi]] <- sigma_d_g
      W[cols, cols] <- toeplitz(sigma_d_g) / N_g
    }
    sigma_d_list_out <- sigma_d_list
    
    # Optimal-GMM normal equations via the generalized inverse, in reduced space.
    red <- gmm_reduced_solve(red_fac, W, tol_rel = eig_tol)
    if (is.null(red)) break
    QtAQ     <- red$QtAQ
    QtAD     <- red$QtAD
    beta_new <- tryCatch(solve(QtAQ, QtAD), error = function(e) NULL)
    if (is.null(beta_new)) break
    
    QtAQ_out <- QtAQ
    eff_ok   <- TRUE
    beta_hat <- beta_new
    
    if (max(abs(beta_hat - beta_old)) < tol) break
  }
  
  if (eff_ok) {
    V_beta <- tryCatch(
      solve(QtAQ_out),
      error = function(e) ginv(as.matrix(QtAQ_out))
    )
  } else {
    V_beta <- ginv(crossprod(Q_H))   # degenerate fallback (not reached in practice)
  }
  
  SE_catt <- sqrt(pmax(0, diag(V_beta)))
  ATT_EW  <- as.numeric(crossprod(w_EW, beta_hat))
  ATT_CW  <- as.numeric(crossprod(w_CW, beta_hat))
  SE_EW   <- sqrt(max(0, as.numeric(crossprod(w_EW, V_beta %*% w_EW))))
  SE_CW   <- sqrt(max(0, as.numeric(crossprod(w_CW, V_beta %*% w_CW))))
  
  names(sigma_d_list_out) <- as.character(all_cohorts)
  
  list(
    beta_hat    = beta_hat,
    SE_catt     = SE_catt,
    catt_out    = data.frame(g = cells$g, t = cells$t,
                             beta_hat = beta_hat, SE = SE_catt),
    ATT_EW      = ATT_EW,
    ATT_CW      = ATT_CW,
    SE_EW       = SE_EW,
    SE_CW       = SE_CW,
    w_EW        = w_EW,
    w_CW        = w_CW,
    cells       = cells,
    eff_ok      = eff_ok,
    n_iter      = n_iter,
    N_beta      = N_beta,
    N_2x2       = N_2x2,
    Q_H_rank    = Q_H_rank,
    variance    = variance,
    sigma_d     = sigma_d_list_out
  )
}


estimate_twfe <- function(dt) {
  tryCatch(
    coef(feols(Y ~ D | unit + time, data = dt))["D"],
    error = function(e) NA_real_
  )
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
    dt_g <- copy(dt)
    dt_g[, first_treat := fifelse(g == 0, Inf, as.numeric(g))]
    coef(did2s(data = as.data.frame(dt_g), yname = "Y",
               first_stage = ~ 0 | unit + time, second_stage = ~ i(D, ref = 0),
               treatment = "D", cluster_var = "unit", verbose = FALSE))["D::1"]
  }, error = function(e) NA_real_)
}

estimate_flex_twfe <- function(dt) {
  tryCatch({
    dt_f <- copy(dt)
    dt_f[, treat_gt := fifelse(D == 1, g * 100L + as.integer(time), 0L)]
    mod <- feols(Y ~ i(treat_gt, ref = 0) | unit + time, data = dt_f)
    mean(coef(mod), na.rm = TRUE)
  }, error = function(e) NA_real_)
}


# =====================================================================
# GLS Flex TWFE I / II / III  (feasible, iterated)
# Flexible TWFE regression (unit FE + time FE + cohort x period treatment
# interactions) estimated by feasible GLS with a BLOCK-DIAGONAL (across
# units) within-unit error covariance whose structure mirrors GMM I/II/III:
#   I   "toeplitz_homo"   : common Toeplitz Sigma (pooled autocovariances)
#   II  "toeplitz_hetero" : cohort-specific Toeplitz Sigma^(g)
#   III "full_hetero"     : cohort-specific UNRESTRICTED Sigma^(g) (full
#                           T x T sample covariance; drops stationarity)
# Sigma is estimated from the regression residuals and the fit is iterated
# to convergence. The design matrix (unit/time/treatment structure) is
# fixed across simulations, so build_flex_design() is called once and only
# the outcome vector is swapped each replication.
# =====================================================================

build_flex_design <- function(dt) {
  dd <- copy(dt)
  setorder(dd, unit, time)
  dd[, treat_gt := fifelse(D == 1L, g * 100L + as.integer(time), 0L)]
  dd[, unit_f := factor(unit)]
  dd[, time_f := factor(time)]
  dd[, tg_f   := relevel(factor(treat_gt), ref = "0")]
  
  Z  <- model.matrix(~ unit_f + time_f + tg_f, data = dd)
  cn <- colnames(Z)
  
  tg_cols <- grep("^tg_f", cn)
  tg_keys <- as.integer(sub("^tg_f", "", cn[tg_cols]))
  catt_g  <- tg_keys %/% 100L
  catt_t  <- tg_keys %%  100L
  
  N_units <- length(unique(dd$unit))
  TT      <- length(unique(dd$time))
  
  row_cohort       <- as.integer(dd$g)
  unit_cohort_sort <- row_cohort[seq(1L, nrow(dd), by = TT)]  # cohort per unit (R_mat row order)
  groups           <- sort(unique(row_cohort))
  rows_by_grp      <- lapply(groups, function(g) which(row_cohort == g))
  units_by_grp     <- lapply(groups, function(g) which(unit_cohort_sort == g))
  Ng_by_grp        <- vapply(units_by_grp, length, integer(1L))
  
  # aggregation weights over CATT cells (match GMM EW / CW conventions)
  start_yrs        <- sort(unique(catt_g))
  K                <- length(start_yrs)
  T_g_count        <- vapply(start_yrs, function(g) sum(catt_g == g), integer(1L))
  names(T_g_count) <- as.character(start_yrs)
  Ng_treated       <- setNames(
    vapply(start_yrs, function(g) Ng_by_grp[match(g, groups)], integer(1L)),
    as.character(start_yrs))
  N_total_post     <- sum(Ng_treated[as.character(start_yrs)] *
                            T_g_count[as.character(start_yrs)])
  
  w_EW <- numeric(length(tg_cols)); w_CW <- numeric(length(tg_cols))
  for (j in seq_along(tg_cols)) {
    gc      <- as.character(catt_g[j])
    w_EW[j] <- (1 / K) / T_g_count[gc]
    w_CW[j] <- Ng_treated[gc] / N_total_post
  }
  w_EW <- w_EW / sum(w_EW)
  w_CW <- w_CW / sum(w_CW)
  
  XtX     <- crossprod(Z)
  XtX_inv <- tryCatch(solve(XtX), error = function(e) ginv(XtX))
  
  list(Z = Z, tg_cols = tg_cols, catt_g = catt_g, catt_t = catt_t,
       N_units = N_units, TT = TT, groups = groups,
       rows_by_grp = rows_by_grp, units_by_grp = units_by_grp,
       Ng_by_grp = Ng_by_grp, w_EW = w_EW, w_CW = w_CW,
       XtX_inv = XtX_inv)
}


.chol_pd <- function(S) {
  S  <- (S + t(S)) / 2
  ok <- tryCatch({ chol(S); TRUE }, error = function(e) FALSE)
  if (!ok) {
    ev <- min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
    diag(S) <- diag(S) + max(0, -ev) + 1e-8
  } else {
    diag(S) <- diag(S) + 1e-8
  }
  chol(S)
}


flex_gls_fit <- function(design, Yvec,
                         cov_type = c("toeplitz_homo", "toeplitz_hetero", "full_hetero"),
                         max_iter = 10L, tol = 1e-6) {
  cov_type <- match.arg(cov_type)
  Z   <- design$Z
  TT  <- design$TT
  Nu  <- design$N_units
  grp <- design$groups
  
  beta <- as.numeric(design$XtX_inv %*% crossprod(Z, Yvec))   # initial OLS
  iter <- 0L
  
  for (iter in seq_len(max_iter)) {
    beta_old <- beta
    r        <- as.numeric(Yvec - Z %*% beta)
    R_mat    <- matrix(r, nrow = Nu, ncol = TT, byrow = TRUE)  # unit x time
    
    sigma_d_pool <- NULL
    if (cov_type == "toeplitz_homo") {
      sigma_d_pool <- numeric(TT)
      for (d in 0L:(TT - 1L))
        sigma_d_pool[d + 1L] <- sum(R_mat[, seq_len(TT - d)] *
                                      R_mat[, (1L + d):TT]) / (Nu * (TT - d))
    }
    
    Yt <- numeric(length(Yvec))
    Zt <- matrix(0.0, nrow = nrow(Z), ncol = ncol(Z))
    
    for (gi in seq_along(grp)) {
      rows_g  <- design$rows_by_grp[[gi]]
      units_g <- design$units_by_grp[[gi]]
      N_g     <- design$Ng_by_grp[gi]
      if (length(rows_g) == 0L) next
      
      if (cov_type == "toeplitz_homo") {
        Sig_g <- toeplitz(sigma_d_pool)
      } else if (cov_type == "toeplitz_hetero") {
        Rg <- R_mat[units_g, , drop = FALSE]
        sd <- numeric(TT)
        for (d in 0L:(TT - 1L))
          sd[d + 1L] <- sum(Rg[, seq_len(TT - d)] * Rg[, (1L + d):TT]) /
          (N_g * (TT - d))
        Sig_g <- toeplitz(sd)
      } else {                                   # full_hetero
        Rg    <- R_mat[units_g, , drop = FALSE]
        Sig_g <- crossprod(Rg) / N_g             # T x T unrestricted
      }
      
      U   <- .chol_pd(Sig_g)
      P_g <- solve(t(U))                         # P_g' P_g = Sigma_g^{-1}
      
      ## transform outcome for this cohort group
      Ymat       <- matrix(Yvec[rows_g], nrow = TT)           # T x N_g
      Yt[rows_g] <- as.numeric(P_g %*% Ymat)
      
      ## transform design (block-diagonal: same P_g per unit)
      A            <- matrix(Z[rows_g, ], nrow = TT)          # T x (N_g * p)
      Zt[rows_g, ] <- matrix(P_g %*% A, nrow = length(rows_g))
    }
    
    XtX  <- crossprod(Zt)
    XtY  <- crossprod(Zt, Yt)
    beta <- tryCatch(as.numeric(solve(XtX, XtY)),
                     error = function(e) as.numeric(ginv(XtX) %*% XtY))
    
    if (max(abs(beta - beta_old)) < tol) break
  }
  
  catt   <- beta[design$tg_cols]
  ATT_EW <- sum(design$w_EW * catt)
  ATT_CW <- sum(design$w_CW * catt)
  list(ATT_EW = ATT_EW, ATT_CW = ATT_CW, catt = catt, n_iter = iter)
}


generate_data <- function(beta_g_vec, r_g_vec) {
  unit_id <- rep(seq_len(N_total), each = T_total)
  time_id <- rep(seq_len(T_total), times = N_total)
  alpha   <- rnorm(N_total)
  lambda  <- rnorm(T_total)
  eps     <- as.vector(t(mvrnorm(N_total, rep(0, T_total), Sigma_ar1)))
  g_vec   <- unit_cohort[unit_id]
  D_vec   <- as.integer(g_vec > 0 & time_id >= g_vec)
  tau_vec <- numeric(N_total * T_total)
  for (c_idx in seq_len(n_cohorts)) {
    g_c  <- treatment_times[c_idx]
    mask <- (g_vec == g_c) & (D_vec == 1)
    tau_vec[mask] <- beta_g_vec[c_idx] * (1 + r_g_vec[c_idx])^(time_id[mask] - g_c)
  }
  dt <- data.table(unit = unit_id, time = time_id,
                   Y = alpha[unit_id] + lambda[time_id] + tau_vec + eps,
                   D = D_vec, g = g_vec)
  dt[, g_inf := fifelse(g == 0, Inf, as.numeric(g))]
  dt[, g_cs  := fifelse(g == 0, 0, as.numeric(g))]
  return(dt)
}

compute_true_att <- function(beta_g_vec, r_g_vec) {
  total_te <- 0; total_obs <- 0
  for (ci in seq_len(n_cohorts)) {
    g_c <- treatment_times[ci]
    for (t in g_c:T_total) {
      total_te  <- total_te  + cohort_size * beta_g_vec[ci] * (1 + r_g_vec[ci])^(t - g_c)
      total_obs <- total_obs + cohort_size
    }
  }
  total_te / total_obs
}


run_simulation <- function(beta_g_vec, r_g_vec, label) {
  true_att <- compute_true_att(beta_g_vec, r_g_vec)
  cat(sprintf("\n=== %s (N=%d) ===\nTrue ATT: %.4f\n", label, N_total, true_att))
  
  est_names <- c("TWFE", "CS", "SA", "Gardner", "Flex_TWFE",
                 "GMM_I", "GMM_II", "GMM_III",
                 "GLS_I", "GLS_II", "GLS_III")
  n_est   <- length(est_names)
  
  # Flex-TWFE design is fixed across sims (deterministic unit/time/treatment
  # structure); build it once here, then only swap the outcome each sim.
  struct_dt <- data.table(
    unit = rep(seq_len(N_total), each = T_total),
    time = rep(seq_len(T_total), times = N_total)
  )
  struct_dt[, g := unit_cohort[unit]]
  struct_dt[, D := as.integer(g > 0 & time >= g)]
  flex_design <- tryCatch(build_flex_design(struct_dt), error = function(e) NULL)
  
  # Generate every dataset SEQUENTIALLY first, so the RNG stream---and hence the
  # full set of simulated datasets---is byte-for-byte identical to the serial
  # implementation (no estimator consumes RNG). Estimation is deterministic, so
  # parallelizing it below does not change any reported value.
  data_list <- lapply(seq_len(n_sims),
                      function(s) generate_data(beta_g_vec, r_g_vec))
  
  # Per-simulation estimation: returns one results row. Contains no RNG.
  estimate_one <- function(dt) {
    att_twfe    <- estimate_twfe(dt)
    att_cs      <- estimate_cs(dt)
    att_sa      <- estimate_sa(dt)
    att_gardner <- estimate_gardner(dt)
    att_flex    <- estimate_flex_twfe(dt)
    
    g1_res <- tryCatch(gmm_staggered_toeplitz(dt, "Y", "time", "unit", "g",
                                              variance = "homo",   has_nt = TRUE), error = function(e) NULL)
    g2_res <- tryCatch(gmm_staggered_toeplitz(dt, "Y", "time", "unit", "g",
                                              variance = "hetero", has_nt = TRUE), error = function(e) NULL)
    hc_res <- tryCatch(gmm_staggered_hetcov(dt, "Y", "time", "unit", "g",
                                            has_nt = TRUE), error = function(e) NULL)
    
    yv   <- dt$Y
    gls1 <- tryCatch(if (!is.null(flex_design)) flex_gls_fit(flex_design, yv, "toeplitz_homo")   else NULL, error = function(e) NULL)
    gls2 <- tryCatch(if (!is.null(flex_design)) flex_gls_fit(flex_design, yv, "toeplitz_hetero") else NULL, error = function(e) NULL)
    gls3 <- tryCatch(if (!is.null(flex_design)) flex_gls_fit(flex_design, yv, "full_hetero")     else NULL, error = function(e) NULL)
    
    data.table(
      TWFE = att_twfe, CS = att_cs, SA = att_sa,
      Gardner = att_gardner, Flex_TWFE = att_flex,
      GMM_I   = if (!is.null(g1_res)) g1_res$ATT_CW else NA_real_,
      GMM_II  = if (!is.null(g2_res)) g2_res$ATT_CW else NA_real_,
      GMM_III = if (!is.null(hc_res)) hc_res$ATT_CW else NA_real_,
      GLS_I   = if (!is.null(gls1)) gls1$ATT_CW else NA_real_,
      GLS_II  = if (!is.null(gls2)) gls2$ATT_CW else NA_real_,
      GLS_III = if (!is.null(gls3)) gls3$ATT_CW else NA_real_
    )
  }
  
  t0 <- proc.time()[3]
  cat(sprintf("Running %d sims on %d core(s)...\n", n_sims, n_cores)); flush(con)
  rows <- if (n_cores > 1L) {
    mclapply(data_list, estimate_one, mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    lapply(data_list, estimate_one)
  }
  cat(sprintf("Estimation finished in %.1f min\n", (proc.time()[3] - t0) / 60)); flush(con)
  
  # Guard against a worker dying (returns a non-data.table): replace with NA row.
  na_row <- data.table(TWFE = NA_real_, CS = NA_real_, SA = NA_real_,
                       Gardner = NA_real_, Flex_TWFE = NA_real_,
                       GMM_I = NA_real_, GMM_II = NA_real_, GMM_III = NA_real_,
                       GLS_I = NA_real_, GLS_II = NA_real_, GLS_III = NA_real_)
  rows <- lapply(rows, function(r) if (is.data.table(r)) r else na_row)
  
  results <- rbindlist(rows)
  results[, sim := seq_len(n_sims)]
  setcolorder(results, c("sim", est_names))
  
  summary_dt <- data.table(Estimator = est_names,
                           Bias = numeric(n_est), Variance = numeric(n_est))
  for (i in seq_len(n_est)) {
    vals <- results[[est_names[i]]]; vals <- vals[!is.na(vals)]
    summary_dt$Bias[i]     <- round(mean(vals) - true_att, 4)
    summary_dt$Variance[i] <- round(var(vals), 4)
  }
  
  cat(sprintf("\nResults (%s):\n", label))
  print(summary_dt)
  
  list(results = results, summary = summary_dt,
       true_att = true_att)
}



beta_het <- c(-16, -12, -10, -9, -2)
r_het    <- c(0.01, 0.04, 0.08, 0.10, 0.07)

## Sweep over AR(1) persistence. For each rho we rebuild the within-unit
## covariance Sigma_ar1 (used by generate_data) and reset the seed so the base
## random draws are identical across rho -- only the serial correlation changes.
res_by_rho        <- vector("list", length(rho_grid))
names(res_by_rho) <- sprintf("rho_%.1f", rho_grid)

for (ri in seq_along(rho_grid)) {
  rho       <- rho_grid[ri]
  Sigma_ar1 <- toeplitz(rho^(0:(T_total - 1L)))   # global; used by generate_data
  set.seed(312844)                                # same base draws across rho
  cat(sprintf("\n################  rho = %.1f  ################\n", rho)); flush(con)
  res_by_rho[[ri]] <- run_simulation(beta_het, r_het,
                                     sprintf("Heterogeneous (rho=%.1f)", rho))
}

## ---- Combined summary across rho: estimators (rows) x rho (cols) ----
est_names <- res_by_rho[[1L]]$summary$Estimator
rho_lab   <- sprintf("rho=%.1f", rho_grid)
bias_tab  <- data.table(Estimator = est_names)
var_tab   <- data.table(Estimator = est_names)
for (ri in seq_along(rho_grid)) {
  s <- res_by_rho[[ri]]$summary
  idx <- match(est_names, s$Estimator)
  bias_tab[[rho_lab[ri]]] <- round(s$Bias[idx], 4)
  var_tab[[ rho_lab[ri]]] <- round(s$Variance[idx], 4)
}

cat("\n==================== SUMMARY ACROSS rho ====================\n")
cat(sprintf("True ATT (constant across rho): %.4f\n", res_by_rho[[1L]]$true_att))
cat("\nBias  (estimator x rho):\n");     print(bias_tab, row.names = FALSE)
cat("\nVariance (estimator x rho):\n");  print(var_tab,  row.names = FALSE)

save(res_by_rho, bias_tab, var_tab, rho_grid,
     file = "/home/rishabh.bijani_ug2023/gmm_sims_01_06/hetcov_dgp2_sweep_het_01_06.RData")

sink(type = "message")
sink()
close(con)

N=50
Variance (estimator x rho):
  Estimator rho=0.0 rho=0.1 rho=0.3 rho=0.5 rho=0.7 rho=0.9
<char>   <num>   <num>   <num>   <num>   <num>   <num>
  TWFE  0.0015  0.0016  0.0026  0.0034  0.0048  0.0060
CS  0.0103  0.0082  0.0087  0.0100  0.0114  0.0114
SA  0.0103  0.0082  0.0087  0.0100  0.0114  0.0114
Gardner  0.0022  0.0023  0.0037  0.0059  0.0087  0.0126
Flex_TWFE  0.0022  0.0023  0.0037  0.0059  0.0087  0.0126
GMM_I  0.0022  0.0023  0.0037  0.0058  0.0078  0.0096
GMM_II  0.0022  0.0023  0.0038  0.0058  0.0080  0.0100
GMM_III  0.0030  0.0031  0.0049  0.0079  0.0114  0.0125
GLS_I  0.0022  0.0023  0.0037  0.0058  0.0078  0.0094
GLS_II  0.0022  0.0024  0.0038  0.0058  0.0081  0.0097
GLS_III  0.0030  0.0031  0.0049  0.0080  0.0115  0.0126


N= 10
Bias  (estimator x rho):
  Estimator rho=0.0 rho=0.1 rho=0.3 rho=0.5 rho=0.7 rho=0.9
<char>   <num>   <num>   <num>   <num>   <num>   <num>
  TWFE  4.0010  4.0077  4.0022  4.0063  4.0065  3.9971
CS -0.0078  0.0029  0.0085 -0.0008  0.0143 -0.0033
SA -0.0078  0.0029  0.0085 -0.0008  0.0143 -0.0033
Gardner -0.0053  0.0053 -0.0025 -0.0043  0.0072 -0.0103
Flex_TWFE -0.0053  0.0053 -0.0025 -0.0043  0.0072 -0.0103
GMM_I -0.0047  0.0046 -0.0019 -0.0054  0.0044 -0.0034
GMM_II -0.0030  0.0047 -0.0018 -0.0060  0.0084 -0.0045
GMM_III -0.0057  0.0052 -0.0045 -0.0042  0.0066 -0.0125
GLS_I -0.0048  0.0045 -0.0019 -0.0056  0.0052 -0.0037
GLS_II  2.1598  1.5177  2.6906  2.1941  2.4454  1.1958
GLS_III -0.0007  0.0024  0.0010 -0.0037  0.0045 -0.0105

Variance (estimator x rho):
  Estimator rho=0.0 rho=0.1 rho=0.3 rho=0.5 rho=0.7 rho=0.9
<char>   <num>   <num>   <num>   <num>   <num>   <num>
  TWFE  0.0069  0.0089  0.0117  0.0166  0.0230  0.0294
CS  0.0401  0.0442  0.0422  0.0480  0.0582  0.0610
SA  0.0401  0.0442  0.0422  0.0480  0.0582  0.0610
Gardner  0.0099  0.0129  0.0187  0.0289  0.0442  0.0648
Flex_TWFE  0.0099  0.0129  0.0187  0.0289  0.0442  0.0648
GMM_I  0.0103  0.0129  0.0185  0.0281  0.0395  0.0489
GMM_II  0.0108  0.0139  0.0191  0.0296  0.0416  0.0514
GMM_III  0.0103  0.0135  0.0204  0.0312  0.0488  0.0742
GLS_I  0.0103  0.0129  0.0185  0.0282  0.0394  0.0481
GLS_II 30.3509 22.4447 36.8005 30.9968 33.8037 18.2023
GLS_III  0.0131  0.0161  0.0234  0.0358  0.0536  0.0782
