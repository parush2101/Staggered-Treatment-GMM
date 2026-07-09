# =============================================================================
# Does the GMM J-test detect violations of parallel trends / no anticipation?
# Derived from hetcov_dgp2_arbitrary_serial.R (same DGP): T=33, cohorts at
# {10,13,16,19,22}, never-treated pool, tau = beta_g (1+r_g)^{t-g}, errors with
# ARBITRARY (mixture-of-AR-kernels) serial correlation.
#
# Estimator: event-study basis (never-treated control, base period t0), r=G(T-1)
# minimal-set moments, iterated homoskedastic-Toeplitz optimal weight. The
# over-identifying restrictions are the pre-treatment placebo moments (t<g);
#   J = phi' Omega^{-1} phi  ~ chi^2_{df},  df = r - N_beta,
# is a joint test of parallel trends + no anticipation.
#
# Experiments:
#   (1) SIZE  under H0 (parallel trends, no anticipation): reject rate ~ nominal?
#   (2) POWER vs differential pre-trends (PT violation).
#   (3) POWER vs anticipation (effect at t=g-1, g-2).
# Base R only.  QUICK=1 for a fast check.
# =============================================================================

set.seed(312844)

T_total   <- 33L
cohorts   <- c(10L,13L,16L,19L,22L); G <- length(cohorts)
beta_het  <- c(-16,-12,-10,-9,-2)
r_het     <- c(0.01,0.04,0.08,0.10,0.07)
t0        <- 7L                      # base pre-period (< min cohort-2, avoids anticipation cells)
n_sims    <- { v <- suppressWarnings(as.integer(Sys.getenv("N_SIMS")))
  if (!is.na(v) && v>=1L) v else if (nzchar(Sys.getenv("QUICK"))) 200L else 1000L }

build_arbitrary_serial_cov <- function(TT, strength=c("low","medium","high")) {
  strength <- match.arg(strength)
  spec <- switch(strength,
    low    = list(rho=c(0.40,0.10,-0.20), w=c(0.45,0.40,0.15)),
    medium = list(rho=c(0.65,0.35,-0.20), w=c(0.55,0.30,0.15)),
    high   = list(rho=c(0.92,0.70,-0.25), w=c(0.60,0.30,0.10)))
  lags <- 0:(TT-1L); S <- matrix(0,TT,TT)
  for (m in seq_along(spec$rho)) S <- S + spec$w[m]*toeplitz(spec$rho[m]^lags)
  d <- sqrt(diag(S)); S <- S/outer(d,d); (S+t(S))/2
}

# ---- moment basis (fixed structure) ----------------------------------------
# PRE_WINDOW=L keeps only placebo (pre-treatment) moments within L periods before
# each g, shrinking the over-identifying df; default Inf keeps all (df=70).
pre_window <- { v <- suppressWarnings(as.numeric(Sys.getenv("PRE_WINDOW")))
  if (!is.na(v) && v>=1) v else Inf }
MOM <- list()
for (g in cohorts) for (t in setdiff(1:T_total, t0)) {
  is_catt <- (t >= g)
  if (is_catt || ((g - t) <= pre_window))
    MOM[[length(MOM)+1L]] <- list(g=g, t=t, catt=is_catt)
}
r        <- length(MOM)                       # G(T-1) = 160
catt_idx <- which(vapply(MOM, `[[`, logical(1), "catt"))
Nbeta    <- length(catt_idx)                  # 90
df       <- r - Nbeta                          # 70
Q <- matrix(0, r, Nbeta); for (j in seq_along(catt_idx)) Q[catt_idx[j], j] <- 1

# ---- one estimator+J evaluation on outcome matrix Y (N x T) ----------------
# (cohort structure passed via globals set in configure())
delta_vec <- function(Y) {
  cm <- matrix(0, G, T_total)
  for (ci in 1:G) cm[ci,] <- colMeans(Y[coh_units[[ci]], , drop=FALSE])
  nm <- colMeans(Y[nev_units, , drop=FALSE])
  d <- numeric(r)
  for (s in 1:r) { g<-MOM[[s]]$g; t<-MOM[[s]]$t; ci<-match(g,cohorts)
    d[s] <- (cm[ci,t]-cm[ci,t0]) - (nm[t]-nm[t0]) }
  d
}
omega <- function(Sig) {                        # Omega = sum_groups N_k W_k Sig W_k'
  Om <- matrix(0, r, r)
  for (k in 1:(G+1L)) { Wk <- matrix(Wc[k,,], r, T_total); Om <- Om + Ngrp[k]*(Wk %*% Sig %*% t(Wk)) }
  (Om+t(Om))/2
}
toeplitz_resid <- function(Y, bhat) {
  th <- matrix(0, N, T_total)
  for (j in seq_along(catt_idx)) { g<-MOM[[catt_idx[j]]]$g; t<-MOM[[catt_idx[j]]]$t
    th[unit_cohort==g, t] <- bhat[j] }
  Yt <- Y - th
  nu <- Yt - rowMeans(Yt) - rep(colMeans(Yt), each=N) + mean(Yt)
  s <- numeric(T_total)
  for (d0 in 0:(T_total-1L)) { a<-nu[,1:(T_total-d0),drop=FALSE]; b<-nu[,(1+d0):T_total,drop=FALSE]; s[d0+1L]<-mean(a*b) }
  toeplitz(s)
}
est_J <- function(Y, iters=3L) {
  d <- delta_vec(Y); bhat <- d[catt_idx]        # A=I start
  Oi <- NULL
  for (j in 1:iters) {
    Om <- omega(toeplitz_resid(Y, bhat)); Oi <- solve(Om)
    QtOi <- crossprod(Q, Oi)                    # t(Q) %*% Oi
    bhat <- as.numeric(solve(QtOi %*% Q, QtOi %*% d))
  }
  phi <- d - Q %*% bhat
  Jstat <- as.numeric(t(phi) %*% Oi %*% phi)
  pchisq(Jstat, df=df, lower.tail=FALSE)         # return p-value
}

# ---- DGP -------------------------------------------------------------------
simulate_Y <- function(viol="none", mag=0) {
  alpha <- rnorm(N); lambda <- rnorm(T_total)
  eps <- matrix(rnorm(N*T_total), N, T_total) %*% L_err
  tau <- matrix(0, N, T_total)
  for (ci in 1:G) { g<-cohorts[ci]; for (t in g:T_total) tau[coh_units[[ci]], t] <- beta_het[ci]*(1+r_het[ci])^(t-g) }
  Y <- matrix(alpha,N,T_total) + matrix(lambda,N,T_total,byrow=TRUE) + tau + eps
  if (viol=="pretrend") {                        # differential linear trend across cohorts (never = 0)
    tt <- (1:T_total) - mean(1:T_total)
    for (ci in 1:G) Y[coh_units[[ci]],] <- Y[coh_units[[ci]],] + matrix(mag*ci*tt, length(coh_units[[ci]]), T_total, byrow=TRUE)
  }
  if (viol=="anticip") {                          # anticipation at t=g-1 and g-2
    for (ci in 1:G) { g<-cohorts[ci]
      Y[coh_units[[ci]], g-1L] <- Y[coh_units[[ci]], g-1L] + mag*beta_het[ci]
      Y[coh_units[[ci]], g-2L] <- Y[coh_units[[ci]], g-2L] + mag*beta_het[ci] }
  }
  Y
}

# ---- configuration (sample size, serial strength) --------------------------
configure <- function(cohort_size, strength) {
  cs <<- cohort_size; nn <<- cohort_size
  unit_cohort <<- c(rep(cohorts, each=cohort_size), rep(0L, nn))
  N <<- length(unit_cohort)
  coh_units <<- lapply(cohorts, function(g) which(unit_cohort==g))
  nev_units <<- which(unit_cohort==0L)
  Ngrp <<- c(rep(cohort_size, G), nn)
  # group weight tensors W_k (r x T): cohorts then never
  Wc <<- array(0, c(G+1L, r, T_total))
  for (ci in 1:G) for (s in 1:r) if (MOM[[s]]$g==cohorts[ci]) {
    Wc[ci, s, MOM[[s]]$t] <<- 1/cohort_size; Wc[ci, s, t0] <<- -1/cohort_size }
  for (s in 1:r) { Wc[G+1L, s, MOM[[s]]$t] <<- -1/nn; Wc[G+1L, s, t0] <<- 1/nn }
  L_err <<- chol(build_arbitrary_serial_cov(T_total, strength))
}
reject_rate <- function(viol="none", mag=0, nsim=n_sims, levels=c(0.05,0.10)) {
  p <- replicate(nsim, est_J(simulate_Y(viol, mag)))
  sapply(levels, function(a) mean(p < a))
}

# ---- run -------------------------------------------------------------------
if (sys.nframe()==0L) {
  cat(sprintf("Moments r=%d  parameters N_beta=%d  df=%d  (chi^2 mean=%d)  n_sims=%d\n\n", r, Nbeta, df, df, n_sims))

  cat("=== (1) SIZE under H0 (parallel trends, no anticipation): rejection rate ===\n")
  cat(sprintf("%-8s %-8s  %8s  %8s\n","N_g","serial","@5%","@10%"))
  for (cz in c(50L,10L)) for (st in c("low","medium","high")) {
    configure(cz, st); rr <- reject_rate("none",0)
    cat(sprintf("%-8d %-8s  %7.1f%%  %7.1f%%\n", cz, st, 100*rr[1], 100*rr[2]))
  }

  cat("\n=== (2) POWER vs differential PRE-TRENDS  (N_g=50, serial=medium) ===\n")
  configure(50L,"medium")
  cat(sprintf("%-12s  %8s\n","slope mag","reject@5%"))
  for (m in c(0,0.02,0.05,0.10,0.20)) cat(sprintf("%-12.3f  %7.1f%%\n", m, 100*reject_rate("pretrend",m)[1]))

  cat("\n=== (3) POWER vs ANTICIPATION  (fraction of impact effect leaked to t=g-1,g-2; N_g=50, serial=medium) ===\n")
  configure(50L,"medium")
  cat(sprintf("%-12s  %8s\n","anticip frac","reject@5%"))
  for (a in c(0,0.05,0.10,0.20,0.40)) cat(sprintf("%-12.3f  %7.1f%%\n", a, 100*reject_rate("anticip",a)[1]))
}
