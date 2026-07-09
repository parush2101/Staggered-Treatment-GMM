# =============================================================================
# Switchers (reversible binary treatment) under ARBITRARY serial correlation.
# OVER-IDENTIFIED switch moments: each switch effect delta is identified by
# several pre-window contrasts (t-1,...,t-K), all valid under no-carryover, so
# the weighting matrix has something to exploit. Four estimators:
#
#   DID_M   : de Chaisemartin & D'Haultfoeuille (2020). Adjacent (k=1) contrast
#             per switch cell, switcher-count (frequency) weighted.  [benchmark]
#   GMM_eq  : over-identified moments, EQUAL weight (A = I).         [+ averaging]
#   GMM_sph : over-identified moments, weight = inverse moment cov implied by
#             SPHERICAL errors (Sigma = I_T).  Switcher analog of Flex TWFE/Gardner.
#   GMM_T   : over-identified moments, weight = inverse moment cov under a
#             homoskedastic-Toeplitz error model.  [+ serial correlation]  <- ours
#
# Multiplicative decomposition of the variance gain, reported each run:
#   var(DID_M)/var(GMM_eq)  = over-identification (averaging in more moments)
#   var(GMM_eq)/var(GMM_sph)= structure present even under spherical errors
#   var(GMM_sph)/var(GMM_T) = PURE serial-correlation channel (the contribution)
#
# No carryover: Y_it responds to the current treatment only. Base R only.
# Full run n_sims=500; quick check: QUICK=1 (or N_SIMS=...).
# =============================================================================

set.seed(312844)

# ---- config ----------------------------------------------------------------
T_total <- 33L
delta   <- -5.0            # true constant instantaneous effect => ATT = delta
sigma   <- 1.0
Kcap    <- 6L              # max clean pre-window length per switch
n_sims  <- { v <- suppressWarnings(as.integer(Sys.getenv("N_SIMS")))
  if (!is.na(v) && v >= 1L) v else if (nzchar(Sys.getenv("QUICK"))) 100L else 500L }
strength <- { s <- Sys.getenv("STRENGTH"); if (nzchar(s)) s else "high" }

# ---- design: LATE switches + PERMANENT stable-control pools -----------------
# Permanent pools are stable over the whole panel, so they are valid controls
# for any pre-window; late switch dates guarantee a long clean pre-window.
n_never  <- 150L  # D = 0 all t  (stable-untreated controls for joiners)
n_always <- 150L  # D = 1 all t  (stable-treated  controls for leavers)
join_at  <- c(20L,23L,26L,29L); n_join  <- 20L   # untreated then switch ON
leave_at <- c(21L,24L,27L,30L); n_leave <- 20L   # treated then switch OFF

build_design <- function() {
  blocks <- list(matrix(0L, n_never, T_total), matrix(1L, n_always, T_total))
  for (g in join_at)  { M <- matrix(0L, n_join,  T_total); M[, g:T_total] <- 1L; blocks[[length(blocks)+1L]] <- M }
  for (g in leave_at) { M <- matrix(1L, n_leave, T_total); M[, g:T_total] <- 0L; blocks[[length(blocks)+1L]] <- M }
  do.call(rbind, blocks)
}
D <- build_design(); N <- nrow(D)
never_idx  <- 1:n_never
always_idx <- (n_never + 1L):(n_never + n_always)
# switch cohorts: (type, switch date, switcher unit indices, control pool)
row0 <- n_never + n_always
switch_cohorts <- list(); r <- row0
for (g in join_at)  { u <- (r+1L):(r+n_join);  r <- r+n_join
  switch_cohorts[[length(switch_cohorts)+1L]] <- list(ty=+1L, t=g, sw=u, ct=never_idx) }
for (g in leave_at) { u <- (r+1L):(r+n_leave); r <- r+n_leave
  switch_cohorts[[length(switch_cohorts)+1L]] <- list(ty=-1L, t=g, sw=u, ct=always_idx) }

alpha  <- rnorm(N, 0, 2)
lambda <- 0.2 * (1:T_total)
tau    <- delta * D                        # no carryover

# ---- arbitrary serial correlation (mixture of AR kernels; from paper code) --
build_arbitrary_serial_cov <- function(TT, strength = c("low","medium","high")) {
  strength <- match.arg(strength)
  spec <- switch(strength,
    low    = list(rho = c(0.40,0.10,-0.20), w = c(0.45,0.40,0.15)),
    medium = list(rho = c(0.65,0.35,-0.20), w = c(0.55,0.30,0.15)),
    high   = list(rho = c(0.92,0.70,-0.25), w = c(0.60,0.30,0.10)))
  lags <- 0:(TT-1L); S <- matrix(0, TT, TT)
  for (m in seq_along(spec$rho)) S <- S + spec$w[m] * toeplitz(spec$rho[m]^lags)
  d <- sqrt(diag(S)); S <- S / outer(d, d); (S + t(S))/2
}
Sigma_err <- sigma^2 * build_arbitrary_serial_cov(T_total, strength)
L_err     <- chol(Sigma_err)

# ---- over-identified moment set: (cohort, pre-lag k) ------------------------
# moment (c,k): (Ybar_sw,t - Ybar_sw,t-k) - (Ybar_ct,t - Ybar_ct,t-k),  E = ty*delta.
MOM <- list()
for (co in switch_cohorts) {
  K <- min(Kcap, co$t - 1L)
  for (k in 1:K) MOM[[length(MOM)+1L]] <- list(ty=co$ty, t=co$t, k=k, sw=co$sw, ct=co$ct)
}
M  <- length(MOM)
Qv <- sapply(MOM, `[[`, "ty")            # incidence +1 joiner / -1 leaver

coefW <- array(0, dim = c(M, N, T_total))
for (s in seq_len(M)) {
  cl <- MOM[[s]]; t <- cl$t; tp <- cl$t - cl$k
  nsw <- length(cl$sw); nct <- length(cl$ct)
  coefW[s, cl$sw, t]  <- coefW[s, cl$sw, t]  + 1/nsw
  coefW[s, cl$sw, tp] <- coefW[s, cl$sw, tp] - 1/nsw
  coefW[s, cl$ct, t]  <- coefW[s, cl$ct, t]  - 1/nct
  coefW[s, cl$ct, tp] <- coefW[s, cl$ct, tp] + 1/nct
}
is_adjacent <- vapply(MOM, function(x) x$k == 1L, logical(1))   # k=1 subset = DID_M moments

moment_vec <- function(Y) { m <- numeric(M); for (s in seq_len(M)) m[s] <- sum(coefW[s,,]*Y); m }

# moment covariance under a within-unit covariance Sig (T x T): sum_i Wi Sig Wi'
switch_omega <- function(Sig) {
  Om <- matrix(0, M, M)
  for (i in seq_len(N)) { Wi <- matrix(coefW[, i, ], nrow = M); Om <- Om + Wi %*% Sig %*% t(Wi) }
  (Om + t(Om))/2
}
Omega_sph <- switch_omega(diag(T_total))   # Sigma = I  (spherical); data-independent

toeplitz_from_resid <- function(Y, betahat) {
  Yt <- Y - betahat * D
  nu <- Yt - rowMeans(Yt) - rep(colMeans(Yt), each = N) + mean(Yt)
  sig <- numeric(T_total)
  for (d in 0:(T_total-1L)) {
    a <- nu[, 1:(T_total-d), drop=FALSE]; b <- nu[, (1+d):T_total, drop=FALSE]
    sig[d+1L] <- mean(a * b)
  }
  toeplitz(sig)
}

# ---- the four estimators ---------------------------------------------------
est_didm <- function(m) {                 # k=1 only, frequency (switcher-count) weight
  num <- 0; den <- 0
  for (s in which(is_adjacent)) {
    w <- length(MOM[[s]]$sw)
    num <- num + w * (MOM[[s]]$ty * m[s]); den <- den + w    # ty*m => +delta
  }
  num / den
}
gls_solve <- function(m, A) as.numeric(solve(t(Qv) %*% A %*% Qv, t(Qv) %*% A %*% m))
est_eq  <- function(m) sum(Qv * m) / sum(Qv * Qv)            # A = I
est_sph <- function(m) gls_solve(m, solve(Omega_sph))       # A = Omega_sph^-1
est_T   <- function(Y, m, iters = 2L) {                     # A = Omega_Toeplitz^-1 (iterated)
  b <- est_eq(m)
  for (j in seq_len(iters)) b <- gls_solve(m, solve(switch_omega(toeplitz_from_resid(Y, b))))
  b
}

# ---- simulation ------------------------------------------------------------
simulate_Y <- function()
  matrix(alpha, N, T_total) + matrix(lambda, N, T_total, byrow=TRUE) + tau +
  matrix(rnorm(N*T_total), N, T_total) %*% L_err

run <- function(n_sims) {
  d <- e <- s <- g <- numeric(n_sims)
  for (r in seq_len(n_sims)) {
    Y <- simulate_Y(); m <- moment_vec(Y)
    d[r] <- est_didm(m); e[r] <- est_eq(m); s[r] <- est_sph(m); g[r] <- est_T(Y, m)
  }
  list(didm=d, eq=e, sph=s, gmmT=g)
}

if (sys.nframe() == 0L) {
  cat(sprintf("N=%d  T=%d  switch cohorts=%d  moments M=%d (K<=%d)  strength=%s  n_sims=%d  true ATT=%.2f\n",
              N, T_total, length(switch_cohorts), M, Kcap, strength, n_sims, delta))
  o <- run(n_sims)
  sm <- function(x) sprintf("mean=%7.4f  bias=%+.4f  sd=%.4f", mean(x), mean(x)-delta, sd(x))
  cat("DID_M   :", sm(o$didm), "\n")
  cat("GMM_eq  :", sm(o$eq),   "\n")
  cat("GMM_sph :", sm(o$sph),  "\n")
  cat("GMM_T   :", sm(o$gmmT), "\n\n")
  r1 <- var(o$didm)/var(o$eq); r2 <- var(o$eq)/var(o$sph); r3 <- var(o$sph)/var(o$gmmT)
  cat("Variance-gain decomposition (each factor >1 => that channel helps):\n")
  cat(sprintf("  over-identification  var(DID_M)/var(GMM_eq)  = %.3f\n", r1))
  cat(sprintf("  spherical structure  var(GMM_eq)/var(GMM_sph)= %.3f\n", r2))
  cat(sprintf("  serial correlation   var(GMM_sph)/var(GMM_T) = %.3f   <- our contribution\n", r3))
  cat(sprintf("  TOTAL  var(DID_M)/var(GMM_T) = %.3f   (check r1*r2*r3 = %.3f)\n",
              var(o$didm)/var(o$gmmT), r1*r2*r3))
}
