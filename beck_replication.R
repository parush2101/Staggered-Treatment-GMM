suppressPackageStartupMessages({
  library(fixest); library(did); library(did2s); library(etwfe)
  library(dplyr); library(MASS)
})

# ── helper: unit-weighted ATT from etwfe ──────────────────────────────────────
flex_att <- function(etwfe_mod, df_cs) {
  coefs <- coef(etwfe_mod); cn <- names(coefs)
  g_v <- as.integer(regmatches(cn, regexpr("(?<=branch_reform::)[0-9]+", cn, perl=TRUE)))
  t_v <- as.integer(regmatches(cn, regexpr("(?<=wrkyr::)[0-9]+",         cn, perl=TRUE)))
  cf  <- data.frame(name=cn, beta=as.numeric(coefs), g=g_v, t=t_v,
                    stringsAsFactors=FALSE)
  csz <- df_cs %>% distinct(statefip, branch_reform) %>%
    count(branch_reform, name="n_units")
  cf  <- merge(cf, csz, by.x="g", by.y="branch_reform", all.x=TRUE)
  cf  <- cf[!is.na(cf$n_units), ]
  W   <- cf$n_units / sum(cf$n_units)
  vcv <- vcov(etwfe_mod)[cf$name, cf$name]
  list(att = sum(cf$beta * W),
       se  = sqrt(as.numeric(t(W) %*% vcv %*% W)))
}

# ── helper: iterated efficient GMM ────────────────────────────────────────────
# Uses not-yet-treated comparisons only; all pre-periods as base periods.
# Omega estimated from autocovariances of two-way-demeaned treatment-adjusted residuals.
# ATT = simple cell-mean (equal weight per identified CATT).
run_gmm <- function(df, n_iter = 4) {
  T_yr  <- 1976:2006
  T_len <- length(T_yr)
  yi    <- function(y) y - 1975L          # year -> column index (1-based)

  # State-cohort structure, ordered by statefip
  sc        <- df %>% distinct(statefip, branch_reform) %>% arrange(statefip)
  N_st      <- nrow(sc)
  unit_coh  <- sc$branch_reform
  # Only cohorts with at least one pre-period in sample (excludes 1976/NY if present)
  treated_g <- sort(unique(sc$branch_reform[sc$branch_reform >= 1977]))
  T_end_id  <- max(treated_g) - 1L        # last year with a not-yet-treated control

  coh_n <- sc %>% filter(branch_reform %in% treated_g) %>%
    count(branch_reform, name = "n")
  get_N <- function(g) coh_n$n[coh_n$branch_reform == g]

  # ── CATT list: (g, t_post) with t_post in [g, T_end_id] ───────────────────
  catt_list <- vector("list", 500L)
  n_catt    <- 0L
  for (g in treated_g) {
    if (g > T_end_id) next
    for (tp in g:T_end_id) {
      n_catt <- n_catt + 1L
      catt_list[[n_catt]] <- c(g, tp)
    }
  }
  catt_list <- catt_list[seq_len(n_catt)]
  catt_g    <- sapply(catt_list, `[`, 1)
  catt_t    <- sapply(catt_list, `[`, 2)

  # ── DiD moments: not-yet-treated, all pre-periods ──────────────────────────
  did_list <- vector("list", 10000L)
  n_did    <- 0L
  for (ci in seq_len(n_catt)) {
    g_c    <- catt_list[[ci]][1]; t_p <- catt_list[[ci]][2]
    notyet <- treated_g[treated_g > t_p]
    if (!length(notyet)) next
    for (t_r in 1976L:(g_c - 1L)) {
      for (g_l in notyet) {
        n_did <- n_did + 1L
        did_list[[n_did]] <- list(ci=ci, gf=g_c, gc=g_l, tp=t_p, tr=t_r)
      }
    }
  }
  did_list <- did_list[seq_len(n_did)]
  cat(sprintf("    %d CATTs, %d DiD moments\n", n_catt, n_did))

  meta_gf <- sapply(did_list, `[[`, "gf")
  meta_gc <- sapply(did_list, `[[`, "gc")
  meta_tp <- sapply(did_list, `[[`, "tp")
  meta_tr <- sapply(did_list, `[[`, "tr")

  # ── Q_H: one +1 per row (not-yet-treated → no bias corrections needed) ─────
  Q_H <- matrix(0, nrow=n_did, ncol=n_catt)
  for (s in seq_len(n_did)) Q_H[s, did_list[[s]]$ci] <- 1L
  Q_H     <- matrix(as.numeric(Q_H), nrow=n_did)
  QtQ_inv <- tryCatch(solve(crossprod(Q_H)), error=function(e) ginv(crossprod(Q_H)))

  # ── C_mat: cohort-pair covariance structure ─────────────────────────────────
  N_f <- sapply(meta_gf, get_N); N_c <- sapply(meta_gc, get_N)
  gg  <- outer(meta_gf, meta_gf, "=="); gc_ <- outer(meta_gf, meta_gc, "==")
  cg  <- outer(meta_gc, meta_gf, "=="); cc  <- outer(meta_gc, meta_gc, "==")
  C_mat <- sweep(gg - gc_, 1, 1/N_f, "*") + sweep(cc - cg, 1, 1/N_c, "*")
  rm(gg, gc_, cg, cc); invisible(gc())

  # ── Lag indices for autocovariance lookup ──────────────────────────────────
  pp <- as.vector(abs(outer(meta_tp, meta_tp, "-")))
  pr <- as.vector(abs(outer(meta_tp, meta_tr, "-")))
  rp <- as.vector(abs(outer(meta_tr, meta_tp, "-")))
  rr <- as.vector(abs(outer(meta_tr, meta_tr, "-")))

  # ── Y matrix (N_st x T_len) and cohort means ───────────────────────────────
  Y_mat <- matrix(NA_real_, nrow=N_st, ncol=T_len)
  for (i in seq_len(N_st)) {
    ri <- df[df$statefip == sc$statefip[i], ]
    Y_mat[i, ] <- ri$ln_gini[order(ri$wrkyr)]
  }
  all_g <- treated_g
  Ybar  <- matrix(NA_real_, nrow=length(all_g), ncol=T_len)
  for (ci in seq_along(all_g)) {
    us <- which(unit_coh == all_g[ci])
    Ybar[ci, ] <- colMeans(Y_mat[us, , drop=FALSE])
  }
  gidx <- function(g) which(all_g == g)

  # ── Delta: observed 2×2 DiD estimates ──────────────────────────────────────
  Delta <- numeric(n_did)
  for (s in seq_len(n_did)) {
    d <- did_list[[s]]
    Delta[s] <- (Ybar[gidx(d$gf), yi(d$tp)] - Ybar[gidx(d$gf), yi(d$tr)]) -
                (Ybar[gidx(d$gc), yi(d$tp)] - Ybar[gidx(d$gc), yi(d$tr)])
  }

  # ── Iterated GMM ────────────────────────────────────────────────────────────
  beta  <- as.numeric(QtQ_inv %*% crossprod(Q_H, Delta))  # OLS GMM start
  Omega <- diag(n_did)

  for (iter in seq_len(n_iter)) {
    beta_old <- beta

    # Subtract current treatment effect estimates from treated cells
    Y_adj <- Y_mat
    for (k in seq_len(n_catt)) {
      us <- which(unit_coh == catt_g[k])
      Y_adj[us, yi(catt_t[k])] <- Y_adj[us, yi(catt_t[k])] - beta[k]
    }

    # Two-way demean
    rm_ <- rowMeans(Y_adj); cm_ <- colMeans(Y_adj); gm_ <- mean(Y_adj)
    res <- Y_adj - outer(rm_, rep(1,T_len)) - outer(rep(1,N_st), cm_) + gm_

    # Autocovariances sigma(d) for d = 0,...,T-1
    sig <- numeric(T_len)
    for (d in 0:(T_len - 1)) {
      r1 <- seq_len(T_len - d)
      sig[d+1] <- sum(res[, r1] * res[, r1+d]) / (N_st * (T_len - d))
    }

    # Build Omega = C_mat * S (element-wise)
    S     <- sig[pp+1] - sig[pr+1] - sig[rp+1] + sig[rr+1]
    Omega <- C_mat * matrix(S, nrow=n_did)
    Omega <- (Omega + t(Omega)) / 2
    diag(Omega) <- diag(Omega) + 1e-6

    OQ <- tryCatch(solve(Omega, Q_H), error=function(e) NULL)
    if (is.null(OQ)) { cat("    solve failed iter", iter, "\n"); break }
    OD   <- solve(Omega, Delta)
    beta <- tryCatch(
      as.numeric(solve(crossprod(Q_H, OQ), crossprod(Q_H, OD))),
      error = function(e) beta_old
    )
    conv <- max(abs(beta - beta_old))
    cat(sprintf("    iter %d: max|Δβ| = %.2e\n", iter, conv))
    if (conv < 1e-8) break
  }

  # ── ATT and SE ──────────────────────────────────────────────────────────────
  w     <- rep(1/n_catt, n_catt)
  att   <- sum(w * beta)
  OQ_f  <- tryCatch(solve(Omega, Q_H), error=function(e) NULL)
  V_eff <- if (!is.null(OQ_f))
    tryCatch(solve(crossprod(Q_H, OQ_f)), error=function(e) ginv(crossprod(Q_H, OQ_f)))
    else ginv(crossprod(Q_H))
  se <- sqrt(as.numeric(t(w) %*% V_eff %*% w))

  list(att=att, se=se, n_catt=n_catt, n_did=n_did)
}

# ── Panel A data ───────────────────────────────────────────────────────────────
dfA          <- read.csv('/Users/parusharora/Downloads/panel_A_beck_replication.csv')
dfA$cohort_sa <- ifelse(dfA$branch_reform < 1976, Inf, as.numeric(dfA$branch_reform))
dfA_37       <- dfA[dfA$branch_reform >= 1976, ]

# ── Panel B data ───────────────────────────────────────────────────────────────
dfB          <- read.csv('/Users/parusharora/panel_B_modern_did_eligible.csv')
dfB$cohort_sa <- dfB$branch_reform

# ═══════════════════════════════════════════════════════════════════════════════
cat("==============================================================\n")
cat("Panel A (49 states, 1976-2006) — Beck et al. (2010)\n")
cat("==============================================================\n\n")

twfe_A <- feols(ln_gini ~ D_branch | statefip + wrkyr, data=dfA, cluster=~statefip)
cat(sprintf("1. TWFE:          ATT = %7.4f  SE = %.4f\n",
            coef(twfe_A)["D_branch"], se(twfe_A)["D_branch"]))

cs_A    <- att_gt("ln_gini","wrkyr","statefip","branch_reform",
                  data=dfA_37, control_group="notyettreated",
                  est_method="dr", bstrap=FALSE, print_details=FALSE)
cs_aggA <- aggte(cs_A, type="group", na.rm=TRUE)
cat(sprintf("2. CS:            ATT = %7.4f  SE = %.4f\n",
            cs_aggA$overall.att, cs_aggA$overall.se))

sa_A    <- feols(ln_gini ~ sunab(cohort_sa, wrkyr) | statefip + wrkyr,
                 data=dfA, cluster=~statefip)
sa_aggA <- aggregate(sa_A, agg="ATT")
cat(sprintf("3. SA:            ATT = %7.4f  SE = %.4f\n",
            sa_aggA[1,"Estimate"], sa_aggA[1,"Std. Error"]))

gard_A <- did2s(dfA, yname="ln_gini",
                first_stage=~0|statefip+wrkyr,
                second_stage=~i(D_branch),
                treatment="D_branch", cluster_var="statefip")
cat(sprintf("4. Gardner:       ATT = %7.4f  SE = %.4f\n",
            coef(gard_A)["D_branch::1"], se(gard_A)["D_branch::1"]))

etwfe_A <- etwfe(fml=ln_gini~1, tvar=wrkyr, gvar=branch_reform,
                 data=dfA_37, vcov=~statefip)
fa <- flex_att(etwfe_A, dfA_37)
cat(sprintf("5. Flexible TWFE: ATT = %7.4f  SE = %.4f\n", fa$att, fa$se))

cat("6. Efficient GMM: building...\n")
gmm_A <- run_gmm(dfA_37, n_iter=4)
cat(sprintf("6. Efficient GMM: ATT = %7.4f  SE = %.4f  (%d CATTs, %d DiDs)\n\n",
            gmm_A$att, gmm_A$se, gmm_A$n_catt, gmm_A$n_did))

# ═══════════════════════════════════════════════════════════════════════════════
cat("==============================================================\n")
cat("Panel B (36 within-sample treated states, 1976-2006)\n")
cat("==============================================================\n\n")

twfe_B <- feols(ln_gini ~ D_branch | statefip + wrkyr, data=dfB, cluster=~statefip)
cat(sprintf("1. TWFE:          ATT = %7.4f  SE = %.4f\n",
            coef(twfe_B)["D_branch"], se(twfe_B)["D_branch"]))

cs_B    <- att_gt("ln_gini","wrkyr","statefip","branch_reform",
                  data=dfB, control_group="notyettreated",
                  est_method="dr", bstrap=FALSE, print_details=FALSE)
cs_aggB <- aggte(cs_B, type="group", na.rm=TRUE)
cat(sprintf("2. CS:            ATT = %7.4f  SE = %.4f\n",
            cs_aggB$overall.att, cs_aggB$overall.se))

sa_B    <- feols(ln_gini ~ sunab(cohort_sa, wrkyr) | statefip + wrkyr,
                 data=dfB, cluster=~statefip)
sa_aggB <- aggregate(sa_B, agg="ATT")
cat(sprintf("3. SA:            ATT = %7.4f  SE = %.4f\n",
            sa_aggB[1,"Estimate"], sa_aggB[1,"Std. Error"]))

gard_B <- did2s(dfB, yname="ln_gini",
                first_stage=~0|statefip+wrkyr,
                second_stage=~i(D_branch),
                treatment="D_branch", cluster_var="statefip")
cat(sprintf("4. Gardner:       ATT = %7.4f  SE = %.4f\n",
            coef(gard_B)["D_branch::1"], se(gard_B)["D_branch::1"]))

etwfe_B <- etwfe(fml=ln_gini~1, tvar=wrkyr, gvar=branch_reform,
                 data=dfB, vcov=~statefip)
fb <- flex_att(etwfe_B, dfB)
cat(sprintf("5. Flexible TWFE: ATT = %7.4f  SE = %.4f\n", fb$att, fb$se))

cat("6. Efficient GMM: building...\n")
gmm_B <- run_gmm(dfB, n_iter=4)
cat(sprintf("6. Efficient GMM: ATT = %7.4f  SE = %.4f  (%d CATTs, %d DiDs)\n\n",
            gmm_B$att, gmm_B$se, gmm_B$n_catt, gmm_B$n_did))
