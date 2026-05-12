###############################################################################
# Beck et al. (2010) Replication — Panel A (Full Sample, 49 States)
#
# Source: Beck, T., Levine, R., & Levkov, A. (2010). Big Bad Banks? The Winners
#   and Losers from Bank Deregulation in the United States. Journal of Finance.
#
# Specification (Arora & Bijani eq. 27):
#   ln(Gini)_st = alpha_s + lambda_t + beta * D_branch_st + eps_st
#
# Panel A: 49 U.S. states, 1976-2006, 1,519 state-year obs
#   Includes 12 always-treated states (deregulated before 1976)
#   Always-treated handling:
#     CS        -> g_cs  = 0    (coded as never-treated control)
#     SA / Flex -> g_inf = Inf  (coded as never-treated reference)
#     Gardner   -> first_treat = Inf (contributes to first-stage FE as control)
#
# Estimators (Table 6, Panel A targets):
#   1. Pooled TWFE           -0.0213  (0.0076)
#   2. Callaway-Sant'Anna    -0.0101  (0.0078)
#   3. Sun-Abraham           -0.0354  (0.0114)
#   4. Gardner / did2s       +0.0195  (0.0067)
#   5. Flexible TWFE         -0.0165  (0.0141)
###############################################################################

library(data.table)
library(fixest)
library(did)
library(did2s)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
dt <- fread("panel_A_beck_replication (1).csv")

cat(sprintf("Panel A: %d obs, %d states, %d years\n",
            nrow(dt), dt[, uniqueN(state)], dt[, uniqueN(wrkyr)]))

first_yr <- min(dt$wrkyr)   # 1976

# Cohort variables — always-treated states have branch_reform < 1976
dt[, g_cs        := fifelse(branch_reform < first_yr, 0L,  as.integer(branch_reform))]
dt[, g_inf       := fifelse(branch_reform < first_yr, Inf, as.numeric(branch_reform))]
dt[, first_treat := fifelse(branch_reform < first_yr, Inf, as.numeric(branch_reform))]

# ===========================================================================
# 1. Pooled TWFE
# ===========================================================================
twfe <- feols(ln_gini ~ D_branch | state + wrkyr,
              data = dt, cluster = ~state)

# ===========================================================================
# 2. Callaway-Sant'Anna (2021)
#    control_group = "notyettreated"; always-treated (g_cs=0) act as controls
# ===========================================================================
cs_out <- att_gt(
  yname         = "ln_gini",
  tname         = "wrkyr",
  idname        = "statefip",
  gname         = "g_cs",
  data          = as.data.frame(dt),
  control_group = "notyettreated",
  est_method    = "dr",
  print_details = FALSE,
  bstrap        = FALSE,
  cband         = FALSE
)
cs_agg <- aggte(cs_out, type = "simple")

# ===========================================================================
# 3. Sun-Abraham (2021)
#    Always-treated coded as g_inf = Inf (never-treated reference in sunab)
# ===========================================================================
sa_mod <- feols(ln_gini ~ sunab(g_inf, wrkyr) | state + wrkyr,
                data = dt, cluster = ~state)
sa_agg <- summary(sa_mod, agg = "ATT")$coeftable[1, ]

# ===========================================================================
# 4. Gardner / did2s (2024)
#    first_stage uses only untreated obs; always-treated (first_treat=Inf)
#    contribute all periods to the first-stage FE as if never-treated
# ===========================================================================
gardner_mod <- did2s(
  data         = as.data.frame(dt),
  yname        = "ln_gini",
  first_stage  = ~ 0 | state + wrkyr,
  second_stage = ~ i(D_branch, ref = 0),
  treatment    = "D_branch",
  cluster_var  = "state",
  verbose      = FALSE
)

# ===========================================================================
# 5. Flexible TWFE (Wooldridge 2025)
#    Fully interacts treatment with cohort x time for within-sample cohorts.
#    treat_gt = g*100 + t when unit is actively treated by a within-sample
#    cohort, 0 otherwise (reference: untreated + always-treated observations).
#    ATT aggregated as equal-weighted mean; SE via delta method.
# ===========================================================================
dt[, treat_gt := fifelse(
  D_branch == 1L & !is.infinite(g_inf),
  as.integer(g_inf) * 100L + as.integer(wrkyr),
  0L
)]

flex_mod  <- feols(ln_gini ~ i(treat_gt, ref = 0) | state + wrkyr,
                   data = dt, cluster = ~state)

flex_coefs <- coef(flex_mod)
flex_vcov  <- vcov(flex_mod)
n_coef     <- length(flex_coefs)
w          <- rep(1 / n_coef, n_coef)
flex_att   <- as.numeric(w %*% flex_coefs)
flex_se    <- as.numeric(sqrt(t(w) %*% flex_vcov %*% w))

# ===========================================================================
# Summary table
# ===========================================================================
paper_att <- c(-0.0213, -0.0101, -0.0354,  0.0195, -0.0165)
paper_se  <- c( 0.0076,  0.0078,  0.0114,  0.0067,  0.0141)

est_names <- c("Pooled TWFE", "CS", "Sun-Abraham", "Gardner", "Flex TWFE")
att_vals  <- c(
  coef(twfe)["D_branch"],
  cs_agg$overall.att,
  sa_agg[1],
  coef(gardner_mod)["D_branch::1"],
  flex_att
)
se_vals <- c(
  se(twfe)["D_branch"],
  cs_agg$overall.se,
  sa_agg[2],
  se(gardner_mod)["D_branch::1"],
  flex_se
)

cat("\n")
cat("=================================================================\n")
cat("  Table 6 Replication — Panel A (49 states, 1,519 obs)\n")
cat("  Outcome: ln(Gini)  |  Treatment: D_branch  |  SE by state\n")
cat("=================================================================\n")
cat(sprintf("  %-14s  %9s  %9s  %9s  %9s\n",
            "Estimator", "ATT", "SE", "Paper ATT", "Paper SE"))
cat("  ", paste(rep("-", 58), collapse = ""), "\n", sep = "")
for (i in seq_along(est_names)) {
  cat(sprintf("  %-14s  %9.4f  %9.4f  %9.4f  %9.4f\n",
              est_names[i], att_vals[i], se_vals[i],
              paper_att[i], paper_se[i]))
}
cat("\n")
