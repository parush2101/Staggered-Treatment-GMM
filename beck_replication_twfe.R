###############################################################################
# Beck et al. (2010) Replication — Panel A (Full Sample, 49 States)
#
# Source: Beck, T., Levine, R., & Levkov, A. (2010). Big Bad Banks? The Winners
#   and Losers from Bank Deregulation in the United States. Journal of Finance.
#
# Specification (Arora & Bijani eq. 27):
#   ln(Gini)_st = alpha_s + lambda_t + beta * D_branch_st + eps_st
#
#   ln(Gini)_st : log Gini coefficient for state s at year t
#   D_branch_st : 1 if state s has deregulated intrastate bank branching by year t
#   alpha_s     : state fixed effects
#   lambda_t    : year fixed effects
#   Standard errors clustered by state
#
# Panel A: full panel of 49 U.S. states, 1976–2006 (1,519 state-year obs)
#   Includes 12 always-treated states (deregulated before 1976)
#
# Target (Table 5 / Table 6, Panel A):
#   beta_hat = -0.0213,  SE = 0.0076  (cluster-robust by state)
###############################################################################

library(data.table)
library(fixest)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
dt <- fread("panel_A_beck_replication (1).csv")

cat(sprintf("Panel A: %d observations, %d states, %d years\n",
            nrow(dt),
            dt[, uniqueN(state)],
            dt[, uniqueN(wrkyr)]))

# ---------------------------------------------------------------------------
# TWFE regression
# ---------------------------------------------------------------------------
twfe_panelA <- feols(
  ln_gini ~ D_branch | state + wrkyr,
  data    = dt,
  cluster = ~state
)

summary(twfe_panelA)

cat(sprintf(
  "\nTWFE coefficient on D_branch: %.4f  (SE = %.4f)\n",
  coef(twfe_panelA)["D_branch"],
  se(twfe_panelA)["D_branch"]
))
cat("Target (Table 5 / Table 6, Panel A): -0.0213  (SE = 0.0076)\n")
