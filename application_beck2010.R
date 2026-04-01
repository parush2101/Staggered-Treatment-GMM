###############################################################################
# Beck (2010) Application: Effect of Branch Banking Deregulation on Inequality
#
# Data:   panel_A_beck_replication.csv
#         49 US states x 31 years (1976-2006)
#         Outcome:   ln_gini (log Gini coefficient)
#         Treatment: D_branch (branch deregulation indicator)
#         Cohort:    branch_reform (year of deregulation)
#
# Estimators:
#   1. TWFE (vanilla)              — single scalar on D_branch
#   2. TWFE Flexible               — cohort x year cell dummies (interacted TWFE)
#   3. Gardner (did2s)             — two-stage DiD, scalar ATT
#   4. Callaway-Sant'Anna (2021)   — att_gt with not-yet-treated control
#   5. Sun-Abraham (2021)          — interaction-weighted estimator via sunab()
#
# Reference: Beck, T., Levine, R., & Levkov, A. (2010). Big Bad Banks?
#   The Winners and Losers from Bank Deregulation in the United States.
#   Journal of Finance, 65(5), 1637-1667.
###############################################################################

library(data.table)
library(fixest)
library(did2s)
library(did)

# ===========================================================================
# 1. Load Data
# ===========================================================================

dt <- fread("panel_A_beck_replication.csv")

# 49 states x 31 years (1976-2006); no never-treated states
# branch_reform: year of deregulation (some pre-sample: 1960, 1970, 1975, 1976)
# Within-sample treated cohorts: branch_reform in {1977, ..., 1999}

cat("States:", uniqueN(dt$state), "\n")
cat("Years: ", paste(range(dt$wrkyr), collapse = "-"), "\n")
cat("Within-sample treated cohorts:",
    paste(sort(unique(dt[branch_reform >= 1977, branch_reform])), collapse = " "), "\n\n")

# ===========================================================================
# 2. TWFE (Vanilla)
# ===========================================================================
# Single scalar coefficient on D_branch; unit + time FE.
# Aggregates CATTs using FWL-residualized weights — can be negative when
# already-treated units act as implicit controls.

m_twfe <- feols(ln_gini ~ D_branch | state + wrkyr, data = dt, cluster = ~state)

att_twfe <- as.numeric(coef(m_twfe)["D_branch"])
se_twfe  <- as.numeric(se(m_twfe)["D_branch"])

# ===========================================================================
# 3. TWFE Flexible (Interacted TWFE)
# ===========================================================================
# Separate CATT per (cohort, post-period) cell via cohort x year interactions.
# ATT reported as simple (cell-weighted) mean over all identified cells.
# Restricted to within-sample treated cohorts (branch_reform >= 1977).

dt_tw <- dt[branch_reform >= 1977]
dt_tw[, cell_id := fifelse(wrkyr >= branch_reform,
                            paste0(branch_reform, "_", wrkyr), "pre")]

m_flex <- feols(ln_gini ~ i(cell_id, ref = "pre") | state + wrkyr,
                data = dt_tw, cluster = ~state)

fc      <- coef(m_flex)
fn      <- names(fc)
keep    <- grepl("cell_id::", fn, fixed = TRUE)
cell_str <- sub("cell_id::", "", fn[keep])
fc_g    <- as.integer(sub("_.*", "", cell_str))
fc_t    <- as.integer(sub(".*_", "", cell_str))
post    <- fc_t >= fc_g

post_names    <- fn[keep][post]
cell_catts    <- fc[post_names]
n_cells       <- sum(post)
att_flex_cell <- mean(cell_catts)

# SE via delta method: ATT = w'beta, w = 1/n => Var(ATT) = w' V w
# V is the full clustered VCV of post-treatment cell CATTs
V_post   <- vcov(m_flex)[post_names, post_names]
w        <- rep(1 / n_cells, n_cells)
se_flex  <- as.numeric(sqrt(t(w) %*% V_post %*% w))

# ===========================================================================
# 4. Gardner (2021) Two-Stage DiD
# ===========================================================================
# Stage 1: estimate unit + time FE on never/not-yet-treated observations only.
# Stage 2: regress residualised outcome on D_branch to recover ATT.
# Corrects for contamination bias from already-treated implicit controls.

m_gard <- did2s(data         = dt,
                yname        = "ln_gini",
                first_stage  = ~0 | state + wrkyr,
                second_stage = ~D_branch,
                treatment    = "D_branch",
                cluster_var  = "state")

att_gard <- as.numeric(coef(m_gard)["D_branch"])
se_gard  <- as.numeric(se(m_gard)["D_branch"])

# ===========================================================================
# 5. Callaway-Sant'Anna (2021)
# ===========================================================================
# att_gt estimates CATT(g,t) for each cohort-period pair using not-yet-treated
# units as the comparison group. Pre-sample-treated states (branch_reform < 1977)
# are recoded as 0 (never-treated) so CS treats them as always-controls.
# ATT aggregated via simple average over all (g,t) post-treatment cells.

dt_cs <- copy(dt)
dt_cs[branch_reform < 1977, branch_reform := 0L]

cs_out <- att_gt(yname         = "ln_gini",
                 tname         = "wrkyr",
                 idname        = "statefip",
                 gname         = "branch_reform",
                 control_group = "notyettreated",
                 data          = as.data.frame(dt_cs),
                 est_method    = "reg",
                 print_details = FALSE)

cs_agg   <- aggte(cs_out, type = "simple")
att_cs   <- cs_agg$overall.att
se_cs    <- cs_agg$overall.se

# ===========================================================================
# 6. Sun-Abraham (2021)
# ===========================================================================
# Interaction-weighted estimator via fixest::sunab(). Pre-sample-treated
# states recoded to 0 (used as never-treated base group).
# ATT aggregated as the standard interaction-weighted average across cohorts.

dt[, g_sa := fifelse(branch_reform >= 1977 & branch_reform <= 2006,
                      branch_reform, 0L)]

m_sa    <- feols(ln_gini ~ sunab(g_sa, wrkyr) | state + wrkyr,
                 data = dt, cluster = ~state)
att_sa  <- as.numeric(aggregate(m_sa, agg = "ATT")[1])
se_sa   <- as.numeric(aggregate(m_sa, agg = "ATT")[2])

# ===========================================================================
# 7. Results Table
# ===========================================================================

sig <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", "")))
pval <- function(est, se) 2 * pnorm(-abs(est / se))

p_twfe <- pval(att_twfe,      se_twfe)
p_flex <- pval(att_flex_cell, se_flex)
p_gard <- pval(att_gard,      se_gard)
p_cs   <- pval(att_cs,        se_cs)
p_sa   <- pval(att_sa,        se_sa)

cat(strrep("=", 76), "\n")
cat("  Beck (2010): Effect of Branch Banking Deregulation on ln(Gini)\n")
cat("  N = 49 states x 31 years (1976-2006)\n")
cat(strrep("=", 76), "\n")
cat(sprintf("  %-30s  %8s  %8s  %7s\n", "Estimator", "ATT", "SE", "p"))
cat(strrep("-", 60), "\n")
cat(sprintf("  %-30s  %8.4f  %8.4f  %7.4f  %s\n",
            "TWFE (vanilla)",            att_twfe,      se_twfe, p_twfe, sig(p_twfe)))
cat(sprintf("  %-30s  %8.4f  %8.4f  %7.4f  %s\n",
            "Sun-Abraham (2021)",        att_sa,        se_sa,   p_sa,   sig(p_sa)))
cat(sprintf("  %-30s  %8.4f  %8.4f  %7.4f  %s\n",
            "Callaway-Sant'Anna (2021)", att_cs,        se_cs,   p_cs,   sig(p_cs)))
cat(sprintf("  %-30s  %8.4f  %8.4f  %7.4f  %s\n",
            "Gardner (2021)",            att_gard,      se_gard, p_gard, sig(p_gard)))
cat(sprintf("  %-30s  %8.4f  %8.4f  %7.4f  %s\n",
            sprintf("TWFE Flexible (%d cells)", n_cells),
            att_flex_cell, se_flex, p_flex, sig(p_flex)))
cat(strrep("=", 76), "\n")
cat("  SE for TWFE Flexible: delta method sqrt(w'Vw), w=1/n, V clustered by state\n")
cat("  CS: not-yet-treated control, simple aggregation; SA: interaction-weighted\n\n")

cat("Cell CATT distribution (TWFE Flexible):\n")
cat(sprintf("  Min=%.4f  p25=%.4f  Median=%.4f  p75=%.4f  Max=%.4f\n\n",
            min(cell_catts), quantile(cell_catts, 0.25),
            median(cell_catts), quantile(cell_catts, 0.75),
            max(cell_catts)))
