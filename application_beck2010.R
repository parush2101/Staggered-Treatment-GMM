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
#   1. TWFE (vanilla)       — single scalar on D_branch
#   2. TWFE Flexible        — cohort x year cell dummies (interacted TWFE)
#   3. Gardner (did2s)      — two-stage DiD, scalar ATT
#
# Reference: Beck, T., Levine, R., & Levkov, A. (2010). Big Bad Banks?
#   The Winners and Losers from Bank Deregulation in the United States.
#   Journal of Finance, 65(5), 1637-1667.
###############################################################################

library(data.table)
library(fixest)
library(did2s)

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

cell_catts    <- fc[keep][post]
n_cells       <- sum(post)
att_flex_cell <- mean(cell_catts)

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
# 5. Results Table
# ===========================================================================

p_twfe <- 2 * pnorm(-abs(att_twfe / se_twfe))
p_gard <- 2 * pnorm(-abs(att_gard / se_gard))

sig <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", "")))

cat(strrep("=", 72), "\n")
cat("  Beck (2010): Effect of Branch Banking Deregulation on ln(Gini)\n")
cat("  N = 49 states x 31 years (1976-2006)\n")
cat(strrep("=", 72), "\n")
cat(sprintf("  %-28s  %9s  %9s  %6s  %s\n",
            "Estimator", "ATT", "SE", "p", "Aggregation"))
cat(strrep("-", 72), "\n")
cat(sprintf("  %-28s  %9.4f  %9.4f  %6.4f  %s\n",
            "TWFE (vanilla)", att_twfe, se_twfe, p_twfe,
            paste0("FWL weights ", sig(p_twfe))))
cat(sprintf("  %-28s  %9.4f  %9s  %6s  %s\n",
            sprintf("TWFE Flexible (%d cells)", n_cells),
            att_flex_cell, "—", "—",
            "simple cell mean"))
cat(sprintf("  %-28s  %9.4f  %9.4f  %6.4f  %s\n",
            "Gardner (2021)", att_gard, se_gard, p_gard,
            paste0("FWL weights ", sig(p_gard))))
cat(strrep("=", 72), "\n\n")

cat("Cell CATT distribution (TWFE Flexible):\n")
cat(sprintf("  Min=%.4f  p25=%.4f  Median=%.4f  p75=%.4f  Max=%.4f\n\n",
            min(cell_catts), quantile(cell_catts, 0.25),
            median(cell_catts), quantile(cell_catts, 0.75),
            max(cell_catts)))

cat("Notes:\n")
cat("  - TWFE and Gardner use the same FWL aggregation weights but differ\n")
cat("    in how fixed effects are estimated: TWFE uses all obs (incl. treated),\n")
cat("    Gardner stage 1 uses only untreated obs. Sign flips: -0.0213 -> +0.0195.\n")
cat("  - TWFE Flexible cell mean avoids negative weights; most CATTs are positive.\n")
cat("  - 8 cohort-1999 post-treatment cells dropped by fixest (collinear with year FE).\n")
