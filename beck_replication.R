suppressPackageStartupMessages({
  library(fixest); library(did); library(did2s); library(etwfe); library(dplyr)
})

df <- read.csv('/Users/parusharora/Downloads/panel_A_beck_replication.csv')

# ── variable prep ──────────────────────────────────────────────────────────────
df$gvar       <- ifelse(df$branch_reform < 1976, 0,   df$branch_reform)
df$cohort_sa  <- ifelse(df$branch_reform < 1976, Inf, as.numeric(df$branch_reform))
df_37         <- df[df$branch_reform >= 1976, ]   # 37 states: exclude always-treated (< 1976)

cat("============================================================\n")
cat("Panel A (full sample, 49 states) — Beck et al. (2010)\n")
cat("============================================================\n\n")

# ── 1. TWFE ───────────────────────────────────────────────────────────────────
twfe <- feols(ln_gini ~ D_branch | statefip + wrkyr, data = df, cluster = ~statefip)
cat(sprintf("1. TWFE:          ATT = %7.4f  SE = %.4f  (paper: -0.0213, 0.0036)\n",
            coef(twfe)["D_branch"], se(twfe)["D_branch"]))

# ── 2. CS (Callaway & Sant'Anna 2021) ─────────────────────────────────────────
# Use 37-state sample (branch_reform >= 1976): exclude always-treated states
# not-yet-treated control group; DR estimator; group-aggregation matches paper
cs_out <- att_gt("ln_gini", "wrkyr", "statefip", "branch_reform",
                 data = df_37, control_group = "notyettreated",
                 est_method = "dr", bstrap = FALSE,
                 print_details = FALSE)
cs_agg <- aggte(cs_out, type = "group", na.rm = TRUE)
cat(sprintf("2. CS:            ATT = %7.4f  SE = %.4f  (paper: -0.0101, 0.0084)\n",
            cs_agg$overall.att, cs_agg$overall.se))

# ── 3. SA (Sun & Abraham 2021) ────────────────────────────────────────────────
sa_mod <- feols(ln_gini ~ sunab(cohort_sa, wrkyr) | statefip + wrkyr,
                data = df, cluster = ~statefip)
sa_agg <- aggregate(sa_mod, agg = "ATT")
cat(sprintf("3. SA:            ATT = %7.4f  SE = %.4f  (paper: -0.030,  0.0217)\n",
            sa_agg[1,"Estimate"], sa_agg[1,"Std. Error"]))

# ── 4. Gardner et al. (2024) ──────────────────────────────────────────────────
# NOTE: did2s v1.0.2 gives +0.0195 vs paper -0.0081.
# Root causes: (a) time FEs unidentified for 1999-2006 (no control obs post-1998),
# (b) 12 always-treated states lack first-stage unit FE, inflating early-year
# time FEs and over-predicting counterfactuals.
gard <- did2s(df, yname = "ln_gini",
              first_stage  = ~0 | statefip + wrkyr,
              second_stage = ~i(D_branch),
              treatment    = "D_branch",
              cluster_var  = "statefip")
cat(sprintf("4. Gardner:       ATT = %7.4f  SE = %.4f  (paper: -0.0081, 0.002)\n",
            coef(gard)["D_branch::1"], se(gard)["D_branch::1"]))

# ── 5. Flexible TWFE (Wooldridge 2025) via etwfe ─────────────────────────────
# Use 37-state sample: all states have branch_reform >= 1976 so etwfe uses
# not-yet-treated as comparison (no gvar=0 group needed).
etwfe_37 <- etwfe(fml = ln_gini ~ 1, tvar = wrkyr, gvar = branch_reform,
                  data = df_37, vcov = ~statefip)

coefs_37 <- coef(etwfe_37)
cn_37    <- names(coefs_37)

# Parse cohort (g) and time (t) from coefficient names like:
# ".Dtreat:branch_reform::1985:wrkyr::1988"
g_37 <- as.integer(regmatches(cn_37, regexpr("(?<=branch_reform::)[0-9]+", cn_37, perl=TRUE)))
t_37 <- as.integer(regmatches(cn_37, regexpr("(?<=wrkyr::)[0-9]+",         cn_37, perl=TRUE)))

cf37 <- data.frame(name = cn_37, beta = as.numeric(coefs_37),
                   g = g_37, t = t_37, stringsAsFactors = FALSE)

cohort_sz_37 <- df_37 %>% distinct(statefip, branch_reform) %>%
  count(branch_reform, name = "n_units")

cf37      <- merge(cf37, cohort_sz_37, by.x = "g", by.y = "branch_reform", all.x = TRUE)
cf37      <- cf37[!is.na(cf37$n_units), ]

total_w   <- sum(cf37$n_units)
att_flex  <- sum(cf37$beta * cf37$n_units) / total_w
W         <- cf37$n_units / total_w
vcv       <- vcov(etwfe_37)[cf37$name, cf37$name]
se_flex   <- sqrt(as.numeric(t(W) %*% vcv %*% W))

cat(sprintf("5. Flexible TWFE: ATT = %7.4f  SE = %.4f  (paper:  0.0052, 0.0095)\n",
            att_flex, se_flex))

cat("\n============================================================\n")
cat("Paper Table 6, Panel A:\n")
cat("  TWFE:          -0.0213  (0.0036)\n")
cat("  CS:            -0.0101  (0.0084)\n")
cat("  SA:            -0.030   (0.0217)\n")
cat("  Gardner:       -0.0081  (0.002)\n")
cat("  Flexible TWFE:  0.0052  (0.0095)\n")
cat("============================================================\n")
