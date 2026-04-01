suppressPackageStartupMessages({
  library(fixest); library(did); library(did2s); library(etwfe); library(dplyr)
})

# ── helpers ───────────────────────────────────────────────────────────────────
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

# ── Panel A data ──────────────────────────────────────────────────────────────
dfA     <- read.csv('/Users/parusharora/Downloads/panel_A_beck_replication.csv')
dfA$cohort_sa <- ifelse(dfA$branch_reform < 1976, Inf, as.numeric(dfA$branch_reform))
dfA_37  <- dfA[dfA$branch_reform >= 1976, ]   # 37 states (drop always-treated)

# ── Panel B data ──────────────────────────────────────────────────────────────
dfB <- read.csv('/Users/parusharora/panel_B_modern_did_eligible.csv')
dfB$cohort_sa <- dfB$branch_reform   # all within-sample; no Inf needed

# ═══════════════════════════════════════════════════════════════════════════════
cat("==============================================================\n")
cat("Panel A (49 states, 1976-2006) — Beck et al. (2010)\n")
cat("==============================================================\n\n")

# 1. TWFE — full 49-state sample
twfe_A <- feols(ln_gini ~ D_branch | statefip + wrkyr, data=dfA, cluster=~statefip)
cat(sprintf("1. TWFE:          ATT = %7.4f  SE = %.4f\n",
            coef(twfe_A)["D_branch"], se(twfe_A)["D_branch"]))

# 2. CS — 37 eligible states (not-yet-treated control, DR, group aggregation)
cs_A   <- att_gt("ln_gini","wrkyr","statefip","branch_reform",
                 data=dfA_37, control_group="notyettreated",
                 est_method="dr", bstrap=FALSE, print_details=FALSE)
cs_aggA <- aggte(cs_A, type="group", na.rm=TRUE)
cat(sprintf("2. CS:            ATT = %7.4f  SE = %.4f\n",
            cs_aggA$overall.att, cs_aggA$overall.se))

# 3. SA — full 49-state sample; Inf for always-treated
sa_A   <- feols(ln_gini ~ sunab(cohort_sa, wrkyr) | statefip + wrkyr,
                data=dfA, cluster=~statefip)
sa_aggA <- aggregate(sa_A, agg="ATT")
cat(sprintf("3. SA:            ATT = %7.4f  SE = %.4f\n",
            sa_aggA[1,"Estimate"], sa_aggA[1,"Std. Error"]))

# 4. Gardner — full 49-state sample
gard_A <- did2s(dfA, yname="ln_gini",
                first_stage=~0|statefip+wrkyr,
                second_stage=~i(D_branch),
                treatment="D_branch", cluster_var="statefip")
cat(sprintf("4. Gardner:       ATT = %7.4f  SE = %.4f\n",
            coef(gard_A)["D_branch::1"], se(gard_A)["D_branch::1"]))

# 5. Flexible TWFE — 37 eligible states
etwfe_A <- etwfe(fml=ln_gini~1, tvar=wrkyr, gvar=branch_reform,
                 data=dfA_37, vcov=~statefip)
fa <- flex_att(etwfe_A, dfA_37)
cat(sprintf("5. Flexible TWFE: ATT = %7.4f  SE = %.4f\n\n", fa$att, fa$se))

# ═══════════════════════════════════════════════════════════════════════════════
cat("==============================================================\n")
cat("Panel B (36 within-sample treated states, 1976-2006)\n")
cat("==============================================================\n\n")

# 1. TWFE
twfe_B <- feols(ln_gini ~ D_branch | statefip + wrkyr, data=dfB, cluster=~statefip)
cat(sprintf("1. TWFE:          ATT = %7.4f  SE = %.4f\n",
            coef(twfe_B)["D_branch"], se(twfe_B)["D_branch"]))

# 2. CS
cs_B    <- att_gt("ln_gini","wrkyr","statefip","branch_reform",
                  data=dfB, control_group="notyettreated",
                  est_method="dr", bstrap=FALSE, print_details=FALSE)
cs_aggB <- aggte(cs_B, type="group", na.rm=TRUE)
cat(sprintf("2. CS:            ATT = %7.4f  SE = %.4f\n",
            cs_aggB$overall.att, cs_aggB$overall.se))

# 3. SA
sa_B    <- feols(ln_gini ~ sunab(cohort_sa, wrkyr) | statefip + wrkyr,
                 data=dfB, cluster=~statefip)
sa_aggB <- aggregate(sa_B, agg="ATT")
cat(sprintf("3. SA:            ATT = %7.4f  SE = %.4f\n",
            sa_aggB[1,"Estimate"], sa_aggB[1,"Std. Error"]))

# 4. Gardner
gard_B <- did2s(dfB, yname="ln_gini",
                first_stage=~0|statefip+wrkyr,
                second_stage=~i(D_branch),
                treatment="D_branch", cluster_var="statefip")
cat(sprintf("4. Gardner:       ATT = %7.4f  SE = %.4f\n",
            coef(gard_B)["D_branch::1"], se(gard_B)["D_branch::1"]))

# 5. Flexible TWFE
etwfe_B <- etwfe(fml=ln_gini~1, tvar=wrkyr, gvar=branch_reform,
                 data=dfB, vcov=~statefip)
fb <- flex_att(etwfe_B, dfB)
cat(sprintf("5. Flexible TWFE: ATT = %7.4f  SE = %.4f\n\n", fb$att, fb$se))
