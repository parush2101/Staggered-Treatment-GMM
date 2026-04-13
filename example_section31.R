###############################################################################
# Section 3.1 Example: Three Cohorts, Three Time Periods
# Arora & Bijani (2026), Figure 1 and Equation (7)
#
# Cohort-time mean outcomes (Eq. 7):
#   Never-treated (g = 0): {Y_1, Y_2, Y_3} = {100, 100, 100}
#   Early-treated (g = 2): {Y_1, Y_2, Y_3} = {110, 130, 125}
#   Late-treated  (g = 3): {Y_1, Y_2, Y_3} = {140, 140, 165}
#
# True CATTs under parallel trends (verified on page 7):
#   beta_22 = 20  (g=2, t=2: impact period)
#   beta_23 = 15  (g=2, t=3: one period after impact)
#   beta_33 = 25  (g=3, t=3: impact period)
#
# DGP: Y_igt = mu_g + alpha_i + tau_gt + epsilon_it
#   mu_g      : cohort baseline intercept (100 / 110 / 140)
#   alpha_i   : unit fixed effect ~ N(0, 1)
#   tau_gt    : true CATT (0 in pre-treatment periods)
#   epsilon_it: idiosyncratic error ~ N(0, 1)
#
# Cohort size: 100 units each (300 total).
# Estimation: Flexible TWFE (Wooldridge 2025).
###############################################################################

library(data.table)
library(fixest)

set.seed(42)

# ===========================================================================
# 1. Parameters
# ===========================================================================

cohort_size <- 100
T_total     <- 3
N_total     <- 3L * cohort_size   # 300

# Cohort membership (g): 0 = never-treated, 2 = treated from t=2, 3 = from t=3
g_cohort <- c(rep(0L, cohort_size), rep(2L, cohort_size), rep(3L, cohort_size))

# Cohort baseline intercepts (match pre-treatment means in Eq. 7)
mu_g <- c("0" = 100, "2" = 110, "3" = 140)

# True CATTs from the paper
true_catt <- c(beta_22 = 20, beta_23 = 15, beta_33 = 25)
true_att  <- mean(true_catt)   # equal-weighted ATT = 20

# ===========================================================================
# 2. Generate Data
# ===========================================================================

unit_id <- rep(seq_len(N_total), each  = T_total)
time_id <- rep(seq_len(T_total), times = N_total)
g_id    <- g_cohort[unit_id]

# Unit fixed effects: alpha_i ~ N(0, 1), one draw per unit
alpha_i <- rnorm(N_total, mean = 0, sd = 1)

# Idiosyncratic errors: epsilon_it ~ N(0, 1)
epsilon  <- rnorm(N_total * T_total, mean = 0, sd = 1)

# Cohort baseline for each observation
baseline <- mu_g[as.character(g_id)]

# Treatment indicator: 1 if unit's cohort has received treatment by time t
D_it <- as.integer((g_id == 2L & time_id >= 2L) |
                   (g_id == 3L & time_id >= 3L))

# True treatment effect tau_gt (= CATT at post-treatment periods, 0 elsewhere)
tau_it <- numeric(N_total * T_total)
tau_it[g_id == 2L & time_id == 2L] <- true_catt["beta_22"]   # 20
tau_it[g_id == 2L & time_id == 3L] <- true_catt["beta_23"]   # 15
tau_it[g_id == 3L & time_id == 3L] <- true_catt["beta_33"]   # 25

# Outcome
Y_it <- baseline + alpha_i[unit_id] + tau_it + epsilon

dt <- data.table(
  unit = unit_id,
  time = time_id,
  g    = g_id,
  D    = D_it,
  tau  = tau_it,
  Y    = Y_it
)

# ===========================================================================
# 3. Verify DGP: cohort-time means should approximate Eq. (7)
# ===========================================================================

means_check <- dcast(
  dt[, .(mean_Y = round(mean(Y), 3)), by = .(g, time)],
  g ~ time, value.var = "mean_Y"
)
setnames(means_check, as.character(1:T_total), paste0("t=", 1:T_total))

cat("=======================================================\n")
cat("  DGP Verification: Cohort-Time Means\n")
cat("=======================================================\n")
cat("  Expected: {100,100,100}, {110,130,125}, {140,140,165}\n\n")
print(means_check)
cat("\n")

# ===========================================================================
# 4. Flexible TWFE Estimation (Wooldridge 2025)
#
#   Y_igt = alpha_i + lambda_t
#           + sum_{g,k>=0} theta_{g,g+k} * 1(G_i=g, t=g+k) + epsilon_it
#
#   Implemented via cohort-time interaction dummies on treated obs.
# ===========================================================================

# Cohort-time cell label: non-empty only for treated (g,t) cells
dt[, cohort_time := fifelse(D == 1L,
                            paste0("g", g, "_t", time),
                            "ref")]

mod_flex <- feols(Y ~ i(cohort_time, ref = "ref") | unit + time, data = dt)

# Extract and label CATT estimates
coef_raw  <- coef(mod_flex)
coef_keys <- gsub("cohort_time::", "", names(coef_raw))
names(coef_raw) <- coef_keys

# Map estimator output keys to CATT names
catt_keys <- c("g2_t2", "g2_t3", "g3_t3")
catt_names <- c("beta_22", "beta_23", "beta_33")
est_catt  <- coef_raw[catt_keys]

# ===========================================================================
# 5. Results Table
# ===========================================================================

cat("=======================================================\n")
cat("  Flexible TWFE: CATT Estimates vs True Values\n")
cat("=======================================================\n")
cat(sprintf("%-10s  %10s  %10s  %10s\n",
            "CATT", "Estimate", "True", "Bias"))
cat(paste(rep("-", 46), collapse = ""), "\n")

for (i in seq_along(catt_keys)) {
  est  <- est_catt[catt_keys[i]]
  tval <- true_catt[catt_names[i]]
  cat(sprintf("%-10s  %10.4f  %10.4f  %10.4f\n",
              catt_keys[i], est, tval, est - tval))
}

cat(paste(rep("-", 46), collapse = ""), "\n")
att_est <- mean(est_catt)
cat(sprintf("%-10s  %10.4f  %10.4f  %10.4f\n",
            "ATT", att_est, true_att, att_est - true_att))
cat("\n")

cat("=======================================================\n")
cat("  Full Model Output\n")
cat("=======================================================\n")
print(summary(mod_flex))
