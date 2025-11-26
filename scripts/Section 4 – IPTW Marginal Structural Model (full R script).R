# =========================
# Section 4 – Advanced Causal Inference
# IPTW MSM + AIPW + G-computation
# + TMLE + Super Learner
# + Causal Forests (GRF)
# + Bayesian g-computation
# =========================

library(survey)
library(tidyverse)
library(broom)
library(gt)
library(ggplot2)

# -------------------------
# 4.0 Output folder
# -------------------------
dir.create("output", showWarnings = FALSE)
dir.create("output/section4_causal", recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 4.1 Build IPTW analysis dataset
#    (complete cases for exposure, outcome, confounders)
# -------------------------

vars_iptw <- c(
  "preterm",          # outcome (factor)
  "unintended_bin",   # exposure (factor: intended / unintended)
  "agecon",           # numeric age at conception
  "race4",            # race/ethnicity
  "parity_cat",       # parity
  "marital_concep",   # marital/cohab at conception
  "educ_cat",         # education
  "POVERTY",          # numeric poverty ratio
  "curr_ins_cat",     # insurance
  "bmi_cat4",         # BMI
  "WGT2022_2023",     # survey weight
  "VEST", "VECL"      # strata / PSUs
)

df_iptw <- nsfg_preg_analytic %>%
  dplyr::select(all_of(vars_iptw)) %>%
  # drop rows with any missing in key variables
  filter(
    !is.na(preterm),
    !is.na(unintended_bin),
    !is.na(agecon),
    !is.na(race4),
    !is.na(parity_cat),
    !is.na(marital_concep),
    !is.na(educ_cat),
    !is.na(POVERTY),
    !is.na(curr_ins_cat),
    !is.na(bmi_cat4),
    !is.na(WGT2022_2023),
    !is.na(VEST),
    !is.na(VECL)
  ) %>%
  mutate(
    # numeric outcome for logistic MSM / DR
    preterm_num    = as.numeric(preterm == "preterm"),
    # ensure exposure reference is "intended"
    unintended_bin = forcats::fct_relevel(unintended_bin, "intended")
  )

# Base survey design (for crude models + balance)
design_base <- svydesign(
  ids     = ~VECL,
  strata  = ~VEST,
  weights = ~WGT2022_2023,
  nest    = TRUE,
  data    = df_iptw
)

# -------------------------
# 4.2 Propensity score model for unintended pregnancy
# -------------------------

ps_model <- glm(
  unintended_bin ~ agecon + race4 + parity_cat +
    marital_concep + educ_cat + POVERTY +
    curr_ins_cat + bmi_cat4,
  data    = df_iptw,
  family  = binomial(),
  weights = WGT2022_2023
)

# Predicted probability of unintended pregnancy (bounded away from 0/1)
df_iptw <- df_iptw %>%
  mutate(
    ps_raw = predict(ps_model, type = "response"),
    ps     = pmin(pmax(ps_raw, 0.01), 0.99)   # numerical stability
  )

# Stabilized weights numerator: marginal P(A=1) and P(A=0), weighted
pA1 <- with(df_iptw,
            sum(WGT2022_2023 * (unintended_bin == "unintended")) /
              sum(WGT2022_2023))
pA0 <- 1 - pA1

df_iptw <- df_iptw %>%
  mutate(
    sw = case_when(
      unintended_bin == "unintended" ~ pA1 / ps,
      unintended_bin == "intended"   ~ pA0 / (1 - ps),
      TRUE ~ NA_real_
    )
  )

# Truncate stabilized weights at 1st and 99th percentiles
q_sw <- quantile(df_iptw$sw, probs = c(0.01, 0.99), na.rm = TRUE)

df_iptw <- df_iptw %>%
  mutate(
    sw_trunc = pmin(pmax(sw, q_sw[1]), q_sw[2]),
    w_comb   = WGT2022_2023 * sw_trunc   # combined survey + IPTW
  )

# Diagnostics table for weights
w_diag <- df_iptw %>%
  summarise(
    n                = n(),
    sw_min           = min(sw, na.rm = TRUE),
    sw_p1            = quantile(sw, 0.01, na.rm = TRUE),
    sw_p50           = median(sw, na.rm = TRUE),
    sw_p99           = quantile(sw, 0.99, na.rm = TRUE),
    sw_max           = max(sw, na.rm = TRUE),
    sw_trunc_min     = min(sw_trunc, na.rm = TRUE),
    sw_trunc_max     = max(sw_trunc, na.rm = TRUE)
  )

tbl4_weights <- w_diag %>%
  gt() %>%
  fmt_number(
    columns = where(is.numeric),
    decimals = 3
  ) %>%
  tab_header(
    title = "Stabilized IPTW Weights: Distribution Before and After Truncation"
  )

gtsave(
  tbl4_weights,
  filename = "output/section4_causal/table4a_weights_diagnostics.html"
)

# Histogram of truncated stabilized weights
weights_hist <- ggplot(df_iptw, aes(x = sw_trunc)) +
  geom_histogram(bins = 40, boundary = 0, closed = "left") +
  xlab("Stabilized IPTW (truncated)") +
  ylab("Count") +
  ggtitle("Distribution of Truncated Stabilized IPTW Weights") +
  theme_minimal()

ggsave(
  filename = "output/section4_causal/figure4a_weights_histogram.png",
  plot     = weights_hist,
  width    = 7,
  height   = 5,
  dpi      = 300
)

# -------------------------
# 4.3 IPTW survey design and MSM
# -------------------------

design_iptw <- svydesign(
  ids     = ~VECL,
  strata  = ~VEST,
  weights = ~w_comb,
  nest    = TRUE,
  data    = df_iptw
)

# Crude survey-weighted model (no IPTW)
msm_crude <- svyglm(
  preterm_num ~ unintended_bin,
  design = design_base,
  family = quasibinomial()
)

# IPTW MSM (causal effect model)
msm_iptw <- svyglm(
  preterm_num ~ unintended_bin,
  design = design_iptw,
  family = quasibinomial()
)

# -------------------------
# 4.4 Effect estimates: crude vs IPTW MSM (OR)
# -------------------------

tidy_crude <- tidy(msm_crude, conf.int = TRUE, exponentiate = TRUE)
tidy_iptw  <- tidy(msm_iptw,  conf.int = TRUE, exponentiate = TRUE)

est_crude <- tidy_crude %>%
  filter(term == "unintended_binunintended") %>%
  mutate(
    Method = "Crude logistic (survey-weighted)",
    Effect = "Odds ratio"
  )

est_iptw <- tidy_iptw %>%
  filter(term == "unintended_binunintended") %>%
  mutate(
    Method = "IPTW MSM (survey-weighted)",
    Effect = "Odds ratio"
  )

effects_logistic <- bind_rows(est_crude, est_iptw) %>%
  transmute(
    Method,
    Effect,
    Estimate = estimate,
    CI_low   = conf.low,
    CI_high  = conf.high,
    p_value  = p.value
  )

tbl4_msm <- effects_logistic %>%
  gt() %>%
  fmt_number(
    columns = c(Estimate, CI_low, CI_high),
    decimals = 2
  ) %>%
  fmt_number(
    columns = c(p_value),
    decimals = 3
  ) %>%
  cols_label(
    Estimate = "Odds ratio (unintended vs intended)",
    CI_low   = "95% CI (lower)",
    CI_high  = "95% CI (upper)",
    p_value  = "p-value"
  ) %>%
  tab_header(
    title = "Unintended Pregnancy and Preterm Birth:",
    subtitle = "Crude vs IPTW Marginal Structural Model (Odds Ratios)"
  )

gtsave(
  tbl4_msm,
  filename = "output/section4_causal/table4b_msm_effects.html"
)

# -------------------------
# 4.5 Covariate balance before vs after IPTW (categorical)
# -------------------------

covars_cat <- c(
  "race4",
  "parity_cat",
  "marital_concep",
  "educ_cat",
  "curr_ins_cat",
  "bmi_cat4"
)

balance_one_cat <- function(var_name, design_obj, design_label) {
  f <- as.formula(paste0("~", var_name, "+unintended_bin"))
  tab <- svytable(f, design_obj)
  df <- as.data.frame(tab)
  names(df)[1:2] <- c("level", "unintended_bin")
  
  df %>%
    group_by(unintended_bin) %>%
    mutate(
      pct = Freq / sum(Freq) * 100
    ) %>%
    ungroup() %>%
    mutate(
      variable  = var_name,
      weighting = design_label
    )
}

bal_cat_base <- bind_rows(
  lapply(covars_cat, balance_one_cat,
         design_obj = design_base,
         design_label = "Baseline")
)

bal_cat_iptw <- bind_rows(
  lapply(covars_cat, balance_one_cat,
         design_obj = design_iptw,
         design_label = "IPTW")
)

bal_cat_all <- bind_rows(bal_cat_base, bal_cat_iptw)

tbl4_balance_cat <- bal_cat_all %>%
  mutate(
    Variable = variable,
    Level    = as.character(level),
    Group    = as.character(unintended_bin)
  ) %>%
  select(Variable, Level, Group, weighting, pct) %>%
  arrange(Variable, Level, Group, weighting) %>%
  gt() %>%
  fmt_number(columns = pct, decimals = 1) %>%
  cols_label(
    Variable  = "Covariate",
    Level     = "Level",
    Group     = "Pregnancy intention",
    weighting = "Weighting scheme",
    pct       = "Weighted % within intention group"
  ) %>%
  tab_header(
    title = "Covariate Balance Before and After IPTW (Categorical Confounders)"
  )

gtsave(
  tbl4_balance_cat,
  filename = "output/section4_causal/table4c_balance_categorical.html"
)

# -------------------------
# 4.6 Covariate balance (continuous: age, poverty)
# -------------------------

covars_cont <- c("agecon", "POVERTY")

balance_one_cont <- function(var_name, design_obj, design_label) {
  f <- as.formula(paste0("~", var_name))
  # means by exposure group
  m <- svyby(f, ~unintended_bin, design_obj, svymean, keep.var = FALSE)
  v <- svyby(f, ~unintended_bin, design_obj, svyvar,  keep.var = FALSE)
  
  # Identify numeric mean and var columns generically
  mean_col <- setdiff(names(m),
                      c("unintended_bin", grep("^se", names(m), value = TRUE)))
  var_col  <- setdiff(names(v),
                      c("unintended_bin", grep("^se", names(v), value = TRUE)))
  
  out_df <- data.frame(
    variable       = var_name,
    weighting      = design_label,
    unintended_bin = m$unintended_bin,
    mean           = as.numeric(m[[mean_col]]),
    var            = as.numeric(v[[var_col]])
  )
  out_df$sd <- sqrt(out_df$var)
  out_df
}

bal_cont_base <- do.call(
  rbind,
  lapply(covars_cont, balance_one_cont,
         design_obj = design_base,
         design_label = "Baseline")
)

bal_cont_iptw <- do.call(
  rbind,
  lapply(covars_cont, balance_one_cont,
         design_obj = design_iptw,
         design_label = "IPTW")
)

bal_cont_all <- rbind(bal_cont_base, bal_cont_iptw)

tbl4_balance_cont <- bal_cont_all %>%
  as_tibble() %>%
  mutate(
    Variable = variable,
    Group    = as.character(unintended_bin)
  ) %>%
  select(Variable, Group, weighting, mean, sd) %>%
  arrange(Variable, Group, weighting) %>%
  gt() %>%
  fmt_number(columns = c(mean, sd), decimals = 2) %>%
  cols_label(
    Variable  = "Covariate",
    Group     = "Pregnancy intention",
    weighting = "Weighting scheme",
    mean      = "Weighted mean",
    sd        = "Weighted SD"
  ) %>%
  tab_header(
    title = "Covariate Balance Before and After IPTW (Continuous Confounders)"
  )

gtsave(
  tbl4_balance_cont,
  filename = "output/section4_causal/table4d_balance_continuous.html"
)

# =====================================================
# 4.7 Advanced: Doubly Robust AIPW + G-computation
# =====================================================

# Helper: compute g-computation and AIPW estimates for a given dataset
compute_estimates <- function(df) {
  # df must contain: preterm_num, unintended_bin, WGT2022_2023,
  # agecon, race4, parity_cat, marital_concep, educ_cat, POVERTY, curr_ins_cat, bmi_cat4
  
  df <- df %>%
    mutate(
      A = ifelse(unintended_bin == "unintended", 1, 0),
      Y = preterm_num
    )
  
  w <- df$WGT2022_2023
  Y <- df$Y
  
  # Propensity model for A
  ps_fit <- glm(
    A ~ agecon + race4 + parity_cat +
      marital_concep + educ_cat + POVERTY +
      curr_ins_cat + bmi_cat4,
    family  = binomial(),
    data    = df,
    weights = w
  )
  ps <- as.numeric(predict(ps_fit, type = "response"))
  ps <- pmin(pmax(ps, 0.01), 0.99)
  
  # Outcome model for Y
  out_fit <- glm(
    Y ~ A + agecon + race4 + parity_cat +
      marital_concep + educ_cat + POVERTY +
      curr_ins_cat + bmi_cat4,
    family  = binomial(),
    data    = df,
    weights = w
  )
  
  # g-computation: predict Y under A=1 and A=0 for all
  df1 <- df; df1$A <- 1
  df0 <- df; df0$A <- 0
  
  m1 <- as.numeric(predict(out_fit, newdata = df1, type = "response"))
  m0 <- as.numeric(predict(out_fit, newdata = df0, type = "response"))
  
  w_sum <- sum(w)
  
  gcomp_psi1 <- sum(w * m1) / w_sum
  gcomp_psi0 <- sum(w * m0) / w_sum
  
  # AIPW / DR
  aipw1 <- m1 + (df$A * (Y - m1)) / ps
  aipw0 <- m0 + ((1 - df$A) * (Y - m0)) / (1 - ps)
  
  aipw_psi1 <- sum(w * aipw1) / w_sum
  aipw_psi0 <- sum(w * aipw0) / w_sum
  
  c(
    gcomp_psi1 = gcomp_psi1,
    gcomp_psi0 = gcomp_psi0,
    aipw_psi1  = aipw_psi1,
    aipw_psi0  = aipw_psi0
  )
}

# Point estimates on full data
est_point <- compute_estimates(df_iptw)

gcomp_psi1 <- est_point["gcomp_psi1"]
gcomp_psi0 <- est_point["gcomp_psi0"]
aipw_psi1  <- est_point["aipw_psi1"]
aipw_psi0  <- est_point["aipw_psi0"]

# Risk differences, risk ratios, odds ratios
gcomp_RD <- gcomp_psi1 - gcomp_psi0
gcomp_RR <- gcomp_psi1 / gcomp_psi0
gcomp_OR <- (gcomp_psi1 / (1 - gcomp_psi1)) /
  (gcomp_psi0 / (1 - gcomp_psi0))

aipw_RD <- aipw_psi1 - aipw_psi0
aipw_RR <- aipw_psi1 / aipw_psi0
aipw_OR <- (aipw_psi1 / (1 - aipw_psi1)) /
  (aipw_psi0 / (1 - aipw_psi0))

# Bootstrap for DR & g-computation CIs
set.seed(2025)
B <- 200
n <- nrow(df_iptw)

boot_mat <- matrix(NA_real_, nrow = B, ncol = 4)
colnames(boot_mat) <- c("gcomp_psi1","gcomp_psi0","aipw_psi1","aipw_psi0")

count <- 0
for (b in 1:B) {
  idx  <- sample.int(n, n, replace = TRUE)
  df_b <- df_iptw[idx, ]
  res_b <- try(compute_estimates(df_b), silent = TRUE)
  if (inherits(res_b, "try-error")) next
  if (any(is.na(res_b))) next
  count <- count + 1
  if (count <= B) {
    boot_mat[count, ] <- res_b
  }
}
if (count == 0) {
  warning("All bootstrap replicates failed; no CIs computed for DR/g-computation.")
} else {
  boot_mat <- boot_mat[1:count, , drop = FALSE]
}

boot_gcomp_RD <- boot_mat[,"gcomp_psi1"] - boot_mat[,"gcomp_psi0"]
boot_gcomp_RR <- boot_mat[,"gcomp_psi1"] / boot_mat[,"gcomp_psi0"]
boot_gcomp_OR <- (boot_mat[,"gcomp_psi1"] / (1 - boot_mat[,"gcomp_psi1"])) /
  (boot_mat[,"gcomp_psi0"] / (1 - boot_mat[,"gcomp_psi0"]))

boot_aipw_RD <- boot_mat[,"aipw_psi1"] - boot_mat[,"aipw_psi0"]
boot_aipw_RR <- boot_mat[,"aipw_psi1"] / boot_mat[,"aipw_psi0"]
boot_aipw_OR <- (boot_mat[,"aipw_psi1"] / (1 - boot_mat[,"aipw_psi1"])) /
  (boot_mat[,"aipw_psi0"] / (1 - boot_mat[,"aipw_psi0"]))

q <- c(0.025, 0.975)

gcomp_RD_CI <- quantile(boot_gcomp_RD, probs = q, na.rm = TRUE)
gcomp_RR_CI <- quantile(boot_gcomp_RR, probs = q, na.rm = TRUE)
gcomp_OR_CI <- quantile(boot_gcomp_OR, probs = q, na.rm = TRUE)

aipw_RD_CI <- quantile(boot_aipw_RD, probs = q, na.rm = TRUE)
aipw_RR_CI <- quantile(boot_aipw_RR, probs = q, na.rm = TRUE)
aipw_OR_CI <- quantile(boot_aipw_OR, probs = q, na.rm = TRUE)

effects_dr <- tribble(
  ~Method,                      ~Effect,          ~Estimate, ~CI_low,             ~CI_high,
  "Outcome regression (g-comp)", "Risk difference", gcomp_RD, gcomp_RD_CI[1],    gcomp_RD_CI[2],
  "Outcome regression (g-comp)", "Risk ratio",      gcomp_RR, gcomp_RR_CI[1],    gcomp_RR_CI[2],
  "Outcome regression (g-comp)", "Odds ratio",      gcomp_OR, gcomp_OR_CI[1],    gcomp_OR_CI[2],
  "AIPW doubly robust",         "Risk difference", aipw_RD,  aipw_RD_CI[1],     aipw_RD_CI[2],
  "AIPW doubly robust",         "Risk ratio",      aipw_RR,  aipw_RR_CI[1],     aipw_RR_CI[2],
  "AIPW doubly robust",         "Odds ratio",      aipw_OR,  aipw_OR_CI[1],     aipw_OR_CI[2]
)

tbl4_dr <- effects_dr %>%
  gt() %>%
  fmt_number(
    columns = c(Estimate, CI_low, CI_high),
    decimals = 3
  ) %>%
  cols_label(
    Method   = "Estimator",
    Effect   = "Effect measure",
    Estimate = "Point estimate",
    CI_low   = "95% CI (lower)",
    CI_high  = "95% CI (upper)"
  ) %>%
  tab_header(
    title = "Unintended Pregnancy → Preterm Birth:",
    subtitle = "Advanced Causal Estimators (G-computation and Doubly Robust AIPW)"
  )

gtsave(
  tbl4_dr,
  filename = "output/section4_causal/table4e_dr_effects.html"
)

# 4.7.1 OR forest plot combining main causal estimators
effects_all_or <- bind_rows(
  effects_logistic,                      # crude & IPTW ORs
  effects_dr %>% filter(Effect == "Odds ratio")
) %>%
  filter(Effect == "Odds ratio") %>%
  filter(
    !is.na(Estimate),
    is.finite(Estimate), is.finite(CI_low), is.finite(CI_high),
    Estimate > 0, CI_low > 0, CI_high > 0
  ) %>%
  mutate(
    Method = factor(
      Method,
      levels = c(
        "Crude logistic (survey-weighted)",
        "IPTW MSM (survey-weighted)",
        "Outcome regression (g-comp)",
        "AIPW doubly robust"
      )
    )
  )

forest_or <- ggplot(
  effects_all_or,
  aes(x = Method, y = Estimate, ymin = CI_low, ymax = CI_high)
) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_pointrange() +
  coord_flip() +
  scale_y_log10() +
  ylab("Odds ratio (log scale)") +
  xlab("") +
  ggtitle("Effect of Unintended Pregnancy on Preterm Birth\nAcross Causal Estimators") +
  theme_minimal()

ggsave(
  filename = "output/section4_causal/figure4b_forest_or_all_methods.png",
  plot     = forest_or,
  width    = 7,
  height   = 5,
  dpi      = 300
)

# =====================================================
# 4.8 TMLE (Targeted Maximum Likelihood Estimation)
#     with Super Learner for PS and outcome
# =====================================================

if (requireNamespace("SuperLearner", quietly = TRUE) &&
    requireNamespace("tmle", quietly = TRUE)) {
  
  library(SuperLearner)
  library(tmle)
  
  df_tmle <- df_iptw %>%
    mutate(
      A = ifelse(unintended_bin == "unintended", 1, 0),
      Y = preterm_num
    )
  
  # Covariate matrix W
  W_tmle <- df_tmle %>%
    dplyr::select(
      agecon, race4, parity_cat, marital_concep,
      educ_cat, POVERTY, curr_ins_cat, bmi_cat4
    )
  
  # Minimal but safe Super Learner library
  SL.lib <- c("SL.mean", "SL.glm")
  
  tmle_fit <- tmle::tmle(
    Y            = df_tmle$Y,
    A            = df_tmle$A,
    W            = W_tmle,
    family       = "binomial",
    Q.SL.library = SL.lib,
    g.SL.library = SL.lib,
    obsWeights   = df_tmle$WGT2022_2023
  )
  
  tmle_RD <- as.numeric(tmle_fit$estimates$ATE$psi)
  tmle_CI <- as.numeric(tmle_fit$estimates$ATE$CI)
  
  tmle_results <- tibble(
    Method   = "TMLE + Super Learner",
    Effect   = "Risk difference",
    Estimate = tmle_RD,
    CI_low   = tmle_CI[1],
    CI_high  = tmle_CI[2]
  )
  
  tbl4_tmle <- tmle_results %>%
    gt() %>%
    fmt_number(columns = c(Estimate, CI_low, CI_high), decimals = 3) %>%
    cols_label(
      Method   = "Estimator",
      Effect   = "Effect measure",
      Estimate = "Point estimate",
      CI_low   = "95% CI (lower)",
      CI_high  = "95% CI (upper)"
    ) %>%
    tab_header(
      title = "TMLE with Super Learner",
      subtitle = "Average causal effect of unintended pregnancy on preterm birth (Risk difference)"
    )
  
  gtsave(
    tbl4_tmle,
    filename = "output/section4_causal/table4f_tmle_superlearner.html"
  )
  
} else {
  message("TMLE section skipped: packages 'SuperLearner' and/or 'tmle' not installed.")
}

# =====================================================
# 4.9 Causal Forests / Generalized Random Forests
# =====================================================

if (requireNamespace("grf", quietly = TRUE)) {
  
  library(grf)
  
  df_cf <- df_iptw %>%
    mutate(
      A = ifelse(unintended_bin == "unintended", 1, 0),
      Y = preterm_num
    )
  
  X_cf <- model.matrix(
    ~ agecon + race4 + parity_cat + marital_concep + educ_cat +
      POVERTY + curr_ins_cat + bmi_cat4,
    data = df_cf
  )[ , -1, drop = FALSE]   # drop intercept
  
  cf <- grf::causal_forest(
    X              = X_cf,
    Y              = df_cf$Y,
    W              = df_cf$A,
    W.hat          = df_cf$ps,           # PS from earlier model
    Y.hat          = NULL,               # let forest learn outcome
    sample.weights = df_cf$WGT2022_2023,
    num.trees      = 2000
  )
  
  ate_cf <- grf::average_treatment_effect(cf, target.sample = "all")
  cf_RD  <- as.numeric(ate_cf["estimate"])
  cf_se  <- as.numeric(ate_cf["std.err"])
  cf_CI  <- cf_RD + c(-1.96, 1.96) * cf_se
  
  cf_results <- tibble(
    Method   = "Causal forest (GRF)",
    Effect   = "Risk difference",
    Estimate = cf_RD,
    CI_low   = cf_CI[1],
    CI_high  = cf_CI[2]
  )
  
  tbl4_cf <- cf_results %>%
    gt() %>%
    fmt_number(columns = c(Estimate, CI_low, CI_high), decimals = 3) %>%
    cols_label(
      Method   = "Estimator",
      Effect   = "Effect measure",
      Estimate = "Point estimate",
      CI_low   = "95% CI (lower)",
      CI_high  = "95% CI (upper)"
    ) %>%
    tab_header(
      title = "Causal Forest (GRF)",
      subtitle = "Average treatment effect of unintended pregnancy on preterm birth (Risk difference)"
    )
  
  gtsave(
    tbl4_cf,
    filename = "output/section4_causal/table4g_causal_forest.html"
  )
  
  # Distribution of individual CATEs
  cate_hat <- predict(cf)$predictions
  
  df_cate <- tibble(
    tau_hat = cate_hat
  )
  
  cate_hist <- ggplot(df_cate, aes(x = tau_hat)) +
    geom_histogram(bins = 40, boundary = 0, closed = "left") +
    xlab("Estimated individual treatment effect (risk difference)") +
    ylab("Count") +
    ggtitle("Distribution of Individual Causal Effects (Causal Forest)") +
    theme_minimal()
  
  ggsave(
    filename = "output/section4_causal/figure4c_cate_histogram_causal_forest.png",
    plot     = cate_hist,
    width    = 7,
    height   = 5,
    dpi      = 300
  )
  
} else {
  message("Causal forest section skipped: package 'grf' not installed.")
}

# =====================================================
# 4.10 Bayesian g-computation with informative priors
# =====================================================

if (requireNamespace("rstanarm", quietly = TRUE)) {
  
  library(rstanarm)
  options(mc.cores = max(1L, parallel::detectCores() - 1L))
  # If rstan is available, enable compiled model caching
  if (requireNamespace("rstan", quietly = TRUE)) {
    rstan::rstan_options(auto_write = TRUE)
  }
  
  df_bayes <- df_iptw %>%
    mutate(
      A = ifelse(unintended_bin == "unintended", 1, 0),
      Y = preterm_num
    )
  
  bayes_fit <- rstanarm::stan_glm(
    Y ~ A + agecon + race4 + parity_cat +
      marital_concep + educ_cat + POVERTY +
      curr_ins_cat + bmi_cat4,
    data   = df_bayes,
    family = binomial(link = "logit"),
    weights = WGT2022_2023,
    prior           = rstanarm::normal(0, 2.5, autoscale = TRUE),
    prior_intercept = rstanarm::normal(0, 5, autoscale = TRUE),
    chains = 2, iter = 2000, seed = 2025, refresh = 0
  )
  
  df1 <- df_bayes; df1$A <- 1
  df0 <- df_bayes; df0$A <- 0
  
  post_p1 <- rstanarm::posterior_epred(bayes_fit, newdata = df1)
  post_p0 <- rstanarm::posterior_epred(bayes_fit, newdata = df0)
  
  w <- df_bayes$WGT2022_2023
  w <- w / sum(w)
  
  psi1_draws <- as.numeric(post_p1 %*% w)
  psi0_draws <- as.numeric(post_p0 %*% w)
  
  eps <- 1e-4
  psi1_draws <- pmin(pmax(psi1_draws, eps), 1 - eps)
  psi0_draws <- pmin(pmax(psi0_draws, eps), 1 - eps)
  
  RD_draws <- psi1_draws - psi0_draws
  RR_draws <- psi1_draws / psi0_draws
  OR_draws <- (psi1_draws / (1 - psi1_draws)) /
    (psi0_draws / (1 - psi0_draws))
  
  bayes_RD <- mean(RD_draws)
  bayes_RR <- mean(RR_draws)
  bayes_OR <- mean(OR_draws)
  
  bayes_RD_CI <- quantile(RD_draws, probs = c(0.025, 0.975))
  bayes_RR_CI <- quantile(RR_draws, probs = c(0.025, 0.975))
  bayes_OR_CI <- quantile(OR_draws, probs = c(0.025, 0.975))
  
  effects_bayes <- tribble(
    ~Method,                    ~Effect,          ~Estimate, ~CI_low,          ~CI_high,
    "Bayesian g-computation",   "Risk difference", bayes_RD, bayes_RD_CI[1],  bayes_RD_CI[2],
    "Bayesian g-computation",   "Risk ratio",      bayes_RR, bayes_RR_CI[1],  bayes_RR_CI[2],
    "Bayesian g-computation",   "Odds ratio",      bayes_OR, bayes_OR_CI[1],  bayes_OR_CI[2]
  )
  
  tbl4_bayes <- effects_bayes %>%
    gt() %>%
    fmt_number(
      columns = c(Estimate, CI_low, CI_high),
      decimals = 3
    ) %>%
    cols_label(
      Method   = "Estimator",
      Effect   = "Effect measure",
      Estimate = "Posterior mean",
      CI_low   = "95% credible interval (lower)",
      CI_high  = "95% credible interval (upper)"
    ) %>%
    tab_header(
      title = "Bayesian g-computation",
      subtitle = "Average causal effect of unintended pregnancy on preterm birth"
    )
  
  gtsave(
    tbl4_bayes,
    filename = "output/section4_causal/table4h_bayesian_gcomp.html"
  )
  
  df_RD <- tibble(RD = RD_draws)
  rd_density <- ggplot(df_RD, aes(x = RD)) +
    geom_density() +
    xlab("Risk difference (unintended vs intended)") +
    ylab("Posterior density") +
    ggtitle("Bayesian g-computation: Posterior of Risk Difference") +
    theme_minimal()
  
  ggsave(
    filename = "output/section4_causal/figure4d_bayes_RD_posterior.png",
    plot     = rd_density,
    width    = 7,
    height   = 5,
    dpi      = 300
  )
  
} else {
  message("Bayesian g-computation section skipped: package 'rstanarm' not installed.")
}
