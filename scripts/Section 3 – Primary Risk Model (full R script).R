# =========================
# Section 3 – Advanced Risk Model
# =========================

library(survey)
library(tidyverse)
library(broom)
library(gt)
library(ggplot2)
library(pROC)      # for AUC / ROC

# -------------------------
# 3.0 Survey design object
# (skip if already created, but safe to re-run)
# -------------------------
nsfg_design <- svydesign(
  ids     = ~VECL,
  strata  = ~VEST,
  weights = ~WGT2022_2023,
  nest    = TRUE,
  data    = nsfg_preg_analytic
)

# -------------------------
# 3.0 Output folder
# -------------------------
dir.create("output", showWarnings = FALSE)
dir.create("output/section3_riskmodel", recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 3.1 Advanced survey-weighted logistic regression
# - splines for age and poverty
# - interaction: unintended × race4
# -------------------------

model_risk_adv <- svyglm(
  preterm ~ unintended_bin * race4 +
    splines::ns(agecon, df = 4) +           # flexible age effect
    parity_cat +
    marital_concep +
    educ_cat +
    splines::ns(POVERTY, df = 4) +          # flexible poverty effect
    curr_ins_cat +
    bmi_cat4,
  design = nsfg_design,
  family = quasibinomial()
)

summary(model_risk_adv)

# -------------------------
# 3.2 Regression OR table (HTML) – focus on interpretable terms
#    (hide spline coefficients)
# -------------------------

tidy_risk_adv <- tidy(
  model_risk_adv,
  conf.int = TRUE,
  exponentiate = TRUE
) %>%
  filter(term != "(Intercept)") %>%
  # drop spline basis terms (they're nuisance, not directly interpretable)
  filter(!grepl("agecon", term),
         !grepl("POVERTY", term)) %>%
  mutate(
    term_label = case_when(
      term == "unintended_binunintended" ~
        "Unintended vs intended (pregnancy intention)",
      grepl("^race4", term) ~
        paste("Race/ethnicity:",
              gsub("race4", "", term)),
      grepl("^parity_cat", term) ~
        paste("Parity:", gsub("parity_cat", "", term)),
      grepl("^marital_concep", term) ~
        paste("Marital status at conception:",
              gsub("marital_concep", "", term)),
      grepl("^educ_cat", term) ~
        paste("Education:", gsub("educ_cat", "", term)),
      grepl("^curr_ins_cat", term) ~
        paste("Insurance:", gsub("curr_ins_cat", "", term)),
      grepl("^bmi_cat4", term) ~
        paste("BMI category:", gsub("bmi_cat4", "", term)),
      grepl("unintended_binunintended:race4", term) ~
        paste("Interaction: unintended × race",
              gsub("unintended_binunintended:race4", "", term)),
      TRUE ~ term
    )
  )

tbl3_risk_adv <- tidy_risk_adv %>%
  transmute(
    Predictor = term_label,
    OR        = estimate,
    CI_low    = conf.low,
    CI_high   = conf.high,
    p_value   = p.value
  ) %>%
  gt() %>%
  fmt_number(
    columns = c(OR, CI_low, CI_high),
    decimals = 2
  ) %>%
  fmt_number(
    columns = c(p_value),
    decimals = 3
  ) %>%
  cols_label(
    CI_low  = "95% CI (lower)",
    CI_high = "95% CI (upper)",
    p_value = "p-value"
  ) %>%
  tab_header(
    title = "Advanced Survey-weighted Logistic Regression:",
    subtitle = "Preterm birth vs pregnancy intention, race interaction, and covariates"
  )

gtsave(
  tbl3_risk_adv,
  filename = "output/section3_riskmodel/table3a_logistic_risk_model_advanced.html"
)

# -------------------------
# 3.3 Forest plot of adjusted ORs (PNG)
# -------------------------

forest_df <- tidy_risk_adv %>%
  filter(!is.na(estimate)) %>%
  arrange(estimate) %>%
  mutate(
    term_label = factor(term_label, levels = term_label)
  )

forest_plot_adv <- ggplot(
  forest_df,
  aes(x = term_label, y = estimate,
      ymin = conf.low, ymax = conf.high)
) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_pointrange() +
  coord_flip() +
  scale_y_log10() +
  ylab("Adjusted Odds Ratio (log scale)") +
  xlab("") +
  ggtitle("Adjusted Odds Ratios for Preterm Birth (Advanced Model)") +
  theme_minimal()

ggsave(
  filename = "output/section3_riskmodel/figure3a_forest_or_advanced.png",
  plot     = forest_plot_adv,
  width    = 7,
  height   = 6,
  dpi      = 300
)

# -------------------------
# Use ONLY complete-case data actually used in the model
# for AUC, calibration, and marginal predictions
# -------------------------

df_model  <- model.frame(model_risk_adv)      # complete cases used in fit
w_model   <- weights(model_risk_adv)          # survey weights for those cases
pred_adv  <- predict(model_risk_adv, type = "response")
preterm_num <- as.numeric(df_model$preterm == "preterm")

# sanity check
stopifnot(
  length(pred_adv) == nrow(df_model),
  length(preterm_num) == nrow(df_model),
  length(w_model) == nrow(df_model)
)

# -------------------------
# 3.4 Model performance: AUC / ROC (PNG + HTML summary)
# -------------------------

roc_obj <- pROC::roc(
  response  = preterm_num,
  predictor = pred_adv,
  weights   = w_model,
  quiet     = TRUE
)

auc_val <- as.numeric(pROC::auc(roc_obj))

# Build ROC curve data frame
roc_df <- tibble(
  tpr = roc_obj$sensitivities,
  fpr = 1 - roc_obj$specificities
)

roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() +
  xlab("False positive rate (1 - specificity)") +
  ylab("True positive rate (sensitivity)") +
  ggtitle(paste0("ROC Curve – Advanced Model (AUC = ",
                 sprintf("%.3f", auc_val), ")")) +
  theme_minimal()

ggsave(
  filename = "output/section3_riskmodel/figure3b_roc_curve_advanced.png",
  plot     = roc_plot,
  width    = 6,
  height   = 6,
  dpi      = 300
)

# AUC summary table
tbl3_auc <- tibble(
  Metric = "AUC (area under ROC curve)",
  Value  = auc_val
) %>%
  gt() %>%
  fmt_number(columns = Value, decimals = 3) %>%
  tab_header(
    title = "Discrimination Performance – Advanced Model"
  )

gtsave(
  tbl3_auc,
  filename = "output/section3_riskmodel/table3b_auc_advanced.html"
)

# -------------------------
# 3.5 Calibration plot (PNG)
#    – survey-weighted observed vs predicted, by decile
# -------------------------

df_perf <- tibble(
  preterm_num = preterm_num,
  pred_adv    = pred_adv,
  w_svy       = w_model
) %>%
  mutate(decile = ntile(pred_adv, 10))

calib_df <- df_perf %>%
  group_by(decile) %>%
  summarise(
    mean_pred = sum(pred_adv * w_svy, na.rm = TRUE) /
      sum(w_svy, na.rm = TRUE),
    obs_risk  = sum(preterm_num * w_svy, na.rm = TRUE) /
      sum(w_svy, na.rm = TRUE),
    .groups = "drop"
  )

calib_plot <- ggplot(calib_df, aes(x = mean_pred, y = obs_risk)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  geom_line() +
  xlab("Mean predicted risk (per decile)") +
  ylab("Observed preterm risk (survey-weighted)") +
  ggtitle("Calibration Plot – Advanced Model") +
  theme_minimal()

ggsave(
  filename = "output/section3_riskmodel/figure3c_calibration_plot_advanced.png",
  plot     = calib_plot,
  width    = 6,
  height   = 6,
  dpi      = 300
)

# -------------------------
# 3.6 Adjusted predicted probabilities (intended vs unintended)
#    – marginal standardization using advanced model
#    – using same complete-case data + weights
# -------------------------

X_adv   <- model.matrix(model_risk_adv)   # design matrix for df_model
beta_adv <- coef(model_risk_adv)

# Identify column for unintended_binunintended
col_unintended <- grep("^unintended_binunintended$", colnames(X_adv))

# Scenario 1: set everyone to intended (0)
X_intended <- X_adv
if (length(col_unintended) == 1) {
  X_intended[, col_unintended] <- 0
}
eta_intended <- drop(X_intended %*% beta_adv)
p_intended   <- plogis(eta_intended)
p_intended_pop <- sum(p_intended * w_model, na.rm = TRUE) / sum(w_model, na.rm = TRUE)

# Scenario 2: set everyone to unintended (1)
X_unintended <- X_adv
if (length(col_unintended) == 1) {
  X_unintended[, col_unintended] <- 1
}
eta_unintended <- drop(X_unintended %*% beta_adv)
p_unintended   <- plogis(eta_unintended)
p_unintended_pop <- sum(p_unintended * w_model, na.rm = TRUE) / sum(w_model, na.rm = TRUE)

adj_pred_df_adv <- tibble(
  intention    = c("Intended pregnancy", "Unintended pregnancy"),
  preterm_prob = c(p_intended_pop, p_unintended_pop)
) %>%
  mutate(
    preterm_pct = 100 * preterm_prob,
    diff_pct    = preterm_pct - first(preterm_pct)
  )

tbl3_pred_adv <- adj_pred_df_adv %>%
  gt() %>%
  fmt_number(
    columns = c(preterm_prob, preterm_pct, diff_pct),
    decimals = 3
  ) %>%
  cols_label(
    intention    = "Pregnancy intention (counterfactual scenario)",
    preterm_prob = "Predicted probability of preterm birth",
    preterm_pct  = "Predicted preterm birth (%)",
    diff_pct     = "Difference vs intended (%)"
  ) %>%
  tab_header(
    title = "Adjusted Predicted Probability of Preterm Birth by Pregnancy Intention",
    subtitle = "Advanced survey-weighted logistic model with splines and interaction"
  )

gtsave(
  tbl3_pred_adv,
  filename = "output/section3_riskmodel/table3c_adjusted_predicted_probabilities_advanced.html"
)
