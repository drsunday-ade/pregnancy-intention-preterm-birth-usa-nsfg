# =========================
# Section 5 – Advanced Risk Prediction & Stratification
# Super Learner + Logistic Baseline
# ROC, Calibration, and High-Risk Targeting
# =========================

library(tidyverse)
library(gt)
library(ggplot2)

# -------------------------
# 5.0 Output folder
# -------------------------
dir.create("output", showWarnings = FALSE)
dir.create("output/section5_prediction", recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 5.1 Prediction dataset
#     (reuse df_iptw from Section 4)
# -------------------------

# df_iptw already contains:
#   preterm_num, unintended_bin, agecon, race4, parity_cat,
#   marital_concep, educ_cat, POVERTY, curr_ins_cat, bmi_cat4,
#   WGT2022_2023, ps, etc.

df_pred <- df_iptw %>%
  mutate(
    preterm_num    = as.numeric(preterm_num),
    unintended_bin = forcats::fct_relevel(unintended_bin, "intended")
  )

Y_vec <- df_pred$preterm_num
w_vec <- df_pred$WGT2022_2023

# Design matrix for prediction (includes intention + covariates)
pred_formula <- ~ unintended_bin + agecon + race4 + parity_cat +
  marital_concep + educ_cat + POVERTY + curr_ins_cat + bmi_cat4

X_mat <- model.matrix(pred_formula, data = df_pred)[, -1, drop = FALSE]  # drop intercept
X_df  <- as.data.frame(X_mat)

# -------------------------
# 5.2 Baseline weighted logistic model
# -------------------------

glm_fit <- glm(
  preterm_num ~ unintended_bin + agecon + race4 + parity_cat +
    marital_concep + educ_cat + POVERTY + curr_ins_cat + bmi_cat4,
  data    = df_pred,
  family  = binomial(),
  weights = WGT2022_2023
)

pred_glm <- as.numeric(predict(glm_fit, type = "response"))

# Bound predictions away from 0/1 (numerical stability for later metrics)
pred_glm <- pmin(pmax(pred_glm, 1e-4), 1 - 1e-4)

# Weighted Brier score
brier_glm <- mean(w_vec * (Y_vec - pred_glm)^2) / mean(w_vec)

# -------------------------
# 5.3 Super Learner ensemble (if available)
# -------------------------

sl_available <- requireNamespace("SuperLearner", quietly = TRUE)

pred_sl  <- NULL
brier_sl <- NA_real_
auc_glm  <- NA_real_
auc_sl   <- NA_real_

if (sl_available) {
  library(SuperLearner)
  
  # Basic but extensible library; add more learners if installed
  SL.lib <- c("SL.mean", "SL.glm")
  if (requireNamespace("randomForest", quietly = TRUE)) {
    SL.lib <- c(SL.lib, "SL.randomForest")
  }
  if (requireNamespace("xgboost", quietly = TRUE)) {
    SL.lib <- c(SL.lib, "SL.xgboost")
  }
  
  # Rescale weights for numerical stability
  w_sl <- w_vec / mean(w_vec)
  
  sl_fit <- SuperLearner::SuperLearner(
    Y          = Y_vec,
    X          = X_df,
    family     = binomial(),
    SL.library = SL.lib,
    obsWeights = w_sl,
    cvControl  = list(V = 10, stratifyCV = TRUE)
  )
  
  # Cross-validated ensemble predictions
  pred_sl <- as.numeric(sl_fit$SL.predict)
  pred_sl <- pmin(pmax(pred_sl, 1e-4), 1 - 1e-4)
  
  # Weighted Brier
  brier_sl <- mean(w_vec * (Y_vec - pred_sl)^2) / mean(w_vec)
}

# -------------------------
# 5.4 AUC (if pROC installed)
# -------------------------

if (requireNamespace("pROC", quietly = TRUE)) {
  library(pROC)
  
  # Logistic ROC / AUC
  roc_glm <- pROC::roc(
    response = Y_vec,
    predictor = pred_glm,
    weights = w_vec,
    quiet = TRUE
  )
  auc_glm <- as.numeric(pROC::auc(roc_glm))
  
  # SuperLearner ROC / AUC (if fitted)
  if (!is.null(pred_sl)) {
    roc_sl <- pROC::roc(
      response = Y_vec,
      predictor = pred_sl,
      weights = w_vec,
      quiet = TRUE
    )
    auc_sl <- as.numeric(pROC::auc(roc_sl))
  }
}

# -------------------------
# 5.5 Performance summary table (AUC + Brier)
# -------------------------

perf_list <- list(
  tibble(
    Model      = "Weighted logistic regression",
    AUC        = auc_glm,
    Brier      = brier_glm
  )
)

if (!is.null(pred_sl)) {
  perf_list <- append(
    perf_list,
    list(
      tibble(
        Model = "Super Learner ensemble",
        AUC   = auc_sl,
        Brier = brier_sl
      )
    )
  )
}

perf_df <- bind_rows(perf_list)

tbl5_perf <- perf_df %>%
  gt() %>%
  fmt_number(
    columns = c(AUC, Brier),
    decimals = 3
  ) %>%
  cols_label(
    Model = "Prediction model",
    AUC   = "AUC (weighted)",
    Brier = "Brier score (weighted)"
  ) %>%
  tab_header(
    title = "Predictive Performance for Preterm Birth",
    subtitle = "Weighted logistic regression vs Super Learner ensemble"
  )

gtsave(
  tbl5_perf,
  filename = "output/section5_prediction/table5a_model_performance.html"
)

# -------------------------
# 5.6 Calibration by risk deciles
# -------------------------

# Choose primary model for calibration: Super Learner if available, else logistic
if (!is.null(pred_sl)) {
  pred_primary   <- pred_sl
  primary_label  <- "Super Learner ensemble"
} else {
  pred_primary   <- pred_glm
  primary_label  <- "Weighted logistic regression"
}

df_calib <- df_pred %>%
  mutate(
    pred_primary = pred_primary,
    decile       = dplyr::ntile(pred_primary, 10)
  ) %>%
  group_by(decile) %>%
  summarise(
    n         = n(),
    wt_sum    = sum(WGT2022_2023),
    pred_mean = sum(WGT2022_2023 * pred_primary) / wt_sum,
    obs_mean  = sum(WGT2022_2023 * preterm_num) / wt_sum,
    .groups   = "drop"
  ) %>%
  mutate(
    abs_diff = obs_mean - pred_mean,
    ratio    = ifelse(pred_mean > 0, obs_mean / pred_mean, NA_real_)
  )

tbl5_calib <- df_calib %>%
  mutate(
    Decile = decile
  ) %>%
  select(
    Decile,
    n,
    wt_sum,
    pred_mean,
    obs_mean,
    abs_diff,
    ratio
  ) %>%
  gt() %>%
  fmt_number(
    columns = c(wt_sum, pred_mean, obs_mean, abs_diff, ratio),
    decimals = 3
  ) %>%
  cols_label(
    Decile    = "Risk decile (lowest → highest)",
    n         = "Unweighted N",
    wt_sum    = "Weighted N",
    pred_mean = "Mean predicted risk",
    obs_mean  = "Observed risk",
    abs_diff  = "Obs − Pred",
    ratio     = "Obs / Pred"
  ) %>%
  tab_header(
    title = paste0(primary_label, ": Calibration by Risk Decile"),
    subtitle = "Weighted observed vs predicted preterm birth risk"
  )

gtsave(
  tbl5_calib,
  filename = "output/section5_prediction/table5b_calibration_deciles.html"
)

# -------------------------
# 5.7 Calibration plot (decile-level)
# -------------------------

calib_plot <- ggplot(df_calib, aes(x = pred_mean, y = obs_mean)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  geom_line() +
  xlim(0, max(df_calib$pred_mean, df_calib$obs_mean, na.rm = TRUE) * 1.05) +
  ylim(0, max(df_calib$pred_mean, df_calib$obs_mean, na.rm = TRUE) * 1.05) +
  xlab("Mean predicted risk of preterm birth") +
  ylab("Weighted observed risk of preterm birth") +
  ggtitle(
    paste0(
      primary_label,
      ": Calibration Plot (Decile-Level)"
    )
  ) +
  theme_minimal()

ggsave(
  filename = "output/section5_prediction/figure5b_calibration_plot.png",
  plot     = calib_plot,
  width    = 7,
  height   = 5,
  dpi      = 300
)

# -------------------------
# 5.8 ROC curve plot (if pROC available)
# -------------------------

if (requireNamespace("pROC", quietly = TRUE)) {
  library(pROC)
  
  # ROC for primary model
  roc_primary <- pROC::roc(
    response  = Y_vec,
    predictor = pred_primary,
    weights   = w_vec,
    quiet     = TRUE
  )
  
  roc_df <- tibble(
    tpr = rev(roc_primary$sensitivities),
    fpr = rev(1 - roc_primary$specificities)
  )
  
  roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    xlab("False positive rate (1 − specificity)") +
    ylab("True positive rate (sensitivity)") +
    ggtitle(
      paste0(
        primary_label,
        ": ROC Curve for Preterm Birth Prediction"
      )
    ) +
    theme_minimal()
  
  ggsave(
    filename = "output/section5_prediction/figure5a_roc_primary_model.png",
    plot     = roc_plot,
    width    = 7,
    height   = 5,
    dpi      = 300
  )
}

# -------------------------
# 5.9 High-risk stratification (top 10% predicted risk)
# -------------------------

thr_90 <- quantile(pred_primary, 0.90, na.rm = TRUE)

df_strat <- df_pred %>%
  mutate(
    pred_primary = pred_primary,
    risk_group = dplyr::if_else(
      pred_primary >= thr_90,
      "Top 10% predicted risk",
      "Lower 90%"
    )
  ) %>%
  group_by(risk_group) %>%
  summarise(
    n         = n(),
    wt_sum    = sum(WGT2022_2023),
    risk_obs  = sum(WGT2022_2023 * preterm_num) / wt_sum,
    risk_pred = sum(WGT2022_2023 * pred_primary) / wt_sum,
    .groups   = "drop"
  )

# Compute contrast (high vs low)
if (nrow(df_strat) == 2) {
  df_high <- df_strat %>% filter(risk_group == "Top 10% predicted risk")
  df_low  <- df_strat %>% filter(risk_group == "Lower 90%")
  
  RD_high_low <- df_high$risk_obs - df_low$risk_obs
  RR_high_low <- df_high$risk_obs / df_low$risk_obs
  
  contrast_row <- tibble(
    risk_group = "Contrast: Top 10% vs Lower 90%",
    n          = NA_integer_,
    wt_sum     = NA_real_,
    risk_obs   = RD_high_low,
    risk_pred  = RR_high_low
  )
  
  df_strat_out <- bind_rows(df_strat, contrast_row)
} else {
  df_strat_out <- df_strat
}

tbl5_strat <- df_strat_out %>%
  gt() %>%
  fmt_number(
    columns = c(wt_sum, risk_obs, risk_pred),
    decimals = 3
  ) %>%
  cols_label(
    risk_group = "Risk group",
    n          = "Unweighted N",
    wt_sum     = "Weighted N",
    risk_obs   = "Observed preterm risk",
    risk_pred  = "Mean predicted risk / contrast"
  ) %>%
  tab_header(
    title = paste0(primary_label, ": High-Risk Stratification"),
    subtitle = "Top 10% of predicted risk vs lower 90%"
  ) %>%
  tab_footnote(
    footnote = "For the contrast row, 'Observed preterm risk' is the risk difference (Top 10% − Lower 90%), and 'Mean predicted risk / contrast' is the risk ratio (Top 10% / Lower 90%).",
    locations = cells_body(
      rows    = risk_group == "Contrast: Top 10% vs Lower 90%",
      columns = c(risk_obs, risk_pred)
    )
  )

gtsave(
  tbl5_strat,
  filename = "output/section5_prediction/table5c_high_risk_stratification.html"
)

# -------------------------
# 5.10 Optional: Random forest variable importance (if available)
# -------------------------

if (requireNamespace("randomForest", quietly = TRUE)) {
  library(randomForest)
  
  rf_fit <- randomForest::randomForest(
    x       = X_df,
    y       = factor(Y_vec, levels = c(0, 1)),
    ntree   = 1000,
    sampsize = pmin(5000, nrow(X_df))
  )
  
  vip <- as.data.frame(randomForest::importance(rf_fit, type = 2))
  vip$Variable <- rownames(vip)
  rownames(vip) <- NULL
  
  vip_tbl <- vip %>%
    arrange(desc(MeanDecreaseGini)) %>%
    mutate(
      Rank = row_number()
    ) %>%
    select(Rank, Variable, MeanDecreaseGini) %>%
    gt() %>%
    fmt_number(
      columns = MeanDecreaseGini,
      decimals = 2
    ) %>%
    cols_label(
      Rank             = "Rank",
      Variable         = "Predictor",
      MeanDecreaseGini = "Random forest importance (MeanDecreaseGini)"
    ) %>%
    tab_header(
      title = "Variable Importance (Random Forest)",
      subtitle = "Exploratory ranking of predictors for preterm birth"
    )
  
  gtsave(
    vip_tbl,
    filename = "output/section5_prediction/table5d_variable_importance_rf.html"
  )
} else {
  message("Random forest variable importance skipped: package 'randomForest' not installed.")
}
