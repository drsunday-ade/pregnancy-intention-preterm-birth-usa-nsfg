# =========================
# Section 7 – Integrated Defence-Mode Summary
# Global synthesis of causal and predictive results
# =========================

library(tidyverse)
library(gt)
library(ggplot2)

# -------------------------
# 7.0 Output folder
# -------------------------
dir.create("output", showWarnings = FALSE)
dir.create("output/section7_defence", recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 7.1 Combined causal-estimator summary table
# -------------------------

summ_list <- list()

# 7.1.1 Crude and IPTW MSM (from Section 4)
if (exists("effects_logistic")) {
  eff_log <- effects_logistic %>%
    mutate(
      Population = "Full sample (survey-weighted MSM)"
    ) %>%
    select(Method, Population, Effect, Estimate, CI_low, CI_high, p_value)
  
  summ_list <- append(summ_list, list(eff_log))
}

# 7.1.2 Logistic variants under trimming / overlap (from Section 6)
if (exists("effects_ps_robust")) {
  eff_ps <- effects_ps_robust %>%
    transmute(
      Method,
      Population = Sample,
      Effect     = "Odds ratio",
      Estimate   = OR,
      CI_low,
      CI_high,
      p_value
    )
  summ_list <- append(summ_list, list(eff_ps))
}

# 7.1.3 Outcome regression and AIPW (from Section 4)
if (exists("effects_dr")) {
  eff_dr <- effects_dr %>%
    mutate(
      Population = "Full sample (survey-weighted, model-based)",
      p_value    = NA_real_
    ) %>%
    select(Method, Population, Effect, Estimate, CI_low, CI_high, p_value)
  
  summ_list <- append(summ_list, list(eff_dr))
}

# 7.1.4 TMLE + Super Learner (from Section 4, optional)
if (exists("tmle_results")) {
  eff_tmle <- tmle_results %>%
    mutate(
      Population = "Full sample (TMLE targeted)",
      p_value    = NA_real_
    ) %>%
    select(Method, Population, Effect, Estimate, CI_low, CI_high, p_value)
  
  summ_list <- append(summ_list, list(eff_tmle))
}

# 7.1.5 Causal forest (from Section 4, optional)
if (exists("cf_results")) {
  eff_cf <- cf_results %>%
    mutate(
      Population = "Overlap population (causal forest)",
      p_value    = NA_real_
    ) %>%
    select(Method, Population, Effect, Estimate, CI_low, CI_high, p_value)
  
  summ_list <- append(summ_list, list(eff_cf))
}

# 7.1.6 Bayesian g-computation (from Section 4, optional)
if (exists("effects_bayes")) {
  eff_bayes <- effects_bayes %>%
    mutate(
      Population = "Full sample (Bayesian g-computation)",
      p_value    = NA_real_
    ) %>%
    select(Method, Population, Effect, Estimate, CI_low, CI_high, p_value)
  
  summ_list <- append(summ_list, list(eff_bayes))
}

if (length(summ_list) > 0) {
  effects_all <- bind_rows(summ_list)
  
  tbl7_causal <- effects_all %>%
    arrange(Effect, Method, Population) %>%
    gt() %>%
    fmt_number(
      columns = c(Estimate, CI_low, CI_high),
      decimals = 3
    ) %>%
    fmt_number(
      columns = where(is.numeric),
      decimals = 3
    ) %>%
    cols_label(
      Method     = "Estimator",
      Population = "Target population / weighting",
      Effect     = "Effect measure",
      Estimate   = "Point estimate",
      CI_low     = "95% CI (lower)",
      CI_high    = "95% CI (upper)",
      p_value    = "p-value (if available)"
    ) %>%
    tab_header(
      title = "Global Summary of Causal Estimators",
      subtitle = "Unintended pregnancy and preterm birth across multiple causal estimands"
    )
  
  gtsave(
    tbl7_causal,
    filename = "output/section7_defence/table7a_causal_estimators_summary.html"
  )
  
} else {
  message("No causal estimators found for Section 7 summary table.")
}

# -------------------------
# 7.2 Forest plot – Odds ratios across estimators
# -------------------------

if (exists("effects_all")) {
  
  effects_or <- effects_all %>%
    filter(Effect == "Odds ratio") %>%
    filter(
      !is.na(Estimate),
      is.finite(Estimate), Estimate > 0,
      is.finite(CI_low),  CI_low  > 0,
      is.finite(CI_high), CI_high > 0
    ) %>%
    mutate(
      Label = paste0(Method, " – ", Population)
    )
  
  if (nrow(effects_or) > 0) {
    effects_or <- effects_or %>%
      arrange(Estimate) %>%
      mutate(
        Label = factor(Label, levels = Label)
      )
    
    forest7_or <- ggplot(
      effects_or,
      aes(x = Label, y = Estimate, ymin = CI_low, ymax = CI_high)
    ) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      geom_pointrange() +
      coord_flip() +
      scale_y_log10() +
      ylab("Odds ratio (log scale)") +
      xlab("") +
      ggtitle("Unintended pregnancy → preterm birth\nOdds ratios across causal estimators") +
      theme_minimal()
    
    ggsave(
      filename = "output/section7_defence/figure7a_forest_or_all_estimators.png",
      plot     = forest7_or,
      width    = 7,
      height   = 5,
      dpi      = 300
    )
  } else {
    message("No valid odds-ratio estimates for Section 7 forest plot.")
  }
  
} else {
  message("effects_all not found; OR forest plot skipped.")
}

# -------------------------
# 7.3 Forest plot – Risk differences across estimators
# -------------------------

if (exists("effects_all")) {
  
  effects_rd <- effects_all %>%
    filter(Effect == "Risk difference") %>%
    filter(
      !is.na(Estimate),
      is.finite(Estimate),
      is.finite(CI_low),
      is.finite(CI_high)
    ) %>%
    mutate(
      Label = paste0(Method, " – ", Population)
    )
  
  if (nrow(effects_rd) > 0) {
    effects_rd <- effects_rd %>%
      arrange(Estimate) %>%
      mutate(
        Label = factor(Label, levels = Label)
      )
    
    forest7_rd <- ggplot(
      effects_rd,
      aes(x = Label, y = Estimate, ymin = CI_low, ymax = CI_high)
    ) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_pointrange() +
      coord_flip() +
      ylab("Risk difference (absolute risk change)") +
      xlab("") +
      ggtitle("Unintended pregnancy → preterm birth\nRisk differences across causal estimators") +
      theme_minimal()
    
    ggsave(
      filename = "output/section7_defence/figure7b_forest_rd_all_estimators.png",
      plot     = forest7_rd,
      width    = 7,
      height   = 5,
      dpi      = 300
    )
  } else {
    message("No valid risk-difference estimates for Section 7 RD forest plot.")
  }
  
} else {
  message("effects_all not found; RD forest plot skipped.")
}

# -------------------------
# 7.4 Predictive performance recap (from Section 5)
# -------------------------

if (exists("perf_df")) {
  
  tbl7_perf <- perf_df %>%
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
      title = "Predictive Performance Recap for Preterm Birth",
      subtitle = "Weighted logistic regression vs Super Learner (if fitted)"
    )
  
  gtsave(
    tbl7_perf,
    filename = "output/section7_defence/table7b_prediction_performance_summary.html"
  )
  
} else {
  message("perf_df not found; predictive performance recap skipped.")
}

# -------------------------
# 7.5 High-risk stratification recap (from Section 5)
# -------------------------

if (exists("df_strat_out")) {
  
  tbl7_strat <- df_strat_out %>%
    gt() %>%
    fmt_number(
      columns = c(wt_sum, risk_obs, risk_pred),
      decimals = 3
    ) %>%
    cols_label(
      risk_group = "Risk group",
      n          = "Unweighted N",
      wt_sum     = "Weighted N",
      risk_obs   = "Observed preterm risk / contrast",
      risk_pred  = "Predicted risk / ratio"
    ) %>%
    tab_header(
      title = "High-Risk Stratification (Top 10% vs Lower 90%)",
      subtitle = "Observed and predicted preterm birth risk"
    )
  
  gtsave(
    tbl7_strat,
    filename = "output/section7_defence/table7c_high_risk_stratification.html"
  )
  
  # Simple plot: observed risk by risk group (excluding contrast row)
  df_strat_plot <- df_strat_out %>%
    filter(risk_group != "Contrast: Top 10% vs Lower 90%")
  
  if (nrow(df_strat_plot) == 2) {
    plot7_highrisk <- ggplot(df_strat_plot, aes(x = risk_group, y = risk_obs)) +
      geom_col() +
      xlab("") +
      ylab("Observed preterm birth risk") +
      ggtitle("Observed preterm birth risk\nTop 10% vs lower 90% predicted risk") +
      theme_minimal()
    
    ggsave(
      filename = "output/section7_defence/figure7c_high_risk_observed_risk.png",
      plot     = plot7_highrisk,
      width    = 7,
      height   = 5,
      dpi      = 300
    )
  } else {
    message("High-risk stratification plot skipped: unexpected number of risk groups.")
  }
  
} else {
  message("df_strat_out not found; high-risk stratification recap skipped.")
}
