# =========================
# Section 6 – Causal Robustness & Sensitivity Analyses
# Positivity, Trimming, Overlap Weighting, E-values
# =========================

library(survey)
library(tidyverse)
library(broom)
library(gt)
library(ggplot2)

# -------------------------
# 6.0 Output folder
# -------------------------
dir.create("output", showWarnings = FALSE)
dir.create("output/section6_robustness", recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 6.1 Safety checks
# -------------------------

if (!exists("df_iptw")) {
  stop("df_iptw not found. Please run Sections 1–4 before Section 6.")
}
if (!("ps" %in% names(df_iptw))) {
  stop("Propensity scores 'ps' not found in df_iptw. Make sure Section 4 ran successfully.")
}
if (!("w_comb" %in% names(df_iptw))) {
  stop("Combined weights 'w_comb' not found in df_iptw. Make sure Section 4 (IPTW) ran successfully.")
}
if (!("preterm_num" %in% names(df_iptw))) {
  stop("'preterm_num' not found in df_iptw. Please ensure outcome coding from previous sections.")
}

# Ensure factors are well-behaved
df_iptw <- df_iptw %>%
  mutate(
    unintended_bin = forcats::fct_relevel(unintended_bin, "intended")
  )

# Base and IPTW survey designs
design_base <- svydesign(
  ids     = ~VECL,
  strata  = ~VEST,
  weights = ~WGT2022_2023,
  nest    = TRUE,
  data    = df_iptw
)

design_iptw <- svydesign(
  ids     = ~VECL,
  strata  = ~VEST,
  weights = ~w_comb,
  nest    = TRUE,
  data    = df_iptw
)

# -------------------------
# 6.2 Propensity-score (PS) diagnostics
# -------------------------

# 6.2.1 PS summary by pregnancy intention
ps_summary <- df_iptw %>%
  group_by(unintended_bin) %>%
  summarise(
    n        = n(),
    wt_sum   = sum(WGT2022_2023),
    ps_min   = min(ps, na.rm = TRUE),
    ps_p25   = quantile(ps, 0.25, na.rm = TRUE),
    ps_med   = median(ps, na.rm = TRUE),
    ps_p75   = quantile(ps, 0.75, na.rm = TRUE),
    ps_max   = max(ps, na.rm = TRUE),
    .groups  = "drop"
  )

tbl6_ps <- ps_summary %>%
  gt() %>%
  fmt_number(
    columns = c(wt_sum, ps_min, ps_p25, ps_med, ps_p75, ps_max),
    decimals = 3
  ) %>%
  cols_label(
    unintended_bin = "Pregnancy intention",
    n              = "Unweighted N",
    wt_sum         = "Weighted N",
    ps_min         = "PS min",
    ps_p25         = "PS 25th pct",
    ps_med         = "PS median",
    ps_p75         = "PS 75th pct",
    ps_max         = "PS max"
  ) %>%
  tab_header(
    title = "Propensity Score Distribution by Pregnancy Intention",
    subtitle = "Unweighted and weighted summaries"
  )

gtsave(
  tbl6_ps,
  filename = "output/section6_robustness/table6a_ps_summary.html"
)

# 6.2.2 Overlap plot: PS histogram by intention
ps_hist <- ggplot(df_iptw, aes(x = ps, fill = unintended_bin)) +
  geom_histogram(
    position = "identity",
    alpha    = 0.4,
    bins     = 30
  ) +
  xlab("Propensity score P(unintended pregnancy | covariates)") +
  ylab("Count") +
  ggtitle("Propensity Score Overlap by Pregnancy Intention") +
  scale_fill_discrete(name = "Pregnancy intention") +
  theme_minimal()

ggsave(
  filename = "output/section6_robustness/figure6a_ps_hist_overlap.png",
  plot     = ps_hist,
  width    = 7,
  height   = 5,
  dpi      = 300
)

# 6.2.3 Overlap plot: PS density by intention
ps_dens <- ggplot(df_iptw, aes(x = ps, colour = unintended_bin)) +
  geom_density(adjust = 1.2, na.rm = TRUE) +
  xlab("Propensity score P(unintended pregnancy | covariates)") +
  ylab("Density") +
  ggtitle("Propensity Score Density by Pregnancy Intention") +
  scale_colour_discrete(name = "Pregnancy intention") +
  theme_minimal()

ggsave(
  filename = "output/section6_robustness/figure6b_ps_density_overlap.png",
  plot     = ps_dens,
  width    = 7,
  height   = 5,
  dpi      = 300
)

# 6.2.4 PS vs truncated IPTW weights (to show extremes)
df_psw <- df_iptw %>%
  filter(!is.na(sw_trunc))

ps_weight_scatter <- ggplot(df_psw, aes(x = ps, y = sw_trunc)) +
  geom_point(alpha = 0.4) +
  xlab("Propensity score") +
  ylab("Truncated stabilized IPTW") +
  ggtitle("Relation Between Propensity Scores and IPTW Weights") +
  theme_minimal()

ggsave(
  filename = "output/section6_robustness/figure6c_ps_vs_truncIPTW.png",
  plot     = ps_weight_scatter,
  width    = 7,
  height   = 5,
  dpi      = 300
)

# -------------------------
# 6.3 Trimmed sample (Crump 0.1–0.9) and overlap weighting
# -------------------------

# 6.3.1 PS trimming per Crump et al. (0.1–0.9)
df_trim <- df_iptw %>%
  filter(ps >= 0.10, ps <= 0.90)

have_trim <- nrow(df_trim) > 0L

if (have_trim) {
  design_trim_base <- svydesign(
    ids     = ~VECL,
    strata  = ~VEST,
    weights = ~WGT2022_2023,
    nest    = TRUE,
    data    = df_trim
  )
  
  design_trim_iptw <- svydesign(
    ids     = ~VECL,
    strata  = ~VEST,
    weights = ~w_comb,
    nest    = TRUE,
    data    = df_trim
  )
}

# 6.3.2 Overlap weights (Li, Morgan, Zaslavsky 2018)
# A = 1 for unintended, 0 for intended
df_iptw <- df_iptw %>%
  mutate(
    A          = ifelse(unintended_bin == "unintended", 1, 0),
    w_overlap  = WGT2022_2023 * if_else(A == 1, 1 - ps, ps)
  )

design_ow <- svydesign(
  ids     = ~VECL,
  strata  = ~VEST,
  weights = ~w_overlap,
  nest    = TRUE,
  data    = df_iptw
)

# -------------------------
# 6.4 Causal effect re-estimation under variants
# -------------------------

# Helper: extract OR for unintended vs intended
extract_or <- function(fit, method_label, sample_label) {
  tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "unintended_binunintended") %>%
    transmute(
      Method  = method_label,
      Sample  = sample_label,
      OR      = estimate,
      CI_low  = conf.low,
      CI_high = conf.high,
      p_value = p.value
    )
}

# 6.4.1 Full-sample crude and IPTW MSM
msm_crude_full <- svyglm(
  preterm_num ~ unintended_bin,
  design = design_base,
  family = quasibinomial()
)

msm_iptw_full <- svyglm(
  preterm_num ~ unintended_bin,
  design = design_iptw,
  family = quasibinomial()
)

est_crude_full <- extract_or(
  msm_crude_full,
  method_label = "Crude logistic (survey-weighted)",
  sample_label = "Full sample"
)

est_iptw_full <- extract_or(
  msm_iptw_full,
  method_label = "IPTW MSM (survey-weighted)",
  sample_label = "Full sample"
)

# 6.4.2 Trimmed-sample models (if available)
if (have_trim) {
  msm_crude_trim <- svyglm(
    preterm_num ~ unintended_bin,
    design = design_trim_base,
    family = quasibinomial()
  )
  
  msm_iptw_trim <- svyglm(
    preterm_num ~ unintended_bin,
    design = design_trim_iptw,
    family = quasibinomial()
  )
  
  est_crude_trim <- extract_or(
    msm_crude_trim,
    method_label = "Crude logistic (survey-weighted)",
    sample_label = "PS-trimmed (0.10–0.90)"
  )
  
  est_iptw_trim <- extract_or(
    msm_iptw_trim,
    method_label = "IPTW MSM (survey-weighted)",
    sample_label = "PS-trimmed (0.10–0.90)"
  )
}

# 6.4.3 Overlap-weighted effect
msm_ow <- svyglm(
  preterm_num ~ unintended_bin,
  design = design_ow,
  family = quasibinomial()
)

est_ow <- extract_or(
  msm_ow,
  method_label = "Overlap weighting (OW)",
  sample_label = "Full overlap sample"
)

# 6.4.4 Combine effects into robustness table
effects_list <- list(
  est_crude_full,
  est_iptw_full,
  est_ow
)

if (have_trim) {
  effects_list <- append(
    effects_list,
    list(est_crude_trim, est_iptw_trim)
  )
}

effects_ps_robust <- bind_rows(effects_list)

tbl6_causal <- effects_ps_robust %>%
  gt() %>%
  fmt_number(
    columns = c(OR, CI_low, CI_high),
    decimals = 3
  ) %>%
  fmt_number(
    columns = p_value,
    decimals = 3
  ) %>%
  cols_label(
    Method  = "Estimator",
    Sample  = "Sample / PS support",
    OR      = "Odds ratio (unintended vs intended)",
    CI_low  = "95% CI (lower)",
    CI_high = "95% CI (upper)",
    p_value = "p-value"
  ) %>%
  tab_header(
    title = "Causal Effect of Unintended Pregnancy on Preterm Birth",
    subtitle = "Robustness across crude, IPTW MSM, trimming, and overlap weighting"
  )

gtsave(
  tbl6_causal,
  filename = "output/section6_robustness/table6b_causal_robustness_logistic.html"
)

# -------------------------
# 6.5 E-values for unmeasured confounding
#     (manual VanderWeele formula; no external package)
# -------------------------

# Helper: E-value for a risk ratio (or OR treated as RR)
evalue_rr <- function(rr) {
  if (!is.finite(rr) || rr <= 0) return(NA_real_)
  if (rr < 1) rr <- 1 / rr
  if (rr <= 1) return(1)
  rr + sqrt(rr * (rr - 1))
}

# Helper: choose CI bound to use for E-value (VanderWeele)
evalue_ci_rr <- function(est, lo, hi) {
  # Work on RR / OR scale; convert to "worse than null" side
  if (!all(is.finite(c(est, lo, hi)))) return(NA_real_)
  if (est >= 1) {
    # Harmful effect; use lower bound if >1
    if (is.na(lo) || lo <= 1) return(NA_real_)
    rr_ci <- lo
  } else {
    # Protective; invert
    if (is.na(hi) || hi >= 1) return(NA_real_)
    rr_ci <- 1 / hi
  }
  evalue_rr(rr_ci)
}

# 1) Choose primary OR – prefer AIPW from Section 4; else IPTW MSM
if (exists("aipw_OR") && exists("aipw_OR_CI")) {
  est_val <- as.numeric(aipw_OR)
  lo_val  <- as.numeric(aipw_OR_CI[1])
  hi_val  <- as.numeric(aipw_OR_CI[2])
  target_label <- "AIPW doubly robust OR"
} else {
  # Fallback: IPTW MSM OR from full sample
  est_row <- est_iptw_full[1, ]
  est_val <- as.numeric(est_row$OR)
  lo_val  <- as.numeric(est_row$CI_low)
  hi_val  <- as.numeric(est_row$CI_high)
  target_label <- "IPTW MSM OR"
}

if (all(is.finite(c(est_val, lo_val, hi_val))) && est_val > 0) {
  
  # Treat OR as approximate RR for E-values
  e_point <- evalue_rr(est_val)
  e_ci    <- evalue_ci_rr(est_val, lo_val, hi_val)
  
  # For transparency, also report which CI bound was used on OR scale
  ci_used_or <- if (est_val >= 1) lo_val else 1 / hi_val
  
  e_df <- tibble(
    Quantity        = c("Point estimate", "Confidence interval bound"),
    OR_used         = c(est_val, ci_used_or),
    E_value         = c(e_point, e_ci)
  )
  
  tbl6_evalues <- e_df %>%
    gt() %>%
    fmt_number(
      columns = c(OR_used, E_value),
      decimals = 3
    ) %>%
    cols_label(
      Quantity = "Quantity",
      OR_used  = "OR used (treated as RR)",
      E_value  = "E-value"
    ) %>%
    tab_header(
      title = "E-values for Unmeasured Confounding",
      subtitle = paste0("Based on ", target_label, " (OR treated as RR)")
    )
  
  gtsave(
    tbl6_evalues,
    filename = "output/section6_robustness/table6c_evalues_unmeasured_confounding.html"
  )
  
} else {
  message("E-value computation skipped: invalid or non-positive effect/CI values.")
}
