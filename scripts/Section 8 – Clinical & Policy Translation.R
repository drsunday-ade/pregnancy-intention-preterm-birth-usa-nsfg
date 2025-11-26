# =========================
# Section 8 – Clinical & Policy Translation
# Decision Curves, Net Benefit, and Risk Thresholds
# =========================

library(tidyverse)
library(gt)
library(ggplot2)

# -------------------------
# 8.0 Output folder
# -------------------------
dir.create("output", showWarnings = FALSE)
dir.create("output/section8_translation", recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 8.1 Safety checks and primary model predictions
# -------------------------

# df_pred should exist from Section 5; if not, reconstruct minimal version from df_iptw
if (!exists("df_pred") && exists("df_iptw")) {
  df_pred <- df_iptw %>%
    mutate(
      preterm_num    = as.numeric(preterm_num),
      unintended_bin = forcats::fct_relevel(unintended_bin, "intended")
    )
}

if (!exists("df_pred")) {
  stop("df_pred not found. Please run Sections 4–5 before Section 8.")
}

# Outcome and weights
if (!("preterm_num" %in% names(df_pred))) {
  stop("'preterm_num' not found in df_pred.")
}
if (!("WGT2022_2023" %in% names(df_pred))) {
  stop("'WGT2022_2023' not found in df_pred.")
}

Y_vec <- as.numeric(df_pred$preterm_num)
w_vec <- df_pred$WGT2022_2023

# Ensure 0/1 coding
Y_vec[is.na(Y_vec)] <- 0

# Ensure we have a primary prediction vector (from Section 5); otherwise rebuild logistic
if (!exists("pred_primary") || length(pred_primary) != nrow(df_pred)) {
  # Fallback: weighted logistic regression using the same formula as Section 5
  glm8_fit <- glm(
    preterm_num ~ unintended_bin + agecon + race4 + parity_cat +
      marital_concep + educ_cat + POVERTY + curr_ins_cat + bmi_cat4,
    data    = df_pred,
    family  = binomial(),
    weights = WGT2022_2023
  )
  pred_primary <- as.numeric(predict(glm8_fit, type = "response"))
  primary_label <- "Weighted logistic regression (fallback)"
} else {
  # Use label from Section 5 if it exists; else generic label
  if (exists("primary_label")) {
    primary_label <- primary_label
  } else {
    primary_label <- "Primary prediction model"
  }
}

# Sanity: constrain predictions to (0,1)
pred_primary <- pmin(pmax(pred_primary, 1e-6), 1 - 1e-6)

# -------------------------
# 8.2 Decision-curve analysis (net benefit)
# -------------------------

# Helper function to compute net benefit at a given threshold
compute_net_benefit <- function(threshold, Y, risk, w) {
  # Treat if predicted risk >= threshold
  treat <- risk >= threshold
  
  w_sum <- sum(w)
  if (w_sum <= 0) return(c(nb_model = NA_real_, nb_all = NA_real_, nb_none = 0))
  
  # Weighted event prevalence
  prev <- sum(w * Y) / w_sum
  
  # Model-based:
  TP <- sum(w * (Y == 1 & treat))
  FP <- sum(w * (Y == 0 & treat))
  
  nb_model <- (TP / w_sum) - (FP / w_sum) * (threshold / (1 - threshold))
  
  # Treat-all strategy:
  # All are treated → TP_all = event count; FP_all = non-events
  TP_all <- sum(w * (Y == 1))
  FP_all <- sum(w * (Y == 0))
  nb_all <- (TP_all / w_sum) - (FP_all / w_sum) * (threshold / (1 - threshold))
  
  # Treat-none net benefit is 0 by definition
  c(nb_model = nb_model, nb_all = nb_all, nb_none = 0)
}

# Threshold grid (restrict to plausible clinical region around preterm risk)
thr_grid <- seq(0.05, 0.30, by = 0.01)

nb_list <- lapply(thr_grid, function(t0) {
  vals <- compute_net_benefit(t0, Y = Y_vec, risk = pred_primary, w = w_vec)
  tibble(
    threshold = t0,
    nb_model  = vals["nb_model"],
    nb_all    = vals["nb_all"],
    nb_none   = vals["nb_none"]
  )
})

df_nb <- bind_rows(nb_list)

# -------------------------
# 8.3 Decision-curve summary table
# -------------------------

# Show selected, interpretable thresholds
sel_thr <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)

df_nb_sel <- df_nb %>%
  filter(threshold %in% sel_thr)

tbl8_decision <- df_nb_sel %>%
  gt() %>%
  fmt_number(
    columns = c(threshold, nb_model, nb_all, nb_none),
    decimals = 3
  ) %>%
  cols_label(
    threshold = "Risk threshold (Pt)",
    nb_model  = paste0("Net benefit: ", primary_label),
    nb_all    = "Net benefit: treat all",
    nb_none   = "Net benefit: treat none (0)"
  ) %>%
  tab_header(
    title = paste0(primary_label, " – Decision-Curve Summary"),
    subtitle = "Net benefit for preterm birth prevention across risk thresholds"
  )

gtsave(
  tbl8_decision,
  filename = "output/section8_translation/table8a_decision_curve_summary.html"
)

# -------------------------
# 8.4 Decision-curve plot
# -------------------------

df_nb_long <- df_nb %>%
  pivot_longer(
    cols = c(nb_model, nb_all, nb_none),
    names_to = "strategy",
    values_to = "net_benefit"
  ) %>%
  mutate(
    strategy = factor(
      strategy,
      levels = c("nb_model", "nb_all", "nb_none"),
      labels = c(primary_label, "Treat all", "Treat none")
    )
  )

dec_curve_plot <- ggplot(df_nb_long, aes(x = threshold, y = net_benefit, colour = strategy)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(size = 1) +
  xlab("Risk threshold (Pt) for intervention") +
  ylab("Net benefit (events averted per patient)") +
  ggtitle(
    paste0(
      primary_label,
      " – Decision Curve for Preterm Birth Risk"
    )
  ) +
  scale_colour_discrete(name = "Strategy") +
  theme_minimal()

ggsave(
  filename = "output/section8_translation/figure8a_decision_curve.png",
  plot     = dec_curve_plot,
  width    = 7,
  height   = 5,
  dpi      = 300
)

# -------------------------
# 8.5 Clinical classification at a reference threshold (Pt = 0.10)
# -------------------------

pt_ref <- 0.10

df_class <- df_pred %>%
  mutate(
    pred_primary = pred_primary,
    treat_flag   = pred_primary >= pt_ref,
    Y            = Y_vec,
    w            = w_vec
  )

# Weighted confusion matrix
TP <- sum(df_class$w * (df_class$Y == 1 & df_class$treat_flag))
FP <- sum(df_class$w * (df_class$Y == 0 & df_class$treat_flag))
FN <- sum(df_class$w * (df_class$Y == 1 & !df_class$treat_flag))
TN <- sum(df_class$w * (df_class$Y == 0 & !df_class$treat_flag))
w_total <- sum(df_class$w)

prev <- (TP + FN) / w_total

sens <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)   # recall
spec <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
ppv  <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
npv  <- ifelse((TN + FN) > 0, TN / (TN + FN), NA_real_)

df_class_summary <- tibble(
  Quantity   = c(
    "Weighted prevalence of preterm birth",
    "Threshold (Pt) for intervention",
    "Sensitivity (recall)",
    "Specificity",
    "Positive predictive value (PPV)",
    "Negative predictive value (NPV)",
    "Weighted TP count",
    "Weighted FP count",
    "Weighted FN count",
    "Weighted TN count"
  ),
  Estimate   = c(
    prev,
    pt_ref,
    sens,
    spec,
    ppv,
    npv,
    TP,
    FP,
    FN,
    TN
  )
)

tbl8_class <- df_class_summary %>%
  gt() %>%
  fmt_number(
    columns = Estimate,
    decimals = 3
  ) %>%
  cols_label(
    Quantity = "Measure",
    Estimate = "Estimate (weighted)"
  ) %>%
  tab_header(
    title = paste0(primary_label, " – Clinical Classification at Pt = 0.10"),
    subtitle = "Weighted performance for targeting high-risk pregnancies"
  )

gtsave(
  tbl8_class,
  filename = "output/section8_translation/table8b_classification_at_0_10.html"
)
