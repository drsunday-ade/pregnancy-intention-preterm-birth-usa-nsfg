# =========================
# Section 2 – Descriptive & Unadjusted Analyses
# =========================

library(survey)
library(tidyverse)
library(gt)

# -------------------------
# 2.0 Survey design object
# -------------------------
nsfg_design <- svydesign(
  ids     = ~VECL,
  strata  = ~VEST,
  weights = ~WGT2022_2023,
  nest    = TRUE,
  data    = nsfg_preg_analytic
)

# -------------------------
# 2.0 Output folders
# -------------------------
dir.create("output", showWarnings = FALSE)
dir.create("output/section2_descriptive", recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 2.1 Descriptive table by pregnancy intention
# -------------------------

# Variables to summarize
vars_cat <- c(
  "preterm",
  "race4",
  "agecon_grp",
  "parity_cat",
  "marital_concep",
  "educ_cat",
  "pov_cat",
  "curr_ins_cat",
  "bmi_cat4"
)

# Nice labels for the table
var_labels <- c(
  preterm        = "Preterm birth outcome",
  race4          = "Race / ethnicity",
  agecon_grp     = "Maternal age at conception (years)",
  parity_cat     = "Parity at interview",
  marital_concep = "Marital/cohabitation status at conception",
  educ_cat       = "Highest education",
  pov_cat        = "Family income (% of poverty line)",
  curr_ins_cat   = "Current health insurance",
  bmi_cat4       = "Pre-pregnancy BMI category"
)

# Function to create one variable’s cross-tab (weighted) by intention
make_one_desc <- function(var_name, design) {
  f <- as.formula(paste0("~", var_name, "+unintended_bin"))
  tab <- svytable(f, design)
  
  df <- as.data.frame(tab)
  # First column is the variable levels, second is unintended_bin
  names(df)[1:2] <- c("level", "unintended_bin")
  
  df %>%
    group_by(unintended_bin) %>%
    mutate(
      pct  = Freq / sum(Freq) * 100,
      cell = sprintf("%.0f (%.1f%%)", Freq, pct)
    ) %>%
    ungroup() %>%
    select(level, unintended_bin, cell) %>%
    mutate(variable = var_name) %>%
    # wide: intended vs unintended
    tidyr::pivot_wider(
      names_from  = unintended_bin,
      values_from = cell,
      values_fill = list(cell = "0 (0.0%)")
    )
}

desc_list <- lapply(vars_cat, make_one_desc, design = nsfg_design)
desc_tbl  <- bind_rows(desc_list)

# Build display table
tbl2_descriptive_df <- desc_tbl %>%
  mutate(
    Characteristic = var_labels[variable],
    Level          = as.character(level)
  ) %>%
  group_by(Characteristic) %>%
  mutate(
    Characteristic = ifelse(row_number() == 1, Characteristic, "")
  ) %>%
  ungroup() %>%
  # columns "intended" and "unintended" come from unintended_bin factor levels
  select(Characteristic, Level, intended, unintended) %>%
  rename(
    `Intended pregnancy`   = intended,
    `Unintended pregnancy` = unintended
  )

tbl2_descriptive <- tbl2_descriptive_df %>%
  gt() %>%
  tab_header(
    title = "Weighted Characteristics of Live-born Singleton Pregnancies by Pregnancy Intention"
  )

gtsave(
  tbl2_descriptive,
  filename = "output/section2_descriptive/table2_descriptive_by_intention.html"
)

# -------------------------
# 2.2 Unadjusted preterm prevalence by intention
# -------------------------

# Create indicator for preterm
nsfg_design_preterm <- update(
  nsfg_design,
  preterm_ind = as.numeric(preterm == "preterm")
)

# Compute weighted prevalence + 95% CI for each intention level
preterm_list <- lapply(
  c("intended", "unintended"),
  function(grp) {
    d_sub <- subset(nsfg_design_preterm, unintended_bin == grp)
    est   <- svymean(~preterm_ind, d_sub)
    ci    <- confint(est)
    
    data.frame(
      intention    = grp,
      preterm_pct  = as.numeric(coef(est)[1]) * 100,
      ci_low_pct   = as.numeric(ci[1, 1]) * 100,
      ci_high_pct  = as.numeric(ci[1, 2]) * 100
    )
  }
)

preterm_by_intention <- bind_rows(preterm_list)

tbl2_preterm <- preterm_by_intention %>%
  mutate(intention = factor(intention, levels = c("intended", "unintended"))) %>%
  arrange(intention) %>%
  gt() %>%
  fmt_number(
    columns  = c(preterm_pct, ci_low_pct, ci_high_pct),
    decimals = 1
  ) %>%
  cols_label(
    intention   = "Pregnancy intention",
    preterm_pct = "Preterm birth (%)",
    ci_low_pct  = "95% CI (lower)",
    ci_high_pct = "95% CI (upper)"
  ) %>%
  tab_header(
    title = "Weighted Prevalence of Preterm Birth by Pregnancy Intention"
  )

gtsave(
  tbl2_preterm,
  filename = "output/section2_descriptive/table2_preterm_prevalence.html"
)
