# =========================
# 0. Setup
# =========================
library(tidyverse)

# Optional: for later survey work (Sections 2–4)
# install.packages("survey")
# library(survey)

# Paths to the CSVs (inside the "data" folder in your project)
preg_path <- "data/NSFG_2022_2023_FemPregPUFData.csv"
resp_path <- "data/NSFG_2022_2023_FemRespPUFData.csv"


# =========================
# 1. Load data
# =========================
preg_raw <- readr::read_csv(preg_path, show_col_types = FALSE)
resp_raw <- readr::read_csv(resp_path, show_col_types = FALSE)

# Keep only what we need from respondent file for now (BMI)
resp_bmi <- resp_raw %>%
  select(CaseID, BMIcat)

# Merge BMIcat into pregnancy file
preg_merged <- preg_raw %>%
  left_join(resp_bmi, by = "CaseID")

# =========================
# 2. Apply inclusion / exclusion criteria
# =========================

preg_cohort <- preg_merged %>%
  # 2.1 Live births only
  filter(OUTCOME == 1) %>%
  # 2.2 Singleton births only
  filter(BORNALIV == 1) %>%
  # 2.3 Valid gestational length categories (1–4)
  filter(!is.na(GEST_LB),
         GEST_LB %in% c(1, 2, 3, 4)) %>%
  # 2.4 Valid intention (exclude "didn't care"/"don't know" for main analysis)
  filter(wantresp %in% c(1, 2, 3, 5))

# (Optional sensitivity analysis sample later)
preg_cohort_recent <- preg_cohort %>%
  filter(RECNT5YRPRG == 1)

# =========================
# 3. Create outcome & exposure
# =========================

preg_cohort <- preg_cohort %>%
  mutate(
    # 3.1 Outcome: preterm vs term/post-term
    preterm = case_when(
      GEST_LB %in% c(1, 2) ~ 1L,  # early or late preterm
      GEST_LB %in% c(3, 4) ~ 0L,  # term or post-term
      TRUE ~ NA_integer_
    ),
    
    # 3.2 Detailed intention categories
    intend_cat = case_when(
      wantresp == 2 ~ "intended_right_time",
      wantresp %in% c(1, 3) ~ "mistimed",
      wantresp == 5 ~ "unwanted",
      wantresp == 4 ~ "didnt_care",
      wantresp == 6 ~ "dont_know",
      TRUE ~ NA_character_
    ),
    
    # 3.3 Binary exposure: unintended vs intended
    unintended_bin = case_when(
      wantresp == 2 ~ 0L,                # intended
      wantresp %in% c(1, 3, 5) ~ 1L,     # mistimed or unwanted
      TRUE ~ NA_integer_
    )
  ) %>%
  # Drop any rows with undefined outcome or exposure
  filter(!is.na(preterm),
         !is.na(unintended_bin))

# =========================
# 4. Recode covariates
# =========================

nsfg_preg_analytic <- preg_cohort %>%
  mutate(
    # 4.1 Outcome & exposure as factors
    preterm = factor(preterm,
                     levels = c(0, 1),
                     labels = c("term_or_postterm", "preterm")),
    unintended_bin = factor(unintended_bin,
                            levels = c(0, 1),
                            labels = c("intended", "unintended")),
    
    # 4.2 Race/ethnicity (HISPRACE2)
    race4 = factor(HISPRACE2,
                   levels = c(1, 2, 3, 4),
                   labels = c("Hispanic",
                              "NH_White",
                              "NH_Black",
                              "NH_Other")),
    
    # 4.3 Parity (PARITY)
    parity_cat = case_when(
      PARITY == 0 ~ "0",
      PARITY == 1 ~ "1",
      PARITY == 2 ~ "2",
      PARITY >= 3 ~ "3+",
      TRUE ~ NA_character_
    ),
    parity_cat = factor(parity_cat,
                        levels = c("0", "1", "2", "3+")),
    
    # 4.4 Maternal age at conception (agecon)
    # You can keep numeric, but define a grouped version too:
    agecon_grp = case_when(
      agecon < 20 ~ "<20",
      agecon >= 20 & agecon <= 24 ~ "20-24",
      agecon >= 25 & agecon <= 29 ~ "25-29",
      agecon >= 30 & agecon <= 34 ~ "30-34",
      agecon >= 35 & agecon <= 39 ~ "35-39",
      agecon >= 40 ~ "40+",
      TRUE ~ NA_character_
    ),
    agecon_grp = factor(agecon_grp,
                        levels = c("<20","20-24","25-29","30-34","35-39","40+")),
    
    # 4.5 Marital/cohabitation status at conception (RMARCON3)
    marital_concep = factor(RMARCON3,
                            levels = c(1, 2, 3),
                            labels = c("married",
                                       "cohabiting",
                                       "neither")),
    
    # 4.6 Education (HIEDUC) collapsed
    educ_cat = case_when(
      HIEDUC %in% c(1, 2, 3, 4) ~ "<=HS",
      HIEDUC %in% c(5, 6, 7)    ~ "SomeCollege/AA",
      HIEDUC %in% c(8, 9, 10, 11) ~ "BAplus",
      TRUE ~ NA_character_
    ),
    educ_cat = factor(educ_cat,
                      levels = c("<=HS", "SomeCollege/AA", "BAplus")),
    
    # 4.7 Poverty ratio (POVERTY) grouped
    pov_cat = case_when(
      POVERTY < 100 ~ "<100%",
      POVERTY >= 100 & POVERTY < 200 ~ "100-199%",
      POVERTY >= 200 & POVERTY < 400 ~ "200-399%",
      POVERTY >= 400 ~ ">=400%",
      TRUE ~ NA_character_
    ),
    pov_cat = factor(pov_cat,
                     levels = c("<100%", "100-199%", "200-399%", ">=400%")),
    
    # 4.8 Insurance coverage (CURR_INS)
    curr_ins_cat = factor(CURR_INS,
                          levels = c(1, 2, 3, 4),
                          labels = c("private_or_Medigap",
                                     "Medicaid_CHIP_state",
                                     "Medicare_military_othergov",
                                     "single_service_or_uninsured")),
    
    # 4.9 BMI category (BMIcat from respondent file)
    bmi_cat4 = factor(BMIcat,
                      levels = c(1, 2, 3, 4, 5),
                      labels = c("underweight",
                                 "normal",
                                 "overweight",
                                 "obese",
                                 "BMI_undef"))
  )

# Optional: check dimensions and a quick glimpse
dim(nsfg_preg_analytic)
dplyr::glimpse(nsfg_preg_analytic, width = 120)
