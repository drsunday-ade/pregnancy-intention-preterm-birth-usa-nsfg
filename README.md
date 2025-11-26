# pregnancy-intention-preterm-birth-usa-nsfg
Reproducible R code, analytic cohort, and outputs for a study of unintended pregnancy and risk of preterm birth in the United States using 2022–23 National Survey of Family Growth data, with survey-weighted causal inference and risk-prediction models.

Unintended pregnancy and preterm birth in the USA: NSFG causal and prediction study

Repository name: pregnancy-intention-preterm-birth-usa-nsfg

This repository contains data-processing and analysis code for the manuscript:

Unintended Pregnancy and Risk of Preterm Birth in the United States:
Causal Inference and Risk Prediction Using National Survey of Family Growth Data

Short title: Unintended pregnancy and preterm birth in the USA: a causal inference and risk prediction study

Authors

Sunday A Adetunji, MD¹˒²*

¹ College of Health, Oregon State University, Corvallis, USA

² Faculty of Clinical Sciences, Obafemi Awolowo University, Ile-Ife, Nigeria

ORCID: https://orcid.org/0000-0001-9321-9957

Corresponding author: adetunjs@oregonstate.edu

Purva Nerandra More, MPH¹

¹ College of Health, Oregon State University, Corvallis, USA

ORCID: https://orcid.org/0009-0006-8262-5857

Study overview

We used nationally representative data from the 2022–2023 National Survey of Family Growth (NSFG) Female Pregnancy file to evaluate whether unintended pregnancy is associated with increased risk of preterm birth among U.S. women aged 15–49 years, and to develop a clinically interpretable risk-prediction model.

Key components of the workflow:

Construction of an analytic cohort of recent singleton live births with complete information on pregnancy intention, gestational age, and key covariates.

Descriptive analyses of pregnancy intention, preterm birth, and social and clinical risk factors.

Causal effect estimation using survey-weighted propensity score and modern causal machine-learning approaches.

Development and internal validation of a prediction model for preterm birth that incorporates pregnancy intention and other risk factors.

Extensive robustness checks, diagnostics, and clinical translation tools (e.g., decision-curve analysis).

All analyses account for the complex survey design and weighting of the NSFG.

Repository structure
pregnancy-intention-preterm-birth-usa-nsfg/
├── data/
│   ├── NSFG_2022_23_FemPreg_PUFData.csv
│   └── NSFG_2022_23_FemPreg_readme.txt
│   (plus any derived/clean analytic datasets created by the scripts)
│
├── scripts/
│   ├── 1.building_the_analytic.R
│   ├── Section 2_Descriptive_Analyses.R
│   ├── Section 3_Risk_Model.R
│   ├── Section 4_Causal_Inference_Full.R
│   ├── Section 5_Prediction_and_Stratification.R
│   ├── Section 6_Sensitivity_Analyses.R
│   ├── Section 7_Defence_Summary.R
│   └── Section 8_Translation_and_Decision_Curve.R
│
└── output/
    ├── section2_descriptive/
    │   └── (descriptive tables and prevalence estimates, HTML + figures)
    ├── section3_riskmodel/
    │   └── (regression models, adjusted risk estimates, ROC and calibration plots)
    ├── section4_causal/
    │   └── (weight diagnostics, causal effect estimates, Bayesian/posterior summaries)
    ├── section5_prediction/
    │   └── (model performance, calibration by deciles, high-risk classification tables)
    ├── section6_robustness/
    │   └── (propensity score overlap, alternative specifications, bounding analyses)
    ├── section7_defence/
    │   └── (summaries for methodological defence and peer review)
    └── section8_translation/
        └── (decision-curve analysis and clinical translation outputs)


data/ contains the original NSFG public-use file and documentation, and any derived analytic datasets.

scripts/ holds the R code used to build the analytic cohort and run each major section of the analysis.

output/ contains reproducible tables (.html) and figures (.png) corresponding to the manuscript’s main and supplementary results.

Data source and availability

This project uses only openly available de-identified human data from the U.S. National Survey of Family Growth (NSFG), distributed by the National Center for Health Statistics (NCHS), Centers for Disease Control and Prevention.

Original data source: NSFG 2022–2023 Female Pregnancy Public-Use Data File (PUF).

How to obtain the data: The NSFG public-use files, documentation, and codebooks can be downloaded free of charge from the NSFG/NCHS website.

In this repository: A copy of the 2022–23 Female Pregnancy PUF .csv file and original NCHS data README are stored in data/ for convenience. Any use of these data must comply with NCHS guidance and users must not attempt to re-identify participants.

No new data were collected for this study.

How to reproduce the analyses
Prerequisites

R version 4.3.x (or later)

Recommended R packages (install if needed):
haven, dplyr, tidyr, stringr, ggplot2,
survey, srvyr, WeightIt, cobalt,
broom, tableone,
SuperLearner, randomForest, glmnet, ranger, xgboost,
and additional packages loaded at the top of each script.

Steps

Clone the repository

git clone https://github.com/drsunday-ade/pregnancy-intention-preterm-birth-usa-nsfg.git
cd pregnancy-intention-preterm-birth-usa-nsfg


Open R / RStudio in the project directory.

Build the analytic cohort

Run:

source("scripts/1.building_the_analytic.R")


This script:

Imports the NSFG Female Pregnancy PUF,

Applies inclusion/exclusion criteria,

Creates derived variables (pregnancy intention, preterm birth, covariates),

Saves the final analytic dataset used by downstream scripts.

Descriptive analyses

source("scripts/Section 2_Descriptive_Analyses.R")


Produces descriptive tables and prevalence estimates saved in output/section2_descriptive/.

Risk model and causal inference

source("scripts/Section 3_Risk_Model.R")
source("scripts/Section 4_Causal_Inference_Full.R")


These scripts fit multivariable regression models, estimate causal effects using survey-weighted and modern causal methods, and save diagnostics and effect estimates to output/section3_riskmodel/ and output/section4_causal/.

Prediction, robustness, and translation

source("scripts/Section 5_Prediction_and_Stratification.R")
source("scripts/Section 6_Sensitivity_Analyses.R")
source("scripts/Section 7_Defence_Summary.R")
source("scripts/Section 8_Translation_and_Decision_Curve.R")


These scripts build and internally validate the prediction model, conduct robustness checks, summarise methodological checks, and generate decision-curve and translation outputs.

The exact numbers reported in the manuscript (e.g., weighted counts ~78·4 million U.S. births) should be reproduced by running the full pipeline on the NSFG 2022–23 Female Pregnancy file.

How to cite

If you use this code or workflow, please cite:

The manuscript (once published) – provisional citation:

Adetunji SA, More PN. Unintended Pregnancy and Risk of Preterm Birth in the United States: Causal Inference and Risk Prediction Using National Survey of Family Growth Data. [Manuscript in review].

The data source:

National Center for Health Statistics. National Survey of Family Growth, 2022–2023 Female Pregnancy Public-Use Data File. Hyattsville, MD: Centers for Disease Control and Prevention.

License

Unless otherwise indicated in individual files, the analysis code in this repository is released under an open-source license (you can specify e.g. MIT or GPL here once you decide). Use of the NSFG data must follow NCHS terms and conditions.

Contact

For questions about the analysis or collaboration enquiries, please contact:

Dr Sunday A Adetunji
College of Health, Oregon State University
Email: adetunjs@oregonstate.edu
