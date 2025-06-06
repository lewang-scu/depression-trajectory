# Project Description

This repository contains statistical analysis scripts used in the study of associations between blood cell markers and depression status based on repeated physical examination records.

## 📊 Dataset Overview

- The dataset includes **repeated measurements** from participants undergoing routine health check-ups.
- Each row represents **one visit** by a participant.

## 🧬 Key Variables

| Variable | Description |
|----------|-------------|
| `dah`    | Unique identifier for each participant. All rows with the same `dah` refer to repeat visits by the same individual. |
| `jcsj`   | Timestamp of the physical examination. |
| `status` | Binary indicator of depression diagnosis status at the next visit:<br>`0` = no depression, `1` = depression. |
| `class`  | Ordinal depression severity at the next visit:<br>Levels: no depression, mild, moderate, severe. |
| Others   | Blood cell parameters, derived ratios (e.g., PLR, MLR, NLR), and health behaviors (e.g., smoking, drinking). |

---

## 📁 Indicator File: `indicator_name.xlsx`

This file contains a list of blood cell biomarkers analyzed in the study.

Columns typically include:

- **Indicator name in Chinese**
- **English translation** *(optional)*
- **Abbreviations** (e.g., PLR = Platelet-to-Lymphocyte Ratio)
- **Units** *(if applicable)*

---

## 🧪 Analysis Overview

The scripts in this repository cover the following statistical models:

- **Poisson regression**:  
  Used to model depression incidence rate; includes dispersion estimation.  
  Robust standard errors are computed using `sandwich::vcovCL()` and `lmtest::coeftest()`.

- **Mixed-effects models**:  
  Account for within-person correlations in repeated measures (`lme4::glmer`).

- **Logistic regression**:  
  Model the association between blood markers and binary depression outcome (`status`).

- **Ordinal logistic regression**:  
  Assess the relationship between biomarkers and depression severity (`class`).

- **Trajectory modeling**:  
  Identify longitudinal patterns of biomarkers over time using latent class modeling (`lcmm::hlme`).

---

## 📂 Scripts Included

- `poisson_regression.R`  
  → Time-updated Poisson models with robust standard errors.

- `poisson_dispersion_analysis.R`  
  → Computes dispersion statistics (Pearson χ² / df).

- `mixed-effect_model.R`  
  → Generalized linear mixed models for repeated outcomes.

- `logistic_and_ordinal_logistic_regression_code.R`  
  → Binary and ordinal logistic regression models.

- `trajectory_analysis.R`  
  → Latent class trajectory models over time.

---

## 🧷 Reproducibility

- All analyses were performed using **R**.
- Robust inference and modeling details are documented within each script.
- Dependencies: `dplyr`, `data.table`, `openxlsx`, `glm2`, `lme4`, `MASS`, `lcmm`, `sandwich`, `lmtest`, etc.

For transparency, please refer to each `.R` file or cite this repository directly in your work.


### Contributors

- **Wang Le** – conceptualization, data preprocessing, statistical modeling  
- **Lin Yifei** – code optimization, reproducibility review  
- **Huang Jin** – project supervision, methodological guidance

