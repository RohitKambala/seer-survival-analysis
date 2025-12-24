# Survival Analysis of Breast Cancer Using SEER Data

## Overview
This project applies survival analysis methods to population-based cancer registry data from the Surveillance, Epidemiology, and End Results (SEER) program to evaluate overall survival patterns among breast cancer patients in the United States.

Using Kaplan–Meier estimation and Cox proportional hazards modeling, the analysis examines how survival varies by stage at diagnosis, sex, and diagnosis period.

## Data Source
- SEER Research Data (public-use)
- Breast cancer cases diagnosed between 2010–2015
- Accessed via SEER*Stat software

⚠️ Raw SEER data are not included due to data use restrictions.

## Methods
- Kaplan–Meier survival curves
- Log-rank tests
- Cox proportional hazards regression
- Proportional hazards diagnostics using Schoenfeld residuals

## Key Findings
- Stage at diagnosis is the strongest predictor of survival.
- Localized disease shows substantially higher survival than regional or distant stages.
- Survival modestly improves for more recent diagnosis periods.
- Sex differences in survival are present but smaller in magnitude.

## Files
- `seer_survival_analysis.Rmd`: Fully reproducible analysis
- `seer_survival_analysis.pdf`: Final rendered report

## Tools
- R
- tidyverse
- survival
- survminer

## Author
Rohit Kambala  
MPH (Biostatistics)
