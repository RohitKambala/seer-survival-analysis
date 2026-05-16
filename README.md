# Survival Analysis of Breast Cancer Using SEER Data

## Overview

This project applies survival analysis methods to population-based cancer 
registry data from the Surveillance, Epidemiology, and End Results (SEER) 
program to evaluate overall survival patterns among breast cancer patients 
in the United States.

Using Kaplan–Meier estimation, Cox proportional hazards modeling, and 
competing risks analysis, the analysis examines how survival varies by 
stage at diagnosis, sex, and diagnosis period. A parallel SAS implementation 
is included alongside the R analysis.

## Data Source

- SEER Research Data (public-use)
- Breast cancer cases diagnosed between 2010–2015
- Accessed via SEER*Stat software

⚠️ Raw SEER data are not included due to data use restrictions.

## Methods

- Kaplan–Meier survival curves and log-rank tests
- Cox proportional hazards regression (unadjusted and multivariable)
- Proportional hazards diagnostics using Schoenfeld residuals
- Competing risks analysis: cumulative incidence functions and Gray's test

## Key Findings

- Stage at diagnosis is the strongest predictor of survival
- Localized disease shows substantially better survival than regional or distant stages
- For localized patients, other-cause mortality exceeds breast cancer death over long follow-up
- Survival modestly improves for more recently diagnosed patients
- Sex differences are present but small and should be interpreted cautiously

## Files

- `seer_survival_analysis.Rmd`: Fully reproducible R analysis
- `seer_survival_analysis.pdf`: Final rendered report
- `seer_survival.sas`: SAS implementation (PROC LIFETEST, PROC PHREG)

## Tools

- R: tidyverse, survival, survminer, cmprsk
- SAS: PROC LIFETEST, PROC PHREG

## Author

Rohit Kambala  
MPH Biostatistics, New York University