/* SEER Breast Cancer -- Survival Analysis
   Rohit Kambala | May 2026
   KM curves and Cox PH in SAS, paralleling my R survival analysis.
   Dataset: seer_breast_survival.csv (351,860 obs, 5 vars) */

/* Import */
PROC IMPORT DATAFILE="/home/u64518804/seer_breast_survival.csv"
    OUT=work.seer DBMS=CSV REPLACE;
    GETNAMES=YES;
RUN;

/* Column names from SEER have spaces and long strings --
   SAS truncates them. Renaming everything to clean short names.
   Creating a binary event indicator: Dead=1, Alive/censored=0.
   Dropping Unknown/unstaged since I'm stratifying by stage. */
DATA work.seer;
    SET work.seer;
    surv_months  = "Survival months"N;
    vital_status = "Vital status recode (study cutof"N;
    stage        = "Summary stage 2000 (1998-2017)"N;
    yr_dx        = "Year of diagnosis"N;
    IF vital_status = 'Dead' THEN event = 1;
    ELSE event = 0;
    IF stage = 'Unknown/unstaged' THEN DELETE;
    DROP "Survival months"N
         "Vital status recode (study cutof"N
         "Summary stage 2000 (1998-2017)"N
         "Year of diagnosis"N
         Sex;
RUN;

/* Verify -- expecting ~105k events, 3 stage levels */
PROC FREQ DATA=work.seer;
    TABLES event stage;
RUN;

/* KM curves by stage.
   NOPRINT suppresses the full KM table -- 346k rows crashes the browser.
   ODS OUTPUT captures the log-rank test as a dataset instead.
   TIME statement syntax: outcome * event_indicator(censored_value) */
PROC LIFETEST DATA=work.seer NOPRINT
              PLOTS=SURVIVAL(ATRISK OUTSIDE);
    TIME surv_months * event(0);
    STRATA stage;
    ODS OUTPUT HomTests = work.km_logrank;
    TITLE 'SEER Breast Cancer -- KM Survival by Stage';
RUN;

PROC PRINT DATA=work.km_logrank;
    TITLE 'Log-Rank Test -- Survival by Stage';
RUN;

/* Cox PH -- stage as categorical, Localized as reference.
   Adding year of diagnosis to adjust for treatment era effects.
   TIES=EFRON matches R coxph() default so HRs are comparable.
   NOSUMMARY reduces output size. */
PROC PHREG DATA=work.seer NOSUMMARY;
    CLASS stage (REF='Localized') / PARAM=REF;
    MODEL surv_months * event(0) = stage yr_dx
                                   / TIES=EFRON RISKLIMITS;
    TITLE 'SEER Breast Cancer -- Cox PH Model';
RUN;
/* Results: Distant HR=9.29 (9.13-9.45), Regional HR=1.67 (1.65-1.69)
   yr_dx HR=0.987 -- slight improvement in survival over time
   All p<0.0001, N=346837 */