/* ============================================================
   SEER Breast Cancer -- Survival Analysis
   Rohit Kambala | May 2026

   Dataset: seer_breast_survival.csv (351,860 obs, 5 vars)
   Goal: KM curves and Cox PH in SAS, paralleling my R
   survival analysis on SEER breast cancer data.

   Analysis structure:
   (1) Import and clean raw SEER export
   (2) Verify event counts and stage distribution
   (3) PROC LIFETEST -- KM curves and log-rank test by stage
   (4) PROC PHREG -- Cox PH model with stage and diagnosis year

   Key results (May 2026):
   Distant HR=9.29 (9.13-9.45), Regional HR=1.67 (1.65-1.69)
   yr_dx HR=0.987 -- ~23% lower hazard over 20 years of data
   Log-rank chi-sq=86,357, p<0.0001
   Distant median survival=28 months vs Localized median=not reached
   All results match R output -- cross-language validation done.

   Clinical interpretation:
   At any given month of follow-up, a Distant stage patient
   is dying at 9.29 times the rate of a Localized patient,
   after adjusting for diagnosis year. That's one of the
   largest HRs in real-world oncology data. yr_dx coefficient
   captures two decades of treatment improvement: targeted
   therapy, HER2 agents, immunotherapy.
   ============================================================ */


/* ---- Step 1: Import ----------------------------------------
   SEER exports come with long column names containing spaces.
   SAS truncates "Vital status recode (study cutoff used)" to
   32 characters at import -- that truncation is expected and
   handled in the DATA step below.
   No GUESSINGROWS needed here because all numeric columns
   imported correctly on default scan.
   ------------------------------------------------------------ */
PROC IMPORT DATAFILE="/home/u64518804/seer_breast_survival.csv"
    OUT=work.seer DBMS=CSV REPLACE;
    GETNAMES=YES;
RUN;


/* ---- Step 2: Clean and recode ------------------------------
   SEER column names have spaces and are too long for clean SAS
   code. Renaming everything to short working names first.

   Event indicator: Dead=1, anything else (Alive or censored)=0.
   This is the standard right-censoring setup. Patients who were
   alive at study cutoff or lost to follow-up are censored --
   they didn't have the event, but they also didn't survive
   the full study period. The survival model handles them
   correctly rather than dropping them.

   Dropping Unknown/unstaged because we're stratifying by stage.
   Keeping them would create a fourth stratum with no clinical
   meaning and would complicate the log-rank test.

   5,023 obs dropped (351,860 minus 346,837).
   ------------------------------------------------------------ */
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


/* ---- Step 3: Verify ----------------------------------------
   Expecting ~105k events (event=1) and exactly 3 stage levels.
   If stage shows more than 3 levels, the DELETE above didn't
   catch all Unknown/unstaged values -- check for trailing
   spaces or alternate spellings in the raw data.
   ------------------------------------------------------------ */
PROC FREQ DATA=work.seer;
    TABLES event stage;
    TITLE 'SEER: event counts and stage distribution';
RUN;


/* ---- Step 4: KM curves and log-rank test -------------------
   PROC LIFETEST produces two things we want:
   (1) The KM survival curve plot by stage
   (2) The log-rank test of whether the three curves differ

   Why ODS EXCLUDE instead of NOPRINT:
   NOPRINT is a blunt instrument -- it suppresses ALL output
   including ODS output objects. That's why km_logrank failed
   to be created in the first version of this code. ODS EXCLUDE
   is surgical: it suppresses only the ProductLimitEstimates
   table (which is 346k rows and crashes the browser) while
   letting everything else through, including HomTests.

   TIME statement: surv_months * event(0) means:
   "time to event is surv_months, event occurred when event=1,
   event=0 means censored." This is the standard SAS syntax
   for right-censored survival data.

   STRATA stage: fits separate KM curves for each stage level
   and runs the log-rank test comparing all three curves.

   PLOTS=SURVIVAL(ATRISK OUTSIDE): adds the at-risk table
   below the x-axis. Standard in oncology publications --
   reviewers expect to see how many patients remain at risk
   at each time point.

   Log-rank test interpretation:
   Null hypothesis: all three survival curves are identical.
   We want to reject this. Chi-sq=86,357, p<0.0001 confirms
   the three curves are massively different.

   Note: log-rank tests equality. It does not quantify the
   size of the difference or adjust for covariates. That's
   what PROC PHREG does below.
   ------------------------------------------------------------ */
ODS EXCLUDE ProductLimitEstimates;

PROC LIFETEST DATA=work.seer
              PLOTS=SURVIVAL(ATRISK OUTSIDE);
    TIME surv_months * event(0);
    STRATA stage;
    ODS OUTPUT HomTests = work.km_logrank;
    TITLE 'SEER Breast Cancer -- KM Survival by Stage';
RUN;

ODS EXCLUDE NONE;

/* Print the log-rank table captured above */
PROC PRINT DATA=work.km_logrank;
    TITLE 'Log-Rank Test -- Survival by Stage';
RUN;


/* ---- Step 5: Cox Proportional Hazards model ----------------
   Cox PH does what the log-rank test cannot:
   (1) Quantifies the difference as a hazard ratio with CI
   (2) Adjusts for covariates (yr_dx here)

   CLASS stage (REF='Localized'): Localized is the reference.
   Every HR is interpreted as "X times the hazard of dying
   compared to a Localized patient diagnosed in the same year."
   Always state the reference category when reporting HRs.

   yr_dx added as a continuous covariate to adjust for
   treatment era effects. Breast cancer survival improved
   substantially from 1990s to 2010s due to better targeted
   therapy, HER2 agents, and immunotherapy. Without adjusting
   for yr_dx, stage HRs would be slightly confounded by the
   uneven distribution of diagnoses across years and stages.

   TIES=EFRON: handles tied survival times. Matches R coxph()
   default so results are directly comparable across languages.
   Efron is more accurate than Breslow for datasets with many
   ties, which is common in administrative databases like SEER.

   NOSUMMARY: suppresses the observation summary table to
   reduce output size. The FREQ table above already confirms N.

   RISKLIMITS: adds 95% confidence intervals to the HR table.
   In a pharma submission you always report HRs with CIs.

   Key results:
   Distant HR=9.29 (9.13, 9.45) -- at any given month, Distant
   patients die at 9.29x the rate of Localized patients.
   Regional HR=1.67 (1.65, 1.69) -- serious but far less severe.
   yr_dx HR=0.987 -- each year of diagnosis reduces hazard by
   1.3%. Over 20 years: 0.987^20 = 0.77, roughly 23% lower
   hazard for a 2015 diagnosis vs a 1995 diagnosis.
   ------------------------------------------------------------ */
PROC PHREG DATA=work.seer NOSUMMARY;
    CLASS stage (REF='Localized') / PARAM=REF;
    MODEL surv_months * event(0) = stage yr_dx
                                   / TIES=EFRON RISKLIMITS;
    TITLE 'SEER Breast Cancer -- Cox PH Model';
RUN;


/* ---- What to do next ---------------------------------------
   1. Add competing risks analysis (PROC LIFETEST with PLOTS=CIF)
      to account for deaths from causes other than breast cancer.
      In a real SEER analysis, competing risks matter because
      older patients may die of other causes before breast cancer
      kills them. This is the planned Q1 portfolio upgrade.
   2. Add TLF-style output: ODS RTF + PROC REPORT to produce
      a publication-quality HR table matching a clinical paper
      format.
   3. Test proportional hazards assumption -- PROC PHREG with
      ASSESS PH statement. Cox PH assumes the HR between groups
      is constant over time. The Distant curve dropping so fast
      early on is worth checking.
   ------------------------------------------------------------ */
