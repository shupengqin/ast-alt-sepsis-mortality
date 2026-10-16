# AST/ALT ratio and mortality in sepsis

Analysis code for a two-database cohort study of the admission AST/ALT ratio
and short-term mortality in adults with sepsis and no chronic liver disease.
Development data are MIMIC-IV version 3.1. The same quartile cutpoints are
applied in the eICU Collaborative Research Database version 2.0.

The databases are available to credentialed users through PhysioNet. This
repository contains analysis code only. It does not contain patient-level records.

## Requirements

Python 3.11 or later, with pandas, numpy, scipy, patsy, statsmodels,
scikit-learn, lifelines, matplotlib, and psycopg2.

## Connection

The scripts connect to PostgreSQL on localhost, port 5442, as user `postgres`.
Set `PGPASSWORD` in the environment before running a script. The password is
not stored in this repository.

The local database names used by the scripts are `mimiciv31` and `eicu`.

## Scripts

- `deritis_paper_analysis.py` builds the cohorts and the main models.
- `deritis_round2.py` runs the multiple imputation.
- `deritis_period_test.py` estimates the period-by-quartile interaction.
- `deritis_closeout.py` runs the proportional-hazards test, the covariate-by-period sensitivity analysis, and the eICU imputation.
- `deritis_paper/make_figures.py` writes the figures.

`manuscript_full.md` is the manuscript source that corresponds to this code.
