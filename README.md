# AST/ALT ratio and mortality in sepsis

Analysis code for the two-database cohort study of the admission AST/ALT ratio
and short-term mortality in adults with sepsis. MIMIC-IV version 3.1 is the
development database. The same quartile cutpoints are applied in the eICU
Collaborative Research Database version 2.0.

This repository contains the code used for the cohorts, multiple imputation,
time-split models, and figures in the manuscript. It does not contain
patient-level records. The databases are available to credentialed users
through PhysioNet.

## Requirements

Python 3.14.5, with pandas 3.0.4, patsy 1.0.2, SciPy 1.18.0, statsmodels 0.14.6,
scikit-learn 1.9.0, lifelines 0.30.0, matplotlib, and psycopg2.

## Connection

The scripts connect to PostgreSQL on localhost, port 5442, as user `postgres`.
Set `PGPASSWORD` in the environment before running a script. The password is
not stored in this repository. The local database names are `mimiciv31` and `eicu`.

Run the scripts from the repository root. Quartile cutpoints are frozen in
`deritis_paper_analysis.py` and are not recomputed. Five stays whose enzyme
chart time is later than the hospital death time are removed in `deritis_drop5.py`
before the reported primary models.

## Scripts

- `deritis_paper_analysis.py` builds the cohorts.
- `deritis_round2.py` and `deritis_revision.py` define the imputation and Cox helpers.
- `deritis_mi_q.py` fits the primary imputed hazard ratios.
- `deritis_period_q.py` and `deritis_cuts_q.py` fit the day-7, day-3, and day-14 splits.
- `deritis_drop5.py` repeats the primary models after the five timing exclusions.
- `deritis_eicu_cluster_ok.py` fits the eICU model with hospital-clustered standard errors.
- `deritis_rerun_two.py` fits the joint enzyme model and the expanded-covariate model.
- `deritis_ipw_all.py` fits the inverse-probability weighted model.
- `deritis_paper/make_figures.py` writes the figures.

`manuscript_full.md` is the manuscript source that corresponds to this code.
