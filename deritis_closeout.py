"""Schoenfeld test, covariate-by-period sensitivity, and eICU imputation."""

from __future__ import annotations

import os
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import sklearn
from lifelines import CoxPHFitter, NelsonAalenFitter
from lifelines import CoxTimeVaryingFitter
from lifelines.statistics import proportional_hazard_test
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer

from deritis_paper_analysis import load_eicu, load_mimic
from deritis_round2 import add_covariates, rubin

OUT = Path(r"F:\CKB数据库\analysis\deritis_paper")
USE = [
    "age", "male", "sofa", "lactate", "log_bili", "creatinine", "aud",
    "u_CCU", "u_SICU", "u_CVICU", "u_Other",
]


def schoenfeld(df):
    need = ["time_28", "event_28", "age", "male", "sofa", "lactate", "log_bili", "creatinine", "aud", "unit", "quartile"]
    d = df.dropna(subset=need).copy()
    d["q2"] = (d["quartile"] == "Q2").astype(float)
    d["q3"] = (d["quartile"] == "Q3").astype(float)
    d["q4"] = (d["quartile"] == "Q4").astype(float)
    for unit in ("CCU", "SICU", "CVICU", "Other"):
        d[f"u_{unit}"] = (d["unit"] == unit).astype(float)
    cols = [
        "time_28", "event_28", "q2", "q3", "q4", "age", "male", "sofa", "lactate",
        "log_bili", "creatinine", "aud", "u_CCU", "u_SICU", "u_CVICU", "u_Other",
    ]
    use = d[cols]
    model = CoxPHFitter()
    model.fit(use, "time_28", "event_28")
    tested = proportional_hazard_test(model, use, time_transform="rank")
    summary = tested.summary.reset_index()
    summary.to_csv(OUT / "ph_test_current.csv", index=False)
    return {str(row["index"]): float(row["p"]) for _, row in summary.iterrows()}


def varying_fit(model_df):
    base = model_df.copy()
    base["id"] = np.arange(len(base))
    early = base.copy()
    early["start"] = 0.0
    early["stop"] = np.minimum(early["time_28"], 7.0)
    early["event"] = ((early["event_28"] == 1) & (early["time_28"] <= 7)).astype(int)
    early["late"] = 0.0
    late = base[base["time_28"] > 7].copy()
    late["start"] = 7.0
    late["stop"] = late["time_28"]
    late["event"] = late["event_28"].astype(int)
    late["late"] = 1.0
    both = pd.concat([early, late], ignore_index=True)
    both = both[both["stop"] > both["start"]]
    for name in ("q2", "q3", "q4", "age", "sofa", "lactate", "log_bili"):
        both[f"{name}_late"] = both[name] * both["late"]
    cols = [
        "id", "start", "stop", "event",
        "q2", "q3", "q4", "q2_late", "q3_late", "q4_late",
        "age", "age_late", "male", "sofa", "sofa_late", "lactate", "lactate_late",
        "log_bili", "log_bili_late", "creatinine", "aud",
        "u_CCU", "u_SICU", "u_CVICU", "u_Other",
    ]
    fit = CoxTimeVaryingFitter()
    fit.fit(both[cols], id_col="id", event_col="event", start_col="start", stop_col="stop", show_progress=False)
    names = list(fit.params_.index)
    var = np.asarray(fit.variance_matrix_, float)
    i_early = names.index("q4")
    i_ratio = names.index("q4_late")
    coef_early = float(fit.params_["q4"])
    coef_ratio = float(fit.params_["q4_late"])
    var_early = float(var[i_early, i_early])
    var_ratio = float(var[i_ratio, i_ratio])
    cov = float(var[i_early, i_ratio])
    return {
        "q4": (coef_early, var_early),
        "q4_late": (coef_ratio, var_ratio),
        "q4_after": (coef_early + coef_ratio, var_early + var_ratio + 2 * cov),
    }


def sensitivity(df):
    d = df.copy()
    d["q2"] = (d["quartile"] == "Q2").astype(float)
    d["q3"] = (d["quartile"] == "Q3").astype(float)
    d["q4"] = (d["quartile"] == "Q4").astype(float)
    for unit in ("CCU", "SICU", "CVICU", "Other"):
        d[f"u_{unit}"] = (d["unit"] == unit).astype(float)
    work = d.dropna(subset=["time_28", "event_28", "quartile"]).copy()
    naf = NelsonAalenFitter()
    naf.fit(work["time_28"], work["event_28"])
    work["na_hazard"] = np.asarray(naf.cumulative_hazard_at_times(work["time_28"]).to_numpy(), float).ravel()
    imp_cols = USE + ["event_28", "na_hazard"]
    base = work[imp_cols]
    store = {term: {"coef": [], "var": []} for term in ("q4", "q4_late", "q4_after")}
    for seed in range(30):
        if seed % 5 == 0:
            print(f"sensitivity {seed + 1}/30", flush=True)
        imputer = IterativeImputer(random_state=seed, max_iter=12, sample_posterior=True)
        filled = pd.DataFrame(imputer.fit_transform(base), columns=imp_cols, index=base.index)
        for col in ("male", "aud", "u_CCU", "u_SICU", "u_CVICU", "u_Other"):
            filled[col] = (filled[col] >= 0.5).astype(float)
        model_df = filled[USE].copy()
        model_df["time_28"] = work["time_28"]
        model_df["event_28"] = work["event_28"]
        model_df["q2"] = work["q2"]
        model_df["q3"] = work["q3"]
        model_df["q4"] = work["q4"]
        fitted = varying_fit(model_df)
        for term, (coef, var) in fitted.items():
            store[term]["coef"].append(coef)
            store[term]["var"].append(var)
    return {term: rubin(store[term]["coef"], store[term]["var"]) for term in store}


def eicu_mi(eicu):
    d = eicu.copy()
    unit = d["unittype"].fillna("")
    d["unit"] = np.select(
        [
            unit.isin(["CSICU", "Cardiac ICU", "CCU-CTICU", "CTICU"]),
            unit.str.contains("MICU|Med-Surg", case=False),
            unit.str.contains("CCU|Cardiac", case=False),
            unit.str.contains("SICU|Surgical", case=False),
        ],
        ["Cardiac", "MICU", "CCU", "SICU"],
        default="Other",
    )
    for name, label in (("u_MICU", "MICU"), ("u_CCU", "CCU"), ("u_SICU", "SICU"), ("u_Cardiac", "Cardiac")):
        d[name] = (d["unit"] == label).astype(float)
    miss = {col: float(d[col].isna().mean()) for col in ("age", "apachescore", "lactate", "log_bili", "creatinine", "y")}
    imp_cols = ["age", "male", "apachescore", "lactate", "log_bili", "creatinine", "y", "aud", "u_MICU", "u_CCU", "u_SICU", "u_Cardiac"]
    work = d.dropna(subset=["y", "quartile"]).copy()
    work["q2"] = (work["quartile"] == "Q2").astype(float)
    work["q3"] = (work["quartile"] == "Q3").astype(float)
    work["q4"] = (work["quartile"] == "Q4").astype(float)
    base = work[imp_cols]
    store = {term: {"coef": [], "var": []} for term in ("q2", "q3", "q4")}
    import statsmodels.api as sm
    for seed in range(30):
        if seed % 10 == 0:
            print(f"eicu {seed + 1}/30", flush=True)
        imputer = IterativeImputer(random_state=seed, max_iter=12, sample_posterior=True)
        filled = pd.DataFrame(imputer.fit_transform(base), columns=imp_cols, index=base.index)
        for col in ("male", "y", "aud", "u_MICU", "u_CCU", "u_SICU", "u_Cardiac"):
            filled[col] = (filled[col] >= 0.5).astype(float)
        X = filled[["age", "male", "apachescore", "lactate", "log_bili", "creatinine", "aud", "u_MICU", "u_CCU", "u_SICU", "u_Cardiac"]].copy()
        X["q2"] = work["q2"]
        X["q3"] = work["q3"]
        X["q4"] = work["q4"]
        X = sm.add_constant(X, has_constant="add")
        keep = [col for col in X.columns if float(X[col].std()) > 1e-8]
        X = X[keep]
        try:
            res = sm.Logit(work["y"], X).fit(disp=False, maxiter=200)
        except Exception:
            X = X.drop(columns=[c for c in X.columns if str(c).startswith("u_")], errors="ignore")
            res = sm.Logit(work["y"], X).fit(disp=False, maxiter=200)
        cov = np.asarray(res.cov_params(), float)
        names = list(X.columns)
        for term in ("q2", "q3", "q4"):
            idx = names.index(term)
            store[term]["coef"].append(float(res.params[term]))
            store[term]["var"].append(float(cov[idx, idx]))
    return {
        "missing": miss,
        "n": int(len(work)),
        "events": int(work["y"].sum()),
        "terms": {term: rubin(store[term]["coef"], store[term]["var"]) for term in store},
    }


def main():
    mimic, cuts, _ = load_mimic()
    mimic = add_covariates(mimic)
    print("schoenfeld", flush=True)
    ph = schoenfeld(mimic)
    print("sensitivity", flush=True)
    sens = sensitivity(mimic)
    eicu = load_eicu(cuts)
    # alcohol flag is created inside eicu_with_unit; repeat the small query here
    import psycopg2
    conn = psycopg2.connect(host="localhost", port=5442, user="postgres", password=os.environ["PGPASSWORD"], dbname="eicu")
    conn.autocommit = True
    with conn.cursor() as cur:
        cur.execute("SET statement_timeout = 0")
    alcohol = pd.read_sql_query(
        """
        SELECT DISTINCT patientunitstayid FROM diagnosis
        WHERE lower(diagnosisstring) LIKE '%alcohol%'
           OR icd9code LIKE '303%' OR icd9code LIKE '305.0%' OR icd9code LIKE '291%'
        """,
        conn,
    )
    conn.close()
    eicu["aud"] = eicu["patientunitstayid"].isin(set(alcohol["patientunitstayid"])).astype(float)
    print("eicu imputation", flush=True)
    eicu_out = eicu_mi(eicu)
    payload = {
        "versions": {"scipy": scipy.__version__, "sklearn": sklearn.__version__},
        "schoenfeld_p": ph,
        "covariate_time_sensitivity": sens,
        "eicu_mi": eicu_out,
        "before_enzyme": json.loads((OUT / "before_enzyme.json").read_text(encoding="utf-8")),
    }
    (OUT / "closeout.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
