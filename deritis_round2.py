"""Landmark follow-up, richer covariates, multiple imputation, and interactions."""

from __future__ import annotations

import os
import json
from pathlib import Path

import numpy as np
import pandas as pd
import psycopg2
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer

from deritis_paper_analysis import cox_fit, load_eicu, load_mimic

OUT = Path(r"F:\CKB数据库\analysis\deritis_paper")
CORE = "q2 + q3 + q4 + age + male + sofa + lactate + log_bili + creatinine + aud + C(unit) - 1"
TERMS = ["q2", "q3", "q4"]


def connect():
    conn = psycopg2.connect(host="localhost", port=5442, user="postgres", password=os.environ["PGPASSWORD"], dbname="mimiciv31")
    conn.autocommit = True
    with conn.cursor() as cur:
        cur.execute("SET statement_timeout = 0")
    return conn


def add_covariates(df):
    conn = connect()
    stays = df[["stay_id", "hadm_id"]].drop_duplicates()
    charlson = pd.read_sql_query(
        "SELECT hadm_id, charlson_comorbidity_index AS charlson FROM mimiciv_derived.charlson",
        conn,
    )
    weight = pd.read_sql_query(
        "SELECT stay_id, weight FROM mimiciv_derived.first_day_weight",
        conn,
    )
    albumin = pd.read_sql_query(
        """
        SELECT DISTINCT ON (b.stay_id) b.stay_id, c.albumin
        FROM mimiciv_derived.icustay_detail b
        JOIN mimiciv_derived.chemistry c
          ON c.hadm_id = b.hadm_id
         AND c.charttime >= b.icu_intime
         AND c.charttime < b.icu_intime + interval '24 hour'
         AND c.albumin IS NOT NULL
        WHERE b.first_icu_stay
        ORDER BY b.stay_id, c.charttime
        """,
        conn,
    )
    ck = pd.read_sql_query(
        """
        SELECT DISTINCT ON (b.stay_id) b.stay_id, e.ck_cpk AS ck_first
        FROM mimiciv_derived.icustay_detail b
        JOIN mimiciv_derived.enzyme e
          ON e.hadm_id = b.hadm_id
         AND e.charttime >= b.icu_intime
         AND e.charttime < b.icu_intime + interval '24 hour'
         AND e.ck_cpk IS NOT NULL
        WHERE b.first_icu_stay
        ORDER BY b.stay_id, e.charttime
        """,
        conn,
    )
    vent = pd.read_sql_query(
        """
        SELECT DISTINCT d.stay_id
        FROM mimiciv_derived.icustay_detail d
        JOIN mimiciv_derived.ventilation v ON v.stay_id = d.stay_id
        WHERE d.first_icu_stay
          AND v.ventilation_status = 'InvasiveVent'
          AND v.starttime < d.icu_intime + interval '24 hour'
          AND v.endtime > d.icu_intime
        """,
        conn,
    )
    vaso = pd.read_sql_query(
        """
        SELECT DISTINCT d.stay_id
        FROM mimiciv_derived.icustay_detail d
        JOIN mimiciv_derived.norepinephrine n ON n.stay_id = d.stay_id
        WHERE d.first_icu_stay
          AND n.starttime < d.icu_intime + interval '24 hour'
          AND n.endtime > d.icu_intime
        """,
        conn,
    )
    rrt = pd.read_sql_query(
        """
        SELECT DISTINCT b.stay_id
        FROM mimiciv_derived.icustay_detail b
        JOIN mimiciv_derived.rrt r ON r.stay_id = b.stay_id
        WHERE b.first_icu_stay
          AND r.charttime >= b.icu_intime
          AND r.charttime < b.icu_intime + interval '24 hour'
          AND r.dialysis_active = 1
        """,
        conn,
    )
    conn.close()
    out = df.merge(charlson, on="hadm_id", how="left").merge(weight, on="stay_id", how="left")
    out = out.merge(albumin, on="stay_id", how="left").merge(ck, on="stay_id", how="left")
    out["vent"] = out["stay_id"].isin(set(vent["stay_id"])).astype(float)
    out["vaso"] = out["stay_id"].isin(set(vaso["stay_id"])).astype(float)
    out["rrt"] = out["stay_id"].isin(set(rrt["stay_id"])).astype(float)
    out["aud"] = out["aud"].fillna(0).astype(float)
    out["log_ck"] = np.log(pd.to_numeric(out["ck_first"], errors="coerce").clip(lower=1))
    return out


def pack(result):
    return {
        "n": result["n"],
        "events": result["events"],
        "terms": result["terms"],
    }


def rubin(estimates, variances):
    estimates = np.asarray(estimates, float)
    variances = np.asarray(variances, float)
    m = len(estimates)
    qbar = float(estimates.mean())
    ubar = float(variances.mean())
    between = float(estimates.var(ddof=1)) if m > 1 else 0.0
    total = ubar + (1 + 1 / m) * between
    se = float(np.sqrt(total))
    z = qbar / se if se else np.nan
    from scipy import stats

    return {
        "hr": float(np.exp(qbar)),
        "hr_low": float(np.exp(qbar - 1.96 * se)),
        "hr_high": float(np.exp(qbar + 1.96 * se)),
        "p": float(2 * stats.norm.sf(abs(z))) if se else None,
    }


def multiple_imputation(df, extra_cols, m=30, split=False):
    from lifelines import NelsonAalenFitter

    d = df.copy()
    d["q2"] = (d["quartile"] == "Q2").astype(float)
    d["q3"] = (d["quartile"] == "Q3").astype(float)
    d["q4"] = (d["quartile"] == "Q4").astype(float)
    for unit in ["CCU", "SICU", "CVICU", "Other"]:
        d[f"u_{unit}"] = (d["unit"] == unit).astype(float)
    use_cols = [
        "age", "male", "sofa", "lactate", "log_bili", "creatinine", "aud",
        "u_CCU", "u_SICU", "u_CVICU", "u_Other",
    ] + extra_cols
    work = d.dropna(subset=["time_28", "event_28", "quartile"]).copy()
    naf = NelsonAalenFitter()
    naf.fit(work["time_28"], work["event_28"])
    work["na_hazard"] = np.asarray(naf.cumulative_hazard_at_times(work["time_28"]).to_numpy(), float).ravel()
    imp_cols = use_cols + ["event_28", "na_hazard"]
    base = work[imp_cols]
    store = {term: {"coef": [], "var": []} for term in TERMS}
    early_store = {term: {"coef": [], "var": []} for term in TERMS}
    late_store = {term: {"coef": [], "var": []} for term in TERMS}
    for seed in range(m):
        if seed % 5 == 0:
            print(f"impute {seed + 1}/{m}", flush=True)
        imputer = IterativeImputer(random_state=seed, max_iter=12, sample_posterior=True)
        filled = pd.DataFrame(imputer.fit_transform(base), columns=imp_cols, index=base.index)
        for col in ["male", "aud", "vaso", "vent", "rrt", "u_CCU", "u_SICU", "u_CVICU", "u_Other"]:
            if col in filled:
                filled[col] = (filled[col] >= 0.5).astype(float)
        model_df = filled[use_cols].copy()
        model_df["time_28"] = work["time_28"]
        model_df["event_28"] = work["event_28"]
        model_df["q2"] = work["q2"]
        model_df["q3"] = work["q3"]
        model_df["q4"] = work["q4"]
        formula = "q2 + q3 + q4 + " + " + ".join(use_cols) + " - 1"
        _fitted, res, X, _use = cox_fit(model_df, formula, TERMS)
        _collect(store, res, X)
        if split:
            early = model_df.copy()
            early["time"] = np.minimum(early["time_28"], 7.0)
            early["event"] = ((early["event_28"] == 1) & (early["time_28"] <= 7)).astype(int)
            _f, res_e, Xe, _ = cox_fit(early, formula, TERMS, time_col="time", event_col="event")
            _collect(early_store, res_e, Xe)
            late = model_df[model_df["time_28"] > 7].copy()
            late["time"] = late["time_28"] - 7.0
            late["event"] = late["event_28"]
            _f, res_l, Xl, _ = cox_fit(late, formula, TERMS, time_col="time", event_col="event")
            _collect(late_store, res_l, Xl)
    out = {
        "m": m,
        "n": int(len(work)),
        "events": int(work["event_28"].sum()),
        "included_outcome": True,
        "terms": {term: rubin(store[term]["coef"], store[term]["var"]) for term in TERMS},
    }
    if split:
        out["day_0_7"] = {term: rubin(early_store[term]["coef"], early_store[term]["var"]) for term in TERMS}
        out["day_7_28"] = {term: rubin(late_store[term]["coef"], late_store[term]["var"]) for term in TERMS}
    return out


def _collect(store, res, X):
    cov = np.asarray(res.cov_params(), float)
    names = list(X.columns)
    for term in TERMS:
        idx = names.index(term)
        store[term]["coef"].append(float(res.params[idx]))
        store[term]["var"].append(float(cov[idx, idx]))


def interaction(df, flag, name):
    part = df[df["quartile"].isin(["Q1", "Q4"])].copy()
    part["q4"] = (part["quartile"] == "Q4").astype(float)
    part["flag"] = part[flag].astype(float)
    formula = "q4 + flag + q4:flag + age + sofa + lactate + log_bili + creatinine + aud - 1"
    if name not in ("Male", "Female"):
        formula = "q4 + flag + q4:flag + age + male + sofa + lactate + log_bili + creatinine + aud - 1"
    fitted, _, X, _ = cox_fit(part, formula, ["q4", "q4:flag"])
    return {"subgroup": name, "n": fitted["n"], "events": fitted["events"], "terms": fitted["terms"]}


def eicu_with_unit(eicu):
    import statsmodels.formula.api as smf

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
    conn = psycopg2.connect(host="localhost", port=5442, user="postgres", password=os.environ["PGPASSWORD"], dbname="eicu")
    conn.autocommit = True
    with conn.cursor() as cur:
        cur.execute("SET statement_timeout = 0")
    alcohol = pd.read_sql_query(
        """
        SELECT DISTINCT patientunitstayid
        FROM diagnosis
        WHERE lower(diagnosisstring) LIKE '%alcohol%'
           OR icd9code LIKE '303%' OR icd9code LIKE '305.0%' OR icd9code LIKE '291%'
        """,
        conn,
    )
    conn.close()
    d["aud"] = d["patientunitstayid"].isin(set(alcohol["patientunitstayid"])).astype(float)
    need = ["y", "quartile", "age", "male", "apachescore", "lactate", "log_bili", "creatinine", "unit", "aud"]
    d = d.dropna(subset=need).copy()
    d["q"] = pd.Categorical(d["quartile"], ["Q1", "Q2", "Q3", "Q4"])
    model = smf.logit(
        "y ~ C(q, Treatment(reference='Q1')) + age + male + apachescore + lactate + log_bili + creatinine + aud + C(unit)",
        data=d,
    ).fit(disp=False)
    conf = model.conf_int()
    terms = {}
    for term in model.params.index:
        if not str(term).startswith("C(q"):
            continue
        terms[str(term)] = {
            "or": float(np.exp(model.params[term])),
            "or_low": float(np.exp(conf.loc[term, 0])),
            "or_high": float(np.exp(conf.loc[term, 1])),
            "p": float(model.pvalues[term]),
        }
    no_aud = d[d["aud"] == 0]
    model2 = smf.logit(
        "y ~ C(q, Treatment(reference='Q1')) + age + male + apachescore + lactate + log_bili + creatinine + C(unit)",
        data=no_aud,
    ).fit(disp=False)
    conf2 = model2.conf_int()
    sens = {}
    for term in model2.params.index:
        if "Q4" not in str(term):
            continue
        sens = {
            "or": float(np.exp(model2.params[term])),
            "or_low": float(np.exp(conf2.loc[term, 0])),
            "or_high": float(np.exp(conf2.loc[term, 1])),
            "p": float(model2.pvalues[term]),
            "n": int(len(no_aud)),
            "events": int(no_aud["y"].sum()),
        }
    return {
        "n": int(len(d)),
        "events": int(d["y"].sum()),
        "aud_prevalence": float(d["aud"].mean()),
        "terms": terms,
        "exclude_aud_q4": sens,
    }


def main():
    mimic, cuts, n_before = load_mimic()
    mimic = add_covariates(mimic)
    mimic["q2"] = (mimic["quartile"] == "Q2").astype(float)
    mimic["q3"] = (mimic["quartile"] == "Q3").astype(float)
    mimic["q4"] = (mimic["quartile"] == "Q4").astype(float)
    missing = {
        col: float(mimic[col].isna().mean())
        for col in ["sofa", "lactate", "log_bili", "creatinine", "charlson", "weight", "albumin", "ck_first", "aud"]
    }
    primary, _, _, _ = cox_fit(mimic, CORE, TERMS)
    expanded_formula = (
        "q2 + q3 + q4 + age + male + sofa + lactate + log_bili + creatinine + aud + "
        "charlson + weight + albumin + vaso + vent + rrt + C(unit) - 1"
    )
    expanded, _, _, _ = cox_fit(mimic, expanded_formula, TERMS)
    ck_formula = expanded_formula.replace(" - 1", " + log_ck - 1")
    ck_model, _, _, _ = cox_fit(mimic, ck_formula, TERMS)
    no_liverish = mimic[(mimic["aud"] == 0) & (mimic["viral_hep"] == 0) & (mimic["fatty_liver"] == 0)]
    sens, _, _, _ = cox_fit(no_liverish, CORE.replace(" + aud", ""), TERMS)
    early = mimic.copy()
    early["time"] = np.minimum(early["time_28"], 7.0)
    early["event"] = ((early["event_28"] == 1) & (early["time_28"] <= 7)).astype(int)
    late = mimic[mimic["time_28"] > 7].copy()
    late["time"] = late["time_28"] - 7.0
    late["event"] = late["event_28"]
    split = {
        "day_0_7": pack(cox_fit(early, CORE, TERMS, time_col="time", event_col="event")[0]),
        "day_7_28": pack(cox_fit(late, CORE, TERMS, time_col="time", event_col="event")[0]),
    }
    mimic["male_flag"] = mimic["male"]
    mimic["age65"] = (mimic["age"] >= 65).astype(float)
    mimic["lact4"] = (mimic["lactate"] >= 4).astype(float)
    mimic["cvicu_flag"] = (mimic["unit"] == "CVICU").astype(float)
    interactions = [
        interaction(mimic, "male_flag", "Male"),
        interaction(mimic, "age65", "Age >= 65"),
        interaction(mimic, "lact4", "Lactate >= 4"),
        interaction(mimic, "cvicu_flag", "Cardiac surgical ICU"),
    ]
    print("imputation", flush=True)
    mi_core = multiple_imputation(mimic, [], m=30, split=True)
    mi_expanded = multiple_imputation(
        mimic, ["charlson", "weight", "albumin", "vaso", "vent", "rrt"], m=30, split=False
    )
    eicu = load_eicu(cuts)
    eicu_out = eicu_with_unit(eicu)
    crude = (
        mimic.groupby("quartile", observed=True)["event_28"]
        .agg(n="size", deaths="sum", rate="mean")
        .reset_index()
        .to_dict(orient="records")
    )
    payload = {
        "n_eligible": n_before,
        "n_with_ratio": int(len(mimic)),
        "cuts": cuts,
        "aud_n": int(mimic["aud"].sum()),
        "viral_n": int(mimic["viral_hep"].sum()),
        "fatty_n": int(mimic["fatty_liver"].sum()),
        "missing_fraction": missing,
        "crude": crude,
        "primary_plus_aud": pack(primary),
        "expanded": pack(expanded),
        "expanded_plus_ck": pack(ck_model),
        "exclude_aud_hepatitis_fatty": pack(sens),
        "time_split": split,
        "interactions": interactions,
        "mi_core": mi_core,
        "mi_expanded": mi_expanded,
        "eicu": eicu_out,
    }
    (OUT / "round2.json").write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    print(json.dumps(payload, indent=2, default=str)[:12000])


if __name__ == "__main__":
    main()
