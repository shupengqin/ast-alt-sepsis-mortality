"""Pre-submission checks: landmark, time cuts, hospital death, selection, missingness.

Writes analysis/deritis_paper/revision_results.json and supplementary tables.
The quartile cutpoints stay frozen.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import psycopg2
from lifelines import CoxTimeVaryingFitter, NelsonAalenFitter
from lifelines.utils import concordance_index
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from scipy import stats

from deritis_paper_analysis import (
    assign_quartiles,
    connect,
    cox_fit,
    fit_cox,
    load_eicu,
    load_mimic,
    prep_ratio,
)
from deritis_round2 import add_covariates, rubin

OUT = Path(Path(__file__).resolve().parent / "deritis_paper")
USE = [
    "age", "male", "sofa", "lactate", "log_bili", "creatinine", "aud",
    "u_CCU", "u_SICU", "u_CVICU", "u_Other",
]
FORMULA = "q2 + q3 + q4 + " + " + ".join(USE) + " - 1"
M = 30


def dump(payload):
    (OUT / "revision_results.json").write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8"
    )


def fmt_p(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return ""
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


def add_quartiles(df):
    d = df.copy()
    d["q2"] = (d["quartile"] == "Q2").astype(float)
    d["q3"] = (d["quartile"] == "Q3").astype(float)
    d["q4"] = (d["quartile"] == "Q4").astype(float)
    for unit in ("CCU", "SICU", "CVICU", "Other"):
        d[f"u_{unit}"] = (d["unit"] == unit).astype(float)
    return d


def absolute_risks(df):
    rows = []
    pieces = [("Overall", df)] + [(q, df[df["quartile"] == q]) for q in ("Q1", "Q2", "Q3", "Q4")]
    for name, sub in pieces:
        early = (sub["event_28"] == 1) & (sub["time_28"] <= 7)
        at_risk = sub["time_28"] > 7
        late = (sub["event_28"] == 1) & (sub["time_28"] > 7)
        rows.append(
            {
                "group": name,
                "n": int(len(sub)),
                "deaths_day_0_7": int(early.sum()),
                "pct_day_0_7": float(early.mean()),
                "alive_at_day_7": int(at_risk.sum()),
                "deaths_day_7_28": int(late.sum()),
                "pct_day_7_28_among_alive": float(late.sum() / at_risk.sum()) if at_risk.sum() else None,
            }
        )
    return rows


def measurement_breakdown(full):
    both = full["ast_raw"].notna() & full["alt_raw"].notna()
    usable = full["ratio"].notna()
    return {
        "eligible": int(len(full)),
        "both_enzymes_in_24h": int(both.sum()),
        "usable_after_bounds": int(usable.sum()),
        "both_but_outside_bounds": int((both & ~usable).sum()),
        "no_paired_result": int((~both).sum()),
    }


def reason_query():
    conn = connect("mimiciv31")
    tab = pd.read_sql_query(
        """
        WITH liver AS (
            SELECT hadm_id,
                MAX(CASE
                    WHEN icd_version = 9 AND replace(icd_code,'.','') ~ '^(5712|5715|5716|5718|5722)' THEN 1
                    WHEN icd_version = 10 AND replace(icd_code,'.','') ~ '^(K70|K74|K717)' THEN 1
                    ELSE 0 END) AS chronic_liver
            FROM mimiciv_hosp.diagnoses_icd GROUP BY hadm_id
        ),
        eligible AS (
            SELECT d.stay_id, d.hadm_id, d.icu_intime
            FROM mimiciv_derived.icustay_detail d
            JOIN mimiciv_derived.sepsis3 s ON s.stay_id = d.stay_id AND s.sepsis3
            LEFT JOIN liver l ON l.hadm_id = d.hadm_id
            WHERE d.first_icu_stay AND d.admission_age >= 18
              AND COALESCE(l.chronic_liver, 0) = 0
        ),
        flags AS (
            SELECT e2.stay_id,
                MAX(CASE WHEN z.ast IS NOT NULL THEN 1 ELSE 0 END) AS any_ast,
                MAX(CASE WHEN z.alt IS NOT NULL THEN 1 ELSE 0 END) AS any_alt
            FROM eligible e2
            LEFT JOIN mimiciv_derived.enzyme z
              ON z.hadm_id = e2.hadm_id
             AND z.charttime >= e2.icu_intime
             AND z.charttime < e2.icu_intime + interval '24 hour'
            GROUP BY e2.stay_id
        )
        SELECT
            SUM(CASE WHEN any_ast = 0 AND any_alt = 0 THEN 1 ELSE 0 END) AS neither,
            SUM(CASE WHEN any_ast = 1 AND any_alt = 0 THEN 1 ELSE 0 END) AS ast_only,
            SUM(CASE WHEN any_ast = 0 AND any_alt = 1 THEN 1 ELSE 0 END) AS alt_only,
            SUM(CASE WHEN any_ast = 1 AND any_alt = 1 THEN 1 ELSE 0 END) AS both_present
        FROM flags
        """,
        conn,
    )
    conn.close()
    return {k: int(tab.loc[0, k]) for k in tab.columns}


def attach_hospital(full):
    conn = connect("mimiciv31")
    adm = pd.read_sql_query(
        """
        SELECT hadm_id, hospital_expire_flag, deathtime
        FROM mimiciv_hosp.admissions
        """,
        conn,
    )
    conn.close()
    adm["deathtime"] = pd.to_datetime(adm["deathtime"], errors="coerce")
    out = full.merge(adm, on="hadm_id", how="left")
    out["hospital_death"] = out["hospital_expire_flag"].fillna(0).astype(int)
    return out


def landmark_frame(df):
    death = df["deathtime"].where(df["deathtime"].notna(), df["dod"])
    landmark = df["icu_intime"] + pd.Timedelta(hours=24)
    died_before = death.notna() & (death < landmark)
    keep = df.loc[~died_before].copy()
    death_k = death.loc[keep.index]
    landmark_k = landmark.loc[keep.index]
    follow = (death_k - landmark_k).dt.total_seconds() / 86400.0
    event = death_k.notna() & (follow >= 0) & (follow <= 28)
    keep["time_lm"] = np.where(event, np.maximum(follow, 0.5), 28.0)
    keep["event_lm"] = event.astype(int)
    return keep, int(died_before.sum())


def period_fit(model_df, cut, time_col="time_28", event_col="event_28"):
    base = model_df.copy()
    base["id"] = np.arange(len(base))
    early = base.copy()
    early["start"] = 0.0
    early["stop"] = np.minimum(early[time_col], float(cut))
    early["event"] = ((early[event_col] == 1) & (early[time_col] <= cut)).astype(int)
    early["late"] = 0.0
    late = base[base[time_col] > cut].copy()
    late["start"] = float(cut)
    late["stop"] = late[time_col]
    late["event"] = late[event_col].astype(int)
    late["late"] = 1.0
    both = pd.concat([early, late], ignore_index=True)
    both = both[both["stop"] > both["start"]]
    for q in ("q2", "q3", "q4"):
        both[f"{q}_late"] = both[q] * both["late"]
    cols = [
        "id", "start", "stop", "event", "q2", "q3", "q4",
        "q2_late", "q3_late", "q4_late",
    ] + USE
    fit = CoxTimeVaryingFitter()
    fit.fit(both[cols], id_col="id", event_col="event", start_col="start", stop_col="stop", show_progress=False)
    names = list(fit.params_.index)
    var = np.asarray(fit.variance_matrix_, float)
    out = {}
    for q in ("q2", "q3", "q4"):
        i_early = names.index(q)
        i_ratio = names.index(f"{q}_late")
        coef_e = float(fit.params_[q])
        coef_r = float(fit.params_[f"{q}_late"])
        var_e = float(var[i_early, i_early])
        var_r = float(var[i_ratio, i_ratio])
        cov = float(var[i_early, i_ratio])
        out[f"{q}_early"] = (coef_e, var_e)
        out[f"{q}_late_ratio"] = (coef_r, var_r)
        out[f"{q}_after"] = (coef_e + coef_r, var_e + var_r + 2 * cov)
    return out


def pooled_from_store(store):
    return {term: rubin(store[term]["coef"], store[term]["var"]) for term in store}


def new_store(keys):
    return {k: {"coef": [], "var": []} for k in keys}


def push(store, mapping):
    for term, (coef, var) in mapping.items():
        store[term]["coef"].append(coef)
        store[term]["var"].append(var)


def cox_terms(df, formula, time_col="time_28", event_col="event_28"):
    fitted, res, X, use = cox_fit(df, formula, [], time_col=time_col, event_col=event_col)
    cov = np.asarray(res.cov_params(), float)
    found = {}
    for i, term in enumerate(X.columns):
        coef = float(res.params[i])
        se = float(np.sqrt(cov[i, i]))
        found[term] = {
            "hr": float(np.exp(coef)),
            "hr_low": float(np.exp(coef - 1.96 * se)),
            "hr_high": float(np.exp(coef + 1.96 * se)),
            "p": float(2 * stats.norm.sf(abs(coef / se))) if se else None,
            "coef": coef,
            "var": float(cov[i, i]),
        }
    return {"n": fitted["n"], "events": int(use[event_col].sum()), "terms": found}, res, X, use


def c_index_pair(df):
    base = " + ".join(USE) + " - 1"
    ext = "log_ratio_z + " + base
    _, res_b, Xb, use_b = cox_terms(df, base)
    same = df.loc[use_b.index]
    _, res_e, Xe, use_e = cox_terms(same, ext)
    use_b = use_b.loc[use_e.index]
    Xb = Xb.loc[use_e.index]
    def c_of(res, X, use, time_col, event_col):
        lp = np.asarray(X, float) @ np.asarray(res.params, float)
        return float(concordance_index(use[time_col], -lp, use[event_col]))
    out = {}
    for label, time_col, event_col in (
        ("day_28", "time_28", "event_28"),
        ("day_7", "time_7", "event_7"),
    ):
        part = use_e.copy()
        part["time_7"] = np.minimum(part["time_28"], 7.0)
        part["event_7"] = ((part["event_28"] == 1) & (part["time_28"] <= 7)).astype(int)
        # rebuild linear predictors on this frame's original columns
        cb = c_of(res_b, Xb.loc[part.index], part, "time_28" if label == "day_28" else "time_7", "event_28" if label == "day_28" else "event_7")
        ce = c_of(res_e, Xe.loc[part.index], part, "time_28" if label == "day_28" else "time_7", "event_28" if label == "day_28" else "event_7")
        out[label] = {"base": cb, "plus_ratio": ce, "delta": ce - cb}
    return out


def prepare_imputation(df, time_col, event_col):
    work = df.dropna(subset=[time_col, event_col, "quartile"]).copy()
    naf = NelsonAalenFitter()
    naf.fit(work[time_col], work[event_col])
    work["na_hazard"] = np.asarray(naf.cumulative_hazard_at_times(work[time_col]).to_numpy(), float).ravel()
    work["y_imp"] = work[event_col].astype(float)
    return work


def fill_once(work, imp_cols, seed):
    imputer = IterativeImputer(random_state=seed, max_iter=12, sample_posterior=True)
    filled = pd.DataFrame(imputer.fit_transform(work[imp_cols]), columns=imp_cols, index=work.index)
    for col in ("male", "aud", "u_CCU", "u_SICU", "u_CVICU", "u_Other"):
        if col in filled:
            filled[col] = (filled[col] >= 0.5).astype(float)
    model = filled[USE].copy()
    for col in ("q2", "q3", "q4", "log_ast_z", "log_alt_z", "log_ratio_z"):
        if col in work:
            model[col] = work[col]
    return model


def threshold_models(full, cuts):
    specs = {
        "primary_5_5000_and_5_2000": (5, 5000, 5, 2000),
        "wide_1_10000": (1, 10000, 1, 10000),
        "tight_10_1000": (10, 1000, 10, 1000),
    }
    out = {}
    masks = {
        name: full["ast_raw"].between(alo, ahi) & full["alt_raw"].between(blo, bhi)
        for name, (alo, ahi, blo, bhi) in specs.items()
    }
    masks["positive_only"] = full["ast_raw"].gt(0) & full["alt_raw"].gt(0)
    for name, ok in masks.items():
        part = full.loc[ok].copy()
        part["ratio"] = part["ast_raw"] / part["alt_raw"]
        part["quartile"] = assign_quartiles(part["ratio"], cuts)
        part["log_ratio"] = np.log(part["ratio"])
        part["log_ratio_z"] = (part["log_ratio"] - part["log_ratio"].mean()) / part["log_ratio"].std(ddof=0)
        part = add_quartiles(part)
        fitted, _, _, _ = cox_terms(part, FORMULA)
        out[name] = {
            "n": fitted["n"],
            "events": fitted["events"],
            "q4": fitted["terms"].get("q4"),
            "quartile_n": {q: int((part["quartile"] == q).sum()) for q in ("Q1", "Q2", "Q3", "Q4")},
        }
        print(name, out[name]["n"], out[name]["q4"], flush=True)
    return out


def ipw(full, cuts):
    d = full.copy()
    d["measured"] = d["ratio"].notna().astype(int)
    for col, med_name in (("sofa", "sofa"), ("lactate", "lactate"), ("creatinine", "creatinine"), ("log_bili", "log_bili")):
        d[f"{col}_miss"] = d[col].isna().astype(float)
        d[col] = d[col].fillna(d[col].median())
    d["aud"] = d["aud"].fillna(0).astype(float)
    for unit in ("CCU", "SICU", "CVICU", "Other"):
        d[f"u_{unit}"] = (d["unit"] == unit).astype(float)
    import statsmodels.api as sm
    xcols = ["age", "male", "sofa", "sofa_miss", "lactate", "lactate_miss", "creatinine", "creatinine_miss", "log_bili", "log_bili_miss", "aud", "u_CCU", "u_SICU", "u_CVICU", "u_Other"]
    design = sm.add_constant(d[xcols], has_constant="add")
    keep = [col for col in design.columns if col == "const" or float(design[col].std()) > 1e-8]
    design = design[keep]
    fit = sm.Logit(d["measured"], design).fit(disp=False, maxiter=200, method="bfgs")
    p = np.clip(fit.predict(design), 0.05, 0.95)
    d["p_sel"] = p
    marg = float(d["measured"].mean())
    measured = d[d["measured"] == 1].copy()
    measured["w"] = marg / measured["p_sel"]
    lo, hi = measured["w"].quantile(0.01), measured["w"].quantile(0.99)
    measured["w"] = measured["w"].clip(lo, hi)
    need = ["time_28", "event_28", "q2", "q3", "q4", "age", "male", "sofa", "lactate", "log_bili", "creatinine", "aud", "w"] + [f"u_{u}" for u in ("CCU", "SICU", "CVICU", "Other")]
    # restore observed covariates for the outcome model; do not use the filled values
    observed = full.loc[measured.index, ["sofa", "lactate", "log_bili", "creatinine", "ratio", "time_28", "event_28", "age", "male", "aud", "unit"]].copy()
    observed["quartile"] = assign_quartiles(observed["ratio"], cuts)
    observed = add_quartiles(observed)
    observed["w"] = measured["w"]
    cc = observed.dropna(subset=need).copy()
    from lifelines import CoxPHFitter
    cols = ["time_28", "event_28", "q2", "q3", "q4"] + USE + ["w"]
    weighted = CoxPHFitter()
    weighted.fit(cc[cols], "time_28", "event_28", weights_col="w", robust=True, show_progress=False)
    plain = CoxPHFitter()
    plain.fit(cc[[c for c in cols if c != "w"]], "time_28", "event_28", show_progress=False)
    def grab(model):
        terms = {}
        for term in ("q2", "q3", "q4"):
            coef = float(model.params_[term])
            se = float(model.standard_errors_[term])
            terms[term] = {
                "hr": float(np.exp(coef)),
                "hr_low": float(np.exp(coef - 1.96 * se)),
                "hr_high": float(np.exp(coef + 1.96 * se)),
                "p": float(2 * stats.norm.sf(abs(coef / se))),
            }
        return terms
    ess = float((cc["w"].sum() ** 2) / (cc["w"] ** 2).sum())
    return {
        "selection_n": int(len(d)),
        "selection_measured": int(d["measured"].sum()),
        "outcome_n": int(len(cc)),
        "outcome_events": int(cc["event_28"].sum()),
        "weight_median": float(cc["w"].median()),
        "weight_p99": float(hi),
        "ess": ess,
        "weighted": grab(weighted),
        "same_rows_unweighted": grab(plain),
    }


def missing_by_quartile(df, cols):
    rows = []
    pieces = [("Overall", df)] + [(q, df[df["quartile"] == q]) for q in ("Q1", "Q2", "Q3", "Q4")]
    for col in cols:
        for name, sub in pieces:
            s = pd.to_numeric(sub[col], errors="coerce")
            rows.append(
                {
                    "variable": col,
                    "group": name,
                    "n": int(len(sub)),
                    "n_missing": int(s.isna().sum()),
                    "pct_missing": float(s.isna().mean()),
                    "observed_median": None if s.dropna().empty else float(s.median()),
                    "observed_q1": None if s.dropna().empty else float(s.quantile(0.25)),
                    "observed_q3": None if s.dropna().empty else float(s.quantile(0.75)),
                }
            )
    return rows


def spline_grid(df, exposure):
    from patsy import build_design_matrices, dmatrix
    formula = f"cr({exposure}, df=4) + " + " + ".join(USE) + " - 1"
    X = dmatrix(formula, df, return_type="dataframe")
    info = X.design_info
    use = df.loc[X.index].dropna(subset=["time_28", "event_28", exposure])
    X = X.loc[use.index]
    res = fit_cox(use["time_28"], use["event_28"], X)
    grid = np.linspace(use[exposure].quantile(0.02), use[exposure].quantile(0.98), 40)
    ref = float(use[exposure].median())
    base = use.iloc[0:1]
    pred = pd.concat([base.assign(**{exposure: v}) for v in list(grid) + [ref]], ignore_index=True)
    mat = build_design_matrices([info], pred, return_type="dataframe")[0]
    mat = mat.reindex(columns=list(info.column_names), fill_value=0.0)
    beta = np.asarray(res.params, float)
    covm = np.asarray(res.cov_params(), float)
    ref_x = np.asarray(mat.iloc[-1], float)
    rows = []
    for i, value in enumerate(grid):
        delta = np.asarray(mat.iloc[i], float) - ref_x
        est = float(delta @ beta)
        se = float(np.sqrt(max(float(delta @ covm @ delta), 0)))
        rows.append(
            {
                "exposure": exposure,
                "value": float(np.exp(value)),
                "hr": float(np.exp(est)),
                "hr_low": float(np.exp(est - 1.96 * se)),
                "hr_high": float(np.exp(est + 1.96 * se)),
            }
        )
    return rows


def mimic_mi(work):
    imp_cols = ["q2", "q3", "q4"] + USE + ["y_imp", "na_hazard"]
    keys = ["q2", "q3", "q4", "log_ast_z", "log_alt_z"] + USE
    store = new_store(keys)
    period_keys = []
    for cut in (3, 7, 14):
        for q in ("q2", "q3", "q4"):
            period_keys += [f"c{cut}_{q}_early", f"c{cut}_{q}_late_ratio", f"c{cut}_{q}_after"]
    period_store = new_store(period_keys)
    c_rows = []
    lactate_medians = []
    joint = "log_ast_z + log_alt_z + " + " + ".join(USE) + " - 1"
    for seed in range(M):
        print(f"mimic mi {seed + 1}/{M}", flush=True)
        model = fill_once(work, imp_cols, seed)
        model["time_28"] = work["time_28"]
        model["event_28"] = work["event_28"]
        fitted, _, _, _ = cox_terms(model, FORMULA)
        for term, rec in fitted["terms"].items():
            if term in store:
                store[term]["coef"].append(rec["coef"])
                store[term]["var"].append(rec["var"])
        joint_fit, _, _, _ = cox_terms(model, joint)
        for term in ("log_ast_z", "log_alt_z"):
            rec = joint_fit["terms"][term]
            store[term]["coef"].append(rec["coef"])
            store[term]["var"].append(rec["var"])
        for cut in (3, 7, 14):
            fitted_p = period_fit(model, cut)
            mapped = {f"c{cut}_{k}": v for k, v in fitted_p.items()}
            push(period_store, mapped)
        model["log_ratio_z"] = work["log_ratio_z"]
        c_rows.append(c_index_pair(model))
        lactate_medians.append(float(model["lactate"].median()))
    return {
        "n": int(len(work)),
        "events": int(work["event_28"].sum()),
        "terms": pooled_from_store(store),
        "periods": pooled_from_store(period_store),
        "c_index": {
            "day_28_base_mean": float(np.mean([r["day_28"]["base"] for r in c_rows])),
            "day_28_plus_mean": float(np.mean([r["day_28"]["plus_ratio"] for r in c_rows])),
            "day_28_delta_mean": float(np.mean([r["day_28"]["delta"] for r in c_rows])),
            "day_28_delta_min": float(np.min([r["day_28"]["delta"] for r in c_rows])),
            "day_28_delta_max": float(np.max([r["day_28"]["delta"] for r in c_rows])),
            "day_7_base_mean": float(np.mean([r["day_7"]["base"] for r in c_rows])),
            "day_7_plus_mean": float(np.mean([r["day_7"]["plus_ratio"] for r in c_rows])),
            "day_7_delta_mean": float(np.mean([r["day_7"]["delta"] for r in c_rows])),
            "day_7_delta_min": float(np.min([r["day_7"]["delta"] for r in c_rows])),
            "day_7_delta_max": float(np.max([r["day_7"]["delta"] for r in c_rows])),
        },
        "imputed_lactate_median_mean": float(np.mean(lactate_medians)),
    }


def landmark_mi(lm):
    work = prepare_imputation(lm, "time_lm", "event_lm")
    imp_cols = ["q2", "q3", "q4"] + USE + ["y_imp", "na_hazard"]
    store = new_store(["q2", "q3", "q4"])
    for seed in range(M):
        print(f"landmark mi {seed + 1}/{M}", flush=True)
        model = fill_once(work, imp_cols, seed)
        model["time_lm"] = work["time_lm"]
        model["event_lm"] = work["event_lm"]
        fitted, _, _, _ = cox_terms(model, FORMULA, time_col="time_lm", event_col="event_lm")
        for term in ("q2", "q3", "q4"):
            store[term]["coef"].append(fitted["terms"][term]["coef"])
            store[term]["var"].append(fitted["terms"][term]["var"])
    return {"n": int(len(work)), "events": int(work["event_lm"].sum()), "terms": pooled_from_store(store)}


def hospital_mi(df):
    d = df.dropna(subset=["hospital_death", "quartile"]).copy()
    d["y_imp"] = d["hospital_death"].astype(float)
    d["na_hazard"] = 0.0
    imp_cols = USE + ["y_imp"]
    store = new_store(["q2", "q3", "q4"])
    import statsmodels.api as sm
    for seed in range(M):
        print(f"hospital mi {seed + 1}/{M}", flush=True)
        model = fill_once(d, imp_cols, seed)
        X = model[USE + ["q2", "q3", "q4"]].copy()
        X = sm.add_constant(X, has_constant="add")
        res = sm.Logit(d["hospital_death"], X).fit(disp=False, maxiter=200)
        cov = np.asarray(res.cov_params(), float)
        names = list(X.columns)
        for term in ("q2", "q3", "q4"):
            idx = names.index(term)
            store[term]["coef"].append(float(res.params[term]))
            store[term]["var"].append(float(cov[idx, idx]))
    crude = []
    for q in ("Q1", "Q2", "Q3", "Q4"):
        sub = d[d["quartile"] == q]
        crude.append({"quartile": q, "n": int(len(sub)), "deaths": int(sub["hospital_death"].sum()), "rate": float(sub["hospital_death"].mean())})
    return {
        "n": int(len(d)),
        "events": int(d["hospital_death"].sum()),
        "crude": crude,
        "terms": pooled_from_store(store),
    }


def eicu_prepare(eicu):
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
    conn = connect("eicu")
    alcohol = pd.read_sql_query(
        """
        SELECT DISTINCT patientunitstayid FROM diagnosis
        WHERE lower(diagnosisstring) LIKE '%alcohol%'
           OR icd9code LIKE '303%' OR icd9code LIKE '305.0%' OR icd9code LIKE '291%'
        """,
        conn,
    )
    hosp = pd.read_sql_query(
        "SELECT patientunitstayid, hospitalid FROM patient",
        conn,
    )
    conn.close()
    d["aud"] = d["patientunitstayid"].isin(set(alcohol["patientunitstayid"])).astype(float)
    d = d.merge(hosp.drop_duplicates("patientunitstayid"), on="patientunitstayid", how="left")
    for name, label in (("u_MICU", "MICU"), ("u_CCU", "CCU"), ("u_SICU", "SICU"), ("u_Cardiac", "Cardiac")):
        d[name] = (d["unit"] == label).astype(float)
    d["q2"] = (d["quartile"] == "Q2").astype(float)
    d["q3"] = (d["quartile"] == "Q3").astype(float)
    d["q4"] = (d["quartile"] == "Q4").astype(float)
    return d


def eicu_cluster_mi(eicu):
    import statsmodels.api as sm
    d = eicu.dropna(subset=["y", "quartile", "hospitalid"]).copy()
    euse = ["age", "male", "apachescore", "lactate", "log_bili", "creatinine", "aud", "u_MICU", "u_CCU", "u_SICU", "u_Cardiac"]
    imp_cols = euse + ["y"]
    store = new_store(["q2", "q3", "q4"])
    plain = new_store(["q2", "q3", "q4"])
    for seed in range(M):
        print(f"eicu cluster {seed + 1}/{M}", flush=True)
        imputer = IterativeImputer(random_state=seed, max_iter=12, sample_posterior=True)
        filled = pd.DataFrame(imputer.fit_transform(d[imp_cols]), columns=imp_cols, index=d.index)
        for col in ("male", "y", "aud", "u_MICU", "u_CCU", "u_SICU", "u_Cardiac"):
            filled[col] = (filled[col] >= 0.5).astype(float)
        X = filled[euse].copy()
        X["q2"] = d["q2"]
        X["q3"] = d["q3"]
        X["q4"] = d["q4"]
        X = sm.add_constant(X, has_constant="add")
        keep = [c for c in X.columns if float(X[c].std()) > 1e-8 or c == "const"]
        X = X[keep]
        y = d["y"]
        groups = d["hospitalid"]
        try:
            res_c = sm.Logit(y, X).fit(disp=False, maxiter=200, cov_type="cluster", cov_kwds={"groups": groups})
        except Exception:
            X = X.drop(columns=[c for c in X.columns if str(c).startswith("u_")], errors="ignore")
            res_c = sm.Logit(y, X).fit(disp=False, maxiter=200, cov_type="cluster", cov_kwds={"groups": groups})
        res_p = sm.Logit(y, X).fit(disp=False, maxiter=200)
        for target, res in ((store, res_c), (plain, res_p)):
            cov = np.asarray(res.cov_params(), float)
            names = list(X.columns)
            for term in ("q2", "q3", "q4"):
                idx = names.index(term)
                target[term]["coef"].append(float(np.asarray(res.params)[idx]))
                target[term]["var"].append(float(cov[idx, idx]))
    return {
        "n": int(len(d)),
        "events": int(d["y"].sum()),
        "hospitals": int(d["hospitalid"].nunique()),
        "cluster": pooled_from_store(store),
        "independent": pooled_from_store(plain),
    }


def main():
    import pickle
    cache = OUT / "_revision_cache.pkl"
    print("load", flush=True)
    if cache.exists():
        mimic, cuts, n_eligible, full = pickle.loads(cache.read_bytes())
        print("cache", len(mimic), len(full), flush=True)
    else:
        mimic, cuts, n_eligible, full = load_mimic(return_frame=True)
        mimic = add_covariates(mimic)
        full = attach_hospital(full)
        mimic = mimic.merge(
            full[["stay_id", "hospital_death", "deathtime", "hospital_expire_flag"]],
            on="stay_id",
            how="left",
        )
        mimic = add_quartiles(mimic)
        cache.write_bytes(pickle.dumps((mimic, cuts, n_eligible, full), protocol=4))
    payload = {"cuts": cuts, "n_eligible": n_eligible}
    payload["measurement"] = measurement_breakdown(full)
    print("reasons", flush=True)
    payload["measurement"].update(reason_query())
    payload["absolute_risk"] = absolute_risks(mimic)
    miss_cols = ["lactate", "log_bili", "creatinine", "sofa", "albumin", "weight", "charlson"]
    payload["missing_by_quartile"] = missing_by_quartile(mimic, [c for c in miss_cols if c in mimic.columns])
    pd.DataFrame(payload["missing_by_quartile"]).to_csv(OUT / "supplement_missingness.csv", index=False)
    pd.DataFrame(payload["absolute_risk"]).to_csv(OUT / "supplement_absolute_risk.csv", index=False)
    print("threshold and ipw", flush=True)
    payload["threshold"] = threshold_models(full, cuts)
    payload["ipw"] = ipw(full, cuts)
    print("landmark cohort", flush=True)
    lm, n_before = landmark_frame(mimic)
    payload["landmark_excluded_died_within_24h"] = n_before
    lm_cc, _, _, _ = cox_terms(lm, FORMULA, time_col="time_lm", event_col="event_lm")
    payload["landmark_complete_case"] = {
        "n": lm_cc["n"],
        "events": lm_cc["events"],
        "q4": lm_cc["terms"].get("q4"),
        "q3": lm_cc["terms"].get("q3"),
        "q2": lm_cc["terms"].get("q2"),
    }
    print("complete-case periods", flush=True)
    cc = mimic.dropna(subset=["time_28", "event_28", "quartile"] + USE).copy()
    payload["period_complete_case"] = {}
    for cut in (3, 7, 14):
        fitted = period_fit(cc, cut)
        payload["period_complete_case"][str(cut)] = {
            k: {"hr": float(np.exp(v[0])), "hr_low": float(np.exp(v[0] - 1.96 * np.sqrt(v[1]))), "hr_high": float(np.exp(v[0] + 1.96 * np.sqrt(v[1])))}
            for k, v in fitted.items() if k.startswith("q4")
        }
        print(cut, payload["period_complete_case"][str(cut)], flush=True)
    print("splines", flush=True)
    spline_rows = spline_grid(mimic, "log_ast") + spline_grid(mimic, "log_alt")
    pd.DataFrame(spline_rows).to_csv(OUT / "supplement_enzyme_splines.csv", index=False)
    coef_rows = []
    primary_cc, _, _, _ = cox_terms(mimic, FORMULA)
    for term, rec in primary_cc["terms"].items():
        coef_rows.append({"model": "complete_case_28day", "term": term, **{k: rec[k] for k in ("hr", "hr_low", "hr_high", "p")}})
    pd.DataFrame(coef_rows).to_csv(OUT / "supplement_model_coefficients.csv", index=False)
    dump(payload)
    print("imputation", flush=True)
    work = prepare_imputation(mimic, "time_28", "event_28")
    payload["mimic_mi"] = mimic_mi(work)
    dump(payload)
    payload["landmark_mi"] = landmark_mi(lm)
    dump(payload)
    payload["hospital_mi"] = hospital_mi(mimic)
    dump(payload)
    print("eicu", flush=True)
    eicu = load_eicu(cuts)
    eicu = eicu_prepare(eicu)
    payload["eicu_missing_by_quartile"] = missing_by_quartile(
        eicu, ["lactate", "log_bili", "creatinine", "apachescore"]
    )
    pd.DataFrame(payload["eicu_missing_by_quartile"]).to_csv(OUT / "supplement_missingness_eicu.csv", index=False)
    payload["eicu_cluster"] = eicu_cluster_mi(eicu)
    dump(payload)
    print("done", flush=True)


if __name__ == "__main__":
    main()
