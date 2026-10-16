"""AST/ALT ratio and mortality in sepsis without chronic liver disease.

Development cohort: MIMIC-IV, all care units, including the cardiac surgical ICU.
Validation cohort: eICU, including cardiac and cardiothoracic units.
"""

from __future__ import annotations

import os
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psycopg2
from patsy import build_design_matrices, dmatrix
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg

OUT = Path(r"F:\CKB数据库\analysis\deritis_paper")
OUT.mkdir(parents=True, exist_ok=True)
COV = "age + male + sofa + lactate + log_bili + creatinine + C(unit)"
QTERMS = ["C(quartile)[T.Q2]", "C(quartile)[T.Q3]", "C(quartile)[T.Q4]"]


def _iqr_row(series):
    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        return ""
    return f"{s.median():.1f} ({s.quantile(0.25):.1f}–{s.quantile(0.75):.1f})"


def _write_missingness(df):
    rows = []
    for label, part in (
        ("Usable AST and ALT", df[df["ratio"].notna()]),
        ("No usable pair", df[df["ratio"].isna()]),
    ):
        rows.append(
            {
                "group": label,
                "n": int(len(part)),
                "age": _iqr_row(part["age"]),
                "male_n_pct": f"{int((part['male'] == 1).sum())} ({(part['male'] == 1).mean() * 100:.1f})",
                "sofa": _iqr_row(part["sofa"]),
                "lactate": _iqr_row(part["lactate"]),
                "creatinine": _iqr_row(part["creatinine"]),
                "cvicu_n_pct": f"{int(part['cvicu'].sum())} ({part['cvicu'].mean() * 100:.1f})",
                "death28_n_pct": f"{int(part['event_icu'].sum())} ({part['event_icu'].mean() * 100:.1f})",
            }
        )
    pd.DataFrame(rows).to_csv(OUT / "missingness_mimic.csv", index=False)


def connect(db):
    conn = psycopg2.connect(
        host="localhost", port=5442, user="postgres", password=os.environ["PGPASSWORD"], dbname=db
    )
    conn.autocommit = True
    with conn.cursor() as cur:
        cur.execute("SET statement_timeout = 0")
        cur.execute("SET work_mem = '256MB'")
    return conn


def prep_ratio(df):
    d = df.copy()
    ok = d["ast"].between(5, 5000) & d["alt"].between(5, 2000)
    d.loc[~ok, ["ast", "alt"]] = np.nan
    d["ratio"] = d["ast"] / d["alt"]
    d["log_ratio"] = np.log(d["ratio"])
    d["log_ast"] = np.log(d["ast"])
    d["log_bili"] = np.log(pd.to_numeric(d["bili"], errors="coerce").clip(lower=0.1))
    return d


def assign_quartiles(ratio, cuts):
    bins = [-np.inf, cuts[0], cuts[1], cuts[2], np.inf]
    return pd.cut(ratio, bins=bins, labels=["Q1", "Q2", "Q3", "Q4"], right=True)


def design(df, formula):
    return dmatrix(formula, df, return_type="dataframe")


def fit_cox(time, event, X):
    model = PHReg(np.asarray(time, float), np.asarray(X, float), status=np.asarray(event, int))
    return model.fit(disp=False)


def hr_of(res, names, term):
    idx = list(names).index(term)
    coef = float(res.params[idx])
    se = float(np.sqrt(np.asarray(res.cov_params())[idx, idx]))
    lo, hi = coef - 1.96 * se, coef + 1.96 * se
    return {
        "hr": float(np.exp(coef)),
        "hr_low": float(np.exp(lo)),
        "hr_high": float(np.exp(hi)),
        "p": float(2 * stats.norm.sf(abs(coef / se))) if se else None,
    }


def cox_fit(df, formula, terms, time_col="time_28", event_col="event_28"):
    X = design(df, formula)
    use = df.loc[X.index].dropna(subset=[time_col, event_col])
    X = X.loc[use.index]
    res = fit_cox(use[time_col], use[event_col], X)
    found = {}
    for term in terms:
        if term in X.columns:
            found[term] = hr_of(res, X.columns, term)
    return {"n": int(len(use)), "events": int(use[event_col].sum()), "terms": found, "llf": float(res.llf)}, res, X, use


def iqr(s):
    s = pd.to_numeric(s, errors="coerce").dropna()
    return f"{s.median():.1f} ({s.quantile(0.25):.1f}-{s.quantile(0.75):.1f})"


def npct(mask):
    m = pd.Series(mask).dropna()
    return f"{int(m.sum())} ({100 * m.mean():.1f}%)"


def km_curve(time, event):
    order = np.argsort(np.asarray(time, float))
    time = np.asarray(time, float)[order]
    event = np.asarray(event, int)[order]
    at_risk = len(time)
    surv = 1.0
    var = 0.0
    times, probs = [0.0], [1.0]
    i = 0
    while i < len(time):
        t = time[i]
        d = c = 0
        while i < len(time) and time[i] == t:
            d += int(event[i] == 1)
            c += int(event[i] == 0)
            i += 1
        if d and at_risk:
            surv *= 1 - d / at_risk
            if at_risk > d:
                var += d / (at_risk * (at_risk - d))
            times.append(float(t))
            probs.append(float(surv))
        at_risk -= d + c
    return np.array(times), np.array(probs)


def logrank_q4(df):
    part = df[df["quartile"].isin(["Q1", "Q4"])].copy()
    part["g"] = (part["quartile"] == "Q4").astype(int)
    o1 = e1 = v = 0.0
    for t in np.sort(part.loc[part["event_28"] == 1, "time_28"].unique()):
        risk = part[part["time_28"] >= t]
        deaths = part[(part["time_28"] == t) & (part["event_28"] == 1)]
        n, n1, d, d1 = len(risk), int((risk["g"] == 1).sum()), len(deaths), int((deaths["g"] == 1).sum())
        if n <= 1 or n1 in (0, n):
            continue
        o1 += d1
        e1 += n1 * d / n
        v += n1 * (n - n1) * d * (n - d) / (n ** 2 * (n - 1))
    chi = (o1 - e1) ** 2 / v if v else np.nan
    return {"chi2": float(chi), "p": float(stats.chi2.sf(chi, 1))}


def load_mimic():
    conn = connect("mimiciv31")
    df = pd.read_sql_query(
        """
        WITH liver AS (
            SELECT hadm_id,
                MAX(CASE
                    WHEN icd_version = 9 AND replace(icd_code,'.','') ~ '^(5712|5715|5716|5718|5722)' THEN 1
                    WHEN icd_version = 10 AND replace(icd_code,'.','') ~ '^(K70|K74|K717)' THEN 1
                    ELSE 0 END) AS chronic_liver,
                MAX(CASE
                    WHEN icd_version = 9 AND replace(icd_code,'.','') ~ '^(303|3050)' THEN 1
                    WHEN icd_version = 10 AND replace(icd_code,'.','') ~ '^F10' THEN 1
                    ELSE 0 END) AS aud,
                MAX(CASE
                    WHEN icd_version = 9 AND replace(icd_code,'.','') ~ '^070' THEN 1
                    WHEN icd_version = 10 AND replace(icd_code,'.','') ~ '^B1[5-9]' THEN 1
                    ELSE 0 END) AS viral_hep,
                MAX(CASE
                    WHEN icd_version = 10 AND replace(icd_code,'.','') ~ '^(K760|K758)' THEN 1
                    ELSE 0 END) AS fatty_liver
            FROM mimiciv_hosp.diagnoses_icd GROUP BY hadm_id
        )
        SELECT d.stay_id, d.hadm_id, d.admission_age AS age, d.gender,
               d.icu_intime, d.dod, i.first_careunit, sf.sofa,
               COALESCE(l.chronic_liver, 0) AS chronic_liver,
               COALESCE(l.aud, 0) AS aud,
               COALESCE(l.viral_hep, 0) AS viral_hep,
               COALESCE(l.fatty_liver, 0) AS fatty_liver
        FROM mimiciv_derived.icustay_detail d
        JOIN mimiciv_derived.sepsis3 s ON s.stay_id = d.stay_id AND s.sepsis3
        JOIN mimiciv_icu.icustays i ON i.stay_id = d.stay_id
        LEFT JOIN mimiciv_derived.first_day_sofa sf ON sf.stay_id = d.stay_id
        LEFT JOIN liver l ON l.hadm_id = d.hadm_id
        WHERE d.first_icu_stay AND d.admission_age >= 18
          AND COALESCE(l.chronic_liver, 0) = 0
        """,
        conn,
    )
    enzymes = pd.read_sql_query(
        """
        SELECT DISTINCT ON (b.stay_id) b.stay_id, e.charttime AS enzyme_time,
               e.ast, e.alt, e.bilirubin_total AS bili, e.ck_cpk AS ck
        FROM mimiciv_derived.icustay_detail b
        JOIN mimiciv_derived.enzyme e
          ON e.hadm_id = b.hadm_id
         AND e.charttime >= b.icu_intime
         AND e.charttime < b.icu_intime + interval '24 hour'
         AND e.ast IS NOT NULL AND e.alt IS NOT NULL
        WHERE b.first_icu_stay
        ORDER BY b.stay_id, e.charttime, e.specimen_id
        """,
        conn,
    )
    creat = pd.read_sql_query(
        """
        SELECT DISTINCT ON (stay_id) stay_id, creatinine
        FROM (
            SELECT b.stay_id, c.creatinine, c.charttime
            FROM mimiciv_derived.icustay_detail b
            JOIN mimiciv_derived.chemistry c
              ON c.hadm_id = b.hadm_id
             AND c.charttime >= b.icu_intime
             AND c.charttime < b.icu_intime + interval '24 hour'
             AND c.creatinine IS NOT NULL
            WHERE b.first_icu_stay AND b.admission_age >= 18
        ) z
        ORDER BY stay_id, charttime
        """,
        conn,
    )
    lact = pd.read_sql_query(
        """
        SELECT DISTINCT ON (stay_id) stay_id, lactate
        FROM (
            SELECT b.stay_id, bg.lactate, bg.charttime
            FROM mimiciv_derived.icustay_detail b
            JOIN mimiciv_derived.bg bg
              ON bg.hadm_id = b.hadm_id
             AND bg.charttime >= b.icu_intime
             AND bg.charttime < b.icu_intime + interval '24 hour'
             AND bg.lactate IS NOT NULL
            WHERE b.first_icu_stay AND b.admission_age >= 18
        ) z
        ORDER BY stay_id, charttime
        """,
        conn,
    )
    labs = creat.merge(lact, on="stay_id", how="outer")
    flow = pd.read_sql_query(
        """
        WITH liver AS (
            SELECT hadm_id,
                MAX(CASE
                    WHEN icd_version = 9 AND replace(icd_code,'.','') ~ '^(5712|5715|5716|5718|5722)' THEN 1
                    WHEN icd_version = 10 AND replace(icd_code,'.','') ~ '^(K70|K74|K717)' THEN 1
                    ELSE 0 END) AS chronic_liver
            FROM mimiciv_hosp.diagnoses_icd GROUP BY hadm_id
        )
        SELECT COUNT(*) AS n_sepsis,
               SUM(CASE WHEN COALESCE(l.chronic_liver, 0) = 1 THEN 1 ELSE 0 END) AS n_liver
        FROM mimiciv_derived.icustay_detail d
        JOIN mimiciv_derived.sepsis3 s ON s.stay_id = d.stay_id AND s.sepsis3
        LEFT JOIN liver l ON l.hadm_id = d.hadm_id
        WHERE d.first_icu_stay AND d.admission_age >= 18
        """,
        conn,
    )
    conn.close()
    df = df.merge(enzymes, on="stay_id", how="left").merge(labs, on="stay_id", how="left")
    df = prep_ratio(df)
    df["icu_intime"] = pd.to_datetime(df["icu_intime"])
    df["dod"] = pd.to_datetime(df["dod"], errors="coerce")
    df["enzyme_time"] = pd.to_datetime(df["enzyme_time"], errors="coerce")
    icu_days = (df["dod"].dt.normalize() - df["icu_intime"].dt.normalize()).dt.days
    df["event_icu"] = ((icu_days.notna()) & (icu_days >= 0) & (icu_days <= 28)).astype(int)
    enzyme_days = (df["dod"].dt.normalize() - df["enzyme_time"].dt.normalize()).dt.days
    before_enzyme = df["enzyme_time"].notna() & enzyme_days.notna() & (enzyme_days < 0)
    (OUT / "before_enzyme.json").write_text(
        json.dumps({"n_death_before_enzyme": int(before_enzyme.sum())}),
        encoding="utf-8",
    )
    df.loc[before_enzyme, "ratio"] = np.nan
    days = enzyme_days.where(df["enzyme_time"].notna(), icu_days)
    df["event_28"] = ((days.notna()) & (days >= 0) & (days <= 28)).astype(int)
    df["time_28"] = np.where(df["event_28"] == 1, np.maximum(days, 0.5), 28.0)
    df["male"] = (df["gender"] == "M").astype(int)
    unit = df["first_careunit"].fillna("")
    df["cvicu"] = unit.str.contains("CVICU", case=False)
    df["unit"] = np.select(
        [
            df["cvicu"],
            unit.str.contains("Medical Intensive Care Unit", case=False),
            unit.str.contains("Coronary", case=False),
            unit.str.contains("Surgical|Trauma", case=False),
        ],
        ["CVICU", "MICU", "CCU", "SICU"],
        default="Other",
    )
    _write_missingness(df)
    exposed = df[df["ratio"].notna()].copy()
    # Frozen MIMIC-IV 25th, 50th, and 75th percentiles. Do not recompute.
    cuts = [1.0357142857142858, 1.4444444444444444, 2.0294117647058822]
    exposed["quartile"] = assign_quartiles(exposed["ratio"], cuts)
    exposed["log_ratio_z"] = (exposed["log_ratio"] - exposed["log_ratio"].mean()) / exposed["log_ratio"].std(ddof=0)
    exposed["log_ast_z"] = (exposed["log_ast"] - exposed["log_ast"].mean()) / exposed["log_ast"].std(ddof=0)
    exposed["log_alt"] = np.log(exposed["alt"])
    exposed["log_alt_z"] = (exposed["log_alt"] - exposed["log_alt"].mean()) / exposed["log_alt"].std(ddof=0)
    (OUT / "mimic_flow.json").write_text(
        json.dumps(
            {
                "n_sepsis": int(flow.loc[0, "n_sepsis"]),
                "n_liver": int(flow.loc[0, "n_liver"]),
                "n_eligible": int(len(df)),
                "n_with_ratio": int(len(exposed)),
            }
        ),
        encoding="utf-8",
    )
    return exposed, cuts, int(len(df))


def load_eicu(cuts):
    conn = connect("eicu")
    df = pd.read_sql_query(
        """
        WITH sepsis AS (
            SELECT DISTINCT patientunitstayid FROM diagnosis
            WHERE lower(diagnosisstring) LIKE '%sepsis%' OR lower(diagnosisstring) LIKE '%septic%'
        ),
        liver AS (
            SELECT DISTINCT patientunitstayid FROM diagnosis
            WHERE lower(diagnosisstring) LIKE '%cirrhos%'
               OR lower(diagnosisstring) LIKE '%chronic liver%'
               OR icd9code LIKE '571.2%' OR icd9code LIKE '571.5%' OR icd9code LIKE '571.6%'
        ),
        labs AS (
            SELECT patientunitstayid, lower(labname) AS labname, labresult, labresultoffset
            FROM lab
            WHERE labresultoffset >= 0 AND labresultoffset < 1440 AND labresult IS NOT NULL
              AND lower(labname) IN ('ast (sgot)','alt (sgpt)','total bilirubin','lactate','creatinine')
        ),
        first_lab AS (
            SELECT DISTINCT ON (patientunitstayid, labname) patientunitstayid, labname, labresult
            FROM labs ORDER BY patientunitstayid, labname, labresultoffset, labresult
        ),
        wide AS (
            SELECT patientunitstayid,
                MAX(labresult) FILTER (WHERE labname='ast (sgot)') AS ast,
                MAX(labresult) FILTER (WHERE labname='alt (sgpt)') AS alt,
                MAX(labresult) FILTER (WHERE labname='total bilirubin') AS bili,
                MAX(labresult) FILTER (WHERE labname='lactate') AS lactate,
                MAX(labresult) FILTER (WHERE labname='creatinine') AS creatinine
            FROM first_lab GROUP BY patientunitstayid
        )
        SELECT p.patientunitstayid, p.age, p.gender, p.unittype, p.hospitaldischargestatus, a.apachescore,
               w.ast, w.alt, w.bili, w.lactate, w.creatinine
        FROM patient p
        JOIN sepsis s ON s.patientunitstayid = p.patientunitstayid
        LEFT JOIN wide w ON w.patientunitstayid = p.patientunitstayid
        LEFT JOIN liver ON liver.patientunitstayid = p.patientunitstayid
        LEFT JOIN (
            SELECT patientunitstayid, MAX(apachescore) AS apachescore
            FROM apachepatientresult GROUP BY patientunitstayid
        ) a ON a.patientunitstayid = p.patientunitstayid
        WHERE p.unitvisitnumber = 1
          AND liver.patientunitstayid IS NULL
        """,
        conn,
    )
    conn.close()
    raw = prep_ratio(df)
    raw["y"] = (raw["hospitaldischargestatus"] == "Expired").astype(int)
    raw["cardiac"] = raw["unittype"].isin(["CSICU", "Cardiac ICU", "CCU-CTICU", "CTICU"])
    miss_rows = []
    for label, part in (
        ("Usable AST and ALT", raw[raw["ratio"].notna()]),
        ("No usable pair", raw[raw["ratio"].isna()]),
    ):
        miss_rows.append(
            {
                "group": label,
                "n": int(len(part)),
                "hospital_death_n_pct": f"{int(part['y'].sum())} ({part['y'].mean() * 100:.1f})",
                "cardiac_n_pct": f"{int(part['cardiac'].sum())} ({part['cardiac'].mean() * 100:.1f})",
            }
        )
    pd.DataFrame(miss_rows).to_csv(OUT / "missingness_eicu.csv", index=False)
    (OUT / "eicu_flow.json").write_text(
        json.dumps({"eligible": int(len(raw)), "with_ratio": int(raw["ratio"].notna().sum())}),
        encoding="utf-8",
    )
    df = raw[raw["ratio"].notna()].copy()
    df["age"] = pd.to_numeric(df["age"].replace({"> 89": "90"}), errors="coerce")
    df["male"] = (df["gender"] == "Male").astype(int)
    df["y"] = (df["hospitaldischargestatus"] == "Expired").astype(int)
    df["quartile"] = assign_quartiles(df["ratio"], cuts)
    df["log_ast_z"] = (df["log_ast"] - df["log_ast"].mean()) / df["log_ast"].std(ddof=0)
    return df


def table_by_quartile(df, death_col, severity_col, severity_name):
    rows = []
    pieces = [("Overall", df)] + [(q, df[df["quartile"] == q]) for q in ["Q1", "Q2", "Q3", "Q4"]]

    def add(label, fn):
        rec = {"variable": label}
        for name, sub in pieces:
            rec[name] = fn(sub)
        rows.append(rec)

    add("N", lambda s: str(len(s)))
    add("Age, years", lambda s: iqr(s["age"]))
    add("Male", lambda s: npct(s["male"] == 1))
    add(severity_name, lambda s: iqr(s[severity_col]))
    add("Lactate, mmol/L", lambda s: iqr(s["lactate"]))
    add("Bilirubin, mg/dL", lambda s: iqr(s["bili"]))
    add("Creatinine, mg/dL", lambda s: iqr(s["creatinine"]))
    add("AST, U/L", lambda s: iqr(s["ast"]))
    add("ALT, U/L", lambda s: iqr(s["alt"]))
    add("AST/ALT", lambda s: iqr(s["ratio"]))
    add("Death", lambda s: npct(s[death_col] == 1))
    return pd.DataFrame(rows)


def mimic_models(df):
    d = df.copy()
    d["q2"] = (d["quartile"] == "Q2").astype(float)
    d["q3"] = (d["quartile"] == "Q3").astype(float)
    d["q4"] = (d["quartile"] == "Q4").astype(float)
    quart, _, _, _ = cox_fit(d, f"q2 + q3 + q4 + {COV} - 1", ["q2", "q3", "q4"])
    per, _, _, _ = cox_fit(df, f"log_ratio_z + {COV} - 1", ["log_ratio_z"])
    enzymes, _, _, _ = cox_fit(df, f"log_ast_z + log_alt_z + {COV} - 1", ["log_ast_z", "log_alt_z"])
    Xs = design(df, f"cr(log_ratio, df=4) + {COV} - 1")
    info = Xs.design_info
    same = df.loc[Xs.index].dropna(subset=["time_28", "event_28"])
    Xs = Xs.loc[same.index]
    Xl = design(same, f"log_ratio + {COV} - 1")
    same = same.loc[Xl.index]
    Xs = Xs.loc[same.index]
    res_s = fit_cox(same["time_28"], same["event_28"], Xs)
    res_l = fit_cox(same["time_28"], same["event_28"], Xl)
    df_diff = int(Xs.shape[1] - Xl.shape[1])
    chi = 2 * (float(res_s.llf) - float(res_l.llf))
    return {
        "quartiles": quart,
        "per_sd": per,
        "ast_and_alt": enzymes,
        "nonlinear": {"chi2": chi, "df": df_diff, "p": float(stats.chi2.sf(chi, df_diff)), "n": int(len(same))},
        "logrank_q4_vs_q1": logrank_q4(df),
    }, res_s, info, same


def rcs_plot(df, res, design_info, path):
    grid = np.linspace(df["log_ratio"].quantile(0.02), df["log_ratio"].quantile(0.98), 70)
    ref = float(np.log(df["ratio"].median()))
    base = df.iloc[0:1].copy()
    pred = pd.concat([base.assign(log_ratio=v, ratio=np.exp(v)) for v in list(grid) + [ref]], ignore_index=True)
    mat = build_design_matrices([design_info], pred, return_type="dataframe")[0]
    mat = mat.reindex(columns=design_info.column_names, fill_value=0.0)
    beta = np.asarray(res.params, float)
    covm = np.asarray(res.cov_params(), float)
    ref_x = np.asarray(mat.iloc[-1], float)
    hrs, lo, hi = [], [], []
    for i in range(len(grid)):
        delta = np.asarray(mat.iloc[i], float) - ref_x
        est = float(delta @ beta)
        se = float(np.sqrt(max(delta @ covm @ delta, 0)))
        hrs.append(np.exp(est))
        lo.append(np.exp(est - 1.96 * se))
        hi.append(np.exp(est + 1.96 * se))
    ratios = np.exp(grid)
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    ax.plot(ratios, hrs, color="#1f4e79", lw=2)
    ax.fill_between(ratios, lo, hi, color="#1f4e79", alpha=0.15)
    ax.axhline(1, color="#666666", lw=1, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("AST/ALT ratio")
    ax.set_ylabel("Hazard ratio for 28-day death")
    ax.set_title("MIMIC-IV adjusted spline")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def km_plot(df, path):
    colors = {"Q1": "#4c78a8", "Q2": "#f2c14e", "Q3": "#e07a3d", "Q4": "#b33a3a"}
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for q, color in colors.items():
        sub = df[df["quartile"] == q]
        t, p = km_curve(sub["time_28"], sub["event_28"])
        ax.step(t, p, where="post", color=color, lw=2, label=f"{q} (n={len(sub)})")
    ax.set_xlim(0, 28)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Days from ICU admission")
    ax.set_ylabel("Survival probability")
    ax.set_title("28-day survival by AST/ALT quartile")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def subgroup_table(df):
    df = df.copy()
    bands = {
        "All": df,
        "Male": df[df["male"] == 1],
        "Female": df[df["male"] == 0],
        "Age < 65": df[df["age"] < 65],
        "Age >= 65": df[df["age"] >= 65],
        "Lactate < 4": df[df["lactate"] < 4],
        "Lactate >= 4": df[df["lactate"] >= 4],
        "MICU": df[df["unit"] == "MICU"],
        "CVICU": df[df["unit"] == "CVICU"],
        "Non-CVICU": df[df["unit"] != "CVICU"],
        "Other ICU": df[(df["unit"] != "MICU") & (df["unit"] != "CVICU")],
    }
    rows = []
    for name, sub in bands.items():
        part = sub[sub["quartile"].isin(["Q1", "Q4"])].copy()
        part["q4"] = (part["quartile"] == "Q4").astype(int)
        cov = "age + sofa + lactate + log_bili + creatinine"
        if name not in ("Male", "Female"):
            cov += " + male"
        try:
            result, _, _, _ = cox_fit(part, f"q4 + {cov} - 1", ["q4"])
            row = {"subgroup": name, "n": result["n"], "events": result["events"]}
            row.update(result["terms"].get("q4", {}))
        except Exception as exc:
            row = {"subgroup": name, "error": str(exc)}
        rows.append(row)
    return pd.DataFrame(rows)


def eicu_logit(df):
    import statsmodels.formula.api as smf

    d = df.dropna(subset=["y", "quartile", "age", "male", "apachescore", "lactate", "log_bili", "creatinine"]).copy()
    d["q"] = pd.Categorical(d["quartile"], ["Q1", "Q2", "Q3", "Q4"])
    model = smf.logit(
        "y ~ C(q, Treatment(reference='Q1')) + age + male + apachescore + lactate + log_bili + creatinine",
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
    rates = d.groupby("q", observed=True)["y"].agg(n="size", death="mean").reset_index()
    return {"n": int(len(d)), "events": int(d["y"].sum()), "terms": terms, "rates": rates.to_dict(orient="records")}


def proportional_hazards(df):
    from lifelines import CoxPHFitter
    from lifelines.statistics import proportional_hazard_test

    need = ["time_28", "event_28", "age", "male", "sofa", "lactate", "log_bili", "creatinine", "unit", "quartile"]
    d = df.dropna(subset=need).copy()
    d["q2"] = (d["quartile"] == "Q2").astype(float)
    d["q3"] = (d["quartile"] == "Q3").astype(float)
    d["q4"] = (d["quartile"] == "Q4").astype(float)
    for unit in ["CCU", "SICU", "CVICU", "Other"]:
        d[f"u_{unit}"] = (d["unit"] == unit).astype(float)
    cols = [
        "time_28", "event_28", "q2", "q3", "q4", "age", "male", "sofa", "lactate",
        "log_bili", "creatinine", "u_CCU", "u_SICU", "u_CVICU", "u_Other",
    ]
    use = d[cols]
    model = CoxPHFitter()
    model.fit(use, "time_28", "event_28")
    tested = proportional_hazard_test(model, use, time_transform="rank")
    summary = tested.summary.reset_index()
    summary.to_csv(OUT / "ph_test.csv", index=False)
    return summary


def time_split_cox(df):
    """Separate Cox models for days 0-7 and days 7-28."""
    need = ["time_28", "event_28", "age", "male", "sofa", "lactate", "log_bili", "creatinine", "unit", "quartile"]
    d = df.dropna(subset=need).copy()
    d["q2"] = (d["quartile"] == "Q2").astype(float)
    d["q3"] = (d["quartile"] == "Q3").astype(float)
    d["q4"] = (d["quartile"] == "Q4").astype(float)
    early = d.copy()
    early["time"] = np.minimum(early["time_28"], 7.0)
    early["event"] = ((early["event_28"] == 1) & (early["time_28"] <= 7)).astype(int)
    late = d[d["time_28"] > 7].copy()
    late["time"] = late["time_28"] - 7.0
    late["event"] = late["event_28"]
    out = {}
    for label, part in (("day_0_7", early), ("day_7_28", late)):
        fitted, _, _, _ = cox_fit(part, f"q2 + q3 + q4 + {COV} - 1", ["q2", "q3", "q4"], time_col="time", event_col="event")
        out[label] = fitted
    return out


def death_logit(df):
    import statsmodels.formula.api as smf

    d = df.dropna(
        subset=["event_28", "quartile", "age", "male", "sofa", "lactate", "log_bili", "creatinine", "unit"]
    ).copy()
    d["q"] = pd.Categorical(d["quartile"], ["Q1", "Q2", "Q3", "Q4"])
    model = smf.logit(
        "event_28 ~ C(q, Treatment(reference='Q1')) + age + male + sofa + lactate + log_bili + creatinine + C(unit)",
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
    return {"n": int(len(d)), "events": int(d["event_28"].sum()), "terms": terms}


def main():
    mimic, cuts, n_before = load_mimic()
    summary, res_s, design_info, same = mimic_models(mimic)
    table_by_quartile(mimic, "event_28", "sofa", "SOFA").to_csv(OUT / "table1_mimic.csv", index=False)
    km_plot(mimic, OUT / "km_quartile.png")
    rcs_plot(same, res_s, design_info, OUT / "rcs_ratio.png")
    subgroup_table(mimic).to_csv(OUT / "subgroups_q4_vs_q1.csv", index=False)
    eicu = load_eicu(cuts)
    eicu_out = eicu_logit(eicu)
    eicu_show = eicu.copy()
    eicu_show["sofa"] = eicu_show["apachescore"]
    table_by_quartile(eicu_show, "y", "apachescore", "APACHE").to_csv(OUT / "table1_eicu.csv", index=False)
    ph = proportional_hazards(mimic)
    logit28 = death_logit(mimic)
    split = time_split_cox(mimic)
    counts = mimic["quartile"].value_counts().sort_index().to_dict()
    payload = {
        "n_eligible_before_enzymes": n_before,
        "n_with_ratio": int(len(mimic)),
        "quartile_counts": {str(k): int(v) for k, v in counts.items()},
        "quartile_cuts_ratio": cuts,
        "mimic": summary,
        "mimic_28d_logit": logit28,
        "time_split": split,
        "eicu": eicu_out,
        "eicu_n_with_ratio": int(len(eicu)),
        "eicu_quartile_counts": {str(k): int(v) for k, v in eicu["quartile"].value_counts().sort_index().items()},
        "ph_rank_p": {str(row["index"]): float(row["p"]) for _, row in ph.iterrows()},
    }
    (OUT / "results.json").write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    print(json.dumps(payload, indent=2, default=str)[:9000])


if __name__ == "__main__":
    main()
