"""Hospital-clustered standard errors for the imputed eICU hospital-death model."""

from pathlib import Path
import json

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer

from deritis_paper_analysis import connect, load_eicu
from deritis_round2 import rubin

CUTS = [1.0357142857142858, 1.4444444444444444, 2.0294117647058822]


def main():
    eicu = load_eicu(CUTS)
    conn = connect("eicu")
    alcohol = pd.read_sql_query(
        """
        SELECT DISTINCT patientunitstayid FROM diagnosis
        WHERE lower(diagnosisstring) LIKE '%alcohol%'
           OR icd9code LIKE '303%' OR icd9code LIKE '305.0%' OR icd9code LIKE '291%'
        """,
        conn,
    )
    hosp = pd.read_sql_query("SELECT patientunitstayid, hospitalid FROM patient", conn)
    conn.close()
    eicu["aud"] = eicu["patientunitstayid"].isin(set(alcohol["patientunitstayid"])).astype(float)
    hosp_map = hosp.drop_duplicates("patientunitstayid").set_index("patientunitstayid")["hospitalid"]
    eicu["hospitalid"] = eicu["patientunitstayid"].map(hosp_map)
    d = eicu.dropna(subset=["y", "quartile", "hospitalid"]).copy()
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
    d["q2"] = (d["quartile"] == "Q2").astype(float)
    d["q3"] = (d["quartile"] == "Q3").astype(float)
    d["q4"] = (d["quartile"] == "Q4").astype(float)
    imp_cols = ["q2", "q3", "q4", "age", "male", "apachescore", "lactate", "log_bili", "creatinine", "y", "aud", "u_MICU", "u_CCU", "u_SICU", "u_Cardiac"]
    plain = {k: {"coef": [], "var": []} for k in ("q2", "q3", "q4")}
    clus = {k: {"coef": [], "var": []} for k in ("q2", "q3", "q4")}
    for seed in range(30):
        print("seed", seed + 1, flush=True)
        imputer = IterativeImputer(random_state=seed, max_iter=12, sample_posterior=True)
        filled = pd.DataFrame(imputer.fit_transform(d[imp_cols]), columns=imp_cols, index=d.index)
        for col in ("male", "y", "aud", "u_MICU", "u_CCU", "u_SICU", "u_Cardiac"):
            filled[col] = (filled[col] >= 0.5).astype(float)
        X = filled[["age", "male", "apachescore", "lactate", "log_bili", "creatinine", "aud", "u_MICU", "u_CCU", "u_SICU", "u_Cardiac"]].copy()
        X["q2"] = d["q2"]
        X["q3"] = d["q3"]
        X["q4"] = d["q4"]
        X = sm.add_constant(X, has_constant="add")
        X = X[[c for c in X.columns if c == "const" or float(X[c].std()) > 1e-8]]
        y = d["y"]
        res = sm.Logit(y, X).fit(disp=False, maxiter=200)
        res_c = sm.Logit(y, X).fit(disp=False, maxiter=200, cov_type="cluster", cov_kwds={"groups": d["hospitalid"]})
        for term in ("q2", "q3", "q4"):
            plain[term]["coef"].append(float(res.params[term]))
            plain[term]["var"].append(float(res.cov_params().loc[term, term]))
            clus[term]["coef"].append(float(res_c.params[term]))
            clus[term]["var"].append(float(res_c.cov_params().loc[term, term]))
        if seed == 0:
            print("seed0", {t: round(float(np.exp(res.params[t])), 3) for t in ("q2", "q3", "q4")})
    payload = {
        "n": int(len(d)),
        "events": int(d["y"].sum()),
        "hospitals": int(d["hospitalid"].nunique()),
        "independent": {k: rubin(v["coef"], v["var"]) for k, v in plain.items()},
        "cluster": {k: rubin(v["coef"], v["var"]) for k, v in clus.items()},
    }
    path = Path(__file__).resolve().parent / "deritis_paper" / "eicu_cluster_q.json"
    open(path, "w", encoding="utf-8").write(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
