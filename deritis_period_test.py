"""Test whether the quartile hazard differs before and after day 7.

Fits a period-by-quartile interaction on each of 30 imputed datasets
and pools the interaction with Rubin's rules. Also refits the eICU
logistic model on the frozen cutpoints.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import NelsonAalenFitter
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from statsmodels.duration.hazard_regression import PHReg

from deritis_paper_analysis import load_eicu, load_mimic
from deritis_round2 import add_covariates, eicu_with_unit, rubin

OUT = Path(r"F:\CKB数据库\analysis\deritis_paper")
USE = [
    "age", "male", "sofa", "lactate", "log_bili", "creatinine", "aud",
    "u_CCU", "u_SICU", "u_CVICU", "u_Other",
]


def period_fit(model_df):
    from lifelines import CoxTimeVaryingFitter

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


def main():
    mimic, cuts, _ = load_mimic()
    mimic = add_covariates(mimic)
    d = mimic.copy()
    d["q2"] = (d["quartile"] == "Q2").astype(float)
    d["q3"] = (d["quartile"] == "Q3").astype(float)
    d["q4"] = (d["quartile"] == "Q4").astype(float)
    for unit in ("CCU", "SICU", "CVICU", "Other"):
        d[f"u_{unit}"] = (d["unit"] == unit).astype(float)
    work = d.dropna(subset=["time_28", "event_28", "quartile"]).copy()
    naf = NelsonAalenFitter()
    naf.fit(work["time_28"], work["event_28"])
    work["na_hazard"] = np.asarray(
        naf.cumulative_hazard_at_times(work["time_28"]).to_numpy(), float
    ).ravel()
    imp_cols = USE + ["event_28", "na_hazard"]
    base = work[imp_cols]
    store = {term: {"coef": [], "var": []} for term in ("q4", "q4_late", "q4_after")}
    m = 30
    for seed in range(m):
        if seed % 5 == 0:
            print(f"impute {seed + 1}/{m}", flush=True)
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
        fitted = period_fit(model_df)
        for term, (coef, var) in fitted.items():
            store[term]["coef"].append(coef)
            store[term]["var"].append(var)
    pooled = {term: rubin(store[term]["coef"], store[term]["var"]) for term in store}
    product = pooled["q4"]["hr"] * pooled["q4_late"]["hr"]
    payload = {"period_interaction": pooled, "product_of_reported_hrs": product, "cuts": cuts}
    (OUT / "period_test.json").write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
