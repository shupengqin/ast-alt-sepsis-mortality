"""Primary hazard ratios with quartiles inside the imputation model.

Also refits the day-7 time split in the 24-hour landmark cohort.
"""

from pathlib import Path
import json
import pickle

import numpy as np

from deritis_revision import (
    FORMULA,
    USE,
    cox_terms,
    fill_once,
    landmark_frame,
    period_fit,
    prepare_imputation,
    rubin,
)

CACHE = Path(__file__).resolve().parent / "deritis_paper" / "_revision_cache.pkl"
OUT = Path(__file__).resolve().parent / "deritis_paper" / "mi_quartile.json"


def pool(store):
    return {k: rubin(v["coef"], v["var"]) for k, v in store.items()}


def main():
    mimic, cuts, n_eligible, full = pickle.loads(open(CACHE, "rb").read())
    work = prepare_imputation(mimic, "time_28", "event_28")
    imp_cols = ["q2", "q3", "q4"] + USE + ["y_imp", "na_hazard"]
    primary = {k: {"coef": [], "var": []} for k in ("q2", "q3", "q4")}
    for seed in range(30):
        print("primary", seed + 1, flush=True)
        model = fill_once(work, imp_cols, seed)
        model["time_28"] = work["time_28"]
        model["event_28"] = work["event_28"]
        fitted, _, _, _ = cox_terms(model, FORMULA)
        for term in primary:
            primary[term]["coef"].append(fitted["terms"][term]["coef"])
            primary[term]["var"].append(fitted["terms"][term]["var"])
    lm, n_excl = landmark_frame(mimic)
    lm_work = prepare_imputation(lm, "time_28", "event_28")
    period = {k: {"coef": [], "var": []} for k in ("q4_early", "q4_late_ratio", "q4_after")}
    for seed in range(30):
        print("landmark", seed + 1, flush=True)
        model = fill_once(lm_work, imp_cols, seed)
        model["time_28"] = lm_work["time_28"]
        model["event_28"] = lm_work["event_28"]
        fitted = period_fit(model, 7)
        for term in period:
            period[term]["coef"].append(fitted[term][0])
            period[term]["var"].append(fitted[term][1])
    payload = {
        "primary_n": int(len(work)),
        "primary_events": int(work["event_28"].sum()),
        "primary": pool(primary),
        "landmark_excluded": n_excl,
        "landmark_n": int(len(lm_work)),
        "landmark_events": int(lm_work["event_28"].sum()),
        "landmark_period": pool(period),
        "cuts": cuts,
        "n_eligible": n_eligible,
    }
    open(OUT, "w", encoding="utf-8").write(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
