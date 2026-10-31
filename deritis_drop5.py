"""Drop stays whose enzyme chart time is after hospital death, then refit the primary models."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from deritis_revision import USE, attach_hospital, cox_terms, fill_once, period_fit, prepare_imputation
from deritis_round2 import rubin

CACHE = Path(Path(__file__).resolve().parent / "deritis_paper" / "_revision_cache.pkl")
OUT = Path(Path(__file__).resolve().parent / "deritis_paper" / "drop5.json")


def pool(store):
    return {k: rubin(v["coef"], v["var"]) for k, v in store.items()}


def fit_store(work, formula, terms, time_col="time_28", event_col="event_28", label="fit"):
    imp_cols = ["q2", "q3", "q4"] + USE + ["y_imp", "na_hazard"]
    store = {k: {"coef": [], "var": []} for k in terms}
    for seed in range(30):
        print(label, seed + 1, flush=True)
        model = fill_once(work, imp_cols, seed)
        model[time_col] = work[time_col]
        model[event_col] = work[event_col]
        for col in ("r_z", "m_z", "log_ast_z", "log_alt_z"):
            if col in work.columns:
                model[col] = work[col]
        if terms[0].startswith("q4_"):
            fitted = period_fit(model, 7, time_col=time_col, event_col=event_col)
            for term in terms:
                store[term]["coef"].append(fitted[term][0])
                store[term]["var"].append(fitted[term][1])
        else:
            fitted, _, _, _ = cox_terms(model, formula, time_col=time_col, event_col=event_col)
            for term in terms:
                store[term]["coef"].append(fitted["terms"][term]["coef"])
                store[term]["var"].append(fitted["terms"][term]["var"])
    return pool(store)


def main():
    mimic, cuts, n_eligible, full = pd.read_pickle(CACHE)
    if "deathtime" not in mimic.columns:
        mimic = attach_hospital(mimic)
    enzyme = pd.to_datetime(mimic["enzyme_time"])
    death = pd.to_datetime(mimic["deathtime"], errors="coerce")
    dod = pd.to_datetime(mimic["dod"], errors="coerce")
    bad = death.notna() & enzyme.notna() & (enzyme > death)
    print("dropped", int(bad.sum()), mimic.loc[bad, "quartile"].value_counts().to_dict(), flush=True)
    keep = mimic.loc[~bad].copy()
    same = dod.notna() & enzyme.notna() & (dod.dt.normalize() == enzyme.dt.normalize())
    same = same.loc[keep.index]
    counts = {}
    for q in ("Q1", "Q2", "Q3", "Q4"):
        sub = keep[keep["quartile"] == q]
        early = (sub["event_28"] == 1) & (sub["time_28"] <= 7)
        alive7 = sub["time_28"] > 7
        late = alive7 & (sub["event_28"] == 1)
        counts[q] = {
            "n": int(len(sub)),
            "deaths": int(sub["event_28"].sum()),
            "week1": int(early.sum()),
            "alive7": int(alive7.sum()),
            "late_deaths": int(late.sum()),
        }
    hours = (enzyme.loc[keep.index] - pd.to_datetime(keep["icu_intime"])).dt.total_seconds() / 3600.0
    desc = {
        "n": int(len(keep)),
        "deaths": int(keep["event_28"].sum()),
        "within_6h": int((hours <= 6).sum()),
        "quartiles": counts,
        "same_day_remaining": int(same.sum()),
        "same_day_by_q": keep.loc[same, "quartile"].value_counts().to_dict(),
    }
    print(json.dumps(desc), flush=True)

    work = prepare_imputation(keep, "time_28", "event_28")
    primary_f = "q2 + q3 + q4 + " + " + ".join(USE) + " - 1"
    primary = fit_store(work, primary_f, ["q2", "q3", "q4"], label="primary")
    period = fit_store(
        work, primary_f, ["q4_early", "q4_late_ratio", "q4_after"], label="period"
    )
    level = keep.copy()
    level["r_z"] = level["log_ast"] - level["log_alt"]
    level["m_z"] = (level["log_ast"] + level["log_alt"]) / 2.0
    level["r_z"] = (level["r_z"] - level["r_z"].mean()) / level["r_z"].std(ddof=0)
    level["m_z"] = (level["m_z"] - level["m_z"].mean()) / level["m_z"].std(ddof=0)
    level["log_ast_z"] = (level["log_ast"] - level["log_ast"].mean()) / level["log_ast"].std(ddof=0)
    level["log_alt_z"] = (level["log_alt"] - level["log_alt"].mean()) / level["log_alt"].std(ddof=0)
    work_l = prepare_imputation(level, "time_28", "event_28")
    enzymes = fit_store(
        work_l,
        "log_ast_z + log_alt_z + " + " + ".join(USE) + " - 1",
        ["log_ast_z", "log_alt_z"],
        label="joint",
    )
    ratio_level = fit_store(
        work_l,
        "r_z + m_z + " + " + ".join(USE) + " - 1",
        ["r_z", "m_z"],
        label="level",
    )
    reduced = keep.loc[~same].copy()
    work_s = prepare_imputation(reduced, "time_28", "event_28")
    same_day = fit_store(
        work_s, primary_f, ["q4_early", "q4_late_ratio", "q4_after"], label="sameday"
    )
    payload = {
        "desc": desc,
        "primary": primary,
        "period": period,
        "joint": enzymes,
        "ratio_level": ratio_level,
        "same_day_n": int(len(reduced)),
        "same_day_period": same_day,
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2), flush=True)


if __name__ == "__main__":
    main()
