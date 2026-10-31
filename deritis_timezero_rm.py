"""Time-zero audit and the ratio-versus-level reparameterization.

Does not recompute quartile cutpoints and does not refit the 24-hour exclusion.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from deritis_paper_analysis import cox_fit
from deritis_revision import USE, attach_hospital, cox_terms, fill_once, prepare_imputation
from deritis_round2 import rubin

CACHE = Path(Path(__file__).resolve().parent / "deritis_paper" / "_revision_cache.pkl")
OUT = Path(Path(__file__).resolve().parent / "deritis_paper" / "timezero_rm.json")


def _summ(hours):
    s = pd.Series(hours).dropna()
    qs = s.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).to_dict()
    return {
        "n": int(s.shape[0]),
        "median": float(s.median()),
        "q1": float(qs[0.25]),
        "q3": float(qs[0.75]),
        "p10": float(qs[0.1]),
        "p90": float(qs[0.9]),
        "min": float(s.min()),
        "max": float(s.max()),
        "within_6h": int((s <= 6).sum()),
        "within_12h": int((s <= 12).sum()),
    }


def pool_term(store, term):
    return rubin(store[term]["coef"], store[term]["var"])


def main():
    mimic, cuts, n_eligible, full = pd.read_pickle(CACHE)
    if "deathtime" not in mimic.columns:
        mimic = attach_hospital(mimic)
    enzyme = pd.to_datetime(mimic["enzyme_time"])
    icu = pd.to_datetime(mimic["icu_intime"])
    dod = pd.to_datetime(mimic["dod"], errors="coerce")
    death = pd.to_datetime(mimic["deathtime"], errors="coerce")
    hours = (enzyme - icu).dt.total_seconds() / 3600.0
    same_day = dod.notna() & enzyme.notna() & (dod.dt.normalize() == enzyme.dt.normalize())
    deaths_28 = mimic["event_28"] == 1
    timed = death.notna()
    lab_after_death = timed & enzyme.notna() & (enzyme > death)
    lab_after_same_day = same_day & lab_after_death
    same_day_timed = same_day & timed
    same_day_date_only = same_day & ~timed
    by_q = {}
    for q in ("Q1", "Q2", "Q3", "Q4"):
        m = mimic["quartile"] == q
        by_q[q] = {
            "n": int(m.sum()),
            "same_day_deaths": int((m & same_day).sum()),
            "deaths_28": int((m & deaths_28).sum()),
        }
    audit = {
        "n": int(len(mimic)),
        "hours_icu_to_enzyme": _summ(hours),
        "same_calendar_day_deaths": int(same_day.sum()),
        "same_day_by_quartile": by_q,
        "deaths_28": int(deaths_28.sum()),
        "deaths_28_with_deathtime": int((deaths_28 & timed).sum()),
        "deaths_28_date_only": int((deaths_28 & ~timed).sum()),
        "lab_after_deathtime": int(lab_after_death.sum()),
        "same_day_with_deathtime": int(same_day_timed.sum()),
        "same_day_lab_after_deathtime": int(lab_after_same_day.sum()),
        "same_day_date_only": int(same_day_date_only.sum()),
        "hours_negative": int((hours < 0).sum()),
    }
    print(json.dumps(audit, indent=2), flush=True)

    keep = mimic.loc[~same_day].copy()
    keep["time_7"] = np.minimum(keep["time_28"], 7.0)
    keep["event_7"] = ((keep["event_28"] == 1) & (keep["time_28"] <= 7)).astype(int)
    work = prepare_imputation(keep, "time_7", "event_7")
    imp_cols = ["q2", "q3", "q4"] + USE + ["y_imp", "na_hazard"]
    early_formula = "q2 + q3 + q4 + " + " + ".join(USE) + " - 1"
    early_store = {k: {"coef": [], "var": []} for k in ("q2", "q3", "q4")}
    for seed in range(30):
        print("early", seed + 1, flush=True)
        model = fill_once(work, imp_cols, seed)
        model["time_7"] = work["time_7"]
        model["event_7"] = work["event_7"]
        fitted, _, _, _ = cox_terms(model, early_formula, time_col="time_7", event_col="event_7")
        for term in early_store:
            early_store[term]["coef"].append(fitted["terms"][term]["coef"])
            early_store[term]["var"].append(fitted["terms"][term]["var"])
    early = {k: pool_term(early_store, k) for k in early_store}
    early_meta = {
        "n": int(len(keep)),
        "excluded_same_day_deaths": int(same_day.sum()),
        "events_day_1_to_7": int(keep["event_7"].sum()),
        "q4": early["q4"],
    }
    print("EARLY", json.dumps(early_meta, indent=2), flush=True)

    level = mimic.copy()
    level["r"] = level["log_ast"] - level["log_alt"]
    level["m"] = (level["log_ast"] + level["log_alt"]) / 2.0
    level["r_z"] = (level["r"] - level["r"].mean()) / level["r"].std(ddof=0)
    level["m_z"] = (level["m"] - level["m"].mean()) / level["m"].std(ddof=0)
    work_l = prepare_imputation(level, "time_28", "event_28")
    imp_l = ["q2", "q3", "q4"] + USE + ["y_imp", "na_hazard"]
    formula = "r_z + m_z + " + " + ".join(USE) + " - 1"
    store = {k: {"coef": [], "var": []} for k in ("r_z", "m_z")}
    for seed in range(30):
        print("level", seed + 1, flush=True)
        model = fill_once(work_l, imp_l, seed)
        model["r_z"] = work_l["r_z"]
        model["m_z"] = work_l["m_z"]
        model["time_28"] = work_l["time_28"]
        model["event_28"] = work_l["event_28"]
        fitted, _, _, _ = cox_terms(model, formula)
        for term in store:
            store[term]["coef"].append(fitted["terms"][term]["coef"])
            store[term]["var"].append(fitted["terms"][term]["var"])
    level_mi = {k: pool_term(store, k) for k in store}

    cc_cols = ["r_z", "m_z", "time_28", "event_28"] + USE
    cc = level.dropna(subset=cc_cols)
    base_f = "r_z + " + " + ".join(USE) + " - 1"
    full_f = "r_z + m_z + " + " + ".join(USE) + " - 1"
    base, _, _, _ = cox_fit(cc, base_f, ["r_z"])
    both, _, _, _ = cox_fit(cc, full_f, ["r_z", "m_z"])
    from scipy import stats

    lrt = float(stats.chi2.sf(2 * (both["llf"] - base["llf"]), 1))
    payload = {
        "audit": audit,
        "exclude_same_day": early_meta,
        "ratio_and_level_imputed": level_mi,
        "complete_case_lrt_p": lrt,
        "complete_case_n": base["n"],
        "complete_case_m": both["terms"]["m_z"],
        "complete_case_r": both["terms"]["r_z"],
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2), flush=True)


if __name__ == "__main__":
    main()
