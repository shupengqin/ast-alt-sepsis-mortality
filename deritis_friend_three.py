"""24-hour landmark with a restarted clock, joint-enzyme C-index, and IPW balance."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from deritis_revision import (
    USE,
    fill_once,
    period_fit,
    prepare_imputation,
    rubin,
    cox_terms,
)
from lifelines.utils import concordance_index

CACHE = Path(Path(__file__).resolve().parent / "deritis_paper" / "_revision_cache.pkl")
OUT = Path(Path(__file__).resolve().parent / "deritis_paper" / "friend_three.json")


def pool(store):
    return {k: rubin(v["coef"], v["var"]) for k, v in store.items()}


def smd(x1, x0, w1=None):
    x1 = np.asarray(x1, float)
    x0 = np.asarray(x0, float)
    if w1 is None:
        m1, v1 = np.nanmean(x1), np.nanvar(x1, ddof=1)
    else:
        w = np.asarray(w1, float)
        w = w / w.sum()
        m1 = np.nansum(w * x1)
        v1 = np.nansum(w * (x1 - m1) ** 2) / (1 - np.nansum(w ** 2))
    m0, v0 = np.nanmean(x0), np.nanvar(x0, ddof=1)
    denom = np.sqrt((v1 + v0) / 2)
    return float((m1 - m0) / denom) if denom else None


def main():
    mimic, cuts, n_eligible, full = pd.read_pickle(CACHE)
    enzyme = pd.to_datetime(mimic["enzyme_time"])
    death_clock = pd.to_datetime(mimic["deathtime"], errors="coerce")
    chart_after = death_clock.notna() & enzyme.notna() & (enzyme > death_clock)
    cohort = mimic.loc[~chart_after].copy()
    icu = pd.to_datetime(cohort["icu_intime"])
    death = death_clock.loc[cohort.index].where(death_clock.loc[cohort.index].notna(), pd.to_datetime(cohort["dod"], errors="coerce"))
    landmark = icu + pd.Timedelta(hours=24)
    alive = death.isna() | (death >= landmark)
    lm = cohort.loc[alive].copy()
    death_k = death.loc[lm.index]
    follow = (death_k - landmark.loc[lm.index]).dt.total_seconds() / 86400.0
    horizon = 27.0
    event = death_k.notna() & follow.notna() & (follow >= 0) & (follow <= horizon)
    lm["time_lm"] = np.where(event, np.maximum(follow, 0.5), horizon)
    lm["event_lm"] = event.astype(int)
    early_dead = int(((lm["event_lm"] == 1) & (lm["time_lm"] <= 6)).sum())
    late_dead = int(((lm["event_lm"] == 1) & (lm["time_lm"] > 6)).sum())
    print("landmark n", len(lm), "excluded", int((~alive).sum()), "early deaths", early_dead, "late deaths", late_dead, flush=True)

    work = prepare_imputation(lm, "time_lm", "event_lm")
    imp_cols = ["q2", "q3", "q4"] + USE + ["y_imp", "na_hazard"]
    store = {k: {"coef": [], "var": []} for k in ("q4_early", "q4_late_ratio", "q4_after")}
    for seed in range(30):
        print("landmark", seed + 1, flush=True)
        model = fill_once(work, imp_cols, seed)
        model["time_lm"] = work["time_lm"]
        model["event_lm"] = work["event_lm"]
        fitted = period_fit(model, 6.0, time_col="time_lm", event_col="event_lm")
        for term in store:
            store[term]["coef"].append(fitted[term][0])
            store[term]["var"].append(fitted[term][1])
    landmark_est = pool(store)

    base_f = " + ".join(USE) + " - 1"
    ratio_f = "log_ratio_z + " + base_f
    both_f = "log_ast_z + log_alt_z + " + base_f
    work_c = prepare_imputation(cohort, "time_28", "event_28")
    for col in ("log_ratio_z", "log_ast_z", "log_alt_z"):
        work_c[col] = cohort[col]
    rows = []
    imp_c = ["q2", "q3", "q4"] + USE + ["y_imp", "na_hazard"]
    for seed in range(30):
        print("cindex", seed + 1, flush=True)
        model = fill_once(work_c, imp_c, seed)
        model["time_28"] = work_c["time_28"]
        model["event_28"] = work_c["event_28"]
        for col in ("log_ratio_z", "log_ast_z", "log_alt_z"):
            model[col] = work_c[col]
        fits = {}
        for name, formula in (("base", base_f), ("ratio", ratio_f), ("both", both_f)):
            fitted, res, X, use = cox_terms(model, formula)
            fits[name] = (res, X, use)
        common = fits["base"][2].index.intersection(fits["ratio"][2].index).intersection(fits["both"][2].index)
        part = model.loc[common].copy()
        part["time_7"] = np.minimum(part["time_28"], 7.0)
        part["event_7"] = ((part["event_28"] == 1) & (part["time_28"] <= 7)).astype(int)
        rec = {}
        for horizon, tcol, ecol in (("day_28", "time_28", "event_28"), ("day_7", "time_7", "event_7")):
            rec[horizon] = {}
            for name, (res, X, use) in fits.items():
                xx = X.loc[common]
                lp = np.asarray(xx, float) @ np.asarray(res.params, float)
                rec[horizon][name] = float(concordance_index(part[tcol], -lp, part[ecol]))
        rows.append(rec)
    cindex = {}
    for horizon in ("day_28", "day_7"):
        base = float(np.mean([r[horizon]["base"] for r in rows]))
        ratio = float(np.mean([r[horizon]["ratio"] for r in rows]))
        both = float(np.mean([r[horizon]["both"] for r in rows]))
        cindex[horizon] = {"base": base, "ratio": ratio, "both": both, "delta_ratio": ratio - base, "delta_both": both - base}

    d = full.copy()
    d["measured"] = d["ratio"].notna().astype(int)
    d["sofa_miss"] = d["sofa"].isna().astype(float)
    d["sofa_fill"] = d["sofa"].fillna(d["sofa"].median())
    d["aud"] = d["aud"].fillna(0).astype(float)
    d["male"] = d["male"].fillna(0).astype(float)
    for unit in ("CCU", "SICU", "CVICU", "Other"):
        d[f"u_{unit}"] = (d["unit"] == unit).astype(float)
    xcols = ["age", "male", "sofa_fill", "sofa_miss", "aud", "u_CCU", "u_SICU", "u_CVICU", "u_Other"]
    X = d[xcols].astype(float)
    X = X[[c for c in X.columns if float(X[c].std()) > 1e-8]]
    clf = LogisticRegression(max_iter=1000)
    clf.fit(X, d["measured"])
    p = np.clip(clf.predict_proba(X)[:, 1], 0.02, 0.98)
    marg = float(d["measured"].mean())
    d["w_raw"] = np.where(d["measured"] == 1, marg / p, np.nan)
    keep_ids = set(cohort["stay_id"])
    measured = d.loc[(d["measured"] == 1) & d["stay_id"].isin(keep_ids)].copy()
    lo, hi = measured["w_raw"].quantile([0.01, 0.99])
    measured["w"] = measured["w_raw"].clip(lo, hi)
    w = measured["w"]
    ess = float((w.sum() ** 2) / np.sum(np.square(w)))
    qs = w.quantile([0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1]).to_dict()
    unmeasured = d.loc[d["measured"] == 0]
    balance = []
    for label, col in (
        ("Age", "age"),
        ("Male", "male"),
        ("SOFA", "sofa_fill"),
        ("SOFA missing", "sofa_miss"),
        ("Alcohol use disorder", "aud"),
        ("CCU", "u_CCU"),
        ("SICU", "u_SICU"),
        ("CVICU", "u_CVICU"),
        ("Other unit", "u_Other"),
    ):
        balance.append({
            "variable": label,
            "smd_unweighted": smd(measured[col], unmeasured[col]),
            "smd_weighted": smd(measured[col], unmeasured[col], measured["w"]),
        })
    payload = {
        "landmark": {
            "n": int(len(lm)),
            "excluded_before_24h": int((~alive).sum()),
            "early_deaths": early_dead,
            "late_deaths": late_dead,
            "estimates": landmark_est,
        },
        "cindex": cindex,
        "ipw": {
            "n": int(len(measured)),
            "ess": ess,
            "clip": [float(lo), float(hi)],
            "quantiles": {str(k): float(v) for k, v in qs.items()},
            "balance": balance,
        },
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
