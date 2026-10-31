"""Day-3 and day-14 splits with quartiles inside the imputation model."""

from pathlib import Path
import json
import pickle

from deritis_revision import USE, period_fit, prepare_imputation, rubin

CACHE = Path(__file__).resolve().parent / "deritis_paper" / "_revision_cache.pkl"
OUT = Path(__file__).resolve().parent / "deritis_paper" / "period_cuts_q.json"


def main():
    mimic, cuts, n_eligible, full = pickle.loads(open(CACHE, "rb").read())
    work = prepare_imputation(mimic, "time_28", "event_28")
    from deritis_revision import fill_once
    imp_cols = ["q2", "q3", "q4"] + USE + ["y_imp", "na_hazard"]
    keys = []
    for cut in (3, 14):
        keys += [f"c{cut}_early", f"c{cut}_ratio", f"c{cut}_after"]
    store = {k: {"coef": [], "var": []} for k in keys}
    for seed in range(30):
        print("cut", seed + 1, flush=True)
        model = fill_once(work, imp_cols, seed)
        model["time_28"] = work["time_28"]
        model["event_28"] = work["event_28"]
        for cut in (3, 14):
            fitted = period_fit(model, cut)
            store[f"c{cut}_early"]["coef"].append(fitted["q4_early"][0])
            store[f"c{cut}_early"]["var"].append(fitted["q4_early"][1])
            store[f"c{cut}_ratio"]["coef"].append(fitted["q4_late_ratio"][0])
            store[f"c{cut}_ratio"]["var"].append(fitted["q4_late_ratio"][1])
            store[f"c{cut}_after"]["coef"].append(fitted["q4_after"][0])
            store[f"c{cut}_after"]["var"].append(fitted["q4_after"][1])
    payload = {k: rubin(v["coef"], v["var"]) for k, v in store.items()}
    open(OUT, "w", encoding="utf-8").write(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
