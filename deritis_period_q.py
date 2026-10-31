"""Day-7 time split with quartiles inside the imputation model."""

from pathlib import Path
import json

from deritis_mi_q import CACHE, fill_once, period_fit, prepare_imputation, rubin
import pickle

from deritis_revision import USE

OUT = Path(__file__).resolve().parent / "deritis_paper" / "period_quartile.json"


def main():
    mimic, cuts, n_eligible, full = pickle.loads(open(CACHE, "rb").read())
    work = prepare_imputation(mimic, "time_28", "event_28")
    imp_cols = ["q2", "q3", "q4"] + USE + ["y_imp", "na_hazard"]
    store = {k: {"coef": [], "var": []} for k in ("q4_early", "q4_late_ratio", "q4_after")}
    for seed in range(30):
        print("period", seed + 1, flush=True)
        model = fill_once(work, imp_cols, seed)
        model["time_28"] = work["time_28"]
        model["event_28"] = work["event_28"]
        fitted = period_fit(model, 7)
        for term in store:
            store[term]["coef"].append(fitted[term][0])
            store[term]["var"].append(fitted[term][1])
    payload = {k: rubin(v["coef"], v["var"]) for k, v in store.items()}
    open(OUT, "w", encoding="utf-8").write(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
