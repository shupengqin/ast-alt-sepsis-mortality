"""Re-estimate two reported sensitivities with quartiles inside the imputation model."""

from pathlib import Path
import json
import pickle

from deritis_closeout import varying_fit
from deritis_revision import USE, cox_terms, fill_once, prepare_imputation, rubin

CACHE = Path(__file__).resolve().parent / "deritis_paper" / "_revision_cache.pkl"
OUT = Path(__file__).resolve().parent / "deritis_paper" / "rerun_two.json"
EXTRA = ["charlson", "weight", "albumin", "vaso", "vent", "rrt"]


def pool(store):
    return {k: rubin(v["coef"], v["var"]) for k, v in store.items()}


def main():
    mimic, cuts, n_eligible, full = pickle.loads(open(CACHE, "rb").read())
    work = prepare_imputation(mimic, "time_28", "event_28")
    imp_core = ["q2", "q3", "q4"] + USE + ["y_imp", "na_hazard"]
    sens = {k: {"coef": [], "var": []} for k in ("q4", "q4_late", "q4_after")}
    joint = {k: {"coef": [], "var": []} for k in ("log_ast_z", "log_alt_z")}
    joint_formula = "log_ast_z + log_alt_z + " + " + ".join(USE) + " - 1"
    for seed in range(30):
        print("core", seed + 1, flush=True)
        model = fill_once(work, imp_core, seed)
        model["time_28"] = work["time_28"]
        model["event_28"] = work["event_28"]
        fitted = varying_fit(model)
        for term, (coef, var) in fitted.items():
            sens[term]["coef"].append(coef)
            sens[term]["var"].append(var)
        joint_fit, _, _, _ = cox_terms(model, joint_formula)
        for term in joint:
            joint[term]["coef"].append(joint_fit["terms"][term]["coef"])
            joint[term]["var"].append(joint_fit["terms"][term]["var"])
    imp_exp = imp_core + EXTRA
    expanded = {k: {"coef": [], "var": []} for k in ("q4",)}
    formula = "q2 + q3 + q4 + " + " + ".join(USE + EXTRA) + " - 1"
    from sklearn.experimental import enable_iterative_imputer
    from sklearn.impute import IterativeImputer
    import pandas as pd
    for seed in range(30):
        print("expanded", seed + 1, flush=True)
        imputer = IterativeImputer(random_state=seed, max_iter=12, sample_posterior=True)
        filled = pd.DataFrame(imputer.fit_transform(work[imp_exp]), columns=imp_exp, index=work.index)
        for col in ("male", "aud", "vaso", "vent", "rrt", "u_CCU", "u_SICU", "u_CVICU", "u_Other"):
            filled[col] = (filled[col] >= 0.5).astype(float)
        model = filled[USE + EXTRA].copy()
        model["q2"] = work["q2"]
        model["q3"] = work["q3"]
        model["q4"] = work["q4"]
        model["time_28"] = work["time_28"]
        model["event_28"] = work["event_28"]
        fitted, _, _, _ = cox_terms(model, formula)
        expanded["q4"]["coef"].append(fitted["terms"]["q4"]["coef"])
        expanded["q4"]["var"].append(fitted["terms"]["q4"]["var"])
    payload = {
        "covariate_time": pool(sens),
        "joint": pool(joint),
        "expanded_q4": pool(expanded)["q4"],
    }
    open(OUT, "w", encoding="utf-8").write(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
