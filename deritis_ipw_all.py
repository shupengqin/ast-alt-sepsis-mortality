"""IPW on every measured stay, without requiring lactate."""

from pathlib import Path
import pickle

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy import stats
from sklearn.linear_model import LogisticRegression

from deritis_revision import add_quartiles

OUT = Path(__file__).resolve().parent / "deritis_paper"
mimic, cuts, n_eligible, full = pickle.loads(open(OUT / "_revision_cache.pkl", "rb").read())
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
clf = LogisticRegression(max_iter=1000, C=1.0)
clf.fit(X, d["measured"])
p = np.clip(clf.predict_proba(X)[:, 1], 0.02, 0.98)
marg = float(d["measured"].mean())
p_hat = pd.Series(p, index=d.index)
src = full.loc[full["ratio"].notna(), ["stay_id"]].copy()
src["w"] = (marg / p_hat).loc[src.index].to_numpy()
measured = mimic.merge(src[["stay_id", "w"]], on="stay_id", how="left")
lo, hi = measured["w"].quantile([0.01, 0.99])
measured["w"] = measured["w"].clip(lo, hi)
cols = ["time_28", "event_28", "q2", "q3", "q4", "age", "male", "sofa", "aud", "u_CCU", "u_SICU", "u_CVICU", "u_Other", "w"]
cc = measured.dropna(subset=cols).copy()
fit_w = CoxPHFitter().fit(cc[cols], "time_28", "event_28", weights_col="w", robust=True, show_progress=False)
fit_u = CoxPHFitter().fit(cc[[c for c in cols if c != "w"]], "time_28", "event_28", show_progress=False)

def grab(model):
    out = {}
    for term in ("q2", "q3", "q4"):
        coef = float(model.params_[term])
        se = float(model.standard_errors_[term])
        out[term] = (
            round(float(np.exp(coef)), 3),
            round(float(np.exp(coef - 1.96 * se)), 3),
            round(float(np.exp(coef + 1.96 * se)), 3),
        )
    return out

ess = float((cc["w"].sum() ** 2) / np.sum(np.square(cc["w"])))
print("n", len(cc), "events", int(cc.event_28.sum()), "ess", round(ess, 1))
print("w", round(float(cc.w.median()), 3), round(float(lo), 3), round(float(hi), 3))
print("p meas", round(float(p[d.measured == 1].mean()), 3), "p unmeas", round(float(p[d.measured == 0].mean()), 3))
print("weighted", grab(fit_w))
print("unweighted", grab(fit_u))
