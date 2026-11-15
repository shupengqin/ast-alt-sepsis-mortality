"""Complete-case time curves for Figure 4A.

The interval at a given day uses the delta method. For a contrast x of the
coefficients, the variance is x' V x, which includes the covariance between
the level and the time terms. It is not the variance of either coefficient alone.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxTimeVaryingFitter

from deritis_revision import USE

ROOT = Path(__file__).resolve().parent
mimic = pd.read_pickle(ROOT / "deritis_paper" / "_revision_cache.pkl")[0]
enzyme = pd.to_datetime(mimic["enzyme_time"])
death = pd.to_datetime(mimic["deathtime"], errors="coerce")
bad = death.notna() & enzyme.notna() & (enzyme > death)
cc = mimic.loc[~bad].dropna(subset=["time_28", "event_28", "q2", "q3", "q4"] + USE).reset_index(drop=True)
times = np.sort(cc.loc[cc["event_28"] == 1, "time_28"].unique())
grid = np.unique(np.concatenate(([0.0], times)))
pieces = []
for start, stop in zip(grid[:-1], grid[1:]):
    sub = cc.loc[cc["time_28"] > start, ["time_28", "event_28", "q2", "q3", "q4"] + USE].copy()
    sub["start"] = start
    sub["stop"] = np.minimum(sub["time_28"].to_numpy(), stop)
    sub = sub.loc[sub["stop"] > start]
    sub["event_i"] = ((sub["time_28"] <= stop) & (sub["event_28"] == 1)).astype(int)
    sub["id"] = sub.index
    logt = np.log(stop)
    sub["q2_logt"] = sub["q2"] * logt
    sub["q3_logt"] = sub["q3"] * logt
    sub["q4_logt"] = sub["q4"] * logt
    sub["q4_t"] = sub["q4"] * stop
    sub["q4_hinge"] = sub["q4"] * max(stop - 7.0, 0.0)
    pieces.append(sub)
long = pd.concat(pieces, ignore_index=True)
common = ["id", "start", "stop", "event_i", "q2", "q3", "q4"] + USE
log_fit = CoxTimeVaryingFitter().fit(
    long[common + ["q2_logt", "q3_logt", "q4_logt"]],
    id_col="id", event_col="event_i", start_col="start", stop_col="stop", show_progress=False,
)
lin_fit = CoxTimeVaryingFitter().fit(
    long[common + ["q4_t", "q4_hinge"]],
    id_col="id", event_col="event_i", start_col="start", stop_col="stop", show_progress=False,
)


def curve(fit, day, kind):
    names = list(fit.params_.index)
    var = np.asarray(fit.variance_matrix_, float)
    x = np.zeros(len(names))
    x[names.index("q4")] = 1.0
    if kind == "log":
        x[names.index("q4_logt")] = np.log(day)
    else:
        x[names.index("q4_t")] = day
        x[names.index("q4_hinge")] = max(day - 7.0, 0.0)
    est = float(x @ np.asarray(fit.params_, float))
    se = float(np.sqrt(max(float(x @ var @ x), 0.0)))
    return float(np.exp(est)), float(np.exp(est - 1.96 * se)), float(np.exp(est + 1.96 * se))


# Reported Table 5 values. The check fails if the band ignores covariance.
hr, lo, hi = curve(lin_fit, 28, "piece")
assert abs(hr - 1.283) < 0.01 and abs(lo - 0.974) < 0.01 and abs(hi - 1.692) < 0.01, (hr, lo, hi)
rows = []
for day in np.linspace(1, 28, 55):
    for kind, label in (("log", "Log time"), ("piece", "Piecewise linear")):
        fit = log_fit if kind == "log" else lin_fit
        hr, lo, hi = curve(fit, day, kind)
        rows.append((label, day, hr, lo, hi))
dest = ROOT / "deritis_paper" / "figures" / "source_time_curves.csv"
pd.DataFrame(rows, columns=["model", "day", "hr", "lo", "hi"]).to_csv(dest, index=False)
print("covariance check passed", round(hr, 3), round(lo, 3), round(hi, 3))
print("wrote", dest)
