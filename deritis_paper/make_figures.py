"""Five multi-panel figures for the AST/ALT sepsis manuscript."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psycopg2
from matplotlib.patches import FancyBboxPatch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.size": 7,
    "axes.linewidth": 0.6,
})

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from deritis_paper.fig_style import (
    EICU_C,
    MIMIC_C,
    OUT,
    Q_COLORS,
    forest,
    panel_letter,
    save_pub,
)
from deritis_paper_analysis import COV, design, fit_cox, load_eicu, load_mimic

QUARTILES = ["Q1", "Q2", "Q3", "Q4"]


def mimic_flow_counts():
    conn = psycopg2.connect(host="localhost", port=5442, user="postgres", password=os.environ["PGPASSWORD"], dbname="mimiciv31")
    conn.autocommit = True
    with conn.cursor() as cur:
        cur.execute("SET statement_timeout = 0")
    df = pd.read_sql_query(
        """
        WITH liver AS (
            SELECT hadm_id,
                MAX(CASE
                    WHEN icd_version = 9 AND replace(icd_code,'.','') ~ '^(5712|5715|5716|5718|5722)' THEN 1
                    WHEN icd_version = 10 AND replace(icd_code,'.','') ~ '^(K70|K74|K717)' THEN 1
                    ELSE 0 END) AS chronic_liver
            FROM mimiciv_hosp.diagnoses_icd GROUP BY hadm_id
        )
        SELECT COALESCE(l.chronic_liver, 0) AS chronic_liver,
               CASE WHEN i.first_careunit ILIKE '%CVICU%' THEN 1 ELSE 0 END AS cvicu
        FROM mimiciv_derived.icustay_detail d
        JOIN mimiciv_derived.sepsis3 s ON s.stay_id = d.stay_id AND s.sepsis3
        JOIN mimiciv_icu.icustays i ON i.stay_id = d.stay_id
        LEFT JOIN liver l ON l.hadm_id = d.hadm_id
        WHERE d.first_icu_stay AND d.admission_age >= 18
        """,
        conn,
    )
    conn.close()
    n0 = int(len(df))
    n_liver = int(df["chronic_liver"].sum())
    rest = df[df["chronic_liver"] == 0]
    return {
        "sepsis_first": n0,
        "chronic_liver": n_liver,
        "cvicu": int(rest["cvicu"].sum()),
        "eligible": int((rest["cvicu"] == 0).sum()),
    }


def eicu_flow_counts():
    conn = psycopg2.connect(host="localhost", port=5442, user="postgres", password=os.environ["PGPASSWORD"], dbname="eicu")
    conn.autocommit = True
    with conn.cursor() as cur:
        cur.execute("SET statement_timeout = 0")
    df = pd.read_sql_query(
        """
        WITH sepsis AS (
            SELECT DISTINCT patientunitstayid FROM diagnosis
            WHERE lower(diagnosisstring) LIKE '%sepsis%' OR lower(diagnosisstring) LIKE '%septic%'
        ),
        liver AS (
            SELECT DISTINCT patientunitstayid FROM diagnosis
            WHERE lower(diagnosisstring) LIKE '%cirrhos%'
               OR lower(diagnosisstring) LIKE '%chronic liver%'
               OR icd9code LIKE '571.2%' OR icd9code LIKE '571.5%' OR icd9code LIKE '571.6%'
        )
        SELECT CASE WHEN liver.patientunitstayid IS NULL THEN 0 ELSE 1 END AS chronic_liver,
               CASE WHEN COALESCE(p.unittype,'') IN ('CSICU','Cardiac ICU','CCU-CTICU','CTICU')
                    THEN 1 ELSE 0 END AS cardiac
        FROM patient p
        JOIN sepsis s ON s.patientunitstayid = p.patientunitstayid
        LEFT JOIN liver ON liver.patientunitstayid = p.patientunitstayid
        WHERE p.unitvisitnumber = 1
        """,
        conn,
    )
    conn.close()
    n0 = int(len(df))
    n_liver = int(df["chronic_liver"].sum())
    rest = df[df["chronic_liver"] == 0]
    return {
        "sepsis_first": n0,
        "chronic_liver": n_liver,
        "cardiac": int(rest["cardiac"].sum()),
        "eligible": int((rest["cardiac"] == 0).sum()),
    }


def qcounts(df):
    return {q: int((df["quartile"] == q).sum()) for q in QUARTILES}


def mortality_table(df, ycol):
    rows = []
    for q in QUARTILES:
        sub = df[df["quartile"] == q]
        rows.append({"quartile": q, "n": int(len(sub)), "deaths": int(sub[ycol].sum()),
                     "rate": float(sub[ycol].mean()) if len(sub) else np.nan})
    return pd.DataFrame(rows)


def km_xy(sub):
    t = sub["time_28"].to_numpy(float)
    e = sub["event_28"].to_numpy(int)
    order = np.argsort(t)
    t, e = t[order], e[order]
    at = len(t)
    surv = 1.0
    xs, ys = [0.0], [1.0]
    i = 0
    while i < len(t):
        tt = t[i]
        d = c = 0
        while i < len(t) and t[i] == tt:
            d += int(e[i] == 1)
            c += int(e[i] == 0)
            i += 1
        if d and at:
            surv *= 1 - d / at
            xs.append(float(tt))
            ys.append(float(surv))
        at -= d + c
    return xs, ys


def design_matrix(df, cols):
    return pd.get_dummies(df[cols], drop_first=True).astype(float)


def fit_prob(df, ycol, cols):
    d = df.dropna(subset=[ycol] + cols).copy()
    y = d[ycol].astype(int).to_numpy()
    X = design_matrix(d, cols)
    clf = LogisticRegression(max_iter=500, solver="lbfgs")
    clf.fit(X, y)
    return y, clf.predict_proba(X)[:, 1]


def nested_metrics(df, ycol, clinical, extra, n_boot=80):
    d = df.dropna(subset=[ycol] + clinical + extra).copy()
    y = d[ycol].astype(int).to_numpy()
    specs = {
        "Clinical": clinical,
        "+ AST/ALT": clinical + ["log_ratio_z"],
        "+ AST level": clinical + ["log_ratio_z", "log_ast_z"],
    }
    point = {}
    for name, cols in specs.items():
        yy, pp = fit_prob(d, ycol, cols)
        point[name] = {
            "auroc": float(roc_auc_score(yy, pp)),
            "auprc": float(average_precision_score(yy, pp)),
        }
    rng = np.random.default_rng(20261015)
    keys = ["ratio_vs_clinical", "both_vs_clinical", "both_vs_ratio"]
    store = {k: {"auroc": [], "auprc": []} for k in keys}
    n = len(d)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        part = d.iloc[idx]
        if part[ycol].nunique() < 2:
            continue
        scores = {}
        ok = True
        for name, cols in specs.items():
            try:
                yy, pp = fit_prob(part, ycol, cols)
                scores[name] = (roc_auc_score(yy, pp), average_precision_score(yy, pp))
            except Exception:
                ok = False
                break
        if not ok:
            continue
        pairs = {
            "ratio_vs_clinical": ("+ AST/ALT", "Clinical"),
            "both_vs_clinical": ("+ AST level", "Clinical"),
            "both_vs_ratio": ("+ AST level", "+ AST/ALT"),
        }
        for key, (a, b) in pairs.items():
            store[key]["auroc"].append(scores[a][0] - scores[b][0])
            store[key]["auprc"].append(scores[a][1] - scores[b][1])
    summary = {"point": point, "n": int(n), "events": int(y.sum())}
    for key, metrics in store.items():
        summary[key] = {}
        for metric, vals in metrics.items():
            arr = np.asarray(vals, float)
            summary[key][metric] = {
                "mean": float(np.mean(arr)),
                "low": float(np.quantile(arr, 0.025)),
                "high": float(np.quantile(arr, 0.975)),
            }
    return summary


def calibration_curve(df, ycol, cols, n_bins=8):
    yy, pp = fit_prob(df, ycol, cols)
    tab = pd.DataFrame({"y": yy, "p": pp})
    tab["bin"] = pd.qcut(tab["p"], n_bins, duplicates="drop")
    return tab.groupby("bin", observed=True).agg(pred=("p", "mean"), obs=("y", "mean"))


def slice_rates(df, ycol, mask):
    sub = df.loc[mask.fillna(False)]
    return [float(sub.loc[sub["quartile"] == q, ycol].mean()) * 100 if (sub["quartile"] == q).any() else np.nan
            for q in QUARTILES]


def draw_column(ax, title, header, fill, steps, groups, group_title):
    ax.set_xlim(0, 10)
    ax.set_ylim(2.3, 12.2)
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.25, 11.2), 9.5, 0.75, boxstyle="square,pad=0",
                                facecolor="#243447", edgecolor="#243447", linewidth=0))
    ax.text(5.0, 11.57, title, ha="center", va="center", color="white", fontsize=8, fontweight="bold")
    y = 10.2
    for i, (main, excl) in enumerate(steps):
        ax.add_patch(FancyBboxPatch((0.3, y - 1.2), 5.15, 1.2, boxstyle="square,pad=0",
                                    facecolor="#F4F6F8", edgecolor="#243447", linewidth=0.6))
        ax.text(2.87, y - 0.6, main, ha="center", va="center", fontsize=6.2, linespacing=1.15)
        if excl:
            ax.annotate("", xy=(5.85, y - 0.6), xytext=(5.45, y - 0.6),
                        arrowprops=dict(arrowstyle="-|>", color="black", lw=0.6))
            ax.add_patch(FancyBboxPatch((5.95, y - 1.05), 3.75, 0.95, boxstyle="square,pad=0",
                                        facecolor="white", edgecolor="black", linewidth=0.6))
            ax.text(7.82, y - 0.57, excl, ha="center", va="center", fontsize=5.8, linespacing=1.1)
        if i < len(steps) - 1:
            ax.annotate("", xy=(2.87, y - 1.42), xytext=(2.87, y - 1.2),
                        arrowprops=dict(arrowstyle="-|>", color="black", lw=0.6))
        y -= 1.75
    ax.add_patch(FancyBboxPatch((0.3, 2.55), 9.4, 2.35, boxstyle="square,pad=0",
                                facecolor="white", edgecolor="black", linewidth=0.7))
    ax.text(5.0, 4.55, group_title, ha="center", va="center", fontsize=6.4, fontweight="bold")
    ax.text(5.0, 3.85, "Cutpoints fixed in MIMIC-IV: 1.0357, 1.4444, 2.0294",
            ha="center", va="center", fontsize=6.2)
    ax.text(5.0, 3.15,
            "    ".join(f"{lab} n={groups[lab]:,}" for lab in QUARTILES),
            ha="center", va="center", fontsize=6.3)


def figure1(m, e, mq, eq):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.6))
    m_base = m["sepsis_first"] - m["chronic_liver"]
    e_base = e["sepsis_first"] - e["chronic_liver"]
    m_steps = [
        (f"First ICU stay with sepsis\nn = {m['sepsis_first']:,}",
         f"Excluded: n = {m['chronic_liver']:,}\nChronic liver disease"),
        (f"No chronic liver disease\nAll care units kept\nn = {m_base:,}",
         f"Excluded: n = {m_base - sum(mq.values()):,}\nNo usable AST/ALT in 24 h"),
        (f"Analysis cohort\nn = {sum(mq.values()):,}", None),
    ]
    e_steps = [
        (f"First ICU stay with sepsis\nn = {e['sepsis_first']:,}",
         f"Excluded: n = {e['chronic_liver']:,}\nChronic liver disease"),
        (f"No chronic liver disease\nCardiac units kept\nn = {e_base:,}",
         f"Excluded: n = {e_base - sum(eq.values()):,}\nNo usable AST/ALT in 24 h"),
        (f"Analysis cohort\nn = {sum(eq.values()):,}", None),
    ]
    draw_column(axes[0], "MIMIC-IV", "#243447", "#F4F6F8",
                m_steps, mq, "Quartile sample")
    draw_column(axes[1], "eICU-CRD", "#243447", "#F4F6F8",
                e_steps, eq, "Same cutpoints")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.96, bottom=0.04, wspace=0.08)
    save_pub(fig, "Figure_1")


def figure2(mimic, m_mort, e_mort):
    fig = plt.figure(figsize=(7.2, 6.0))
    gs = fig.add_gridspec(2, 2, hspace=0.48, wspace=0.38)
    ax = fig.add_subplot(gs[0, 0])
    panel_letter(ax, "A")
    for q in QUARTILES:
        xs, ys = km_xy(mimic[mimic["quartile"] == q])
        ax.step(xs, ys, where="post", color=Q_COLORS[q], lw=1.5,
                label=f"{q} (n={(mimic['quartile']==q).sum():,})")
    ax.set_xlim(0, 28)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Days from enzyme measurement")
    ax.set_ylabel("Survival probability")
    ax.set_title("MIMIC-IV, all with a ratio", fontsize=8, pad=10)
    ax.legend(loc="lower left", fontsize=5.4)

    ax = fig.add_subplot(gs[0, 1])
    panel_letter(ax, "B")
    x = np.arange(4)
    w = 0.36
    mr = m_mort["rate"].to_numpy() * 100
    er = e_mort["rate"].to_numpy() * 100
    ax.bar(x - w / 2, mr, w, color=MIMIC_C, label="MIMIC-IV, 28-day death")
    ax.bar(x + w / 2, er, w, color=EICU_C, label="eICU, hospital death")
    ax.set_xticks(x)
    ax.set_xticklabels(QUARTILES)
    ax.set_ylabel("Death (%)")
    ax.set_ylim(0, 52)
    ax.set_title("Crude death, two endpoints", fontsize=8, pad=10)
    ax.legend(loc="upper left", fontsize=5.2)

    ax = fig.add_subplot(gs[1, 0])
    panel_letter(ax, "C")
    forest(
        ax,
        ["MIMIC Q2, HR", "eICU Q2, OR", "MIMIC Q3, HR", "eICU Q3, OR", "MIMIC Q4, HR", "eICU Q4, OR"],
        [1.05, 0.87, 1.19, 1.05, 1.45, 1.64],
        [0.96, 0.76, 1.09, 0.91, 1.33, 1.43],
        [1.15, 1.01, 1.30, 1.21, 1.59, 1.88],
        [MIMIC_C, EICU_C] * 3,
        "Adjusted associations",
        "HR (MIMIC-IV) or OR (eICU)",
    )

    ax = fig.add_subplot(gs[1, 1])
    panel_letter(ax, "D")
    rcs_panel(ax, mimic)
    fig.savefig(OUT / "Figure_2.pdf", bbox_inches="tight")
    fig.savefig(OUT / "Figure_2.svg", bbox_inches="tight")
    fig.savefig(OUT / "Figure_2.tiff", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "Figure_2.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def labeled_bars(ax, tab, ylabel, title):
    rates = tab["rate"].to_numpy() * 100
    ns = tab["n"].to_numpy()
    ax.bar(np.arange(4), rates, color=[Q_COLORS[q] for q in QUARTILES])
    for i, (r, n) in enumerate(zip(rates, ns)):
        ax.text(i, r + 1.2, f"{r:.1f}%\nn={n:,}", ha="center", va="bottom", fontsize=5.4)
    ax.set_xticks(np.arange(4))
    ax.set_xticklabels(QUARTILES)
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, 58)
    ax.set_title(title, fontsize=8, pad=10)


def rcs_panel(ax, mimic):
    from patsy import build_design_matrices
    formula = f"cr(log_ratio, df=4) + {COV} - 1"
    X = design(mimic, formula)
    info = X.design_info
    same = mimic.loc[X.index].dropna(subset=["time_28", "event_28"])
    X = X.loc[same.index]
    res = fit_cox(same["time_28"], same["event_28"], X)
    grid = np.linspace(same["log_ratio"].quantile(0.02), same["log_ratio"].quantile(0.98), 60)
    ref = float(np.log(same["ratio"].median()))
    base = same.iloc[0:1]
    pred = pd.concat([base.assign(log_ratio=v, ratio=float(np.exp(v))) for v in list(grid) + [ref]], ignore_index=True)
    mat = build_design_matrices([info], pred, return_type="dataframe")[0]
    mat = mat.reindex(columns=list(info.column_names), fill_value=0.0)
    beta = np.asarray(res.params, float)
    covm = np.asarray(res.cov_params(), float)
    ref_x = np.asarray(mat.iloc[-1], float)
    hrs, lo, hi = [], [], []
    for i in range(len(grid)):
        delta = np.asarray(mat.iloc[i], float) - ref_x
        est = float(delta @ beta)
        se = float(np.sqrt(max(float(delta @ covm @ delta), 0)))
        hrs.append(np.exp(est))
        lo.append(np.exp(est - 1.96 * se))
        hi.append(np.exp(est + 1.96 * se))
    ratios = np.exp(grid)
    ax.plot(ratios, hrs, color=MIMIC_C, lw=1.6)
    ax.fill_between(ratios, lo, hi, color=MIMIC_C, alpha=0.15, linewidth=0)
    ax.axhline(1, color="#888888", lw=0.6, ls="--")
    ax.set_xscale("log")
    from matplotlib.ticker import FixedLocator, FuncFormatter
    ax.xaxis.set_major_locator(FixedLocator([0.5, 1, 2, 4]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v:g}"))
    ax.xaxis.set_minor_locator(FixedLocator([]))
    ax.set_xlabel("AST/ALT ratio")
    ax.set_ylabel("Hazard ratio")
    ax.set_title("Adjusted spline", fontsize=8, pad=10)


def incr_forest(ax, mimic_m, eicu_m, metric, title):
    labels = ["+ AST/ALT vs clinical", "+ AST level vs clinical", "AST level vs ratio"]
    keys = ["ratio_vs_clinical", "both_vs_clinical", "both_vs_ratio"]
    y = np.arange(len(labels))[::-1]
    ax.axvline(0, color="#888888", lw=0.6, ls="--")
    for i, key in enumerate(keys):
        for block, col, shift in ((mimic_m, MIMIC_C, 0.12), (eicu_m, EICU_C, -0.12)):
            rec = block[key][metric]
            yy = y[i] + shift
            ax.plot([rec["low"], rec["high"]], [yy, yy], color=col, lw=1.15)
            ax.plot(rec["mean"], yy, "o", color=col, ms=4)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel(f"Difference in apparent {metric.upper()}")
    ax.set_title(title, fontsize=8, pad=10)


def figure3(m_met, e_met):
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4))
    names = ["Clinical", "+ AST/ALT", "+ AST level"]
    ax = axes[0, 0]
    panel_letter(ax, "A")
    ax.plot(np.arange(3), [m_met["point"][k]["auroc"] for k in names], "o-", color=MIMIC_C, label="MIMIC-IV")
    ax.plot(np.arange(3), [e_met["point"][k]["auroc"] for k in names], "s-", color=EICU_C, label="eICU-CRD")
    ax.set_xticks(np.arange(3))
    ax.set_xticklabels(names, rotation=12, ha="right")
    ax.set_ylabel("Apparent AUROC")
    ax.set_ylim(0.65, 0.90)
    ax.set_title("Nested models", fontsize=8, pad=10)
    ax.legend()

    ax = axes[0, 1]
    panel_letter(ax, "B")
    ax.plot(np.arange(3), [m_met["point"][k]["auprc"] for k in names], "o-", color=MIMIC_C, label="MIMIC-IV")
    ax.plot(np.arange(3), [e_met["point"][k]["auprc"] for k in names], "s-", color=EICU_C, label="eICU-CRD")
    ax.set_xticks(np.arange(3))
    ax.set_xticklabels(names, rotation=12, ha="right")
    ax.set_ylabel("Apparent AUPRC")
    ax.set_ylim(0.40, 0.70)
    ax.set_title("Precision-recall", fontsize=8, pad=10)
    ax.legend()

    ax = axes[1, 0]
    panel_letter(ax, "C")
    incr_forest(ax, m_met, e_met, "auroc", "Incremental AUROC")
    ax.plot([], [], "o", color=MIMIC_C, label="MIMIC-IV")
    ax.plot([], [], "o", color=EICU_C, label="eICU-CRD")
    ax.legend(loc="lower right")

    ax = axes[1, 1]
    panel_letter(ax, "D")
    incr_forest(ax, m_met, e_met, "auprc", "Incremental AUPRC")
    fig.tight_layout()
    save_pub(fig, "Figure_3")


def figure4(curves):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3))
    for ax, title, letter, base, full, col in [
        (axes[0], "MIMIC-IV", "A", curves["m_base"], curves["m_full"], MIMIC_C),
        (axes[1], "eICU-CRD", "B", curves["e_base"], curves["e_full"], EICU_C),
    ]:
        panel_letter(ax, letter)
        ax.plot([0, 1], [0, 1], ls="--", color="#888888", lw=0.7)
        ax.plot(base["pred"], base["obs"], "o-", color="#666666", ms=3.4, label="Clinical")
        ax.plot(full["pred"], full["obs"], "o-", color=col, ms=3.4, label="Clinical + AST/ALT")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel("Mean predicted mortality")
        ax.set_ylabel("Observed mortality")
        ax.set_title(title, fontsize=8, pad=10)
        ax.legend(loc="upper left")
    fig.tight_layout()
    save_pub(fig, "Figure_4")


def slice_rates(df, ycol, mask):
    sub = df.loc[mask.fillna(False)]
    out = []
    for q in QUARTILES:
        part = sub[sub["quartile"] == q]
        out.append(float(part[ycol].mean()) * 100 if len(part) else np.nan)
    return out


def figure5(mimic, subgroup):
    fig = plt.figure(figsize=(7.2, 3.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15], wspace=0.55)
    ax = fig.add_subplot(gs[0, 0])
    panel_letter(ax, "A")
    x = np.arange(4)
    w = 0.36
    low = slice_rates(mimic, "event_28", mimic["lactate"] < 4)
    high = slice_rates(mimic, "event_28", mimic["lactate"] >= 4)
    ax.bar(x - w / 2, low, w, color="#8AA0B4", label="Lactate < 4")
    ax.bar(x + w / 2, high, w, color="#8E5A8A", label="Lactate >= 4")
    ax.set_xticks(x)
    ax.set_xticklabels(QUARTILES)
    ax.set_ylim(0, 75)
    ax.set_ylabel("28-day mortality (%)")
    ax.set_title("MIMIC-IV, by lactate", fontsize=8, pad=10)
    ax.legend()

    ax = fig.add_subplot(gs[0, 1])
    panel_letter(ax, "B")
    rows = [
        ("Men vs women", 1.03, 0.84, 1.27),
        ("Age >= 65", 0.80, 0.64, 1.00),
        ("Lactate >= 4", 1.09, 0.85, 1.41),
        ("Cardiac surgical ICU", 0.61, 0.36, 1.03),
    ][::-1]
    y = np.arange(len(rows))
    ax.axvline(1, color="#888888", lw=0.6, ls="--")
    for yi, (label, hr, lo, hi) in enumerate(rows):
        ax.plot([lo, hi], [yi, yi], color="#243447", lw=1.05)
        ax.plot(hr, yi, "o", color="#243447", ms=3.6)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlim(0.25, 1.70)
    ax.set_xlabel("Interaction HR")
    ax.set_title("Effect modification", fontsize=8, pad=10)
    save_pub(fig, "Figure_3")


def main():
    print("flow", flush=True)
    m_flow = mimic_flow_counts()
    e_flow = eicu_flow_counts()
    print(m_flow, e_flow, flush=True)
    mimic, cuts, _ = load_mimic()
    eicu = load_eicu(cuts)
    eicu["log_ratio_z"] = (eicu["log_ratio"] - eicu["log_ratio"].mean()) / eicu["log_ratio"].std(ddof=0)
    mq, eq = qcounts(mimic), qcounts(eicu)
    m_mort = mortality_table(mimic, "event_28")
    e_mort = mortality_table(eicu, "y")
    m_mort.to_csv(OUT / "source_mimic_quartile_mortality.csv", index=False)
    e_mort.to_csv(OUT / "source_eicu_quartile_mortality.csv", index=False)
    figure1(m_flow, e_flow, mq, eq)
    figure2(mimic, m_mort, e_mort)
    sub = pd.read_csv(r"F:\CKB数据库\analysis\deritis_paper\subgroups_q4_vs_q1.csv")
    figure5(mimic, sub)
    print("done", OUT)


if __name__ == "__main__":
    main()
