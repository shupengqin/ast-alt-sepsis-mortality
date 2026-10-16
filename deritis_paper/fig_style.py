"""Shared style for the AST/ALT multi-panel figures."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(r"F:\CKB数据库\analysis\deritis_paper\figures")
OUT.mkdir(parents=True, exist_ok=True)

Q_COLORS = {
    "Q1": "#3B6EA5",
    "Q2": "#4C9A8A",
    "Q3": "#E07A3D",
    "Q4": "#8E5A8A",
}
MIMIC_C = "#3B6EA5"
EICU_C = "#E07A3D"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.6,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.frameon": False,
    "legend.fontsize": 6,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})


def save_pub(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.tiff", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def panel_letter(ax, letter):
    ax.text(
        -0.08, 1.08, letter, transform=ax.transAxes,
        fontsize=10, fontweight="bold", va="bottom", ha="left", color="black",
    )


def forest(ax, labels, estimates, lows, highs, colors, title, xlabel):
    from matplotlib.ticker import FixedLocator, FuncFormatter

    y = np.arange(len(labels))[::-1]
    ax.axvline(1, color="#888888", lw=0.6, ls="--")
    for i, (est, lo, hi, col) in enumerate(zip(estimates, lows, highs, colors)):
        ax.plot([lo, hi], [y[i], y[i]], color=col, lw=1.2, solid_capstyle="butt")
        ax.plot(est, y[i], "o", color=col, ms=4.2, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator([1, 2, 3]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v:g}"))
    ax.xaxis.set_minor_locator(FixedLocator([]))
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=8, pad=10)
    ax.set_xlim(0.65, 3.2)
