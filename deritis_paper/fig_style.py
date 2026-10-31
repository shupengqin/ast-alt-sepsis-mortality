"""Shared style for the AST/ALT multi-panel figures."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(parents=True, exist_ok=True)

Q_COLORS = {
    "Q1": "#3B6EA5",
    "Q2": "#4C9A8A",
    "Q3": "#E07A3D",
    "Q4": "#8E5A8A",
}
MIMIC_C = "#3B6EA5"
EICU_C = "#1B4F72"

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
    common = dict(bbox_inches="tight", facecolor="white", transparent=False)
    fig.savefig(OUT / f"{stem}.pdf", **common)
    fig.savefig(OUT / f"{stem}.svg", **common)
    fig.savefig(OUT / f"{stem}.tiff", dpi=600, **common)
    png = OUT / f"{stem}.png"
    fig.savefig(png, dpi=300, **common)
    plt.close(fig)
    try:
        from PIL import Image
        for path, dpi in ((png, 300), (OUT / f"{stem}.tiff", 600)):
            image = Image.open(path)
            if image.mode == "RGB":
                continue
            canvas = Image.new("RGB", image.size, (255, 255, 255))
            if image.mode == "RGBA":
                canvas.paste(image, mask=image.split()[-1])
            else:
                canvas.paste(image.convert("RGB"))
            canvas.save(path, dpi=(dpi, dpi))
    except Exception:
        pass


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
        if abs(hi - lo) < 1e-8:
            ax.plot(est, y[i], "o", markerfacecolor="white", markeredgecolor=col, ms=4.2, zorder=3)
        else:
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
