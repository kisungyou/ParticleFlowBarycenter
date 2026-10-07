"""Shared PGF export using the manuscript's pdfLaTeX/Computer Modern fonts.
Figures are constructed at their final display width, with no later font scaling.
All plot labels use 12-point type; tick and legend text use 10-point type.
"""
from pathlib import Path
import os
import tempfile
if "MPLCONFIGDIR" not in os.environ:
    _mpl_cache = tempfile.TemporaryDirectory(prefix="otto-figures-")
    os.environ["MPLCONFIGDIR"] = _mpl_cache.name
import matplotlib as mpl
mpl.use("pgf")
mpl.rcParams.update({
    "pgf.texsystem": "pdflatex", "pgf.rcfonts": False,
    "pgf.preamble": r"\usepackage{amsmath,amssymb}",
    "font.family": "serif", "font.serif": [], "font.size": 12,
    "text.usetex": True, "axes.labelsize": 12, "axes.titlesize": 12,
    "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 10,
    "legend.title_fontsize": 10, "axes.linewidth": .65,
    "lines.linewidth": 1.15, "lines.markersize": 4.5,
    "xtick.major.width": .65, "ytick.major.width": .65,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": False, "savefig.pad_inches": .01, "savefig.dpi": 600,
    "figure.facecolor": "white", "axes.facecolor": "white",
})
# The 12pt article class has 390pt text width, increased by one inch.
TEXTWIDTH_IN = (390 + 72.27) / 72.27
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#6C5700", "#333333"]
OUTPUT = Path(os.environ.get("OTTO_FIGURE_OUTPUT", Path(__file__).resolve().parent / "build")).resolve()
def save_figure(fig, stem):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    # Fixed canvas retains the intended final font sizes when included in TeX.
    fig.savefig(OUTPUT / (stem + ".pgf"))
    fig.savefig(OUTPUT / (stem + ".pdf"))
    return OUTPUT / (stem + ".pdf")
