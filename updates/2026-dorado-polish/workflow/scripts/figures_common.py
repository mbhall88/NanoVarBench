"""Shared loading and styling for the figure and table scripts (#14).

Imported by fig1_best_f1_depth.py, fig2_pr_curves.py, fig3_per_sample.py, table1_runtime.py
and table_s1_per_sample.py, which all read the aggregated tables (results.tsv, pr_curves.tsv,
depth.tsv, benchmarks.tsv) and nothing else.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

VAR_TYPES = ["SNP", "INDEL"]  # ALL is in results.tsv but not in the figures
READ_MODELS = ["hac", "sup"]
# Okabe-Ito colours, which hold up for colour-blind readers. Arms A-C all run Clair3 so they
# are shades of blue and grey; Dorado (Arm D) is the orange that has to stand out.
COLOURS = {"A": "#7f7f7f", "B": "#56b4e9", "C": "#0072b2", "D": "#d55e00"}
MARKERS = {"A": "o", "B": "s", "C": "^", "D": "D"}
# Arms A-C overlap almost exactly at most Depths, so each gets a different line width and
# the thinner lines on top stay visible through the wider ones under them.
LINEWIDTHS = {"A": 5.0, "B": 3.2, "C": 1.8, "D": 1.8}
ZORDERS = {"A": 2, "B": 3, "C": 4, "D": 5}
DND_COLOUR = "#b2182b"
DND_URL = "https://github.com/nanoporetech/dorado/issues/1599"

DEPTH_NOTE = (
    "Depth is a per-position cap applied with rasusa aln, not a random genome-wide subsample, "
    "so low-Depth recall is not directly comparable with the eLife paper."
)

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "axes.facecolor": "white",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": "#e5e5e5",
        "grid.linewidth": 0.6,
        "axes.axisbelow": True,
        "legend.frameon": False,
        "svg.fonttype": "none",  # keep text as text in the SVGs
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


def arm_label(arm, labels):
    """The Arm's name in CONTEXT.md, e.g. "Arm D (Dorado)"."""
    return f"Arm {arm} ({labels[arm]})"


def read_results(path):
    """results.tsv with typed columns, restricted to the SNP and INDEL rows the figures use."""
    df = pd.read_csv(path, sep="\t", dtype={"dnd_sample": str})
    df["dnd_sample"] = df["dnd_sample"].str.lower() == "true"
    return df[df["var_type"].isin(VAR_TYPES)]


def present(values, wanted):
    """The wanted values that occur in the data, in the wanted order."""
    return [v for v in wanted if v in set(values)]


def save(fig, outputs, dpi):
    """Write the figure to every output path (PNG and SVG), by the path's extension."""
    for path in outputs:
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def species_short(species):
    """Escherichia coli -> E. coli."""
    genus, *rest = species.split()
    return f"{genus[0]}. {' '.join(rest)}"
