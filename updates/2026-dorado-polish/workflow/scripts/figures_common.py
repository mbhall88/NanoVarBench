"""Shared loading and styling for the figure and table scripts (#14).

Imported by fig1_best_f1_depth.py, fig2_pr_curves.py, fig3_per_sample.py, table1_runtime.py
and table_s1_per_sample.py, which all read the aggregated tables (results.tsv, pr_curves.tsv,
depth.tsv, benchmarks.tsv) and nothing else.

When the AF filter analysis (#20) is enabled, Figures 1-3 and Table S1 also draw one extra
series from its tables (#28): Clair3 with the AF filter at one AF threshold, which is an extra
analysis and not an Arm. `af` below is the rule's figures.af_filter params, {"arm", "threshold"},
or None when the series is off, in which case nothing here changes the output.
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
# The AF filter series sits in the Arm dictionaries under its own key, apart from any Arm
# letter. It is bluish green, the Okabe-Ito colour that isn't a blue, grey or orange.
AF_SERIES = "AF"
COLOURS[AF_SERIES] = "#009e73"
MARKERS[AF_SERIES] = "P"
# Arms A-C overlap almost exactly at most Depths, so each gets a different line width and
# the thinner lines on top stay visible through the wider ones under them.
LINEWIDTHS = {"A": 5.0, "B": 3.2, "C": 1.8, "D": 1.8, AF_SERIES: 1.8}
ZORDERS = {"A": 2, "B": 3, "C": 4, "D": 5, AF_SERIES: 6}
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


def arm_label(arm, labels, af=None):
    """The Arm's name in CONTEXT.md, e.g. "Arm D (Dorado)". The AF filter series is "Arm C + AF
    filter (0.65), extra analysis": the threshold is there, and it says it isn't an Arm."""
    if arm == AF_SERIES:
        return f"Arm {af['arm']} + AF filter ({af['threshold']}), extra analysis"
    return f"Arm {arm} ({labels[arm]})"


def af_note(af):
    """One sentence for a figure's note saying what the AF filter series is."""
    return (
        f"Arm {af['arm']} + AF filter ({af['threshold']}) is an extra analysis, not an Arm: "
        f"Arm {af['arm']}'s Clair3 run diploid, each het call resolved by FORMAT/AF >= {af['threshold']}."
    )


def arm_order(arms):
    """The Arms in letter order, then the AF filter series, if it is among them."""
    return sorted(a for a in arms if a != AF_SERIES) + [a for a in arms if a == AF_SERIES]


def read_results(path, af_path=None, af=None):
    """results.tsv with typed columns, restricted to the SNP and INDEL rows the figures use.
    With the AF filter series (`af`, from clair3_af_filter.tsv at `af_path`), its rows are added
    with the Arm AF_SERIES: the same columns as results.tsv's, Best F1 and Default-PASS."""
    df = pd.read_csv(path, sep="\t", dtype={"dnd_sample": str})
    if af:
        df = pd.concat([df, read_af_series(af_path, af)], ignore_index=True)
    df["dnd_sample"] = df["dnd_sample"].astype(str).str.lower() == "true"
    return df[df["var_type"].isin(VAR_TYPES)]


def read_af_series(path, af, as_text=False):
    """The AF filter's rows for the series' Arm and AF threshold from clair3_af_filter.tsv
    (its `calls == af_filter` rows), as results.tsv-style rows with the Arm AF_SERIES. With
    `as_text`, every column is left as the file wrote it."""
    dtype = str if as_text else {"dnd_sample": str, "af_threshold": str}
    df = pd.read_csv(path, sep="\t", dtype=dtype)
    df = df[(df["calls"] == "af_filter") & (df["arm"] == af["arm"]) & (df["af_threshold"] == af["threshold"])]
    if df.empty:
        raise ValueError(f"no AF filter rows for Arm {af['arm']} at {af['threshold']} in {path}")
    return df.drop(columns=["calls", "af_threshold", "clair3_options"]).assign(arm=AF_SERIES)


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
