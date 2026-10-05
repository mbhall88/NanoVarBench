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
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator  # noqa: E402

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

# F1, precision and recall go on a logit axis, which spreads out the differences close to 1
# that a linear axis squashes together. A perfect score has no logit, so it is drawn at
# 1 - LOGIT_CLIP: Q60, the same cap as f1_qscore's.
LOGIT_CLIP = 1e-6
LOGIT_TICKS = [
    0.5, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.998, 0.999, 0.9995, 0.9998, 0.9999, 0.99999,
    1 - LOGIT_CLIP,
]  # fmt: skip
LOGIT_NOTE = "F1 is on a logit scale."

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


def logit_clip(values):
    """Scores clipped into the logit scale's domain, so a perfect score can be drawn."""
    return np.clip(np.asarray(values, dtype=float), LOGIT_CLIP, 1 - LOGIT_CLIP)


def logit_axis(ax, which, values, max_ticks=6, perfect_column=False):
    """Put the "x" or "y" axis on a logit scale, limited to the values with a little padding,
    with ticks from LOGIT_TICKS labelled as decimals. Returns the function that maps scores to
    where they are drawn.

    A perfect score is drawn at 1 - LOGIT_CLIP and labelled 1. With `perfect_column`, it is
    instead drawn in a column of its own just past the panel's best imperfect score, after a
    dotted line, so one perfect score doesn't squash the rest of the panel: the axis is broken
    there, and the distance to the column means nothing."""
    v = np.asarray(values, dtype=float)
    logit = lambda x: np.log(x / (1 - x))  # noqa: E731
    expit = lambda z: 1 / (1 + np.exp(-z))  # noqa: E731
    if perfect_column and (v >= 1).any():
        imperfect = logit(logit_clip(v[v < 1])) if (v < 1).any() else logit(np.array([0.9999]))
        lo, top = imperfect.min(), imperfect.max()
        gap = max(0.15 * (top - lo), 0.4)
        perfect = expit(top + gap)
        ax.axvline(expit(top + gap / 2), color="#999999", lw=0.7, ls=":", zorder=1) if which == "x" else ax.axhline(
            expit(top + gap / 2), color="#999999", lw=0.7, ls=":", zorder=1
        )  # fmt: skip
        z = np.append(imperfect, top + gap)
        tick_max = expit(top + gap / 3)
    else:
        perfect = 1 - LOGIT_CLIP
        z = logit(logit_clip(v))
        tick_max = 1
    pad = 0.06 * (z.max() - z.min()) + 0.1
    lim = (expit(z.min() - pad), expit(z.max() + pad))
    getattr(ax, f"set_{which}scale")("logit")
    getattr(ax, f"set_{which}lim")(lim)
    ticks = [t for t in LOGIT_TICKS[:-1] if lim[0] <= t <= min(lim[1], tick_max)]
    while len(ticks) > max_ticks:  # thin from the bottom, keeping the ticks closest to 1
        ticks = ticks[::-2][::-1]
    labels = [f"{t:g}" for t in ticks]
    if lim[0] <= perfect <= lim[1]:
        ticks, labels = ticks + [perfect], labels + ["1"]
    axis = ax.xaxis if which == "x" else ax.yaxis
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(FixedFormatter(labels))
    axis.set_minor_locator(NullLocator())
    return lambda x: np.where(np.asarray(x, dtype=float) >= 1, perfect, logit_clip(x))
