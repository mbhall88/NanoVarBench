"""Figure 2: precision-recall curves at the Depths in the config's figures.pr_depths (#14).

One panel per variant type (SNP, INDEL; rows) x Read model x Depth (columns), with a curve
per Arm from the QUAL sweep (pr_curves.tsv). The Default-PASS point (the PASS-only score,
no QUAL threshold) is marked on each curve. Each curve pools the Samples in the run: truth
and query counts are summed over the Samples at each QUAL threshold, so it is the PR curve of
all their variants together. A Sample whose calls end below a threshold contributes no calls
there: all its truth variants count as missed. The axes are zoomed to the curves' upper
right corner, so each panel has its own limits. With the AF filter analysis enabled (#28) there
is a further curve and Default-PASS point for Clair3 with the AF filter, from its own curves
table (the AF filter's QUAL sweep for the one threshold the figure shows).
"""

import sys
from pathlib import Path

sys.stderr = open(snakemake.log[0], "w")
sys.path.insert(0, str(Path(__file__).parent))

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

from figures_common import (  # noqa: E402
    AF_SERIES, COLOURS, DEPTH_NOTE, MARKERS, READ_MODELS, VAR_TYPES, ZORDERS, af_note, arm_label,
    arm_order, present, read_results, save,
)  # fmt: skip

cfg = snakemake.params
af = cfg.af_filter
COUNTS = ["truth_tp", "truth_fn", "query_tp", "query_fp"]

results = read_results(snakemake.input.results, snakemake.input.af_results[0] if af else None, af)
curves = pd.read_csv(snakemake.input.pr_curves, sep="\t")
if af:
    af_curves = pd.read_csv(snakemake.input.af_curves[0], sep="\t", dtype={"af_threshold": str})
    af_curves = af_curves[
        (af_curves["arm"] == af["arm"]) & (af_curves["af_threshold"] == af["threshold"])
    ]
    curves = pd.concat([curves, af_curves.drop(columns="af_threshold").assign(arm=AF_SERIES)])
curves = curves[curves["var_type"].isin(VAR_TYPES)]
read_models = present(results["read_model"], READ_MODELS)
arms = arm_order(results["arm"].unique())
depths = [d for d in cfg.depths if d in set(results["depth"])]
for d in cfg.depths:
    if d not in depths:
        print(f"skipping {d}x: not in the results", file=sys.stderr)
if not depths:
    raise ValueError(f"none of the Depths {list(cfg.depths)} are in the results")


def pr(counts):
    """Precision and recall from summed counts; precision is undefined with no calls."""
    calls = counts["query_tp"] + counts["query_fp"]
    truth = counts["truth_tp"] + counts["truth_fn"]
    return pd.DataFrame(
        {
            "precision": (counts["query_tp"] / calls).where(calls > 0),
            "recall": counts["truth_tp"] / truth,
        }
    )


def pooled_curve(grp):
    """Sum each Sample's counts at every QUAL threshold from 0 to the largest in the group.
    Beyond a Sample's last threshold it has no calls: no TP or FP, and all its truth missed."""
    thresholds = range(int(grp["min_qual"].max()) + 1)
    total = 0
    for _, g in grp.groupby("sample"):
        g = g.set_index("min_qual").reindex(thresholds)
        g["truth_fn"] = g["truth_fn"].fillna(g["truth_total"].max())
        total = total + g[COUNTS].fillna(0)
    return pr(total).dropna()


def pooled_point(rows):
    return pr(rows[COUNTS].sum().to_frame().T).iloc[0]


fig, axes = plt.subplots(
    len(VAR_TYPES), len(read_models) * len(depths),
    figsize=(2.9 * len(read_models) * len(depths) + 0.4, 6.4), squeeze=False,
)  # fmt: skip
for i, var_type in enumerate(VAR_TYPES):
    for j, (read_model, depth) in enumerate((rm, d) for rm in read_models for d in depths):
        ax = axes[i, j]
        xs, ys = [], []  # points the zoom has to include
        for arm in arms:
            sel = (curves["read_model"] == read_model) & (curves["depth"] == depth)
            sel &= (curves["arm"] == arm) & (curves["var_type"] == var_type)
            curve = pooled_curve(curves[sel])
            ax.plot(
                curve["recall"], curve["precision"], color=COLOURS[arm], lw=1.6,
                zorder=ZORDERS[arm], clip_on=True,
            )  # fmt: skip
            keep = (results["read_model"] == read_model) & (results["depth"] == depth)
            keep &= (results["arm"] == arm) & (results["var_type"] == var_type)
            point = pooled_point(results[keep & (results["scoring_mode"] == "default_pass")])
            best = pooled_point(results[keep & (results["scoring_mode"] == "sweep_best")])
            ax.plot(
                [point["recall"]], [point["precision"]], marker=MARKERS[arm], ms=6.5,
                color=COLOURS[arm], mec="black", mew=0.8, ls="none", zorder=ZORDERS[arm] + 5,
            )  # fmt: skip
            # Zoom to the Default-PASS point, the Best F1 point and the unthresholded end.
            xs += [point["recall"], best["recall"], curve["recall"].iloc[0]]
            ys += [point["precision"], best["precision"], curve["precision"].iloc[0]]
        pad_x = 0.25 * (max(xs) - min(xs)) + 0.002
        pad_y = 0.25 * (max(ys) - min(ys)) + 0.002
        # A little room above 1 so a curve along precision 1 isn't clipped; no tick is drawn there.
        xlo, ylo = min(xs) - pad_x, min(ys) - pad_y
        ax.set_xlim(xlo, min(max(xs) + pad_x, 1.0 + 0.03 * (1.0 - xlo)))
        ax.set_ylim(ylo, min(max(ys) + pad_y, 1.0 + 0.03 * (1.0 - ylo)))
        xlim, ylim = ax.get_xlim(), ax.get_ylim()  # set_ticks can widen the limits
        for axis, set_ticks in ((ax.xaxis, ax.set_xticks), (ax.yaxis, ax.set_yticks)):
            axis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))
            set_ticks([t for t in axis.get_majorticklocs() if t <= 1.0 + 1e-9])
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_title(f"{var_type}, {read_model} reads, {depth}x", fontsize=9)
        ax.set_xlabel("Recall")
        if j == 0:
            ax.set_ylabel("Precision")

handles = [
    Line2D([], [], color=COLOURS[a], lw=1.8, label=arm_label(a, cfg.arm_labels, af)) for a in arms
]  # fmt: skip
handles.append(
    Line2D([], [], color="white", mfc="#bbbbbb", mec="black", marker="o", ms=6.5, label="Default-PASS score (one marker per Arm)")
)  # fmt: skip
fig.legend(
    handles=handles, loc="lower center", ncol=3 if af else len(handles),
    bbox_to_anchor=(0.5, -0.045 if af else -0.01),
)
fig.suptitle("Precision-recall curves (QUAL sweep)", y=0.99, fontsize=11)
note = f"Curves pool all {results['sample'].nunique()} Samples. Each panel is zoomed to its own range. {DEPTH_NOTE}"
if af:
    note += "\n" + af_note(af)
fig.text(0.5, -0.045 if not af else -0.08, note, ha="center", va="top", fontsize=7.5, color="#444444")
fig.tight_layout(rect=(0, 0.05, 1, 0.97))
save(fig, [snakemake.output.png, snakemake.output.svg], cfg.dpi)
print(f"arms={arms} depths={depths} read_models={read_models}", file=sys.stderr)
