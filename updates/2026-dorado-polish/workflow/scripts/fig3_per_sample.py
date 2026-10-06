"""Figure 3: per-Sample Best F1 at every Depth, with the dnd Samples highlighted (#14).

A dot plot with a row per Sample and a dot per Arm (Best F1, the `sweep_best` scoring mode).
Columns are Depths and rows are variant type x Read model. The dnd Samples, where Dorado's
bacterial model is reported to make systematic errors (dorado#1599), are shaded and their
names are in red. Each panel has its own x limits, because the Depths differ so much in F1,
and F1 is on a logit scale.
With the AF filter analysis enabled (#28) there is a further dot per Sample for Clair3 with the
AF filter.
"""

import sys
from pathlib import Path

sys.stderr = open(snakemake.log[0], "w")
sys.path.insert(0, str(Path(__file__).parent))

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from figures_common import (  # noqa: E402
    COLOURS, DND_COLOUR, MARKERS, READ_MODELS, VAR_TYPES, arm_label, arm_order, logit_axis,
    present, read_results, save, species_short,
)  # fmt: skip

cfg = snakemake.params
af = cfg.af_filter
results = read_results(snakemake.input.results, snakemake.input.af_results[0] if af else None, af)
best = results[results["scoring_mode"] == "sweep_best"]
read_models = present(best["read_model"], READ_MODELS)
arms = arm_order(best["arm"].unique())
depths = sorted(best["depth"].unique())

# dnd Samples first, then the rest by species; the first row is at the top.
samples = (
    best[["sample", "species", "dnd_sample"]]
    .drop_duplicates()
    .sort_values(["dnd_sample", "species"], ascending=[False, True])
    .reset_index(drop=True)
)
ypos = {sample: i for i, sample in enumerate(samples["sample"])}
step = 0.2 if len(arms) <= 4 else 0.16  # keep a Sample's dots inside its row
dodge = {arm: (k - (len(arms) - 1) / 2) * step for k, arm in enumerate(arms)}

rows = [(vt, rm) for rm in read_models for vt in VAR_TYPES]
fig, axes = plt.subplots(
    len(rows), len(depths), figsize=(3.0 * len(depths) + 1.6, 3.6 * len(rows)),
    sharey=True, squeeze=False,
)  # fmt: skip
for i, (var_type, read_model) in enumerate(rows):
    for j, depth in enumerate(depths):
        ax = axes[i, j]
        panel = best[
            (best["var_type"] == var_type) & (best["read_model"] == read_model) & (best["depth"] == depth)
        ]  # fmt: skip
        for _, s in samples[samples["dnd_sample"]].iterrows():
            ax.axhspan(ypos[s["sample"]] - 0.5, ypos[s["sample"]] + 0.5, color=DND_COLOUR, alpha=0.09, lw=0)
        drawn_at = logit_axis(ax, "x", panel["f1"], max_ticks=4, perfect_column=True)
        for arm in arms:
            a = panel[panel["arm"] == arm]
            ax.scatter(
                drawn_at(a["f1"]), [ypos[s] + dodge[arm] for s in a["sample"]], s=16, marker=MARKERS[arm],
                color=COLOURS[arm], edgecolor="white", linewidth=0.3, zorder=3,
            )  # fmt: skip
        ax.set_ylim(len(samples) - 0.5, -0.5)
        ax.set_title(f"{var_type}, {read_model} reads, {depth}x", fontsize=9)
        ax.set_xlabel("Best F1")
        ax.grid(axis="y", visible=False)
        ax.tick_params(axis="x", labelsize=8)
        if j == 0:
            ax.set_yticks(range(len(samples)))
            ax.set_yticklabels(
                [species_short(sp) for sp in samples["species"]], fontstyle="italic", fontsize=8
            )  # fmt: skip
            for label, dnd in zip(ax.get_yticklabels(), samples["dnd_sample"]):
                if dnd:
                    label.set_color(DND_COLOUR)
                    label.set_fontweight("bold")

handles = [
    Line2D([], [], color=COLOURS[a], marker=MARKERS[a], ls="none", ms=6, label=arm_label(a, cfg.arm_labels, af))
    for a in arms
]  # fmt: skip
handles.append(Patch(color=DND_COLOUR, alpha=0.25, label=f"dnd Sample (dorado#1599)"))
fig.legend(
    handles=handles, loc="lower center", ncol=3 if af else len(handles),
    bbox_to_anchor=(0.5, 0),
)
fig.tight_layout(rect=(0, 0.035 if af else 0.025, 1, 1))
save(fig, [snakemake.output.png, snakemake.output.svg], cfg.dpi)
print(f"arms={arms} depths={depths} read_models={read_models} samples={len(samples)}", file=sys.stderr)
