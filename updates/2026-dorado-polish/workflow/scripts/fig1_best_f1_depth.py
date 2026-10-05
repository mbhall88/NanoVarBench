"""Figure 1: Best F1 against Depth (#14).

One panel per variant type (SNP, INDEL; rows) and Read model (hac, sup; columns). Each Arm
has a solid line for its Best F1 (the `sweep_best` scoring mode) and Arms C and D (the
config's figures.default_pass_arms) also a dashed line for the Default-PASS score. Lines are
medians over the Samples in the run, so Arm D's Default-PASS score sits apart from its Best
F1 and can be compared with Clair3's. With the AF filter analysis enabled (#28) there is a
further series, Clair3 with the AF filter, which gets both lines like Arms C and D. F1 is on a
logit scale, which spreads out the differences close to 1.
"""

import sys
from pathlib import Path

sys.stderr = open(snakemake.log[0], "w")
sys.path.insert(0, str(Path(__file__).parent))

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FixedLocator, NullLocator  # noqa: E402

from figures_common import (  # noqa: E402
    AF_SERIES, COLOURS, DEPTH_NOTE, LINEWIDTHS, LOGIT_NOTE, MARKERS, READ_MODELS, VAR_TYPES,
    ZORDERS, af_note, arm_label, arm_order, logit_axis, logit_clip, present, read_results, save,
)  # fmt: skip

cfg = snakemake.params
af = cfg.af_filter
results = read_results(snakemake.input.results, snakemake.input.af_results[0] if af else None, af)
read_models = present(results["read_model"], READ_MODELS)
arms = arm_order(results["arm"].unique())
depths = sorted(results["depth"].unique())
n_samples = results["sample"].nunique()

# Median over Samples of each Arm's Best F1 / Default-PASS F1 at each Depth.
medians = (
    results.groupby(["scoring_mode", "var_type", "read_model", "arm", "depth"])["f1"]
    .median()
    .reset_index()
)

fig, axes = plt.subplots(
    len(VAR_TYPES), len(read_models), figsize=(3.9 * len(read_models) + 0.4, 6.6),
    sharex=True, sharey="row", squeeze=False,
)  # fmt: skip
for i, var_type in enumerate(VAR_TYPES):
    for j, read_model in enumerate(read_models):
        ax = axes[i, j]
        panel = medians[(medians["var_type"] == var_type) & (medians["read_model"] == read_model)]
        for arm in arms:
            best = panel[(panel["scoring_mode"] == "sweep_best") & (panel["arm"] == arm)]
            ax.plot(
                best["depth"], logit_clip(best["f1"]), color=COLOURS[arm], marker=MARKERS[arm], ms=4.5,
                lw=LINEWIDTHS[arm], alpha=0.9 if arm != "A" else 0.45, zorder=ZORDERS[arm],
                solid_capstyle="round",
            )  # fmt: skip
            if arm in cfg.default_pass_arms or arm == AF_SERIES:
                dp = panel[(panel["scoring_mode"] == "default_pass") & (panel["arm"] == arm)]
                ax.plot(
                    dp["depth"], logit_clip(dp["f1"]), color=COLOURS[arm], marker=MARKERS[arm], ms=4.5,
                    mfc="white", lw=1.6, ls=(0, (4, 2)), zorder=ZORDERS[arm] + 5,
                )  # fmt: skip
        # One logit y range per row (sharey), from the lines drawn in both Read models.
        row = medians[(medians["var_type"] == var_type) & medians["arm"].isin(arms)]
        drawn = (row["scoring_mode"] == "sweep_best") | row["arm"].isin([*cfg.default_pass_arms, AF_SERIES])
        logit_axis(ax, "y", row.loc[drawn, "f1"])
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(FixedLocator(depths))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_xticklabels([f"{d}x" for d in depths])
        if i == 0:
            ax.set_title(f"{read_model} reads")
        if i == len(VAR_TYPES) - 1:
            ax.set_xlabel("Depth")
        if j == 0:
            ax.set_ylabel(f"{var_type} F1")

handles = [
    Line2D([], [], color=COLOURS[a], marker=MARKERS[a], ms=5, lw=2.2, label=arm_label(a, cfg.arm_labels, af))
    for a in arms
]  # fmt: skip
handles += [
    Line2D([], [], color="black", lw=1.8, label="Best F1 (QUAL sweep)"),
    Line2D(
        [], [], color="black", lw=1.6, ls=(0, (4, 2)), marker="o", mfc="white", ms=4.5,
        label="Default-PASS score (Arms " + " and ".join(cfg.default_pass_arms) + (", and the AF filter)" if af else ")"),
    ),
]  # fmt: skip
fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.06 if af else 0.045))
fig.suptitle("Best F1 against Depth", y=0.98, fontsize=11)
note = f"Lines are medians over {n_samples} Samples. {LOGIT_NOTE} {DEPTH_NOTE}"
if af:
    note += "\n" + af_note(af)
fig.text(0.5, 0.012, note, ha="center", va="bottom", fontsize=7.5, color="#444444", wrap=True)
fig.tight_layout(rect=(0, 0.14 if af else 0.11, 1, 0.97))
save(fig, [snakemake.output.png, snakemake.output.svg], cfg.dpi)
print(f"arms={arms} depths={depths} read_models={read_models} samples={n_samples}", file=sys.stderr)
