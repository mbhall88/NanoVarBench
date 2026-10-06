"""Figure 4: wall time and peak memory of variant calling against Depth.

Built from benchmarks.tsv, like Table 1. Each Arm's calling step (Clair3 for Arms A-C, dorado
polish on the GPU for Arm D) is a line of medians over the Samples and Read models at each
Depth, with bars for the range. Dorado's timing-only `--device cpu` re-run is an open marker at
the Depths it ran at. Alignment isn't drawn: Table 1 has it. Both y axes are log scale, and
points at the same Depth are nudged apart so their bars don't overlap. Peak memory is the job's
peak RSS, i.e. host memory: GPU memory isn't measured.
"""

import sys
from pathlib import Path

sys.stderr = open(snakemake.log[0], "w")
sys.path.insert(0, str(Path(__file__).parent))

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator  # noqa: E402

from figures_common import COLOURS, MARKERS, ZORDERS, arm_label, save  # noqa: E402

cfg = snakemake.params
bench = pd.read_csv(snakemake.input.benchmarks, sep="\t", dtype={"tool_version": str})
bench["timing_only"] = bench["timing_only"].astype(str).str.lower() == "true"
call = bench[bench["step"] == "call"].assign(max_rss_gb=lambda d: d["max_rss_mb"] / 1000)
depths = sorted(call["depth"].unique())
arms = sorted(call["arm"].unique())

TOOLS = {"clair3": "Clair3", "dorado": "dorado polish"}
METRICS = [("wall_time_s", "Wall time (s)"), ("max_rss_gb", "Peak RAM (GB)")]
TICKS = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]

# One series per Arm and device: the Arm's own calling run, plus Dorado's CPU re-run.
series = []
for (arm, device, timing_only), grp in call.groupby(["arm", "device", "timing_only"]):
    stats = grp.groupby("depth")[[m for m, _ in METRICS]].agg(["median", "min", "max"])
    tool = f"{TOOLS[grp['tool'].iloc[0]]} {grp['tool_version'].iloc[0]}"
    label = f"{arm_label(arm, cfg.arm_labels)}: {tool}, {device.upper()}"
    series.append({"arm": arm, "timing_only": timing_only, "stats": stats, "label": label})
series.sort(key=lambda s: (s["arm"], s["timing_only"]))
nudge = {i: 1.07 ** (i - (len(series) - 1) / 2) for i in range(len(series))}

fig, axes = plt.subplots(1, len(METRICS), figsize=(9.2, 3.9), squeeze=False)
for ax, (metric, ylabel) in zip(axes[0], METRICS):
    values = []
    for i, s in enumerate(series):
        st = s["stats"][metric].dropna()
        x = [d * nudge[i] for d in st.index]
        open_marker = s["timing_only"]
        ax.errorbar(
            x, st["median"], yerr=[st["median"] - st["min"], st["max"] - st["median"]],
            color=COLOURS[s["arm"]], marker=MARKERS[s["arm"]], ms=5.5,
            mfc="white" if open_marker else COLOURS[s["arm"]], ls="none" if open_marker else "-",
            lw=1.6, elinewidth=1, capsize=2.5, zorder=ZORDERS[s["arm"]],
        )  # fmt: skip
        values += [*st["min"], *st["max"]]
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator(depths))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{d}x" for d in depths]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlim(min(depths) / 1.25, max(depths) * 1.25)
    ax.set_yscale("log")
    lo, hi = min(values) / 1.3, max(values) * 1.3
    ax.set_ylim(lo, hi)
    ticks = [t for t in TICKS if lo <= t <= hi]
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in ticks]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("Depth")
    ax.set_ylabel(ylabel)

handles = [
    Line2D(
        [], [], color=COLOURS[s["arm"]], marker=MARKERS[s["arm"]], ms=5.5, lw=1.6,
        ls="none" if s["timing_only"] else "-", mfc="white" if s["timing_only"] else COLOURS[s["arm"]],
        label=s["label"] + (" (timing only)" if s["timing_only"] else ""),
    )
    for s in series
]  # fmt: skip
fig.legend(handles=handles, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0))
fig.tight_layout(rect=(0, 0.2, 1, 1))
save(fig, [snakemake.output.png, snakemake.output.svg], cfg.dpi)
print(f"series={[s['label'] for s in series]} depths={depths}", file=sys.stderr)
