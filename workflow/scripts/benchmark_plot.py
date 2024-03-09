import sys

sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# colourblind-friendly palette from colour universal design (CUD)
# https://jfly.uni-koeln.de/color/
# https://nanx.me/oneclust/reference/cud.html
named_colors = {
    "black": "#000000",
    "orange": "#e69f00",
    "skyblue": "#56b4e9",
    "bluish green": "#009e73",
    "yellow": "#f0e442",
    "blue": "#0072b2",
    "vermilion": "#d55e00",
    "reddish purple": "#cc79a7",
}
cud_palette = list(named_colors.values())
sns.set_theme(style="whitegrid")


def cud(n: int = len(cud_palette), start: int = 0) -> list[str]:
    remainder = cud_palette[:start]
    palette = cud_palette[start:] + remainder
    return palette[:n]


gsize = {}
for p in Path("../results/truth").rglob("*.fai"):
    sample = p.parts[-2]
    size = sum(int(l.split("\t")[1]) for l in p.read_text().splitlines())
    gsize[sample] = size
frames = []
for p in Path("../results/benchmark/call/mutref").rglob("*.tsv"):
    df = pd.read_csv(p, sep="\t")
    sample = p.stem
    model = p.parts[-2].split("_")[-1].split("@")[0]
    mode = p.parts[-4]
    dp = int(p.parts[-5][:-1])
    bp = dp * gsize[sample]
    caller = p.parts[-6]
    df["sample"] = sample
    df["model"] = model
    df["mode"] = mode
    df["depth"] = dp
    df["caller"] = caller
    df["bp"] = bp
    # use rate which is sec/Mbp
    df["rate"] = df["cpu_time"] / df["bp"] * 1e6

    frames.append(df)
df = pd.concat(frames)
df.head()
y = "caller"
hue = y
order = sorted(df[hue].unique())
# move longshot to the end
order.remove("longshot")
order.append("longshot")
palette = cud(n=len(df[hue].unique()))
fig, axes = plt.subplots(nrows=2, figsize=(10, 10), dpi=300, sharey=True)

# plot memory
x = "max_rss"
mem_ax = axes[0]
violin_alpha = 0.2
strip_alpha = 0.2
orient = "h"

kwargs = dict(
    data=df,
    x=x,
    y=y,
    hue=hue,
    order=order,
    hue_order=order,
    palette=palette,
    orient=orient,
    dodge=False,
)
# sns.violinplot(**kwargs, cut=0, inner="quartile", ax=mem_ax)
sns.boxenplot(**kwargs, ax=mem_ax, fill=None, showfliers=False)
# for violin in mem_ax.collections:
# violin.set_facecolor(to_rgba(violin.get_facecolor(), alpha=violin_alpha))

sns.stripplot(
    **kwargs, alpha=strip_alpha, edgecolor="gray", linewidth=0.5, ax=mem_ax, jitter=0.2
)
mem_ax.set_xscale("log")
ticks = [
    (100, "100MB"),
    (500, "500MB"),
    (1000, "1GB"),
    (2000, "2GB"),
    # (3000, "3GB"),
    (4000, "4GB"),
    (8000, "8GB"),
]
FS = 12
mem_ax.set_xticks([t[0] for t in ticks])
mem_ax.set_xticklabels([t[1] for t in ticks], fontsize=FS)
mem_ax.set_xlabel("Max. RAM usage", fontsize=FS)
mem_ax.set_ylabel("")
mem_ax.tick_params(axis="both", which="major", labelsize=FS)

# plot time (rate)
x = "rate"
rt_ax = axes[1]
kwargs["x"] = x

sns.boxenplot(**kwargs, ax=rt_ax, fill=None, showfliers=False)
# sns.violinplot(**kwargs, cut=0, inner="quartile", ax=rt_ax)
# for violin in rt_ax.collections:
#     violin.set_facecolor(to_rgba(violin.get_facecolor(), alpha=violin_alpha))

sns.stripplot(**kwargs, alpha=strip_alpha, edgecolor="gray", linewidth=1, ax=rt_ax)

rt_ax.set_xscale("log")
ticks = [
    (0.1, "0.1s/Mbp"),
    (1, "1s/Mbp"),
    (10, "10s/Mbp"),
    (60, "1m/Mbp"),
    (600, "10m/Mbp"),
]
rt_ax.set_xticks([t[0] for t in ticks])
rt_ax.set_xticklabels([t[1] for t in ticks], fontsize=FS)
rt_ax.set_xlabel("CPU time", fontsize=FS)
rt_ax.set_ylabel("")
rt_ax.tick_params(axis="both", which="major", labelsize=FS)
rt_ax.set_ylabel("")
rt_ax.tick_params(axis="both", which="major", labelsize=FS)
