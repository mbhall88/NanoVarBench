import sys

sys.stderr = open(snakemake.log[0], "w")
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

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


frames = []
modes = ["simplex", "duplex"]
models = ["fast", "hac", "sup"]
for p in snakemake.input.csvs:
    path = Path(p)
    model = path.parts[-1].split("_")[4].split("@")[0]
    mode = path.parts[-3]
    if model == "fast" and mode == "duplex":
        continue

    df = pd.read_csv(path)
    df["mode"] = mode
    df["model"] = model
    frames.append(df)

df = pd.concat(frames)
# rename filename column to sample
df.rename(columns={"filename": "sample"}, inplace=True)
# make the mode and model columns the first and second columns
col = df.pop("model")
df.insert(0, col.name, col)
col = df.pop("mode")
df.insert(0, col.name, col)
df.sort_values(
    by=["mode", "model", "sample"], ascending=[False, True, True], inplace=True
)

df.to_csv(snakemake.output.csv, index=False)


# functions to convert between identity and score
def identity2qscore(x):
    a = np.array(x)
    return -10 * np.log10(1 - a)


def qscore2identity(x):
    a = np.array(x)
    return 1 - np.power(10, -a / 10)


fig, ax = plt.subplots(1, 1, figsize=(7, 4), dpi=300, sharex=False, sharey=True)
palette = cud(n=len(df["model"].unique()))
for i, x in enumerate(["median_identity"]):
    if i == 0:
        legend = "brief"
    else:
        legend = False
    # ax = axes[i]
    if x == "median_identity" and df[x].max() > 1:
        df[x] = df[x] / 100

    sns.boxplot(
        data=df,
        y="mode",
        x=x,
        hue="model",
        ax=ax,
        showfliers=False,
        fill=None,
        gap=0.2,
        legend=legend,
        palette=palette,
    )
    sns.stripplot(
        data=df,
        y="mode",
        x=x,
        hue="model",
        ax=ax,
        dodge=True,
        linewidth=1,
        alpha=0.5,
        edgecolor="black",
        legend=False,
        palette=palette,
    )

    ax.set_ylabel("Read type")
    if x == "median_identity":
        ax.set_xscale("logit")
        xticks = [10, 15, 20, 25, 30, 35, 40]
        xticklabels = [f"{qscore2identity(xval):.2%}\n (Q{xval})" for xval in xticks]
        ax.set_xticks(qscore2identity(xticks))
        ax.set_xticklabels(xticklabels)
        ax.set_xlabel("median read identity (Qscore)")

    if legend:
        handles, labels = ax.get_legend_handles_labels()
        for h in handles:
            h.set_linewidth(3)
        ax.legend(
            handles=handles,
            labels=labels,
            title="Model",
            framealpha=1.0,
            fancybox=True,
            shadow=True,
        )

fig.tight_layout()
fig.savefig(snakemake.output.pdf)
