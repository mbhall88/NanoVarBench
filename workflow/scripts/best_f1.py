import sys

sys.stderr = open(snakemake.log[0], "w")
from pathlib import Path

import matplotlib.pyplot as plt
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


def cud(n: int = len(cud_palette), start: int = 0) -> list[str]:
    remainder = cud_palette[:start]
    palette = cud_palette[start:] + remainder
    return palette[:n]


sns.set_theme(style="whitegrid")
no_indels = snakemake.params.no_indels

output = {}
for path in snakemake.output.figures:
    p = Path(path)
    metric = p.stem
    if metric == "f1":
        metric = "F1_SCORE"
    elif metric == "precision":
        metric = "PREC"
    elif metric == "recall":
        metric = "RECALL"
    else:
        raise ValueError(f"Unknown metric: {metric}")
    output[metric] = p

tsvs = snakemake.input.pr
tsvs.extend(snakemake.input.illumina_pr)

frames = []

for p in map(Path, tsvs):
    df = pd.read_csv(p, sep="\t")
    df["sample"] = p.parent.name

    if "illumina" in str(p):
        df["caller"] = "illumina"
        df["depth"] = "illumina"
        df["mode"] = "simplex"
        df["version"] = "illumina"
        df["model"] = "illumina"
    else:
        df["caller"] = p.parts[-7]
        depth = int(p.parts[-6][:-1])
        df["depth"] = depth
        df["mode"] = p.parts[-5]
        df["version"] = p.parts[-4]
        model = p.parts[-3].split("_")[-1].split("@")[0]
        df["model"] = model
    frames.append(df)
    # duplicate Illumina as duplex also
    if "illumina" in str(p):
        df2 = df.copy()
        df2["caller"] = "illumina"
        df2["depth"] = "illumina"
        df2["mode"] = "duplex"
        df2["version"] = "illumina"
        df2["model"] = "illumina"
        frames.append(df2)

pr_df = pd.concat(frames)
pr_df.reset_index(inplace=True, drop=True)

metrics = ["F1_SCORE", "PREC", "RECALL"]
x = "caller"
hue = "model"
cols = "mode"
ncols = pr_df[cols].nunique()
rows = ("SNP", "INDEL")
nrows = len(rows)
dataix = pr_df.groupby([x, hue, cols, "VAR_TYPE", "sample"])["F1_SCORE"].idxmax()
data = pr_df.iloc[dataix]
order = sorted(set(data[x]))
# move illumina to end
order.remove("illumina")
order.append("illumina")
hue_order = sorted(set(data[hue]))
# move illumina to end
hue_order.remove("illumina")
hue_order.append("illumina")
# map each caller in hue order to a colour in cud
pal = {c: cud()[i] for i, c in enumerate(hue_order)}

for y in metrics:
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(10, 10),
        dpi=300,
        sharex=True,
        sharey=True,
    )
    i = 0
    for vartype in rows:
        for mode in reversed(data[cols].unique()):
            legend = i == 0
            ax = axes.flatten()[i]
            # transform precision so logit shows values of 1.0 - i.e. change values of 1.0 to 0.99999
            df = data.query("VAR_TYPE == @vartype and mode == @mode")
            if mode == "duplex":
                df = df.query("model != 'fast'")
            if vartype == "INDEL":
                df = df.query("caller not in @no_indels")

            cap = 0.99999
            df.loc[:, y] = df[y].apply(lambda v: cap if v > cap else v)
            yticks = [0.01, 0.1, 0.5, 0.8, 0.9, 0.99, 0.999, 0.9999, cap]
            yticklabels = [f"{yval:.2%}" for yval in yticks]
            sns.boxplot(
                data=df,
                x=x,
                y=y,
                order=order,
                hue=hue,
                ax=ax,
                palette=pal,
                fill=False,
                fliersize=0,
                legend=legend,
                gap=0.2,
            )
            sns.stripplot(
                data=df,
                x=x,
                y=y,
                order=order,
                hue=hue,
                ax=ax,
                palette=pal,
                alpha=0.5,
                dodge=True,
                legend=False,
                linewidth=0.5,
                edgecolor="black",
            )

            ax.set_yscale("logit", nonpositive="clip")
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticklabels)
            # make ylabels more readable
            ylabel = {
                "F1_SCORE": "F1 score",
                "PREC": "Precision",
                "RECALL": "Recall",
            }[y]
            ax.set_ylabel(f"{vartype} {ylabel}")
            # rotate the xlabels
            ax.tick_params(axis="x", labelrotation=90, labelsize=12)
            ax.set_xlabel("")

            if legend:
                handles, labels = ax.get_legend_handles_labels()

                for h in handles:
                    h.set_linewidth(3)

                ax.legend(
                    handles=handles,
                    labels=labels,
                    framealpha=1.0,
                    fancybox=True,
                    shadow=True,
                    title=hue.capitalize(),
                )

                legend = False

            if i < ncols:
                ax.set_title(mode)

            i += 1

    fig.tight_layout()
    fig.savefig(output[y])

dataix = pr_df.groupby([x, hue, cols, "VAR_TYPE", "sample"])["F1_SCORE"].idxmax()
data = pr_df.iloc[dataix]
data = data.query(
    "VAR_TYPE not in ('ALL', 'SV') and not (mode == 'duplex' and model == 'fast')"
)
# make the mode and model columns the first and second columns
col = data.pop("model")
data.insert(0, col.name, col)
col = data.pop("mode")
data.insert(0, col.name, col)
col = data.pop("sample")
data.insert(0, col.name, col)
col = data.pop("caller")
data.insert(0, col.name, col)
data.sort_values(
    by=["mode", "model", "sample", "caller"],
    ascending=[False, True, True, True],
    inplace=True,
)
data.to_csv(snakemake.output.csv, index=False)
