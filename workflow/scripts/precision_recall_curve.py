import sys

sys.stderr = open(snakemake.log[0], "w")
import math
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
sns.set_theme(style="whitegrid")
no_indels = snakemake.params.no_indels
MODEL = snakemake.wildcards.model.split("_")[-1].split("@")[0]


def cud(n: int = len(cud_palette), start: int = 0) -> list[str]:
    remainder = cud_palette[:start]
    palette = cud_palette[start:] + remainder
    return palette[:n]


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
        mode = p.parts[-5]
        df["mode"] = mode
        df["version"] = p.parts[-4]
        model = p.parts[-3].split("_")[-1].split("@")[0]
        assert model == MODEL, (p, model, snakemake.wildcards.model)
        if mode == "duplex" and model == "fast":
            continue
        df["model"] = model
    frames.append(df)
    # duplicate Illumina as duplex also
    if "illumina" in str(p) and MODEL != "fast":
        df2 = df.copy()
        df2["caller"] = "illumina"
        df2["depth"] = "illumina"
        df2["mode"] = "duplex"
        df2["version"] = "illumina"
        df2["model"] = "illumina"
        frames.append(df2)

pr_df = pd.concat(frames)
pr_df.reset_index(inplace=True)
model = MODEL

samples = set(pr_df["sample"])
metrics = []
for vartype in ("SNP", "INDEL"):
    for mode in ("simplex", "duplex"):
        for caller in pr_df["caller"].unique():
            if caller in no_indels and vartype == "INDEL":
                continue

            if caller != "illumina":
                data = pr_df.query(
                    f"caller == '{caller}' and VAR_TYPE == '{vartype}' and mode == '{mode}' and model == '{model}'"
                )
            else:
                data = pr_df.query(
                    f"caller == '{caller}' and VAR_TYPE == '{vartype}' and mode == '{mode}' and model == 'illumina'"
                )
            for q in sorted(set(data["MIN_QUAL"])):
                # check if every sample has this MIN_QUAL
                subdf = data.query(f"MIN_QUAL == {q}")
                if set(subdf["sample"]) == samples:
                    tps = subdf["TRUTH_TP"].sum()
                    fps = subdf["QUERY_FP"].sum()
                    fns = subdf["TRUTH_FN"].sum()
                    precision = tps / (tps + fps)
                    recall = tps / (tps + fns)
                    f1 = 2 * (precision * recall) / (precision + recall)
                    metrics.append(
                        (caller, q, precision, recall, f1, vartype, mode, model)
                    )


aggdf = pd.DataFrame(
    metrics,
    columns=["caller", "QUAL", "precision", "recall", "f1", "vartype", "mode", "model"],
)

vartypes = ["SNP", "INDEL"]
modes = ["simplex"]

if MODEL != "fast":
    modes.append("duplex")

fig, axes = plt.subplots(
    nrows=len(vartypes),
    ncols=len(modes),
    figsize=(10, 10),
    dpi=300,
    sharex=True,
    sharey=True,
)
x = "recall"
y = "precision"
hue = "caller"

hue_order = sorted(set(aggdf[hue]))
# move illumina to end
hue_order.remove("illumina")
hue_order.append("illumina")
# map each caller in hue order to a colour in cud
pal = {c: cud()[i] for i, c in enumerate(hue_order)}

i = 0
legend = True
for vartype in vartypes:
    for mode in modes:
        ax = axes.flatten()[i]
        data = aggdf.query("vartype == @vartype and mode == @mode and model == @model")
        cap = 0.99999
        data.loc[:, y] = data[y].apply(lambda v: cap if v > cap else v)

        if vartype == "INDEL":
            data = data.query("caller not in @no_indels")

        sns.lineplot(
            data=data,
            x=x,
            y=y,
            hue=hue,
            hue_order=sorted(pal),
            ax=ax,
            palette=pal,
            alpha=0.9,
            linewidth=2,
            legend=legend,
        )

        if legend:
            handles, labels = ax.get_legend_handles_labels()
            ax.legend().remove()
            legend = False

        ax.set_yscale("logit", nonpositive="clip")
        yticks = [0.8, 0.9, 0.95, 0.99, 0.999, 0.9999, cap]
        yticklabels = [f"{yval:.2%}" for yval in yticks]
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)

        xticks = [0, 0.25, 0.5, 0.75, 1.0]
        xticklabels = [f"{xval:.2%}" for xval in xticks]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        ax.set_title(f"{vartype} {mode}")
        i += 1

# make legend lines thicker
for h in handles:
    h.set_linewidth(3)

plt.tight_layout()
leg_cols = math.ceil((len(aggdf[hue].unique()) + 1) / 2)
fig.legend(
    handles=handles,
    labels=labels,
    loc="upper center",
    bbox_to_anchor=(0.52, 0.59),
    ncol=leg_cols,
    title="",
    framealpha=1.0,
    fancybox=True,
    shadow=True,
)

fig.savefig(snakemake.output.pdf)
