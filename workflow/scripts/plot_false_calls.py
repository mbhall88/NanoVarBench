import sys

sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import seaborn as sns
from matplotlib.colors import SymLogNorm

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
HOM_THRESH = 5


def cud(n: int = len(cud_palette), start: int = 0) -> list[str]:
    remainder = cud_palette[:start]
    palette = cud_palette[start:] + remainder
    return palette[:n]


palette = {
    "TP": named_colors["blue"],
    "FP": named_colors["orange"],
    "FN": named_colors["vermilion"],
}

df = pd.read_csv(snakemake.input.table, low_memory=False)
# create new column, indel_len, which is the absolute value of the difference in length between the ref and alt columns
df["indel_len"] = abs(df["ref"].str.len() - df["alt"].str.len())

frames = []
tsvs = snakemake.input.pr_wo_repeats
tsvs.extend(snakemake.input.pr_w_repeats)
tsvs.extend(snakemake.input.pr_illumina_wo_repeats)
tsvs.extend(snakemake.input.pr_illumina_w_repeats)

for p in map(Path, tsvs):
    _df = pd.read_csv(p, sep="\t")
    _df["sample"] = p.parent.name
    mask = "without" in p.name
    if "illumina" in str(p):
        caller = "illumina"
    else:
        caller = p.parts[-7]
        if caller != "clair3":
            continue
        dp = int(p.parts[-6][:-1])
        if dp < 100:
            continue
        mode = p.parts[-5]
        if mode == "duplex":
            continue
        model = p.parts[-3].split("_")[-1].split("@")[0]
        if model != "sup":
            continue
    _df["mask_repeats"] = mask
    _df["caller"] = caller
    frames.append(_df)

mask_df = pd.concat(frames)
mask_df.reset_index(inplace=True, drop=True)
dataix = (
    mask_df.query("VAR_TYPE == 'ALL'")
    .groupby(["sample", "caller", "mask_repeats"])["F1_SCORE"]
    .idxmax()
)
best_mask = mask_df.iloc[dataix]


def plot_fp_heatmaps():
    for caller in df["caller"].unique():
        if caller in ["illumina", "longshot"]:
            continue
        # create a heatmap of the indel_len vs. homlen columns
        max_indel_len = df.query(
            "vartype=='INDEL' and mode == 'simplex' and caller==@caller"
        )["indel_len"].max()
        max_homlen = df.query(
            "vartype=='INDEL' and mode == 'simplex' and caller==@caller"
        )["homlen"].max()
        # set vmin and vmax to the min and max of the number of FP variants
        vmax = (
            df.query(
                "vartype=='INDEL' and mode == 'simplex' and decision == 'FP' and caller == @caller"
            )
            .groupby(["caller", "model", "homlen", "indel_len"])
            .count()["sample"]
            .max()
        )

        fig, axes = plt.subplots(
            figsize=(10, 10), dpi=300, nrows=2, ncols=2, sharex=True
        )

        for i, model in enumerate(["fast", "hac", "sup", "illumina"]):
            ax = axes.flatten()[i]
            if model == "illumina":
                data = df.query("caller == 'illumina'")
            else:
                data = df.query("caller == @caller")
                data = data.query("model == @model and mode == 'simplex'")
            data = data.query("vartype == 'INDEL' and decision == 'FP'")

            # data = data.query(f"homlen >= {HOM_THRESH}")
            x = "indel_len"
            y = "homlen"
            # turn data into a pivot table where the cells are the counts of the number of variants
            # in each bin
            pivot = data.pivot_table(index=x, columns=y, aggfunc="size", fill_value=0)
            total_fps = data.query("decision == 'FP'").shape[0]
            # calculate the percentage of the total number of FP variants in each bin
            # pivot = pivot / total_fps * 100
            # fill in the missing columns and rows
            for j in range(1, max_indel_len + 1):
                if j not in pivot.index:
                    pivot.loc[j] = 0
            for j in range(1, max_homlen + 1):
                if j not in pivot.columns:
                    pivot[j] = 0
            pivot.sort_index(axis=0, inplace=True)
            pivot.sort_index(axis=1, inplace=True)

            # create a heatmap of the pivot table
            mask = pivot == 0

            fmt = "d"
            cbar_ticks = [0, 1, 10, 50, 100, 500, 1000, 2000]
            cbar_kws = {
                "ticks": cbar_ticks,
                "format": mticker.FixedFormatter(cbar_ticks),
                "label": "Number of FP indels" if i in [1, 3] else "",
            }
            sns.heatmap(
                pivot,
                ax=ax,
                cmap="coolwarm",
                cbar_kws=cbar_kws,
                # square=True,
                # annot=True,
                linecolor="black",
                linewidths=0.5,
                cbar=True,
                fmt=fmt,
                # mask=mask,
                # vmin=vmin,
                # vmax=vmax,
                norm=SymLogNorm(vmax=vmax, linthresh=1),
            )
            ax.set_title(model)
            if i < 2:
                ax.set_xlabel("")
            else:
                ax.set_xlabel("Homopolymer length (bp)")
            if i in [1, 3]:
                ax.set_ylabel("")
            else:
                ax.set_ylabel("Indel length (bp)")

            # turn off gridlines
            ax.grid(False)

            # add a vertical dashed line at the HOM_THRESH and align to left

            ax.axvline(
                HOM_THRESH - 1,
                color="red",
                linestyle="-",
                snap=True,
                label="Homopolymer threshold",
            )
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
            # rotate xticklabels
            ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

        fig.tight_layout()
        handles, labels = ax.get_legend_handles_labels()
        for h in handles:
            h.set_linewidth(2)
        fig.legend(
            handles=handles,
            labels=labels,
            loc="center",
            fancybox=True,
            shadow=True,
            bbox_to_anchor=(0.51, 0.515),
        )
        outfile = ""
        for f in snakemake.output.fp_pdfs:
            if caller in f:
                outfile = f
                break
        if not outfile:
            raise ValueError(f"Could not find output file for caller {caller}")
        fig.savefig(outfile)


def plot_fn_heatmaps():
    for caller in df["caller"].unique():
        if caller in ["illumina", "longshot"]:
            continue
        # create a heatmap of the indel_len vs. homlen columns
        max_indel_len = df.query(
            "vartype=='INDEL' and mode == 'simplex' and caller==@caller"
        )["indel_len"].max()
        max_homlen = df.query(
            "vartype=='INDEL' and mode == 'simplex' and caller==@caller"
        )["homlen"].max()
        # set vmin and vmax to the min and max of the number of FP variants
        vmax = (
            df.query(
                "vartype=='INDEL' and mode == 'simplex' and decision == 'FN' and caller == @caller"
            )
            .groupby(["caller", "model", "homlen", "indel_len"])
            .count()["sample"]
            .max()
        )

        fig, axes = plt.subplots(
            figsize=(10, 10), dpi=300, nrows=2, ncols=2, sharex=True
        )

        for i, model in enumerate(["fast", "hac", "sup", "illumina"]):
            ax = axes.flatten()[i]
            if model == "illumina":
                data = df.query("caller == 'illumina'")
            else:
                data = df.query("caller == @caller")
                data = data.query("model == @model and mode == 'simplex'")
            data = data.query("vartype == 'INDEL' and decision == 'FN'")

            # data = data.query(f"homlen >= {HOM_THRESH}")
            x = "indel_len"
            y = "homlen"
            # turn data into a pivot table where the cells are the counts of the number of variants
            # in each bin
            pivot = data.pivot_table(index=x, columns=y, aggfunc="size", fill_value=0)
            total_fps = data.query("decision == 'FN'").shape[0]
            # calculate the percentage of the total number of FP variants in each bin
            # pivot = pivot / total_fps * 100
            # fill in the missing columns and rows
            for j in range(1, max_indel_len + 1):
                if j not in pivot.index:
                    pivot.loc[j] = 0
            for j in range(1, max_homlen + 1):
                if j not in pivot.columns:
                    pivot[j] = 0
            pivot.sort_index(axis=0, inplace=True)
            pivot.sort_index(axis=1, inplace=True)

            # create a heatmap of the pivot table
            mask = pivot == 0

            fmt = "d"
            cbar_ticks = [0, 1, 10, 50, 100, 500, 1000, 2000]
            cbar_kws = {
                "ticks": cbar_ticks,
                "format": mticker.FixedFormatter(cbar_ticks),
                "label": "Number of FN indels" if i in [1, 3] else "",
            }
            sns.heatmap(
                pivot,
                ax=ax,
                cmap="coolwarm",
                cbar_kws=cbar_kws,
                linecolor="black",
                linewidths=0.5,
                cbar=True,
                fmt=fmt,
                norm=SymLogNorm(vmax=vmax, linthresh=1),
            )
            ax.set_title(model)
            if i < 2:
                ax.set_xlabel("")
            else:
                ax.set_xlabel("Homopolymer length (bp)")
            if i in [1, 3]:
                ax.set_ylabel("")
            else:
                ax.set_ylabel("Indel length (bp)")

            # turn off gridlines
            ax.grid(False)

            # add a vertical dashed line at the HOM_THRESH and align to left

            ax.axvline(
                HOM_THRESH - 1,
                color="red",
                linestyle="-",
                snap=True,
                label="Homopolymer threshold",
            )
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
            # rotate xticklabels
            ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

        fig.tight_layout()
        handles, labels = ax.get_legend_handles_labels()
        for h in handles:
            h.set_linewidth(2)
        fig.legend(
            handles=handles,
            labels=labels,
            loc="center",
            fancybox=True,
            shadow=True,
            bbox_to_anchor=(0.51, 0.515),
        )
        outfile = ""
        for f in snakemake.output.fn_pdfs:
            if caller in f:
                outfile = f
                break
        if not outfile:
            raise ValueError(f"Could not find output file for caller {caller}")
        fig.savefig(outfile)


def plot_density_and_distance():
    fig = plt.figure(dpi=300, figsize=(10, 10))

    gs = fig.add_gridspec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[0:, -1])

    ax1.sharey(ax2)
    ax1.sharex(ax2)

    for caller, ax in [("illumina", ax1), ("clair3", ax2)]:
        palette = {
            "TP": named_colors["blue"],
            "FP": named_colors["orange"],
            "FN": named_colors["vermilion"],
        }
        if caller == "illumina":
            qry = f"caller == '{caller}'"
        else:
            qry = f"caller == '{caller}' and model != 'fast'"
        data = df.query(qry)
        bins = 40
        sns.histplot(
            data=data,
            x="density",
            hue="decision",
            bins=bins,
            stat="percent",
            common_norm=False,
            element="step",
            # fill=False,
            palette=palette,
            ax=ax,
            legend=True,
        )

        ax.set_title(caller)

        ax.set_xlabel("Variants in 100bp window around call")

    hue_order = sorted(set(df["caller"]))
    # move illumina to end
    hue_order.remove("illumina")
    hue_order.append("illumina")
    # map each caller in hue order to a colour in cud
    palette = {c: cud()[i] for i, c in enumerate(hue_order)}

    cap = 0.99999

    ax = ax3
    vals = []
    for b in [True, False]:
        for caller in ["illumina", "clair3"]:
            if caller == "illumina":
                qry = f"caller == '{caller}'"
            else:
                qry = f"caller == '{caller}' and  model == 'sup' and mode == 'simplex'"
            data = df.query(qry)
            threshold = int(b)
            tp = data.query(f"decision == 'TP' and dist >= {threshold}").shape[0]
            fp = data.query(f"decision == 'FP' and dist >= {threshold}").shape[0]
            fn = data.query(f"decision == 'FN' and dist >= {threshold}").shape[0]
            precision = tp / (tp + fp) if tp + fp > 0 else 1
            recall = tp / (tp + fn) if tp + fn > 0 else 1
            vals.append(
                {
                    "caller": caller,
                    "metric": "recall",
                    "value": recall,
                    "exclude_repetitive": b,
                }
            )
            vals.append(
                {
                    "caller": caller,
                    "metric": "precision",
                    "value": precision,
                    "exclude_repetitive": b,
                }
            )

    x = "mask_repeats"
    y = "F1_SCORE"
    hue = "caller"

    best_mask.loc[:, y] = best_mask[y].apply(lambda v: cap if v > cap else v)

    sns.stripplot(
        data=best_mask,
        x=x,
        y=y,
        hue=hue,
        ax=ax,
        palette=palette,
        dodge=True,
        linewidth=0.5,
        alpha=0.5,
        legend=False,
    )
    sns.boxplot(
        data=best_mask,
        x=x,
        y=y,
        hue=hue,
        ax=ax,
        palette=palette,
        showfliers=False,
        dodge=True,
        fill=False,
    )

    ax.set_ylabel("F1 score")
    ax.set_yscale("logit")
    yticks = [0.8, 0.9, 0.99, 0.995, 0.999, 0.9995, 0.9999, cap]
    yticklabels = [f"{yval:.2%}" for yval in yticks]
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    ax.set_xlabel("Repeats masked")

    fig.tight_layout()
    # add a, b, c, d labels
    for i, ax in enumerate([ax1, ax2, ax3]):
        ax.text(
            -0.05,
            1.03,
            chr(97 + i),
            transform=ax.transAxes,
            fontsize=14,
            fontweight="bold",
            va="top",
        )

    fig.savefig(snakemake.output.density_pdf)


plot_density_and_distance()
plot_fp_heatmaps()
plot_fn_heatmaps()
