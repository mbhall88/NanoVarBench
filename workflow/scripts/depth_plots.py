import sys

sys.stderr = open(snakemake.log[0], "w")

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
import seaborn as sns

DUPLEX_DP_CAP = snakemake.params.duplex_depth_cap
NO_INDELS = snakemake.params.no_indels
CAP = 0.99999
bbox = (0.52, 0.72)

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


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2.0, n - 1)
    return m, m - h, m + h


def plot(
    vartype: str,
    df: pd.DataFrame,
    hue: str,
    col: str,
    illumina_df: pd.DataFrame,
    illumina: dict,
    x: str,
    y: str,
    yticks: list[str],
    yticklabels: list[str],
    hue_order: list[str],
):
    data = df.query(f"VAR_TYPE == '{vartype}' and caller != 'longshot'")
    palette = cud(n=len(data[hue].unique()))
    # palette = "colorblind"
    row_order = ["F1 Score", "Recall", "Precision"]
    ncols = len(data.query("mode == 'simplex'")[col].unique()) + len(
        data.query("mode == 'duplex'")[col].unique()
    )
    nrows = len(row_order)
    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols, figsize=(10, 10), dpi=300, sharey=True, sharex=False
    )
    ill_colour = "red"
    legend = True
    i = 0
    for metric in row_order:
        ill_mean, ill_ci_low, ill_ci_high = illumina[vartype][metric]
        for mode in ["simplex", "duplex"]:
            for model in data[col].unique():
                if model == "fast" and mode == "duplex":
                    continue
                d = data.query(
                    f"{col} == '{model}' and metric == '{metric}' and mode == '{mode}'"
                )
                # transform precision so logit shows values of 1.0 - i.e. change values of 1.0 to 0.99999
                cap = CAP
                d.loc[:, y] = d[y].apply(lambda v: cap if v == 1.0 else v)
                ax = axes.flatten()[i]
                sns.pointplot(
                    ax=ax,
                    data=d,
                    x=x,
                    y=y,
                    hue=hue,
                    palette=palette,
                    units="sample",
                    dodge=True,
                    legend=legend,
                    hue_order=hue_order,
                    alpha=0.6,
                    err_kws={"linewidth": 1.5, "alpha": 0.6},
                )

                # add horizontal line for Illumina with fill_between for Std
                ax.axhline(ill_mean, color=ill_colour, linestyle="--", alpha=0.5)
                # get the x-axis limits
                xlim = ax.get_xlim()
                ax.fill_between(
                    x=xlim,
                    y1=ill_ci_low,
                    y2=ill_ci_high,
                    color=ill_colour,
                    alpha=0.1,
                )

                if i < ncols:
                    ax.set_title(f"{model} {mode}")

                if i < ncols * 2:
                    # remove x-axis labels and ticks
                    ax.set(xlabel="", xticklabels=[])
                else:
                    ax.set_xlabel("Depth")

                ax.set_ylabel(metric)

                ax.set_yscale("logit", nonpositive="clip")

                ax.set_yticks(yticks)
                ax.set_yticklabels(yticklabels)
                if legend:
                    handles, labels = ax.get_legend_handles_labels()
                    # add Illumina to legend
                    handles.append(
                        plt.Line2D([0], [0], color=ill_colour, linestyle="--")
                    )
                    labels.append("Illumina")
                    legend = False
                    ax.legend().remove()

                # plot the number of sample at each depth
                ax2 = ax.twinx()
                xs = []
                ys = []
                for depth in d[x].unique():
                    count = d.query(f"{x} == {depth}")["sample"].nunique()
                    xs.append(depth)
                    ys.append(count)

                sns.barplot(x=xs, y=ys, ax=ax2, color="black", alpha=0.1)
                ax2.grid(False)
                if (i + 1) % ncols == 0:
                    ax2.set_ylabel("Number of samples")
                else:
                    ax2.set(yticks=[], ylabel="", yticklabels=[])

                # get rid of all minor y-axis ticks
                ax.yaxis.set_minor_locator(plt.NullLocator())
                i += 1

    plt.tight_layout()
    leg_cols = math.ceil((len(data[hue].unique()) + 1) / 2)
    fig.legend(
        handles=handles,
        labels=labels,
        loc="upper center",
        bbox_to_anchor=bbox,
        ncol=leg_cols,
        title="",
        framealpha=1.0,
        fancybox=True,
        shadow=True,
    )
    return fig


def main():

    # read in the QC files to get full depth of each sample
    frames = []
    for p in map(Path, snakemake.input.postfilter_stats):
        df = pd.read_csv(p)
        df["model"] = df["model"].apply(lambda x: x.split("_")[-1].split("@")[0])
        # rename filename column to sample and mean_coverage to depth
        df = df.rename(columns={"filename": "sample", "mean_coverage": "depth"})
        df["mode"] = p.parts[-3]
        frames.append(df)

    depth_df = pd.concat(frames)
    depth_df.set_index(["sample", "mode", "model"], inplace=True, verify_integrity=True)

    illumina_frames = []
    frames = []
    for p in map(Path, snakemake.input.ont_pr):
        df = pd.read_csv(p, sep="\t")
        sample = p.parent.name
        df["sample"] = sample

        df["caller"] = p.parts[-7]
        dp = int(p.parts[-6][:-1])
        mode = p.parts[-5]
        model = p.parts[-3].split("_")[-1].split("@")[0]
        full_dp = depth_df.loc[(sample, mode, model), "depth"]
        if full_dp < (dp * 0.9):
            continue
        if mode == "duplex" and model == "fast":
            continue
        if mode == "duplex" and dp > DUPLEX_DP_CAP:
            continue
        df["depth"] = dp
        df["mode"] = mode
        df["version"] = p.parts[-4]
        df["model"] = model
        frames.append(df)

    for p in map(Path, snakemake.input.illumina_pr):
        df = pd.read_csv(p, sep="\t")
        sample = p.parent.name
        df["sample"] = sample

        df["caller"] = "illumina"
        illumina_frames.append(df)

    illumina_df = pd.concat(illumina_frames)
    illumina_df.reset_index(inplace=True, drop=True)

    pr_df = pd.concat(frames)
    pr_df.reset_index(inplace=True, drop=True)

    # drop all columns except the following
    keep_cols = [
        "MIN_QUAL",
        "PREC",
        "RECALL",
        "F1_SCORE",
        "VAR_TYPE",
        "sample",
        "caller",
        "depth",
        "mode",
        "version",
        "model",
    ]
    pr_df = pr_df[keep_cols]

    x = "depth"
    y = "value"
    hue = "caller"
    col = "model"

    best_ix = pr_df.groupby([x, hue, col, "sample", "VAR_TYPE", "mode"])[
        "F1_SCORE"
    ].idxmax()
    best_df = pr_df.iloc[best_ix]
    best_df = best_df.melt(
        id_vars=keep_cols[4:],
        value_vars=["PREC", "RECALL", "F1_SCORE"],
        var_name="metric",
        value_name="value",
    )
    # rename the values to be more descriptive
    best_df["metric"] = best_df["metric"].replace(
        {"PREC": "Precision", "RECALL": "Recall", "F1_SCORE": "F1 Score"}
    )

    # gather the best Illumina F1 scores for each sample and then
    # take the mean and std of F1, Recall, and Precision
    ix = illumina_df.groupby(["sample", "VAR_TYPE"])["F1_SCORE"].idxmax()
    illumina_best_df = illumina_df.iloc[ix]
    illumina_best_df = illumina_best_df.melt(
        id_vars=["sample", "VAR_TYPE"],
        value_vars=["PREC", "RECALL", "F1_SCORE"],
        var_name="metric",
        value_name="value",
    )
    illumina_best_df["metric"] = illumina_best_df["metric"].replace(
        {"PREC": "Precision", "RECALL": "Recall", "F1_SCORE": "F1 Score"}
    )
    illumina_snp = illumina_best_df.query("VAR_TYPE == 'SNP'")
    illumina_snp_f1 = illumina_snp.query("metric == 'F1 Score'")["value"]
    illumina_snp_f1_mean, illumina_snp_f1_low, illumina_snp_f1_high = (
        mean_confidence_interval(illumina_snp_f1)
    )
    illumina_snp_recall = illumina_snp.query("metric == 'Recall'")["value"]
    illumina_snp_recall_mean, illumina_snp_recall_low, illumina_snp_recall_high = (
        mean_confidence_interval(illumina_snp_recall)
    )
    illumina_snp_precision = illumina_snp.query("metric == 'Precision'")["value"]
    (
        illumina_snp_precision_mean,
        illumina_snp_precision_low,
        illumina_snp_precision_high,
    ) = mean_confidence_interval(illumina_snp_precision)
    illumina_indel = illumina_best_df.query("VAR_TYPE == 'INDEL'")
    illumina_indel_f1 = illumina_indel.query("metric == 'F1 Score'")["value"]
    illumina_indel_f1_mean, illumina_indel_f1_low, illumina_indel_f1_high = (
        mean_confidence_interval(illumina_indel_f1)
    )
    illumina_indel_recall = illumina_indel.query("metric == 'Recall'")["value"]
    (
        illumina_indel_recall_mean,
        illumina_indel_recall_low,
        illumina_indel_recall_high,
    ) = mean_confidence_interval(illumina_indel_recall)
    illumina_indel_precision = illumina_indel.query("metric == 'Precision'")["value"]
    (
        illumina_indel_precision_mean,
        illumina_indel_precision_low,
        illumina_indel_precision_high,
    ) = mean_confidence_interval(illumina_indel_precision)
    illumina = {
        "SNP": {
            "F1 Score": (
                illumina_snp_f1_mean,
                illumina_snp_f1_low,
                illumina_snp_f1_high,
            ),
            "Recall": (
                illumina_snp_recall_mean,
                illumina_snp_recall_low,
                illumina_snp_recall_high,
            ),
            "Precision": (
                illumina_snp_precision_mean,
                illumina_snp_precision_low,
                illumina_snp_precision_high,
            ),
        },
        "INDEL": {
            "F1 Score": (
                illumina_indel_f1_mean,
                illumina_indel_f1_low,
                illumina_indel_f1_high,
            ),
            "Recall": (
                illumina_indel_recall_mean,
                illumina_indel_recall_low,
                illumina_indel_recall_high,
            ),
            "Precision": (
                illumina_indel_precision_mean,
                illumina_indel_precision_low,
                illumina_indel_precision_high,
            ),
        },
    }

    hue_order = sorted(pr_df["caller"].unique())

    cap = CAP
    yticks = [0.01, 0.1, 0.25, 0.5, 0.9, 0.99, 0.999, 0.9999, cap]
    yticklabels = [f"{yval:.2%}" for yval in yticks]

    vartype = "SNP"
    snp_fig = plot(
        vartype=vartype,
        df=best_df,
        hue=hue,
        col=col,
        illumina_df=illumina_df,
        illumina=illumina,
        x=x,
        y=y,
        yticks=yticks,
        yticklabels=yticklabels,
        hue_order=hue_order,
    )
    snp_fig.savefig(snakemake.output.snp_fig)

    # move no indel callers to end
    for caller in NO_INDELS:
        hue_order.remove(caller)
        hue_order.append(caller)

    vartype = "INDEL"
    indel_fig = plot(
        vartype=vartype,
        df=best_df.query("caller not in @NO_INDELS"),
        hue=hue,
        col=col,
        illumina_df=illumina_df,
        illumina=illumina,
        x=x,
        y=y,
        yticks=yticks,
        yticklabels=yticklabels,
        hue_order=hue_order,
    )


if __name__ == "__main__":
    main()
