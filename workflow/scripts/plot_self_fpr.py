import sys

sys.stderr = open(snakemake.log[0], "w")

import math

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.cbook import boxplot_stats

SAMPLE2SPECIES = snakemake.params.sample2species

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


def plot_fpr(df, y, title=""):
    row = "depth"
    nrows = len(df[row].unique())
    col = "mode"
    ncols = len(df[col].unique())
    hue = "model"
    x = "caller"
    palette = cud(n=len(df[hue].unique()))
    annotated_species = set()
    fig, axes = plt.subplots(
        ncols=ncols,
        nrows=nrows,
        figsize=(10, 10),
        dpi=300,
        squeeze=False,
        sharey=True,
        tight_layout=True,
    )
    legend_drawn = False
    for i, depth in enumerate(df[row].unique()):
        depth_df = df.query(f"{row} == '{depth}'")
        for j, mode in enumerate(depth_df[col].unique()):
            mode_df = depth_df.query(f"{col} == '{mode}'")
            ax = axes[i, j]
            sns.boxplot(
                x=x,
                y=y,
                data=mode_df,
                ax=ax,
                hue=hue,
                fill=False,
                fliersize=0,
                palette=palette,
                gap=0.2,
            )
            strip = sns.stripplot(
                x=x,
                y=y,
                data=mode_df,
                ax=ax,
                hue=hue,
                dodge=True,
                alpha=0.75,
                legend=False,
                linewidth=0.5,
                edgecolor="black",
                palette=palette,
            )

            for collection_i, collection in enumerate(strip.collections):
                caller_i = math.floor(collection_i / len(df[hue].unique()))
                caller = mode_df[x].unique()[caller_i]
                model_i = collection_i % len(df[hue].unique())
                offsets = collection.get_offsets()
                xs = []
                ys = []
                for xval, yval in offsets:
                    xs.append(xval)
                    ys.append(yval)

                stats = boxplot_stats(ys)[0]
                whishi = stats["whishi"]
                outliers = [y for y in stats["fliers"] if y > whishi]
                for ix, yval in enumerate(ys):
                    if yval in outliers:
                        # annotate with the name of sample that is an outlier
                        caller_df = mode_df.query(f"{x} == '{caller}'")
                        # get the index of the y value closest to the outlier
                        s_ix = (caller_df[y] - yval).abs().idxmin()
                        sample = caller_df["sample"].loc[s_ix]
                        species = SAMPLE2SPECIES[sample]
                        annotated_species.add(species)
                        s = species.split("_")[0][0] + species.split("_")[1][0]
                        ax.annotate(
                            s,
                            (xs[ix], yval),
                            xycoords="data",
                            textcoords="offset points",
                            xytext=(3, 50),
                            ha="left",
                            va="bottom",
                            fontsize=6,
                            arrowprops=dict(
                                arrowstyle="->",
                                lw=1,
                                color=cud()[model_i],
                                relpos=(0, 1),
                            ),
                        )

            ax.set_title(f"{mode}")
            if "FPR" in y:
                ax.set_ylabel("FPR")
                ax.set_yscale("symlog", linthresh=0.0000001)
                yticks = [
                    0,
                    0.00000001,
                    0.0000001,
                    0.000001,
                    0.00001,
                    0.0001,
                    0.001,
                    0.01,
                ]
            elif "FP/Mbp" in y:
                ax.set_ylabel("FP/Mbp")
                ax.set_yscale("symlog", linthresh=1)
                yticks = [0, 1, 10, 100, 1000, 10000]
            else:
                ax.set_ylabel(y)
            ax.set_xlabel("Caller")

            ax.set_yticks(yticks)
            ax.set_yticklabels(yticks)
            # rotate x labels
            for tick in ax.get_xticklabels():
                tick.set_rotation(90)

            if legend_drawn:
                ax.get_legend().remove()
            else:
                ax.legend(loc="upper right", title=hue)
                legend_drawn = True

    # add a legend for the species initials
    if annotated_species:
        annotated_species = sorted(list(annotated_species))
        abbrevs = [s.split("_")[0][0] + s.split("_")[1][0] for s in annotated_species]
        handles = [
            plt.Line2D(
                [0], [0], marker="o", color="w", markerfacecolor="w", markersize=0
            )
        ] * len(annotated_species)
        labels = list(
            [
                f"{a} = $\it{{{s.split('_')[0]}}}$ $\it{{{s.split('_')[1]}}}$"
                for a, s in zip(abbrevs, annotated_species)
            ]
        )
        axes[0, 1].legend(
            handles,
            labels,
            loc="best",
            # title="Species",
            fontsize=6,
        )

    fig.suptitle(title, fontsize=16)

    return fig, axes


def main():
    df = pd.read_csv(snakemake.input.csv)
    df.sort_values(by=["sample", "caller"], inplace=True)
    # shorten model column from dna_r10.4.1_e8.2_400bps_fast@v4.3.0 to fast
    df["model"] = df["model"].apply(lambda x: x.split("_")[-1].split("@")[0])
    df["seqlen"] = df["FP"] + df["TN"]
    df["SNP_FPR"] = df["SNP"] / df["seqlen"]
    df["INDEL_FPR"] = df["INDEL"] / df["seqlen"]
    df["SNP_FP/Mbp"] = df["SNP"] / df["seqlen"] * 1e6
    df["INDEL_FP/Mbp"] = df["INDEL"] / df["seqlen"] * 1e6
    df["FP/Mbp"] = df["FP"] / df["seqlen"] * 1e6

    fpr_rate = df.groupby(
        ["version", "mode", "model", "depth", "caller", "sample"]
    ).apply(lambda a: a[:])
    fpr_rate.to_csv(snakemake.output.csv, index=False)

    # SNP FPR plot
    snp_fig, snp_axes = plot_fpr(df, "SNP_FP/Mbp", "SNP False positives per Mbp")
    snp_fig.savefig(snakemake.output.snps_png)

    # INDEL FPR plot
    no_indels = snakemake.params.no_indels
    data = df.query(f"caller not in {no_indels}")
    indel_fig, indel_axes = plot_fpr(
        data, "INDEL_FP/Mbp", "Indel False positives per Mbp"
    )
    indel_fig.savefig(snakemake.output.indels_png)


main()
