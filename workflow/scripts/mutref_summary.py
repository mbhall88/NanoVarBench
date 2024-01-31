"""This script takes a list of log files from the mutref pipeline and summarizes the 
results into CSV and Latex tables
"""
import sys

sys.stderr = open(snakemake.log[0], "w")

import re
from pathlib import Path

import pandas as pd


def main():
    # Parse log files
    data = []
    for log in snakemake.input.logs:
        sample = Path(log).stem
        species = snakemake.params.samples2species[sample].replace("_", " ")
        with open(log) as f:
            contents = f.read()

        # example line: [2024-01-25 12:05:58] SUCCESS  | Downloaded 10000 genomes
        m = re.search(r"Downloaded (\d+) genomes", contents)
        if m:
            n_genomes = int(m.group(1))
        else:
            raise Exception(f"Could not find number of genomes in {log}")

        # example line: [2024-01-25 12:07:55] SUCCESS  | Using GCF_018995165.1_ASM1899516v1_genomic.fna.gz, which has a mash distance of 0.00393593 and 8530/10000 shared hashes
        pattern = r"Using (?P<acc>GC[AF]_\d+\.\d+).*\.fna\.gz, which has a mash distance of (?P<dist>\d+\.\d+) and (?P<shared_hashes>\d+)/(?P<total_hashes>\d+) shared hashes"
        m = re.search(pattern, contents)
        if m:
            acc = m.group("acc")
            dist = float(m.group("dist"))
            shared_hashes = int(m.group("shared_hashes"))
            total_hashes = int(m.group("total_hashes"))
        else:
            raise Exception(f"Could not find donor information in {log}")

        # example line: [2024-01-25 12:08:55] INFO     | # substitutions: 2217
        # but must confirm the line contains INFO so as not to confuse with minimap2 stderr
        pattern = r"\[.*\] INFO     \| # substitutions: (?P<snps>\d+)"
        m = re.search(pattern, contents)
        if m:
            snps = int(m.group("snps"))
        else:
            raise Exception(f"Could not find substitutions in {log}")

        # example lines: [2024-01-25 12:08:55] INFO     | # 1bp insertions: 63
        pattern = r"\[.*\] INFO     \| # 1bp insertions: (?P<insertions>\d+)"
        m = re.search(pattern, contents)
        if m:
            insertions_1bp = int(m.group("insertions"))
        else:
            raise Exception(f"Could not find insertions in {log}")

        # example lines: [2024-01-25 12:08:55] INFO     | # 2bp insertions: 63
        pattern = r"\[.*\] INFO     \| # 2bp insertions: (?P<insertions>\d+)"
        m = re.search(pattern, contents)
        if m:
            insertions_2bp = int(m.group("insertions"))
        else:
            raise Exception(f"Could not find insertions in {log}")

        # example lines: [2024-01-25 12:08:55] INFO     | # [3,50) insertions: 11
        pattern = r"\[.*\] INFO     \| # \[3,50\) insertions: (?P<insertions>\d+)"
        m = re.search(pattern, contents)
        if m:
            insertions_long = int(m.group("insertions"))
        else:
            raise Exception(f"Could not find insertions in {log}")

        # example lines: [2024-01-25 12:08:55] INFO     | # 1bp deletions: 63
        pattern = r"\[.*\] INFO     \| # 1bp deletions: (?P<deletions>\d+)"
        m = re.search(pattern, contents)
        if m:
            deletions_1bp = int(m.group("deletions"))
        else:
            raise Exception(f"Could not find deletions in {log}")

        # example lines: [2024-01-25 12:08:55] INFO     | # 2bp deletions: 63
        pattern = r"\[.*\] INFO     \| # 2bp deletions: (?P<deletions>\d+)"
        m = re.search(pattern, contents)
        if m:
            deletions_2bp = int(m.group("deletions"))
        else:
            raise Exception(f"Could not find deletions in {log}")

        # example lines: [2024-01-25 12:08:55] INFO     | # [3,50) deletions: 11
        pattern = r"\[.*\] INFO     \| # \[3,50\) deletions: (?P<deletions>\d+)"
        m = re.search(pattern, contents)
        if m:
            deletions_long = int(m.group("deletions"))
        else:
            raise Exception(f"Could not find deletions in {log}")

        deletions = deletions_1bp + deletions_2bp + deletions_long
        insertions = insertions_1bp + insertions_2bp + insertions_long
        n_variants = snps + insertions + deletions

        data.append(
            (
                sample,
                species,
                n_genomes,
                acc,
                dist,
                shared_hashes,
                total_hashes,
                snps,
                insertions_1bp,
                insertions_2bp,
                insertions_long,
                insertions,
                deletions_1bp,
                deletions_2bp,
                deletions_long,
                deletions,
                n_variants,
            )
        )

    # Create dataframe
    df = pd.DataFrame(
        data,
        columns=[
            "sample",
            "species",
            "n_genomes",
            "acc",
            "dist",
            "shared_hashes",
            "total_hashes",
            "snps",
            "insertions_1bp",
            "insertions_2bp",
            "insertions_long",
            "insertions",
            "deletions_1bp",
            "deletions_2bp",
            "deletions_long",
            "deletions",
            "n_variants",
        ],
    )

    # Write CSV
    df.to_csv(snakemake.output.csv, index=False)

    latex_cols = [
        "sample",
        "species",
        "dist",
        "snps",
        "insertions",
        "deletions",
        "n_variants",
    ]
    latex_df = df[latex_cols].copy()
    rename_cols = {
        "sample": "Sample",
        "species": "Species",
        "dist": "Distance",
        "snps": "SNPs",
        "insertions": "Insertions",
        "deletions": "Deletions",
        "n_variants": "Total variants",
    }
    latex_df = latex_df.rename(columns=rename_cols)
    latex_df = latex_df.sort_values("Total variants", ascending=False)

    # Write Latex
    caption = (
        "Summary of the Mash distance and number of variants found between each sample and "
        "its donor genome."
    )
    label = "tab:mutref_summary"
    position = "ht"
    latex = latex_df.to_latex(
        index=False,
        escape=True,
        na_rep="-",
        float_format="%.4f",
        caption=caption,
        label=label,
        position=position,
        formatters={
            "Species": lambda s: f"\\textit{{{s}}}",
        },
    )
    # add fullwidth environment
    lines = latex.splitlines()
    lines.insert(1, "\\begin{fullwidth}")
    lines.insert(-1, "\\end{fullwidth}")
    latex = "\n".join(lines)
    with open(snakemake.output.latex, "w") as f:
        print(latex, file=f)


if __name__ == "__main__":
    main()
