"""PR curves of the AF filter series Figure 2 draws (#28; figures.af_filter).

The QUAL sweep's precision-recall curve of Clair3 with the AF filter at one AF threshold, for
each Read set: one row per QUAL threshold and variant type, with pr_curves.tsv's columns
(keyed by sample, Read model, Depth and Arm, the Arm being the Clair3 Arm whose alignment,
Calling model and options the calls used) and af_threshold. Nothing here goes into
results.tsv or pr_curves.tsv.
"""

import csv
import math
import sys

sys.stderr = open(snakemake.log[0], "w")

VAR_TYPES = ["SNP", "INDEL", "ALL"]
F1_RESOLUTION = 1e-6  # vcfdist writes F1 to 6 decimal places


def f1_qscore(f1):
    """-10*log10(1 - F1), capped at Q60 for F1 = 1 (as in pr_curves.tsv)."""
    return f"{-10 * math.log10(max(1 - float(f1), F1_RESOLUTION)):.6f}"


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


rows = []
for read_set, path in zip(snakemake.params.read_sets, snakemake.input.files):
    for r in read_tsv(path):
        if r["VAR_TYPE"] in VAR_TYPES:
            rows.append(
                {
                    "sample": read_set["sample"],
                    "read_model": read_set["read_model"],
                    "depth": read_set["depth"],
                    "arm": snakemake.params.arm,
                    "af_threshold": snakemake.params.af_threshold,
                    "var_type": r["VAR_TYPE"],
                    "min_qual": r["MIN_QUAL"],
                    "precision": r["PREC"],
                    "recall": r["RECALL"],
                    "f1": r["F1_SCORE"],
                    "f1_qscore": f1_qscore(r["F1_SCORE"]),
                    "truth_total": r["TRUTH_TOTAL"],
                    "truth_tp": r["TRUTH_TP"],
                    "truth_fn": r["TRUTH_FN"],
                    "query_total": r["QUERY_TOTAL"],
                    "query_tp": r["QUERY_TP"],
                    "query_fp": r["QUERY_FP"],
                }
            )

with open(snakemake.output.tsv, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
print(f"rows={len(rows)}", file=sys.stderr)
