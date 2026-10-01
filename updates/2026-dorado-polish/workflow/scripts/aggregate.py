"""Aggregate vcfdist summaries, actual depth and caller info into tidy tables.

results.tsv has one row per Sample x Arm x Read model x Depth x variant type x scoring
mode (the #6 contract):
  sweep_best    vcfdist's THRESHOLD == BEST row from the QUAL sweep (Best F1)
  default_pass  the THRESHOLD == NONE row from the PASS-only run (Default-PASS score)
Both use vcfdist's precision, recall, F1 and counts. F1 Q-score is computed here as
-10*log10(1 - F1) from the F1 column, which has 6 decimal places, so a perfect F1 is capped
at Q60 (1 - F1 is floored at 1e-6). vcfdist's own F1_QSCORE column writes 100 for F1 = 1.

depth.tsv has the actual per-contig depth of each Read set, and pr_curves.tsv the QUAL
sweep's precision-recall curve for each Read set and Arm: one row per QUAL threshold, keyed
like results.tsv (without scoring_mode, since a curve is the sweep) and with its column names.
"""

import csv
import math
import sys

sys.stderr = open(snakemake.log[0], "w")

VAR_TYPES = ["SNP", "INDEL", "ALL"]
MODES = {"sweep_best": ("sweep_summary", "BEST"), "default_pass": ("pass_summary", "NONE")}
F1_RESOLUTION = 1e-6  # vcfdist writes F1 to 6 decimal places


def f1_qscore(f1):
    """-10*log10(1 - F1), capped at Q60 for F1 = 1 (the F1 column's resolution)."""
    return f"{-10 * math.log10(max(1 - float(f1), F1_RESOLUTION)):.6f}"


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_coverage(path):
    """samtools coverage output -> list of per-contig dicts."""
    with open(path) as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]
    return rows


samples = {row["sample"]: row for row in read_tsv(snakemake.input.samples)}
arms = snakemake.params.arms
read_models = snakemake.params.read_models

results, depths, curves = [], [], []
seen_read_sets = set()
for c in snakemake.params.combos:
    sample, read_model, depth, arm = c["sample"], c["read_model"], int(c["depth"]), c["arm"]
    keys = {"sample": sample, "read_model": read_model, "depth": depth, "arm": arm}

    coverage = read_coverage(c["coverage"])
    total_len = sum(int(r["endpos"]) - int(r["startpos"]) + 1 for r in coverage)
    actual = (
        sum((int(r["endpos"]) - int(r["startpos"]) + 1) * float(r["meandepth"]) for r in coverage)
        / total_len
    )
    read_set = (sample, read_model, depth)
    if read_set not in seen_read_sets:  # depth is a Read set property, shared by its Arms
        seen_read_sets.add(read_set)
        for r in coverage:
            depths.append(
                {
                    "sample": sample,
                    "read_model": read_model,
                    "depth": depth,
                    "contig": r["rname"],
                    "length": int(r["endpos"]) - int(r["startpos"]) + 1,
                    "num_reads": r["numreads"],
                    "breadth_pct": r["coverage"],
                    "mean_depth": r["meandepth"],
                }
            )

    info = read_tsv(c["caller_info"])[0]
    arm_cfg = arms[arm]
    sample_row = samples[sample]
    common = {
        "sample": sample,
        "species": sample_row["species"],
        "dnd_sample": sample_row["dnd"],
        "read_model": read_model,
        "basecall_model": read_models[read_model],
        "depth": depth,
        "actual_depth": f"{actual:.2f}",
        "actual_depth_by_contig": ";".join(
            f"{r['rname']}={float(r['meandepth']):.2f}" for r in coverage
        ),
        "arm": arm,
        "aligner": arm_cfg["aligner"],
        "aligner_version": arm_cfg["aligner_version"],
        "preset": arm_cfg["preset"],
        "caller": info["caller"],
        "caller_version": info["caller_version"],
        "calling_model": info["calling_model"],
        # Dorado records its weights file's SHA256; Clair3 lists each model file's, as
        # name=sha256;... Container digests are recorded for tools that run in one.
        "calling_model_sha256": info.get("calling_model_sha256")
        or info.get("calling_model_weights_sha256", ""),
        "container": info.get("container", ""),
        "hardware": info["hardware"],
    }

    for mode, (summary_key, threshold) in MODES.items():
        rows = [r for r in read_tsv(c[summary_key]) if r["THRESHOLD"] == threshold]
        by_type = {r["VAR_TYPE"]: r for r in rows}
        for var_type in VAR_TYPES:
            r = by_type[var_type]
            results.append(
                {
                    **common,
                    "scoring_mode": mode,
                    "var_type": var_type,
                    "min_qual": r["MIN_QUAL"],
                    "precision": r["PREC"],
                    "recall": r["RECALL"],
                    "f1": r["F1_SCORE"],
                    "f1_qscore": f1_qscore(r["F1_SCORE"]),
                    "truth_tp": r["TRUTH_TP"],
                    "query_tp": r["QUERY_TP"],
                    "truth_fn": r["TRUTH_FN"],
                    "query_fp": r["QUERY_FP"],
                }
            )

    for r in read_tsv(c["sweep_pr"]):
        if r["VAR_TYPE"] in VAR_TYPES:
            curves.append(
                {
                    **keys,
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


def write(path, rows):
    with open(path, "w", newline="") as fh:
        if not rows:
            return
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


write(snakemake.output.results, results)
write(snakemake.output.depth, depths)
write(snakemake.output.pr_curves, curves)
print(f"results={len(results)} depth={len(depths)} pr_curves={len(curves)}", file=sys.stderr)
