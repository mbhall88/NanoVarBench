"""Tables for the smallvar pilot (#29). Nothing here goes into results.tsv or benchmarks.tsv.

smallvar_pilot.tsv has one row per Sample x Read model x Depth x calls x Arm x variant type x
scoring mode:
  calls = smallvar_haploid  dorado smallvar with the configured model forced in
                            (--model-override) and the whole genome hemizygous, on the Arm's
                            alignment
  calls = main              the Arm's own calls, scored as in results.tsv
  calls = af_filter         the AF filter (#20) on the Arm's diploid Clair3 calls, at
                            af_threshold
basecall_model_mismatch says how the calling model's basecall model differs from the Read
model's: version (hac v4.3.0 reads, hac v6.0.0 smallvar model), tier_and_version (sup v4.3.0
reads) or none (the Arms' own Calling models, which match their reads). The scoring modes are
those of results.tsv: sweep_best (Best F1, the QUAL sweep's BEST row) and default_pass (the
Default-PASS score, the PASS-only run's NONE row).

smallvar_pilot_runtime.tsv has one row per smallvar job: Snakemake's benchmark of the job,
which only runs `dorado smallvar`. hardware and threads are config values, as in
benchmarks.tsv.
"""

import csv
import math
import sys

sys.stderr = open(snakemake.log[0], "w")

VAR_TYPES = ["SNP", "INDEL", "ALL"]
MODES = {"sweep_best": ("sweep_summary", "BEST"), "default_pass": ("pass_summary", "NONE")}
F1_RESOLUTION = 1e-6  # vcfdist writes F1 to 6 decimal places
# Snakemake benchmark column -> our column, as in benchmarks.tsv
BENCHMARK_COLUMNS = {
    "s": "wall_time_s",
    "max_rss": "max_rss_mb",
    "max_vms": "max_vms_mb",
    "cpu_time": "cpu_time_s",
    "mean_load": "mean_cpu_pct",
}


def f1_qscore(f1):
    """-10*log10(1 - F1), capped at Q60 for F1 = 1 (as in results.tsv)."""
    return f"{-10 * math.log10(max(1 - float(f1), F1_RESOLUTION)):.6f}"


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(path, table):
    with open(path, "w", newline="") as fh:
        if not table:
            return
        writer = csv.DictWriter(fh, fieldnames=list(table[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)


samples = {row["sample"]: row for row in read_tsv(snakemake.input.samples)}
read_basecall_models = snakemake.params.read_basecall_models

rows = []
for s in snakemake.params.score_sets:
    sample_row = samples[s["sample"]]
    common = {
        "sample": s["sample"],
        "species": sample_row["species"],
        "dnd_sample": sample_row["dnd"],
        "read_model": s["read_model"],
        "basecall_model": read_basecall_models[s["read_model"]],
        "depth": int(s["depth"]),
        "calls": s["calls"],
        "arm": s["arm"],
        "af_threshold": s["af_threshold"],
        "calling_model": s["calling_model"],
        "basecall_model_mismatch": s["basecall_model_mismatch"],
    }
    for mode, (summary_key, threshold) in MODES.items():
        by_type = {r["VAR_TYPE"]: r for r in read_tsv(s[summary_key]) if r["THRESHOLD"] == threshold}
        for var_type in VAR_TYPES:
            r = by_type[var_type]
            rows.append(
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

runtime = []
for run in snakemake.params.runs:
    bench = read_tsv(run["benchmark"])[0]
    runtime.append(
        {
            "sample": run["sample"],
            "read_model": run["read_model"],
            "depth": int(run["depth"]),
            "arm": run["arm"],
            "tool": "dorado smallvar",
            "tool_version": snakemake.params.dorado_version[run["arm"]],
            "calling_model": snakemake.params.model,
            "basecall_model_mismatch": snakemake.params.mismatch[run["read_model"]],
            "device": snakemake.params.device,
            "hardware": snakemake.params.hardware,
            "threads": snakemake.params.threads,
            **{ours: bench[theirs] for theirs, ours in BENCHMARK_COLUMNS.items()},
        }
    )

write(snakemake.output.tsv, rows)
write(snakemake.output.runtime, runtime)
print(f"rows={len(rows)} runtime={len(runtime)}", file=sys.stderr)
