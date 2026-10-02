"""Tables for the AF filter analysis (#20). Nothing here goes into results.tsv.

clair3_af_filter.tsv has one row per Sample x Read model x Depth x Arm x calls x AF threshold
x variant type x scoring mode:
  calls = main       the Arm's own calls, scored as in results.tsv (Clair3 --haploid_precise,
                     or the reference Arm's caller); af_threshold is empty
  calls = af_filter  Clair3 without a haploid mode, HETs resolved at af_threshold
The scoring modes are those of results.tsv: sweep_best (Best F1, the QUAL sweep's BEST row)
and default_pass (the Default-PASS score, the PASS-only run's NONE row).

clair3_af_filter_summary.tsv summarises it over the Samples, per Arm x calls x AF threshold x
Read model x Depth x variant type x scoring mode:
  median_f1, min_f1, max_f1         over the Samples
  truth_fn, query_fp                summed over the Samples
  median_f1_diff_vs_main            median per-Sample F1 minus the same Arm's main F1
  n_better, n_tied, n_worse         Samples whose F1 is above, equal to or below the Arm's main
  n_best                            Samples where this AF threshold has the highest F1 of all
                                    the thresholds (ties count for each)
  median_loss_vs_best_threshold     median per-Sample F1 lost by using this AF threshold
                                    rather than that Sample's best one
  gap_to_reference_closed           (median_f1 - main median_f1) / (reference Arm median_f1 -
                                    main median_f1), when the reference Arm's median is higher
The comparison columns are empty on the main rows.
"""

import csv
import math
import statistics
import sys
from collections import defaultdict

sys.stderr = open(snakemake.log[0], "w")

VAR_TYPES = ["SNP", "INDEL", "ALL"]
MODES = {"sweep_best": ("sweep_summary", "BEST"), "default_pass": ("pass_summary", "NONE")}
F1_RESOLUTION = 1e-6  # vcfdist writes F1 to 6 decimal places


def f1_qscore(f1):
    """-10*log10(1 - F1), capped at Q60 for F1 = 1 (as in results.tsv)."""
    return f"{-10 * math.log10(max(1 - float(f1), F1_RESOLUTION)):.6f}"


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


samples = {row["sample"]: row for row in read_tsv(snakemake.input.samples)}
reference_arm = snakemake.params.reference_arm

rows = []
for s in snakemake.params.score_sets:
    sample_row = samples[s["sample"]]
    common = {
        "sample": s["sample"],
        "species": sample_row["species"],
        "dnd_sample": sample_row["dnd"],
        "read_model": s["read_model"],
        "depth": int(s["depth"]),
        "arm": s["arm"],
        "calls": s["calls"],
        "af_threshold": s["af_threshold"],
        "clair3_options": snakemake.params.clair3_options if s["calls"] == "af_filter" else "",
    }
    for mode, (summary_key, threshold) in MODES.items():
        by_type = {
            r["VAR_TYPE"]: r for r in read_tsv(s[summary_key]) if r["THRESHOLD"] == threshold
        }
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

# --- Summary over the Samples ---------------------------------------------------------
# F1 per Sample, keyed by everything else.
f1 = {}
groups = defaultdict(list)
for r in rows:
    group = (r["arm"], r["calls"], r["af_threshold"], r["read_model"], r["depth"], r["var_type"], r["scoring_mode"])
    groups[group].append(r)
    f1[group + (r["sample"],)] = float(r["f1"])
thresholds = sorted({r["af_threshold"] for r in rows if r["calls"] == "af_filter"})


def fmt(x):
    return "" if x is None else f"{x:.6f}"


summary = []
for group, members in groups.items():
    arm, calls, af, rm, depth, var_type, mode = group
    values = [float(r["f1"]) for r in members]
    median = statistics.median(values)
    out = {
        "arm": arm,
        "calls": calls,
        "af_threshold": af,
        "read_model": rm,
        "depth": depth,
        "var_type": var_type,
        "scoring_mode": mode,
        "n_samples": len(members),
        "median_f1": fmt(median),
        "min_f1": fmt(min(values)),
        "max_f1": fmt(max(values)),
        "truth_fn": sum(int(r["truth_fn"]) for r in members),
        "query_fp": sum(int(r["query_fp"]) for r in members),
        "median_f1_diff_vs_main": "",
        "n_better": "",
        "n_tied": "",
        "n_worse": "",
        "n_best": "",
        "median_loss_vs_best_threshold": "",
        "gap_to_reference_closed": "",
    }
    if calls == "af_filter":
        main_key = (arm, "main", "", rm, depth, var_type, mode)
        diffs, losses, n_best = [], [], 0
        for r in members:
            s = r["sample"]
            diffs.append(float(r["f1"]) - f1[main_key + (s,)])
            best = max(f1[(arm, "af_filter", t, rm, depth, var_type, mode, s)] for t in thresholds)
            losses.append(best - float(r["f1"]))
            n_best += float(r["f1"]) == best
        main_median = statistics.median(f1[main_key + (r["sample"],)] for r in members)
        out.update(
            median_f1_diff_vs_main=fmt(statistics.median(diffs)),
            n_better=sum(d > 0 for d in diffs),
            n_tied=sum(d == 0 for d in diffs),
            n_worse=sum(d < 0 for d in diffs),
            n_best=n_best,
            median_loss_vs_best_threshold=fmt(statistics.median(losses)),
        )
        ref_key = (reference_arm, "main", "", rm, depth, var_type, mode)
        if reference_arm and ref_key in groups:
            ref_median = statistics.median(float(r["f1"]) for r in groups[ref_key])
            if ref_median > main_median:
                out["gap_to_reference_closed"] = fmt((median - main_median) / (ref_median - main_median))
    summary.append(out)

order = {"main": 0, "af_filter": 1}
summary.sort(
    key=lambda r: (r["read_model"], r["var_type"], r["scoring_mode"], r["depth"], r["arm"], order[r["calls"]], r["af_threshold"])
)


def write(path, table):
    with open(path, "w", newline="") as fh:
        if not table:
            return
        writer = csv.DictWriter(fh, fieldnames=list(table[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)


write(snakemake.output.tsv, rows)
write(snakemake.output.summary, summary)
print(f"rows={len(rows)} summary={len(summary)}", file=sys.stderr)
