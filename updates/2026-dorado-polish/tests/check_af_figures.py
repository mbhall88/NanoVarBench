"""Check Seam 1's AF filter series in the figures and Table S1 (#28), with the analysis on (the
run's own results/) and off (af_off/results/, rendered from the same aggregated tables).

On: Table S1 has one more row per Read set, labelled as the AF filter variant and an extra
analysis, with the scores clair3_af_filter.tsv has for that Arm and threshold; the Arms' rows
are the same as with the analysis off; the PR curves table is that Arm's QUAL sweep at the
threshold; and each figure names the series. Off: no AF filter rows, nothing about it in the
figures, and the Arms' rows are as in the run.

    python3 check_af_figures.py <OUTDIR> --arm C --threshold 0.65 \
        --read-models hac sup --depths 5 10 25 50
"""

import argparse
import csv
import sys
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("outdir", type=Path)
ap.add_argument("--arm", required=True)
ap.add_argument("--threshold", required=True)
ap.add_argument("--read-models", nargs="+", required=True)
ap.add_argument("--depths", nargs="+", type=int, required=True)
args = ap.parse_args()
THRESHOLD = f"{float(args.threshold):.2f}"
READ_SETS = {(rm, str(d)) for rm in args.read_models for d in args.depths}
ON, OFF = args.outdir / "results", args.outdir / "af_off/results"
FIGURES = ("fig1_best_f1_depth", "fig2_pr_curves", "fig3_per_sample_best_f1")
failures = []


def check(ok, message):
    print(("PASS " if ok else "FAIL ") + message)
    if not ok:
        failures.append(message)


def read_csv(path, delimiter=","):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def is_af(row):
    return "AF filter" in row["Arm"]


on = read_csv(ON / "tables/table_s1_per_sample.csv")
off = read_csv(OFF / "tables/table_s1_per_sample.csv")

# Table S1, analysis on: one AF filter row per Read set, labelled, with the AF table's scores.
af_rows = [r for r in on if is_af(r)]
check(
    {(r["Read model"], r["Depth (x)"]) for r in af_rows} == READ_SETS and len(af_rows) == len(READ_SETS),
    f"Table S1 has one AF filter row per Read set ({len(af_rows)} rows)",
)
label = af_rows[0]["Arm"] if af_rows else ""
check(
    f"AF filter ({THRESHOLD})" in label and "extra analysis" in label and label.startswith(args.arm),
    f"the AF filter rows are labelled as an extra analysis with the threshold: {label!r}",
)
table = read_csv(ON / "tables/clair3_af_filter.tsv", "\t")
want = {
    (r["read_model"], r["depth"], r["var_type"], r["scoring_mode"]): r
    for r in table
    if r["calls"] == "af_filter" and r["arm"] == args.arm and r["af_threshold"] == THRESHOLD
}
bad = []
for r in af_rows:
    for var_type in ("SNP", "INDEL"):
        for mode, name in (("sweep_best", "Best"), ("default_pass", "Default-PASS")):
            src = want[(r["Read model"], r["Depth (x)"], var_type, mode)]
            for column, metric in (("precision", "precision"), ("recall", "recall"), ("f1", "F1")):
                if r[f"{var_type} {name} {metric}"] != src[column]:
                    bad.append((r["Read model"], r["Depth (x)"], var_type, name, metric))
check(not bad, f"the AF filter rows have the AF filter table's scores at {THRESHOLD}: {bad[:3]}")
main = [r for r in on if not is_af(r)]
check(
    all(r["Actual depth (x)"] for r in af_rows)
    and all(
        r["Actual depth (x)"] == next(m for m in main if (m["Read model"], m["Depth (x)"]) == (r["Read model"], r["Depth (x)"]))["Actual depth (x)"]
        for r in af_rows
    ),
    "the AF filter rows have their Read set's actual depth",
)

# Table S1, analysis off: the same Arms' rows and none for the AF filter.
check(not any(is_af(r) for r in off), "Table S1 with the analysis off has no AF filter rows")
check(main == off, f"the Arms' rows are the same with the analysis on and off ({len(off)} rows)")

# The PR curves table: the Arm's QUAL sweep at the threshold, for every Read set.
curves = read_csv(ON / "tables/clair3_af_filter_pr_curves.tsv", "\t")
check(
    {r["arm"] for r in curves} == {args.arm} and {r["af_threshold"] for r in curves} == {THRESHOLD},
    f"the PR curves table is Arm {args.arm} at {THRESHOLD}",
)
check(
    {(r["read_model"], r["depth"]) for r in curves} == READ_SETS
    and {r["var_type"] for r in curves} == {"SNP", "INDEL", "ALL"},
    "the PR curves table has every Read set and variant type",
)
# Each Best F1 row is a point on the curve at its QUAL threshold, as for the Arms' curves.
point = {(r["read_model"], r["depth"], r["var_type"], r["min_qual"]): r["f1"] for r in curves}
check(
    all(point.get((rm, d, v, s["min_qual"])) == s["f1"] for (rm, d, v, m), s in want.items() if m == "sweep_best"),
    "every AF filter Best F1 is a point on its PR curve",
)

# The figures: the series is named in each figure when on, and nowhere when off.
for name in FIGURES:
    on_svg = (ON / f"figures/{name}.svg").read_text()
    off_svg = (OFF / f"figures/{name}.svg").read_text()
    check(
        f"{args.arm} + AF filter ({THRESHOLD}), extra analysis" in on_svg,
        f"{name} names the AF filter series, as an extra analysis",
    )
    check("AF filter" not in off_svg, f"{name} with the analysis off doesn't mention it")
    check((OFF / f"figures/{name}.png").stat().st_size > 0, f"{name}.png is rendered with the analysis off")

if failures:
    sys.exit(f"{len(failures)} check(s) failed")
print("AF figure checks passed")
