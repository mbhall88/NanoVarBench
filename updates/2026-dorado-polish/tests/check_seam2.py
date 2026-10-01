"""Check Seam 2's outputs: the aggregated tables built from tests/seam2/'s hand-written inputs.

    python3 check_seam2.py <OUTDIR>

Only the output tables are read. Every expected value is a literal from the hand-written
inputs (see tests/seam2.sh for which Read set and Arm gets which).
"""

import csv
import sys
from pathlib import Path

outdir = Path(sys.argv[1])
tables = outdir / "results/tables"
WORKED = "ATCC_25922__202309"
ARMS = ("A", "B", "C", "D")
DEPTHS = ("25", "50")
VAR_TYPES = ("SNP", "INDEL", "ALL")
MODES = ("sweep_best", "default_pass")


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


failures = []


def check(ok, message):
    print(("PASS " if ok else "FAIL ") + message)
    if not ok:
        failures.append(message)


SAMPLES = [r["sample"] for r in read_tsv(Path(__file__).parents[1] / "config/samples.tsv")]
check(len(SAMPLES) == 14, f"config/samples.tsv lists the 14 Samples: {SAMPLES}")

results = read_tsv(tables / "results.tsv")
RESULT_KEYS = ("sample", "read_model", "depth", "arm", "var_type", "scoring_mode")
by_key = {tuple(r[k] for k in RESULT_KEYS): r for r in results}
check(len(by_key) == len(results), f"results.tsv has one row per key ({len(results)} rows)")


def row(arm, var_type, mode, sample=WORKED, depth="50"):
    return by_key[(sample, "hac", depth, arm, var_type, mode)]


def fields(r, *names):
    return tuple(r[n] for n in names)


# 1. One row per Sample x Read model x Depth x Arm x variant type x scoring mode, with
# sweep-best and default-PASS rows for every Arm, and no SV rows.
want = {
    (s, "hac", d, a, t, m)
    for s in SAMPLES
    for d in DEPTHS
    for a in ARMS
    for t in VAR_TYPES
    for m in MODES
}
check(set(by_key) == want, f"results.tsv keys: {len(by_key)} rows, want {len(want)}")

# 2. sweep_best is vcfdist's THRESHOLD == BEST row from the QUAL sweep (Best F1), not the
# unthresholded NONE row: in the Arm D worked example BEST is at QUAL>=12 (SNP, ALL) and
# QUAL>=7 (INDEL).
COUNTS = ("min_qual", "truth_tp", "query_tp", "truth_fn", "query_fp", "precision", "recall", "f1")
for var_type, want_row in {
    "SNP": ("12", "990", "990", "10", "10", "0.990000", "0.990000", "0.990000"),
    "INDEL": ("7", "100", "100", "0", "0", "1.000000", "1.000000", "1.000000"),
    "ALL": ("12", "1090", "1090", "10", "10", "0.990909", "0.990909", "0.990909"),
}.items():
    got = fields(row("D", var_type, "sweep_best"), *COUNTS)
    check(got == want_row, f"Arm D {var_type} sweep_best is the sweep's BEST row: {got}")

# 3. default_pass is the PASS-only run's unthresholded NONE row (the Default-PASS score), not
# its BEST row (QUAL>=15 for SNP and ALL) and not anything from the QUAL sweep.
for var_type, want_row in {
    "SNP": ("0", "985", "985", "15", "2", "0.997974", "0.985000", "0.991444"),
    "INDEL": ("0", "98", "98", "2", "0", "1.000000", "0.980000", "0.989899"),
    "ALL": ("0", "1083", "1083", "17", "2", "0.998157", "0.984545", "0.991304"),
}.items():
    got = fields(row("D", var_type, "default_pass"), *COUNTS)
    check(got == want_row, f"Arm D {var_type} default_pass is the PASS run's NONE row: {got}")

# 4. Every Arm keeps its own sweep and default-PASS rows. Clair3 sets FILTER too, so Arm A's
# Default-PASS score drops its LowQual calls (worked_A) and differs from its Best F1.
for var_type, sweep, pas in (
    ("SNP", ("0", "960", "40", "0.979592"), ("0", "955", "45", "0.976982")),
    ("INDEL", ("0", "98", "2", "0.989899"), ("0", "97", "3", "0.984772")),
    ("ALL", ("0", "1058", "42", "0.980538"), ("0", "1052", "48", "0.977695")),
):
    for mode, want_row in (("sweep_best", sweep), ("default_pass", pas)):
        got = fields(row("A", var_type, mode), "min_qual", "truth_tp", "truth_fn", "f1")
        check(got == want_row, f"Arm A {var_type} {mode}: {got}")
# ... and every other Read set and Arm got the filler's rows.
for (s, _, d, a, t, m), r in by_key.items():
    if (s, d) == (WORKED, "50") and a in ("A", "D"):
        continue
    want_tp = {
        "sweep_best": {"SNP": "900", "INDEL": "80", "ALL": "980"},
        "default_pass": {"SNP": "890", "INDEL": "78", "ALL": "968"},
    }[m][t]
    if r["truth_tp"] != want_tp:
        check(False, f"{s} {d}x Arm {a} {t} {m} truth_tp {r['truth_tp']}, want {want_tp}")

# 5. F1 Q-score is -10*log10(1 - F1), from the F1 column. F1 has 6 decimal places, so a
# perfect F1 (1 - F1 = 0) is capped at Q60, the column's resolution (vcfdist writes 100).
QSCORE_TOL = 1e-4
for arm, var_type, mode, want_q in (
    ("D", "SNP", "sweep_best", 20.0),  # F1 0.990000
    ("D", "INDEL", "sweep_best", 60.0),  # F1 1.000000: capped
    ("D", "ALL", "sweep_best", 20.413883),  # F1 0.990909
    ("D", "SNP", "default_pass", 20.677292),  # F1 0.991444
    ("A", "INDEL", "sweep_best", 19.956356),  # F1 0.989899
):
    got = row(arm, var_type, mode)["f1_qscore"]
    check(
        abs(float(got) - want_q) < QSCORE_TOL,
        f"Arm {arm} {var_type} {mode} f1_qscore {got} == {want_q}",
    )

# 6. Actual depth is joined from each Read set's own coverage: the worked Read set has a
# 4,000 bp chromosome at 50x and a 1,000 bp plasmid at 75x, so its length-weighted mean is
# 55x. Every other Read set has one 5,000 bp contig at 26x. Every Arm of a Read set shares it.
for (s, _, d, a, t, m), r in by_key.items():
    worked = (s, d) == (WORKED, "50")
    want_depth = (
        ("55.00", "chromosome=50.00;plasmid=75.00") if worked else ("26.00", "chromosome=26.00")
    )
    got = fields(r, "actual_depth", "actual_depth_by_contig")
    if got != want_depth or r["depth"] != d:
        check(False, f"{s} {d}x Arm {a} {t} {m} actual depth {got}, want {want_depth}")
check(
    fields(row("B", "SNP", "sweep_best"), "depth", "actual_depth") == ("50", "55.00"),
    "the worked Read set's Depth (the cap) sits next to its actual depth, 55.00",
)
check(
    fields(row("D", "SNP", "sweep_best", depth="25"), "depth", "actual_depth") == ("25", "26.00"),
    "ATCC_25922 hac 25x has its own actual depth, 26.00",
)
# depth.tsv has one row per contig per Read set (not per Arm).
depth = read_tsv(tables / "depth.tsv")
got = sorted((r["sample"], r["depth"], r["contig"], r["length"], r["mean_depth"]) for r in depth)
want_depth_rows = sorted(
    [(WORKED, "50", "chromosome", "4000", "50.0"), (WORKED, "50", "plasmid", "1000", "75.0")]
    + [
        (s, d, "chromosome", "5000", "26.0")
        for s in SAMPLES
        for d in DEPTHS
        if (s, d) != (WORKED, "50")
    ]
)
check(got == want_depth_rows, f"depth.tsv has each Read set's contigs once ({len(depth)} rows)")

# 7. The dnd flag: S. enterica (ATCC_10708) and V. parahaemolyticus (ATCC_17802) are the dnd
# samples (dorado#1599); the other 12 Samples aren't.
DND = {"ATCC_10708__202309": "Salmonella enterica", "ATCC_17802__202309": "Vibrio parahaemolyticus"}
for s in SAMPLES:
    rows = [r for k, r in by_key.items() if k[0] == s]
    flags = {r["dnd_sample"] for r in rows}
    want_flag = "true" if s in DND else "false"
    check(
        flags == {want_flag},
        f"{s} dnd_sample {sorted(flags)} == {want_flag} on all {len(rows)} rows",
    )
    if s in DND:
        check({r["species"] for r in rows} == {DND[s]}, f"{s} is {DND[s]}")

# 8. pr_curves.tsv: the QUAL sweep's precision and recall at each threshold, keyed like
# results.tsv and with its column names, so figure code can join the two.
curves = read_tsv(tables / "pr_curves.tsv")
CURVE_COLUMNS = [
    "sample", "read_model", "depth", "arm", "var_type", "min_qual",
    "precision", "recall", "f1", "f1_qscore",
    "truth_total", "truth_tp", "truth_fn", "query_total", "query_tp", "query_fp",
]
got = list(curves[0]) if curves else []
check(got == CURVE_COLUMNS, f"pr_curves.tsv columns: {got}")
check(
    {fields(r, *RESULT_KEYS[:4]) for r in curves} == {k[:4] for k in by_key},
    "pr_curves.tsv has a curve for every Sample x Read model x Depth x Arm in results.tsv",
)
check(
    {r["var_type"] for r in curves} == set(VAR_TYPES), "pr_curves.tsv has SNP, INDEL and ALL only"
)
CURVE_KEYS = (*RESULT_KEYS[:4], "var_type", "min_qual")
curve = {fields(r, *CURVE_KEYS): r for r in curves}
check(len(curve) == len(curves), f"pr_curves.tsv has one row per key and threshold ({len(curves)})")
for var_type in VAR_TYPES:
    got = [
        r["min_qual"]
        for r in curves
        if fields(r, "sample", "depth", "arm", "var_type") == (WORKED, "50", "D", var_type)
    ]
    check(got == ["0", "7", "12", "60"], f"Arm D {var_type} curve thresholds {got}")
for var_type, min_qual, want_q in (
    ("INDEL", "7", 60.0),  # F1 1.000000: capped
    ("SNP", "60", 4.771217),  # F1 0.666667
    ("ALL", "7", 18.672602),  # F1 0.986425
):
    got = curve.get((WORKED, "hac", "50", "D", var_type, min_qual), {}).get("f1_qscore", "")
    check(
        got != "" and abs(float(got) - want_q) < QSCORE_TOL,
        f"Arm D {var_type} curve at QUAL>={min_qual} f1_qscore {got!r} == {want_q}",
    )
# Each Best F1 row is the point on its curve at the best threshold.
JOINED = ("precision", "recall", "f1", "f1_qscore", "truth_tp", "truth_fn", "query_tp", "query_fp")
NO_POINT = dict.fromkeys(JOINED)
mismatched = [
    k
    for k, r in by_key.items()
    if k[5] == "sweep_best"
    and fields(curve.get((*k[:5], r["min_qual"]), NO_POINT), *JOINED) != fields(r, *JOINED)
]
check(not mismatched, f"each sweep_best row is its curve's point at min_qual: {mismatched}")

if failures:
    sys.exit(f"{len(failures)} check(s) failed")
print("Seam 2 passed")
