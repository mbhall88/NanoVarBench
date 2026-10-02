"""Check Seam 1's AF filter analysis (#20): Clair3 ran without a haploid mode, the AF filter
resolved every HET by its AF and left the rest alone, its tables carry each threshold next to
the Arm's own --haploid_precise calls and the reference Arm, and results.tsv is unchanged. With
--reuse-dry-run, also check that a fresh work_dir reading this run's (input_work_dir)
schedules only the AF filter's own jobs.

    python3 check_af_filter.py <OUTDIR> --arms C --reference-arm D --thresholds 0.50 0.80 \
        --read-models hac sup --depths 5 10 25 50 [--reuse-dry-run dry_run.txt]
"""

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("outdir", type=Path)
ap.add_argument("--arms", nargs="+", required=True)
ap.add_argument("--reference-arm", required=True)
ap.add_argument("--thresholds", nargs="+", required=True)
ap.add_argument("--read-models", nargs="+", required=True)
ap.add_argument("--depths", nargs="+", type=int, required=True)
ap.add_argument("--reuse-dry-run", type=Path, help="snakemake -n output of the reuse run")
args = ap.parse_args()
outdir = args.outdir
SAMPLE = "ATCC_25922_fixture"
THRESHOLDS = sorted(f"{float(t):.2f}" for t in args.thresholds)
READ_SETS = [(rm, d) for rm in args.read_models for d in sorted(args.depths)]
TOP = max(args.depths)
VAR_TYPES = ("SNP", "INDEL", "ALL")
MODES = ("sweep_best", "default_pass")
failures = []


def check(ok, message):
    print(("PASS " if ok else "FAIL ") + message)
    if not ok:
        failures.append(message)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_vcf(path):
    """(header lines, records), each record a dict with its key, GT alleles and FORMAT values."""
    header, records = [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line.rstrip("\n"))
                continue
            f = line.rstrip("\n").split("\t")
            fmt = dict(zip(f[8].split(":"), f[9].split(":")))
            gt = [int(a) if a != "." else -1 for a in re.split(r"[/|]", fmt["GT"])]
            records.append({"key": (f[0], f[1], f[3], f[4]), "gt": gt, "raw_gt": fmt["GT"], "fmt": fmt})
    return header, records


def call_dir(rm, d, arm):
    return outdir / f"work/call_clair3_af/{SAMPLE}/{rm}/{d}x/{arm}"


# 1. Clair3 ran without a haploid mode, with the Arm's other options, so it wrote diploid
# genotypes; at the fixture's repeat it calls HETs.
n_het = 0
for rm, d in READ_SETS:
    for arm in args.arms:
        header, records = read_vcf(call_dir(rm, d, arm) / "variants.vcf.gz")
        cmdline = next((h for h in header if h.startswith("##cmdline=")), "")
        check(
            "--haploid" not in cmdline and "--no_phasing_for_fa" in cmdline and "--enable_long_indel" in cmdline,
            f"{rm} {d}x Arm {arm}: Clair3 ran without a haploid mode, with the Arm's other options",
        )
        check(
            records and all(len(r["gt"]) == 2 for r in records),
            f"{rm} {d}x Arm {arm}: {len(records)} records, all with diploid genotypes",
        )
        n_het += sum(len(set(r["gt"])) == 2 for r in records)
check(n_het > 0, f"Clair3 called {n_het} HETs over the Read sets")


# 2. The AF filter (af_filter): every record is kept, HETs become homozygous for the ALT with
# the highest AF when it is >= the threshold and homozygous REF otherwise, and homozygous calls
# are unchanged.
def expected_gt(r, t):
    alleles = sorted({a for a in r["gt"] if a > 0})
    if len(set(r["gt"])) < 2:
        return r["gt"]
    afs = [float(x) for x in r["fmt"]["AF"].split(",")]
    best = max(alleles, key=lambda a: (afs[a - 1], -a))
    return [best, best] if afs[best - 1] >= float(t) else [0, 0]


for rm, d in READ_SETS:
    for arm in args.arms:
        _, raw = read_vcf(call_dir(rm, d, arm) / "variants.vcf.gz")
        n_filtered = []
        for t in THRESHOLDS:
            header, out = read_vcf(call_dir(rm, d, arm) / f"af{t}.vcf.gz")
            check(
                [r["key"] for r in out] == [r["key"] for r in raw],
                f"{rm} {d}x Arm {arm} AF {t}: the AF filter keeps every record",
            )
            wrong = [
                (r["key"], r["raw_gt"], o["raw_gt"])
                for r, o in zip(raw, out)
                if o["gt"] != expected_gt(r, t)
            ]
            check(not wrong, f"{rm} {d}x Arm {arm} AF {t}: every HET resolved by AF, homs unchanged: {wrong[:3]}")
            check(any(f"min_af={float(t)}" in h for h in header), f"{rm} {d}x Arm {arm} AF {t}: header records the threshold")
            vcf = outdir / f"results/clair3_af_filter/calls/{SAMPLE}/{rm}/{d}x/{arm}.af{t}.filter.vcf.gz"
            _, filtered = read_vcf(vcf)
            check(
                all(r["raw_gt"] in ("1", "2", "3") for r in filtered),
                f"{rm} {d}x Arm {arm} AF {t}: the Filter chain left only haploid ALT calls",
            )
            n_filtered.append(len(filtered))
        # A higher threshold only turns more HETs into REF.
        check(
            n_filtered == sorted(n_filtered, reverse=True),
            f"{rm} {d}x Arm {arm}: filtered records fall as the threshold rises: {n_filtered}",
        )

# 3. clair3_af_filter.tsv: each Read set x variant type x scoring mode has the AF Arm's main
# calls, every threshold, and the reference Arm. The main rows are the scores results.tsv has.
table = read_tsv(outdir / "results/tables/clair3_af_filter.tsv")
got = sorted((r["read_model"], int(r["depth"]), r["arm"], r["calls"], r["af_threshold"], r["var_type"], r["scoring_mode"]) for r in table)
methods = [(a, "main", "") for a in args.arms] + [(a, "af_filter", t) for a in args.arms for t in THRESHOLDS]
if args.reference_arm not in args.arms:
    methods.append((args.reference_arm, "main", ""))
want = sorted((rm, d, *m, v, mode) for rm, d in READ_SETS for m in methods for v in VAR_TYPES for mode in MODES)
check(got == want, f"clair3_af_filter.tsv has {len(got)} rows, want {len(want)}")
check({r["sample"] for r in table} == {SAMPLE}, "clair3_af_filter.tsv is for the fixture Sample")

results = read_tsv(outdir / "results/tables/results.tsv")
res = {(r["read_model"], int(r["depth"]), r["arm"], r["var_type"], r["scoring_mode"]): r for r in results}
fields = ("precision", "recall", "f1", "f1_qscore", "min_qual", "truth_tp", "query_tp", "truth_fn", "query_fp")
mismatched = [
    (r["read_model"], r["depth"], r["arm"], r["var_type"], r["scoring_mode"])
    for r in table
    if r["calls"] == "main"
    and any(r[f] != res[(r["read_model"], int(r["depth"]), r["arm"], r["var_type"], r["scoring_mode"])][f] for f in fields)
]
check(not mismatched, f"the main rows are the Arms' results.tsv scores: {mismatched[:3]}")
check(
    all(
        "--haploid" not in r["clair3_options"] and bool(r["clair3_options"]) == (r["calls"] == "af_filter")
        for r in table
    ),
    "the af_filter rows record Clair3's options, without a haploid mode",
)

# results.tsv is untouched by the analysis: Arms A-D only, and none of its columns.
check({r["arm"] for r in results} == {"A", "B", "C", "D"}, "results.tsv has only Arms A-D")
check("calls" not in results[0] and "af_threshold" not in results[0], "results.tsv has no AF filter columns")

# 4. The summary: with one Sample the median is that Sample's F1, the comparisons are against
# the same Arm's main calls, and the main rows have none.
summary = read_tsv(outdir / "results/tables/clair3_af_filter_summary.tsv")
by_key = {(r["read_model"], int(r["depth"]), r["arm"], r["calls"], r["af_threshold"], r["var_type"], r["scoring_mode"]): r for r in table}
check(len(summary) == len(table), f"the summary has a row per table row (one Sample): {len(summary)}")
bad = []
for s in summary:
    key = (s["read_model"], int(s["depth"]), s["arm"], s["calls"], s["af_threshold"], s["var_type"], s["scoring_mode"])
    r = by_key[key]
    f1 = float(r["f1"])
    ok = s["n_samples"] == "1" and abs(float(s["median_f1"]) - f1) < 1e-9 and s["truth_fn"] == r["truth_fn"]
    if s["calls"] == "main":
        ok &= s["median_f1_diff_vs_main"] == "" and s["n_best"] == ""
    else:
        main = float(by_key[key[:3] + ("main", "") + key[5:]]["f1"])
        best = max(float(by_key[key[:3] + ("af_filter", t) + key[5:]]["f1"]) for t in THRESHOLDS)
        ok &= abs(float(s["median_f1_diff_vs_main"]) - (f1 - main)) < 1e-9
        ok &= int(s["n_better"]) + int(s["n_tied"]) + int(s["n_worse"]) == 1
        ok &= s["n_best"] == str(int(f1 == best))
        ok &= abs(float(s["median_loss_vs_best_threshold"]) - (best - f1)) < 1e-9
    if not ok:
        bad.append(key)
check(not bad, f"summary medians, sums and comparisons follow the per-Sample table: {bad[:3]}")

# 5. The analysis' point (#6): at the fixture's repeat Clair3 calls the true ALTs as HETs, which
# --haploid_precise drops. At the lowest threshold the AF filter keeps more true SNPs than the
# Arm's own calls at the largest Depth.
t0 = THRESHOLDS[0]
for rm in args.read_models:
    for arm in args.arms:
        main = int(by_key[(rm, TOP, arm, "main", "", "SNP", "sweep_best")]["truth_tp"])
        af = int(by_key[(rm, TOP, arm, "af_filter", t0, "SNP", "sweep_best")]["truth_tp"])
        print(f"NOTE {rm} {TOP}x Arm {arm}: SNP TPs {main} with --haploid_precise, {af} with the AF filter at {t0}")
        check(af > main, f"{rm} {TOP}x Arm {arm}: the AF filter at {t0} recovers SNPs --haploid_precise misses")

# 6. Reusing a finished run (input_work_dir): a fresh work_dir schedules only the AF filter's
# own jobs, none of the Read sets, alignments, Arms' calls or scores.
if args.reuse_dry_run:
    text = args.reuse_dry_run.read_text()
    stats = text[text.index("Job stats:"):].splitlines()
    jobs = {}
    for line in stats[3:]:
        parts = line.split()
        if len(parts) != 2 or parts[0] == "total":
            break
        jobs[parts[0]] = int(parts[1])
    allowed = {"call_clair3_af", "af_filter", "filter_calls_af", "vcfdist_af", "clair3_af_filter_tables", "clair3_af_filter"}
    n_runs = len(READ_SETS) * len(args.arms)
    check(set(jobs) <= allowed, f"the reuse run schedules only the AF filter's jobs: {jobs}")
    check(
        jobs.get("call_clair3_af") == n_runs and jobs.get("vcfdist_af") == 2 * n_runs * len(THRESHOLDS),
        f"the reuse run calls each Read set once per Arm and scores every threshold twice: {jobs}",
    )

if failures:
    sys.exit(f"{len(failures)} check(s) failed")
print("AF filter checks passed")
