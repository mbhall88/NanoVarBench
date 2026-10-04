"""Check Seam 1's smallvar pilot (#29): dorado smallvar ran with the model forced in and the
whole genome hemizygous, wrote haploid GT:GQ calls (no AF or allele depths, so the AF filter
can't apply), its table carries every smallvar row's Basecall-model mismatch next to the
compared Arms and the AF filter, and results.tsv and benchmarks.tsv are unchanged. With
--reuse-dry-run, also check that a fresh work_dir reading this run's (input_work_dir)
schedules only the pilot's own jobs.

    python3 check_smallvar_pilot.py <OUTDIR> --arms D --compare-arms C D --af-arm C \
        --af-threshold 0.50 --read-models hac sup --depths 5 10 25 50 \
        [--reuse-dry-run dry_run.txt]
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("outdir", type=Path)
ap.add_argument("--arms", nargs="+", required=True)
ap.add_argument("--compare-arms", nargs="+", required=True)
ap.add_argument("--af-arm", required=True)
ap.add_argument("--af-threshold", required=True)
ap.add_argument("--read-models", nargs="+", required=True)
ap.add_argument("--depths", nargs="+", type=int, required=True)
ap.add_argument("--reuse-dry-run", type=Path, help="snakemake -n output of the reuse run")
args = ap.parse_args()
outdir = args.outdir
SAMPLE = "ATCC_25922_fixture"
MODEL = "dna_r10.4.1_e8.2_400bps_hac@v6.0.0_smallvar@v1.0"
MISMATCH = {"hac": "version", "sup": "tier_and_version"}  # reads are v4.3.0, the model hac v6.0.0
AF_T = f"{float(args.af_threshold):.2f}"
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
    """(header lines, records as split fields)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    header, records = [], []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                header.append(line)
            else:
                records.append(line.split("\t"))
    return header, records


# 1. smallvar's own VCF: the model was forced in (Dorado warns it doesn't match the reads), the
# calls are haploid, and FORMAT has only GT and GQ, with INFO empty: no AF or allele depths.
for rm, d in READ_SETS:
    for arm in args.arms:
        call = outdir / f"work/call_dorado_smallvar/{SAMPLE}/{rm}/{d}x/{arm}"
        log = (outdir / f"work/logs/call_dorado_smallvar/{SAMPLE}.{rm}.{d}x.{arm}.log").read_text()
        check(
            f"{MODEL}/config.toml" in log and "not compatible with the input BAM" in log,
            f"{rm} {d}x Arm {arm}: smallvar used the overridden model and warned of the mismatch",
        )
        header, records = read_vcf(call / "variants.vcf")
        fmt_ids = {h.split("ID=")[1].split(",")[0] for h in header if h.startswith("##FORMAT=")}
        check(fmt_ids == {"GT", "GQ"}, f"{rm} {d}x Arm {arm}: FORMAT declares only GT and GQ: {sorted(fmt_ids)}")
        check(
            all(r[8] == "GT:GQ" and r[7] == "." for r in records),
            f"{rm} {d}x Arm {arm}: every record is GT:GQ with no INFO (no AF, AD or DP)",
        )
        check(
            records and all("/" not in r[9] and "|" not in r[9] for r in records),
            f"{rm} {d}x Arm {arm}: {len(records)} records, all haploid (whole genome hemizygous)",
        )
        _, filtered = read_vcf(outdir / f"results/smallvar_pilot/calls/{SAMPLE}/{rm}/{d}x/{arm}.smallvar.filter.vcf.gz")
        check(
            filtered and all(r[9].split(":")[0] in ("1", "2", "3") for r in filtered),
            f"{rm} {d}x Arm {arm}: the Filter chain left {len(filtered)} haploid ALT calls",
        )

# 2. smallvar_pilot.tsv: each Read set x variant type x scoring mode has smallvar's calls, each
# compared Arm's own calls and the AF filter at its threshold.
table = read_tsv(outdir / "results/tables/smallvar_pilot.tsv")
got = sorted((r["read_model"], int(r["depth"]), r["calls"], r["arm"], r["af_threshold"], r["var_type"], r["scoring_mode"]) for r in table)
methods = (
    [("smallvar_haploid", a, "") for a in args.arms]
    + [("main", a, "") for a in args.compare_arms]
    + [("af_filter", args.af_arm, AF_T)]
)
want = sorted((rm, d, *m, v, mode) for rm, d in READ_SETS for m in methods for v in VAR_TYPES for mode in MODES)
check(got == want, f"smallvar_pilot.tsv has {len(got)} rows, want {len(want)}")
check({r["sample"] for r in table} == {SAMPLE}, "smallvar_pilot.tsv is for the fixture Sample")

# Every smallvar row is labelled with its Basecall-model mismatch; the Arms' rows aren't.
bad = [
    (r["read_model"], r["depth"], r["calls"])
    for r in table
    if (r["calls"] == "smallvar_haploid")
    != (r["basecall_model_mismatch"] == MISMATCH[r["read_model"]] and r["calling_model"] == MODEL)
    or (r["calls"] != "smallvar_haploid" and r["basecall_model_mismatch"] != "none")
]
check(not bad, f"smallvar rows carry the model and their mismatch (hac: version, sup: tier_and_version): {bad[:3]}")

fields = ("precision", "recall", "f1", "f1_qscore", "min_qual", "truth_tp", "query_tp", "truth_fn", "query_fp")
results = read_tsv(outdir / "results/tables/results.tsv")
res = {(r["read_model"], int(r["depth"]), r["arm"], r["var_type"], r["scoring_mode"]): r for r in results}
mismatched = [
    (r["read_model"], r["depth"], r["arm"], r["var_type"], r["scoring_mode"])
    for r in table
    if r["calls"] == "main"
    and any(r[f] != res[(r["read_model"], int(r["depth"]), r["arm"], r["var_type"], r["scoring_mode"])][f] for f in fields)
]
check(not mismatched, f"the main rows are the Arms' results.tsv scores: {mismatched[:3]}")

af = read_tsv(outdir / "results/tables/clair3_af_filter.tsv")
af_rows = {
    (r["read_model"], int(r["depth"]), r["var_type"], r["scoring_mode"]): r
    for r in af
    if r["calls"] == "af_filter" and r["arm"] == args.af_arm and r["af_threshold"] == AF_T
}
mismatched = [
    (r["read_model"], r["depth"], r["var_type"], r["scoring_mode"])
    for r in table
    if r["calls"] == "af_filter"
    and any(r[f] != af_rows[(r["read_model"], int(r["depth"]), r["var_type"], r["scoring_mode"])][f] for f in fields)
]
check(not mismatched, f"the af_filter rows are clair3_af_filter.tsv's at {AF_T}: {mismatched[:3]}")

# results.tsv and benchmarks.tsv are untouched by the pilot.
check({r["arm"] for r in results} == {"A", "B", "C", "D"}, "results.tsv has only Arms A-D")
check("basecall_model_mismatch" not in results[0], "results.tsv has no pilot columns")
bench = read_tsv(outdir / "results/tables/benchmarks.tsv")
check(not any("smallvar" in r["tool"] for r in bench), "benchmarks.tsv has no smallvar rows")

# 3. The runtime table: one row per smallvar job, labelled, with the job's benchmark.
runtime = read_tsv(outdir / "results/tables/smallvar_pilot_runtime.tsv")
check(
    sorted((r["read_model"], int(r["depth"]), r["arm"]) for r in runtime)
    == sorted((rm, d, a) for rm, d in READ_SETS for a in args.arms),
    f"smallvar_pilot_runtime.tsv has a row per Read set and Arm: {len(runtime)}",
)
check(
    all(
        float(r["wall_time_s"]) > 0 and float(r["max_rss_mb"]) > 0 and r["tool"] == "dorado smallvar"
        and r["calling_model"] == MODEL and r["basecall_model_mismatch"] == MISMATCH[r["read_model"]]
        for r in runtime
    ),
    "every runtime row has a wall time, peak RSS, the model and its mismatch",
)

# 4. At the largest Depth smallvar finds the fixture's SNPs, despite the mismatch. How it compares
# with the Arms on real data is the pilot's question, not the seam's.
by_key = {(r["read_model"], int(r["depth"]), r["calls"], r["arm"], r["var_type"], r["scoring_mode"]): r for r in table}
for rm in args.read_models:
    for arm in args.arms:
        sv = int(by_key[(rm, TOP, "smallvar_haploid", arm, "SNP", "sweep_best")]["truth_tp"])
        ref = int(by_key[(rm, TOP, "main", arm, "SNP", "sweep_best")]["truth_tp"]) if arm in args.compare_arms else None
        print(f"NOTE {rm} {TOP}x: SNP TPs {sv} with smallvar, {ref} with Arm {arm}")
        check(sv > 0, f"{rm} {TOP}x Arm {arm}: smallvar calls true SNPs")

# 5. Reusing a finished run (input_work_dir): a fresh work_dir schedules only the pilot's own
# jobs, none of the Read sets, alignments, Arms' calls or scores.
if args.reuse_dry_run:
    text = args.reuse_dry_run.read_text()
    stats = text[text.index("Job stats:"):].splitlines()
    jobs = {}
    for line in stats[3:]:
        parts = line.split()
        if len(parts) != 2 or parts[0] == "total":
            break
        jobs[parts[0]] = int(parts[1])
    allowed = {"call_dorado_smallvar", "filter_calls_smallvar", "vcfdist_smallvar", "smallvar_pilot_tables", "smallvar_pilot"}
    n_runs = len(READ_SETS) * len(args.arms)
    check(set(jobs) <= allowed, f"the reuse run schedules only the pilot's jobs: {jobs}")
    check(
        jobs.get("call_dorado_smallvar") == n_runs and jobs.get("vcfdist_smallvar") == 2 * n_runs,
        f"the reuse run calls each Read set once per Arm and scores it twice: {jobs}",
    )

if failures:
    sys.exit(f"{len(failures)} check(s) failed")
print("smallvar pilot checks passed")
