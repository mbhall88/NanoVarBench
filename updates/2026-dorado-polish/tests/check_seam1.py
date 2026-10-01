"""Check Seam 1's outputs: the aggregated results table and the fixture's known TPs.

    python3 check_seam1.py <OUTDIR> <expected_tp.tsv>
"""

import csv
import sys
from pathlib import Path

outdir, expected_path = Path(sys.argv[1]), Path(sys.argv[2])
KEYS = ("ATCC_25922_fixture", "hac", "25")


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


failures = []


def check(ok, message):
    print(("PASS " if ok else "FAIL ") + message)
    if not ok:
        failures.append(message)


# 1. One row per Arm D x hac x {SNP, INDEL, ALL} x {sweep_best, default_pass}.
results = read_tsv(outdir / "results/tables/results.tsv")
got = sorted((r["arm"], r["read_model"], r["var_type"], r["scoring_mode"]) for r in results)
want = sorted(
    ("D", "hac", t, m) for t in ("SNP", "INDEL", "ALL") for m in ("sweep_best", "default_pass")
)
check(got == want, f"results.tsv rows are Arm D x hac x variant type x scoring mode: {got}")
for r in results:
    check(
        (r["sample"], r["read_model"], r["depth"]) == KEYS,
        f"row keys {r['sample']}/{r['read_model']}/{r['depth']}",
    )

# 2. Dorado's version and the Calling model it loaded are recorded.
for r in results:
    check(r["caller"] == "dorado", f"caller recorded as {r['caller']!r}")
    check(r["caller_version"].startswith("2.1.2"), f"caller_version {r['caller_version']!r}")
    check(
        r["calling_model"].startswith("dna_r10.4.1_e8.2_400bps_polish_bacterial"),
        f"calling_model {r['calling_model']!r}",
    )
    break

# 3. The fixture's known variants are true positives (all records, no QUAL threshold).
# vcfdist's truth.tsv is 0-based and writes indels without the VCF anchor base.
def vcfdist_key(contig, pos, ref, alt):
    pos = int(pos)
    if len(ref) == len(alt):
        return (contig, str(pos - 1), ref, alt)
    return (contig, str(pos), ref[1:], alt[1:])


expected = read_tsv(expected_path)
truth = read_tsv(outdir / "work/score/ATCC_25922_fixture/hac/25x/D/sweep/truth.tsv")
status = {(r["CONTIG"], r["POS"], r["REF"], r["ALT"]): r["ERRTYPE"] for r in truth}
missed = [
    e for e in expected if status.get(vcfdist_key(e["CONTIG"], e["POS"], e["REF"], e["ALT"])) != "TP"
]
check(not missed, f"all {len(expected)} known fixture variants are TPs (missed: {missed})")

# ... and the Default-PASS rows (every PASS record, unthresholded) count them.
by_mode = {(r["var_type"], r["scoring_mode"]): r for r in results}
n_snp = sum(len(e["REF"]) == len(e["ALT"]) for e in expected)
for var_type, n in (("SNP", n_snp), ("INDEL", len(expected) - n_snp), ("ALL", len(expected))):
    tp = int(by_mode[(var_type, "default_pass")]["truth_tp"])
    check(tp >= n, f"{var_type} default_pass truth_tp {tp} >= {n} known TPs")

# 4. Actual depth is reported per contig.
depth = read_tsv(outdir / "results/tables/depth.tsv")
check(len(depth) == 1 and float(depth[0]["mean_depth"]) > 0, f"depth.tsv: {depth}")

# 5. The timing-only CPU run of Dorado is checked against the main calls, not scored.
cmp = read_tsv(outdir / "results/calls/ATCC_25922_fixture/hac/25x/D.cpu_vs_main.tsv")
check(len(cmp) == 1 and cmp[0]["same_calls"] == "yes", f"CPU Dorado calls match the main run: {cmp}")
check((outdir / "results/benchmarks/call_dorado_cpu/ATCC_25922_fixture.hac.25x.D.tsv").exists(),
      "CPU Dorado benchmark recorded")

if failures:
    sys.exit(f"{len(failures)} check(s) failed")
print("Seam 1 passed")
