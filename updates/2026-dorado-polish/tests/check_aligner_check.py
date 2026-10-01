"""Check Seam 1's dorado aligner check (#13): the report exists, covers SNP, INDEL and ALL in
both scoring modes, and the two paths agree closely on the fixture.

    python3 check_aligner_check.py <OUTDIR>

The fixture is tiny (63 truth variants), so one record is a big F1 step: this checks the
check's outputs and that dorado aligner's calls are close to the reheadered minimap2 path's,
not that they are identical.
"""

import csv
import sys
from pathlib import Path

outdir = Path(sys.argv[1])
out = outdir / "results/dorado_aligner_check/ATCC_25922_fixture/hac/25x"
failures = []


def check(ok, message):
    print(("PASS " if ok else "FAIL ") + message)
    if not ok:
        failures.append(message)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


report = out / "D.dorado_aligner_vs_minimap2.md"
check(report.exists() and "Calls differ materially" in report.read_text(), f"{report.name} written")
rows = read_tsv(out / "D.dorado_aligner_vs_minimap2.tsv")

# One F1 row per score and variant type.
f1 = {(r["metric"], r["var_type"]): r for r in rows if r["section"] == "vcfdist" and r["metric"].endswith("_f1")}
want = {(m, t) for m in ("best_f1", "default_pass_f1") for t in ("SNP", "INDEL", "ALL")}
check(set(f1) == want, f"F1 rows for Best and Default-PASS x SNP, INDEL, ALL: {sorted(f1)}")
for key, r in f1.items():
    a, b = float(r["path1"]), float(r["path2"])
    check(0 < a <= 1 and 0 < b <= 1 and abs(a - b) <= 0.05, f"{key} F1 path 1 {a}, path 2 {b}")

# Path 2 is not an Arm: it adds nothing to results.tsv.
results = read_tsv(outdir / "results/tables/results.tsv")
check({r["arm"] for r in results} == {"A", "B", "C", "D"}, "results.tsv has only Arms A-D")

# The two paths' calls are close, and path 2's polish took no --any-bam.
metrics = {(r["section"], r["metric"]): r["value"] for r in rows if r["section"] != "vcfdist"}
n1, n2 = int(metrics[("filtered_vcf", "path1_records")]), int(metrics[("filtered_vcf", "path2_records")])
check(n1 > 0 and n2 > 0, f"both paths made calls: {n1} and {n2} filtered records")
n_differ = int(metrics[("filtered_vcf", "records_differ")])
check(n_differ <= 0.05 * max(n1, n2), f"{n_differ} of {max(n1, n2)} filtered records differ")
log = outdir / "work/logs/call_dorado_aligner_bam/ATCC_25922_fixture.hac.25x.D.log"
check(log.exists(), "dorado polish ran on the dorado aligner BAM")
info = read_tsv(out / "D.dorado_aligner.caller_info.tsv")[0]
info1 = read_tsv(outdir / "results/calls/ATCC_25922_fixture/hac/25x/D.caller_info.tsv")[0]
check(info["calling_model"] == info1["calling_model"], "both paths resolved the same Calling model")

if failures:
    sys.exit(f"{len(failures)} check(s) failed")
print("dorado aligner check passed")
