"""Check Seam 1's outputs: the aggregated results table and the fixture's known TPs.

    python3 check_seam1.py <OUTDIR> <expected_tp.tsv>
"""

import csv
import re
import sys
from pathlib import Path

outdir, expected_path = Path(sys.argv[1]), Path(sys.argv[2])
known_fn_path = expected_path.with_name("clair3_known_fn.tsv")
KEYS = ("ATCC_25922_fixture", "hac", "25")
ARMS = ("A", "B", "C", "D")
# Arm -> (aligner version, preset, caller, caller version prefix). Arms B, C and D share one
# alignment; Arm A has its own.
ARM_SPEC = {
    "A": ("2.26", "map-ont", "clair3", "1.0.5"),
    "B": ("2.31", "lr:hq", "clair3", "1.0.5"),
    "C": ("2.31", "lr:hq", "clair3", "2.0.3"),
    "D": ("2.31", "lr:hq", "dorado", "2.1.2"),
}


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


failures = []


def check(ok, message):
    print(("PASS " if ok else "FAIL ") + message)
    if not ok:
        failures.append(message)


# 1. One row per Arm x hac x {SNP, INDEL, ALL} x {sweep_best, default_pass}, for Arms A-D.
results = read_tsv(outdir / "results/tables/results.tsv")
got = sorted((r["arm"], r["read_model"], r["var_type"], r["scoring_mode"]) for r in results)
want = sorted(
    (a, "hac", t, m)
    for a in ARMS
    for t in ("SNP", "INDEL", "ALL")
    for m in ("sweep_best", "default_pass")
)
check(got == want, f"results.tsv rows are Arm x hac x variant type x scoring mode: {got}")
for r in results:
    check(
        (r["sample"], r["read_model"], r["depth"]) == KEYS,
        f"row keys {r['sample']}/{r['read_model']}/{r['depth']}",
    )

# 2. Each Arm's aligner, caller and Calling model are recorded, with a checksum for the
# Calling model's files and the digest of the container each tool ran in.
SHA256 = re.compile(r"[0-9a-f]{64}")
by_arm = {a: [r for r in results if r["arm"] == a] for a in ARMS}
for arm, (aln_version, preset, caller, caller_version) in ARM_SPEC.items():
    for r in by_arm[arm]:
        check(
            (r["aligner"], r["aligner_version"], r["preset"]) == ("minimap2", aln_version, preset),
            f"Arm {arm} aligner {r['aligner']} {r['aligner_version']} {r['preset']}",
        )
        check(r["caller"] == caller, f"Arm {arm} caller {r['caller']!r}")
        check(
            r["caller_version"].startswith(caller_version),
            f"Arm {arm} caller_version {r['caller_version']!r}",
        )
        check(
            bool(SHA256.search(r["calling_model_sha256"])),
            f"Arm {arm} calling_model_sha256 {r['calling_model_sha256']!r}",
        )
        break
    r = by_arm[arm][0]
    if caller == "clair3":
        # the Calling model follows the Read model: hac reads -> the hac v4.3.0 model
        check(
            r["calling_model"] == "r1041_e82_400bps_hac_v430",
            f"Arm {arm} calling_model {r['calling_model']!r}",
        )
        check(
            "@sha256:" in r["container"],
            f"Arm {arm} container is pinned by digest: {r['container']!r}",
        )
    else:
        check(
            r["calling_model"].startswith("dna_r10.4.1_e8.2_400bps_polish_bacterial"),
            f"Arm {arm} calling_model {r['calling_model']!r}",
        )

# Arms A and B run Clair3 1.0.5 in one container, and Arm C runs 2.0.3 in another, so C's
# Calling-model checksums differ from the TF models of A and B.
sha = {a: by_arm[a][0]["calling_model_sha256"] for a in ARMS}
check(sha["A"] == sha["B"], "Arms A and B use the same bundled TF Calling model")
check(sha["C"] != sha["B"], "Arm C uses the HKU PyTorch Calling model, not B's TF one")
digest = {a: by_arm[a][0]["container"] for a in ("A", "B", "C")}
check(digest["A"] == digest["B"], "Arms A and B ran in the same Clair3 1.0.5 container")
check(digest["C"] != digest["B"], "Arm C ran in the Clair3 2.0.3 container")

# 3. The fixture's known variants are true positives (all records, no QUAL threshold).
# vcfdist's truth.tsv is 0-based and writes indels without the VCF anchor base.
def vcfdist_key(contig, pos, ref, alt):
    pos = int(pos)
    if len(ref) == len(alt):
        return (contig, str(pos - 1), ref, alt)
    return (contig, str(pos), ref[1:], alt[1:])


expected = read_tsv(expected_path)
# Clair3 can't recover all of them on this fixture. Eight of the 62 sit in a dense cluster
# (24 truth variants in 3.3 kb, 20.0-23.3 kb) where 3-4 of ~25 reads keep the reference allele.
# Clair3's full-alignment model (1.0.5 and 2.0.3 alike) calls seven of the eight as
# low-confidence hets (e.g. 22049: 4 ref / 19 alt reads) and `--haploid_precise` drops hets,
# so they never reach the VCF. The eighth, 23283, is called, but its neighbours 23277 and
# 23279 aren't, so vcfdist can't match the cluster. All three Clair3 Arms miss exactly these
# eight with the paper's options, so it isn't a workflow fault. The Clair3 Arms must still get
# every other known variant, and Dorado must get all of them. A Clair3 Arm that starts
# catching some of the eight passes (with a NOTE); one that misses any other variant fails.
known_fn = {(e["CONTIG"], e["POS"], e["REF"], e["ALT"]) for e in read_tsv(known_fn_path)}
assert known_fn <= {(e["CONTIG"], e["POS"], e["REF"], e["ALT"]) for e in expected}


def is_snp(e):
    return len(e["REF"]) == len(e["ALT"])


by_mode = {(r["arm"], r["var_type"], r["scoring_mode"]): r for r in results}
for arm in ARMS:
    truth = read_tsv(outdir / f"work/score/ATCC_25922_fixture/hac/25x/{arm}/sweep/truth.tsv")
    status = {(r["CONTIG"], r["POS"], r["REF"], r["ALT"]): r["ERRTYPE"] for r in truth}
    allowed_fn = known_fn if ARM_SPEC[arm][2] == "clair3" else set()
    required = [e for e in expected if (e["CONTIG"], e["POS"], e["REF"], e["ALT"]) not in allowed_fn]
    missed = [
        e
        for e in required
        if status.get(vcfdist_key(e["CONTIG"], e["POS"], e["REF"], e["ALT"])) != "TP"
    ]
    if allowed_fn:
        n_fn = sum(
            status.get(vcfdist_key(*key)) != "TP" for key in allowed_fn
        )
        print(f"NOTE Arm {arm} misses {n_fn} of the {len(allowed_fn)} known Clair3 FNs")
    check(
        not missed,
        f"Arm {arm}: {len(required)} of the {len(expected)} known fixture variants are TPs "
        f"(missed: {missed})",
    )

    # ... and the Default-PASS rows (every PASS record, unthresholded) count them.
    for var_type, pick in (("SNP", is_snp), ("INDEL", lambda e: not is_snp(e)), ("ALL", bool)):
        n = sum(1 for e in required if pick(e))
        tp = int(by_mode[(arm, var_type, "default_pass")]["truth_tp"])
        check(tp >= n, f"Arm {arm} {var_type} default_pass truth_tp {tp} >= {n} required TPs")

# 4. Actual depth is reported per contig.
depth = read_tsv(outdir / "results/tables/depth.tsv")
check(len(depth) == 1 and float(depth[0]["mean_depth"]) > 0, f"depth.tsv: {depth}")

# 5. The timing-only CPU run of Dorado is checked against the main calls, not scored. Clair3
# has no CPU re-run: it only runs on CPU.
cmp = read_tsv(outdir / "results/calls/ATCC_25922_fixture/hac/25x/D.cpu_vs_main.tsv")
check(len(cmp) == 1 and cmp[0]["same_calls"] == "yes", f"CPU Dorado calls match the main run: {cmp}")
check((outdir / "results/benchmarks/call_dorado_cpu/ATCC_25922_fixture.hac.25x.D.tsv").exists(),
      "CPU Dorado benchmark recorded")
calls = outdir / "results/calls/ATCC_25922_fixture/hac/25x"
for arm in "ABC":
    check(not (calls / f"{arm}.cpu_vs_main.tsv").exists(), f"Arm {arm} has no CPU re-run")
    check(not (outdir / f"results/benchmarks/call_dorado_cpu/ATCC_25922_fixture.hac.25x.{arm}.tsv").exists(),
          f"Arm {arm} has no Dorado CPU benchmark")

# 6. Clair3's runtime is benchmarked for Arms A-C, and the Filter chain output exists for all.
for arm in "ABC":
    check((outdir / f"results/benchmarks/call_clair3/ATCC_25922_fixture.hac.25x.{arm}.tsv").exists(),
          f"Arm {arm} Clair3 benchmark recorded")
for arm in ARMS:
    check((calls / f"{arm}.filter.vcf.gz").exists(), f"Arm {arm} filtered VCF written")

# 7. Arms B, C and D share one lr:hq alignment, and Arm A has its own map-ont one: the
# Read set was aligned exactly twice for calling (plus the primary-only subsampling BAM).
bams = sorted(
    p.name for p in (outdir / "work/align/ATCC_25922_fixture/hac/25x").glob("*.bam")
    if not p.name.endswith(".rg.bam")
)
check(bams == ["minimap2-2.26-map-ont.bam", "minimap2-2.31-lrhq.bam"],
      f"Read set aligned once per aligner/version/preset: {bams}")

if failures:
    sys.exit(f"{len(failures)} check(s) failed")
print("Seam 1 passed")
