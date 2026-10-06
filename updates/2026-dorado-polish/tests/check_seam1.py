"""Check Seam 1's outputs: the aggregated results table, each Read set's actual depth, and the
fixture's known TPs.

    python3 check_seam1.py <OUTDIR> <expected_tp.tsv> --read-models hac sup \
        --depths 5 10 25 50 --cpu-depth 50
"""

import argparse
import csv
import gzip
import os
import re
import sys
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("outdir", type=Path)
ap.add_argument("expected", type=Path)
ap.add_argument("--read-models", nargs="+", required=True)
ap.add_argument("--depths", nargs="+", type=int, required=True)
ap.add_argument("--cpu-depth", type=int, required=True, help="the Depth Dorado is re-run on CPU at")
args = ap.parse_args()
outdir, expected_path = args.outdir, args.expected
READ_MODELS, DEPTHS = args.read_models, sorted(args.depths)
SAMPLE = "ATCC_25922_fixture"
READ_SETS = [(rm, d) for rm in READ_MODELS for d in DEPTHS]
# The known variants are checked at the largest Depth, where the fixture has enough coverage for
# every one of them to be callable.
TOP = DEPTHS[-1]
DEPTH_TOLERANCE = 0.05  # the chromosome's actual depth must be within 5% of the Depth
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


def read_tsv_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


failures = []


def check(ok, message):
    print(("PASS " if ok else "FAIL ") + message)
    if not ok:
        failures.append(message)


# 1. One row per Arm x Read model x Depth x {SNP, INDEL, ALL} x {sweep_best, default_pass}, for
# Arms A-D: hac and sup Read sets are produced independently at every Depth.
results = read_tsv(outdir / "results/tables/results.tsv")
got = sorted(
    (r["arm"], r["read_model"], int(r["depth"]), r["var_type"], r["scoring_mode"]) for r in results
)
want = sorted(
    (a, rm, d, t, m)
    for a in ARMS
    for rm, d in READ_SETS
    for t in ("SNP", "INDEL", "ALL")
    for m in ("sweep_best", "default_pass")
)
check(
    got == want,
    f"results.tsv has {len(got)} rows, want Arm x Read model x Depth x variant type x scoring "
    f"mode ({len(want)})",
)
if got != want:
    print("  missing:", sorted(set(want) - set(got))[:10])
    print("  unexpected:", sorted(set(got) - set(want))[:10])
check({r["sample"] for r in results} == {SAMPLE}, "results.tsv is for the fixture Sample only")

# 2. Each Arm's aligner, caller and Calling model come from the config, so results.tsv carries
# them exactly: the caller's version, the Calling model (which follows the Read model), the
# digest of the container each Clair3 ran in (Dorado isn't containerised) and the hardware.
# The jobs record nothing, so there is no Calling-model checksum column.
check(
    "calling_model_sha256" not in results[0],
    "results.tsv has no calling_model_sha256 column (Clair3's model files are checked at download)",
)
DORADO_MODEL = "dna_r10.4.1_e8.2_400bps_polish_bacterial_methylation_v5.0.0"
CLAIR3_CONTAINER = {
    "1.0.5": "docker://quay.io/mbhall88/clair3@sha256:6a4c352e7d14ebb67bdad4be366dcbf6afa66cc7f0ee0b37a05cebe3f4754735",
    "2.0.3": "docker://hkubal/clair3@sha256:d56df84c2c7c508aaf0623b124d9926d366a58f34236505b26c055f848afa202",
}
by_arm = {
    (a, rm): [r for r in results if r["arm"] == a and r["read_model"] == rm]
    for a in ARMS
    for rm in READ_MODELS
}
for (arm, rm), rows in by_arm.items():
    aln_version, preset, caller, caller_version = ARM_SPEC[arm]
    check(
        len({(r["calling_model"], r["container"], r["hardware"]) for r in rows}) == 1,
        f"Arm {arm} {rm}: one Calling model, container and hardware across Depths",
    )
    r = rows[0]
    check(
        (r["aligner"], r["aligner_version"], r["preset"]) == ("minimap2", aln_version, preset),
        f"Arm {arm} {rm} aligner {r['aligner']} {r['aligner_version']} {r['preset']}",
    )
    check(r["caller"] == caller, f"Arm {arm} {rm} caller {r['caller']!r}")
    check(
        r["caller_version"] == caller_version,
        f"Arm {arm} {rm} caller_version {r['caller_version']!r}",
    )
    check(
        r["basecall_model"] == f"dna_r10.4.1_e8.2_400bps_{rm}@v4.3.0",
        f"Arm {arm} {rm} basecall_model {r['basecall_model']!r}",
    )
    if caller == "clair3":
        check(
            r["calling_model"] == f"r1041_e82_400bps_{rm}_v430",
            f"Arm {arm} {rm} calling_model {r['calling_model']!r}",
        )
        check(
            r["container"] == CLAIR3_CONTAINER[caller_version],
            f"Arm {arm} {rm} container is Clair3 {caller_version}'s pinned digest: {r['container']!r}",
        )
        check(r["hardware"].startswith("CPU: "), f"Arm {arm} {rm} hardware {r['hardware']!r}")
    else:
        check(
            r["calling_model"] == DORADO_MODEL and r["container"] == "",
            f"Arm {arm} {rm} calling_model {r['calling_model']!r}, container {r['container']!r}",
        )
        # The config's Dorado model is the one dorado polish --bacteria resolves for these
        # reads: its log says so (the rule itself no longer parses it).
        log = (outdir / f"work/logs/call_dorado/{SAMPLE}.{rm}.{TOP}x.{arm}.log").read_text()
        resolved = re.findall(r"Resolved model from input data: (\S+)", log)
        check(
            resolved and resolved[-1] == r["calling_model"],
            f"Arm {arm} {rm}: dorado resolved {resolved[-1:]} == the config's {r['calling_model']!r}",
        )

# Arms A and B run Clair3 1.0.5 in one container, and Arm C runs 2.0.3 in another. Every Read
# model has its own Clair3 Calling model.
for rm in READ_MODELS:
    digest = {a: by_arm[(a, rm)][0]["container"] for a in ("A", "B", "C")}
    check(digest["A"] == digest["B"], f"{rm}: Arms A and B ran in the same Clair3 1.0.5 container")
    check(digest["C"] != digest["B"], f"{rm}: Arm C ran in the Clair3 2.0.3 container")
check(
    len({by_arm[("B", rm)][0]["calling_model"] for rm in READ_MODELS}) == len(READ_MODELS),
    "Arm B has a different Calling model for each Read model",
)

# 3. Each Read set's actual depth is reported per contig, next to its Depth, and the
# chromosome's is within 5% of the Depth (the Depth is a per-position cap spread evenly over
# every contig, ADR-0005).
depth_rows = read_tsv(outdir / "results/tables/depth.tsv")
check(
    sorted((r["read_model"], int(r["depth"])) for r in depth_rows) == sorted(READ_SETS),
    f"depth.tsv has one row per Read set (the fixture is a single contig): {len(depth_rows)} rows",
)
for r in depth_rows:
    rm, cap, actual = r["read_model"], int(r["depth"]), float(r["mean_depth"])
    check(
        r["contig"] == "chromosome" and abs(actual - cap) <= DEPTH_TOLERANCE * cap,
        f"{rm} {cap}x: chromosome actual depth {actual:.2f}x is {actual / cap - 1:+.1%} from the "
        f"Depth (limit +/-{DEPTH_TOLERANCE:.0%})",
    )
for r in results:
    by_contig = dict(kv.split("=") for kv in r["actual_depth_by_contig"].split(";"))
    cap, actual = int(r["depth"]), float(r["actual_depth"])
    if not (
        set(by_contig) == {"chromosome"}
        and abs(float(by_contig["chromosome"]) - actual) < 0.01
        and abs(actual - cap) <= DEPTH_TOLERANCE * cap
    ):
        check(False, f"results.tsv {r['arm']} {r['read_model']} {cap}x: actual depth {actual}, "
                     f"by contig {r['actual_depth_by_contig']}")
check(True, "results.tsv carries each Read set's actual depth, matching depth.tsv and the Depth")


# 4. The fixture's known variants are true positives (all records, no QUAL threshold) at the
# largest Depth, bar the few each caller is known to miss in the repeats (known_fn.tsv).
# vcfdist's truth.tsv is 0-based and writes indels without the VCF anchor base.
def vcfdist_key(contig, pos, ref, alt):
    pos = int(pos)
    if len(ref) == len(alt):
        return (contig, str(pos - 1), ref, alt)
    return (contig, str(pos), ref[1:], alt[1:])


expected = read_tsv(expected_path)


def known_fn(read_model, caller):
    """The known variants this caller may miss on this Read model's fixture reads
    (tests/fixture/known_fn.tsv).

    All of them sit in 20.0-23.3 kb, a dense cluster (24 truth variants in 3.3 kb) inside
    ATCC_25922's pair of 7.5 kb, 96.5%-identical repeats, whose reads multi-map and leave mixed
    bases at the variants. That is the cause the eLife paper's Appendix 2 gives for Clair3's
    outlier on this Sample, where most callers missed 45 of 47 SNP FNs in the repeats; Clair3
    (1.0.5 and 2.0.3 alike) calls such sites as low-confidence hets, which `--haploid_precise`
    drops. So it isn't a workflow fault, and the three Clair3 Arms miss exactly the same ones.
    The list is what each caller misses at the largest Depth on the fixture reads, which are
    fixed (a seeded subsample); regenerate it if the fixture changes. A caller that catches
    some of the listed variants passes (with a NOTE); one that misses any other variant fails.
    """
    keys = {
        (e["CONTIG"], e["POS"], e["REF"], e["ALT"])
        for e in read_tsv(expected_path.with_name("known_fn.tsv"))
        if (e["READ_MODEL"], e["CALLER"]) == (read_model, caller)
    }
    assert keys <= {(e["CONTIG"], e["POS"], e["REF"], e["ALT"]) for e in expected}
    return keys


def is_snp(e):
    return len(e["REF"]) == len(e["ALT"])


by_mode = {(r["arm"], r["read_model"], int(r["depth"]), r["var_type"], r["scoring_mode"]): r
           for r in results}
for rm in READ_MODELS:
    for arm in ARMS:
        truth = read_tsv(outdir / f"work/score/{SAMPLE}/{rm}/{TOP}x/{arm}/sweep/truth.tsv")
        status = {(r["CONTIG"], r["POS"], r["REF"], r["ALT"]): r["ERRTYPE"] for r in truth}
        allowed_fn = known_fn(rm, ARM_SPEC[arm][2])
        required = [
            e for e in expected if (e["CONTIG"], e["POS"], e["REF"], e["ALT"]) not in allowed_fn
        ]
        missed = [
            e
            for e in required
            if status.get(vcfdist_key(e["CONTIG"], e["POS"], e["REF"], e["ALT"])) != "TP"
        ]
        if allowed_fn:
            n_fn = sum(status.get(vcfdist_key(*key)) != "TP" for key in allowed_fn)
            print(f"NOTE {rm} {TOP}x Arm {arm} misses {n_fn} of the {len(allowed_fn)} known FNs")
        check(
            not missed,
            f"{rm} {TOP}x Arm {arm}: {len(required)} of the {len(expected)} known fixture "
            f"variants are TPs (missed: {missed})",
        )

        # ... and the Default-PASS rows (every PASS record, unthresholded) count them.
        for var_type, pick in (("SNP", is_snp), ("INDEL", lambda e: not is_snp(e)), ("ALL", bool)):
            n = sum(1 for e in required if pick(e))
            tp = int(by_mode[(arm, rm, TOP, var_type, "default_pass")]["truth_tp"])
            check(tp >= n, f"{rm} {TOP}x Arm {arm} {var_type} default_pass truth_tp {tp} >= {n} required TPs")

# 5. The timing-only CPU run of Dorado is checked against the main calls, not scored, at the
# Depth(s) in cpu_run. Clair3 has no CPU re-run: it only runs on CPU.
for rm, d in READ_SETS:
    calls = outdir / f"results/calls/{SAMPLE}/{rm}/{d}x"
    cpu_bench = lambda arm: outdir / f"results/benchmarks/call_dorado_cpu/{SAMPLE}.{rm}.{d}x.{arm}.tsv"
    if d == args.cpu_depth:
        cmp = read_tsv(calls / "D.cpu_vs_main.tsv")
        check(len(cmp) == 1 and cmp[0]["same_calls"] == "yes",
              f"{rm} {d}x: CPU Dorado calls match the main run: {cmp}")
        check(cpu_bench("D").exists(), f"{rm} {d}x: CPU Dorado benchmark recorded")
    else:
        check(not (calls / "D.cpu_vs_main.tsv").exists(), f"{rm} {d}x: no CPU re-run of Dorado")
    for arm in "ABC":
        check(not (calls / f"{arm}.cpu_vs_main.tsv").exists(), f"{rm} {d}x: Arm {arm} has no CPU re-run")
        check(not cpu_bench(arm).exists(), f"{rm} {d}x: Arm {arm} has no Dorado CPU benchmark")

    # 6. Clair3's runtime is benchmarked for Arms A-C, and the Filter chain output exists for all.
    for arm in "ABC":
        check((outdir / f"results/benchmarks/call_clair3/{SAMPLE}.{rm}.{d}x.{arm}.tsv").exists(),
              f"{rm} {d}x: Arm {arm} Clair3 benchmark recorded")
    for arm in ARMS:
        check((calls / f"{arm}.filter.vcf.gz").exists(), f"{rm} {d}x: Arm {arm} filtered VCF written")

    # 7. Arms B, C and D share one lr:hq alignment, and Arm A has its own map-ont one: each
    # Read set was aligned exactly twice for calling (plus the primary-only subsampling BAM).
    bams = sorted(
        p.name for p in (outdir / f"work/align/{SAMPLE}/{rm}/{d}x").glob("*.bam")
        if not p.name.endswith(".rg.bam")
    )
    check(bams == ["minimap2-2.26-map-ont.bam", "minimap2-2.31-lrhq.bam"],
          f"{rm} {d}x: Read set aligned once per aligner/version/preset: {bams}")

# The VCFs a run commits carry no site paths: the Filter chain swaps the run's directories
# (here all under OUTDIR) for placeholders such as {work_dir}.
site_dirs = {str(outdir), os.path.realpath(outdir)}
vcfs = sorted((outdir / "results").rglob("*.vcf.gz"))
leaky = []
for vcf in vcfs:
    with gzip.open(vcf, "rt") as fh:
        header = "".join(line for line in fh if line.startswith("##"))
    if any(d in header for d in site_dirs):
        leaky.append(str(vcf.relative_to(outdir)))
check(vcfs and not leaky, f"{len(vcfs)} filtered VCF headers carry no site paths (leaky: {leaky})")

# 8. benchmarks.tsv: one alignment and one calling row for every Read set x Arm, plus the
# timing-only Dorado CPU re-run (an extra Arm D calling row) at the CPU Depth only, and sane
# numbers: positive times and memory, the configured threads, and one hardware model per
# device (a config value). Dorado runs on CPU here (the fixture has no GPU), so the main Arm D
# row and the re-run are told apart by timing_only.
bench = read_tsv(outdir / "results/tables/benchmarks.tsv")
got = sorted((r["read_model"], int(r["depth"]), r["arm"], r["step"], r["timing_only"]) for r in bench)
want = sorted(
    [(rm, d, a, step, "false") for rm, d in READ_SETS for a in ARMS for step in ("align", "call")]
    + [(rm, d, "D", "call", "true") for rm, d in READ_SETS if d == args.cpu_depth]
)
check(got == want, f"benchmarks.tsv has {len(want)} rows: align and call per Arm, plus the CPU re-run at {args.cpu_depth}x")
check({r["sample"] for r in bench} == {SAMPLE}, "benchmarks.tsv is for the fixture Sample")
bad = [
    (r["read_model"], r["depth"], r["arm"], r["step"], r["timing_only"])
    for r in bench
    if not (
        float(r["wall_time_s"]) > 0
        and float(r["max_rss_mb"]) > 0
        and float(r["cpu_time_s"]) > 0
        and re.fullmatch(r"\d+", r["threads"])
        and r["hardware"]
    )
]
check(not bad, f"every benchmark row has positive times and memory, threads and hardware: {bad}")
check(
    all(r["hardware"].startswith("CPU: ") and r["device"] == "cpu" for r in bench),
    "all fixture timings are on CPU, recorded as CPU: <model>",
)
check(len({r["hardware"] for r in bench}) == 1, f"one hardware model: {sorted({r['hardware'] for r in bench})}")
check(
    "driver" not in bench[0] and "host" not in bench[0] and "command_wall_time_s" not in bench[0],
    "benchmarks.tsv has no host, driver or command_wall_time_s column (nothing is run in the jobs)",
)
for r in bench:
    if r["timing_only"] == "true" and r["arm"] != "D":
        check(False, f"timing-only row for Arm {r['arm']}")
check(
    {r["tool"] for r in bench if r["step"] == "call"} == {"clair3", "dorado"}
    and {r["tool"] for r in bench if r["step"] == "align"} == {"minimap2"},
    "benchmark rows name their tool",
)
# Dorado's CPU-only re-run never becomes an accuracy row: results.tsv has exactly one Arm D row
# per Read set x variant type x scoring mode.
d_rows = [r for r in results if r["arm"] == "D"]
check(
    len(d_rows) == len(READ_SETS) * 3 * 2,
    f"results.tsv has no extra Arm D rows from the CPU re-run ({len(d_rows)} rows)",
)

# 9. The figures and tables (#14) exist, are non-empty, and the tables have the rows Table 1 and
# Table S1 promise. (What the figures look like is checked by eye on the full run.)
PNG_MAGIC = b"\x89PNG\r\n\x1a\n"
for name in ("fig1_best_f1_depth", "fig2_pr_curves", "fig3_per_sample_best_f1", "fig4_runtime_memory"):
    png = outdir / f"results/figures/{name}.png"
    svg = outdir / f"results/figures/{name}.svg"
    check(png.exists() and png.read_bytes()[:8] == PNG_MAGIC, f"{name}.png is a PNG")
    check(svg.exists() and b"<svg" in svg.read_bytes()[:2000], f"{name}.svg is an SVG")

table1 = read_tsv_csv(outdir / "results/tables/table1_runtime_memory.csv")
tools = {(r["Tool"], r["Version"], r["Device"], r["Timing only"]) for r in table1}
for want_tool in (
    ("minimap2", "2.26", "CPU", "no"),
    ("minimap2", "2.31", "CPU", "no"),
    ("Clair3", "1.0.5", "CPU", "no"),
    ("Clair3", "2.0.3", "CPU", "no"),
    ("dorado polish", "2.1.2", "CPU", "no"),  # the fixture's main Dorado run is on CPU, not GPU
    ("dorado polish", "2.1.2", "CPU", "yes"),
):
    check(want_tool in tools, f"Table 1 has {want_tool}")
check(
    {r["Depth (x)"] for r in table1 if r["Timing only"] == "yes"} == {str(args.cpu_depth)},
    f"Table 1's Dorado CPU row is only at {args.cpu_depth}x",
)
# The shared alignment is one row for Arms B, C and D, not one each.
shared = [r for r in table1 if r["Step"] == "Alignment" and r["Version"] == "2.31"]
check(
    len(shared) == len(DEPTHS) and all(r["Arms"].count(",") == 2 for r in shared),
    "Table 1 counts the alignment shared by Arms B, C and D once",
)
check(
    all(float(r["Wall time median (s)"]) > 0 and float(r["Peak RSS median (MB)"]) > 0 for r in table1),
    "Table 1's times and memory are positive",
)

# The AF filter's rows (the run has the analysis on, #28) are checked by check_af_figures.py;
# here, the Arms' rows alone.
table_s1 = [
    r
    for r in read_tsv_csv(outdir / "results/tables/table_s1_per_sample.csv")
    if "AF filter" not in r["Arm"]
]
check(
    len(table_s1) == len(READ_SETS) * len(ARMS),
    f"Table S1 has a row per Read set x Arm ({len(table_s1)} rows)",
)
check(
    all(r["Actual depth by contig (x)"].startswith("chromosome ") for r in table_s1),
    "Table S1 gives the actual depth of every contig",
)
check(
    {r["Arm"] for r in table_s1} == {"A (paper)", "B (lr:hq post)", "C (current Clair3)", "D (Dorado)"},
    "Table S1 labels the Arms as in CONTEXT.md",
)
by_row = {(r["Read model"], r["Depth (x)"], r["Arm"][0]): r for r in table_s1}
res_row = {
    (r["read_model"], r["depth"], r["arm"], r["var_type"], r["scoring_mode"]): r for r in results
}
check(
    all(
        by_row[(rm, str(d), a)]["INDEL Best F1"] == res_row[(rm, str(d), a, "INDEL", "sweep_best")]["f1"]
        and by_row[(rm, str(d), a)]["INDEL Default-PASS F1"]
        == res_row[(rm, str(d), a, "INDEL", "default_pass")]["f1"]
        for rm, d in READ_SETS
        for a in ARMS
    ),
    "Table S1 keeps Best F1 and the Default-PASS score apart, as in results.tsv",
)

if failures:
    sys.exit(f"{len(failures)} check(s) failed")
print("Seam 1 passed")
