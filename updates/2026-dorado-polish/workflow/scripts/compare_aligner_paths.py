"""Compare Dorado's calls from two alignments of the same Read set (#13).

  path 1  the Arm's minimap2 2.31 lr:hq BAM, RG header added, dorado polish --any-bam
  path 2  a dorado aligner BAM, the same RG header, dorado polish without --any-bam

For Dorado's raw VCFs and for the filtered VCFs that go to vcfdist it counts the records
that differ (CHROM, POS, REF, ALT, FILTER, GT) and the largest QUAL difference, then compares
the vcfdist Best F1 (sweep, THRESHOLD == BEST) and Default-PASS F1 (PASS only, THRESHOLD ==
NONE) for SNP, INDEL and ALL. Writes a long-format TSV and a short Markdown report.
"""

import csv
import gzip
import sys

sys.stderr = open(snakemake.log[0], "w")

wc = snakemake.wildcards
inp = snakemake.input
material = snakemake.params.material
VAR_TYPES = ["SNP", "INDEL", "ALL"]
SCORES = {"best_f1": ("sweep", "BEST"), "default_pass_f1": ("pass", "NONE")}
SHOW = 10  # differing records listed in the report


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_vcf(path):
    """{(CHROM, POS, REF, ALT): (FILTER, GT, QUAL)}. Filtered VCFs are bgzipped, which gzip
    reads."""
    opener = gzip.open if str(path).endswith(".gz") else open
    records = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            key = (f[0], int(f[1]), f[3], f[4])
            if key in records:
                raise ValueError(f"{path}: duplicate record {key}")
            qual = float(f[5]) if f[5] != "." else float("nan")
            records[key] = (f[6], f[9].split(":")[0], qual)
    return records


def compare(path1, path2):
    r1, r2 = read_vcf(path1), read_vcf(path2)
    shared = r1.keys() & r2.keys()
    only1, only2 = sorted(r1.keys() - r2.keys()), sorted(r2.keys() - r1.keys())
    changed = sorted(k for k in shared if r1[k][:2] != r2[k][:2])  # FILTER or GT differ
    dq = {k: abs(r1[k][2] - r2[k][2]) for k in shared}
    worst = max(dq, key=dq.get) if dq else None
    return {
        "path1_records": len(r1),
        "path2_records": len(r2),
        "shared_sites": len(shared),
        "only_path1": len(only1),
        "only_path2": len(only2),
        "filter_or_gt_differs": len(changed),
        "records_differ": len(only1) + len(only2) + len(changed),
        "max_abs_qual_diff": dq[worst] if worst else 0.0,
        "max_qual_diff_at": worst,
        "details": {
            "only_path1": [(k, r1[k]) for k in only1],
            "only_path2": [(k, r2[k]) for k in only2],
            "changed": [(k, r1[k], r2[k]) for k in changed],
        },
    }


raw = compare(inp.raw1, inp.raw2)
filt = compare(inp.filter1, inp.filter2)


def scores(summary_path, threshold):
    rows = {r["VAR_TYPE"]: r for r in read_tsv(summary_path) if r["THRESHOLD"] == threshold}
    return {t: rows[t] for t in VAR_TYPES}


f1 = {}  # (score, var_type) -> (path 1 row, path 2 row)
for name, (mode, threshold) in SCORES.items():
    s1 = scores(inp[f"{mode}1"], threshold)
    s2 = scores(inp[f"{mode}2"], threshold)
    for t in VAR_TYPES:
        f1[(name, t)] = (s1[t], s2[t])

info1, info2 = read_tsv(inp.info1)[0], read_tsv(inp.info2)[0]

# --- TSV: one row per measurement -------------------------------------------------------
rows = []
for label, res in (("raw_vcf", raw), ("filtered_vcf", filt)):
    for k in (
        "path1_records",
        "path2_records",
        "shared_sites",
        "only_path1",
        "only_path2",
        "filter_or_gt_differs",
        "records_differ",
        "max_abs_qual_diff",
    ):
        rows.append({"section": label, "metric": k, "var_type": "", "path1": "", "path2": "",
                     "value": res[k] if k != "max_abs_qual_diff" else f"{res[k]:.4f}"})
for (name, t), (a, b) in f1.items():
    for col, metric in (("F1_SCORE", name), ("PREC", name.replace("f1", "precision")),
                        ("RECALL", name.replace("f1", "recall")), ("TRUTH_FN", name.replace("f1", "fn")),
                        ("QUERY_FP", name.replace("f1", "fp"))):
        va, vb = float(a[col]), float(b[col])
        rows.append({"section": "vcfdist", "metric": metric, "var_type": t,
                     "path1": a[col], "path2": b[col], "value": f"{vb - va:.6f}"})
best_diff = max(abs(float(b["F1_SCORE"]) - float(a["F1_SCORE"]))
                for (name, t), (a, b) in f1.items() if name == "best_f1")
pass_diff = max(abs(float(b["F1_SCORE"]) - float(a["F1_SCORE"]))
                for (name, t), (a, b) in f1.items() if name == "default_pass_f1")
is_material = filt["records_differ"] > material["max_records"] or best_diff > material["max_best_f1_diff"]
rows.append({"section": "summary", "metric": "materially_different", "var_type": "",
             "path1": "", "path2": "", "value": "yes" if is_material else "no"})

with open(snakemake.output.tsv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)


# --- Markdown report --------------------------------------------------------------------
def fmt_key(key):
    chrom, pos, ref, alt = key
    ref = ref if len(ref) <= 12 else ref[:12] + "..."
    alt = alt if len(alt) <= 12 else alt[:12] + "..."
    return f"{chrom}:{pos} {ref}>{alt}"


def detail_lines(res):
    d = res["details"]
    out = []
    for label, items in (("only in path 1", d["only_path1"]), ("only in path 2", d["only_path2"])):
        for key, (flt, gt, qual) in items[:SHOW]:
            out.append(f"  - {label}: {fmt_key(key)} FILTER={flt} GT={gt} QUAL={qual:g}")
        if len(items) > SHOW:
            out.append(f"  - ... and {len(items) - SHOW} more {label}")
    for key, a, b in d["changed"][:SHOW]:
        out.append(f"  - {fmt_key(key)}: path 1 FILTER={a[0]} GT={a[1]}, path 2 FILTER={b[0]} GT={b[1]}")
    if len(d["changed"]) > SHOW:
        out.append(f"  - ... and {len(d['changed']) - SHOW} more FILTER/GT changes")
    return out


def where(res):
    k = res["max_qual_diff_at"]
    return f" (at {fmt_key(k)})" if k and res["max_abs_qual_diff"] > 0 else ""


title = f"{wc.sample} {wc.read_model} {wc.depth}x, Arm {wc.arm}"
L = [
    f"# dorado aligner vs reheadered minimap2: {title}",
    "",
    f"**Calls differ materially: {'YES' if is_material else 'no'}** "
    f"(more than {material['max_records']} filtered records differ, or a Best F1 difference "
    f"above {material['max_best_f1_diff']}).",
    "",
    "Both paths start from the same Read set and add the same single `@RG` line. "
    "Both run `dorado polish --bacteria --vcf --min-depth 2 --ignore-read-groups` on the same "
    "device with the same Calling model.",
    "",
    "- **Path 1:** minimap2 2.31 `lr:hq` (`-aL --cs --MD`), then `dorado polish --any-bam`.",
    "- **Path 2:** `dorado aligner` (default `lr:hq`, bundled minimap2), then `dorado polish` "
    "without `--any-bam`.",
    "",
    f"Dorado: path 1 `{info1['caller_version']}`, path 2 `{info2['caller_version']}`. "
    f"Calling model: path 1 `{info1['calling_model']}`, path 2 `{info2['calling_model']}`. "
    f"Device: path 1 {info1['device']} ({info1['hardware']}), path 2 {info2['device']} "
    f"({info2['hardware']}).",
    "",
    "## Calls",
    "",
    "| VCF | path 1 records | path 2 records | records that differ | only path 1 | only path 2 "
    "| FILTER/GT differs | largest QUAL difference |",
    "|---|---|---|---|---|---|---|---|",
]
for label, res in (("Dorado raw", raw), ("after Filter chain", filt)):
    L.append(
        f"| {label} | {res['path1_records']} | {res['path2_records']} | {res['records_differ']} "
        f"| {res['only_path1']} | {res['only_path2']} | {res['filter_or_gt_differs']} "
        f"| {res['max_abs_qual_diff']:.4f} |"
    )
L += ["", "A record is identified by CHROM, POS, REF and ALT. It differs if it is in only one "
      "VCF or if its FILTER or GT differs. QUAL is compared on records in both VCFs.", ""]
L.append(f"Largest QUAL difference: raw {raw['max_abs_qual_diff']:.4f}{where(raw)}, "
         f"filtered {filt['max_abs_qual_diff']:.4f}{where(filt)}.")
for label, res in (("Dorado raw", raw), ("After the Filter chain", filt)):
    lines = detail_lines(res)
    if lines:
        L += ["", f"{label}, differing records:", *lines]
L += ["", "## vcfdist (v2.6.4, whole genome)", "",
      "Best F1 is the QUAL sweep's best row. Default-PASS F1 scores PASS records only with no "
      "threshold. Difference is path 2 minus path 1.", "",
      "| Score | Type | path 1 F1 | path 2 F1 | difference | path 1 FN / FP | path 2 FN / FP |",
      "|---|---|---|---|---|---|---|"]
for (name, t), (a, b) in f1.items():
    d = float(b["F1_SCORE"]) - float(a["F1_SCORE"])
    L.append(
        f"| {'Best F1' if name == 'best_f1' else 'Default-PASS F1'} | {t} | {a['F1_SCORE']} "
        f"| {b['F1_SCORE']} | {d:+.6f} | {a['TRUTH_FN']} / {a['QUERY_FP']} "
        f"| {b['TRUTH_FN']} / {b['QUERY_FP']} |"
    )
L += ["", f"Largest absolute Best F1 difference: {best_diff:.6f}. "
      f"Largest absolute Default-PASS F1 difference: {pass_diff:.6f}.", ""]

with open(snakemake.output.report, "w") as fh:
    fh.write("\n".join(L))
print(f"filtered records differ={filt['records_differ']} best_f1_diff={best_diff}", file=sys.stderr)
