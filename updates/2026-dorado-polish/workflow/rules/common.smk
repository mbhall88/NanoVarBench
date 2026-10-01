import csv
import re
from pathlib import Path


# Containers, pinned by digest (#6). The tag each digest was taken from is in the comment.
CONTAINERS = {
    # mulled minimap2 2.26 + samtools 1.17, so Arm A's aligner can pipe into samtools sort
    # (tag 7e6194c85b2f194e301c71cdda1c002a754e8cc1-0)
    "minimap2-2.26": "docker://quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4@sha256:450ac9cf2a6a291d7ef0d3cebe29066c67edc642d325f568ec9196b8b83d51f7",
    # mulled minimap2 2.31-r1302 + samtools 1.23.1, so the aligner can pipe into samtools sort
    # (tag b411340b52d82a9c276d87c7a3dcffc880be762f-0)
    "minimap2-2.31": "docker://quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4@sha256:966a1318a02cc3cda1785ccf62a4db2390e88dd48f303befa6c4a0a89a241a49",
    # samtools:1.24--h9dcdb79_1
    "samtools": "docker://quay.io/biocontainers/samtools@sha256:a130447589651ed09252aa95a5e4f4132942cdb54d835d81a04a9a930d656561",
    # bcftools:1.24--h118bc1c_2
    "bcftools": "docker://quay.io/biocontainers/bcftools@sha256:a3e0d3007ffe325c409b398f660840a3e7574d076219c6e82fc994ced87d47c3",
    # seqkit:2.14.0--hb192632_0
    "seqkit": "docker://quay.io/biocontainers/seqkit@sha256:45fb535880be37dfed5be5517111fb8bfdd6234ef36e725b126a6131b1af2ef0",
    # rasusa:5.1.0--hfa8f182_0
    "rasusa": "docker://quay.io/biocontainers/rasusa@sha256:e8c7b92c66abd96fdb861a4b7977c6ac47f1e9df6428d539a6616bb257bb0846",
    # quay.io/mbhall88/clair3:1.0.5, the paper's Clair3 (Arms A and B); bundles the TF
    # r1041_e82_400bps_{hac,sup}_v430 Calling models in /opt/models
    "clair3-1.0.5": "docker://quay.io/mbhall88/clair3@sha256:6a4c352e7d14ebb67bdad4be366dcbf6afa66cc7f0ee0b37a05cebe3f4754735",
    # hkubal/clair3:v2.0.3 (Arm C); the v4.3.0 PyTorch Calling models aren't bundled
    "clair3-2.0.3": "docker://hkubal/clair3@sha256:d56df84c2c7c508aaf0623b124d9926d366a58f34236505b26c055f848afa202",
    # timd1/vcfdist:v2.6.4 (ADR-0003)
    "vcfdist": "docker://timd1/vcfdist@sha256:d8b14a999a290f3b21dd4cde3bf52f2ad814b252823b8a4d9a01b548ae71dee3",
}

WORKFLOW_DIR = Path(workflow.basedir)
SCRIPTS = WORKFLOW_DIR / "scripts"
ENVS = WORKFLOW_DIR / "envs"

WORK = Path(config["work_dir"])
RESULTS = Path(config["results_dir"])
LOGS = WORK / "logs"
READS_DIR = Path(config["reads_dir"])
TRUTH_DIR = Path(config["truth_dir"])


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


SAMPLES = {row["sample"]: row for row in read_tsv(config["samples"])}
RUNS = {(row["sample"], row["read_model"]): row for row in read_tsv(config["runs"])}
RUN_URLS = {row["run"]: row for row in RUNS.values()}

RUN_SAMPLES = config["run"]["samples"]
RUN_READ_MODELS = config["run"]["read_models"]
RUN_DEPTHS = [int(d) for d in config["run"]["depths"]]
RUN_ARMS = config["run"]["arms"]
ARMS = config["arms"]
SCORING_MODES = ["sweep", "pass"]

for s in RUN_SAMPLES:
    if s not in SAMPLES:
        raise ValueError(f"Sample {s} is not in {config['samples']}")
    for rm in RUN_READ_MODELS:
        if (s, rm) not in RUNS:
            raise ValueError(f"No {rm} run for Sample {s} in {config['runs']}")
for a in RUN_ARMS:
    if a not in ARMS:
        raise ValueError(f"Arm {a} is not defined under 'arms' in the config")
CALLERS = ("dorado", "clair3")
for a, spec in ARMS.items():
    if spec["caller"] not in CALLERS:
        raise ValueError(f"Arm {a}: no calling rule for caller {spec['caller']}")
    if spec["caller"] == "clair3" and f"clair3-{spec['caller_version']}" not in CONTAINERS:
        raise ValueError(f"Arm {a}: no container pinned for Clair3 {spec['caller_version']}")


def alignment_id(aligner, version, preset):
    """Name an alignment by aligner, version and preset, so Arms that share all three
    (B, C and D) share one BAM."""
    return f"{aligner}-{version}-{preset.replace(':', '')}"


# Every alignment the config refers to: each Arm's, the subsampling alignment and the one
# actual depth is measured on.
ALIGNMENTS = {}
for spec in [
    *ARMS.values(),
    config["subsample"]["alignment"],
    config["actual_depth"]["alignment"],
]:
    aid = alignment_id(spec["aligner"], str(spec["aligner_version"]), spec["preset"])
    ALIGNMENTS[aid] = {
        "aligner": spec["aligner"],
        "version": str(spec["aligner_version"]),
        "preset": spec["preset"],
    }


def arm_alignment(arm):
    a = ARMS[arm]
    return alignment_id(a["aligner"], str(a["aligner_version"]), a["preset"])


SUBSAMPLE_ALN = alignment_id(
    config["subsample"]["alignment"]["aligner"],
    str(config["subsample"]["alignment"]["aligner_version"]),
    config["subsample"]["alignment"]["preset"],
)
DEPTH_ALN = alignment_id(
    config["actual_depth"]["alignment"]["aligner"],
    str(config["actual_depth"]["alignment"]["aligner_version"]),
    config["actual_depth"]["alignment"]["preset"],
)


def aligner_container(wildcards):
    spec = ALIGNMENTS[wildcards.aln]
    key = f"{spec['aligner']}-{spec['version']}"
    if key not in CONTAINERS:
        raise ValueError(f"No container pinned for {key}")
    return CONTAINERS[key]


def arms_with_caller(caller):
    return [a for a in RUN_ARMS if ARMS[a]["caller"] == caller]


def arm_pattern(caller):
    """Wildcard constraint for the Arms (configured, not just run) that use a caller, so a
    caller's rules never match another caller's Arms."""
    arms = [a for a, spec in ARMS.items() if spec["caller"] == caller]
    return "|".join(map(re.escape, arms)) or "(?!)"


def basecall_model(wildcards):
    return config["read_models"][wildcards.read_model]


def raw_reads(wildcards):
    run = RUNS[(wildcards.sample, wildcards.read_model)]["run"]
    return READS_DIR / wildcards.read_model / f"{run}_1.fastq.gz"


wildcard_constraints:
    sample="|".join(map(re.escape, SAMPLES)),
    read_model="|".join(map(re.escape, config["read_models"])),
    depth=r"\d+",
    arm="|".join(map(re.escape, ARMS)),
    aln="|".join(map(re.escape, ALIGNMENTS)),
    mode="|".join(SCORING_MODES),
    run=r"[A-Za-z0-9]+",


# Path stems shared across rule files.
READ_SET = WORK / "reads/{sample}/{read_model}/{depth}x"
ALIGN = WORK / "align/{sample}/{read_model}/{depth}x"
CALL = WORK / "call/{sample}/{read_model}/{depth}x/{arm}"
SCORE = WORK / "score/{sample}/{read_model}/{depth}x/{arm}/{mode}"
KEYS = "{sample}.{read_model}.{depth}x"
# Every caller writes one of these per Arm: caller, version, Calling model and its
# checksums, and where it ran.
CALLER_INFO = RESULTS / "calls/{sample}/{read_model}/{depth}x/{arm}.caller_info.tsv"
