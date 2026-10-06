# Calling: one rule per caller family (Dorado, Clair3), parameterised by Arm. Each rule only
# runs its tool, since Snakemake's benchmark: directive times the whole job script (#12). What
# ran (caller, version, Calling model, container, hardware) is read from the config by the
# aggregation, not recorded per job.

DORADO = config["dorado"]
DORADO_MODELS_DIR = Path(DORADO["models_dir"])


def dorado_bin(arm):
    version = str(ARMS[arm]["caller_version"])
    try:
        return DORADO["bin"][version]
    except KeyError:
        raise ValueError(f"No dorado binary configured for version {version} (Arm {arm})")


def dorado_model(arm):
    return DORADO["models"][ARMS[arm]["calling_model"]]


DORADO_THREADS = 8
DORADO_DEVICE = "cpu" if DORADO["device"] == "cpu" else "gpu"  # for the benchmark hardware


rule download_dorado_model:
    """Fetch a Dorado model (a polishing model, or the smallvar pilot's) into models_dir, so
    GPU nodes without internet can read it locally. Skipped when the model is already there.
    `dorado download` skips a model whose directory exists, even empty, and Snakemake makes
    the output's directory before the job runs, so it downloads to a temporary directory and
    moves the files in (#32)."""
    output:
        config=DORADO_MODELS_DIR / "{model}/config.toml",
        weights=DORADO_MODELS_DIR / "{model}/weights.pt",
    log:
        LOGS / "download_dorado_model/{model}.log",
    resources:
        mem_mb=2000,
        runtime=30,
    params:
        bin=DORADO["bin"][str(DORADO["download_with"])],
        models_dir=DORADO_MODELS_DIR,
    shell:
        """
        exec 2> {log}
        tmp=$(mktemp -d -p {params.models_dir}); trap 'rm -rf "$tmp"' EXIT
        {params.bin} download --model '{wildcards.model}' --models-directory "$tmp"
        mv "$tmp/{wildcards.model}/config.toml" '{output.config}'
        mv "$tmp/{wildcards.model}/weights.pt" '{output.weights}'
        """


localrules:
    download_dorado_model,


rule call_dorado:
    """dorado polish --vcf on the RG-reheadered minimap2 BAM (ADR-0004), with --min-depth 2
    to match Clair3 (ADR-0002). Only the tool runs here: the job is timed (#12)."""
    input:
        bam=lambda wc: ALIGN / f"{arm_alignment(wc.arm)}.rg.bam",
        bai=lambda wc: ALIGN / f"{arm_alignment(wc.arm)}.rg.bam.bai",
        mutref=rules.extract_truth.output.mutref,
        faidx=rules.index_mutref.output.faidx,
        model=lambda wc: DORADO_MODELS_DIR / dorado_model(wc.arm) / "weights.pt",
    output:
        vcf=CALL / "variants.vcf",
    wildcard_constraints:
        arm=arm_pattern("dorado"),
    log:
        LOGS / "call_dorado/{sample}.{read_model}.{depth}x.{arm}.log",
    benchmark:
        BENCH_DORADO
    threads: DORADO_THREADS
    resources:
        mem_mb=16000,
        runtime=30,
    params:
        bin=lambda wc: dorado_bin(wc.arm),
        model_flag=lambda wc: f"--{ARMS[wc.arm]['calling_model']}",
        models_dir=DORADO_MODELS_DIR,
        device=DORADO["device"],
        min_depth=DORADO["min_depth"],
        any_bam="--any-bam",  # a dorado aligner BAM doesn't need it (#13)
    shell:
        """
        outdir=$(dirname {output.vcf})
        {params.bin} polish {input.bam} {input.mutref} {params.model_flag} --vcf \
            --min-depth {params.min_depth} {params.any_bam} --ignore-read-groups \
            --models-directory {params.models_dir} --threads {threads} \
            --device {params.device} -o "$outdir" -v > "$outdir/stdout.txt" 2> {log}
        """


# Timing-only CPU run of the same Dorado command (#12), on the Depths in dorado.cpu_run.
# It isn't scored: compare_dorado_devices checks its calls against call_dorado's instead.
CPU_RUN = DORADO.get("cpu_run") or {}
DORADO_CPU_THREADS = CPU_RUN.get("threads", 8)
CALL_CPU = WORK / "call_cpu/{sample}/{read_model}/{depth}x/{arm}"


use rule call_dorado as call_dorado_cpu with:
    output:
        vcf=CALL_CPU / "variants.vcf",
    log:
        LOGS / "call_dorado_cpu/{sample}.{read_model}.{depth}x.{arm}.log",
    benchmark:
        BENCH_DORADO_CPU
    threads: DORADO_CPU_THREADS
    resources:
        mem_mb=16000,
        runtime=120,
    params:
        bin=lambda wc: dorado_bin(wc.arm),
        model_flag=lambda wc: f"--{ARMS[wc.arm]['calling_model']}",
        models_dir=DORADO_MODELS_DIR,
        device="cpu",
        min_depth=DORADO["min_depth"],
        any_bam="--any-bam",


rule compare_dorado_devices:
    """Check that the CPU run's calls match call_dorado's: same records (CHROM, POS, REF,
    ALT, FILTER, GT), with QUAL allowed to differ (the GPU runs the model in half
    precision, #7). Writes one row with the record counts and the largest QUAL change."""
    input:
        main=CALL / "variants.vcf",
        cpu=CALL_CPU / "variants.vcf",
    output:
        tsv=RESULTS / "calls/{sample}/{read_model}/{depth}x/{arm}.cpu_vs_main.tsv",
    log:
        LOGS / "compare_dorado_devices/{sample}.{read_model}.{depth}x.{arm}.log",
    resources:
        mem_mb=1000,
        runtime=10,
    shell:
        """
        exec 2> {log}
        tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
        for f in main cpu; do
            vcf={input.main}; [ $f = cpu ] && vcf={input.cpu}
            grep -v '^#' "$vcf" | awk -F'\\t' -v OFS='\\t' \
                '{{split($10, s, ":"); print $1":"$2":"$4":"$5":"$7":"s[1], $6}}' | sort > "$tmp/$f"
        done
        join -t $'\\t' "$tmp/main" "$tmp/cpu" > "$tmp/both"
        n_main=$(wc -l < "$tmp/main"); n_cpu=$(wc -l < "$tmp/cpu"); n_both=$(wc -l < "$tmp/both")
        max_dq=$(awk -F'\\t' '{{d = $2 - $3; if (d < 0) d = -d; if (d > m) m = d}} END {{printf "%.4f", m}}' "$tmp/both")
        same=no; [ "$n_main" -eq "$n_cpu" ] && [ "$n_main" -eq "$n_both" ] && same=yes
        printf 'main_records\\tcpu_records\\tshared_records\\tsame_calls\\tmax_abs_qual_diff\\n%s\\t%s\\t%s\\t%s\\t%s\\n' \
            "$n_main" "$n_cpu" "$n_both" "$same" "$max_dq" > {output.tsv}
        """


# --- Clair3 (Arms A, B and C) ---------------------------------------------------------
# The same options for every Clair3 Arm, matching the paper's rule (workflow/scripts/
# callers/clair3.sh at the repo root) for haploid bacterial calling. --min_coverage is left
# at its default of 2 in both versions, which is the Dorado --min-depth 2 (ADR-0002), and
# the Fine-tuned model is never used (ADR-0001).

CLAIR3 = config["clair3"]
CLAIR3_MODELS_DIR = Path(CLAIR3["models_dir"])
CLAIR3_MODEL_FILES = ["pileup.pt", "full_alignment.pt"]  # an HKU PyTorch model's files
CLAIR3_SOURCES = ("bundled_tf", "hku_pytorch")

for a, spec in ARMS.items():
    if spec["caller"] == "clair3" and spec["calling_model"] not in CLAIR3_SOURCES:
        raise ValueError(
            f"Arm {a}: calling_model must be one of {CLAIR3_SOURCES}, not {spec['calling_model']}"
        )


def clair3_model(wildcards):
    """The Calling model follows the Read model: hac reads use the hac model."""
    return CLAIR3["calling_models"][wildcards.read_model]


def clair3_model_path(wildcards):
    if ARMS[wildcards.arm]["calling_model"] == "bundled_tf":
        return f"{CLAIR3['bundled_models_dir']}/{clair3_model(wildcards)}"  # inside the container
    return str(CLAIR3_MODELS_DIR / clair3_model(wildcards))


def clair3_model_files(wildcards):
    """The downloaded model files an Arm needs; none for the models bundled in the image."""
    if ARMS[wildcards.arm]["calling_model"] == "bundled_tf":
        return []
    return [CLAIR3_MODELS_DIR / clair3_model(wildcards) / f for f in CLAIR3_MODEL_FILES]


def clair3_container(wildcards):
    return CONTAINERS[f"clair3-{ARMS[wildcards.arm]['caller_version']}"]


rule download_clair3_model:
    """Download an HKU-converted PyTorch Calling model (for Clair3 2.x) and check each file
    against the SHA256 pinned in the config. HKU publishes no checksums, so the pins were
    recorded from a download on 2026-10-01 (the files' sizes and hashes were identical on a
    second download). Apptainer bind-mounts models_dir into the Clair3 container."""
    output:
        files=[CLAIR3_MODELS_DIR / "{model}" / f for f in CLAIR3_MODEL_FILES],
    wildcard_constraints:
        model="|".join(map(re.escape, CLAIR3["model_sha256"])),
    log:
        LOGS / "download_clair3_model/{model}.log",
    resources:
        mem_mb=1000,
        runtime=15,
    params:
        url=CLAIR3["models_url"],
        sha256=lambda wc: CLAIR3["model_sha256"][wc.model],
        files=CLAIR3_MODEL_FILES,
    run:
        import hashlib
        import time
        import urllib.request

        outdir = Path(output.files[0]).parent
        outdir.mkdir(parents=True, exist_ok=True)
        with open(log[0], "w") as lg:
            for name in params.files:
                url = f"{params.url}/{wildcards.model}/{name}"
                part = outdir / (name + ".part")
                for attempt in range(3):
                    try:
                        urllib.request.urlretrieve(url, part)
                        break
                    except OSError as err:
                        print(f"{url}: attempt {attempt + 1} failed: {err}", file=lg)
                        if attempt == 2:
                            raise
                        time.sleep(10)
                got = hashlib.sha256(part.read_bytes()).hexdigest()
                print(f"{url} sha256 {got} (want {params.sha256[name]})", file=lg)
                if got != params.sha256[name]:
                    part.unlink()
                    raise ValueError(f"{url}: sha256 {got} != pinned {params.sha256[name]}")
                part.rename(outdir / name)


localrules:
    download_clair3_model,


rule call_clair3:
    """Clair3 on the Arm's alignment: 8 CPU threads, haploid precise (the paper's options).
    Only the tool runs here: the job is timed (#12)."""
    input:
        bam=lambda wc: ALIGN / f"{arm_alignment(wc.arm)}.bam",
        bai=lambda wc: ALIGN / f"{arm_alignment(wc.arm)}.bam.bai",
        mutref=rules.extract_truth.output.mutref,
        faidx=rules.index_mutref.output.faidx,
        model_files=clair3_model_files,
    output:
        vcf=CALL / "variants.vcf.gz",
    wildcard_constraints:
        arm=arm_pattern("clair3"),
    log:
        LOGS / "call_clair3/{sample}.{read_model}.{depth}x.{arm}.log",
    benchmark:
        BENCH_CLAIR3
    threads: CLAIR3["threads"]
    resources:
        mem_mb=16000,
        runtime=120,
    params:
        model_path=clair3_model_path,
        options=" ".join(CLAIR3["options"]),
    container:
        clair3_container
    shell:
        """
        exec &> {log}
        outdir=$(mktemp -d -p $(dirname {output.vcf}) clair3.XXXXXX)
        trap 'rm -rf "$outdir"' EXIT
        [ -r {params.model_path}/pileup.pt ] || [ -r {params.model_path}/pileup.index ] || \
            {{ echo "no Calling model files in {params.model_path}" >&2; exit 1; }}

        /opt/bin/run_clair3.sh --bam_fn={input.bam} --ref_fn={input.mutref} \
            --threads={threads} --model_path={params.model_path} --output="$outdir" \
            --sample_name={wildcards.sample} {params.options}
        mv "$outdir/merge_output.vcf.gz" {output.vcf}
        """


def dorado_cpu_targets():
    depths = [int(d) for d in CPU_RUN.get("depths", [])]
    return [
        str(rules.compare_dorado_devices.output.tsv).format(
            sample=s, read_model=rm, depth=d, arm=a
        )
        for s in RUN_SAMPLES
        for rm in RUN_READ_MODELS
        for d in RUN_DEPTHS
        if d in depths
        for a in arms_with_caller("dorado")
    ]
