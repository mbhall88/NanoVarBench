# Calling: one rule per caller family, parameterised by Arm.

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


rule download_dorado_model:
    """Fetch a Dorado polishing model into models_dir, so GPU nodes without internet can
    read it locally. Skipped when the model is already there."""
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
        {params.bin} download --model {wildcards.model} --models-directory {params.models_dir} 2> {log}
        """


localrules:
    download_dorado_model,


rule call_dorado:
    """dorado polish --vcf on the RG-reheadered minimap2 BAM (ADR-0004), with --min-depth 2
    to match Clair3 (ADR-0002). Records the Dorado version, the Calling model it resolved
    and the device it ran on."""
    input:
        bam=lambda wc: ALIGN / f"{arm_alignment(wc.arm)}.rg.bam",
        bai=lambda wc: ALIGN / f"{arm_alignment(wc.arm)}.rg.bam.bai",
        mutref=rules.extract_truth.output.mutref,
        faidx=rules.index_mutref.output.faidx,
        model=lambda wc: DORADO_MODELS_DIR / dorado_model(wc.arm) / "weights.pt",
    output:
        vcf=CALL / "variants.vcf",
        info=RESULTS / "calls/{sample}/{read_model}/{depth}x/{arm}.caller_info.tsv",
    log:
        LOGS / "call_dorado/{sample}.{read_model}.{depth}x.{arm}.log",
    benchmark:
        RESULTS / "benchmarks/call_dorado/{sample}.{read_model}.{depth}x.{arm}.tsv"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=30,
    params:
        bin=lambda wc: dorado_bin(wc.arm),
        model_flag=lambda wc: f"--{ARMS[wc.arm]['calling_model']}",
        models_dir=DORADO_MODELS_DIR,
        device=DORADO["device"],
        min_depth=DORADO["min_depth"],
    shell:
        """
        outdir=$(dirname {output.vcf})
        {params.bin} polish {input.bam} {input.mutref} {params.model_flag} --vcf \
            --min-depth {params.min_depth} --any-bam --ignore-read-groups \
            --models-directory {params.models_dir} --threads {threads} \
            --device {params.device} -o "$outdir" -v > "$outdir/stdout.txt" 2> {log}

        version=$({params.bin} --version 2>&1 | tail -n 1)
        model=$(sed -n 's/.*Resolved model from input data: \\([^ ]*\\).*/\\1/p' {log} | tail -n 1)
        [ -n "$model" ] || {{ echo "could not find the resolved model in {log}" >&2; exit 1; }}
        weights_sha256=$(sha256sum {params.models_dir}/"$model"/weights.pt | cut -d' ' -f1)
        if [ "{params.device}" = cpu ]; then
            hardware="cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ *//')"
        else
            hardware=$(nvidia-smi --query-gpu=name,driver_version --format=csv,noheader | sed -n 1p)
        fi
        mkdir -p "$(dirname {output.info})"
        printf 'caller\\tcaller_version\\tcalling_model\\tcalling_model_weights_sha256\\tdevice\\thardware\\thost\\n' > {output.info}
        printf 'dorado\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "$version" "$model" "$weights_sha256" \
            "{params.device}" "$hardware" "$(hostname)" >> {output.info}
        """


# Timing-only CPU run of the same Dorado command (#12), on the Depths in dorado.cpu_run.
# It isn't scored: compare_dorado_devices checks its calls against call_dorado's instead.
CPU_RUN = DORADO.get("cpu_run") or {}
CALL_CPU = WORK / "call_cpu/{sample}/{read_model}/{depth}x/{arm}"


use rule call_dorado as call_dorado_cpu with:
    output:
        vcf=CALL_CPU / "variants.vcf",
        info=RESULTS / "calls/{sample}/{read_model}/{depth}x/{arm}.cpu.caller_info.tsv",
    log:
        LOGS / "call_dorado_cpu/{sample}.{read_model}.{depth}x.{arm}.log",
    benchmark:
        RESULTS / "benchmarks/call_dorado_cpu/{sample}.{read_model}.{depth}x.{arm}.tsv"
    threads: CPU_RUN.get("threads", 8)
    resources:
        mem_mb=16000,
        runtime=120,
    params:
        bin=lambda wc: dorado_bin(wc.arm),
        model_flag=lambda wc: f"--{ARMS[wc.arm]['calling_model']}",
        models_dir=DORADO_MODELS_DIR,
        device="cpu",
        min_depth=DORADO["min_depth"],


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
