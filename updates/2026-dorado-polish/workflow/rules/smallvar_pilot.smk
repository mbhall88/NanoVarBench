# Optional pilot (#29): `dorado smallvar`, Dorado's diploid small-variant caller, on a Dorado
# Arm's alignment (the RG-reheadered minimap2 BAM, ADR-0004). Dorado 2.1.2 has smallvar models
# only for hac v5.2.0 and v6.0.0 reads: none for v4.3.0 (our reads) and none for sup. The
# configured model (hac v6.0.0) is forced in with --model-override, which skips Dorado's
# compatibility check, so every result is a Basecall-model mismatch (on sup reads a double
# one: tier and version).
#
#   the Arm (results.tsv): RG BAM -> dorado polish --bacteria --vcf          -> Filter chain -> vcfdist
#   this pilot:            RG BAM -> dorado smallvar, whole genome hemizygous -> Filter chain -> vcfdist
#
# The whole genome (genome.bed) is passed as --hemizygous-regions, so smallvar makes haploid
# calls. Its diploid mode can't take the AF filter (#20): smallvar's VCF has only GT and GQ in
# FORMAT, and its INFO/DP header line is never filled in, so there is no AF or allele depth to
# resolve a het by.
#
# This is not an Arm: it adds nothing to results.tsv or benchmarks.tsv. Its tables are
# <results_dir>/tables/smallvar_pilot{,_runtime}.tsv. `snakemake ... smallvar_pilot` runs
# just this.

SV = config.get("smallvar_pilot") or {}
SV_ARMS = SV.get("arms") or []
SV_MODEL = SV.get("model") or ""
SV_COMPARE_ARMS = SV.get("compare_arms") or []
SV_AF = SV.get("af_filter") or {}
SV_AF_ARM = SV_AF.get("arm")
SV_AF_THRESHOLD = f"{float(SV_AF['threshold']):.2f}" if SV_AF.get("threshold") is not None else None
# The Read sets' BAMs, Mutated references, Truth sets and the compared Arms' scores are read
# from input_work_dir, and the AF filter's scores from af_filter.work_dir, so finished runs can
# be reused without re-making them. Empty = this run's work_dir, where the workflow makes them.
SV_IN = Path(SV.get("input_work_dir") or config["work_dir"])
SV_AF_IN = Path(SV_AF.get("work_dir") or config["work_dir"])
SV_THREADS = 8

if SV_ARMS:
    for a in SV_ARMS:
        if a not in ARMS or ARMS[a]["caller"] != "dorado":
            raise ValueError(f"smallvar_pilot.arms: {a} must be a Dorado Arm")
    if not SV_MODEL:
        raise ValueError("smallvar_pilot.model is empty")
    for a in SV_COMPARE_ARMS:
        if a not in ARMS:
            raise ValueError(f"smallvar_pilot.compare_arms: {a} is not an Arm")
    if SV_AF_ARM and (SV_AF_ARM not in ARMS or ARMS[SV_AF_ARM]["caller"] != "clair3"):
        raise ValueError(f"smallvar_pilot.af_filter.arm: {SV_AF_ARM} must be a Clair3 Arm")
    if SV_AF_ARM and SV_AF_THRESHOLD is None:
        raise ValueError("smallvar_pilot.af_filter.threshold is empty")


def tier_version(model):
    """(tier, version) of a basecall or smallvar model name, e.g. ('hac', 'v6.0.0') for
    dna_r10.4.1_e8.2_400bps_hac@v6.0.0_smallvar@v1.0 (the first <tier>@<version> in it)."""
    m = re.search(r"_(fast|hac|sup)@(v[\d.]+)", model)
    if not m:
        raise ValueError(f"Can't read a basecall model tier and version from {model}")
    return m.groups()


def basecall_model_mismatch(read_model):
    """How the smallvar model's basecall model differs from the Read model's: version (hac
    v4.3.0 reads, hac v6.0.0 model), tier_and_version (sup v4.3.0 reads) or none."""
    reads, model = tier_version(config["read_models"][read_model]), tier_version(SV_MODEL)
    diffs = [name for name, r, m in zip(("tier", "version"), reads, model) if r != m]
    return "_and_".join(diffs) or "none"


SV_IN_ALIGN = SV_IN / "align/{sample}/{read_model}/{depth}x"
SV_IN_TRUTH = SV_IN / "truth/{sample}"
SV_IN_SCORE = SV_IN / "score/{sample}/{read_model}/{depth}x/{arm}/{mode}"
SV_AF_SCORE = SV_AF_IN / "score_clair3_af/{sample}/{read_model}/{depth}x/{arm}/af{af}/{mode}"
# The pilot's calls and scores. They have their own directories so they can't be mistaken for
# an Arm's.
SV_CALL = WORK / "call_dorado_smallvar/{sample}/{read_model}/{depth}x/{arm}"
SV_SCORE = WORK / "score_dorado_smallvar/{sample}/{read_model}/{depth}x/{arm}/{mode}"
SV_OUT = RESULTS / "smallvar_pilot/calls/{sample}/{read_model}/{depth}x"
BENCH_SMALLVAR = RESULTS / "benchmarks/call_dorado_smallvar/{sample}.{read_model}.{depth}x.{arm}.tsv"


rule call_dorado_smallvar:
    """dorado smallvar on the Arm's RG-reheadered BAM, with the smallvar model forced in
    (--model-override) and the whole genome hemizygous, so the calls are haploid. --any-bam
    (hidden, as for polish) accepts a minimap2 BAM, and --min-depth 2 matches the Arms
    (ADR-0002). Only the tool runs here: the job is timed."""
    input:
        bam=lambda wc: SV_IN_ALIGN / f"{arm_alignment(wc.arm)}.rg.bam",
        bai=lambda wc: SV_IN_ALIGN / f"{arm_alignment(wc.arm)}.rg.bam.bai",
        mutref=SV_IN_TRUTH / "mutreference.fna",
        faidx=SV_IN_TRUTH / "mutreference.fna.fai",
        hemizygous=SV_IN_TRUTH / "genome.bed",
        model=DORADO_MODELS_DIR / f"{SV_MODEL}/weights.pt",
    output:
        vcf=SV_CALL / "variants.vcf",
    wildcard_constraints:
        arm=arm_pattern("dorado"),
    log:
        LOGS / "call_dorado_smallvar/{sample}.{read_model}.{depth}x.{arm}.log",
    benchmark:
        BENCH_SMALLVAR
    threads: SV_THREADS
    resources:
        mem_mb=16000,
        runtime=30,
    params:
        bin=lambda wc: dorado_bin(wc.arm),
        model=DORADO_MODELS_DIR / SV_MODEL,
        device=DORADO["device"],
        min_depth=DORADO["min_depth"],
    shell:
        """
        outdir=$(dirname {output.vcf})
        {params.bin} smallvar {input.bam} {input.mutref} --model-override '{params.model}' \
            --hemizygous-regions {input.hemizygous} --min-depth {params.min_depth} \
            --any-bam --ignore-read-groups --threads {threads} --device {params.device} \
            -o "$outdir" -v > "$outdir/stdout.txt" 2> {log}
        """


use rule filter_calls as filter_calls_smallvar with:
    input:
        vcf=rules.call_dorado_smallvar.output.vcf,
        reference=SV_IN_TRUTH / "mutreference.fna",
        faidx=SV_IN_TRUTH / "mutreference.fna.fai",
        filter_script=SCRIPTS / "filter_hets.py",
    output:
        vcf=SV_OUT / "{arm}.smallvar.filter.vcf.gz",
        csi=SV_OUT / "{arm}.smallvar.filter.vcf.gz.csi",
    log:
        LOGS / "filter_calls_smallvar/{sample}.{read_model}.{depth}x.{arm}.log",
    wildcard_constraints:
        arm=arm_pattern("dorado"),
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10,


# The same vcfdist settings as the Arms, resources included: -t and -r come from them.
use rule vcfdist as vcfdist_smallvar with:
    input:
        query=rules.filter_calls_smallvar.output.vcf,
        truth=SV_IN_TRUTH / "truth.vcf.gz",
        truth_index=SV_IN_TRUTH / "truth.vcf.gz.csi",
        mutref=SV_IN_TRUTH / "mutreference.fna",
        faidx=SV_IN_TRUTH / "mutreference.fna.fai",
        bed=SV_IN_TRUTH / "genome.bed",
    output:
        summary=SV_SCORE / "precision-recall-summary.tsv",
        pr=SV_SCORE / "precision-recall.tsv",
        truth=SV_SCORE / "truth.tsv",
        query=SV_SCORE / "query.tsv",
    log:
        LOGS / "vcfdist_smallvar/{sample}.{read_model}.{depth}x.{arm}.{mode}.log",
    wildcard_constraints:
        arm=arm_pattern("dorado"),
    benchmark:
        RESULTS / "benchmarks/vcfdist_smallvar/{sample}.{read_model}.{depth}x.{arm}.{mode}.tsv"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=30,


def sv_calling_model(arm, read_model):
    if ARMS[arm]["caller"] == "dorado":
        return dorado_model(arm)
    return CLAIR3["calling_models"][read_model]


def sv_score_sets():
    """Every set of scores the table compares, per Read set: smallvar on each pilot Arm's
    alignment, each compared Arm's own calls (the scores results.tsv has) and the AF filter's
    calls at its threshold."""
    sets = []
    for s in RUN_SAMPLES:
        for rm in RUN_READ_MODELS:
            for d in RUN_DEPTHS:
                keys = dict(sample=s, read_model=rm, depth=str(d))
                summaries = lambda stem, **kw: {
                    f"{mode}_summary": f"{str(stem).format(**keys, **kw, mode=mode)}/precision-recall-summary.tsv"
                    for mode in SCORING_MODES
                }
                for a in SV_ARMS:
                    sets.append(
                        dict(
                            keys,
                            arm=a,
                            calls="smallvar_haploid",
                            af_threshold="",
                            calling_model=SV_MODEL,
                            basecall_model_mismatch=basecall_model_mismatch(rm),
                            **summaries(SV_SCORE, arm=a),
                        )
                    )
                for a in SV_COMPARE_ARMS:
                    sets.append(
                        dict(
                            keys,
                            arm=a,
                            calls="main",
                            af_threshold="",
                            calling_model=sv_calling_model(a, rm),
                            basecall_model_mismatch="none",
                            **summaries(SV_IN_SCORE, arm=a),
                        )
                    )
                if SV_AF_ARM:
                    sets.append(
                        dict(
                            keys,
                            arm=SV_AF_ARM,
                            calls="af_filter",
                            af_threshold=SV_AF_THRESHOLD,
                            calling_model=sv_calling_model(SV_AF_ARM, rm),
                            basecall_model_mismatch="none",
                            **summaries(SV_AF_SCORE, arm=SV_AF_ARM, af=SV_AF_THRESHOLD),
                        )
                    )
    return sets


def sv_runs():
    """The smallvar jobs whose benchmark files go into the runtime table."""
    return [
        dict(
            sample=s,
            read_model=rm,
            depth=str(d),
            arm=a,
            benchmark=str(BENCH_SMALLVAR).format(sample=s, read_model=rm, depth=d, arm=a),
        )
        for s in RUN_SAMPLES
        for rm in RUN_READ_MODELS
        for d in RUN_DEPTHS
        for a in SV_ARMS
    ]


rule smallvar_pilot_tables:
    """Best F1 and the Default-PASS score of smallvar next to the compared Arms and the AF
    filter, per Read set, and smallvar's wall time and memory. Every smallvar row carries its
    Basecall-model mismatch."""
    input:
        files=[s[f"{m}_summary"] for s in sv_score_sets() for m in SCORING_MODES],
        benchmarks=[r["benchmark"] for r in sv_runs()],
        samples=config["samples"],
    output:
        tsv=RESULTS / "tables/smallvar_pilot.tsv",
        runtime=RESULTS / "tables/smallvar_pilot_runtime.tsv",
    log:
        LOGS / "smallvar_pilot_tables.log",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
    params:
        score_sets=sv_score_sets(),
        runs=sv_runs(),
        model=SV_MODEL,
        mismatch={rm: basecall_model_mismatch(rm) for rm in RUN_READ_MODELS} if SV_ARMS else {},
        read_basecall_models=config["read_models"],
        dorado_version={a: str(ARMS[a]["caller_version"]) for a in SV_ARMS},
        device=DORADO_DEVICE,
        hardware=HARDWARE[DORADO_DEVICE],
        threads=SV_THREADS,
    script:
        "../scripts/smallvar_pilot.py"


localrules:
    smallvar_pilot_tables,


def smallvar_pilot_targets():
    if not SV_ARMS:
        return []
    return [rules.smallvar_pilot_tables.output.tsv, rules.smallvar_pilot_tables.output.runtime]


rule smallvar_pilot:
    """Run the smallvar pilot (#29) without the rest of the workflow:
    snakemake ... smallvar_pilot"""
    input:
        smallvar_pilot_targets(),
    localrule: True
