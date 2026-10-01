# Optional check (#13) that the reheadered-minimap2 path (ADR-0004) doesn't change Dorado's
# calls. For each Read set under `dorado_aligner_check` in the config (empty = off), the Read
# set is aligned a second way, with `dorado aligner` (its default preset is lr:hq, with its
# bundled minimap2). The same @RG line is added, and `dorado polish` runs WITHOUT --any-bam,
# because the BAM has the @PG line from dorado aligner that polish looks for. Both paths'
# calls then go through the Filter chain and vcfdist, and compare_aligner_paths reports how
# they differ.
#
#   path 1 (the Arm's normal run): minimap2 2.31 lr:hq -> @RG -> polish --any-bam
#   path 2 (this check):           dorado aligner      -> @RG -> polish
#
# This is not an Arm: it adds nothing to results.tsv, and its outputs are in
# <results_dir>/dorado_aligner_check/.

CHECK = config.get("dorado_aligner_check") or {}
CHECK_SAMPLES = CHECK.get("samples") or []
CHECK_READ_MODELS = CHECK.get("read_models") or []
CHECK_DEPTHS = [int(d) for d in CHECK.get("depths") or []]
CHECK_ARM = CHECK.get("arm", "D")  # the Dorado Arm whose binary and Calling model are used
CHECK_MATERIAL = CHECK.get("material", {"max_records": 5, "max_best_f1_diff": 0.0005})

if CHECK_SAMPLES or CHECK_READ_MODELS or CHECK_DEPTHS:
    if CHECK_ARM not in ARMS or ARMS[CHECK_ARM]["caller"] != "dorado":
        raise ValueError(f"dorado_aligner_check.arm {CHECK_ARM} must be a Dorado Arm")
    for s in CHECK_SAMPLES:
        if s not in SAMPLES:
            raise ValueError(f"dorado_aligner_check: Sample {s} is not in {config['samples']}")
        for rm in CHECK_READ_MODELS:
            if (s, rm) not in RUNS:
                raise ValueError(f"dorado_aligner_check: no {rm} run for Sample {s}")

# Path 2's BAMs, VCFs and scores. They have their own directories so they can't be mistaken
# for an Arm's (the Arms' BAMs are in ALIGN).
CHECK_ALIGN = WORK / "align_dorado_aligner/{sample}/{read_model}/{depth}x"
CHECK_CALL = WORK / "call_dorado_aligner/{sample}/{read_model}/{depth}x/{arm}"
CHECK_SCORE = WORK / "score_dorado_aligner/{sample}/{read_model}/{depth}x/{arm}/{mode}"
CHECK_OUT = RESULTS / "dorado_aligner_check/{sample}/{read_model}/{depth}x"


rule dorado_align:
    """Align the Read set with dorado aligner (default preset lr:hq, bundled minimap2). It
    writes an unsorted BAM to stdout, with the @PG line (ID:aligner) polish looks for."""
    input:
        fastq=rules.read_set.output.fastq,
        mutref=rules.extract_truth.output.mutref,
    output:
        bam=temp(CHECK_ALIGN / "dorado_aligner.unsorted.bam"),
    log:
        LOGS / "dorado_align/{sample}.{read_model}.{depth}x.log",
    benchmark:
        RESULTS / "benchmarks/dorado_align/{sample}.{read_model}.{depth}x.tsv"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=120,
    params:
        bin=lambda wc: dorado_bin(CHECK_ARM),
    shell:
        """
        {params.bin} aligner {input.mutref} {input.fastq} --threads {threads} \
            > {output.bam} 2> {log}
        """


rule dorado_align_sort:
    """Sort and index dorado aligner's BAM. samtools sort keeps the @PG lines, so this is
    all polish needs: it fails here if the aligner's @PG line is lost."""
    input:
        bam=rules.dorado_align.output.bam,
    output:
        bam=CHECK_ALIGN / "dorado_aligner.bam",
        bai=CHECK_ALIGN / "dorado_aligner.bam.bai",
    log:
        LOGS / "dorado_align_sort/{sample}.{read_model}.{depth}x.log",
    threads: 4
    resources:
        mem_mb=8000,
        runtime=60,
    container:
        CONTAINERS["samtools"]
    shell:
        """
        exec 2> {log}
        samtools sort -@ {threads} -m 1G -T {output.bam}.tmp -o {output.bam} {input.bam}
        samtools index {output.bam}
        samtools view -H {output.bam} | grep '^@PG' >&2
        samtools view -H {output.bam} | grep -q '^@PG.*ID:aligner' \
            || {{ echo "the @PG line from dorado aligner is missing" >&2; exit 1; }}
        """


# The same single @RG line as path 1 (rg_reheader), so the two BAMs differ only in how the
# reads were aligned.
use rule rg_reheader as dorado_aligner_reheader with:
    input:
        bam=CHECK_ALIGN / "dorado_aligner.bam",
    output:
        bam=CHECK_ALIGN / "dorado_aligner.rg.bam",
        bai=CHECK_ALIGN / "dorado_aligner.rg.bam.bai",
    log:
        LOGS / "dorado_aligner_reheader/{sample}.{read_model}.{depth}x.log",


# dorado polish with the same flags as call_dorado (--min-depth 2, --ignore-read-groups, the
# same Calling model and device), but no --any-bam. It needs a GPU like call_dorado does:
# the Bunya profile gives it the same entry.
use rule call_dorado as call_dorado_aligner_bam with:
    input:
        bam=CHECK_ALIGN / "dorado_aligner.rg.bam",
        bai=CHECK_ALIGN / "dorado_aligner.rg.bam.bai",
        mutref=rules.extract_truth.output.mutref,
        faidx=rules.index_mutref.output.faidx,
        model=lambda wc: DORADO_MODELS_DIR / dorado_model(wc.arm) / "weights.pt",
    output:
        vcf=CHECK_CALL / "variants.vcf",
        info=CHECK_OUT / "{arm}.dorado_aligner.caller_info.tsv",
    log:
        LOGS / "call_dorado_aligner_bam/{sample}.{read_model}.{depth}x.{arm}.log",
    benchmark:
        RESULTS / "benchmarks/call_dorado_aligner_bam/{sample}.{read_model}.{depth}x.{arm}.tsv"
    params:
        bin=lambda wc: dorado_bin(wc.arm),
        model_flag=lambda wc: f"--{ARMS[wc.arm]['calling_model']}",
        models_dir=DORADO_MODELS_DIR,
        device=DORADO["device"],
        min_depth=DORADO["min_depth"],
        any_bam="",


use rule filter_calls as filter_calls_dorado_aligner with:
    input:
        vcf=CHECK_CALL / "variants.vcf",
        reference=rules.extract_truth.output.mutref,
        faidx=rules.index_mutref.output.faidx,
        filter_script=SCRIPTS / "filter_hets.py",
    output:
        vcf=CHECK_OUT / "{arm}.dorado_aligner.filter.vcf.gz",
        csi=CHECK_OUT / "{arm}.dorado_aligner.filter.vcf.gz.csi",
    log:
        LOGS / "filter_calls_dorado_aligner/{sample}.{read_model}.{depth}x.{arm}.log",


use rule vcfdist as vcfdist_dorado_aligner with:
    input:
        query=rules.filter_calls_dorado_aligner.output.vcf,
        truth=rules.extract_truth.output.truth,
        truth_index=rules.index_truth.output.index,
        mutref=rules.extract_truth.output.mutref,
        faidx=rules.index_mutref.output.faidx,
        bed=rules.index_mutref.output.bed,
    output:
        summary=CHECK_SCORE / "precision-recall-summary.tsv",
        pr=CHECK_SCORE / "precision-recall.tsv",
        truth=CHECK_SCORE / "truth.tsv",
        query=CHECK_SCORE / "query.tsv",
    log:
        LOGS / "vcfdist_dorado_aligner/{sample}.{read_model}.{depth}x.{arm}.{mode}.log",
    benchmark:
        RESULTS / "benchmarks/vcfdist_dorado_aligner/{sample}.{read_model}.{depth}x.{arm}.{mode}.tsv"


def check_keys(wc):
    return dict(sample=wc.sample, read_model=wc.read_model, depth=wc.depth, arm=wc.arm)


rule compare_aligner_paths:
    """Compare path 1 (Arm's minimap2 BAM, --any-bam) with path 2 (dorado aligner BAM, no
    --any-bam) on the same Read set: records that differ in Dorado's raw and filtered VCFs,
    the largest QUAL difference, and the vcfdist Best F1 and Default-PASS F1 differences."""
    input:
        raw1=lambda wc: str(CALL / "variants.vcf").format(**check_keys(wc)),
        raw2=rules.call_dorado_aligner_bam.output.vcf,
        filter1=lambda wc: str(rules.filter_calls.output.vcf).format(**check_keys(wc)),
        filter2=rules.filter_calls_dorado_aligner.output.vcf,
        info1=lambda wc: str(CALLER_INFO).format(**check_keys(wc)),
        info2=rules.call_dorado_aligner_bam.output.info,
        sweep1=lambda wc: str(rules.vcfdist.output.summary).format(**check_keys(wc), mode="sweep"),
        pass1=lambda wc: str(rules.vcfdist.output.summary).format(**check_keys(wc), mode="pass"),
        sweep2=lambda wc: str(rules.vcfdist_dorado_aligner.output.summary).format(
            **check_keys(wc), mode="sweep"
        ),
        pass2=lambda wc: str(rules.vcfdist_dorado_aligner.output.summary).format(
            **check_keys(wc), mode="pass"
        ),
    output:
        tsv=CHECK_OUT / "{arm}.dorado_aligner_vs_minimap2.tsv",
        report=CHECK_OUT / "{arm}.dorado_aligner_vs_minimap2.md",
    wildcard_constraints:
        arm=arm_pattern("dorado"),
    log:
        LOGS / "compare_aligner_paths/{sample}.{read_model}.{depth}x.{arm}.log",
    resources:
        mem_mb=1000,
        runtime=10,
    params:
        material=CHECK_MATERIAL,
    script:
        "../scripts/compare_aligner_paths.py"


def aligner_check_targets():
    return [
        str(rules.compare_aligner_paths.output.report).format(
            sample=s, read_model=rm, depth=d, arm=CHECK_ARM
        )
        for s in CHECK_SAMPLES
        for rm in CHECK_READ_MODELS
        for d in CHECK_DEPTHS
    ]


rule dorado_aligner_check:
    """Run the dorado aligner check (#13) on the Read sets in `dorado_aligner_check`, without
    the rest of the workflow: snakemake ... dorado_aligner_check"""
    input:
        aligner_check_targets(),
    localrule: True
