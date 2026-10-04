# Optional extra analysis (#20): the AF filter. Clair3 runs without --haploid_precise, so it
# reports its diploid genotypes, and each het call is then made homozygous for its ALT when
# that ALT's allele frequency (FORMAT/AF) is at least the AF threshold, and homozygous REF
# otherwise. The calls then go through the Filter chain and vcfdist unchanged, for every AF
# threshold in the config.
#
#   the Arm (results.tsv): Clair3 --haploid_precise          -> Filter chain -> vcfdist
#   this analysis:         Clair3 (diploid) -> AF threshold  -> Filter chain -> vcfdist
#
# Why diploid rather than --haploid_sensitive: both recover the repeat SNPs --haploid_precise
# drops (#6), but --haploid_sensitive writes 0/1 and 1/1 alike as "1" and drops multi-allelic
# (1/2) sites outright (Clair3's CallVariants.py). Diploid output keeps Clair3's own het/hom
# call, so the AF threshold only resolves the hets and leaves the hom calls alone.
#
# This is not an Arm: it adds nothing to results.tsv. Its tables are
# <results_dir>/tables/clair3_af_filter{,_summary}.tsv, plus the PR curves of the one series
# the figures add (clair3_af_filter_pr_curves.tsv, #28). `snakemake ... clair3_af_filter`
# runs just this.

AF = config.get("clair3_af_filter") or {}
AF_ARMS = AF.get("arms") or []
AF_THRESHOLDS = [f"{float(t):.2f}" for t in AF.get("thresholds") or []]
AF_REFERENCE_ARM = AF.get("reference_arm")
# The Read sets' BAMs, Mutated references, Truth sets and the Arms' own scores are read from
# input_work_dir, so a finished run's can be reused without re-making them. Empty = this run's
# work_dir, where the workflow makes them.
AF_IN = Path(AF.get("input_work_dir") or config["work_dir"])

if AF_ARMS:
    for a in AF_ARMS:
        if a not in ARMS or ARMS[a]["caller"] != "clair3":
            raise ValueError(f"clair3_af_filter.arms: {a} must be a Clair3 Arm")
    if not AF_THRESHOLDS:
        raise ValueError("clair3_af_filter.thresholds is empty")
    for t in AF_THRESHOLDS:
        if not 0 < float(t) <= 1:
            raise ValueError(f"clair3_af_filter.thresholds: {t} is not in (0, 1]")
    if AF_REFERENCE_ARM is not None and AF_REFERENCE_ARM not in ARMS:
        raise ValueError(f"clair3_af_filter.reference_arm {AF_REFERENCE_ARM} is not an Arm")



def figure_af_filter():
    """The AF filter series the figures and Table S1 add (figures.af_filter, #28): None when
    the analysis isn't run for that Arm, else its Arm and AF threshold, the threshold written
    as the tables write it."""
    spec = config["figures"].get("af_filter")
    if not spec or spec["arm"] not in AF_ARMS:
        return None
    threshold = f"{float(spec['threshold']):.2f}"
    if threshold not in AF_THRESHOLDS:
        raise ValueError(
            f"figures.af_filter.threshold {threshold} is not in clair3_af_filter.thresholds"
        )
    return {"arm": spec["arm"], "threshold": threshold}


FIGURE_AF = figure_af_filter()

# The Arms' Clair3 options without the haploid mode, so Clair3 writes 0/1, 1/1 and 1/2.
CLAIR3_AF_OPTIONS = [
    o for o in CLAIR3["options"] if o not in ("--haploid_precise", "--haploid_sensitive")
]

AF_IN_ALIGN = AF_IN / "align/{sample}/{read_model}/{depth}x"
AF_IN_TRUTH = AF_IN / "truth/{sample}"
AF_IN_SCORE = AF_IN / "score/{sample}/{read_model}/{depth}x/{arm}/{mode}"
# This analysis' calls and scores. They have their own directories so they can't be mistaken
# for an Arm's.
AF_CALL = WORK / "call_clair3_af/{sample}/{read_model}/{depth}x/{arm}"
AF_SCORE = WORK / "score_clair3_af/{sample}/{read_model}/{depth}x/{arm}/af{af}/{mode}"
AF_OUT = RESULTS / "clair3_af_filter/calls/{sample}/{read_model}/{depth}x"


# The Arm's Clair3 command on the Arm's alignment, with only the haploid mode dropped. It isn't
# timed for Table 1, so the Bunya profile doesn't pin it to the timed CPU nodes.
use rule call_clair3 as call_clair3_af with:
    input:
        bam=lambda wc: AF_IN_ALIGN / f"{arm_alignment(wc.arm)}.bam",
        bai=lambda wc: AF_IN_ALIGN / f"{arm_alignment(wc.arm)}.bam.bai",
        mutref=AF_IN_TRUTH / "mutreference.fna",
        faidx=AF_IN_TRUTH / "mutreference.fna.fai",
        model_files=clair3_model_files,
    output:
        vcf=AF_CALL / "variants.vcf.gz",
    log:
        LOGS / "call_clair3_af/{sample}.{read_model}.{depth}x.{arm}.log",
    benchmark:
        RESULTS / "benchmarks/call_clair3_af/{sample}.{read_model}.{depth}x.{arm}.tsv"
    threads: CLAIR3["threads"]
    resources:
        mem_mb=16000,
        runtime=120,
    params:
        model_path=clair3_model_path,
        options=" ".join(CLAIR3_AF_OPTIONS),


rule af_filter:
    """Resolve the diploid calls' hets at one AF threshold: ALT if FORMAT/AF >= it, else REF.
    Hom calls are left alone. The Filter chain then runs on the result as on any Arm's calls."""
    input:
        vcf=rules.call_clair3_af.output.vcf,
        script=SCRIPTS / "af_filter.py",
    output:
        vcf=AF_CALL / "af{af}.vcf.gz",
    wildcard_constraints:
        arm=arm_pattern("clair3"),
        af=r"\d\.\d\d",
    log:
        LOGS / "af_filter/{sample}.{read_model}.{depth}x.{arm}.af{af}.log",
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10,
    conda:
        ENVS / "filter.yaml"
    shell:
        "python {input.script} --min-af {wildcards.af} {input.vcf} -o {output.vcf} 2> {log}"


use rule filter_calls as filter_calls_af with:
    input:
        vcf=rules.af_filter.output.vcf,
        reference=AF_IN_TRUTH / "mutreference.fna",
        faidx=AF_IN_TRUTH / "mutreference.fna.fai",
        filter_script=SCRIPTS / "filter_hets.py",
    output:
        vcf=AF_OUT / "{arm}.af{af}.filter.vcf.gz",
        csi=AF_OUT / "{arm}.af{af}.filter.vcf.gz.csi",
    log:
        LOGS / "filter_calls_af/{sample}.{read_model}.{depth}x.{arm}.af{af}.log",
    wildcard_constraints:
        arm=arm_pattern("clair3"),
        af=r"\d\.\d\d",
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10,


# The same vcfdist settings as the Arms, resources included: -t and -r come from them.
use rule vcfdist as vcfdist_af with:
    input:
        query=rules.filter_calls_af.output.vcf,
        truth=AF_IN_TRUTH / "truth.vcf.gz",
        truth_index=AF_IN_TRUTH / "truth.vcf.gz.csi",
        mutref=AF_IN_TRUTH / "mutreference.fna",
        faidx=AF_IN_TRUTH / "mutreference.fna.fai",
        bed=AF_IN_TRUTH / "genome.bed",
    output:
        summary=AF_SCORE / "precision-recall-summary.tsv",
        pr=AF_SCORE / "precision-recall.tsv",
        truth=AF_SCORE / "truth.tsv",
        query=AF_SCORE / "query.tsv",
    log:
        LOGS / "vcfdist_af/{sample}.{read_model}.{depth}x.{arm}.af{af}.{mode}.log",
    wildcard_constraints:
        arm=arm_pattern("clair3"),
        af=r"\d\.\d\d",
    benchmark:
        RESULTS / "benchmarks/vcfdist_af/{sample}.{read_model}.{depth}x.{arm}.af{af}.{mode}.tsv"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=30,


def af_score_sets():
    """Every set of scores the tables compare, per Read set: each AF Arm's own calls
    (--haploid_precise, the scores results.tsv has), the same Arm at every AF threshold, and
    the reference Arm's calls (e.g. Dorado)."""
    sets = []
    for s in RUN_SAMPLES:
        for rm in RUN_READ_MODELS:
            for d in RUN_DEPTHS:
                keys = dict(sample=s, read_model=rm, depth=str(d))
                main = lambda arm, mode: str(AF_IN_SCORE).format(**keys, arm=arm, mode=mode)
                for a in AF_ARMS:
                    sets.append(
                        dict(
                            keys,
                            arm=a,
                            calls="main",
                            af_threshold="",
                            sweep_summary=f"{main(a, 'sweep')}/precision-recall-summary.tsv",
                            pass_summary=f"{main(a, 'pass')}/precision-recall-summary.tsv",
                        )
                    )
                    for t in AF_THRESHOLDS:
                        af = lambda mode: str(AF_SCORE).format(**keys, arm=a, af=t, mode=mode)
                        sets.append(
                            dict(
                                keys,
                                arm=a,
                                calls="af_filter",
                                af_threshold=t,
                                sweep_summary=f"{af('sweep')}/precision-recall-summary.tsv",
                                pass_summary=f"{af('pass')}/precision-recall-summary.tsv",
                            )
                        )
                if AF_REFERENCE_ARM and AF_REFERENCE_ARM not in AF_ARMS:
                    a = AF_REFERENCE_ARM
                    sets.append(
                        dict(
                            keys,
                            arm=a,
                            calls="main",
                            af_threshold="",
                            sweep_summary=f"{main(a, 'sweep')}/precision-recall-summary.tsv",
                            pass_summary=f"{main(a, 'pass')}/precision-recall-summary.tsv",
                        )
                    )
    return sets


rule clair3_af_filter_tables:
    """Best F1 and the Default-PASS score of the AF filter at every AF threshold, next to the
    Arm's own --haploid_precise calls and the reference Arm, per Sample and summarised."""
    input:
        files=[s[k] for s in af_score_sets() for k in ("sweep_summary", "pass_summary")],
        samples=config["samples"],
    output:
        tsv=RESULTS / "tables/clair3_af_filter.tsv",
        summary=RESULTS / "tables/clair3_af_filter_summary.tsv",
    log:
        LOGS / "clair3_af_filter_tables.log",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
    params:
        score_sets=af_score_sets(),
        clair3_options=" ".join(CLAIR3_AF_OPTIONS),
        reference_arm=AF_REFERENCE_ARM or "",
    script:
        "../scripts/clair3_af_filter.py"


localrules:
    clair3_af_filter_tables,


def af_curve_files():
    """The QUAL sweep's precision-recall file of the figures' AF filter series, per Read set."""
    if not FIGURE_AF:
        return []
    return [
        str(AF_SCORE).format(
            sample=s, read_model=rm, depth=d, arm=FIGURE_AF["arm"], af=FIGURE_AF["threshold"],
            mode="sweep",
        )
        + "/precision-recall.tsv"
        for s in RUN_SAMPLES
        for rm in RUN_READ_MODELS
        for d in RUN_DEPTHS
    ]


rule clair3_af_filter_pr_curves:
    """The QUAL sweep's precision-recall curves of the AF filter series that Figure 2 draws
    (figures.af_filter), in pr_curves.tsv's layout. The curves of the other thresholds aren't
    tabulated; the Arms' own are in pr_curves.tsv."""
    input:
        files=af_curve_files(),
    output:
        tsv=RESULTS / "tables/clair3_af_filter_pr_curves.tsv",
    log:
        LOGS / "clair3_af_filter_pr_curves.log",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
    params:
        read_sets=[
            dict(sample=s, read_model=rm, depth=d)
            for s in RUN_SAMPLES
            for rm in RUN_READ_MODELS
            for d in RUN_DEPTHS
        ],
        arm=FIGURE_AF["arm"] if FIGURE_AF else "",
        af_threshold=FIGURE_AF["threshold"] if FIGURE_AF else "",
    script:
        "../scripts/clair3_af_filter_pr_curves.py"


localrules:
    clair3_af_filter_pr_curves,


def af_filter_targets():
    if not AF_ARMS:
        return []
    targets = [rules.clair3_af_filter_tables.output.tsv, rules.clair3_af_filter_tables.output.summary]
    if FIGURE_AF:
        targets.append(rules.clair3_af_filter_pr_curves.output.tsv)
    return targets


rule clair3_af_filter:
    """Run the AF filter analysis (#20) without the rest of the workflow:
    snakemake ... clair3_af_filter"""
    input:
        af_filter_targets(),
    localrule: True
