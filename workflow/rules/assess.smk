rule assess_self_calls:
    input:
        vcf=RESULTS
        / "call/self/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.filter.vcf.gz",
        reference=rules.align_to_self.input.reference,
    output:
        csv=RESULTS
        / "assess/self/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.csv",
        json=RESULTS
        / "assess/self/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.json",
    log:
        LOGS
        / "assess_self_calls/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.log",
    resources:
        runtime="2m",
        mem_mb=500,
    conda:
        ENVS / "assess_self_calls.yaml"
    script:
        SCRIPTS / "assess_self_calls.py"


rule combine_self_calls:
    input:
        csvs=expand(
            RESULTS
            / "assess/self/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.csv",
            caller=CALLERS,
            depth=[MAX_DEPTH],
            mode=MODES,
            version=VERSIONS,
            model=MODELS,
            sample=SAMPLES,
        ),
    output:
        csv=RESULTS / "assess/self/self_calls.csv",
    log:
        log=LOGS / "combine_self_calls.log",
    resources:
        runtime="2m",
        mem_mb=500,
    container:
        "docker://python:3.11-slim"
    script:
        SCRIPTS / "combine_csvs.py"


rule plot_self_fpr:
    input:
        csv=rules.combine_self_calls.output.csv,
    output:
        snps_png=FIGURES / "assess/self/self_fpr.snps.png",
        indels_png=FIGURES / "assess/self/self_fpr.indels.png",
        csv=TABLES / "assess/self/self_fpr.csv",
    log:
        LOGS / "plot_self_fpr.log",
    resources:
        runtime="5m",
        mem_mb=500,
    params:
        sample2species={sample: infer_species(sample) for sample in SAMPLES},
        no_indels=config.get("no_indels", []),
    conda:
        ENVS / "plot_self_fpr.yaml"
    script:
        SCRIPTS / "plot_self_fpr.py"


rule assess_mutref_calls:
    input:
        query_vcf=RESULTS
        / "call/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.filter.vcf.gz",
        truth_vcf=rules.create_mutref.output.truth_vcf,
        mutref=rules.create_mutref.output.mutref,
        faidx=rules.faidx_mutref.output.faidx,
    output:
        pr_summary=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.precision-recall-summary.tsv",
        pr=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.precision-recall.tsv",
        summary=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.summary.vcf",
        query=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.query.tsv",
        truth=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.truth.tsv",
        dist_summary=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.distance-summary.tsv",
        dist=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.distance.tsv",
    log:
        LOGS
        / "assess_mutref_calls/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.log",
    resources:
        runtime="5m",
        mem_mb=500,
    container:
        "docker://timd1/vcfdist:v2.3.3"
    params:
        opts=f"--largest-variant {config['truth']['max_indel']} --credit-threshold 1.0 -d",
        prefix=lambda wildcards, output: Path(output.pr).with_suffix("").with_suffix(""),
    shell:
        """
        exec 2> {log}
        echo "Calculated maximum QUAL score..." 1>&2
        MAX_QUAL=$(bgzip -dc {input.query_vcf} | grep -v '^#' | cut -f 6 | sort -gr | sed -n '1p')
        echo "MAX_QUAL=$MAX_QUAL" 1>&2
        tmpbed=$(mktemp -u).bed
        echo "Writing BED file that covers the whole reference..." 1>&2
        awk '{{print $1"\t"0"\t"$2}}' {input.faidx} > $tmpbed
        cat $tmpbed 1>&2
        echo "Running vcfdist..." 1>&2
        vcfdist {input.query_vcf} {input.truth_vcf} {input.mutref} {params.opts} \
            -p {params.prefix}. -b "$tmpbed" -mx $MAX_QUAL
        """
