rule stats_prefilter:
    input:
        reads=infer_original_reads_path,
    output:
        stats=RESULTS / "QC/stats/prefilter/{version}/{model}/{sample}.tsv",
    log:
        LOGS / "stats_prefilter/{version}/{model}/{sample}.log",
    resources:
        mem_mb=GB,
        runtime="1h",
    benchmark:
        BENCH / "stats_prefilter/{version}/{model}/{sample}.tsv"
    container:
        "docker://quay.io/biocontainers/seqkit:2.6.1--h9ee0642_0"
    shell:
        "seqkit stats -Ta {input.reads} > {output.stats} 2> {log}"


rule combine_stats_prefilter:
    input:
        expand(
            RESULTS / "QC/stats/prefilter/{version}/{model}/{sample}.tsv",
            version=VERSIONS,
            model=MODELS,
            sample=SAMPLES,
        ),
    output:
        stats=RESULTS / "QC/stats/prefilter/{version}/{model}/combined.tsv",
        md_stats=RESULTS / "QC/stats/prefilter/{version}/{model}/combined.md",
    log:
        LOGS / "combine_stats_prefilter/{version}/{model}/combined.log",
    resources:
        mem_mb=int(0.5 * GB),
        runtime="2m",
    container:
        "docker://quay.io/biocontainers/csvtk:0.29.0--h9ee0642_0"
    shell:
        """
        csvtk concat -tT -o {output.stats} {input} 2> {log}
        csvtk csv2md -t -o {output.md_stats} {output.stats} 2>> {log}
        """


rule filter_reads:
    input:
        reads=infer_original_reads_path,
    output:
        reads=RESULTS / "QC/filter/{version}/{model}/{sample}.fq.gz",
    log:
        LOGS / "filter_reads/{version}/{model}/{sample}.log",
    resources:
        mem_mb=GB,
        runtime="1h",
    benchmark:
        BENCH / "filter_reads/{version}/{model}/{sample}.tsv"
    container:
        "docker://quay.io/biocontainers/seqkit:2.6.1--h9ee0642_0"
    params:
        min_length=config["QC"]["min_length"],
        min_qscore=config["QC"]["min_qscore"],
    shell:
        "seqkit seq -m {params.min_length} -Q {params.min_qscore} -o {output.reads} {input.reads} 2> {log}"


use rule stats_prefilter as stats_postfilter with:
    input:
        reads=rules.filter_reads.output.reads,
    output:
        stats=RESULTS / "QC/stats/postfilter/{version}/{model}/{sample}.tsv",
    log:
        LOGS / "stats_postfilter/{version}/{model}/{sample}.log",
    benchmark:
        BENCH / "stats_postfilter/{version}/{model}/{sample}.tsv"


use rule combine_stats_prefilter as combine_stats_postfilter with:
    input:
        expand(
            RESULTS / "QC/stats/postfilter/{version}/{model}/{sample}.tsv",
            version=VERSIONS,
            model=MODELS,
            sample=SAMPLES,
        ),
    output:
        stats=RESULTS / "QC/stats/postfilter/{version}/{model}/combined.tsv",
        md_stats=RESULTS / "QC/stats/postfilter/{version}/{model}/combined.md",
    log:
        LOGS / "combine_stats_postfilter/{version}/{model}/combined.log",


rule faidx_reference:
    input:
        reference=infer_reference_genome,
    output:
        faidx=RESULTS / "reference/faidx/{sample}.fai",
    log:
        LOGS / "faidx_reference/{sample}.log",
    resources:
        mem_mb=500,
        runtime="5m",
    container:
        "docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0"
    shell:
        "samtools faidx --fai-idx {output.faidx} {input.reference} 2> {log}"


rule downsample_reads:
    input:
        reads=rules.filter_reads.output.reads,
        faidx=rules.faidx_reference.output.faidx,
    output:
        reads=RESULTS
        / "QC/downsample/{depth}x/{version}/{model}/{sample}.{depth}x.fq.gz",
    log:
        LOGS / "downsample_reads/{depth}x/{version}/{model}/{sample}.log",
    resources:
        mem_mb=GB,
        runtime="1h",
    benchmark:
        BENCH / "downsample_reads/{depth}x/{version}/{model}/{sample}.tsv"
    container:
        "docker://quay.io/mbhall88/rasusa:0.7.1"
    params:
        seed=20240102,
    shell:
        "rasusa -i {input.reads} -o {output.reads} -c {wildcards.depth} -g {input.faidx} -s {params.seed} 2> {log}"
