rule faidx_reference:
    input:
        reference=infer_reference_genome,
    output:
        fasta=RESULTS / "reference/{sample}.fa",
        faidx=RESULTS / "reference/{sample}.fa.fai",
    log:
        LOGS / "faidx_reference/{sample}.log",
    resources:
        mem_mb=500,
        runtime="5m",
    container:
        "docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0"
    shell:
        """
        cp {input.reference} {output.fasta} 2> {log}
        samtools faidx {output.fasta} 2>> {log}
        """


rule stats_prefilter:
    input:
        reads=rules.basecall.output.reads,
        reference=rules.faidx_reference.output.fasta,
    output:
        stats=RESULTS / "QC/stats/prefilter/{mode}/{version}/{model}/{sample}.txt",
    log:
        LOGS / "stats_prefilter/{mode}/{version}/{model}/{sample}.log",
    threads: 4
    resources:
        mem_mb=8 * GB,
        runtime="30m",
    benchmark:
        BENCH / "stats_prefilter/{mode}/{version}/{model}/{sample}.tsv"
    conda:
        ENVS / "qc_stats.yaml"
    params:
        mm2_opts="-ax map-ont --secondary=no",
        samtools_opts="-bh -F 0x900",
        cramino_opts="--hist",
    shell:
        """
        tmpbam=$(mktemp -u --suffix .bam) 2> {log}
        (minimap2 -t {threads} {params.mm2_opts} {input.reference} {input.reads} | \
        samtools view {params.samtools_opts} -o "$tmpbam") 2>> {log}
        cramino -t {threads} {params.cramino_opts} "$tmpbam" > {output.stats} 2>> {log}
        """


rule combine_stats_prefilter:
    input:
        stats=expand(
            RESULTS / "QC/stats/prefilter/{{mode}}/{{version}}/{{model}}/{sample}.txt",
            sample=SAMPLES,
        ),
    output:
        stats=RESULTS / "QC/stats/prefilter/{mode}/{version}/{model}.csv",
    log:
        LOGS / "combine_stats_prefilter/{mode}/{version}/{model}.log",
    resources:
        runtime="1m",
    container:
        "docker://python:3.11-slim"
    script:
        SCRIPTS / "combine_stats.py"


rule filter_reads:
    input:
        reads=rules.basecall.output.reads,
    output:
        reads=temp(RESULTS / "QC/filter/{mode}/{version}/{model}/{sample}.fq.gz"),
    log:
        LOGS / "filter_reads/{mode}/{version}/{model}/{sample}.log",
    resources:
        mem_mb=GB,
        runtime="1h",
    benchmark:
        BENCH / "filter_reads/{mode}/{version}/{model}/{sample}.tsv"
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
        reference=rules.faidx_reference.output.fasta,
    output:
        stats=RESULTS / "QC/stats/postfilter/{mode}/{version}/{model}/{sample}.txt",
    log:
        LOGS / "stats_postfilter/{mode}/{version}/{model}/{sample}.log",
    benchmark:
        BENCH / "stats_postfilter/{mode}/{version}/{model}/{sample}.tsv"


use rule combine_stats_prefilter as combine_stats_postfilter with:
    input:
        stats=expand(
            RESULTS
            / "QC/stats/postfilter/{{mode}}/{{version}}/{{model}}/{sample}.txt",
            sample=SAMPLES,
        ),
    output:
        stats=RESULTS / "QC/stats/postfilter/{mode}/{version}/{model}.csv",
    log:
        LOGS / "combine_stats_postfilter/{mode}/{version}/{model}.log",


rule downsample_reads:
    input:
        reads=rules.filter_reads.output.reads,
        faidx=rules.faidx_reference.output.faidx,
    output:
        reads=RESULTS
        / "QC/downsample/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.fq.gz",
    log:
        LOGS / "downsample_reads/{depth}x/{mode}/{version}/{model}/{sample}.log",
    resources:
        mem_mb=GB,
        runtime="15m",
    benchmark:
        BENCH / "downsample_reads/{depth}x/{mode}/{version}/{model}/{sample}.tsv"
    container:
        "docker://quay.io/mbhall88/rasusa:0.8.0"
    params:
        seed=20240102,
    shell:
        "rasusa -i {input.reads} -o {output.reads} -c {wildcards.depth} -g {input.faidx} -s {params.seed} 2> {log}"


use rule stats_prefilter as stats_downsample with:
    input:
        reads=rules.downsample_reads.output.reads,
        reference=rules.faidx_reference.output.fasta,
    output:
        stats=RESULTS
        / "QC/stats/downsample/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.txt",
    log:
        LOGS / "stats_downsample/{depth}x/{mode}/{version}/{model}/{sample}.log",
    benchmark:
        BENCH / "stats_downsample/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.tsv"


use rule combine_stats_prefilter as combine_stats_downsample with:
    input:
        stats=expand(
            RESULTS
            / "QC/stats/downsample/{{depth}}x/{{mode}}/{{version}}/{{model}}/{sample}.{{depth}}x.txt",
            sample=SAMPLES,
        ),
    output:
        stats=RESULTS / "QC/stats/downsample/{depth}x/{mode}/{version}/{model}.csv",
    log:
        LOGS / "combine_stats_downsample/{depth}x/{mode}/{version}/{model}.log",


rule preprocess_illumina:
    input:
        r1=infer_illumina_reads_1,
        r2=infer_illumina_reads_2,
    output:
        r1=RESULTS / "preprocess/illumina/{sample}_1.fq.gz",
        r2=RESULTS / "preprocess/illumina/{sample}_2.fq.gz",
        json=RESULTS / "preprocess/illumina/{sample}.json",
    log:
        LOGS / "preprocess_illumina/{sample}.log",
    resources:
        mem_mb=8 * GB,
        runtime="1h",
    threads: 4
    container:
        "docker://quay.io/biocontainers/fastp:0.23.4--hadf994f_2"
    shadow:
        "shallow"
    params:
        opts="-l 30 --cut_tail --dedup --detect_adapter_for_pe",
    shell:
        """
        fastp {params.opts} -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
            -w {threads} -j {output.json} &> {log}
        """
