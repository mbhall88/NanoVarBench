rule QC:
    input:
        reads=infer_original_reads_path,
    output:
        reads=RESULTS / "QC/subsample/{version}/{model}/{sample}.fq.gz",
    log:
        LOGS / "QC/{version}/{model}/{sample}.log",
    resources:
        mem_mb=GB,
        runtime="1h",
    benchmark:
        BENCH / "QC/{version}/{model}/{sample}.tsv",
    conda:
        ENVS / "qc.yaml",  # todo: write env file