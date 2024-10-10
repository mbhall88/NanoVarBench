rule align_to_mutref:
    input:
        reads=rules.downsample_reads.output.reads,
        reference=rules.create_mutref.output.mutref,
    output:
        alignment=RESULTS
        / "align/mutref/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.bam",
        index=RESULTS
        / "align/mutref/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.bam.bai",
    log:
        LOGS / "align_to_mutref/{depth}x/{mode}/{version}/{model}/{sample}.log",
    benchmark:
        BENCH / "align_to_mutref/{depth}x/{mode}/{version}/{model}/{sample}.tsv"
    threads: 4
    resources:
        mem_mb=4 * GB,
        runtime="5m",
    conda:
        ENVS / "align.yaml"
    params:
        preset="map-ont",
        opts="-aL --cs --MD",
    shell:
        """
        (minimap2 {params.opts} -t {threads} -x {params.preset} {input.reference} {input.reads} | \
          samtools sort -@ {threads} -o {output.alignment}) 2> {log}
        samtools index {output.alignment} 2>> {log}
        """