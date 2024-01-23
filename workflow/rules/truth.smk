rule create_mutref:
    input:
        genome=infer_reference_genome,
        script=SCRIPTS / "create_mutref.py",
    output:
        mutref=RESULTS / "truth/{sample}/mutreference.fna",
        truth_vcf=RESULTS / "truth/{sample}/truth.vcf.gz",
        stats=RESULTS / "truth/{sample}/vcfstats.txt",
    log:
        LOGS / "create_mutref/{sample}.log",
    resources:
        mem_mb=2 * GB,
        runtime="6h",
    threads: 8
    conda:
        ENVS / "create_mutref.yaml"
    params:
        taxid=infer_taxid,
        distance=config["mash_distance"],
        outdir=lambda wildcards, output: Path(output.mutref).parent,
        max_indel=config["max_indel"],
        asm_lvl='""',
        taxonomy="ncbi",
        flags="-v --force",
    shadow:
        "shallow"
    shell:
        """
        python {input.script} \
            {params.flags} \
            -s {params.taxid} \
            -d {params.distance} \
            -t {threads} \
            -o {params.outdir} \
            -I {params.max_indel} \
            -l {params.asm_lvl} \
            -T {params.taxonomy} \
            {input.genome} 2> {log}
        """
