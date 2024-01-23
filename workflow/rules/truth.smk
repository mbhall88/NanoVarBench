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
        mem_mb=1 * GB,
        runtime="1h",
    threads: 4
    conda:
        ENVS / "create_mutref.yaml"
    params:
        species=infer_species,
        distance=config["mash_distance"],
        outdir=lambda wildcards, output: Path(output.mutref).parent,
        max_indel=config["max_indel"],
        flags="-v --force",
    shadow:
        "shallow"
    shell:
        """
        python {input.script} \
            {params.flags} \
            -s {params.species} \
            -d {params.distance} \
            -t {threads} \
            -o {params.outdir} \
            -I {params.max_indel} \
            {input.genome} 2> {log}
        """
