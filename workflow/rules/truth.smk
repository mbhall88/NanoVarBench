truth_config = config["truth"]


rule create_mutref:
    input:
        genome=infer_reference_genome,
        script=SCRIPTS / "create_mutref.py",
    output:
        mutref=RESULTS / "truth/{sample}/mutreference.fna",
        truth_vcf=RESULTS / "truth/{sample}/truth.vcf.gz",
        stats=RESULTS / "truth/{sample}/vcfstats.txt",
        mutdonor=RESULTS / "truth/{sample}/mutdonor.fna",
    log:
        LOGS / "create_mutref/{sample}.log",
    resources:
        mem_mb=8 * GB,
        runtime="12h",
    threads: 8
    conda:
        ENVS / "create_mutref.yaml"
    params:
        taxid=infer_taxid,
        distance=truth_config["mash_distance"],
        max_distance=truth_config["max_mash_distance"],
        min_distance=truth_config["min_mash_distance"],
        max_asm=truth_config["max_assemblies"],
        outdir=lambda wildcards, output: Path(output.mutref).parent,
        max_indel=truth_config["max_indel"],
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
            -M {params.max_distance} \
            -m {params.min_distance} \
            -t {threads} \
            -o {params.outdir} \
            -I {params.max_indel} \
            -l {params.asm_lvl} \
            -T {params.taxonomy} \
            -A {params.max_asm} \
            {input.genome} 2> {log}
        """

rule plot_synteny:
    input:
        reference=rules.create_mutref.input.genome,
        donor=rules.create_mutref.output.mutdonor,
        script=SCRIPTS / "plot_synteny.py",
    output:
        minimap2_plot=FIGURES / "plot_synteny/{sample}/minimap2.plotsr.png",
        nucmer_plot=FIGURES / "plot_synteny/{sample}/nucmer.plotsr.png",
    log:
        LOGS / "plot_synteny/{sample}.log",
    threads: 2
    resources:
        mem_mb=2 * GB,
        runtime="20m",
    conda:
        ENVS / "plot_synteny.yaml"
    params:
        outdir=lambda wildcards, output: Path(output.minimap2_plot).parent,
        flags="-v",
    shell:
        """
        python {input.script} \
            {params.flags} \
            -t {threads} \
            -o {params.outdir} \
            --ref {input.reference} \
            --donor {input.donor} 2> {log}
        """