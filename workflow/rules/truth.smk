truth_config = config["truth"]


# the create_mutref script can do this download, but doing it manually means troubleshooting is faster
rule download_genomes:
    output:
        outdir=directory(RESULTS / "truth/genomes/{sample}/"),
    log:
        LOGS / "download_genomes/{sample}.log",
    threads: 32
    resources:
        mem_mb=4 * GB,
        runtime="4h",
    conda:
        ENVS / "download_genomes.yaml"
    shadow:
        "shallow"
    params:
        taxid=infer_taxid,
        max_asm=truth_config["max_assemblies"],
        taxonomy="ncbi",
        opts="-d refseq -g bacteria -f genomic.fna.gz -m -a",
    shell:
        """
        genome_updater.sh {params.opts} -o {output.outdir} -M {params.taxonomy} \
            -T {params.taxid} -A species:{params.max_asm} -m -a -t {threads} -l "" 2>{log}
        """


rule create_mutref:
    input:
        genome=rules.faidx_reference.output.fasta,
        script=SCRIPTS / "create_mutref.py",
        genomes=rules.download_genomes.output.outdir,
    output:
        mutref=RESULTS / "truth/{sample}/mutreference.fna",
        truth_vcf=RESULTS / "truth/{sample}/truth.vcf.gz",
        stats=RESULTS / "truth/{sample}/vcfstats.txt",
        mutdonor=RESULTS / "truth/{sample}/mutdonor.fna",
    log:
        LOGS / "create_mutref/{sample}.log",
    resources:
        mem_mb=lambda wildcards, attempt: 4 * GB * attempt,
        runtime="1h",
    threads: 8
    conda:
        ENVS / "create_mutref.yaml"
    params:
        taxid=infer_taxid,
        ani=truth_config["ANI"],
        max_ani=truth_config["max_ANI"],
        min_ani=truth_config["min_ANI"],
        max_asm=truth_config["max_assemblies"],
        max_contam=get_max_contamination,
        min_completeness=get_min_completeness,
        min_comp_perc=get_min_completeness_percentile,
        outdir=lambda wildcards, output: Path(output.mutref).parent,
        max_indel=truth_config["max_indel"],
        asm_lvl='""',
        taxonomy="ncbi",
        flags="--remove-overlaps -vv --force",  # -O see https://github.com/mbhall88/NanoVarBench/issues/5
    shadow:
        "shallow"
    shell:
        """
        python {input.script} \
            {params.flags} \
            -s {params.taxid} \
            -a {params.ani} \
            -M {params.max_ani} \
            -m {params.min_ani} \
            -C {params.min_completeness} \
            -P {params.min_comp_perc} \
            -X {params.max_contam} \
            -t {threads} \
            -o {params.outdir} \
            -I {params.max_indel} \
            -l {params.asm_lvl} \
            -T {params.taxonomy} \
            -A {params.max_asm} \
            -f {input.genomes} \
            {input.genome} 2> {log}
        """


rule faidx_mutref:
    input:
        reference=rules.create_mutref.output.mutref,
    output:
        faidx=RESULTS / "truth/{sample}/mutreference.fna.fai",
    log:
        LOGS / "faidx_mutref/{sample}.log",
    resources:
        mem_mb=500,
        runtime="5m",
    container:
        "docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0"
    shell:
        "samtools faidx {input.reference} 2> {log}"


rule mutref_summary:
    input:
        logs=expand(LOGS / "create_mutref/{sample}.log", sample=SAMPLES),
        mutrefs=expand(RESULTS / "truth/{sample}/mutreference.fna", sample=SAMPLES),
    output:
        csv=TABLES / "mutref_summary.csv",
        latex=TABLES / "mutref_summary.tex",
    log:
        LOGS / "mutref_summary.log",
    resources:
        mem_mb=500,
        runtime="5m",
    conda:
        ENVS / "mutref_summary.yaml"
    params:
        samples2species={sample: infer_species(sample) for sample in SAMPLES},
    script:
        SCRIPTS / "mutref_summary.py"


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
    shadow:
        "shallow"
    shell:
        """
        python {input.script} \
            {params.flags} \
            -t {threads} \
            -o {params.outdir} \
            --ref {input.reference} \
            --donor {input.donor} 2> {log}
        """


rule make_full_bed:
    input:
        faidx=rules.faidx_mutref.output.faidx,
    output:
        bed=RESULTS / "truth/{sample}/{sample}.bed",
    log:
        LOGS / "make_full_bed/{sample}.log",
    resources:
        mem_mb=500,
        runtime="1m",
    shell:
        """
        awk '{{print $1"\t"0"\t"$2}}' {input.faidx} | sort -k1,1 -k2,2n > {output.bed} 2> {log}
        """


rule identify_repetitive_regions:
    input:
        reference=rules.create_mutref.output.mutref,
        faidx=rules.faidx_mutref.output.faidx,
    output:
        repetitive_bed=RESULTS / "truth/{sample}/{sample}.repetitive_regions.bed",
        unique_bed=RESULTS / "truth/{sample}/{sample}.unique_regions.bed",
    log:
        LOGS / "identify_repetitive_regions/{sample}.log",
    resources:
        mem_mb=500,
        runtime="5m",
    threads: 2
    conda:
        ENVS / "identify_repetitive_regions.yaml"
    shadow:
        "shallow"
    params:
        min_exact_match=20,
        min_aln_identity=60,
    shell:
        """
        exec 2>{log}
        delta=$(mktemp -u).delta
        nucmer --maxmatch --nosimplify --delta "$delta" \
            -l {params.min_exact_match} {input.reference} {input.reference}
        show-coords -rTH -I {params.min_aln_identity} "$delta" | 
            awk '{{if ($1 != $3 && $2 != $4) print $0}}' | 
            awk '{{print $8"\t"$1"\t"$2}}' |
            sort -k1,1 -k2,2n |
            bedtools merge -i - > {output.repetitive_bed}
        total_bases=$(grep -v '^>' {input.reference} | tr -d '\n' | wc -c)
        total_repeats=$(awk '{{sum+=$3-$2}} END {{print sum}}' {output.repetitive_bed})
        echo "Total bases: $total_bases" >> {log}
        echo "Total repeats: $total_repeats" >> {log}
        perc_repeats=$(printf %.2f $(echo "$total_repeats / $total_bases * 100" | bc -l))
        echo "Percentage of genome that is repetitive: $perc_repeats%" >> {log}
        bedtools complement -i {output.repetitive_bed} -g {input.faidx} > {output.unique_bed}
        """
