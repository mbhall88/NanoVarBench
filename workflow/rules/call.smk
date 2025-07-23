REPEAT = config.get("repeat", 1)

rule call_mutref_illumina:
    input:
        r1=rules.preprocess_illumina.output.r1,
        r2=rules.preprocess_illumina.output.r2,
        reference=rules.faidx_mutref.input.reference,
    output:
        vcf=RESULTS / "call/mutref/illumina/{sample}/{sample}.raw.vcf.gz",
        alignment=RESULTS / "call/mutref/illumina/{sample}/{sample}.bam",
    log:
        LOGS / "call_mutref_illumina/{sample}.log",
    benchmark:
        repeat(
            BENCH / "call_mutref_illumina/{sample}.tsv",
            REPEAT,
        )
    threads: 4
    resources:
        mem_mb=12 * GB,
        runtime="4h",
    params:
        opts="--force --prefix {sample} --mapqual 0 --mincov 2 --minqual 0",
        outdir=lambda wildcards, output: Path(output.vcf).parent,
        raw_vcf=lambda wildcards, output: Path(output.vcf).with_suffix(""),
    container:
        "docker://quay.io/biocontainers/snippy:4.6.0--hdfd78af_3"
    shell:
        """
        snippy {params.opts} --cpus {threads} --reference {input.reference} \
            --R1 {input.r1} --R2 {input.r2} --outdir {params.outdir} 2> {log}
        bcftools view -o {output.vcf} {params.raw_vcf} 2>> {log}
        bcftools index -f {output.vcf} 2>> {log}
        """


caller = "bcftools"


rule call_mutref_bcftools:
    input:
        alignment=rules.align_to_mutref.output.alignment,
        reference=rules.align_to_mutref.input.reference,
        faidx=rules.faidx_mutref.output.faidx,
    output:
        vcf=RESULTS
        / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.{{depth}}x.{caller}.vcf.gz",
    log:
        LOGS
        / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.log",
    benchmark:
        repeat(
            BENCH
            / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.tsv",
            REPEAT,
        )
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 4**attempt * GB,
        runtime="6h",
    container:
        "docker://quay.io/biocontainers/bcftools:1.19--h8b25389_0"
    shadow:
        "shallow"
    script:
        "../scripts/callers/bcftools.sh"


caller = "clair3"


rule call_mutref_clair3:
    input:
        alignment=rules.align_to_mutref.output.alignment,
        reference=rules.align_to_mutref.input.reference,
        faidx=rules.faidx_mutref.output.faidx,
    output:
        vcf=RESULTS
        / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.{{depth}}x.{caller}.vcf.gz",
    log:
        LOGS
        / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.log",
    benchmark:
        repeat(
            BENCH
            / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.tsv",
            REPEAT,
        )
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 4**attempt * GB,
        runtime="6h",
    container:
        "docker://quay.io/mbhall88/clair3:1.0.10"
    shadow:
        "shallow"
    script:
        "../scripts/callers/clair3.sh"


caller = "deepvariant"


rule call_mutref_deepvariant:
    input:
        alignment=rules.align_to_mutref.output.alignment,
        reference=rules.align_to_mutref.input.reference,
        faidx=rules.faidx_mutref.output.faidx,
    output:
        vcf=RESULTS
        / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.{{depth}}x.{caller}.vcf.gz",
    log:
        LOGS
        / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.log",
    benchmark:
        repeat(
            BENCH
            / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.tsv",
            REPEAT,
        )
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 4**attempt * GB,
        runtime="2d",
    container:
        "docker://google/deepvariant:1.9.0"
    shadow:
        "shallow"
    script:
        "../scripts/callers/deepvariant.sh"


caller = "freebayes"


rule call_mutref_freebayes:
    input:
        alignment=rules.align_to_mutref.output.alignment,
        reference=rules.align_to_mutref.input.reference,
        faidx=rules.faidx_mutref.output.faidx,
    output:
        vcf=RESULTS
        / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.{{depth}}x.{caller}.vcf.gz",
    log:
        LOGS
        / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.log",
    benchmark:
        repeat(
            BENCH
            / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.tsv",
            REPEAT,
        )
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16 * GB,
        runtime="3d",
    conda:
        ENVS / f"{caller}.yaml"
    shadow:
        "shallow"
    script:
        "../scripts/callers/freebayes.sh"


caller = "longshot"


rule call_mutref_longshot:
    input:
        alignment=rules.align_to_mutref.output.alignment,
        reference=rules.align_to_mutref.input.reference,
        faidx=rules.faidx_mutref.output.faidx,
    output:
        vcf=RESULTS
        / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.{{depth}}x.{caller}.vcf.gz",
    log:
        LOGS
        / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.log",
    benchmark:
        repeat(
            BENCH
            / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.tsv",
            REPEAT,
        )
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 4**attempt * GB,
        runtime="1d",
    conda:
        ENVS / f"{caller}.yaml"
    shadow:
        "shallow"
    script:
        "../scripts/callers/longshot.sh"


caller = "medaka"


rule call_mutref_medaka:
    input:
        reads=rules.align_to_mutref.input.reads,
        reference=rules.align_to_mutref.input.reference,
        faidx=rules.faidx_mutref.output.faidx,
    output:
        vcf=RESULTS
        / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.{{depth}}x.{caller}.vcf.gz",
    log:
        LOGS
        / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.log",
    benchmark:
        repeat(
            BENCH
            / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.tsv",
            REPEAT,
        )
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 4**attempt * GB,
        runtime="1d",
    container:
        "docker://quay.io/biocontainers/medaka:2.0.1--py310he807b20_0"
    shadow:
        "shallow"
    script:
        "../scripts/callers/medaka.sh"


caller = "nanocaller"


rule call_mutref_nanocaller:
    input:
        alignment=rules.align_to_mutref.output.alignment,
        reference=rules.align_to_mutref.input.reference,
        faidx=rules.faidx_mutref.output.faidx,
    output:
        vcf=RESULTS
        / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.{{depth}}x.{caller}.vcf.gz",
    log:
        LOGS
        / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.log",
    benchmark:
        repeat(
            BENCH
            / f"call/mutref/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.tsv",
            REPEAT,
        )
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 4**attempt * GB,
        runtime="6h",
    container:
        "docker://genomicslab/nanocaller:3.4.1"
    shadow:
        "shallow"
    script:
        "../scripts/callers/nanocaller.sh"


rule filter_variants:
    input:
        vcf=RESULTS
        / "call/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.vcf.gz",
        reference=rules.faidx_mutref.input.reference,
        faidx=rules.faidx_mutref.output.faidx,
        filter_script=SCRIPTS / "filter_hets.py",
    output:
        vcf=RESULTS
        / "call/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.filter.vcf.gz",
    log:
        LOGS
        / "filter_variants/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.log",
    resources:
        mem_mb=int(0.5 * GB),
        runtime="3m",
    conda:
        ENVS / "filter_variants.yaml"
    params:
        max_indel=config["truth"].get("max_indel", 50),
    shell:
        """
        exec 2> {log}
        contigs=$(mktemp -u).contigs.txt
        header=$(mktemp -u).header.txt
        trap 'rm -f $contigs $header' EXIT 

        # bcftools reheader only adds contigs that appear in the VCF, we want all contigs
        awk '{{print "##contig=<ID="$1",length="$2">"}}' {input.faidx} > "$contigs"  # make contig lines with all contigs
        (bcftools view -h {input.vcf} |                                            # output VCF header
            grep -v "^##contig=" |                                                 # remove contig lines
            sed -e "3r $contigs") > "$header"                                      # add contig lines after 3rd line

        (bcftools reheader -h "$header" {input.vcf} |                       # replace VCF header with new header containing all contigs
            python {input.filter_script} |                                  # make heterozygous calls homozygous for allele with most depth
            bcftools view -i 'GT="alt"' |                                   # remove non-alt alleles 
            bcftools view -e 'ALT="."' |                                    # remove sites with no alt allele (NanoCaller bug)
            bcftools norm -f {input.reference} -a -c e -m - |               # normalise and left-align indels
            bcftools norm -aD |                                             # remove duplicates after normalisation
            #bcftools norm -D -N |                                            # remove duplicates without normalisation
            bcftools filter -e 'abs(ILEN)>{params.max_indel} || ALT="*"' |  # remove long indels or sites with unobserved alleles
            bcftools +setGT - -- -t a -n c:M |                              # make genotypes haploid e.g., 1/1 -> 1
            bcftools sort |                                                 # sort VCF
            bcftools view -i 'GT="A"' -o {output.vcf})                      # remove non-alt alleles and write index
        
        # make sure the VCF is not empty
        counts=$(bcftools +counts {output.vcf})

        echo "$counts" 1>&2

        # extract the total number of sites
        num_sites=$(echo "$counts" | grep "Number of sites" | awk '{{print $4}}')

        # Check if the number of sites is 0
        if [[ "$num_sites" -eq 0 ]]; then
            echo "Error: Number of sites is 0" 1>&2
            exit 1
        fi
        
        bcftools index -f {output.vcf}
        """


use rule filter_variants as filter_variants_illumina with:
    input:
        vcf=rules.call_mutref_illumina.output.vcf,
        reference=rules.faidx_mutref.input.reference,
        faidx=rules.faidx_mutref.output.faidx,
        filter_script=SCRIPTS / "filter_hets.py",
    output:
        vcf=RESULTS / "call/{ref}/illumina/{sample}/{sample}.filter.vcf.gz",
    log:
        LOGS / "filter_variants_illumina/{ref}/{sample}.log",
