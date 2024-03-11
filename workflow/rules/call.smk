REPEAT = config.get("repeat", 1)


rule call_self_illumina:
    input:
        r1=rules.preprocess_illumina.output.r1,
        r2=rules.preprocess_illumina.output.r2,
        reference=rules.faidx_reference.output.fasta,
    output:
        vcf=RESULTS / "call/self/illumina/{sample}/{sample}.vcf.gz",
        alignment=RESULTS / "call/self/illumina/{sample}/{sample}.bam",
    log:
        LOGS / "call_self_illumina/{sample}.log",
    benchmark:
        repeat(
            BENCH / "call_self_illumina/{sample}.tsv",
            REPEAT,
        )
    threads: 4
    resources:
        mem_mb=12 * GB,
        runtime="4h",
    params:
        opts="--force --prefix {sample}",
        outdir=lambda wildcards, output: Path(output.vcf).parent,
    container:
        "docker://quay.io/biocontainers/snippy:4.6.0--hdfd78af_3"
    shell:
        """
        snippy {params.opts} --cpus {threads} --reference {input.reference} \
            --R1 {input.r1} --R2 {input.r2} --outdir {params.outdir} 2> {log}
        """


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


rule call_self_bcftools:
    input:
        alignment=rules.align_to_self.output.alignment,
        reference=rules.align_to_self.input.reference,
        faidx=rules.faidx_reference.output.faidx,
    output:
        vcf=RESULTS
        / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.{{depth}}x.{caller}.vcf.gz",
    log:
        LOGS
        / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.log",
    benchmark:
        repeat(
            BENCH
            / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.tsv",
            REPEAT,
        )
    threads: 4
    resources:
        mem_mb=16 * GB,
        runtime="6h",
    container:
        "docker://quay.io/biocontainers/bcftools:1.19--h8b25389_0"
    shadow:
        "shallow"
    script:
        "../scripts/callers/bcftools.sh"


use rule call_self_bcftools as call_mutref_bcftools with:
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


caller = "clair3"


rule call_self_clair3:
    input:
        alignment=rules.align_to_self.output.alignment,
        reference=rules.align_to_self.input.reference,
        faidx=rules.faidx_reference.output.faidx,
    output:
        vcf=RESULTS
        / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.{{depth}}x.{caller}.vcf.gz",
    log:
        LOGS
        / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.log",
    benchmark:
        repeat(
            BENCH
            / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.tsv",
            REPEAT,
        )
    threads: 4
    resources:
        mem_mb=8 * GB,
        runtime="6h",
    container:
        "docker://quay.io/mbhall88/clair3:1.0.5"
    shadow:
        "shallow"
    script:
        "../scripts/callers/clair3.sh"


use rule call_self_clair3 as call_mutref_clair3 with:
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


caller = "deepvariant"


rule call_self_deepvariant:
    input:
        alignment=rules.align_to_self.output.alignment,
        reference=rules.align_to_self.input.reference,
        faidx=rules.faidx_reference.output.faidx,
    output:
        vcf=RESULTS
        / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.{{depth}}x.{caller}.vcf.gz",
    log:
        LOGS
        / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.log",
    benchmark:
        repeat(
            BENCH
            / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.tsv",
            REPEAT,
        )
    threads: 4
    resources:
        mem_mb=8 * GB,
        runtime="1d",
    container:
        "docker://google/deepvariant:1.6.0"
    shadow:
        "shallow"
    script:
        "../scripts/callers/deepvariant.sh"


use rule call_self_deepvariant as call_mutref_deepvariant with:
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


caller = "freebayes"


rule call_self_freebayes:
    input:
        alignment=rules.align_to_self.output.alignment,
        reference=rules.align_to_self.input.reference,
        faidx=rules.faidx_reference.output.faidx,
    output:
        vcf=RESULTS
        / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.{{depth}}x.{caller}.vcf.gz",
    log:
        LOGS
        / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.log",
    benchmark:
        repeat(
            BENCH
            / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.tsv",
            REPEAT,
        )
    threads: 8
    resources:
        mem_mb=16 * GB,
        runtime="1d",
    conda:
        ENVS / f"{caller}.yaml"
    shadow:
        "shallow"
    script:
        "../scripts/callers/freebayes.sh"


use rule call_self_freebayes as call_mutref_freebayes with:
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


caller = "longshot"


rule call_self_longshot:
    input:
        alignment=rules.align_to_self.output.alignment,
        reference=rules.align_to_self.input.reference,
        faidx=rules.faidx_reference.output.faidx,
    output:
        vcf=RESULTS
        / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.{{depth}}x.{caller}.vcf.gz",
    log:
        LOGS
        / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.log",
    benchmark:
        repeat(
            BENCH
            / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.tsv",
            REPEAT,
        )
    threads: 4
    resources:
        mem_mb=8 * GB,
        runtime="6h",
    conda:
        ENVS / f"{caller}.yaml"
    shadow:
        "shallow"
    script:
        "../scripts/callers/longshot.sh"


use rule call_self_longshot as call_mutref_longshot with:
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


caller = "nanocaller"


rule call_self_nanocaller:
    input:
        alignment=rules.align_to_self.output.alignment,
        reference=rules.align_to_self.input.reference,
        faidx=rules.faidx_reference.output.faidx,
    output:
        vcf=RESULTS
        / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.{{depth}}x.{caller}.vcf.gz",
    log:
        LOGS
        / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.log",
    benchmark:
        repeat(
            BENCH
            / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.tsv",
            REPEAT,
        )
    threads: 4
    resources:
        mem_mb=8 * GB,
        runtime="6h",
    container:
        "docker://genomicslab/nanocaller:3.4.1"
    shadow:
        "shallow"
    script:
        "../scripts/callers/nanocaller.sh"


use rule call_self_nanocaller as call_mutref_nanocaller with:
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


rule filter_variants:
    input:
        vcf=RESULTS
        / "call/{ref}/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.vcf.gz",
        reference=infer_vcf_reference,
        faidx=infer_vcf_reference_faidx,
        filter_script=SCRIPTS / "filter_hets.py",
    output:
        vcf=RESULTS
        / "call/{ref}/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.filter.vcf.gz",
    log:
        LOGS
        / "filter_variants/{ref}/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.log",
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
            bcftools filter -e 'abs(ILEN)>{params.max_indel} || ALT="*"' |  # remove long indels or sites with unobserved alleles
            bcftools +setGT - -- -t a -n c:M |                              # make genotypes haploid e.g., 1/1 -> 1
            bcftools sort |                                                 # sort VCF
            bcftools view -i 'GT="A"' -o {output.vcf})                      # remove non-alt alleles and write index
        bcftools index -f {output.vcf}
        """


use rule filter_variants as filter_variants_illumina with:
    input:
        vcf=infer_illumina_vcf_to_filter,
        reference=infer_vcf_reference,
        faidx=infer_vcf_reference_faidx,
        filter_script=SCRIPTS / "filter_hets.py",
    output:
        vcf=RESULTS / "call/{ref}/illumina/{sample}/{sample}.filter.vcf.gz",
    log:
        LOGS / "filter_variants_illumina/{ref}/{sample}.log",
