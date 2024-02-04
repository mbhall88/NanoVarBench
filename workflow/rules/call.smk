REPEAT = config.get("repeat", 1)

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
        mem_mb=8 * GB,
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
    threads: 4
    resources:
        mem_mb=8 * GB,
        runtime="12h",
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


rule filter_variants:
    input:
        vcf=RESULTS
        / "call/{ref}/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.vcf.gz",
        reference=infer_vcf_reference,
        faidx=infer_vcf_reference_faidx,
    output:
        vcf=RESULTS
        / "call/{ref}/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.filter.vcf.gz",
    log:
        LOGS
        / "filter_variants/{ref}/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.log",
    resources:
        mem_mb=int(0.5 * GB),
        runtime="3m",
    container:
        "docker://quay.io/biocontainers/bcftools:1.19--h8b25389_0"
    params:
        max_indel=config["truth"].get("max_indel", 50),
    shell:
        """
        (bcftools reheader -f {input.faidx} {input.vcf} |        # add contigs to header
            bcftools view -i 'GT="alt"' |                        # remove non-alt alleles
            bcftools norm -f {input.reference} -a -c e -m - |    # normalise and left-align indels
            bcftools norm -aD |                                  # remove duplicates after normalisation
            bcftools filter -e 'abs(ILEN)>{params.max_indel}' |  # remove long indels
            sed 's|1/1|1|g' |                                    # make genotypes haploid e.g., 1/1 -> 1
            bcftools view -o {output.vcf} --write-index) 2> {log}
        """
