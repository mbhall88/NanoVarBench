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
        mem_mb=8 * GB,
        runtime="6h",
    container:
        "docker://quay.io/biocontainers/bcftools:1.19--h8b25389_0"
    shadow:
        "shallow"
    script:
        "../scripts/callers/bcftools.sh"


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
        mem_mb=8 * GB,
        runtime="6h",
    container:
        "docker://quay.io/mbhall88/clair3:1.0.5"
    shadow:
        "shallow"
    script:
        "../scripts/callers/clair3.sh"


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
        runtime="6h",
    container:
        "docker://google/deepvariant:1.6.0"
    shadow:
        "shallow"
    script:
        "../scripts/callers/deepvariant.sh"


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
        mem_mb=8 * GB,
        runtime="6h",
    container:
        "docker://google/deepvariant:1.6.0"
    shadow:
        "shallow"
    script:
        "../scripts/callers/deepvariant.sh"


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
        runtime="6h",
    conda:
        ENVS / f"{caller}.yaml"
    shadow:
        "shallow"
    script:
        "../scripts/callers/freebayes.sh"


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
    threads: 4
    resources:
        mem_mb=8 * GB,
        runtime="6h",
    conda:
        ENVS / f"{caller}.yaml"
    shadow:
        "shallow"
    script:
        "../scripts/callers/freebayes.sh"


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
        mem_mb=8 * GB,
        runtime="6h",
    conda:
        ENVS / f"{caller}.yaml"
    shadow:
        "shallow"
    script:
        "../scripts/callers/longshot.sh"
