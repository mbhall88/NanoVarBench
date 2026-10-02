# Each Arm's alignment of the Read set, the RG reheader Dorado needs, and actual depth.


rule align:
    """Align a Read set. Arms with the same aligner, version and preset share this BAM."""
    input:
        fastq=rules.read_set.output.fastq,
        mutref=rules.extract_truth.output.mutref,
    output:
        bam=ALIGN / "{aln}.bam",
        bai=ALIGN / "{aln}.bam.bai",
        info=ALIGN_INFO,
    log:
        LOGS / "align/{sample}.{read_model}.{depth}x.{aln}.log",
    benchmark:
        BENCH_ALIGN
    threads: 8
    resources:
        mem_mb=8000,
        runtime=60,
    params:
        preset=lambda wc: ALIGNMENTS[wc.aln]["preset"],
    container:
        aligner_container
    shell:
        """
        exec 2> {log}
        minimap2 -t {threads} -aL --cs --MD -x {params.preset} {input.mutref} {input.fastq} \
            | samtools sort -@ 2 -T {output.bam}.tmp -o {output.bam} -
        samtools index {output.bam}

        hardware="cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ *//')"
        mkdir -p "$(dirname {output.info})"
        printf 'device\\thardware\\thost\\tthreads\\n' > {output.info}
        printf 'cpu\\t%s\\t%s\\t%s\\n' "$hardware" "$(cat /proc/sys/kernel/hostname)" "{threads}" >> {output.info}
        """


rule rg_reheader:
    """Add one @RG header line carrying the Read model's basecall_model, which Dorado reads
    from the header (ADR-0004). Per-read RG:Z tags aren't needed (#7)."""
    input:
        bam=rules.align.output.bam,
    output:
        bam=ALIGN / "{aln}.rg.bam",
        bai=ALIGN / "{aln}.rg.bam.bai",
    log:
        LOGS / "rg_reheader/{sample}.{read_model}.{depth}x.{aln}.log",
    resources:
        mem_mb=2000,
        runtime=20,
    params:
        rg=lambda wc: "\\t".join(
            [
                "@RG",
                f"ID:{RUNS[(wc.sample, wc.read_model)]['run']}_{basecall_model(wc)}",
                f"SM:{wc.sample}",
                "PL:ONT",
                f"DS:basecall_model={basecall_model(wc)}",
            ]
        ),
    container:
        CONTAINERS["samtools"]
    shell:
        """
        exec 2> {log}
        header=$(mktemp); trap 'rm -f "$header"' EXIT
        samtools view -H {input.bam} | grep -v '^@RG' > "$header"
        printf '{params.rg}\\n' >> "$header"
        samtools reheader "$header" {input.bam} > {output.bam}
        samtools index {output.bam}
        samtools view -H {output.bam} | grep '^@RG' >&2
        """


rule actual_depth:
    """Per-contig depth of a Read set, measured on the actual_depth alignment (Arm B's)."""
    input:
        bam=ALIGN / f"{DEPTH_ALN}.bam",
        bai=ALIGN / f"{DEPTH_ALN}.bam.bai",
    output:
        tsv=RESULTS / "depth/{sample}.{read_model}.{depth}x.coverage.tsv",
    log:
        LOGS / "actual_depth/{sample}.{read_model}.{depth}x.log",
    resources:
        mem_mb=2000,
        runtime=20,
    container:
        CONTAINERS["samtools"]
    shell:
        "samtools coverage {input.bam} > {output.tsv} 2> {log}"
