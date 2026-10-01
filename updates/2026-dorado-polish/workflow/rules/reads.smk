# Read QC and depth-capped subsampling into Read sets (ADR-0005).


rule qc_reads:
    """Keep reads >= min_length bp with mean quality >= min_qual (the paper's QC)."""
    input:
        fastq=raw_reads,
    output:
        fastq=WORK / "reads/{sample}/{read_model}/qc.fq.gz",
    log:
        LOGS / "qc_reads/{sample}.{read_model}.log",
    threads: 8
    resources:
        mem_mb=4000,
        runtime=120,
    params:
        min_length=config["qc"]["min_length"],
        min_qual=config["qc"]["min_qual"],
    container:
        CONTAINERS["seqkit"]
    shell:
        """
        seqkit seq -j {threads} --min-len {params.min_length} --min-qual {params.min_qual} \
            -o {output.fastq} {input.fastq} 2> {log}
        """


rule subsample_align:
    """Primary-only alignment of all QC'd reads, the input to rasusa aln."""
    input:
        fastq=rules.qc_reads.output.fastq,
        mutref=rules.extract_truth.output.mutref,
    output:
        bam=temp(WORK / "reads/{sample}/{read_model}/subsample.{aln}.primary.bam"),
        bai=temp(WORK / "reads/{sample}/{read_model}/subsample.{aln}.primary.bam.bai"),
    log:
        LOGS / "subsample_align/{sample}.{read_model}.{aln}.log",
    threads: 32
    resources:
        mem_mb=32000,
        runtime=240,
    params:
        preset=lambda wc: ALIGNMENTS[wc.aln]["preset"],
    container:
        aligner_container
    shell:
        """
        exec 2> {log}
        minimap2 -t {threads} -a -x {params.preset} --secondary=no {input.mutref} {input.fastq} \
            | samtools view -u -F 0x904 - \
            | samtools sort -@ 4 -m 2G -T {output.bam}.tmp -o {output.bam} -
        samtools index {output.bam}
        """


rule subsample_read_ids:
    """Cap per-position depth at the Depth on every contig and list the reads kept."""
    input:
        bam=WORK / f"reads/{{sample}}/{{read_model}}/subsample.{SUBSAMPLE_ALN}.primary.bam",
        bai=WORK / f"reads/{{sample}}/{{read_model}}/subsample.{SUBSAMPLE_ALN}.primary.bam.bai",
    output:
        ids=READ_SET / "read_ids.txt",
    log:
        LOGS / "subsample_read_ids/{sample}.{read_model}.{depth}x.log",
    resources:
        mem_mb=8000,
        runtime=60,
    params:
        seed=config["subsample"]["seed"],
    container:
        CONTAINERS["rasusa"]
    shell:
        """
        exec 2> {log}
        rasusa aln -c {wildcards.depth} -s {params.seed} -O sam {input.bam} \
            | grep -v '^@' | cut -f1 | sort -u > {output.ids}
        echo "reads kept: $(wc -l < {output.ids})" >&2
        """


rule read_set:
    """Pull the kept reads out of the QC'd FASTQ. This FASTQ is the Read set every Arm
    realigns."""
    input:
        fastq=rules.qc_reads.output.fastq,
        ids=rules.subsample_read_ids.output.ids,
    output:
        fastq=READ_SET / "reads.fq.gz",
    log:
        LOGS / "read_set/{sample}.{read_model}.{depth}x.log",
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    container:
        CONTAINERS["seqkit"]
    shell:
        """
        exec 2> {log}
        seqkit grep -j {threads} -f {input.ids} -o {output.fastq} {input.fastq}
        n_ids=$(wc -l < {input.ids})
        n_reads=$(seqkit stats -T {output.fastq} | awk -F'\\t' 'NR==2{{print $4}}')
        echo "ids=$n_ids reads=$n_reads" >&2
        [ "$n_ids" -eq "$n_reads" ]
        """
