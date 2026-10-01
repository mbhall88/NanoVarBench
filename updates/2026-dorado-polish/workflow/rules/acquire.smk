# Read and Truth set acquisition. Downloads land in reads_dir and truth_dir and are only
# fetched when missing. Every download is checked against the expected MD5 (and FASTQs
# with gzip -t) by this workflow, never by the downloader: kingfisher 0.5.0's own
# --check-md5sums always reports "MD5sum OK" (#6).


rule download_reads:
    output:
        fastq=READS_DIR / "{read_model}/{run}_1.fastq.gz",
    log:
        LOGS / "download_reads/{read_model}/{run}.log",
    params:
        url=lambda wc: RUN_URLS[wc.run]["fastq_url"],
        md5=lambda wc: RUN_URLS[wc.run]["fastq_md5"],
        method=config["download"]["method"],
        kingfisher=config["download"].get("kingfisher", "kingfisher"),
        attempts=config["download"].get("attempts", 3),
    resources:
        mem_mb=4000,
        runtime=12 * 60,
    shell:
        """
        exec 2> {log}
        out={output.fastq}
        verify() {{ echo "{params.md5}  $out" | md5sum -c - && gzip -t "$out"; }}
        for attempt in $(seq 1 {params.attempts}); do
            rm -f "$out"
            if [ "{params.method}" = kingfisher-ascp ]; then
                {params.kingfisher} get -r {wildcards.run} -m ena-ascp --output-directory "$(dirname "$out")" || true
            else
                curl -fsSL --retry 3 -o "$out" "{params.url}" || true
            fi
            if [ -s "$out" ] && verify; then exit 0; fi
            echo "attempt $attempt failed verification" >&2
            sleep 30
        done
        rm -f "$out"
        exit 1
        """


rule download_truth:
    output:
        tarball=TRUTH_DIR / "{sample}.tar.gz",
    log:
        LOGS / "download_truth/{sample}.log",
    params:
        url=lambda wc: SAMPLES[wc.sample]["truth_url"],
        md5=lambda wc: SAMPLES[wc.sample]["truth_md5"],
    resources:
        mem_mb=1000,
        runtime=30,
    shell:
        """
        exec 2> {log}
        curl -fsSL --retry 3 -o {output.tarball} "{params.url}"
        echo "{params.md5}  {output.tarball}" | md5sum -c -
        """


rule extract_truth:
    """Unpack the Mutated reference and Truth set from the Zenodo 10867171 tarball."""
    input:
        tarball=rules.download_truth.output.tarball,
    output:
        mutref=WORK / "truth/{sample}/mutreference.fna",
        truth=WORK / "truth/{sample}/truth.vcf.gz",
    log:
        LOGS / "extract_truth/{sample}.log",
    resources:
        mem_mb=1000,
        runtime=10,
    shell:
        """
        exec 2> {log}
        tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
        tar -xzf {input.tarball} -C "$tmp" {wildcards.sample}/mutreference.fna {wildcards.sample}/truth.vcf.gz
        mv "$tmp/{wildcards.sample}/mutreference.fna" {output.mutref}
        mv "$tmp/{wildcards.sample}/truth.vcf.gz" {output.truth}
        """


rule index_mutref:
    input:
        mutref=rules.extract_truth.output.mutref,
    output:
        faidx=WORK / "truth/{sample}/mutreference.fna.fai",
        bed=WORK / "truth/{sample}/genome.bed",
    log:
        LOGS / "index_mutref/{sample}.log",
    resources:
        mem_mb=1000,
        runtime=10,
    container:
        CONTAINERS["samtools"]
    shell:
        """
        exec 2> {log}
        samtools faidx {input.mutref}
        awk 'BEGIN{{OFS="\\t"}} {{print $1, 0, $2}}' {output.faidx} | sort -k1,1 -k2,2n > {output.bed}
        """


rule index_truth:
    input:
        truth=rules.extract_truth.output.truth,
    output:
        index=WORK / "truth/{sample}/truth.vcf.gz.csi",
    log:
        LOGS / "index_truth/{sample}.log",
    resources:
        mem_mb=1000,
        runtime=10,
    container:
        CONTAINERS["bcftools"]
    shell:
        "bcftools index -f {input.truth} 2> {log}"
