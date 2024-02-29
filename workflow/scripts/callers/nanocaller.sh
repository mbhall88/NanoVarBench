#!/usr/bin/env bash
set -euo pipefail

# send all stderr and stdout from this script to the log file
exec &>"${snakemake_log[0]}"

aln="${snakemake_input[alignment]}"
ref="${snakemake_input[reference]}"
finalvcf="${snakemake_output[vcf]}"
sample="${snakemake_wildcards[sample]}"
threads="${snakemake[threads]}"
tmpdir=$(mktemp -d)

NanoCaller --sample "$sample" --prefix "$sample" --cpu "$threads" --haploid_genome \
    --bam "$aln" --ref "$ref" --preset ont --output "$tmpdir" --mincov 2

mv "${tmpdir}/${sample}.vcf.gz" "$finalvcf"
