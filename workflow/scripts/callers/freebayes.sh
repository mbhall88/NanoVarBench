#!/usr/bin/env bash
set -euo pipefail

# send all stderr and stdout from this script to the log file
exec &>"${snakemake_log[0]}"

aln="${snakemake_input[alignment]}"
ref="${snakemake_input[reference]}"
finalvcf="${snakemake_output[vcf]}"
sample="${snakemake_wildcards[sample]}"
chunk_size=500000
threads="${snakemake[threads]}"

freebayes-parallel <(fasta_generate_regions.py "$ref" "$chunk_size") "$threads" \
    --haplotype-length -1 -m 10 -q 10 -p 1 --min-coverage 2 -f "$ref" "$aln" |
    bcftools view -o "$finalvcf"
