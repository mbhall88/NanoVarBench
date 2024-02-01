#!/usr/bin/env bash
set -euo pipefail

# send all stderr and stdout from this script to the log file
exec &>"${snakemake_log[0]}"

aln="${snakemake_input[alignment]}"
ref="${snakemake_input[reference]}"
finalvcf="${snakemake_output[vcf]}"
sample="${snakemake_wildcards[sample]}"
outvcf=$(mktemp -u).vcf

longshot -AFn \
    -b "$aln" \
    -f "$ref" \
    -o "$outvcf" \
    -I 50 \
    -s "$sample"

bcftools view -e 'GT="het"' -o "$finalvcf" "$outvcf"