#!/usr/bin/env bash
set -euo pipefail

exec 2>"${snakemake_log[0]}" # send all stderr from this script to the log file

reads="${snakemake_input[reads]}"
ref="${snakemake_input[reference]}"
outvcf="${snakemake_output[vcf]}"
model="${snakemake_wildcards[model]}"
version="${snakemake_wildcards[version]}"

# Example model is dna_r10.4.1_e8.2_400bps_fast
model="${model#dna_}"          # Remove the "dna_" prefix
model="${model//./}"           # Remove all periods

# Check if the model ends with "fast" and replace it with "hac"
if [[ "$model" == *fast ]]; then
    model="${model%fast}hac"   # Replace "fast" with "hac"
fi

# Construct the model_name
model_name="${model}_variant_${version}"

tempdir=$(mktemp -d)

medaka_variant -f -m "$model_name" -i "$reads" -r "$ref" -o "$tempdir" \
    -t "${snakemake[threads]}"

bcftools view -o "$outvcf" "$tempdir/"medaka.annotated.vcf