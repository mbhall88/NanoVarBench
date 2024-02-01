#!/usr/bin/env bash
set -euo pipefail

exec 2>"${snakemake_log[0]}" # send all stderr from this script to the log file

aln="${snakemake_input[alignment]}"
ref="${snakemake_input[reference]}"
outvcf="${snakemake_output[vcf]}"
sample="${snakemake_wildcards[sample]}"
tmpvcf=$(mktemp -u).vcf.gz

run_deepvariant --model_type ONT_R104 \
    --output_vcf "$tmpvcf" \
    --num_shards "${snakemake[threads]}" \
    --reads "$aln" \
    --ref "$ref" \
    --sample_name "$sample" \
    --novcf_stats_report \
    --noruntime_report

bcftools view -i 'GT="AA"' -o "$outvcf" "$tmpvcf"