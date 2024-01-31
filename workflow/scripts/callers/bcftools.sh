#!/usr/bin/env bash
set -euo pipefail

exec 2>"${snakemake_log[0]}" # send all stderr from this script to the log file

aln="${snakemake_input[alignment]}"
ref="${snakemake_input[reference]}"
outvcf="${snakemake_output[vcf]}"

bcftools mpileup -f "$ref" \
    -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    -Ou --threads "${snakemake[threads]}" -Q 10 -x -M 10000 -h 100 "$aln" |
    bcftools call --ploidy 1 --threads "${snakemake[threads]}" -m --prior 0.005 \
        -o "$outvcf" --variants-only
