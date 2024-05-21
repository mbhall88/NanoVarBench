#!/usr/bin/env bash
set -euo pipefail

exec 2>"${snakemake_log[0]}" # send all stderr from this script to the log file

reads="${snakemake_input[reads]}"
ref="${snakemake_input[reference]}"
outvcf="${snakemake_output[vcf]}"
model="${snakemake_wildcards[model]}"

# if the model is sup, use r1041_e82_400bps_sup_variant_v4.3.0 otherwise r1041_e82_400bps_hac_variant_v4.3.0
if [ "$model" == "sup" ]; then
    model_name="r1041_e82_400bps_sup_variant_v4.3.0"
else
    model_name="r1041_e82_400bps_hac_variant_v4.3.0"
fi

tempdir=$(mktemp -d)

medaka_haploid_variant -m "$model_name" -i "$reads" -r "$ref" -o "$tempdir" \
    -t "${snakemake[threads]}"

bcftools view -o "$outvcf" "$tempdir/"medaka.annotated.vcf