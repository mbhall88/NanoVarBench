#!/usr/bin/env bash
set -euo pipefail

# send all stderr and stdout from this script to the log file
exec &>"${snakemake_log[0]}"

aln="${snakemake_input[alignment]}"
ref="${snakemake_input[reference]}"
outvcf="${snakemake_output[vcf]}"
model="${snakemake_wildcards[model]}"
version="${snakemake_wildcards[version]}"
# need to map model name e.g. dna_r10.4.1_e8.2_400bps_sup to clair model name e.g. r1041_e82_400bps_sup_v430
# this is a bit hacky but it works
model_name=$(echo "$model" | sed -E 's/.*dna_(.*)@.*/\1/')
model_name=$(echo "$model_name" | sed -E 's/\.//g')
model_name="${model_name}_${version}"
model_name=$(echo "$model_name" | sed -E 's/\.//g')

# if the model contains _fast@ then we need to use the hac model
if [[ "$model" == *_fast@* ]]; then
    model_name=$(echo "$model_name" | sed -E 's/_fast/_hac/')
fi

model_path="/opt/models/${model_name}"
tmpoutdir=$(mktemp -d)
sample="${snakemake_wildcards[sample]}"

trap 'rm -rf "$tmpoutdir"' EXIT

run_clair3.sh \
    --bam_fn="$aln" \
    --ref_fn="$ref" \
    --threads="${snakemake[threads]}" \
    --platform="ont" \
    --model_path="$model_path" \
    --output="$tmpoutdir" \
    --sample_name="$sample" \
    --include_all_ctgs \
    --haploid_precise \
    --no_phasing_for_fa \
    --enable_long_indel

mv "${tmpoutdir}/merge_output.vcf.gz" "$outvcf"
