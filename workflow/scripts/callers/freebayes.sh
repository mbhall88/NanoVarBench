#!/usr/bin/env bash
set -euo pipefail

# send all stderr and stdout from this script to the log file
exec &>"${snakemake_log[0]}"

aln="${snakemake_input[alignment]}"
ref="${snakemake_input[reference]}"
finalvcf="${snakemake_output[vcf]}"
sample="${snakemake_wildcards[sample]}"

# freebayes -f %(ref_file)s -F %(af_hard)s -r {region} --haplotype-length -1 
# %(calling_params)s %(bam_file)s | annotate_maaf.py | bcftools +fill-tags | 
# bcftools view -c 1 | bcftools norm -f %(ref_file)s -Oz -o %(prefix)s.{region_safe}.vcf.gz" % vars(self)

freebayes --haplotype-length -1 -m 10 -q 10 -p 1 --min-coverage 2 -f "$ref" "$aln" | 
    bcftools view -o "$finalvcf"