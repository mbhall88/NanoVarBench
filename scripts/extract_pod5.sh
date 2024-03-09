#!/usr/bin/env bash
set -euo pipefail


# create the mapping file that maps the reads ids to their expected output files
# mapping="pod5_mapping.csv"
# echo "target,read_id" > "$mapping"

# loop over every fastq file in the directory
# extract the sample name from the file name
# fq=$1
# sample=$(basename $fq .fq.gz)
# read_ids=$2
# mapping="$sample.mapping.csv"
# extract the read ids from the fastq file and create a two column CSV where the first column
# is the target pod5 name, which is the sample name with a batch number that increases 
# every 4000 reads and the second column is the read id
# seqkit seq -in "$fq" | awk -v sample="$sample" 'BEGIN {i=0} {if (i%4000==0) {batch++} print sample "_batch" batch ".pod5," $1; i++}' > "$mapping"
# echo "target,read_id" > "$sample.mapping.csv"
# rapidgzip -dc $fq | 
    # python /data/scratch/projects/punim2009/NanoVarBench/scripts/extract_read_ids.py $read_ids |
    # awk -v sample="$sample" 'BEGIN {i=0} {if (i%4000==0) {batch++} print sample "_batch" batch ".pod5," $1; i++}' >> "$sample.mapping.csv"

# print suvvessful message to stderr
# echo "Extracted read ids from $fq" >&2

POD5_DIR="$1"
MAPPING="$2"
OUTDIR="$3"

if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi

pod5 subset -o "$OUTDIR" -r -f -t 8 --csv "$MAPPING" "$POD5_DIR"