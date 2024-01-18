#!/usr/bin/env bash
set -euo pipefail

# This script combines demuxed duplex reads for the ONT_202309 runs.

function print_usage_and_exit() {
    echo "Usage: $0 [-r run_dir] [-o outdir]" >&2
    echo "  -r run_dir: demux directory" >&2
    echo "  -o outdir: output directory" >&2
    echo "  -h: print this help message" >&2
    exit 0
}

function parse_params() {
    while getopts 'm:r:o:h' flag; do
        case "${flag}" in
        h) print_usage_and_exit ;;
        r) run_dir=$(realpath "${OPTARG}") ;;
        o) final_dir=$(realpath "${OPTARG}") ;;
        *) error "Unexpected option ${flag}" ;;
        esac
    done
}

function safe_cat() {
    if [ -f "$1" ]; then
        cat "$1"
    else
        echo "WARNING: $(realpath "$1") does not exist" >&2
    fi
}

parse_params "$@"

cd "$run_dir" || exit 1

# get a list of all unique models which are directory names under each directory in the current directory
models=$(find . -mindepth 2 -maxdepth 2 -type d -printf '%f\n' | sort -u)

# loop over models
for model in $models; do
    mkdir -p "$final_dir/$model"
    {
        safe_cat ONT_230922A/"$model"/barcode73.fastq
        safe_cat ONT_230922B/"$model"/barcode65.fastq
        safe_cat ONT_230922C/"$model"/barcode01.fastq
        safe_cat ONT_230926_HW/"$model"/barcode01.fastq
        safe_cat ONT_230922A/"$model"/barcode81.fastq
        safe_cat ONT_230922B/"$model"/barcode73.fastq
        safe_cat ONT_230922C/"$model"/barcode09.fastq
        safe_cat ONT_230926_HW/"$model"/barcode09.fastq
    } >"$final_dir/$model"/ATCC_10708__202309.fq.gz

    {
        safe_cat ONT_230922A/"$model"/barcode74.fastq
        safe_cat ONT_230922B/"$model"/barcode66.fastq
        safe_cat ONT_230922C/"$model"/barcode02.fastq
        safe_cat ONT_230926_HW/"$model"/barcode02.fastq
        safe_cat ONT_230922A/"$model"/barcode82.fastq
        safe_cat ONT_230922B/"$model"/barcode74.fastq
        safe_cat ONT_230922C/"$model"/barcode10.fastq
        safe_cat ONT_230926_HW/"$model"/barcode10.fastq
    } >"$final_dir/$model"/ATCC_14035__202309.fq.gz

    {
        safe_cat ONT_230922A/"$model"/barcode89.fastq
        safe_cat ONT_230922B/"$model"/barcode81.fastq
        safe_cat ONT_230922C/"$model"/barcode17.fastq
        safe_cat ONT_230926_HW/"$model"/barcode17.fastq
        safe_cat ONT_230922A/"$model"/barcode90.fastq
        safe_cat ONT_230922B/"$model"/barcode82.fastq
        safe_cat ONT_230922C/"$model"/barcode18.fastq
        safe_cat ONT_230926_HW/"$model"/barcode18.fastq
    } >"$final_dir/$model"/ATCC_17802__202309.fq.gz

    {
        safe_cat ONT_230922A/"$model"/barcode79.fastq
        safe_cat ONT_230922B/"$model"/barcode71.fastq
        safe_cat ONT_230922C/"$model"/barcode07.fastq
        safe_cat ONT_230926_HW/"$model"/barcode07.fastq
        safe_cat ONT_230922A/"$model"/barcode87.fastq
        safe_cat ONT_230922B/"$model"/barcode79.fastq
        safe_cat ONT_230922C/"$model"/barcode15.fastq
        safe_cat ONT_230926_HW/"$model"/barcode15.fastq
    } >"$final_dir/$model"/ATCC_19119__202309.fq.gz

    {
        safe_cat ONT_230922A/"$model"/barcode75.fastq
        safe_cat ONT_230922B/"$model"/barcode67.fastq
        safe_cat ONT_230922C/"$model"/barcode03.fastq
        safe_cat ONT_230926_HW/"$model"/barcode03.fastq
        safe_cat ONT_230922A/"$model"/barcode83.fastq
        safe_cat ONT_230922B/"$model"/barcode75.fastq
        safe_cat ONT_230922C/"$model"/barcode11.fastq
        safe_cat ONT_230926_HW/"$model"/barcode11.fastq
    } >"$final_dir/$model"/ATCC_25922__202309.fq.gz

    {
        safe_cat ONT_230922A/"$model"/barcode77.fastq
        safe_cat ONT_230922B/"$model"/barcode69.fastq
        safe_cat ONT_230922C/"$model"/barcode05.fastq
        safe_cat ONT_230926_HW/"$model"/barcode05.fastq
        safe_cat ONT_230922A/"$model"/barcode85.fastq
        safe_cat ONT_230922B/"$model"/barcode77.fastq
        safe_cat ONT_230922C/"$model"/barcode13.fastq
        safe_cat ONT_230926_HW/"$model"/barcode13.fastq
    } >"$final_dir/$model"/ATCC_33560__202309.fq.gz

    {
        safe_cat ONT_230922A/"$model"/barcode76.fastq
        safe_cat ONT_230922B/"$model"/barcode68.fastq
        safe_cat ONT_230922C/"$model"/barcode04.fastq
        safe_cat ONT_230926_HW/"$model"/barcode04.fastq
        safe_cat ONT_230922A/"$model"/barcode84.fastq
        safe_cat ONT_230922B/"$model"/barcode76.fastq
        safe_cat ONT_230922C/"$model"/barcode12.fastq
        safe_cat ONT_230926_HW/"$model"/barcode12.fastq
    } >"$final_dir/$model"/ATCC_35221__202309.fq.gz

    {
        safe_cat ONT_230922A/"$model"/barcode78.fastq
        safe_cat ONT_230922B/"$model"/barcode70.fastq
        safe_cat ONT_230922C/"$model"/barcode06.fastq
        safe_cat ONT_230926_HW/"$model"/barcode06.fastq
        safe_cat ONT_230922A/"$model"/barcode86.fastq
        safe_cat ONT_230922B/"$model"/barcode78.fastq
        safe_cat ONT_230922C/"$model"/barcode14.fastq
        safe_cat ONT_230926_HW/"$model"/barcode14.fastq
    } >"$final_dir/$model"/ATCC_35897__202309.fq.gz

    {
        safe_cat ONT_230922A/"$model"/barcode80.fastq
        safe_cat ONT_230922B/"$model"/barcode72.fastq
        safe_cat ONT_230922C/"$model"/barcode08.fastq
        safe_cat ONT_230926_HW/"$model"/barcode08.fastq
        safe_cat ONT_230922A/"$model"/barcode88.fastq
        safe_cat ONT_230922B/"$model"/barcode80.fastq
        safe_cat ONT_230922C/"$model"/barcode16.fastq
        safe_cat ONT_230926_HW/"$model"/barcode16.fastq
        safe_cat ONT_230926_HW/"$model"/barcode19.fastq
    } >"$final_dir/$model"/ATCC_BAA-679__202309.fq.gz

done
