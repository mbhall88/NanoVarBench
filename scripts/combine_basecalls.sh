#!/usr/bin/env bash
set -euo pipefail


function print_usage_and_exit() {
    echo "Usage: $0 [-m model] [-d model_dir] [-r run_dir] [-o outdir]" >&2
    echo "  -m model: dorado model name" >&2
    echo "  -r run_dir: run directory (with pod5)" >&2
    echo "  -o outdir: output directory" >&2
    echo "  -h: print this help message" >&2
    exit 0
}

function parse_params() {
    while getopts 'm:r:o:h' flag; do
        case "${flag}" in
        h) print_usage_and_exit ;;
        m) model="${OPTARG}" ;;
        r) run_dir="${OPTARG}" ;;
        o) final_dir="${OPTARG}" ;;
        *) error "Unexpected option ${flag}" ;;
        esac
    done
}

parse_params "$@"

run_dir=$(realpath "$run_dir")
final_dir=$(realpath "$final_dir")

mkdir -p "$final_dir/$model"
cd "$run_dir" || exit 1

{
    cat ONT_230922A/"$model"/demultiplexed/barcode73.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode65.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode01.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode01.fastq.gz
    cat ONT_230922A/"$model"/demultiplexed/barcode81.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode73.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode09.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode09.fastq.gz
} >"$final_dir/$model"/ATCC_10708__202309.fq.gz

{
    cat ONT_230922A/"$model"/demultiplexed/barcode74.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode66.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode02.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode02.fastq.gz
    cat ONT_230922A/"$model"/demultiplexed/barcode82.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode74.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode10.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode10.fastq.gz
} >"$final_dir/$model"/ATCC_14035__202309.fq.gz

{
    cat ONT_230922A/"$model"/demultiplexed/barcode89.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode81.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode17.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode17.fastq.gz
    cat ONT_230922A/"$model"/demultiplexed/barcode90.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode82.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode18.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode18.fastq.gz
} >"$final_dir/$model"/ATCC_17802__202309.fq.gz

{
    cat ONT_230922A/"$model"/demultiplexed/barcode79.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode71.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode07.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode07.fastq.gz
    cat ONT_230922A/"$model"/demultiplexed/barcode87.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode79.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode15.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode15.fastq.gz
} >"$final_dir/$model"/ATCC_19119__202309.fq.gz

{
    cat ONT_230922A/"$model"/demultiplexed/barcode75.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode67.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode03.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode03.fastq.gz
    cat ONT_230922A/"$model"/demultiplexed/barcode83.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode75.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode11.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode11.fastq.gz
} >"$final_dir/$model"/ATCC_25922__202309.fq.gz

{
    cat ONT_230922A/"$model"/demultiplexed/barcode77.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode69.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode05.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode05.fastq.gz
    cat ONT_230922A/"$model"/demultiplexed/barcode85.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode77.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode13.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode13.fastq.gz
} >"$final_dir/$model"/ATCC_33560__202309.fq.gz

{
    cat ONT_230922A/"$model"/demultiplexed/barcode76.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode68.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode04.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode04.fastq.gz
    cat ONT_230922A/"$model"/demultiplexed/barcode84.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode76.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode12.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode12.fastq.gz
} >"$final_dir/$model"/ATCC_35221__202309.fq.gz

{
    cat ONT_230922A/"$model"/demultiplexed/barcode78.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode70.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode06.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode06.fastq.gz
    cat ONT_230922A/"$model"/demultiplexed/barcode86.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode78.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode14.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode14.fastq.gz
} >"$final_dir/$model"/ATCC_35897__202309.fq.gz

{
    cat ONT_230922A/"$model"/demultiplexed/barcode80.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode72.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode08.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode08.fastq.gz
    cat ONT_230922A/"$model"/demultiplexed/barcode88.fastq.gz
    cat ONT_230922B/"$model"/demultiplexed/barcode80.fastq.gz
    cat ONT_230922C/"$model"/demultiplexed/barcode16.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode16.fastq.gz
    cat ONT_230926_HW/"$model"/demultiplexed/barcode19.fastq.gz
} >"$final_dir/$model"/ATCC_BAA-679__202309.fq.gz
