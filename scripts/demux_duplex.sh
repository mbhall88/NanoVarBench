#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
DEMUX_SCRIPT="$SCRIPT_DIR"/dupdemux.py

function print_usage_and_exit() {
    echo "Usage: $0 [-r run_dir] [-o outdir]" >&2
    echo "  -r run_dir: directory containing duplex BAMs to demux" >&2
    echo "  -o outdir: output directory" >&2
    echo "  -h: print this help message" >&2
    exit 0
}

function parse_params() {
    while getopts 'm:r:o:h' flag; do
        case "${flag}" in
        h) print_usage_and_exit ;;
        r) run_dir="$(realpath "${OPTARG}")" ;;
        o) final_dir="$(realpath "${OPTARG}")" ;;
        *) error "Unexpected option ${flag}" ;;
        esac
    done
}

parse_params "$@"

for bamfile in "$run_dir"/*.bam; do
    sample=$(basename "$bamfile" .bam)
    # split sample on . with the first part being the run and the second part the model
    IFS='.' read -r -a array <<<"$sample"
    run=${array[0]}
    model=${array[1]}
    model="dna_r10.4.1_e8.2_400bps_${model}@v4.3.0"
    outdir="$final_dir"/"$run"/"$model"
    mkdir -p "$outdir"
    summary="/data/gpfs/projects/punim2009/data/ont/basecall/${run}/${model}/sequencing_summary.txt"
    
    echo "Demultiplexing $bamfile"
    python "$DEMUX_SCRIPT" -dfU -o "$outdir" "$summary" "$bamfile" 2> "$outdir"/"demux_$sample".log &
done