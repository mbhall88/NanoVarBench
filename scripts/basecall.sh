#!/usr/bin/env bash
set -euxo pipefail

dorado=dorado

function print_usage_and_exit() {
    echo "Usage: $0 [-m model] [-d model_dir] [-r run_dir] [-k barcode_kit] [-D dorado_exec] [-o outdir]" >&2
    echo "  -m model: dorado model name" >&2
    echo "  -d model_dir: dorado model directory" >&2
    echo "  -r run_dir: run directory (with pod5)" >&2
    echo "  -o outdir: output directory" >&2
    echo "  -k barcode_kit: barcode kit name" >&2
    echo "  -D dorado_exec: path to dorado executable" >&2
    echo "  -h: print this help message" >&2
    exit 0
}

function parse_params() {
    while getopts 'm:d:r:o:D:k:o:h' flag; do
        case "${flag}" in
        h) print_usage_and_exit ;;
        m) model="${OPTARG}" ;;
        d) model_dir=$(realpath "${OPTARG}") ;;
        k) kit="${OPTARG}" ;;
        r) run_dir=$(realpath "${OPTARG}") ;;
        D) dorado="${OPTARG}" ;;
        o) outdir=$(realpath "${OPTARG}") ;;
        *) error "Unexpected option ${flag}" ;;
        esac
    done
}

parse_params "$@"

mkdir -p "$outdir"
cd "$outdir" || exit 1

printf "\n\n\n\n%s %s\n--------------------------------------------------------------------------------\n" "$run_dir" "$model" >&2

mkdir -p "$model"
"$dorado" basecaller --recursive "$model_dir"/"$model" "$run_dir" --kit-name "$kit" >"$model"/reads.bam

cd "$model" || exit 1
"$dorado" summary reads.bam >sequencing_summary.txt
"$dorado" demux --output-dir demultiplexed --no-classify reads.bam

cd demultiplexed || exit 1
for b in *_barcode*.bam; do
    name=$(echo "$b" | grep -oP "barcode\d+")
    samtools fastq -T '*' "$b" | tr '\t' ' ' | pigz -p8 >"$name".fastq.gz
done
samtools fastq -T '*' unclassified.bam | tr '\t' ' ' | pigz -p8 >unclassified.fastq.gz
rm ./*.bam
