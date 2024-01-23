#!/usr/bin/env bash
set -euxo pipefail

dorado=dorado
kit=""
subcommand="basecall"

function print_usage_and_exit() {
    echo "Usage: $0 [-m model] [-d model_dir] [-r run_dir] [-k barcode_kit] [-D dorado_exec] [-o outdir] [-x]" >&2
    echo "  -m model: dorado model name" >&2
    echo "  -d model_dir: dorado model directory" >&2
    echo "  -r run_dir: run directory (with pod5)" >&2
    echo "  -o outdir: output directory" >&2
    echo "  -k barcode_kit: barcode kit name" >&2
    echo "  -D dorado_exec: path to dorado executable" >&2
    echo "  -x: duplex basecalling" >&2
    echo "  -h: print this help message" >&2
    exit 0
}

function parse_params() {
    while getopts 'm:d:r:o:D:k:o:xh' flag; do
        case "${flag}" in
        h) print_usage_and_exit ;;
        m) model="${OPTARG}" ;;
        d) model_dir=$(realpath "${OPTARG}") ;;
        k) kit="${OPTARG}" ;;
        r) run_dir=$(realpath "${OPTARG}") ;;
        D) dorado="${OPTARG}" ;;
        o) outdir=$(realpath "${OPTARG}") ;;
        x) subcommand="duplex" ;;
        *) error "Unexpected option ${flag}" ;;
        esac
    done

    # don't allow duplex basecalling with a kit
    if [[ "$subcommand" == "duplex" && -n "$kit" ]]; then
        echo "Duplex basecalling does not support a barcode kit" >&2
        exit 1
    fi

}

parse_params "$@"

mkdir -p "$outdir"
cd "$outdir" || exit 1

printf "\n\n\n\n%s %s\n--------------------------------------------------------------------------------\n" "$run_dir" "$model" >&2

mkdir -p "$model"
opts=(--recursive "$model_dir"/"$model" "$run_dir")
if [[ -n "$kit" ]]; then
    opts+=(--kit-name "$kit")
fi
"$dorado" "$subcommand" "${opts[@]}" >"$model"/reads.bam

cd "$model" || exit 1
"$dorado" summary reads.bam >sequencing_summary.txt
if [[ -n "$kit" ]]; then
    "$dorado" demux --output-dir demultiplexed --no-classify reads.bam

    cd demultiplexed || exit 1
    for b in *_barcode*.bam; do
        name=$(echo "$b" | grep -oP "barcode\d+")
        samtools fastq -T '*' "$b" | tr '\t' ' ' | pigz -p8 >"$name".fastq.gz
    done
    samtools fastq -T '*' unclassified.bam | tr '\t' ' ' | pigz -p8 >unclassified.fastq.gz
else
    if [[ "$subcommand" == "duplex" ]]; then
        samtools view -h -d 'dx:1' reads.bam | samtools fastq -T '*' - | tr '\t' ' ' | pigz -p8 >reads.fastq.gz
    else
        samtools fastq -T '*' reads.bam | tr '\t' ' ' | pigz -p8 >reads.fastq.gz
    fi
fi

rm ./*.bam
