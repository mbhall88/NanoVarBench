#!/usr/bin/env bash
set -euo pipefail

dorado=dorado
num_gpus=1
memory=8G
time=12h

function print_usage_and_exit() {
    echo "Usage: $0 [-m model] [-d model_dir] [-r run_dir] [-o outdir] [-D dorado_exec] [-o outdir]" >&2
    echo "  -m model: dorado model name" >&2
    echo "  -d model_dir: dorado model directory" >&2
    echo "  -r run_dir: run directory (with pod5)" >&2
    echo "  -D dorado_exec: path to dorado executable [default: $dorado]" >&2
    echo "  -o outdir: output directory" >&2
    echo "  -s basecall_script: path to basecall script" >&2
    echo "  -g num_gpus: number of GPUs to use [default: $num_gpus]" >&2
    echo "  -M memory: memory to use [default: $memory]" >&2
    echo "  -t time: time to use [default: $time]" >&2
    echo "  -h: print this help message" >&2
    exit 0
}

function parse_params() {
    while getopts 'm:d:r:D:s:g:M:t:o:h' flag; do
        case "${flag}" in
        h) print_usage_and_exit ;;
        m) model="${OPTARG}" ;;
        d) model_dir="$(realpath "${OPTARG}")" ;;
        r) run_dir="$(realpath "${OPTARG}")" ;;
        D) dorado="${OPTARG}" ;;
        s) basecall_script="$(realpath "${OPTARG}")" ;;
        g) num_gpus="${OPTARG}" ;;
        M) memory="${OPTARG}" ;;
        t) time="${OPTARG}" ;;
        o) outdir="$(realpath "${OPTARG}")" ;;
        *) error "Unexpected option ${flag}" ;;
        esac
    done
}

parse_params "$@"

kit=SQK-NBD114-96
for run in ONT_230922A ONT_230922B; do
    ssubmit -t "$time" -m "$memory" "basecall_${run}_${model}" "bash $basecall_script -m $model -d $model_dir -r $run_dir/$run -o $outdir/$run -D $dorado -k $kit" -- --gres=gpu:"$num_gpus" -p gpu-a100
done

kit=SQK-RBK114-96
for run in ONT_230922C ONT_230926_HW; do
    ssubmit -t "$time" -m "$memory" "basecall_${run}_${model}" "bash $basecall_script -m $model -d $model_dir -r $run_dir/$run -D $dorado -o $outdir/$run -k $kit" -- --gres=gpu:"$num_gpus" -p gpu-a100
done

