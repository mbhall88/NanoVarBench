#!/usr/bin/env bash
# Seam 2: run the workflow's aggregation on hand-written vcfdist, depth and benchmark inputs
# (tests/seam2/), then check the output tables. No GPU, Slurm, containers or
# downloads; it takes a few seconds. Needs snakemake on PATH.
#
#   [OUTDIR=dir] tests/seam2.sh [extra snakemake args]
#
# It runs every Sample in config/samples.tsv (all 14) x hac x 25 and 50x x Arms A-D. Each
# Read set and Arm gets the files a real run leaves where the aggregate rule reads them:
#   ATCC_25922__202309 hac 50x Arm D  tests/seam2/worked_D (the worked example)
#   ATCC_25922__202309 hac 50x Arm A  tests/seam2/worked_A
#   everything else                   tests/seam2/filler
# and ATCC_25922__202309 hac 50x gets depth/worked.coverage.tsv (the rest get
# depth/filler.coverage.tsv), and the same split for benchmarks/{worked,filler}/ (one file
# per step and Arm, plus the Dorado CPU re-run at 50x). Only the aggregate and benchmarks
# rules may run (--allowed-rules), so Snakemake never tries to rebuild the hand-written
# inputs. Which caller, version, Calling model, container and hardware each row carries comes
# from the config, so there are no per-job caller files; the test config sets
# benchmark_hardware to made-up models to show it is the config's. OUTDIR defaults to a new
# temporary directory.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
UPDATE=$(dirname "$HERE")
IN=$HERE/seam2
OUTDIR=${OUTDIR:-$(mktemp -d)}
OUTDIR=$(mkdir -p "$OUTDIR" && cd "$OUTDIR" && pwd)
echo "Seam 2 output: $OUTDIR"

SAMPLES=$(awk -F'\t' 'NR > 1 {print $1}' "$UPDATE/config/samples.tsv")
DEPTHS="25 50"
ARMS="A B C D"
READ_MODEL=hac
WORKED=ATCC_25922__202309
CPU_DEPTH=50  # dorado.cpu_run.depths in config/config.yaml
ALN_A=minimap2-2.26-map-ont
ALN_SHARED=minimap2-2.31-lrhq

# make_inputs DIR: hand-written inputs for every Read set and Arm, where a run leaves them.
make_inputs() {
  local out=$1 s d a src bsrc
  for s in $SAMPLES; do
    for d in $DEPTHS; do
      cov=filler
      [[ $s == "$WORKED" && $d == 50 ]] && cov=worked
      mkdir -p "$out/results/depth" "$out/results/benchmarks/align" \
        "$out/results/benchmarks/call_clair3" "$out/results/benchmarks/call_dorado" \
        "$out/results/benchmarks/call_dorado_cpu"
      cp "$IN/depth/$cov.coverage.tsv" "$out/results/depth/$s.$READ_MODEL.${d}x.coverage.tsv"
      bsrc=filler
      [[ $s == "$WORKED" && $d == 50 ]] && bsrc=worked
      bk=$s.$READ_MODEL.${d}x
      cp "$IN/benchmarks/$bsrc/align_A.tsv" "$out/results/benchmarks/align/$bk.$ALN_A.tsv"
      cp "$IN/benchmarks/$bsrc/align_B.tsv" "$out/results/benchmarks/align/$bk.$ALN_SHARED.tsv"
      for a in $ARMS; do
        src=filler
        [[ $s == "$WORKED" && $d == 50 && $a == D ]] && src=worked_D
        [[ $s == "$WORKED" && $d == 50 && $a == A ]] && src=worked_A
        score=$out/work/score/$s/$READ_MODEL/${d}x/$a
        mkdir -p "$score/sweep" "$score/pass"
        cp "$IN/$src/sweep/precision-recall-summary.tsv" "$IN/$src/sweep/precision-recall.tsv" "$score/sweep/"
        cp "$IN/$src/pass/precision-recall-summary.tsv" "$score/pass/"
        caller=call_clair3; [[ $a == D ]] && caller=call_dorado
        cp "$IN/benchmarks/$bsrc/call_$a.tsv" "$out/results/benchmarks/$caller/$bk.$a.tsv"
        if [[ $a == D && $d == "$CPU_DEPTH" ]]; then  # the timing-only CPU re-run
          cp "$IN/benchmarks/$bsrc/call_D_cpu.tsv" "$out/results/benchmarks/call_dorado_cpu/$bk.D.tsv"
        fi
      done
    done
  done

  mkdir -p "$out/seam2_config"
  cat > "$out/seam2_config/config.yaml" <<YAML
samples: $UPDATE/config/samples.tsv
runs: $UPDATE/config/runs.tsv
work_dir: $out/work
results_dir: $out/results
benchmark_hardware: {gpu: "Test GPU", cpu: "Test CPU"}
run:
  samples: [$(echo $SAMPLES | sed 's/ /, /g')]
  read_models: [$READ_MODEL]
  depths: [$(echo $DEPTHS | sed 's/ /, /g')]
  arms: [$(echo $ARMS | sed 's/ /, /g')]
YAML
}

make_inputs "$OUTDIR"

cd "$UPDATE"
# (targets first: --allowed-rules and --configfile swallow the arguments after them)
snakemake aggregate benchmarks -s workflow/Snakefile --cores 1 --show-failed-logs "$@" \
  --configfile "$OUTDIR/seam2_config/config.yaml" --allowed-rules aggregate benchmarks

python3 "$HERE/check_seam2.py" "$OUTDIR"
