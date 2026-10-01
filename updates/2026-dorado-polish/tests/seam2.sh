#!/usr/bin/env bash
# Seam 2: run the workflow's aggregation on hand-written vcfdist, depth and caller inputs
# (tests/seam2/), then check the output tables. No GPU, Slurm, containers or downloads;
# it takes a few seconds. Needs snakemake on PATH.
#
#   [OUTDIR=dir] tests/seam2.sh [extra snakemake args]
#
# It runs every Sample in config/samples.tsv (all 14) x hac x 25 and 50x x Arms A-D. Each
# Read set and Arm gets the files a real run leaves where the aggregate rule reads them:
#   ATCC_25922__202309 hac 50x Arm D  tests/seam2/worked_D (the worked example)
#   ATCC_25922__202309 hac 50x Arm A  tests/seam2/worked_A
#   everything else                   tests/seam2/filler
# and ATCC_25922__202309 hac 50x gets depth/worked.coverage.tsv (the rest get
# depth/filler.coverage.tsv). Only the aggregate rule may run (--allowed-rules), so
# Snakemake never tries to rebuild the hand-written inputs. OUTDIR defaults to a new
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

for s in $SAMPLES; do
  for d in $DEPTHS; do
    cov=filler
    [[ $s == "$WORKED" && $d == 50 ]] && cov=worked
    mkdir -p "$OUTDIR/results/depth"
    cp "$IN/depth/$cov.coverage.tsv" "$OUTDIR/results/depth/$s.$READ_MODEL.${d}x.coverage.tsv"
    for a in $ARMS; do
      src=filler
      [[ $s == "$WORKED" && $d == 50 && $a == D ]] && src=worked_D
      [[ $s == "$WORKED" && $d == 50 && $a == A ]] && src=worked_A
      score=$OUTDIR/work/score/$s/$READ_MODEL/${d}x/$a
      mkdir -p "$score/sweep" "$score/pass" "$OUTDIR/results/calls/$s/$READ_MODEL/${d}x"
      cp "$IN/$src/sweep/precision-recall-summary.tsv" "$IN/$src/sweep/precision-recall.tsv" "$score/sweep/"
      cp "$IN/$src/pass/precision-recall-summary.tsv" "$score/pass/"
      cp "$IN/caller_info/$a.caller_info.tsv" "$OUTDIR/results/calls/$s/$READ_MODEL/${d}x/"
    done
  done
done

mkdir -p "$OUTDIR/seam2_config"
cat > "$OUTDIR/seam2_config/config.yaml" <<YAML
samples: $UPDATE/config/samples.tsv
runs: $UPDATE/config/runs.tsv
work_dir: $OUTDIR/work
results_dir: $OUTDIR/results
run:
  samples: [$(echo $SAMPLES | sed 's/ /, /g')]
  read_models: [$READ_MODEL]
  depths: [$(echo $DEPTHS | sed 's/ /, /g')]
  arms: [$(echo $ARMS | sed 's/ /, /g')]
YAML

cd "$UPDATE"
snakemake -s workflow/Snakefile --configfile "$OUTDIR/seam2_config/config.yaml" \
  --cores 1 --allowed-rules aggregate --show-failed-logs "$@" aggregate

python3 "$HERE/check_seam2.py" "$OUTDIR"
