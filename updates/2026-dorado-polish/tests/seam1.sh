#!/usr/bin/env bash
# Seam 1: run the whole workflow on the tiny fixture (tests/fixture/), with Dorado on CPU,
# then check the outputs. Needs snakemake (with apptainer and conda) on PATH.
#
#   DORADO=/path/to/dorado-2.1.2/bin/dorado [DORADO_MODELS_DIR=dir] \
#       [CLAIR3_MODELS_DIR=dir] [OUTDIR=dir] tests/seam1.sh [extra snakemake args]
#
# It covers Arms A-D. DORADO_MODELS_DIR defaults to $OUTDIR/models, which the workflow fills
# with `dorado download` (needs internet). CLAIR3_MODELS_DIR defaults to
# $OUTDIR/clair3_models, which the workflow fills by downloading the HKU PyTorch Calling
# models for Arm C and checking their SHA256s (needs internet). OUTDIR defaults to a new
# temporary directory.
# Conda envs and container images are cached in .snakemake/ and shared with real runs.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
UPDATE=$(dirname "$HERE")
FIXTURE=$HERE/fixture
: "${DORADO:?set DORADO to the dorado 2.1.2 binary}"
OUTDIR=${OUTDIR:-$(mktemp -d)}
OUTDIR=$(mkdir -p "$OUTDIR" && cd "$OUTDIR" && pwd)
MODELS=${DORADO_MODELS_DIR:-$OUTDIR/models}
CLAIR3_MODELS=${CLAIR3_MODELS_DIR:-$OUTDIR/clair3_models}
CLAIR3_MODELS=$(mkdir -p "$CLAIR3_MODELS" && cd "$CLAIR3_MODELS" && pwd)
CORES=${CORES:-8}
echo "Seam 1 output: $OUTDIR"

md5() { awk -v f="$1" '$2 == f {print $1}' "$FIXTURE/md5.txt"; }
mkdir -p "$OUTDIR/seam1_config"
printf 'sample\tspecies\tdnd\ttruth_url\ttruth_md5\n%s\t%s\t%s\t%s\t%s\n' \
  ATCC_25922_fixture "Escherichia coli" false \
  "file://$FIXTURE/ATCC_25922_fixture.tar.gz" "$(md5 ATCC_25922_fixture.tar.gz)" \
  > "$OUTDIR/seam1_config/samples.tsv"
printf 'sample\tread_model\trun\tfastq_url\tfastq_md5\tfastq_bytes\n%s\t%s\t%s\t%s\t%s\t%s\n' \
  ATCC_25922_fixture hac FIXTURE1 "file://$FIXTURE/reads.fastq.gz" "$(md5 reads.fastq.gz)" \
  "$(stat -c %s "$FIXTURE/reads.fastq.gz")" > "$OUTDIR/seam1_config/runs.tsv"
cat > "$OUTDIR/seam1_config/config.yaml" <<YAML
samples: $OUTDIR/seam1_config/samples.tsv
runs: $OUTDIR/seam1_config/runs.tsv
reads_dir: $OUTDIR/downloads/reads
truth_dir: $OUTDIR/downloads/truth
download: {method: curl}
work_dir: $OUTDIR/work
results_dir: $OUTDIR/results
run:
  samples: [ATCC_25922_fixture]
  read_models: [hac]
  depths: [25]
  arms: [A, B, C, D]
dorado:
  bin: {"2.1.2": $DORADO}
  models_dir: $MODELS
  device: cpu
  cpu_run: {depths: [25], threads: 2}
clair3:
  models_dir: $CLAIR3_MODELS
YAML

cd "$UPDATE"
snakemake -s workflow/Snakefile --configfile "$OUTDIR/seam1_config/config.yaml" \
  --cores "$CORES" --software-deployment-method apptainer conda \
  --apptainer-args "--bind $UPDATE,$OUTDIR,$CLAIR3_MODELS" \
  --conda-prefix "$UPDATE/.snakemake/conda" --apptainer-prefix "$UPDATE/.snakemake/singularity" \
  --show-failed-logs "$@"

python3 "$HERE/check_seam1.py" "$OUTDIR" "$FIXTURE/expected_tp.tsv"
