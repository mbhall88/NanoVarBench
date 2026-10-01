#!/usr/bin/env bash
# Seam 1: run the whole workflow on the tiny fixture (tests/fixture/), with Dorado on CPU,
# then check the outputs. Needs snakemake (with apptainer and conda) on PATH.
#
#   DORADO=/path/to/dorado-2.1.2/bin/dorado [DORADO_MODELS_DIR=dir] \
#       [CLAIR3_MODELS_DIR=dir] [OUTDIR=dir] tests/seam1.sh [extra snakemake args]
#
# It covers Arms A-D on a hac and a sup Read set at every Depth in DEPTHS (default
# "5 10 25 50", the config default), and checks that the chromosome's actual depth is within
# 5% of each Depth. READ_MODELS (default "hac sup") and DEPTHS pick a subset for a quicker run.
# The Dorado CPU re-run is only done at the largest Depth, and the dorado aligner check
# (#13) on hac at 25x when both are in the run.
# DORADO_MODELS_DIR defaults to $OUTDIR/models, which the workflow fills with
# `dorado download` (needs internet). CLAIR3_MODELS_DIR defaults to
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
READ_MODELS=${READ_MODELS:-hac sup}
DEPTHS=${DEPTHS:-5 10 25 50}
CPU_DEPTH=$(tr ' ' '\n' <<< "$DEPTHS" | sort -n | tail -1)
# The dorado aligner check (#13) runs on hac at 25x, so only when both are in this run.
ALIGNER_CHECK=
if [[ " $READ_MODELS " == *" hac "* && " $DEPTHS " == *" 25 "* ]]; then
  ALIGNER_CHECK="{samples: [ATCC_25922_fixture], read_models: [hac], depths: [25]}"
fi
yaml_list() { local out=; for x in "$@"; do out+="${out:+, }$x"; done; echo "[$out]"; }
echo "Seam 1 output: $OUTDIR"

md5() { awk -v f="$1" '$2 == f {print $1}' "$FIXTURE/md5.txt"; }
mkdir -p "$OUTDIR/seam1_config"
printf 'sample\tspecies\tdnd\ttruth_url\ttruth_md5\n%s\t%s\t%s\t%s\t%s\n' \
  ATCC_25922_fixture "Escherichia coli" false \
  "file://$FIXTURE/ATCC_25922_fixture.tar.gz" "$(md5 ATCC_25922_fixture.tar.gz)" \
  > "$OUTDIR/seam1_config/samples.tsv"
printf 'sample\tread_model\trun\tfastq_url\tfastq_md5\tfastq_bytes\n' > "$OUTDIR/seam1_config/runs.tsv"
for rm in $READ_MODELS; do
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' ATCC_25922_fixture "$rm" "FIXTURE${rm^^}" \
    "file://$FIXTURE/reads.$rm.fastq.gz" "$(md5 reads.$rm.fastq.gz)" \
    "$(stat -c %s "$FIXTURE/reads.$rm.fastq.gz")" >> "$OUTDIR/seam1_config/runs.tsv"
done
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
  read_models: $(yaml_list $READ_MODELS)
  depths: $(yaml_list $DEPTHS)
  arms: [A, B, C, D]
dorado_aligner_check: ${ALIGNER_CHECK:-{\}}
dorado:
  bin: {"2.1.2": $DORADO}
  models_dir: $MODELS
  device: cpu
  cpu_run: {depths: [$CPU_DEPTH], threads: 2}
clair3:
  models_dir: $CLAIR3_MODELS
YAML

cd "$UPDATE"
snakemake -s workflow/Snakefile --configfile "$OUTDIR/seam1_config/config.yaml" \
  --cores "$CORES" --software-deployment-method apptainer conda \
  --apptainer-args "--bind $UPDATE,$OUTDIR,$CLAIR3_MODELS" \
  --conda-prefix "$UPDATE/.snakemake/conda" --apptainer-prefix "$UPDATE/.snakemake/singularity" \
  --show-failed-logs "$@"

python3 "$HERE/check_seam1.py" "$OUTDIR" "$FIXTURE/expected_tp.tsv" \
  --read-models $READ_MODELS --depths $DEPTHS --cpu-depth "$CPU_DEPTH"
if [ -n "$ALIGNER_CHECK" ]; then python3 "$HERE/check_aligner_check.py" "$OUTDIR"; fi
