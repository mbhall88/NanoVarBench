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
# The AF filter analysis (#20) runs for Arm C on every Read set at AF_THRESHOLDS (default
# "0.5 0.65 0.8"), and a dry-run with a fresh work_dir reading this run's (input_work_dir) checks
# that only its own jobs are scheduled. The figures and Table S1 then have its series at
# AF_FIGURE_THRESHOLD (default 0.65, figures.af_filter), and the same figures and table rendered
# from this run's aggregated tables with the analysis off are checked to be without it (#28).
# The smallvar pilot (#29) runs for Arm D on every Read set, on CPU like the rest of Dorado here,
# with the AF filter at the lowest AF threshold beside it, and a reuse dry-run like the AF
# filter's checks that only its own jobs are scheduled.
# DORADO_MODELS_DIR defaults to $OUTDIR/models, which the workflow fills with
# `dorado download` (needs internet); the seam then checks both models arrived whole (#32). CLAIR3_MODELS_DIR defaults to
# $OUTDIR/clair3_models, which the workflow fills by downloading the HKU PyTorch Calling
# models for Arm C and checking their SHA256s (needs internet). OUTDIR defaults to a new
# temporary directory.
# It also renders the figures and tables (#14) from the fixture's aggregated tables, and
# check_seam1.py confirms they exist.
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
AF_THRESHOLDS=${AF_THRESHOLDS:-0.5 0.65 0.8}
AF_FIGURE_THRESHOLD=${AF_FIGURE_THRESHOLD:-0.65}
SMALLVAR_AF=$(tr ' ' '\n' <<< "$AF_THRESHOLDS" | sort -g | head -1)
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
clair3_af_filter: {arms: [C], thresholds: $(yaml_list $AF_THRESHOLDS), reference_arm: D}
figures: {af_filter: {arm: C, threshold: $AF_FIGURE_THRESHOLD}}
smallvar_pilot: {arms: [D], compare_arms: [C, D], af_filter: {arm: C, threshold: $SMALLVAR_AF}}
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

# With no DORADO_MODELS_DIR, the models dir started empty (#32): the polishing and smallvar
# models must each have been downloaded whole, with no temporary download left behind.
if [ -z "${DORADO_MODELS_DIR:-}" ]; then
  n=0
  for d in "$MODELS"/*/; do
    for f in config.toml weights.pt; do
      [ -s "$d$f" ] || { echo "FAIL fresh models dir: $d$f missing or empty" >&2; exit 1; }
    done
    case $(basename "$d") in tmp.*) echo "FAIL fresh models dir: leftover $d" >&2; exit 1;; esac
    n=$((n + 1))
  done
  [ "$n" -eq 2 ] || { echo "FAIL fresh models dir: $n models in $MODELS, want 2" >&2; exit 1; }
  echo "fresh models dir checks passed"
fi

# The AF filter reusing this run's Read sets, alignments and scores from a fresh work_dir: a
# dry-run, whose job counts check_af_filter.py reads.
cat > "$OUTDIR/seam1_config/af_reuse.yaml" <<YAML
work_dir: $OUTDIR/af_reuse/work
results_dir: $OUTDIR/af_reuse/results
clair3_af_filter: {arms: [C], thresholds: $(yaml_list $AF_THRESHOLDS), reference_arm: D, input_work_dir: $OUTDIR/work}
YAML
snakemake clair3_af_filter -n -s workflow/Snakefile \
  --configfile "$OUTDIR/seam1_config/config.yaml" "$OUTDIR/seam1_config/af_reuse.yaml" \
  --software-deployment-method apptainer conda \
  --conda-prefix "$UPDATE/.snakemake/conda" --apptainer-prefix "$UPDATE/.snakemake/singularity" \
  > "$OUTDIR/af_reuse.dry_run.txt"
python3 "$HERE/check_af_filter.py" "$OUTDIR" --arms C --reference-arm D --thresholds $AF_THRESHOLDS \
  --read-models $READ_MODELS --depths $DEPTHS --reuse-dry-run "$OUTDIR/af_reuse.dry_run.txt"

# The figures and Table S1 with the AF filter off: this run's aggregated tables in a fresh
# results_dir, the analysis disabled, and only the figure and table rules allowed to run.
mkdir -p "$OUTDIR/af_off/results/tables"
cp "$OUTDIR"/results/tables/{results,pr_curves,depth}.tsv "$OUTDIR/af_off/results/tables/"
cat > "$OUTDIR/seam1_config/af_off.yaml" <<YAML
work_dir: $OUTDIR/af_off/work
results_dir: $OUTDIR/af_off/results
clair3_af_filter: {arms: []}
YAML
OFF=$OUTDIR/af_off/results
snakemake "$OFF/figures/fig1_best_f1_depth.png" "$OFF/figures/fig1_best_f1_depth.svg" \
  "$OFF/figures/fig2_pr_curves.png" "$OFF/figures/fig2_pr_curves.svg" \
  "$OFF/figures/fig3_per_sample_best_f1.png" "$OFF/figures/fig3_per_sample_best_f1.svg" \
  "$OFF/tables/table_s1_per_sample.csv" -s workflow/Snakefile \
  --configfile "$OUTDIR/seam1_config/config.yaml" "$OUTDIR/seam1_config/af_off.yaml" \
  --cores "$CORES" --software-deployment-method conda \
  --conda-prefix "$UPDATE/.snakemake/conda" \
  --allowed-rules fig1_best_f1_depth fig2_pr_curves fig3_per_sample_best_f1 table_s1_per_sample
python3 "$HERE/check_af_figures.py" "$OUTDIR" --arm C --threshold "$AF_FIGURE_THRESHOLD" \
  --read-models $READ_MODELS --depths $DEPTHS

# The smallvar pilot reusing this run's Read sets, alignments and scores (and its AF filter
# scores) from a fresh work_dir: a dry-run, whose job counts check_smallvar_pilot.py reads.
cat > "$OUTDIR/seam1_config/smallvar_reuse.yaml" <<YAML
work_dir: $OUTDIR/smallvar_reuse/work
results_dir: $OUTDIR/smallvar_reuse/results
smallvar_pilot:
  arms: [D]
  compare_arms: [C, D]
  af_filter: {arm: C, threshold: $SMALLVAR_AF, work_dir: $OUTDIR/work}
  input_work_dir: $OUTDIR/work
YAML
snakemake smallvar_pilot -n -s workflow/Snakefile \
  --configfile "$OUTDIR/seam1_config/config.yaml" "$OUTDIR/seam1_config/smallvar_reuse.yaml" \
  --software-deployment-method apptainer conda \
  --conda-prefix "$UPDATE/.snakemake/conda" --apptainer-prefix "$UPDATE/.snakemake/singularity" \
  > "$OUTDIR/smallvar_reuse.dry_run.txt"
python3 "$HERE/check_smallvar_pilot.py" "$OUTDIR" --arms D --compare-arms C D --af-arm C \
  --af-threshold "$SMALLVAR_AF" --read-models $READ_MODELS --depths $DEPTHS \
  --reuse-dry-run "$OUTDIR/smallvar_reuse.dry_run.txt"
