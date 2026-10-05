# 2026 update: Dorado polish vs Clair3

A self-contained Snakemake workflow that benchmarks `dorado polish --bacteria --vcf` against
Clair3 on the NanoVarBench data. The paper's workflow at the repository root is untouched.
[PLAN.md](PLAN.md) has the plan, [CONTEXT.md](CONTEXT.md) the glossary, and
[docs/adr/](docs/adr/) the decisions.

The workflow runs all four Arms end to end: A, B and C with Clair3 (#9) and D with Dorado
(#8). For each Sample x Read model x Depth it:

1. downloads the reads and Truth set if they are missing, verifying their MD5s itself;
2. QCs the reads (seqkit, >=1000 bp, Q>=10);
3. caps depth per contig with `rasusa aln` on a primary-only minimap2 2.31 `lr:hq`
   alignment, then pulls the kept reads into the Read set FASTQ (ADR-0005);
4. aligns the Read set for each Arm, and for Dorado adds the `@RG` header line (ADR-0004).
   Arms with the same aligner, version and preset share one BAM: B, C and D use the
   minimap2 2.31 `lr:hq` one, and Arm A has its own minimap2 2.26 `map-ont` BAM;
5. calls with Clair3 (Arms A and B: 1.0.5 with its bundled TF Calling models; Arm C: 2.0.3
   with HKU's converted PyTorch models) or `dorado polish` (Arm D), runs the NanoVarBench
   Filter chain, and scores with vcfdist 2.6.4 (QUAL sweep and PASS only, ADR-0003);
6. writes `results/tables/results.tsv` (one row per Sample x Arm x Read model x Depth x
   variant type x scoring mode), `depth.tsv` (actual per-contig depth), `pr_curves.tsv`
   and `versions.txt`. The scoring modes are `sweep_best` (Best F1, vcfdist's
   `THRESHOLD == BEST` row from the QUAL sweep) and `default_pass` (the Default-PASS score:
   PASS records only, no QUAL threshold), and both are reported for every Arm, since Clair3
   sets FILTER as well as Dorado. `f1_qscore` is −10·log10(1 − F1), capped at Q60 for a
   perfect F1. `pr_curves.tsv` has the QUAL sweep's precision and recall at each threshold,
   keyed and named like `results.tsv`. `dnd_sample` flags the dnd samples (`dnd` in
   [config/samples.tsv](config/samples.tsv)). Each row also carries its Arm's `caller`,
   `caller_version`, `calling_model`, `container` (the pinned digest, for Clair3) and
   `hardware`. They come from the config (the Arm, `dorado.models`, `clair3.calling_models`,
   the pinned `CONTAINERS` and `benchmark_hardware`), not from files the jobs write. The
   jobs run only their tool, so they can be timed cleanly (see Benchmarks);
7. re-runs Dorado with `--device cpu` on the Depths in `dorado.cpu_run` (50x) for the
   runtime comparison (#12). These runs are timing only: they aren't scored and never appear
   in `results.tsv`, and `results/calls/.../<arm>.cpu_vs_main.tsv` checks their calls match
   the main run's. Clair3 only runs on CPU, so it has no re-run;
8. writes `results/tables/benchmarks.tsv` (see Benchmarks below);
9. draws the post's figures and tables from those aggregated tables alone (see Figures and
   tables below).

Arms are config entities (`arms:` in [config/config.yaml](config/config.yaml)); `run:`
picks the Samples, Read models, Depths and Arms.

**Results:** the full run's aggregated tables, the post's figures and tables, and the configs
that made them are in [final/](final/) (#15).

## Benchmarks (#12)

`results/tables/benchmarks.tsv` is Table 1: one row per Sample x Read model x Depth x Arm x
step x `timing_only`, built by the `benchmarks` rule from Snakemake's `benchmark:` files
(`results/benchmarks/`). The timed rules (`align`, `call_dorado`, `call_dorado_cpu`,
`call_clair3`) only run their tool, so `wall_time_s` is the tool's: no version calls, weight
hashing, `nvidia-smi` or `hostname` run in a job, and there is no separate timing of the
command.

- **Steps.** `align` is the Arm's alignment (minimap2 with `samtools sort` and `index`). Arms
  with the same aligner, version and preset share one BAM, so B, C and D carry the same
  alignment row and Arm A has its own; don't add alignment rows across Arms. `call` is the
  calling step alone: Clair3 (Arms A-C, 8 threads) or `dorado polish` (Arm D, on the GPU).
  Everything else is left out: read QC, subsampling, the `@RG` reheader, the Filter chain and
  scoring. Arm D at the Depths in `dorado.cpu_run` (50x) has a second `call` row with
  `timing_only = true`: `dorado polish --device cpu` on 8 threads, the same as Clair3. It
  isn't scored and isn't in `results.tsv`.
- **Columns.** The keys, `step`, `timing_only`, `tool`, `tool_version`, `alignment`, `device`
  (`gpu` or `cpu`), `hardware` (the GPU model, or `CPU: <model>`), `threads`, then
  `wall_time_s`, `max_rss_mb`, `max_vms_mb`, `max_uss_mb`, `max_pss_mb`, `cpu_time_s`,
  `mean_cpu_pct`, `io_in_mb` and `io_out_mb`. Times are in seconds and memory in MB, as in
  Snakemake's benchmark files. `max_rss_mb` is the whole job's.
- **Hardware is a config value.** `hardware` is `benchmark_hardware` in
  [config/config.yaml](config/config.yaml) for the row's `device` (`gpu: "NVIDIA H100 80GB
  HBM3"`, `cpu: "AMD EPYC 9745"`), and `threads` is the rule's configured threads. Nothing
  is measured in a job, so there are no host or driver columns. The hardware was confirmed
  once, outside the workflow: `nvidia-smi` on bun118 and bun119 gave NVIDIA H100 80GB HBM3
  with driver 610.43.02, and the `epyc5` nodes are AMD EPYC 9745. The Bunya profile pins it.
- **One piece of hardware.** `call_dorado` is pinned to one H100 variant, the **H100 SXM**
  (`gpu_sxm`, "NVIDIA H100 80GB HBM3"), never an A100 or a MIG slice. `gpu_cuda`'s H100s are
  PCIe, a different part with different speed, so letting `call_dorado` pick either would mix
  hardware. SXM was also the partition with free H100s when the benchmark was run. The timed
  CPU steps (alignment, Clair3 and the Dorado CPU re-run) are pinned to the `epyc5` nodes
  (AMD EPYC 9745). The profile's partition and constraint are what guarantee this; the
  workflow doesn't check it. If you change them, change `benchmark_hardware` to match. Nodes
  are shared, so CPU timings include whatever else ran on the node.
- **Re-timing.** To re-time the calling steps, force just those rules (and `align` to time
  the alignment too). The targets go before `--config`, which takes the rest of the line:

  ```sh
  snakemake -s workflow/Snakefile --workflow-profile profiles/bunya --configfile config/local.yaml \
      --forcerun align call_dorado call_dorado_cpu call_clair3 \
      --config 'run={samples: [ATCC_25922__202309], read_models: [hac, sup], depths: [5, 10, 25, 50], arms: [A, B, C, D]}'
  ```

## Figures and tables (#14)

The `fig*` and `table*` rules ([workflow/rules/figures.smk](workflow/rules/figures.smk)) read
only `results/tables/` (`results.tsv`, `pr_curves.tsv`, `depth.tsv`, `benchmarks.tsv`), never a
VCF, BAM or score file, so they can be rebuilt from the aggregated tables without `work_dir`.
They run in a conda env ([workflow/envs/plot.yaml](workflow/envs/plot.yaml): pandas and
matplotlib) and are part of `rule all`. Arms are labelled as in CONTEXT.md ("Arm D (Dorado)",
from `figures.arm_labels` in the config). Figures are written as PNG and SVG to
`results/figures/`, tables as CSV to `results/tables/`. Nothing here is committed with the
code: the rendered outputs go in once, with the final aggregated tables (#15).

- **Figure 1** (`fig1_best_f1_depth`): Best F1 against Depth, with SNP and INDEL rows and hac
  and sup columns. Solid lines are each Arm's Best F1 (`sweep_best`), medians over the
  Samples. Dashed lines are the Default-PASS score of Arms C and D
  (`figures.default_pass_arms`), so Arm D's default is shown apart from its Best F1 and beside
  Clair3's. F1 is on a logit scale, which spreads out the differences close to 1 that a
  linear axis squashes together. The figure's note says Depth is a per-position `rasusa aln`
  cap, not a random genome-wide subsample.
- **Figure 2** (`fig2_pr_curves`): precision-recall curves from the QUAL sweep at the Depths in
  `figures.pr_depths` (10 and 50x; a Depth missing from the run is skipped), per variant type
  and Read model, one curve per Arm. The Samples are pooled: truth and query counts are summed
  at each QUAL threshold, so a curve is the PR curve of all their variants together. The
  Default-PASS score is a marker on each curve, pooled the same way. Each panel is zoomed to
  its own range.
- **Figure 3** (`fig3_per_sample_best_f1`): every Sample's Best F1 at every Depth, a dot per
  Arm, with the dnd Samples shaded and their names in red (dorado#1599). F1 is on a logit
  scale. A perfect score has no logit, so perfect scores are drawn in a column of their own
  after a dotted line, just past the panel's best imperfect score: the axis is broken there.
- **The AF filter series** (#28): with the AF filter analysis enabled (see below), Figures 1-3
  and Table S1 gain a series for Clair3 with the AF filter, labelled "Arm C + AF filter (0.65),
  extra analysis" since it is not one of the Arms. The Arm and threshold are
  `figures.af_filter` in the config (`arm: C`, `threshold: 0.65`); the threshold must be one of
  `clair3_af_filter.thresholds`, and the series is drawn when `clair3_af_filter.arms` includes
  the Arm. Figure 1 gives it a solid Best F1 line and a dashed Default-PASS line like Arms C and
  D, Figure 2 its PR curve (from `tables/clair3_af_filter_pr_curves.tsv`) and Default-PASS
  point, and Figure 3 a dot per Sample. It is bluish green with a plus marker, and each
  figure's note says it is an extra analysis. In Table S1 it is one more row per Read set,
  after the Arms' rows, with the Arm column "C + AF filter (0.65), extra analysis, not an Arm".
  The series is read from `tables/clair3_af_filter.tsv`, so the figure and table rules need
  that table (and the PR curves table) next to the main ones. With the analysis off, the
  figures and tables are as without it.
- **Table 1** (`tables/table1_runtime_memory.csv`): runtime and memory, from
  `benchmarks.tsv`. One row per step, tool and Arm(s) at each Depth, with the median and range
  of wall time and peak RSS over the Samples and Read models. The alignment shared by Arms B,
  C and D is one row (not counted three times), and Dorado on CPU (timing only, 50x) is a row
  beside Dorado on GPU. Peak RSS is host memory: GPU memory isn't measured.
- **Table S1** (`tables/table_s1_per_sample.csv`): a row per Sample x Read model x Depth x Arm
  (and the AF filter series, when it is on),
  with the Read set's actual depth, overall and per contig (from `depth.tsv`), and for SNP
  and INDEL the Best F1 (with the QUAL threshold it was reached at) and the Default-PASS score,
  each with precision and recall. Its headers are readable, for the site's interactive
  `csv-table`.

To render the figures and tables from a copy of the aggregated tables, put them in a
`results_dir`'s `tables/` and ask for just these rules. Snakemake doesn't re-run the upstream
jobs when the tables already exist:

```sh
snakemake results/figures/fig1_best_f1_depth.png results/figures/fig2_pr_curves.png \
    results/figures/fig3_per_sample_best_f1.png results/tables/table1_runtime_memory.csv \
    results/tables/table_s1_per_sample.csv -s workflow/Snakefile --cores 1 \
    --software-deployment-method conda --allowed-rules fig1_best_f1_depth fig2_pr_curves \
    fig3_per_sample_best_f1 table1_runtime_memory table_s1_per_sample
```

With the AF filter series on, also copy `clair3_af_filter.tsv` and
`clair3_af_filter_pr_curves.tsv` into `tables/`, or set `clair3_af_filter.arms: []` to render
without it. Make the PR curves table, if the run lacks it, from the AF filter's score files
in `work_dir`: `snakemake results/tables/clair3_af_filter_pr_curves.tsv --allowed-rules
clair3_af_filter_pr_curves`.

## Clair3 (Arms A, B and C)

All three use the same options as the paper's rule, for haploid bacterial calling:
`--platform=ont --include_all_ctgs --haploid_precise --no_phasing_for_fa --enable_long_indel`
on 8 CPU threads. `--min_coverage` stays at its default of 2 in both versions, which matches
Dorado's `--min-depth 2` (ADR-0002). The Calling model follows the Read model (hac reads use
`r1041_e82_400bps_hac_v430`), and the bacteria Fine-tuned model is never used (ADR-0001).

- **Arms A and B** run `quay.io/mbhall88/clair3:1.0.5`, pinned by digest, with the TF models
  bundled at `/opt/models`.
- **Arm C** runs `hkubal/clair3:v2.0.3`, pinned by digest, with HKU's converted PyTorch
  models from <https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/>. The
  `download_clair3_model` rule fetches them to `clair3.models_dir` and checks each file against
  the SHA256 pinned in `config/config.yaml`. HKU doesn't publish checksums, so the pins are
  from a download on 2026-10-01. Apptainer bind-mounts the directory (the Bunya profile binds
  `/scratch`), so it must be somewhere Apptainer can see.

## dorado aligner check (#13, optional)

Arm D feeds `dorado polish` a minimap2 BAM with a hand-added `@RG` line and `--any-bam`
(ADR-0004). To check that this doesn't change Dorado's calls, list Read sets under
`dorado_aligner_check` in the config (empty = off, and it isn't an Arm: it adds nothing to
`results.tsv`). For each one the check aligns the Read set with `dorado aligner` 2.1.2
(default preset `lr:hq`, its bundled minimap2), sorts and indexes the BAM, adds the same
`@RG` line, and runs `dorado polish` **without** `--any-bam`, with the other flags unchanged.
Nothing extra is needed to keep the `@PG` line polish looks for (`ID:aligner`): `samtools sort`
and `samtools reheader` both keep it, and `dorado_align_sort` fails if it is missing. Both
paths' calls go through the Filter chain and vcfdist, and `compare_aligner_paths` writes
`results/dorado_aligner_check/<sample>/<read_model>/<depth>x/<arm>.dorado_aligner_vs_minimap2.{md,tsv}`:
the records that differ (CHROM, POS, REF, ALT, FILTER, GT), the largest QUAL difference and
the Best F1 and Default-PASS F1 differences for SNP, INDEL and ALL, flagged "materially
different" above the `material` limits in the config. Its `dorado polish` job needs the same
GPU as `call_dorado`, so the Bunya profile gives it the same resources.

```sh
# only this check, not the full table (the target goes before --config, which takes many values)
snakemake dorado_aligner_check -s workflow/Snakefile --workflow-profile profiles/bunya \
    --configfile config/local.yaml \
    --config 'dorado_aligner_check={samples: [ATCC_25922__202309], read_models: [hac, sup], depths: [50]}'
```

(With the check configured, a plain `snakemake` run does it too.) Seam 1 runs it on the
fixture.

## AF filter (#20, optional)

`--haploid_precise` drops every het call, and in repeats Clair3 calls the true ALT as a
low-confidence het when reads from the other copy mix in (#6). The AF filter tests whether
keeping those hets and resolving them by allele frequency closes Clair3's SNP gap to Dorado.
It isn't an Arm and changes nothing in `results.tsv`. For each Clair3 Arm listed under
`clair3_af_filter.arms` (empty = off), on every Read set in `run`:

1. Clair3 runs on the Arm's alignment with the Arm's version, Calling model and options, minus
   `--haploid_precise`, so it writes diploid genotypes (`call_clair3_af`).
2. For each AF threshold in `clair3_af_filter.thresholds` (0.5 to 0.8 in steps of 0.05),
   `af_filter` ([workflow/scripts/af_filter.py](workflow/scripts/af_filter.py)) makes each het
   (0/1, or 1/2) homozygous for the ALT in its genotype with the highest `FORMAT/AF` when that
   AF is at least the threshold, and homozygous REF otherwise. Hom calls are left alone.
3. The result goes through the Filter chain and vcfdist with the Arms' settings (QUAL sweep and
   PASS only).

It uses Clair3's diploid output, not `--haploid_sensitive`. Both recover the repeat SNPs (#6),
but `--haploid_sensitive` writes 0/1 and 1/1 alike as `1` and drops multi-allelic (1/2) sites
(Clair3's `CallVariants.py`), so the filter couldn't tell a het from a hom call.

`clair3_af_filter_tables` writes two tables, and `clair3_af_filter_pr_curves` a third, the QUAL
sweep's PR curves for the one Arm and threshold the figures draw (`figures.af_filter`, #28; the
other thresholds' curves aren't tabulated):

- `results/tables/clair3_af_filter.tsv`: one row per Sample x Read model x Depth x Arm x
  `calls` x `af_threshold` x variant type x scoring mode, with the columns of `results.tsv`'s
  scores. `calls` is `af_filter` at each threshold, or `main` for the Arm's own
  `--haploid_precise` calls and for `reference_arm` (Arm D), copied from their scores so the
  table stands alone. `scoring_mode` is `sweep_best` (Best F1) or `default_pass`.
- `results/tables/clair3_af_filter_summary.tsv`: the same over the Samples. For each Arm,
  `calls`, threshold, Read model, Depth, variant type and scoring mode: the median, min and
  max F1, the summed FN and FP, and for the `af_filter` rows the median F1 change from the
  Arm's own calls (with how many Samples are better, tied and worse), `n_best` (Samples where
  this threshold is the best one), `median_loss_vs_best_threshold` (what using this one
  threshold costs against each Sample's best), and `gap_to_reference_closed`: the share of the
  gap in median F1 between the Arm's own calls and the reference Arm that this closes, when
  the reference Arm is ahead.
- `results/tables/clair3_af_filter_pr_curves.tsv`: `pr_curves.tsv`'s layout (one row per
  Sample x Read model x Depth x QUAL threshold x variant type) plus `af_threshold`, for the
  figures' Arm and threshold, so Figure 2 can draw the series.

To run it on a finished run's Read sets without re-making them, set
`clair3_af_filter.input_work_dir` to that run's `work_dir`, and give this run a new
`work_dir` and `results_dir`. Its BAMs, Mutated references, Truth sets and the Arms' scores are
then read from there, so only the analysis' own jobs run. Check with `-n` first:

```sh
# my_af_filter.yaml sets work_dir, results_dir, run, and
#   clair3_af_filter: {arms: [A, C], input_work_dir: /path/to/full/run/work}
snakemake clair3_af_filter -s workflow/Snakefile --workflow-profile profiles/bunya \
    --configfile config/local.yaml my_af_filter.yaml -n
```

The Bunya profile groups each threshold's `af_filter`, `filter_calls_af` and two
`vcfdist_af` jobs into one Slurm job, since each takes seconds. `call_clair3_af` isn't timed,
so it isn't pinned to the timed CPU nodes. Seam 1 runs it for Arm C on the fixture, then
checks that a fresh `work_dir` reading the fixture run's schedules only these jobs.

## smallvar pilot (#29, optional)

`dorado smallvar` is Dorado's diploid small-variant caller. This pilot asks whether it is worth
a place in the post. Dorado 2.1.2 ships smallvar models only for hac v5.2.0 and v6.0.0 reads
(`dna_r10.4.1_e8.2_400bps_hac@v{5.2.0,6.0.0}_smallvar@v1.0`): none for v4.3.0, which our reads
were basecalled with, and none for sup. The pilot forces in `smallvar_pilot.model` (hac v6.0.0)
with `--model-override`, which skips Dorado's compatibility check (Dorado warns "Variant calling
model is not compatible with the input BAM. This may produce inferior results"). **Every result
is a Basecall-model mismatch**: of version on hac reads, and of tier and version on sup reads.
The tables label each smallvar row (`basecall_model_mismatch`: `version` or
`tier_and_version`). It isn't an Arm and changes nothing in `results.tsv` or `benchmarks.tsv`.
For each Dorado Arm listed under `smallvar_pilot.arms` (empty = off), on every Read set in `run`:

1. `call_dorado_smallvar` runs `dorado smallvar` on the Arm's RG-reheadered BAM (ADR-0004), with
   the whole genome (`genome.bed`) as `--hemizygous-regions`, so it makes haploid calls. It
   needs the hidden `--any-bam` flag, as polish does, and gets `--min-depth 2` (ADR-0002) and
   `--ignore-read-groups`. The rest are Dorado's defaults, including `--pass-qual-filter 3`.
   Only the tool runs in the job, which is timed.
2. The calls go through the Filter chain and vcfdist with the Arms' settings (QUAL sweep and
   PASS only).

**No AF filter.** smallvar's VCF has only `GT` and `GQ` in FORMAT. INFO declares `DP` but
leaves every record's INFO empty, and `--gvcf` adds only reference blocks (`END`, `LEN`).
With no AF and no allele depths, a het can't be resolved by AF, so the pilot has no diploid +
AF filter mode. In the diploid mode the Filter chain's `filter_hets.py` would turn every het
into REF (it has no AD or AC to go on), which is what `--hemizygous-regions` gives anyway. On
a 200 kb test (ATCC_25922 hac 50x) the diploid run's hom calls were exactly the haploid run's
records, and its 3 hets were all false.

`smallvar_pilot_tables` writes two tables:

- `results/tables/smallvar_pilot.tsv`: one row per Sample x Read model x Depth x `calls` x Arm
  x variant type x scoring mode, with the columns of `results.tsv`'s scores. `calls` is
  `smallvar_haploid`, `main` for the Arms in `smallvar_pilot.compare_arms` (C and D), copied
  from their scores, or `af_filter` for Arm C's AF filter (#20) at `af_filter.threshold`
  (0.65), copied from its scores. `calling_model` and `basecall_model_mismatch` say which model
  made each row's calls and how it differs from the reads' basecall model (`none` for the
  Arms).
- `results/tables/smallvar_pilot_runtime.tsv`: one row per smallvar job, with the job's wall
  time and memory and, as in `benchmarks.tsv`, the device, hardware and threads from the config.
  The Bunya profile runs it on the same H100 SXM as `call_dorado`.

To run it on a finished run's Read sets without re-making them, set
`smallvar_pilot.input_work_dir` to that run's `work_dir` (BAMs, Mutated references, Truth sets
and the Arms' scores) and `smallvar_pilot.af_filter.work_dir` to the AF filter run's (its
scores), and give this run a new `work_dir` and `results_dir`. Check with `-n` first:

```sh
# my_smallvar.yaml sets work_dir, results_dir, run, and
#   smallvar_pilot: {arms: [D], input_work_dir: /path/to/full/run/work,
#                    af_filter: {arm: C, threshold: 0.65, work_dir: /path/to/af/run/work}}
snakemake smallvar_pilot -s workflow/Snakefile --workflow-profile profiles/bunya \
    --configfile config/local.yaml my_smallvar.yaml -n
```

`download_smallvar_model` fetches the model into `dorado.models_dir`. The Bunya profile groups
each Read set's Filter chain and two vcfdist runs into one Slurm job. Seam 1 runs the pilot for
Arm D on the fixture, on CPU.

## Requirements

Snakemake 9 with the Slurm executor plugin, Apptainer, conda/mamba, and the
[Dorado 2.1.2](https://cdn.oxfordnanoportal.com/software/analysis/dorado-2.1.2-linux-x64.tar.gz)
binary. Every other tool runs from a container pinned by digest
([workflow/rules/common.smk](workflow/rules/common.smk)), except the Filter chain's
bcftools + cyvcf2, which use a conda env ([workflow/envs/filter.yaml](workflow/envs/filter.yaml)),
and the figure and table scripts ([workflow/envs/plot.yaml](workflow/envs/plot.yaml)).

## Run

Paths in [config/config.yaml](config/config.yaml) default to directories inside this one
(`data/`, `resources/`, `work/`, all gitignored). To point them elsewhere, copy
[config/local.example.yaml](config/local.example.yaml) to `config/local.yaml` (gitignored)
and pass it with `--configfile`. From this directory:

```sh
conda activate snakemake   # activate it; don't just put its bin/ on PATH
snakemake -s workflow/Snakefile --workflow-profile profiles/bunya --configfile config/local.yaml
```

If Snakemake's env is on PATH but not activated, conda's stacked activation leaves it ahead
of the rule's env, so the Filter chain runs the wrong Python (`No module named 'cyvcf2'`).

`profiles/bunya` submits to Slurm on Bunya: CPU rules (including `call_dorado_cpu`) go to
`general`, with the timed ones on the `epyc5` nodes, and `call_dorado` goes to one full H100
SXM on `gpu_sxm` (see Benchmarks). Large intermediates go to `work_dir`, and small results go
to `results_dir`. Neither is committed: the full run's final aggregated tables and figures
are in [final/](final/) (#15).

## Test

Seam 1 runs the whole workflow on a 40 kb fixture (tests/fixture/) for Arms A-D, on a hac
and a sup Read set at every Depth in the config (5, 10, 25 and 50x), with Dorado on CPU at
the largest Depth. It checks the results table, that each Arm's rows carry the config's
caller version, Calling model (for Dorado, the model `dorado polish` logged as resolved),
container digest and hardware, that the chromosome's actual depth is within 5% of
the Depth for every Read set, that the fixture's known variants are true positives at
50x, that `benchmarks.tsv` exists and is sane, and that the figures and tables exist (PNG and
SVG figures, and Table 1 and Table S1 with the rows they should have). It takes several minutes, plus a download
of the Clair3 images (about 6 GB) and HKU's models on a first run:

```sh
DORADO=/path/to/dorado DORADO_MODELS_DIR=/path/to/models tests/seam1.sh
# optionally CLAIR3_MODELS_DIR=/path/to/clair3_models to reuse downloaded models
# optionally READ_MODELS="hac" DEPTHS="25 50" for a quicker run on a subset
```

The fixture's reads (`reads.hac.fastq.gz`, `reads.sup.fastq.gz`, about 1.8 MB each) are the
real ATCC_25922 reads for the window, thinned to an even ~52x; `tests/fixture/make_fixture.sh`
records how they were made.

Seam 1 also runs the AF filter for Arm C at `AF_THRESHOLDS` (default 0.5, 0.65 and 0.8), and
`tests/check_af_filter.py` checks it: Clair3 ran without a haploid mode, every het was resolved
by its AF and every hom call left alone, the tables' main rows match `results.tsv`, the
summary follows the per-Sample table, the AF filter recovers repeat SNPs that
`--haploid_precise` drops at 50x, and reusing the run's `work_dir` schedules only the AF
filter's jobs. The run's figures and Table S1 then have the AF filter series at
`AF_FIGURE_THRESHOLD` (default 0.65), and seam 1 renders them again from the same aggregated
tables with the analysis off. `tests/check_af_figures.py` checks that Table S1 has one labelled
AF filter row per Read set with the AF table's scores, that the Arms' rows are the same on and
off, that the PR curves table is the Arm's sweep at the threshold, and that each figure names
the series when on and never when off.

Seam 1 also runs the smallvar pilot for Arm D, beside Arm C's AF filter at the lowest
`AF_THRESHOLDS`, and `tests/check_smallvar_pilot.py` checks it: smallvar used the overridden
model and warned of the mismatch, its calls are haploid with only GT and GQ (no AF or allele
depths), every smallvar row carries its Basecall-model mismatch, the Arm and AF filter rows
match `results.tsv` and `clair3_af_filter.tsv`, neither `results.tsv` nor `benchmarks.tsv`
has smallvar rows, the runtime table has a row per Read set, and reusing the run's `work_dir`
schedules only the pilot's jobs.

Seam 2 tests the aggregation on its own. It writes hand-written vcfdist summaries and
precision-recall curves, `samtools coverage` output and benchmark files
(tests/seam2/) where a run leaves them, for all 14 Samples x hac x 25 and 50x x Arms A-D. It
then runs only the `aggregate` and `benchmarks` rules and checks the output tables:

- `sweep_best` rows are the sweep's BEST row, and `default_pass` rows are the PASS-only
  run's unthresholded NONE row, kept separate for every Arm;
- F1 Q-score, including the Q60 cap;
- actual depth for each Read set;
- the dnd flag for all 14 Samples;
- each PR curve's keys, and that every Best F1 row is a point on its curve;
- `results.tsv`'s columns (no `calling_model_sha256`) and each Arm's config-derived caller,
  version, Calling model, container and hardware;
- `benchmarks.tsv`: its keys and units, the config's hardware and threads on every row, the
  worked Read set's literal timings, separate alignment rows (shared by B, C and D), the
  Dorado CPU rows present at 50x only, and no CPU run in `results.tsv`.

It takes a few seconds and needs only Snakemake (no GPU, Slurm, containers or downloads):

```sh
tests/seam2.sh   # OUTDIR=dir to keep the output somewhere other than a new temp dir
```
