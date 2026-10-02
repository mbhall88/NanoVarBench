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
  Clair3's. The figure's note says Depth is a per-position `rasusa aln` cap, not a random
  genome-wide subsample.
- **Figure 2** (`fig2_pr_curves`): precision-recall curves from the QUAL sweep at the Depths in
  `figures.pr_depths` (10 and 50x; a Depth missing from the run is skipped), per variant type
  and Read model, one curve per Arm. The Samples are pooled: truth and query counts are summed
  at each QUAL threshold, so a curve is the PR curve of all their variants together. The
  Default-PASS score is a marker on each curve, pooled the same way. Each panel is zoomed to
  its own range.
- **Figure 3** (`fig3_per_sample_best_f1`): every Sample's Best F1 at every Depth, a dot per
  Arm, with the dnd Samples shaded and their names in red (dorado#1599).
- **Table 1** (`tables/table1_runtime_memory.csv`): runtime and memory, from
  `benchmarks.tsv`. One row per step, tool and Arm(s) at each Depth, with the median and range
  of wall time and peak RSS over the Samples and Read models. The alignment shared by Arms B,
  C and D is one row (not counted three times), and Dorado on CPU (timing only, 50x) is a row
  beside Dorado on GPU. Peak RSS is host memory: GPU memory isn't measured.
- **Table S1** (`tables/table_s1_per_sample.csv`): a row per Sample x Read model x Depth x Arm,
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
`general`, and `call_dorado` goes to one full H100 on `gpu_cuda` or `gpu_sxm`. Large
intermediates go to `work_dir`, and small results go to `results_dir`. Neither is committed:
only the full run's final aggregated tables are, at the end (#15).

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
