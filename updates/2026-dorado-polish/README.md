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
   and `versions.txt`. Each Arm's `calls/.../<arm>.caller_info.tsv` records the caller
   version, Calling model, the SHA256 of its files, and (for Clair3) the container digest,
   and `results.tsv` carries them too;
7. re-runs Dorado with `--device cpu` on the Depths in `dorado.cpu_run` (50x) for the
   runtime comparison (#12). These runs are timing only: they aren't scored, and
   `results/calls/.../<arm>.cpu_vs_main.tsv` checks their calls match the main run's.
   Clair3 only runs on CPU, so it has no re-run.

Arms are config entities (`arms:` in [config/config.yaml](config/config.yaml)); `run:`
picks the Samples, Read models, Depths and Arms.

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

## Requirements

Snakemake 9 with the Slurm executor plugin, Apptainer, conda/mamba, and the
[Dorado 2.1.2](https://cdn.oxfordnanoportal.com/software/analysis/dorado-2.1.2-linux-x64.tar.gz)
binary. Every other tool runs from a container pinned by digest
([workflow/rules/common.smk](workflow/rules/common.smk)), except the Filter chain's
bcftools + cyvcf2, which use a conda env ([workflow/envs/filter.yaml](workflow/envs/filter.yaml)).

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
intermediates go to `work_dir`, and small results go to `results_dir`, which is committed.

## Test

Seam 1 runs the whole workflow on a 40 kb fixture (tests/fixture/) for Arms A-D, with
Dorado on CPU, then checks the results table, that each Arm records its versions, Calling
model checksums and container digest, and that the fixture's known variants are true
positives. It takes a few minutes, plus a download of the Clair3 images (about 6 GB) and
HKU's models on a first run:

```sh
DORADO=/path/to/dorado DORADO_MODELS_DIR=/path/to/models tests/seam1.sh
# optionally CLAIR3_MODELS_DIR=/path/to/clair3_models to reuse downloaded models
```
