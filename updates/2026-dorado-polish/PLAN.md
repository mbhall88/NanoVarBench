# Plan: Dorado polish vs Clair3 (NanoVarBench update)

Agreed 2026-09-30. Terms are defined in [CONTEXT.md](CONTEXT.md); the reasoning behind key decisions is in [docs/adr/](docs/adr/).

## Question

Is `dorado polish --bacteria --vcf` as good as Clair3 for bacterial ONT variant calling, and what does it cost in runtime and resources?

## Venue

- **Post:** on mbhall88.github.io, as a page bundle in the site repo.
  - Cite eLife (`cite "10.7554/eLife.98300"`) near the top, with a "please cite" box.
  - Give the post its own Zenodo DOI.
  - Link to the workflow rather than embedding the scripts.
- **Code:** a self-contained Snakemake workflow in `NanoVarBench/updates/2026-dorado-polish/`. Reuse the Zenodo truth sets, the filter chain and `filter_hets.py`; leave the paper's `workflow/` untouched.
- **NanoVarBench PR:**
  - Add a dated entry to the README's updates section that links to the post and asks people to cite the paper.
  - Add a `CITATION.cff`.
  - Fix the BibTeX licence (it says CC-BY-SA, but Europe PMC lists CC BY). Check against the eLife article page first.
- **Ryan Wick:** thanked in the acknowledgements, and reviews the draft before publishing.
- These docs live in that subdirectory (moved there in #8).

## Data

- **Samples:** all 14 samples in NanoVarBench `config/accessions.csv`.
- **Reads:** simplex hac and sup, basecalled with v4.3.0, from SRA. No rebasecalling ([ADR-0001](docs/adr/0001-v430-reads-and-no-finetuned-clair3.md)).
- **Download integrity:** use kingfisher `-m ena-ascp` (Aspera). Then check every FASTQ against ENA's `fastq_md5` and with `gzip -t`, retrying on failure. kingfisher 0.4.1's own `--check-md5sums` always logs "MD5sum OK" because `check_md5sum()` returns the hexdigest before comparing. Its `ena-ftp` route (aria2c `-x8`) produced crc32-corrupt files when ENA refused parallel connections. On 2026-10-01, Aspera fetched 7 hac read sets (0.7–6.4 GB) in 10–25 min each, all verified. ENA FTP ran at 0.2–0.6 MB/s per stream.
- **Read QC:** apply to hac and sup separately: seqkit keeping reads ≥1000 bp and Q≥10.
- **Subsampling** (ADR-0005): hac and sup are subsampled independently at Depths of 5, 10, 25 and 50x. 100x is dropped because KPC2 and AJ292 can't reach it. For each Read model:
  1. Align the QC'd reads with minimap2 2.31 `lr:hq` to the Mutated reference, keeping primary alignments only (`-F 0x904`).
  2. Run `rasusa aln -c <Depth> -s 20240102`, which caps per-position depth at the Depth evenly across every contig, so plasmids no longer take the depth budget.
  3. Extract the read IDs and pull those reads from the QC'd FASTQ. That FASTQ is the Read set, and every Arm realigns it.
  4. Report the actual per-contig mean depth of each Read set, measured after Arm B's realignment. Plasmids come back somewhat above the Depth (#7 test: 57–75x at a 50x cap), because whole reads bring their supplementary parts back with them.
  - **Caveat for the post:** capped depth is much more uniform than the paper's random genome-wide subsampling (at 5x, the paper's method left about 12% of positions below 2.5x). So low-Depth recall isn't directly comparable with the eLife figures, including Arm A. It's the same for every Arm here, though.
  - **Depth check:** per-contig depth for all 14 Samples, and the tests of alternative approaches, are in `data/depth_check/`.
- **Truth sets:** Zenodo 10867171. Everything is scored on the whole genome, with no repeat-excluded set.

## Arms

See [ADR-0002](docs/adr/0002-four-arm-design.md).

| Arm | Aligner | Caller | Calling model |
|---|---|---|---|
| A | minimap2 2.26 `-aL --cs --MD -x map-ont` | Clair3 1.0.5 | TF v4.3.0 hac/sup |
| B | minimap2 2.31 `-aL --cs --MD -x lr:hq` | Clair3 1.0.5 | TF v4.3.0 hac/sup |
| C | minimap2 2.31 `lr:hq` (same BAM as B) | Clair3 2.0.3 | HKU-converted PyTorch v4.3.0 hac/sup |
| D | minimap2 2.31 `lr:hq` + RG reheader | Dorado polish 2.1.2 | `--bacteria` |

- **Clair3 options (all Clair3 Arms):** `--platform=ont --include_all_ctgs --haploid_precise --no_phasing_for_fa --enable_long_indel`, with `--min_coverage` at its default of 2. 8 threads, on CPU.
- **Clair3 models for Arm C:** download from https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/r1041_e82_400bps_{hac,sup}_v430/ and mount them into the container, since they aren't bundled in `hkubal/clair3:v2.0.3`. Record each file's SHA256.
- **Dorado:**
  - **Command:** `dorado polish --bacteria --vcf --min-depth 2 --any-bam --ignore-read-groups`, with the RG added by `samtools addreplacerg` (`DS:basecall_model=dna_r10.4.1_e8.2_400bps_{hac,sup}@v4.3.0`). See [ADR-0004](docs/adr/0004-rg-header-and-any-bam.md).
  - **Hardware:** full GPUs only (A100 or H100), never MIG slices. `gpu_cuda` jobs need `--qos=gpu`. Full A100s can be booked for days, so accuracy jobs may use any full GPU, since the calls don't depend on device (#7: GPU and CPU gave identical calls, QUAL differing by less than 0.9). Table 1 timings use one pinned GPU model, A100 or H100, to be decided after the #7 timing runs. Also run `--device cpu` at 50x for timing (87 s on 32 threads in #7).
  - **Slurm sizing** (Bunya docs, `UQ-RCC/hpc-docs` guides/Slurm-Tips.md): Dorado used about 1.9 GB RSS and 5–7 s of GPU time per 50x Read set, so request small (4 CPUs, 8 GB) and bundle many Read sets into each GPU job with a Snakemake `group`. `--constraint=cuda80gb --batch=cuda80gb` with `--gres=gpu:1` accepts any full 80 GB A100 or H100. The `debug` QoS (priority 30, ≤1 h, 2 running and 20 submitted) and `short` QoS (priority 20, ≤12 h) cut queue times compared with `gpu` (priority 10).
  - **In Snakemake**, use the Slurm executor plugin's resources, not `ssubmit`. Set `slurm_account` in the profile's `default-resources`. The Dorado rule uses `slurm_partition="gpu_cuda"`, `gres="gpu:h100:1"` (or `a100`), `slurm_extra="--qos=gpu"` (there's no per-rule QoS resource, and `--slurm-qos` is global), and small `runtime` and `mem_mb`. CPU rules use `general` with `normal`. Bundle Read sets with `group:` plus `--group-components`, checking this works with the executor.
  - **Header check:** run one sample from a `dorado aligner` BAM to confirm the calls match.
- **Left out:** Clair3's fine-tuned model, `dorado smallvar`, the newer polishing models, and duplex reads. All are listed as caveats in the post.

## Evaluation

- **Filter chain:** the same NanoVarBench chain for every Arm. It isn't a pass-through: `norm -a` splits Dorado's MNP and complex records into separate records (#7: 4,718 became 4,885). Guard `filter_hets.py` against hets that lack AD or AC.
  1. Reheader contigs.
  2. `filter_hets.py`.
  3. Keep alt genotypes.
  4. `bcftools norm -a -c e -m -`, then `norm -aD`.
  5. Drop indels with |ILEN|>50 and `*` alleles.
  6. `+setGT` to haploid.
- **Scoring:** vcfdist v2.6.4 with `--largest-variant 50 --credit-threshold 1.0 -mx <max QUAL> -b <whole-genome bed>`. Keep all records (QUAL sweep) and never pass `-s`. Dorado caps QUAL at 60, and in #7 the best threshold was QUAL≥0 for every variant type. A flat sweep is a finding in itself, for Ryan's question about whether QUAL separates bad calls. See [ADR-0003](docs/adr/0003-vcfdist-2-6-4.md).
- **Metrics:** Best F1, precision and recall per variant type (SNP, INDEL, ALL), plus F1 Q-score, PR curves, and the Default-PASS score for every Arm (Clair3 sets FILTER too, so Dorado's default is compared with Clair3's).
- **Runtime:** wall time and peak RAM from Snakemake's `benchmark:` output.

## Outputs

**Figures and tables:**
1. **Figure 1:** Best F1 against Depth. SNP and INDEL panels, hac and sup, with one line per Arm.
2. **Figure 2:** PR curves at 50x.
3. **Figure 3:** per-sample Best F1 at 50x, with the dnd samples highlighted (link dorado#1599).
4. **Table 1:** runtime and memory.
5. **Table S1:** interactive per-sample table with actual per-contig depths.
6. **Left out:** the FP/FN characterisation, unless Dorado and Clair3 differ in an interesting way.

**Post outline:**
1. TL;DR
2. Background: Ryan's question, the eLife paper and cite box, what `polish --vcf` is
3. Methods, including the RG/`--any-bam` callout box
4. Results
5. Why not the fine-tuned Clair3 model
6. Caveats
7. Conclusion
8. Acknowledgements

## Software

| Tool | Source |
|---|---|
| Dorado 2.1.2 | https://cdn.oxfordnanoportal.com/software/analysis/dorado-2.1.2-linux-x64.tar.gz (CUDA 12.8; driver ≥525.105). Pre-download the bacterial model with `dorado download`; pass `--models-directory`. |
| Clair3 1.0.5 | `docker://quay.io/mbhall88/clair3:1.0.5` (includes `/opt/models/r1041_e82_400bps_{hac,sup}_v430`) |
| Clair3 2.0.3 | `docker://hkubal/clair3:v2.0.3` |
| minimap2 | mulled with samtools so the aligner can pipe into `samtools sort`: 2.26 + samtools 1.17 (`mulled-v2-66534bcb…:7e6194c8…-0`), 2.31 + samtools 1.23.1 (`mulled-v2-66534bcb…:b411340b…-0`) |
| samtools | `quay.io/biocontainers/samtools:1.24--h9dcdb79_1` |
| bcftools | `quay.io/biocontainers/bcftools:1.24--h118bc1c_2` |
| rasusa | `quay.io/biocontainers/rasusa:5.1.0--hfa8f182_0` |
| seqkit | `quay.io/biocontainers/seqkit:2.14.0--hb192632_0` |
| vcfdist | `docker://timd1/vcfdist:v2.6.4` |

Pin every container by digest in the workflow.

## Storage

- **Compute:** on HPC scratch (`work_dir`, set in the gitignored `config/local.yaml`). Scratch is purged: the earlier lr:hq working directory was lost this way.
- **Git:** commit the small results as each stage finishes: the tables, depth, benchmark and version TSVs, caller info and figures. Filtered VCFs aren't committed (decided 2026-10-02): they'd add about 20 MB of binary files that every rerun rewrites, and the workflow regenerates them.
- **Zenodo:** deposit the final set of filtered VCFs and summaries, from the full run, along with the post DOI.
- **Not kept:** BAMs and reads, which can be regenerated from SRA.

## Out of scope

- The possible medaka model-selection bug on NanoVarBench `main`.
