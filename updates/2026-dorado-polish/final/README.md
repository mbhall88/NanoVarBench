# Final results (#15)

The aggregated tables, figures and post tables from the full run, made by this workflow and
committed once. Per-run outputs (VCFs, vcfdist summaries, BAMs, reads) aren't in git; the
filtered VCFs and vcfdist summaries go to Zenodo (#17).

If you use these results, please cite the NanoVarBench paper:
Hall MB et al. (2024) Benchmarking reveals superiority of deep learning variant callers on
bacterial nanopore sequence data. *eLife* 13:RP98300. https://doi.org/10.7554/eLife.98300

## What was run

| Run | Config | Read sets | Tables |
|---|---|---|---|
| Full run: Arms A-D | [config/full_run.yaml](config/full_run.yaml) | 14 Samples x hac/sup x 5, 10, 25, 50x (112) | `results`, `pr_curves`, `depth`, `benchmarks`, `versions` |
| AF filter analysis (#20), Arms A and C | [config/af_filter.yaml](config/af_filter.yaml) | the same 112 | `clair3_af_filter*` |
| smallvar pilot (#29), Arm D's alignments | [config/smallvar_pilot.yaml](config/smallvar_pilot.yaml) | ATCC_25922, KPC2, AMtb_1 x hac/sup x 10, 50x (12) | `smallvar_pilot*` |

Each config was used with `config/local.yaml` (work and results directories, the Dorado
binary and models) and `profiles/bunya`. Every job in all three runs succeeded. Dorado ran on
one full H100 SXM ("NVIDIA H100 80GB HBM3"); the timed CPU steps ran on AMD EPYC 9745 nodes.
Tool versions and pinned container digests are in [tables/versions.txt](tables/versions.txt).

## Tables

Columns are described in the update's [README](../README.md). Files over 1 MB are gzipped.

| File | Contents |
|---|---|
| `tables/results.tsv` | Best F1 (`sweep_best`) and Default-PASS (`default_pass`) scores: one row per Sample x Read model x Depth x Arm x variant type x scoring mode |
| `tables/pr_curves.tsv.gz` | precision and recall at every QUAL threshold of the sweep |
| `tables/depth.tsv` | each Read set's actual depth, overall and per contig |
| `tables/benchmarks.tsv` | wall time, CPU time and peak RSS of every timed step |
| `tables/versions.txt` | tool versions and container digests |
| `tables/clair3_af_filter.tsv.gz` | the AF filter's scores per Sample, Arm and AF threshold (0.50-0.80), beside the Arm's `--haploid_precise` scores and Arm D's |
| `tables/clair3_af_filter_summary.tsv` | its medians over the Samples per threshold |
| `tables/clair3_af_filter_pr_curves.tsv.gz` | its PR curves, for Arm C at 0.65 (Figure 2) |
| `tables/smallvar_pilot.tsv` | smallvar's scores beside Arms C and D and Arm C + AF filter (0.65), on the pilot's Read sets; every smallvar row is a Basecall-model mismatch |
| `tables/smallvar_pilot_runtime.tsv` | smallvar's wall time and peak RSS on the H100 |
| `tables/table1_runtime_memory.csv` | Table 1: runtime and memory |
| `tables/table_s1_per_sample.csv` | Table S1: per-Sample results with actual depths, for the post's interactive table |

## Figures

PNG and SVG, with the Arm C + AF filter (0.65) series from #28:

- `figures/fig1_best_f1_depth`: Best F1 against Depth, with Default-PASS dashed
- `figures/fig2_pr_curves`: PR curves at 10x and 50x, the Samples pooled
- `figures/fig3_per_sample_best_f1`: per-Sample Best F1 at every Depth, the dnd Samples shaded

To re-render them, and Tables 1 and S1, put the (gunzipped) tables in a `results_dir`'s
`tables/` and run only the figure and table rules, as in the update's README.
