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

PNG and SVG. The images carry no titles or notes; these are their captions.

**Figure 1** (`figures/fig1_best_f1_depth`). Median F1 over the 14 Samples against Depth, for
SNPs (top) and indels (bottom) with hac (left) and sup (right) reads. Solid lines are Best F1
(the best QUAL threshold); dashed lines with open markers are the Default-PASS score (PASS
records only), for Arms C and D and the AF filter. F1 is on a logit scale, which spreads out
the differences close to 1. Points are nudged sideways so they don't hide each other; each is
at the Depth below it. AF is allele frequency (FORMAT/AF). Arm C + AF filter (0.65) is Arm C's
Clair3 run diploid, with each het call made the ALT when its AF is ≥ 0.65 and REF otherwise;
it isn't one of the Arms.

**Figure 2** (`figures/fig2_pr_curves`). Precision-recall curves over QUAL thresholds at 10x
and 50x, for SNPs (top) and indels (bottom). Each curve pools the 14 Samples, summing their
truth and query counts at each threshold. Markers show each Arm's Default-PASS score, pooled
the same way. Each panel is zoomed to its own range. Arm C + AF filter (0.65) is as in
Figure 1.

**Figure 3** (`figures/fig3_per_sample_best_f1`). Best F1 for every Sample at every Depth,
for SNPs and indels with hac and sup reads. F1 is on a logit scale and each panel has its own
x-axis. A perfect score (F1 = 1) has no place on a logit scale, so perfect scores are drawn in
a column of their own after the dotted line, just past the panel's best imperfect score. The
shaded Samples, *S. enterica* and *V. parahaemolyticus*, carry dnd phosphorothioate systems,
where Dorado's bacterial model is reported to make systematic errors with hac v6.0 reads
([dorado#1599](https://github.com/nanoporetech/dorado/issues/1599)). Arm C + AF filter (0.65)
is as in Figure 1.

**Figure 4** (`figures/fig4_runtime_memory`). Wall time (left) and peak RAM (right) of variant
calling against Depth, both on log scales. Points are medians over the 28 Read sets (14
Samples, hac and sup) and bars show the range. Clair3 (Arms A-C) ran on 8 threads of an AMD
EPYC 9745 and `dorado polish` (Arm D) on one NVIDIA H100 80GB HBM3 with 8 threads; the open
marker is Dorado's timing-only re-run on 8 CPU threads at 50x. Alignment (under 30 s with minimap2)
isn't shown; it is in Table 1. Peak RAM is host memory: Dorado's GPU memory isn't measured.

To re-render them, and Tables 1 and S1, put the (gunzipped) tables in a `results_dir`'s
`tables/` and run only the figure and table rules, as in the update's README.
