# dorado aligner vs reheadered minimap2: ATCC_25922__202309 hac 50x, Arm D

**Calls differ materially: no** (more than 5 filtered records differ, or a Best F1 difference above 0.0005).

Both paths start from the same Read set and add the same single `@RG` line. Both run `dorado polish --bacteria --vcf --min-depth 2 --ignore-read-groups` on the same device with the same Calling model.

- **Path 1:** minimap2 2.31 `lr:hq` (`-aL --cs --MD`), then `dorado polish --any-bam`.
- **Path 2:** `dorado aligner` (default `lr:hq`, bundled minimap2), then `dorado polish` without `--any-bam`.

Dorado: path 1 `2.1.2+8b8fc5d`, path 2 `2.1.2+8b8fc5d`. Calling model: path 1 `dna_r10.4.1_e8.2_400bps_polish_bacterial_methylation_v5.0.0`, path 2 `dna_r10.4.1_e8.2_400bps_polish_bacterial_methylation_v5.0.0`. Device: path 1 cuda:all (NVIDIA H100 80GB HBM3, 610.43.02), path 2 cuda:all (NVIDIA H100 80GB HBM3, 610.43.02).

## Calls

| VCF | path 1 records | path 2 records | records that differ | only path 1 | only path 2 | FILTER/GT differs | largest QUAL difference |
|---|---|---|---|---|---|---|---|
| Dorado raw | 4714 | 4714 | 0 | 0 | 0 | 0 | 3.2740 |
| after Filter chain | 4880 | 4880 | 0 | 0 | 0 | 0 | 3.2740 |

A record is identified by CHROM, POS, REF and ALT. It differs if it is in only one VCF or if its FILTER or GT differs. QUAL is compared on records in both VCFs.

Largest QUAL difference: raw 3.2740 (at chromosome:967809 G>A), filtered 3.2740 (at chromosome:967809 G>A). 916 of 4714 raw records have a different QUAL, 17 by more than 1. (For scale: GPU and CPU runs of the same BAM differ by up to 0.87, because the GPU runs the model in half precision, #7.)

## vcfdist (v2.6.4, whole genome)

Best F1 is the QUAL sweep's best row. Default-PASS F1 scores PASS records only with no threshold. Difference is path 2 minus path 1.

| Score | Type | path 1 F1 | path 2 F1 | difference | path 1 FN / FP | path 2 FN / FP |
|---|---|---|---|---|---|---|
| Best F1 | SNP | 0.999117 | 0.999117 | +0.000000 | 7 / 1 | 7 / 1 |
| Best F1 | INDEL | 0.990241 | 0.990241 | +0.000000 | 6 / 1 | 6 / 1 |
| Best F1 | ALL | 0.998465 | 0.998465 | +0.000000 | 13 / 2 | 13 / 2 |
| Default-PASS F1 | SNP | 0.999117 | 0.999117 | +0.000000 | 7 / 1 | 7 / 1 |
| Default-PASS F1 | INDEL | 0.990241 | 0.990241 | +0.000000 | 6 / 1 | 6 / 1 |
| Default-PASS F1 | ALL | 0.998465 | 0.998465 | +0.000000 | 13 / 2 | 13 / 2 |

Largest absolute Best F1 difference: 0.000000. Largest absolute Default-PASS F1 difference: 0.000000.
