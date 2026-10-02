# Feed Dorado a reheadered minimap2 BAM via `--any-bam`

`dorado polish` requires an `@PG` line from `dorado aligner` and an `@RG` line carrying `DS:basecall_model=...`. The SRA FASTQs carry neither. We align with minimap2 2.31 `lr:hq`, add the read group with `samtools addreplacerg` (`basecall_model=dna_r10.4.1_e8.2_400bps_{hac,sup}@v4.3.0`), and run polish with the hidden `--any-bam` flag plus `--ignore-read-groups`. That way Clair3 (Arm C) and Dorado (Arm D) see byte-identical alignments. One sample is also run from a `dorado aligner` BAM to confirm that the calls match.

## Considered Options

- **`dorado aligner` for all Arms.** This avoids the hidden flag, but still needs the RG added by hand, and it ties the Clair3 Arms to Dorado's bundled minimap2.
- **Rebasecall to BAM.** Rejected in ADR-0001.

## Confirmed on real data (#7, ATCC_25922 hac 50x)

Dorado 2.1.2 accepts the reheadered BAM. The minimum it needs is `--any-bam` plus a single `@RG` header line whose DS contains `basecall_model=dna_r10.4.1_e8.2_400bps_hac@v4.3.0`. Tagging each read with `RG:Z` isn't needed, and with only one read group neither is `--ignore-read-groups`: both gave identical VCF records. We keep `--ignore-read-groups` as a harmless safeguard. Without the workaround, Dorado fails with "Input BAM file was not aligned using Dorado." or "Input BAM file has no basecaller models listed in the header."

## Confirmed against `dorado aligner` (#13, ATCC_25922 at 50x, hac and sup)

Aligning the same Read set with `dorado aligner` 2.1.2 (default `lr:hq`, bundled minimap2 2.31), sorting it, adding the same `@RG` line and running `dorado polish` without `--any-bam` gives the same calls as the reheadered minimap2 path. In the filtered VCFs 0 of about 4,880 records differ in CHROM, POS, REF, ALT, FILTER or GT, and the Best F1 and Default-PASS F1 differences are 0 for SNP, INDEL and ALL, for both Read models. QUAL is not identical: it differs on 8-19% of records, by at most 3.3, so a QUAL threshold could move a record. `dorado aligner` writes the `@PG ID:aligner` line polish looks for, and `samtools sort` and `samtools reheader` keep it, so nothing extra is needed beyond sorting (its stdout BAM is unsorted). The reports are in `results/dorado_aligner_check/`.
