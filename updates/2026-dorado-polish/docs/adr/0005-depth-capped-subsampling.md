# Cap depth evenly per contig with `rasusa aln`; subsample hac and sup independently; drop 100x

The paper subsampled reads randomly to a genome-wide mean depth (`rasusa reads`). A depth check across all 14 Samples found this undershoots the chromosome wherever non-chromosome bases take up the base budget. Plasmids at up to about 75× chromosome copy number cost ATCC_25922 about 12% (50x gave 43.8x on the chromosome). Unmapped, apparently cross-barcode chimeric bases cost ATCC_10708 about 21%. We now align QC'd reads first, then cap per-position depth with `rasusa aln -c <Depth> -s 20240102` on primary alignments, so every contig sits at the Depth. The selected reads form the Read set that every Arm realigns.

hac and sup are QC'd and subsampled independently, as in the paper. We dropped the planned matched hac/sup read-ID design because only 81 of the first 200 sup read UUIDs for ATCC_25922 appear in the hac run at all. 100x is dropped because KPC2 (about 64x on the chromosome in total) and AJ292 (about 97x) can't reach it.

## Considered Options

- **Keep genome-wide `rasusa reads` (the paper's method).** Comparable with the eLife numbers, but chromosome depth drifts by up to 21% below target.
- **Per-Sample corrected target for `rasusa reads`.** Hits the chromosome mean and keeps natural random depth variation, but plasmids stay at their natural copy number (e.g. about 3,900x).
- **Genome size = chromosome length.** Made the undershoot worse (43.1x at 50x).

## Consequences

- Depth is far more uniform than random subsampling. At 5x, random subsampling leaves about 12% of positions below 2.5x, and capping leaves almost none. So low-Depth recall will look better than the eLife figures, including for Arm A. The post must say this.
- The subsampling alignment (minimap2 2.31 `lr:hq`) is separate from each Arm's own alignment. After realignment, plasmids come out somewhat above the cap (57–75x at 50x) because whole reads bring their supplementary parts back.
- Reads longer than average are slightly favoured (N50 4.8 kb → 5.3 kb in the ATCC_25922 test).
