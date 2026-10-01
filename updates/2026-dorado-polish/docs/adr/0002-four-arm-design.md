# Rerun the old Clair3 configurations as Arms A and B instead of quoting published numbers

The obvious comparison is just Dorado against the latest Clair3. We add two baseline Arms: A (minimap2 2.26 `map-ont`, Clair3 1.0.5) reproduces the eLife paper, and B (minimap2 2.31 `lr:hq`, Clair3 1.0.5) reproduces the lr:hq post. That way A→B→C→D changes one thing at a time (preset, then Clair3 version and model conversion, then caller). We rerun them rather than quoting published numbers because the original VCFs were lost in a scratch purge, and new subsamples plus a new vcfdist version (ADR-0003) would make quoted numbers incomparable.

The old tool versions in Arms A and B are pinned on purpose. Don't "upgrade" them to latest.

## Consequences

- B→C changes the Clair3 version and the TF→PyTorch model conversion together. v2 can't load the TF weights, so the two effects can't be separated.
- Arm D runs with `--min-depth 2`, not Dorado's default of 0, to match Clair3's `--min_coverage` default of 2, which the paper also used. Otherwise the Arms would differ in their depth threshold as well as their caller.
