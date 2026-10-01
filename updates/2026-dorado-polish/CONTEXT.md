# Dorado polish vs Clair3 (NanoVarBench update)

A follow-up to the NanoVarBench eLife paper (doi:10.7554/eLife.98300) asking whether `dorado polish --vcf` can replace Clair3 for variant calling from bacterial ONT reads against a haploid reference.

## Language

### Data

**Sample**:
One of the 14 NanoVarBench bacterial isolates, each with its own assembly and ONT reads.
_Avoid_: isolate, genome, strain (when meaning the benchmark unit)

**Read model**:
The basecalling model tier and version of a read set, here always simplex hac or sup at v4.3.0.
_Avoid_: basecaller, accuracy, model (unqualified; clashes with calling models)

**Read set**:
The reads for one Sample × Read model × Depth combination.

**Depth**:
The per-position depth cap a Read set was subsampled to with `rasusa aln` (5, 10, 25 or 50x), applied evenly across every contig. A Read set's actual depth is its measured per-contig mean after realignment, which sits close to the cap on chromosomes and somewhat above it on plasmids.
_Avoid_: coverage (except when quoting a tool's flag); target depth for genome-wide rasusa read subsampling (not used here)

**Mutated reference**:
A sample's own assembly with a donor genome's variants applied. Reads are aligned to it, and calls are scored against its truth set.
_Avoid_: mutref (in prose), donor reference

**Truth set**:
The known variants separating a sample's Mutated reference from its reads, taken from Zenodo 10867171.
_Avoid_: truth VCF, gold standard

### Comparison

**Arm**:
One aligner + caller + calling model combination (A–D) that is run over every Read set.
_Avoid_: condition, method, pipeline

**Arm A (paper)**:
minimap2 2.26 `map-ont` with Clair3 1.0.5, reproducing the eLife configuration.

**Arm B (lr:hq post)**:
minimap2 2.31 `lr:hq` with Clair3 1.0.5, reproducing the minimap2 preset blog post.

**Arm C (current Clair3)**:
minimap2 2.31 `lr:hq` with Clair3 2.0.3 and the HKU-converted v4.3.0 models.

**Arm D (Dorado)**:
minimap2 2.31 `lr:hq` with `dorado polish --bacteria --vcf` 2.1.2.

**Calling model**:
The trained network a variant caller uses (Clair3's `r1041_e82_400bps_{hac,sup}_v430`, or Dorado's bacterial polishing model), as distinct from the Read model.
_Avoid_: model (unqualified)

**Fine-tuned model**:
Clair3's `r1041_e82_400bps_sup_v430_bacteria_finetuned`, trained on 12 of the 14 samples and therefore excluded.

**dnd samples**:
The samples with phosphorothioate (dnd) systems, S. enterica and V. parahaemolyticus, where Dorado's bacterial model is reported to make systematic errors (dorado#1599).

### Evaluation

**Filter chain**:
The NanoVarBench post-call normalisation applied to every Arm's calls before scoring.
_Avoid_: filtering (unqualified; clashes with QUAL/PASS filtering)

**QUAL sweep**:
Scoring the calls at every QUAL threshold, with all records kept regardless of FILTER.

**Best F1**:
The F1 at the QUAL threshold that maximises it for a given variant type, i.e. vcfdist's `THRESHOLD == BEST` row.
_Avoid_: F1 (unqualified, when the threshold matters)

**Default-PASS score**:
Arm D scored using only records Dorado marks PASS, i.e. what a user gets without tuning a threshold.
