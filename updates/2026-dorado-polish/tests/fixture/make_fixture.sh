# Build the Seam 1 fixture: a 40 kb window of ATCC_25922__202309's Mutated reference,
# the matching slice of its Truth set (positions shifted to the window), and, for each Read
# model (hac and sup), reads for that window. Kept for provenance; the outputs are committed,
# so the test never runs this.
#
#   make_fixture.sh <truth dir with ATCC_25922__202309/> <hac BAM> <sup BAM> <apptainer image
#                   cache dir>
#
# The BAMs hold QC'd (seqkit, >=1000 bp, Q>=10) reads of each Read model, aligned with
# minimap2 2.31 lr:hq to the whole Mutated reference, at well over 100x in the window (a random
# subsample of the QC'd reads is enough). Set FIXTURE_OUT to write somewhere other than this
# directory.
#
# Reads: every primary alignment overlapping the window is clipped to the window
# (clip_reads.py), so depth is flat to the window's ends. Reads lying wholly inside the window
# would leave depth ramping up from both ends, so a 50x cap could never be reached. The
# clipped reads are QC'd again (>=1000 bp, Q>=10), aligned to the window, and thinned with
# `rasusa aln -c 52 -s 20240102` to an even ~52x, so the fixture holds little more than the
# 50x Seam 1 asks of it (~1.8 MB each) and stays small enough to commit. After the workflow's
# own cap and realignment, the window's chromosome depth comes out at ~98% of the Depth at
# 5, 10, 25 and 50x.
set -euo pipefail
TRUTH_DIR=$1
HAC_BAM=$2
SUP_BAM=$3
export APPTAINER_CACHEDIR=$4
OUT=${FIXTURE_OUT:-$(cd "$(dirname "$0")" && pwd)}
HERE=$(cd "$(dirname "$0")" && pwd)
SAMPLE=ATCC_25922__202309
NAME=ATCC_25922_fixture
CONTIG=chromosome
START=3948001   # 1-based, inclusive
END=3988000
LEN=$((END - START + 1))
POOL_DEPTH=52
SEED=20240102
SAMTOOLS=docker://quay.io/biocontainers/samtools@sha256:a130447589651ed09252aa95a5e4f4132942cdb54d835d81a04a9a930d656561
BCFTOOLS=docker://quay.io/biocontainers/bcftools@sha256:a3e0d3007ffe325c409b398f660840a3e7574d076219c6e82fc994ced87d47c3
SEQKIT=docker://quay.io/biocontainers/seqkit@sha256:45fb535880be37dfed5be5517111fb8bfdd6234ef36e725b126a6131b1af2ef0
RASUSA=docker://quay.io/biocontainers/rasusa@sha256:e8c7b92c66abd96fdb861a4b7977c6ac47f1e9df6428d539a6616bb257bb0846
MINIMAP2=docker://quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4@sha256:966a1318a02cc3cda1785ccf62a4db2390e88dd48f303befa6c4a0a89a241a49
ax() { local img=$1; shift; apptainer exec -B /scratch "$img" "$@"; }

mkdir -p "$OUT"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/$NAME"

# Mutated reference window, renamed to the original contig name.
ax $SAMTOOLS samtools faidx "$TRUTH_DIR/$SAMPLE/mutreference.fna" "$CONTIG:$START-$END" \
  | sed "1s/.*/>$CONTIG/" > "$tmp/$NAME/mutreference.fna"

# Truth set slice, shifted into window coordinates.
{
  ax $BCFTOOLS bcftools view -h "$TRUTH_DIR/$SAMPLE/truth.vcf.gz" | grep -v '^##contig=' | sed '$d'
  echo "##contig=<ID=$CONTIG,length=$LEN>"
  ax $BCFTOOLS bcftools view -h "$TRUTH_DIR/$SAMPLE/truth.vcf.gz" | tail -1
  ax $BCFTOOLS bcftools view -H -r "$CONTIG:$START-$END" "$TRUTH_DIR/$SAMPLE/truth.vcf.gz" \
    | awk -v off=$((START - 1)) 'BEGIN{OFS="\t"} {$2 = $2 - off; print}'
} | ax $BCFTOOLS bcftools view -Oz -o "$tmp/$NAME/truth.vcf.gz"
tar -C "$tmp" -czf "$OUT/$NAME.tar.gz" "$NAME"

# Reads for one Read model, back in their sequenced orientation.
make_reads() {
  local rm=$1 bam=$2 w="$tmp/$1"
  mkdir -p "$w"
  ax $SAMTOOLS samtools view -F 0x904 "$bam" "$CONTIG:$START-$END" \
    | python3 "$HERE/clip_reads.py" $START $END > "$w/clipped.fq"
  ax $SEQKIT seqkit seq --min-len 1000 --min-qual 10 -o "$w/clipped.qc.fq" "$w/clipped.fq"
  ax $MINIMAP2 sh -c "minimap2 -t 4 -a -x lr:hq --secondary=no $tmp/$NAME/mutreference.fna $w/clipped.qc.fq \
      | samtools view -u -F 0x904 - | samtools sort -o $w/primary.bam - && samtools index $w/primary.bam"
  ax $RASUSA rasusa aln -c $POOL_DEPTH -s $SEED -O sam "$w/primary.bam" \
    | grep -v '^@' | cut -f1 | sort -u > "$w/ids.txt"
  ax $SEQKIT seqkit grep -f "$w/ids.txt" "$w/clipped.qc.fq" -o "$w/pool.fq"
  gzip -9n -c "$w/pool.fq" > "$OUT/reads.$rm.fastq.gz"
  echo "$rm: $(wc -l < "$w/ids.txt") reads"
}
make_reads hac "$HAC_BAM"
make_reads sup "$SUP_BAM"

# Truth variants the test expects every Arm to call as TPs (bar the Clair3 FNs listed in
# clair3_known_fn.<Read model>.tsv): those at least 8 kb from either end of the window.
ax $BCFTOOLS bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$tmp/$NAME/truth.vcf.gz" \
  | awk -v len=$LEN 'BEGIN{OFS="\t"; print "CONTIG", "POS", "REF", "ALT"} $2 > 8000 && $2 < len - 8000' \
  > "$OUT/expected_tp.tsv"

( cd "$OUT" && md5sum "$NAME.tar.gz" reads.hac.fastq.gz reads.sup.fastq.gz > md5.txt )
echo "expected TPs: $(( $(wc -l < "$OUT/expected_tp.tsv") - 1 ))"
cat "$OUT/md5.txt"
