#!/usr/bin/env bash
# Build the Seam 1 fixture: a 40 kb window of ATCC_25922__202309's Mutated reference,
# the matching slice of its Truth set (positions shifted to the window), and the hac reads
# whose primary alignment lies wholly inside the window. Kept for provenance; the outputs
# are committed, so the test never runs this.
#
#   make_fixture.sh <truth dir with ATCC_25922__202309/> <Read set BAM aligned to the
#                   Mutated reference> <apptainer image cache dir>
#
# The fixture was built from data/depth_check/ATCC_25922__202309/rasusa_aln/
# c50_s20240102.realn.bam (an ADR-0005 50x Read set realigned with minimap2 2.31 lr:hq).
set -euo pipefail
TRUTH_DIR=$1
BAM=$2
export APPTAINER_CACHEDIR=$3
OUT=$(cd "$(dirname "$0")" && pwd)
SAMPLE=ATCC_25922__202309
NAME=ATCC_25922_fixture
CONTIG=chromosome
START=3948001   # 1-based, inclusive
END=3988000
LEN=$((END - START + 1))
SAMTOOLS=docker://quay.io/biocontainers/samtools@sha256:a130447589651ed09252aa95a5e4f4132942cdb54d835d81a04a9a930d656561
BCFTOOLS=docker://quay.io/biocontainers/bcftools@sha256:a3e0d3007ffe325c409b398f660840a3e7574d076219c6e82fc994ced87d47c3
ax() { local img=$1; shift; apptainer exec -B /scratch "$img" "$@"; }

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

# Reads whose primary alignment sits wholly inside the window, back in their sequenced
# orientation.
ax $SAMTOOLS samtools view -F 0x904 "$BAM" "$CONTIG:$START-$END" \
  | awk -v s=$START -v e=$END '{
      len = 0; c = $6
      while (match(c, /[0-9]+[MDN=X]/)) { len += substr(c, RSTART, RLENGTH - 1); c = substr(c, RSTART + RLENGTH) }
      if ($4 >= s && $4 + len - 1 <= e) print $1 }' \
  | sort -u > "$tmp/ids.txt"
ax $SAMTOOLS samtools view -b -N "$tmp/ids.txt" -F 0x900 -o "$tmp/reads.bam" "$BAM"
ax $SAMTOOLS samtools fastq "$tmp/reads.bam" 2>/dev/null | gzip -9n > "$OUT/reads.fastq.gz"

# Truth variants the test expects Arm D to call as TPs: those at least 8 kb from either end
# of the window, where the contained reads give full depth.
ax $BCFTOOLS bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$tmp/$NAME/truth.vcf.gz" \
  | awk -v len=$LEN 'BEGIN{OFS="\t"; print "CONTIG", "POS", "REF", "ALT"} $2 > 8000 && $2 < len - 8000' \
  > "$OUT/expected_tp.tsv"

( cd "$OUT" && md5sum "$NAME.tar.gz" reads.fastq.gz > md5.txt )
echo "reads: $(wc -l < "$tmp/ids.txt")  expected TPs: $(( $(wc -l < "$OUT/expected_tp.tsv") - 1 ))"
cat "$OUT/md5.txt"
