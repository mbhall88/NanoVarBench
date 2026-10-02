"""Resolve Clair3's diploid HET calls by allele frequency (the AF filter, #20).

A HET (0/1, or 1/2 at a multi-allelic site) is made homozygous for the ALT in its genotype
with the highest FORMAT/AF when that AF is >= --min-af, and homozygous REF otherwise, so the
Filter chain's next steps drop it. Homozygous calls are left alone: they are what
--haploid_precise keeps. A HET without an AF is made homozygous REF and reported on stderr,
like filter_hets.py does for a HET without AD or AC.

    python af_filter.py --min-af 0.6 in.vcf.gz -o out.vcf.gz
"""

import argparse
import math
import sys

import cyvcf2


def resolve(alleles, afs, min_af):
    """The homozygous allele a HET with these (non-missing) alleles becomes: the ALT with
    the highest AF if that AF is >= min_af, else REF (0). Ties go to the first ALT listed."""
    best, best_af = 0, None
    for a in sorted(set(alleles)):
        if a == 0:
            continue
        af = afs[a - 1] if a - 1 < len(afs) else None
        if af is None or math.isnan(af) or af < 0:  # cyvcf2 writes missing values as nan or < 0
            continue
        if best_af is None or af > best_af:
            best, best_af = a, af
    if best_af is None:
        return None
    return best if best_af >= min_af else 0


def main(args):
    vcf = cyvcf2.VCF(args.vcf)
    vcf.add_to_header(
        f"##af_filter=<min_af={args.min_af},Description=\"HETs made homozygous for the ALT with "
        f"the highest FORMAT/AF when it is >= {args.min_af}, else homozygous REF\">"
    )
    out = cyvcf2.Writer(args.output, vcf, mode="wz")
    counts = {"het_to_alt": 0, "het_to_ref": 0, "het_no_af": 0, "not_het": 0}
    for variant in vcf:
        gt = variant.genotypes[0]
        alleles = [a for a in gt[:-1] if a >= 0]  # the last item is the phased flag
        if len(set(alleles)) < 2:
            counts["not_het"] += 1
            out.write_record(variant)
            continue
        afs = []
        if "AF" in variant.FORMAT:
            afs = [float(x) for x in variant.format("AF")[0]]
        allele = resolve(alleles, afs, args.min_af)
        if allele is None:
            print(
                f"No AF for HET at {variant.CHROM}:{variant.POS}; making it homozygous REF",
                file=sys.stderr,
            )
            counts["het_no_af"] += 1
            allele = 0
        else:
            counts["het_to_alt" if allele else "het_to_ref"] += 1
        variant.genotypes = [[allele, allele, False]]
        out.write_record(variant)
    out.close()
    print(" ".join(f"{k}={v}" for k, v in counts.items()), file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf", help="Clair3's VCF, called without a haploid mode")
    parser.add_argument("--min-af", type=float, required=True, help="the AF threshold")
    parser.add_argument("-o", "--output", required=True, help="bgzipped VCF to write")
    main(parser.parse_args())
