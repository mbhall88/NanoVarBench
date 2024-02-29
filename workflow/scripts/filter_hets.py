"""This script takes a VCF file and forces HETs to be homozygous for the allele with 
the highest depth
"""

import argparse

import cyvcf2


def main(args):
    vcf = cyvcf2.VCF(args.vcf, gts012=True)
    vcf_out = cyvcf2.Writer("-", vcf, mode="w")
    for variant in vcf:
        v_type = variant.gt_types[0]
        if v_type == 1:  # HET
            use_ref = True
            allelic_counts = variant.INFO.get("AC")
            if "AD" in variant.FORMAT:
                allele_depths = variant.format("AD")[0]
                ref_depth = allele_depths[0]
                alt_depth = allele_depths[1]
                if alt_depth > ref_depth:
                    use_ref = False
            elif allelic_counts is not None:
                ref_count = allelic_counts[0]
                alt_count = allelic_counts[1]
                if alt_count > ref_count:
                    use_ref = False
            else:
                raise KeyError(f"Could not find allele counts for variant {variant}")

            if not use_ref:
                variant.genotypes[0] = [1, 1, False]
            else:
                variant.genotypes[0] = [0, 0, False]
            variant.genotypes = variant.genotypes

        vcf_out.write_record(variant)

    vcf_out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "vcf",
        help="Input VCF file",
        default="-",
        type=argparse.FileType("r"),
        nargs="?",
    )
    args = parser.parse_args()
    main(args)
