"""This script takes a VCF file and calculates the variant density in a given window 
size around each variant (i.e. if -w 100 then 50bp is taken up and downstream). That 
density value is then added as INFO/DEN tag and written to stdout"""

import argparse
import sys

import cyvcf2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf")
    parser.add_argument("-w", "--window-size", type=int, default=1000)
    args = parser.parse_args()

    vcffile = args.vcf
    window_size = args.window_size

    vcf = cyvcf2.VCF(vcffile)
    vcf2 = cyvcf2.VCF(vcffile)
    vcf.add_info_to_header(
        {
            "ID": "DEN",
            "Description": f"Number of variants in a {window_size}bp window around this variant",
            "Type": "Integer",
            "Number": "1",
        }
    )

    writer = cyvcf2.Writer("-", vcf)

    for variant in vcf:
        chrom = variant.CHROM
        pos = variant.POS

        start = pos - (window_size // 2)
        end = pos + (window_size // 2)

        seqix = vcf.seqnames.index(chrom)
        seqlen = vcf.seqlens[seqix]

        # Account for if start and end are within window size of the start or end of the chromosome
        if start < 0:
            diff = abs(start)
            # add the different to end
            end += diff
            start = 0
        if end > seqlen:
            diff = end - seqlen
            # subtract the difference from start
            start -= diff
            end = seqlen

            if start < 0:
                print(
                    "WARNING: {chrom} is shorter than window size, using entire chromosome",
                    file=sys.stderr,
                )
                start = 0

        den = 0
        region = f"{chrom}:{start}-{end}"
        for v in vcf2(region):
            if v.POS != pos:
                den += 1

        variant.INFO["DEN"] = den
        writer.write_record(variant)

    writer.close()
    vcf.close()
    vcf2.close()


if __name__ == "__main__":
    main()
