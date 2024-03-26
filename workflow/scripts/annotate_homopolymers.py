"""This script takes a VCF and a reference fasta file and annotates homopolymers in the 
VCF with the length of the homopolymer. The homopolymer length is added as INFO/HPL tag
"""

import argparse
import sys
from typing import Dict

import cyvcf2


def fasta_to_dict(fasta: str) -> Dict[str, str]:
    """Convert a fasta file to a dictionary"""
    d = {}
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                name = line.split()[0][1:]
                d[name] = ""
            else:
                d[name] += line.strip()
    return d


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf")
    parser.add_argument("reference")
    args = parser.parse_args()

    seqs = fasta_to_dict(args.reference)
    vcf = cyvcf2.VCF(args.vcf)

    # add HPL to header
    vcf.add_info_to_header(
        {
            "ID": "HPL",
            "Description": "Length of the homopolymer the variant is within",
            "Type": "Integer",
            "Number": "1",
        }
    )

    writer = cyvcf2.Writer("-", vcf)
    writer.write_header()

    for variant in vcf:
        chrom = variant.CHROM
        pos = variant.POS
        ref = variant.REF
        alt = variant.ALT[0]
        start = pos - 1
        if pos == 6111:
            print(variant, file=sys.stderr)
            print(ref, alt, file=sys.stderr)
            print(seqs[chrom][start : start + 10], file=sys.stderr)

        if variant.is_deletion:
            # check that the first and last base of the REF are different
            # homopolymer deletions will often have a REF allele such as CAAA>C or GC>G or CAAG>C
            # we want to start counting length from the start of the run of the last base in REF that is the same
            # but different from the first base in REF
            last_base = ref[-1]
            if pos == 6111:
                print(last_base, file=sys.stderr)
            ix_in_ref = len(ref)
            if pos == 6111:
                print(ix_in_ref, file=sys.stderr)
            for base in ref[1:][::-1]:
                if pos == 6111:
                    print((base, last_base, ix_in_ref), file=sys.stderr)
                if base != last_base:
                    break
                ix_in_ref -= 1

            if pos == 6111:
                print(ix_in_ref, file=sys.stderr)
            start = start + ix_in_ref

            if pos == 6111:
                print(start, file=sys.stderr)

            hom_start_base = last_base
            hom_len = 0
            if pos == 6111:
                print(seqs[chrom][start : start + 10], file=sys.stderr)

            for base in seqs[chrom][start:]:
                if pos == 6111:
                    print((base, hom_start_base, hom_len), file=sys.stderr)
                if base != hom_start_base:
                    break
                hom_len += 1

            if pos == 6111:
                print(hom_len, file=sys.stderr)

            variant.INFO["HPL"] = hom_len

        elif variant.is_indel:  # must be an insertion
            if len(ref) > len(alt):
                print(variant, file=sys.stderr)
                raise ValueError("Expected insertion, got deletion")
            # do the same, but for the ALT allele
            last_base = alt[-1]
            ix_in_alt = len(alt) - 1
            hom_len = 0
            for base in alt[1:][-1]:
                if base != last_base:
                    break
                ix_in_alt -= 1
                hom_len += 1

            start += 1
            hom_start_base = last_base
            for base in seqs[chrom][start:]:
                if base != hom_start_base:
                    break
                hom_len += 1

            variant.INFO["HPL"] = hom_len

        else:
            variant.INFO["HPL"] = 0

        writer.write_record(variant)

    writer.close()
    vcf.close()


if __name__ == "__main__":
    main()
