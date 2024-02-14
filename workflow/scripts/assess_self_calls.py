import sys

sys.stderr = open(snakemake.log[0], "w")
import json
from collections import Counter, defaultdict
from pathlib import Path
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


def canonical_kmer(kmer: str) -> str:
    """Return the canonical kmer"""
    return min(kmer, kmer[::-1])


def main():
    seqs = fasta_to_dict(snakemake.input.reference)
    vcf = cyvcf2.VCF(snakemake.input.vcf)
    seqlens = sum(vcf.seqlens)
    variants = {"SNP": 0, "INS": 0, "DEL": 0}
    substitution_types = defaultdict(int)
    sub_contexts = defaultdict(int)
    hom_indels = defaultdict(list)

    for variant in vcf:
        if variant.is_snp:
            variants["SNP"] += 1
            substitution_type = f"{variant.REF}>{variant.ALT[0]}"
            substitution_types[substitution_type] += 1
            # extract 5-mer with this base in the middle
            start = variant.start - 2
            end = variant.start + 3
            kmer = seqs[variant.CHROM][start:end]
            assert kmer[2] == variant.REF
            kmer = canonical_kmer(kmer)
            sub_contexts[kmer] += 1

        elif variant.is_indel:
            alt = variant.ALT[0]
            if len(variant.REF) > len(alt):
                variants["DEL"] += 1
                # check if homopolymer
                if len(alt) == 1 and variant.REF[0] == alt:
                    deletion = variant.REF[1:]
                    unique_bases = set(deletion)
                    if len(unique_bases) == 1:
                        hom_indels["hom_del"].append(len(deletion))
            else:
                variants["INS"] += 1
                # check if homopolymer
                if len(variant.REF) == 1 and variant.REF == alt[0]:
                    insertion = alt[1:]
                    unique_bases = set(insertion)
                    if len(unique_bases) == 1:
                        hom_indels["hom_ins"].append(len(insertion))
        else:
            raise ValueError(f"Unknown variant type {variant}")

    n_fps = sum(variants.values())
    n_tns = seqlens - n_fps
    fpr = n_fps / seqlens
    snp_tns = seqlens - variants["SNP"]
    indel_tns = seqlens - (variants["INS"] + variants["DEL"])
    snp_fpr = variants["SNP"] / seqlens
    indel_fpr = (variants["INS"] + variants["DEL"]) / seqlens
    with open(snakemake.output.csv, "w") as output, open(
        snakemake.output.json, "w"
    ) as output_json:
        print(
            ",".join(
                [
                    "SNP",
                    "INS",
                    "DEL",
                    "INDEL",
                    "FP",
                    "SNP_TN",
                    "INDEL_TN",
                    "TN",
                    "SNP_FPR",
                    "INDEL_FPR",
                    "FPR",
                    "caller",
                    "sample",
                    "depth",
                    "mode",
                    "version",
                    "model",
                    "2bp_hom_del",
                    "3bp_hom_del",
                    "4bp_hom_del",
                    "5+bp_hom_del",
                    "hom_del",
                    "2bp_hom_ins",
                    "3bp_hom_ins",
                    "4bp_hom_ins",
                    "5+bp_hom_ins",
                    "hom_ins",
                    "hom_indels",
                ]
            ),
            file=output,
        )
        p = Path(snakemake.input.vcf)
        caller = p.parts[-6]
        depth = p.parts[-5]
        mode = p.parts[-4]
        version = p.parts[-3]
        model = p.parts[-2]
        sample = p.name.split(".")[0]

        hom_dels = Counter(hom_indels["hom_del"])
        hom_ins = Counter(hom_indels["hom_ins"])

        prop_2bp_hom_del = (
            hom_dels.get(2, 0) / variants["DEL"] if variants["DEL"] != 0 else 0
        )
        prop_3bp_hom_del = (
            hom_dels.get(3, 0) / variants["DEL"] if variants["DEL"] != 0 else 0
        )
        prop_4bp_hom_del = (
            hom_dels.get(4, 0) / variants["DEL"] if variants["DEL"] != 0 else 0
        )
        prop_5bp_or_more_hom_del = (
            (sum([v for k, v in hom_dels.items() if k >= 5]) / variants["DEL"])
            if variants["DEL"] != 0
            else 0
        )
        prop_hom_del = (
            (sum([v for k, v in hom_dels.items()]) / variants["DEL"])
            if variants["DEL"] != 0
            else 0
        )
        prop_2bp_hom_ins = (
            hom_ins.get(2, 0) / variants["INS"] if variants["INS"] != 0 else 0
        )
        prop_3bp_hom_ins = (
            hom_ins.get(3, 0) / variants["INS"] if variants["INS"] != 0 else 0
        )
        prop_4bp_hom_ins = (
            hom_ins.get(4, 0) / variants["INS"] if variants["INS"] != 0 else 0
        )
        prop_5bp_or_more_hom_ins = (
            (sum([v for k, v in hom_ins.items() if k >= 5]) / variants["INS"])
            if variants["INS"] != 0
            else 0
        )
        prop_hom_ins = (
            (sum([v for k, v in hom_ins.items()]) / variants["INS"])
            if variants["INS"] != 0
            else 0
        )
        prop_hom_indels = (
            (
                sum([v for k, v in hom_dels.items()])
                + sum([v for k, v in hom_ins.items()])
            )
            / (variants["INS"] + variants["DEL"])
            if (variants["INS"] + variants["DEL"]) != 0
            else 0
        )

        print(
            ",".join(
                [
                    str(variants["SNP"]),
                    str(variants["INS"]),
                    str(variants["DEL"]),
                    str(variants["INS"] + variants["DEL"]),
                    str(n_fps),
                    str(snp_tns),
                    str(indel_tns),
                    str(n_tns),
                    str(snp_fpr),
                    str(indel_fpr),
                    str(fpr),
                    caller,
                    sample,
                    depth,
                    mode,
                    version,
                    model,
                    str(prop_2bp_hom_del),
                    str(prop_3bp_hom_del),
                    str(prop_4bp_hom_del),
                    str(prop_5bp_or_more_hom_del),
                    str(prop_hom_del),
                    str(prop_2bp_hom_ins),
                    str(prop_3bp_hom_ins),
                    str(prop_4bp_hom_ins),
                    str(prop_5bp_or_more_hom_ins),
                    str(prop_hom_ins),
                    str(prop_hom_indels),
                ]
            ),
            file=output,
        )

        json.dump(
            {
                "SNP": variants["SNP"],
                "INS": variants["INS"],
                "DEL": variants["DEL"],
                "INDEL": variants["INS"] + variants["DEL"],
                "FP": n_fps,
                "SNP_TN": snp_tns,
                "INDEL_TN": indel_tns,
                "TN": n_tns,
                "SNP_FPR": snp_fpr,
                "INDEL_FPR": indel_fpr,
                "FPR": fpr,
                "caller": caller,
                "sample": sample,
                "depth": depth,
                "mode": mode,
                "version": version,
                "model": model,
                "substition_types": substitution_types,
                "homopolymer_deletions": hom_dels,
                "homopolymer_insertions": hom_ins,
                "2bp_hom_del": prop_2bp_hom_del,
                "3bp_hom_del": prop_3bp_hom_del,
                "4bp_hom_del": prop_4bp_hom_del,
                "5+bp_hom_del": prop_5bp_or_more_hom_del,
                "hom_del": prop_hom_del,
                "2bp_hom_ins": prop_2bp_hom_ins,
                "3bp_hom_ins": prop_3bp_hom_ins,
                "4bp_hom_ins": prop_4bp_hom_ins,
                "5+bp_hom_ins": prop_5bp_or_more_hom_ins,
                "hom_ins": prop_hom_ins,
                "hom_indels": prop_hom_indels,
            },
            output_json,
            indent=4,
            sort_keys=True,
        )


if __name__ == "__main__":
    main()
