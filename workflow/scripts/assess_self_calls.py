import sys

sys.stderr = open(snakemake.log[0], "w")
from pathlib import Path

import cyvcf2


def main():
    vcf = cyvcf2.VCF(snakemake.input.vcf)
    seqlens = sum(vcf.seqlens)
    variants = {"SNP": 0, "INS": 0, "DEL": 0}

    for variant in vcf:
        if variant.is_snp:
            variants["SNP"] += 1
        elif variant.is_indel:
            if len(variant.REF) > len(variant.ALT[0]):
                variants["DEL"] += 1
            else:
                variants["INS"] += 1
        else:
            raise ValueError(f"Unknown variant type {variant}")

    n_fps = sum(variants.values())
    n_tns = seqlens - n_fps
    fpr = n_fps / seqlens
    snp_tns = seqlens - variants["SNP"]
    indel_tns = seqlens - (variants["INS"] + variants["DEL"])
    snp_fpr = variants["SNP"] / seqlens
    indel_fpr = (variants["INS"] + variants["DEL"]) / seqlens
    with open(snakemake.output.csv, "w") as output:
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
                ]
            ),
            file=output,
        )


if __name__ == "__main__":
    main()
