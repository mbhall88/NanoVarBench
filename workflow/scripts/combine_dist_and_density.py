import sys

sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path

import cyvcf2


def extract_decision(record: cyvcf2.Variant) -> str:
    """Extract the decision from the FORMAT/BD fields.
    It can take the values TP, FP, FN, or .
    """
    bd = record.format("BD")
    # bd is a numpy array, convert to list
    bd = bd.tolist()
    truth_decision = bd[0]
    query_decision = bd[1]

    if truth_decision == ".":
        if query_decision == ".":
            print(record, file=sys.stderr)
            raise ValueError("Both truth and query decisions are missing")
        return query_decision
    else:
        return truth_decision


def main():
    vcfs = snakemake.input.vcfs
    vcfs.extend(snakemake.input.illumina_vcfs)

    with open(snakemake.output.csv, "w") as f:
        print(
            "sample,caller,mode,model,chrom,pos,ref,alt,decision,dist,density,homlen,vartype",
            file=f,
        )
        for path in map(Path, vcfs):
            vcf = cyvcf2.VCF(str(path))
            sample = path.parent.name
            if "illumina" in str(path):
                caller = "illumina"
                mode = ""
                model = ""
            else:
                caller = path.parts[-7]
                mode = path.parts[-5]
                model = path.parts[-3].split("_")[-1].split("@")[0]

            for record in vcf:
                decision = extract_decision(record)
                dist = record.INFO.get("DIST")
                if dist is None:
                    print(record, file=sys.stderr)
                    print(path, file=sys.stderr)
                    raise ValueError("No DIST field found in the record")
                density = record.INFO.get("DEN")
                if density is None:
                    print(record, file=sys.stderr)
                    print(path, file=sys.stderr)
                    raise ValueError("No DEN field found in the record")

                hom_len = record.INFO.get("HPL")
                if hom_len is None:
                    print(record, file=sys.stderr)
                    print(path, file=sys.stderr)
                    raise ValueError("No HPL field found in the record")

                if record.is_indel:
                    vtype = "INDEL"
                elif record.is_snp:
                    vtype = "SNP"
                else:
                    raise ValueError(f"Unknown variant type {record}")

                print(
                    ",".join(
                        map(
                            str,
                            [
                                sample,
                                caller,
                                mode,
                                model,
                                record.CHROM,
                                record.POS,
                                record.REF,
                                record.ALT[0],
                                decision,
                                dist,
                                density,
                                hom_len,
                                vtype,
                            ],
                        )
                    ),
                    file=f,
                )


main()
