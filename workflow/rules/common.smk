from pathlib import Path


def infer_original_reads_path(wildcards):
    reads_dir = Path(pep.get_sample(wildcards.sample)["reads_dir"])
    return (
        reads_dir
        / f"{wildcards.mode}/{wildcards.version}/{wildcards.model}/{wildcards.sample}.fq.gz"
    )


def infer_reference_genome(wildcards):
    return Path(pep.get_sample(wildcards.sample)["reference_path"])


def infer_species(sample):
    return pep.get_sample(sample)["species"]


def infer_taxid(wildcards):
    return pep.get_sample(wildcards.sample)["taxid"]


def infer_vcf_reference(wildcards):
    if wildcards.ref == "self":
        return str(RESULTS / f"reference/{wildcards.sample}.fa")
    elif wildcards.ref == "mutref":
        return str(RESULTS / f"truth/{wildcards.sample}/mutreference.fna")


def infer_vcf_reference_faidx(wildcards):
    return infer_vcf_reference(wildcards) + ".fai"
