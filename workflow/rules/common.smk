from pathlib import Path


def infer_original_reads_path(wildcards):
    reads_dir = Path(pep.get_sample(wildcards.sample)["reads_dir"])
    return reads_dir / f"{wildcards.version}/{wildcards.model}/{wildcards.sample}.fq.gz"


def infer_reference_genome(wildcards):
    return Path(pep.get_sample(wildcards.sample)["reference_path"])
