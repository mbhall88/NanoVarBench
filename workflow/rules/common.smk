from pathlib import Path


def infer_original_reads_path(wildcards):
    reads_dir = Path(pep.get_sample(wildcards.sample)["reads_dir"])
    return (
        reads_dir
        / f"{wildcards.mode}/{wildcards.version}/{wildcards.model}/{wildcards.sample}.fq.gz"
    )


def infer_reference_genome(wildcards):
    return Path(pep.get_sample(wildcards.sample)["reference_path"])


def infer_illumina_reads_1(wildcards):
    return Path(pep.get_sample(wildcards.sample)["illumina_1"])


def infer_illumina_reads_2(wildcards):
    return Path(pep.get_sample(wildcards.sample)["illumina_2"])


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


def infer_illumina_vcf_to_filter(wildcards):
    if wildcards.ref == "self":
        return (
            RESULTS / f"call/self/illumina/{wildcards.sample}/{wildcards.sample}.vcf.gz"
        )
    elif wildcards.ref == "mutref":
        return (
            RESULTS
            / f"call/mutref/illumina/{wildcards.sample}/{wildcards.sample}.raw.vcf.gz"
        )
    else:
        raise ValueError(f"Unknown ref: {wildcards.ref}")


def get_min_completeness(wildcards):
    default = truth_config["min_completeness"]
    sample_specific_truth_config = config.get("sample_specific_truth", {}).get(
        wildcards.sample
    )

    if sample_specific_truth_config is not None:
        return sample_specific_truth_config.get("min_completeness", default)

    return default


def get_min_completeness_percentile(wildcards):
    default = truth_config["min_completeness_percentile"]
    sample_specific_truth_config = config.get("sample_specific_truth", {}).get(
        wildcards.sample
    )

    if sample_specific_truth_config is not None:
        return sample_specific_truth_config.get("min_completeness_percentile", default)

    return default


def get_max_contamination(wildcards):
    default = truth_config["max_contamination"]
    sample_specific_truth_config = config.get("sample_specific_truth", {}).get(
        wildcards.sample
    )

    if sample_specific_truth_config is not None:
        return sample_specific_truth_config.get("max_contamination", default)

    return default
