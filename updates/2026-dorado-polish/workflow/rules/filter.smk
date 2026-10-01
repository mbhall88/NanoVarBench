# The NanoVarBench Filter chain, copied from the paper's `filter_variants` rule
# (workflow/rules/call.smk at the repo root) and applied unchanged to every Arm's calls.


# Each caller's raw calls for an Arm (bcftools reads both plain and bgzipped VCFs).
RAW_CALLS = {"dorado": CALL / "variants.vcf", "clair3": CALL / "variants.vcf.gz"}


def raw_calls(wildcards):
    caller = ARMS[wildcards.arm]["caller"]
    try:
        return RAW_CALLS[caller]
    except KeyError:
        raise ValueError(f"No calling rule for caller {caller} (Arm {wildcards.arm})")


rule filter_calls:
    input:
        vcf=raw_calls,
        reference=rules.extract_truth.output.mutref,
        faidx=rules.index_mutref.output.faidx,
        filter_script=SCRIPTS / "filter_hets.py",
    output:
        vcf=RESULTS / "calls/{sample}/{read_model}/{depth}x/{arm}.filter.vcf.gz",
        csi=RESULTS / "calls/{sample}/{read_model}/{depth}x/{arm}.filter.vcf.gz.csi",
    log:
        LOGS / "filter_calls/{sample}.{read_model}.{depth}x.{arm}.log",
    resources:
        mem_mb=1000,
        runtime=10,
    conda:
        ENVS / "filter.yaml"
    params:
        max_indel=config["max_indel"],
    shell:
        """
        exec 2> {log}
        contigs=$(mktemp -u).contigs.txt
        header=$(mktemp -u).header.txt
        trap 'rm -f $contigs $header' EXIT

        # bcftools reheader only adds contigs that appear in the VCF, we want all contigs
        awk '{{print "##contig=<ID="$1",length="$2">"}}' {input.faidx} > "$contigs"  # make contig lines with all contigs
        (bcftools view -h {input.vcf} |                                            # output VCF header
            grep -v "^##contig=" |                                                 # remove contig lines
            sed -e "3r $contigs") > "$header"                                      # add contig lines after 3rd line

        (bcftools reheader -h "$header" {input.vcf} |                       # replace VCF header with new header containing all contigs
            python {input.filter_script} |                                  # make heterozygous calls homozygous for allele with most depth
            bcftools view -i 'GT="alt"' |                                   # remove non-alt alleles
            bcftools view -e 'ALT="."' |                                    # remove sites with no alt allele (NanoCaller bug)
            bcftools norm -f {input.reference} -a -c e -m - |               # normalise and left-align indels
            bcftools norm -aD |                                             # remove duplicates after normalisation
            bcftools filter -e 'abs(ILEN)>{params.max_indel} || ALT="*"' |  # remove long indels or sites with unobserved alleles
            bcftools +setGT - -- -t a -n c:M |                              # make genotypes haploid e.g., 1/1 -> 1
            bcftools sort |                                                 # sort VCF
            bcftools view -i 'GT="A"' -o {output.vcf})                      # remove non-alt alleles and write index
        bcftools index -f {output.vcf}
        """
