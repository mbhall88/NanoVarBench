# Scoring with vcfdist v2.6.4 over the whole genome (ADR-0003). Never pass -s: in v2.6+ it
# is --max-supercluster-size, not --smallest-variant.
#   sweep: every record regardless of FILTER, swept over QUAL (Best F1 is the BEST row)
#   pass:  PASS records only, unthresholded (the Default-PASS score is the NONE row)


rule vcfdist:
    input:
        query=rules.filter_calls.output.vcf,
        truth=rules.extract_truth.output.truth,
        truth_index=rules.index_truth.output.index,
        mutref=rules.extract_truth.output.mutref,
        faidx=rules.index_mutref.output.faidx,
        bed=rules.index_mutref.output.bed,
    output:
        summary=SCORE / "precision-recall-summary.tsv",
        pr=SCORE / "precision-recall.tsv",
        truth=SCORE / "truth.tsv",
        query=SCORE / "query.tsv",
    log:
        LOGS / "vcfdist/{sample}.{read_model}.{depth}x.{arm}.{mode}.log",
    benchmark:
        RESULTS / "benchmarks/vcfdist/{sample}.{read_model}.{depth}x.{arm}.{mode}.tsv"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=30,
    params:
        opts=(
            f"--largest-variant {config['vcfdist']['largest_variant']} "
            f"--credit-threshold {config['vcfdist']['credit_threshold']}"
        ),
        filter=lambda wc: "-f PASS" if wc.mode == "pass" else "",
        max_ram_gb=lambda wc, resources: max(1, int(resources.mem_mb / 1000) - 1),
    container:
        CONTAINERS["vcfdist"]
    shell:
        """
        exec 2> {log}
        # -mx takes an integer, so round the highest QUAL in the calls up.
        MAX_QUAL=$(zcat {input.query} | grep -v '^#' | cut -f 6 | sort -gr | sed -n '1p' \
            | awk '{{c = int($1); if (c < $1) c++; print c}}')
        echo "MAX_QUAL=$MAX_QUAL" >&2
        vcfdist {input.query} {input.truth} {input.mutref} {params.opts} \
            -mx "$MAX_QUAL" -b {input.bed} {params.filter} \
            -t {threads} -r {params.max_ram_gb} -p $(dirname {output.summary})/
        """
