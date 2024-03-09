rule assess_self_calls:
    input:
        vcf=RESULTS
        / "call/self/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.filter.vcf.gz",
        reference=rules.align_to_self.input.reference,
    output:
        csv=RESULTS
        / "assess/self/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.csv",
        json=RESULTS
        / "assess/self/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.json",
    log:
        LOGS
        / "assess_self_calls/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.log",
    resources:
        runtime="2m",
        mem_mb=500,
    conda:
        ENVS / "assess_self_calls.yaml"
    script:
        SCRIPTS / "assess_self_calls.py"


use rule assess_self_calls as assess_self_calls_illumina with:
    input:
        vcf=RESULTS / "call/self/illumina/{sample}/{sample}.filter.vcf.gz",
        reference=rules.call_self_illumina.input.reference,
    output:
        csv=RESULTS / "assess/self/illumina/{sample}.csv",
        json=RESULTS / "assess/self/illumina/{sample}.json",
    log:
        LOGS / "assess_self_calls_illumina/{sample}.log",


rule combine_self_calls_csvs:
    input:
        csvs=expand(
            RESULTS
            / "assess/self/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.csv",
            caller=CALLERS,
            depth=[MAX_DEPTH],
            mode=MODES,
            version=VERSIONS,
            model=MODELS,
            sample=SAMPLES,
        ),
    output:
        csv=RESULTS / "assess/self/self_calls.csv",
    log:
        LOGS / "combine_self_calls_csvs.log",
    resources:
        runtime="2m",
        mem_mb=500,
    container:
        "docker://python:3.11-slim"
    script:
        SCRIPTS / "combine_csvs.py"


use rule combine_self_calls_csvs as combine_self_calls_csvs_illumina with:
    input:
        csvs=expand(
            RESULTS / "assess/self/illumina/{sample}.csv",
            sample=SAMPLES,
        ),
    output:
        csv=RESULTS / "assess/self/self_calls_illumina.csv",
    log:
        LOGS / "combine_self_calls_csvs_illumina.log",


rule combine_self_calls_jsons:
    input:
        jsons=expand(
            RESULTS
            / "assess/self/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.json",
            caller=CALLERS,
            depth=[MAX_DEPTH],
            mode=MODES,
            version=VERSIONS,
            model=MODELS,
            sample=SAMPLES,
        ),
    output:
        json=RESULTS / "assess/self/self_calls.json",
    log:
        LOGS / "combine_self_calls_jsons.log",
    resources:
        runtime="2m",
        mem_mb=2 * GB,
    container:
        "docker://python:3.11-slim"
    script:
        SCRIPTS / "combine_jsons.py"


use rule combine_self_calls_jsons as combine_self_calls_jsons_illumina with:
    input:
        jsons=expand(
            RESULTS / "assess/self/illumina/{sample}.json",
            sample=SAMPLES,
        ),
    output:
        json=RESULTS / "assess/self/self_calls_illumina.json",
    log:
        LOGS / "combine_self_calls_jsons_illumina.log",


rule plot_self_fpr:
    input:
        csv=rules.combine_self_calls_csvs.output.csv,
        illumina_csv=rules.combine_self_calls_csvs_illumina.output.csv,
    output:
        snps_png=FIGURES / "assess/self/self_fpr.snps.png",
        indels_png=FIGURES / "assess/self/self_fpr.indels.png",
        csv=TABLES / "assess/self/self_fpr.csv",
    log:
        LOGS / "plot_self_fpr.log",
    resources:
        runtime="5m",
        mem_mb=500,
    params:
        sample2species={sample: infer_species(sample) for sample in SAMPLES},
        no_indels=config.get("no_indels", []),
    conda:
        ENVS / "plot_self_fpr.yaml"
    script:
        SCRIPTS / "plot_self_fpr.py"


rule assess_mutref_calls:
    input:
        query_vcf=RESULTS
        / "call/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.filter.vcf.gz",
        truth_vcf=rules.create_mutref.output.truth_vcf,
        mutref=rules.create_mutref.output.mutref,
        faidx=rules.faidx_mutref.output.faidx,
        bed=rules.make_full_bed.output.bed,
    output:
        pr_summary=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.precision-recall-summary.tsv",
        pr=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.precision-recall.tsv",
        summary=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.summary.vcf",
        query=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.query.tsv",
        truth=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.truth.tsv",
        dist_summary=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.distance-summary.tsv",
        dist=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.distance.tsv",
    log:
        LOGS
        / "assess_mutref_calls/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.log",
    resources:
        runtime="5m",
        mem_mb=500,
    container:
        "docker://timd1/vcfdist:v2.3.3"
    params:
        opts=f"--largest-variant {config['truth']['max_indel']} --credit-threshold 1.0 -d",
        prefix=lambda wildcards, output: Path(output.pr).with_suffix("").with_suffix(""),
    shell:
        """
        exec 2> {log}
        echo "Calculated maximum QUAL score..." 1>&2
        MAX_QUAL=$(bgzip -dc {input.query_vcf} | grep -v '^#' | cut -f 6 | sort -gr | sed -n '1p')
        echo "MAX_QUAL=$MAX_QUAL" 1>&2
        echo "Running vcfdist..." 1>&2
        vcfdist {input.query_vcf} {input.truth_vcf} {input.mutref} {params.opts} \
            -p {params.prefix}. -b {input.bed} -mx $MAX_QUAL
        """


use rule assess_mutref_calls as assess_mutref_calls_illumina with:
    input:
        query_vcf=RESULTS / "call/mutref/illumina/{sample}/{sample}.filter.vcf.gz",
        truth_vcf=rules.create_mutref.output.truth_vcf,
        mutref=rules.create_mutref.output.mutref,
        faidx=rules.faidx_mutref.output.faidx,
        bed=rules.make_full_bed.output.bed,
    output:
        pr_summary=RESULTS
        / "assess/mutref/illumina/{sample}/{sample}.precision-recall-summary.tsv",
        pr=RESULTS / "assess/mutref/illumina/{sample}/{sample}.precision-recall.tsv",
        summary=RESULTS / "assess/mutref/illumina/{sample}/{sample}.summary.vcf",
        query=RESULTS / "assess/mutref/illumina/{sample}/{sample}.query.tsv",
        truth=RESULTS / "assess/mutref/illumina/{sample}/{sample}.truth.tsv",
        dist_summary=RESULTS
        / "assess/mutref/illumina/{sample}/{sample}.distance-summary.tsv",
        dist=RESULTS / "assess/mutref/illumina/{sample}/{sample}.distance.tsv",
    log:
        LOGS / "assess_mutref_calls_illumina/{sample}.log",


rule assess_mutref_calls_without_repetitive_regions:
    input:
        query_vcf=RESULTS
        / "call/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.filter.vcf.gz",
        truth_vcf=rules.create_mutref.output.truth_vcf,
        mutref=rules.create_mutref.output.mutref,
        faidx=rules.faidx_mutref.output.faidx,
        bed=rules.identify_repetitive_regions.output.unique_bed,
    output:
        pr_summary=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.without_repetitive_regions.precision-recall-summary.tsv",
        pr=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.without_repetitive_regions.precision-recall.tsv",
        summary=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.without_repetitive_regions.summary.vcf",
        query=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.without_repetitive_regions.query.tsv",
        truth=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.without_repetitive_regions.truth.tsv",
        dist_summary=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.without_repetitive_regions.distance-summary.tsv",
        dist=RESULTS
        / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.without_repetitive_regions.distance.tsv",
    log:
        LOGS
        / "assess_mutref_calls_without_repetitive_regions/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.log",
    resources:
        runtime="5m",
        mem_mb=500,
    container:
        "docker://timd1/vcfdist:v2.3.3"
    params:
        opts=f"--largest-variant {config['truth']['max_indel']} --credit-threshold 1.0 -d",
        prefix=lambda wildcards, output: Path(output.pr).with_suffix("").with_suffix(""),
    shell:
        """
        exec 2> {log}
        echo "Calculated maximum QUAL score..." 1>&2
        MAX_QUAL=$(bgzip -dc {input.query_vcf} | grep -v '^#' | cut -f 6 | sort -gr | sed -n '1p')
        echo "MAX_QUAL=$MAX_QUAL" 1>&2
        echo "Running vcfdist..." 1>&2
        vcfdist {input.query_vcf} {input.truth_vcf} {input.mutref} {params.opts} \
            -p {params.prefix}. -b {input.bed} -mx $MAX_QUAL
        """


use rule assess_mutref_calls_without_repetitive_regions as assess_mutref_calls_illumina_without_repetitive_regions with:
    input:
        query_vcf=RESULTS / "call/mutref/illumina/{sample}/{sample}.filter.vcf.gz",
        truth_vcf=rules.create_mutref.output.truth_vcf,
        mutref=rules.create_mutref.output.mutref,
        faidx=rules.faidx_mutref.output.faidx,
        bed=rules.identify_repetitive_regions.output.unique_bed,
    output:
        pr_summary=RESULTS
        / "assess/mutref/illumina/{sample}/{sample}.without_repetitive_regions.precision-recall-summary.tsv",
        pr=RESULTS
        / "assess/mutref/illumina/{sample}/{sample}.without_repetitive_regions.precision-recall.tsv",
        summary=RESULTS
        / "assess/mutref/illumina/{sample}/{sample}.without_repetitive_regions.summary.vcf",
        query=RESULTS
        / "assess/mutref/illumina/{sample}/{sample}.without_repetitive_regions.query.tsv",
        truth=RESULTS
        / "assess/mutref/illumina/{sample}/{sample}.without_repetitive_regions.truth.tsv",
        dist_summary=RESULTS
        / "assess/mutref/illumina/{sample}/{sample}.without_repetitive_regions.distance-summary.tsv",
        dist=RESULTS
        / "assess/mutref/illumina/{sample}/{sample}.without_repetitive_regions.distance.tsv",
    log:
        LOGS / "assess_mutref_calls_illumina_without_repetitive_regions/{sample}.log",


rule depth_plots:
    input:
        postfilter_stats=expand(
            RESULTS / "QC/stats/postfilter/{mode}/{version}/{model}.csv",
            mode=MODES,
            version=VERSIONS,
            model=MODELS,
        ),
        ont_pr=expand(
            RESULTS
            / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.precision-recall.tsv",
            caller=CALLERS,
            depth=DEPTHS,
            mode=MODES,
            version=VERSIONS,
            model=MODELS,
            sample=SAMPLES,
        ),
        illumina_pr=expand(
            RESULTS / "assess/mutref/illumina/{sample}/{sample}.precision-recall.tsv",
            sample=SAMPLES,
        ),
    output:
        snp_fig=FIGURES / "depth_plots.snp.pdf",
        indel_fig=FIGURES / "depth_plots.indel.pdf",
    log:
        LOGS / "depth_plots.log",
    resources:
        runtime="10m",
        mem_mb=4 * GB,
    params:
        duplex_depth_cap=50,
        no_indels=config.get("no_indels", []),
    conda:
        ENVS / "depth_plots.yaml"
    script:
        SCRIPTS / "depth_plots.py"


rule precision_recall_curve:
    input:
        pr=expand(
            RESULTS
            / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{{model}}/{sample}/{sample}.precision-recall.tsv",
            caller=CALLERS,
            depth=[MAX_DEPTH],
            mode=MODES,
            version=VERSIONS,
            sample=SAMPLES,
        ),
        illumina_pr=expand(
            RESULTS / "assess/mutref/illumina/{sample}/{sample}.precision-recall.tsv",
            sample=SAMPLES,
        ),
    output:
        pdf=FIGURES / "precision_recall_curve.{model}.pdf",
    log:
        LOGS / "precision_recall_curve/{model}.log",
    resources:
        runtime="30m",
        mem_mb=8 * GB,
    conda:
        ENVS / "precision_recall_curve.yaml"
    params:
        no_indels=config.get("no_indels", []),
    script:
        SCRIPTS / "precision_recall_curve.py"


rule best_f1:
    input:
        pr=expand(
            RESULTS
            / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.precision-recall.tsv",
            caller=CALLERS,
            depth=[MAX_DEPTH],
            mode=MODES,
            version=VERSIONS,
            sample=SAMPLES,
            model=MODELS,
        ),
        illumina_pr=expand(
            RESULTS / "assess/mutref/illumina/{sample}/{sample}.precision-recall.tsv",
            sample=SAMPLES,
        ),
    output:
        figures=[
            FIGURES / f"best_f1_plots/{metric}.pdf"
            for metric in ["f1", "recall", "precision"]
        ],
    log:
        LOGS / "best_f1.log",
    resources:
        runtime="30m",
        mem_mb=8 * GB,
    conda:
        ENVS / "precision_recall_curve.yaml"
    params:
        no_indels=config.get("no_indels", []),
    script:
        SCRIPTS / "best_f1.py"


rule read_summary:
    input:
        csvs=expand(
            RESULTS / "QC/stats/prefilter/{mode}/{version}/{model}.csv",
            mode=MODES,
            version=VERSIONS,
            model=MODELS,
        ),
    output:
        csv=TABLES / "read_summary.csv",
        pdf=FIGURES / "read_identity.pdf",
    log:
        LOGS / "read_summary.log",
    resources:
        runtime="5m",
        mem_mb=500,
    conda:
        ENVS / "precision_recall_curve.yaml"
    script:
        SCRIPTS / "identity_plot.py"