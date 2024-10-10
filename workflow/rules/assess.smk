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
    log:
        LOGS
        / "assess_mutref_calls/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.log",
    resources:
        runtime="20m",
        mem_mb=int(2 * GB),
    container:
        "docker://timd1/vcfdist:v2.3.3"
    params:
        opts=f"--largest-variant {config['truth']['max_indel']} --credit-threshold 1.0",
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
    log:
        LOGS / "assess_mutref_calls_illumina/{sample}.log",


rule annotate_dist_to_repeats_and_density:
    input:
        vcf=RESULTS
        / "assess/mutref/{caller}/100x/{mode}/{version}/{model}/{sample}/{sample}.summary.vcf",
        repeats=RESULTS / "truth/{sample}/{sample}.repetitive_regions.bed",
        reference=rules.create_mutref.output.mutref,
        script=SCRIPTS / "variant_density.py",
        hom_script=SCRIPTS / "annotate_homopolymers.py",
    output:
        vcf=RESULTS
        / "assess/mutref/{caller}/100x/{mode}/{version}/{model}/{sample}/{sample}.summary.annotated.vcf.gz",
    log:
        LOGS
        / "annotate_dist_to_repeats_and_density/mutref/{caller}/{mode}/{version}/{model}/{sample}.log",
    resources:
        mem_mb=int(0.5 * GB),
        runtime="3m",
    params:
        window=100,  # window size for density calculation
    conda:
        ENVS / "annotate.yaml"
    shell:
        """
        exec 2> {log}
        tmpinvcf=$(mktemp -u).vcf.gz
        >&2 echo "Annotating VCF with homopolymer length..."
        python {input.hom_script} {input.vcf} {input.reference} | bcftools view -o $tmpinvcf
        bcftools index -f $tmpinvcf

        >&2 echo "Annotating VCF with density..."
        den_vcf=$(mktemp -u).vcf.gz
        python {input.script} -w {params.window} $tmpinvcf | bcftools view -o $den_vcf
        bcftools index -f $den_vcf

        >&2 echo "Annotating VCF with distance to repetitive regions..."
        hdr=$(mktemp -u).hdr
        dist=$(mktemp -u).tab.gz
        echo '##INFO=<ID=DIST,Number=1,Type=Integer,Description="Distance to closest repetitive region">' > $hdr
        bedtools closest -a $den_vcf -b {input.repeats} -d -t first | 
            cut -f1,2,15 |
            bgzip -c > $dist
        tabix -f -s1 -b2 -e2 $dist
        bcftools annotate -a $dist -h $hdr -c CHROM,POS,DIST $den_vcf -o {output.vcf}
        bcftools index -f {output.vcf}

        rm -f $den_vcf $hdr $dist $tmpinvcf
        >&2 echo "Done."
        """


use rule annotate_dist_to_repeats_and_density as annotate_dist_to_repeats_and_density_illumina with:
    input:
        vcf=RESULTS / "assess/mutref/illumina/{sample}/{sample}.summary.vcf",
        repeats=RESULTS / "truth/{sample}/{sample}.repetitive_regions.bed",
        reference=rules.create_mutref.output.mutref,
        script=SCRIPTS / "variant_density.py",
        hom_script=SCRIPTS / "annotate_homopolymers.py",
    output:
        vcf=RESULTS
        / "assess/mutref/illumina/{sample}/{sample}.summary.annotated.vcf.gz",
    log:
        LOGS / "annotate_dist_to_repeats_and_density_illumina/mutref/{sample}.log",


rule combine_annotations:
    input:
        vcfs=expand(
            RESULTS
            / "assess/mutref/{caller}/100x/{mode}/{{version}}/{model}/{sample}/{sample}.summary.annotated.vcf.gz",
            caller=CALLERS,
            mode=MODES,
            model=MODELS,
            sample=SAMPLES,
        ),
        illumina_vcfs=expand(
            RESULTS
            / "assess/mutref/illumina/{sample}/{sample}.summary.annotated.vcf.gz",
            sample=SAMPLES,
        ),
    output:
        csv=RESULTS / "assess/mutref/annotations/{version}_annotations.csv",
    log:
        LOGS / "combine_annotations/{version}.log",
    resources:
        runtime="20m",
        mem_mb=2 * GB,
    container:
        "docker://quay.io/biocontainers/cyvcf2:0.30.28--py310hcf1fb4a_0"
    script:
        SCRIPTS / "combine_dist_and_density.py"


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
    log:
        LOGS
        / "assess_mutref_calls_without_repetitive_regions/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.log",
    resources:
        runtime="20m",
        mem_mb=int(2 * GB),
    container:
        "docker://timd1/vcfdist:v2.3.3"
    params:
        opts=rules.assess_mutref_calls.params.opts,
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
    log:
        LOGS / "assess_mutref_calls_illumina_without_repetitive_regions/{sample}.log",


rule depth_plots:
    input:
        postfilter_stats=expand(
            RESULTS / "QC/stats/postfilter/{mode}/{{version}}/{model}.csv",
            mode=MODES,
            model=MODELS,
        ),
        ont_pr=expand(
            RESULTS
            / "assess/mutref/{caller}/{depth}x/{mode}/{{version}}/{model}/{sample}/{sample}.precision-recall.tsv",
            caller=CALLERS,
            depth=DEPTHS,
            mode=MODES,
            model=MODELS,
            sample=SAMPLES,
        ),
        illumina_pr=expand(
            RESULTS / "assess/mutref/illumina/{sample}/{sample}.precision-recall.tsv",
            sample=SAMPLES,
        ),
    output:
        snp_fig=FIGURES / "depth_plots/{version}.snp.pdf",
        indel_fig=FIGURES / "depth_plots/{version}.indel.pdf",
    log:
        LOGS / "depth_plots/{version}.log",
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
        csv=TABLES / "best_f1.csv",
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


rule benchmark_resources:
    input:
        faidx=expand(RESULTS / "truth/{sample}/mutreference.fna.fai", sample=SAMPLES),
        call_benchmark=expand(
            BENCH
            / "call/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.tsv",
            caller=CALLERS,
            depth=DEPTHS,
            mode=MODES,
            version=VERSIONS,
            model=MODELS,
            sample=SAMPLES,
        ),
        align_benchmark=expand(
            BENCH / "align_to_mutref/{depth}x/{mode}/{version}/{model}/{sample}.tsv",
            sample=SAMPLES,
            depth=DEPTHS,
            mode=MODES,
            version=VERSIONS,
            model=MODELS,
        ),
    output:
        pdf=FIGURES / "benchmark_resources.pdf",
        csv=TABLES / "benchmark_resources.csv",
    log:
        LOGS / "benchmark_resources.log",
    resources:
        runtime="5m",
        mem_mb=500,
    conda:
        ENVS / "precision_recall_curve.yaml"
    script:
        SCRIPTS / "benchmark_plot.py"

rule plot_false_calls:
    input:
        table=rules.combine_annotations.output.csv,
        pr_wo_repeats=expand(
            RESULTS
            / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.without_repetitive_regions.precision-recall.tsv",
            caller=CALLERS,
            depth=[MAX_DEPTH],
            mode=MODES,
            version=VERSIONS,
            model=MODELS,
            sample=SAMPLES,
        ),
        pr_illumina_wo_repeats=expand(
            RESULTS
            / "assess/mutref/illumina/{sample}/{sample}.without_repetitive_regions.precision-recall.tsv",
            sample=SAMPLES,
        ),
        pr_w_repeats=expand(
            RESULTS
            / "assess/mutref/{caller}/{depth}x/{mode}/{version}/{model}/{sample}/{sample}.precision-recall.tsv",
            caller=CALLERS,
            depth=[MAX_DEPTH],
            mode=MODES,
            version=VERSIONS,
            model=MODELS,
            sample=SAMPLES,
        ),
        pr_illumina_w_repeats=expand(
            RESULTS / "assess/mutref/illumina/{sample}/{sample}.precision-recall.tsv",
            sample=SAMPLES,
        ),
    output:
        density_pdf=FIGURES / "false_calls/false_calls.density.pdf",
        fn_pdfs=expand(
            FIGURES / "false_calls/homopolymer.{caller}.fns.pdf",
            caller=[c for c in CALLERS if c not in NO_INDELS],
        ),
        fp_pdfs=expand(
            FIGURES / "false_calls/homopolymer.{caller}.fps.pdf",
            caller=[c for c in CALLERS if c not in NO_INDELS],
        ),
    log:
        LOGS / "plot_false_calls.log",
    resources:
        runtime="10m",
        mem_mb=8 * GB,
    conda:
        ENVS / "precision_recall_curve.yaml"
    script:
        SCRIPTS / "plot_false_calls.py"