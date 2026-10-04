# Figures and tables for the post (#14), made from the aggregated tables alone: no rule here
# reads a VCF, BAM or score file, so they can be rebuilt from results/tables/ without the
# work directory. Rendered outputs are not committed with the code (see README).
#
# When the AF filter analysis (#20) is enabled for the Arm in figures.af_filter, Figures 1-3
# and Table S1 add its series (#28), read from the analysis' own aggregated tables. Otherwise
# these inputs are empty, the scripts get no series, and the outputs are unchanged.

FIGURES = RESULTS / "figures"
TABLES = RESULTS / "tables"
FIGURE_CONFIG = config["figures"]

FIGURE_NAMES = {
    "fig1": "fig1_best_f1_depth",
    "fig2": "fig2_pr_curves",
    "fig3": "fig3_per_sample_best_f1",
}
TABLE1 = TABLES / "table1_runtime_memory.csv"
TABLE_S1 = TABLES / "table_s1_per_sample.csv"
# The AF filter series' inputs: its per-Sample scores and, for Figure 2, its PR curves.
AF_RESULTS = [rules.clair3_af_filter_tables.output.tsv] if FIGURE_AF else []
AF_CURVES = [rules.clair3_af_filter_pr_curves.output.tsv] if FIGURE_AF else []


def figure_files(name):
    return {ext: FIGURES / f"{FIGURE_NAMES[name]}.{ext}" for ext in ("png", "svg")}


def figure_targets():
    return [
        *(path for name in FIGURE_NAMES for path in figure_files(name).values()),
        TABLE1,
        TABLE_S1,
    ]


rule fig1_best_f1_depth:
    """Figure 1: Best F1 against Depth, with Arm C's and D's Default-PASS scores dashed."""
    input:
        results=rules.aggregate.output.results,
        af_results=AF_RESULTS,
    output:
        **figure_files("fig1"),
    log:
        LOGS / "fig1_best_f1_depth.log",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
    conda:
        ENVS / "plot.yaml"
    params:
        arm_labels=FIGURE_CONFIG["arm_labels"],
        default_pass_arms=FIGURE_CONFIG["default_pass_arms"],
        af_filter=FIGURE_AF,
        dpi=FIGURE_CONFIG["dpi"],
    script:
        "../scripts/fig1_best_f1_depth.py"


rule fig2_pr_curves:
    """Figure 2: precision-recall curves from the QUAL sweep, with the Default-PASS points."""
    input:
        results=rules.aggregate.output.results,
        pr_curves=rules.aggregate.output.pr_curves,
        af_results=AF_RESULTS,
        af_curves=AF_CURVES,
    output:
        **figure_files("fig2"),
    log:
        LOGS / "fig2_pr_curves.log",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10,
    conda:
        ENVS / "plot.yaml"
    params:
        arm_labels=FIGURE_CONFIG["arm_labels"],
        depths=FIGURE_CONFIG["pr_depths"],
        af_filter=FIGURE_AF,
        dpi=FIGURE_CONFIG["dpi"],
    script:
        "../scripts/fig2_pr_curves.py"


rule fig3_per_sample_best_f1:
    """Figure 3: every Sample's Best F1 at every Depth, with the dnd Samples highlighted."""
    input:
        results=rules.aggregate.output.results,
        af_results=AF_RESULTS,
    output:
        **figure_files("fig3"),
    log:
        LOGS / "fig3_per_sample_best_f1.log",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
    conda:
        ENVS / "plot.yaml"
    params:
        arm_labels=FIGURE_CONFIG["arm_labels"],
        af_filter=FIGURE_AF,
        dpi=FIGURE_CONFIG["dpi"],
    script:
        "../scripts/fig3_per_sample.py"


rule table1_runtime_memory:
    """Table 1: wall time and peak memory of each alignment and calling step, from benchmarks.tsv."""
    input:
        benchmarks=rules.benchmarks.output.tsv,
    output:
        csv=TABLE1,
    log:
        LOGS / "table1_runtime_memory.log",
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
    conda:
        ENVS / "plot.yaml"
    params:
        arm_labels=FIGURE_CONFIG["arm_labels"],
    script:
        "../scripts/table1_runtime.py"


rule table_s1_per_sample:
    """Table S1: per-Sample results for every Read set and Arm, with actual per-contig depth,
    as a CSV for the site's interactive csv-table."""
    input:
        results=rules.aggregate.output.results,
        depth=rules.aggregate.output.depth,
        af_results=AF_RESULTS,
    output:
        csv=TABLE_S1,
    log:
        LOGS / "table_s1_per_sample.log",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=5,
    conda:
        ENVS / "plot.yaml"
    params:
        arm_labels=FIGURE_CONFIG["arm_labels"],
        af_filter=FIGURE_AF,
    script:
        "../scripts/table_s1_per_sample.py"


localrules:
    fig1_best_f1_depth,
    fig2_pr_curves,
    fig3_per_sample_best_f1,
    table1_runtime_memory,
    table_s1_per_sample,
