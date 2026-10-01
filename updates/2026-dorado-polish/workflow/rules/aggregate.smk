# Aggregation into tidy tables (results contract in #6), plus tool versions.


def combos():
    return [
        dict(sample=s, read_model=rm, depth=d, arm=a)
        for s in RUN_SAMPLES
        for rm in RUN_READ_MODELS
        for d in RUN_DEPTHS
        for a in RUN_ARMS
    ]


def combo_inputs(c):
    keys = dict(c, depth=str(c["depth"]))
    score = lambda mode: str(SCORE).format(**keys, mode=mode)
    return {
        **keys,
        "sweep_summary": f"{score('sweep')}/precision-recall-summary.tsv",
        "pass_summary": f"{score('pass')}/precision-recall-summary.tsv",
        "sweep_pr": f"{score('sweep')}/precision-recall.tsv",
        "caller_info": str(rules.call_dorado.output.info).format(**keys),
        "coverage": str(rules.actual_depth.output.tsv).format(**keys),
    }


rule aggregate:
    input:
        files=[
            path
            for c in combos()
            for key, path in combo_inputs(c).items()
            if key not in ("sample", "read_model", "depth", "arm")
        ],
        samples=config["samples"],
    output:
        results=RESULTS / "tables/results.tsv",
        depth=RESULTS / "tables/depth.tsv",
        pr_curves=RESULTS / "tables/pr_curves.tsv",
    log:
        LOGS / "aggregate.log",
    resources:
        mem_mb=2000,
        runtime=10,
    params:
        combos=[combo_inputs(c) for c in combos()],
        arms=ARMS,
        read_models=config["read_models"],
    script:
        "../scripts/aggregate.py"


# Tools run from a container: the command that prints each one's version.
VERSION_COMMANDS = {
    "minimap2-2.31": "minimap2 --version; samtools --version | sed -n 1p",
    "samtools": "samtools --version | sed -n 1p",
    "bcftools": "bcftools --version | sed -n 1p",
    "seqkit": "seqkit version",
    "rasusa": "rasusa --version",
    "vcfdist": "vcfdist --version",
}


rule tool_version:
    output:
        txt=RESULTS / "versions/{tool}.txt",
    log:
        LOGS / "tool_version/{tool}.log",
    wildcard_constraints:
        tool="|".join(map(re.escape, VERSION_COMMANDS)),
    params:
        cmd=lambda wc: VERSION_COMMANDS[wc.tool],
        image=lambda wc: CONTAINERS[wc.tool],
    resources:
        mem_mb=500,
        runtime=5,
    container:
        lambda wc: CONTAINERS[wc.tool]
    shell:
        """
        exec 2> {log}
        {{ echo "image: {params.image}"; {params.cmd}; }} > {output.txt} 2>&1
        """


rule filter_env_version:
    output:
        txt=RESULTS / "versions/filter_env.txt",
    log:
        LOGS / "tool_version/filter_env.log",
    resources:
        mem_mb=500,
        runtime=5,
    conda:
        ENVS / "filter.yaml"
    shell:
        """
        exec 2> {log}
        {{ bcftools --version | sed -n 1p; python -c 'import cyvcf2, sys; print("cyvcf2", cyvcf2.__version__, "python", sys.version.split()[0])'; }} > {output.txt}
        """


def used_tools():
    tools = {"samtools", "bcftools", "seqkit", "rasusa", "vcfdist"}
    for spec in ALIGNMENTS.values():
        tools.add(f"{spec['aligner']}-{spec['version']}")
    return sorted(tools)


rule versions:
    input:
        tools=expand(RESULTS / "versions/{tool}.txt", tool=used_tools()),
        filter_env=rules.filter_env_version.output.txt,
    output:
        txt=RESULTS / "tables/versions.txt",
    resources:
        mem_mb=500,
        runtime=5,
    shell:
        """
        for f in {input}; do echo "== $(basename $f .txt)"; cat $f; done > {output.txt}
        """


localrules:
    aggregate,
    versions,
    tool_version,
    filter_env_version,
