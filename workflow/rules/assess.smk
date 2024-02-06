rule assess_self_calls:
    input:
        vcf=RESULTS
        / "call/self/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.filter.vcf.gz",
    output:
        csv=RESULTS
        / "assess/self/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.csv",
    log:
        log=RESULTS
        / "assess_self_calls/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.log",
    resources:
        runtime="2m",
        mem_mb=500,
    conda:
        ENVS / "assess_self_calls.yaml"
    script:
        SCRIPTS / "assess_self_calls.py"


rule combine_self_calls:
    input:
        csvs=expand(
            RESULTS
            / "assess/self/{caller}/{depth}x/{mode}/{version}/{model}/{sample}.{depth}x.{caller}.csv",
            caller=CALLERS,
            depth=DEPTHS,
            mode=MODES,
            version=VERSIONS,
            model=MODELS,
            sample=SAMPLES,
        ),
    output:
        csv=RESULTS / "assess/self/self_calls.csv",
    log:
        log=RESULTS / "combine_self_calls.log",
    resources:
        runtime="2m",
        mem_mb=500,
    container:
        "docker://python:3.11-slim"
    script:
        SCRIPTS / "combine_csvs.py"
