rule download_pod5: 
    output:
        pod5=temp(directory(RESULTS / "pod5/{sample}"))
    log:
        LOGS / "download_pod5/{sample}.log"
    resources:
        mem_mb=2_000,
        runtime="1d",
    params:
        opts="",
        url=infer_pod5_url
    shell:
        """
        mkdir -p {output.pod5} 2> {log}
        (wget {params.opts} -O - {params.url} | tar x --directory={output.pod5} --strip-components=1) 2>> {log}
        """


rule download_dorado_model:
    output:
        model=directory(RESULTS / "dorado_models/{model}@{version}")
    
    log:
        LOGS / "download_dorado_model/{model}/{version}.log"
    resources:
        mem_mb=1_000,
        runtime="5m",
    params:
        model_dir=lambda wildcards, output: Path(output.model).parent,
        opts="--overwrite",
        model="{model}@{version}"
    shell:
        "dorado download --model {params.model} --models-directory {params.model_dir} {params.opts} 2> {log}"

rule basecall:
    input:
        pod5=rules.download_pod5.output.pod5,
        model=rules.download_dorado_model.output.model
    output:
        reads=RESULTS / "basecall/{mode}/{version}/{model}/{sample}.fq.gz",
    log:
        LOGS / "basecall/{mode}/{version}/{model}/{sample}.log",
    benchmark:
        BENCH / "basecall/{mode}/{version}/{model}/{sample}.tsv"
    resources:
        mem_mb=20 * GB,
        runtime="12h",
        partition="gpu-a100",
        slurm="gres=gpu:1",
    threads: 2
    params:
        opts="-r",
        model_dir=lambda wildcards, input: Path(input.model).parent,
    conda:
        ENVS / "align.yaml"
    shell:
        """
        if [ {wildcards.mode} = "duplex" ]; then
            subcommand="duplex"
            tag_opt=(-d "dx:1")
        else
            subcommand="basecaller"
            tag_opt=()
        fi

        (dorado $subcommand {params.opts} {input.model} {input.pod5} | \
            samtools fastq -T '*' "${{tag_opt[@]}}" | gzip) 2> {log} > {output.reads}
        """