REPEAT = config.get("repeat", 1)

for caller, caller_config in CALLERS.items():
    if "container" in caller_config:

        rule:
            name:
                f"{caller}_call_self"
            input:
                alignment=rules.align_to_self.output.alignment,
                reference=rules.align_to_self.input.reference,
                faidx=rules.faidx_reference.output.faidx,
            output:
                vcf=RESULTS
                / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.{{depth}}x.{caller}.vcf.gz",
            log:
                LOGS
                / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.log",
            benchmark:
                repeat(
                    BENCH
                    / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.tsv",
                    REPEAT,
                )
            threads: caller_config.get("threads", 4)
            resources:
                mem_mb=caller_config.get("memory", 16 * GB),
                runtime=caller_config.get("runtime", "1h"),
            container:
                caller_config["container"]
            shadow:
                "shallow"
            script:
                SCRIPTS / f"callers/{caller}.sh"

    elif "conda" in caller_config:

        rule:
            name:
                f"{caller}_call_self"
            input:
                alignment=rules.align_to_self.output.alignment,
                reference=rules.align_to_self.input.reference,
                faidx=rules.faidx_reference.output.faidx,
            output:
                vcf=RESULTS
                / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.{{depth}}x.{caller}.vcf.gz",
            log:
                LOGS
                / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.log",
            benchmark:
                repeat(
                    BENCH
                    / f"call/self/{caller}/{{depth}}x/{{mode}}/{{version}}/{{model}}/{{sample}}.tsv",
                    REPEAT,
                )
            threads: caller_config.get("threads", 4)
            resources:
                mem_mb=caller_config.get("memory", 16 * GB),
                runtime=caller_config.get("runtime", "1h"),
            conda:
                caller_config["conda"]
            shadow:
                "shallow"
            script:
                SCRIPTS / f"callers/{caller}.sh"

    else:
        raise KeyError(f"Caller {caller} has no container or conda environment")
