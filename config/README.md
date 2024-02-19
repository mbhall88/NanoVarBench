# Configuring the pipeline

- [Configuring the pipeline](#configuring-the-pipeline)
  - [`config.yaml`](#configyaml)
    - [Variant caller configuration](#variant-caller-configuration)
      - [Example](#example)
  - [`pep/project_config.yaml`](#pepproject_configyaml)
    - [`sample_table`](#sample_table)
    - [`reads_dir`](#reads_dir)
    - [`reference_path`](#reference_path)
    - [`illumina_1` and `illumina_2`](#illumina_1-and-illumina_2)

This directory contains files for configuring the pipeline. Instructions relevant for each files are as follows:

## [`config.yaml`](./config.yaml)

The pipeline configuration file.

- `pepfile` is the path to the [PEP configuration file][pepconfig] (see [here](#pepproject_configyaml)).
- `pepschema` is the path to the [PEP schema][pepschema] file (see [here](../schemas/pep.yaml)).
- `model` is a list of all basecaller model versions you want to run the pipeline for. The version is a key (e.g., `v4.3.0:`) and the values are a list of the [full name for the models](https://github.com/nanoporetech/dorado#available-basecalling-models) you want to analyse.
- `QC` lists parameters for read quality control
  - `min_length` is the minimum length of a read to keep
  - `min_qscore` is the per read minimum mean quality score to keep
  - `depths` indicates to what depths the readset will be randomly subsampled to. By using multiple depths you can compare variant calling metrics across varying read depth. Only enter a single value if you don't want to do this, or set to 10000000 if you don't want subsampling.
- `mode` lists the sequencing modes you want to test - i.e., simplex, duplex
- `max_indel` the maximum indel length that will be allowed in the truth variants and assessed in the called variants.
- `truth` parameters related to the generation of the truth set of variants
  - `mash_distance` The [mash] distance between the reference and donor when generating the truth set. The closest distance will be chosen. Mash distance approximates the sequence divergence (i.e., 1-ANI).
  - `min_mash_distance` The minimum mash distance between the reference and donor when generating the truth set
  - `max_mash_distance` The maximum mash distance between the reference and donor when generating the truth set
  - `max_assemblies` The maximum number of assemblies to download for distance calculation
- `repeat` the number of times to repeat each variant calling rule. This is intended for benchmarking purposes and can be set to 1 for normal use.
- `callers` see the [callers section](#variant-caller-configuration).
- `no_indels` a list of `callers` that do not call indels

### Variant caller configuration

In the interest of making it easier to add/remove variant callers to the benchmark, there are some specific configuration requirements for variant caller addition.

Firstly, the `caller` section of the config is a list of variant callers you will be using. It is just a descriptive name for the caller, and doesn't have the be the actual command used to invoke the tool on the command line.

A Bash script that runs the variant caller must be provided. This is executed as a [Snakemake Bash script][smk-script]. See [the caller rules](../workflow/rules/call.smk) for the files and parameters available in the `snakemake` object in the scrpt. The script **must be named/located at `workflow/scripts/callers/<caller>.sh`**, where `<caller>` is the name of the caller provided in the config file. In addition, the script **must produce a compressed VCF file** (i.e., `.vcf.gz`). You must also provide either a container or conda enviornment file that the rule should be executed in. If using a conda environment, it must be located at `workflow/envs/<caller>.yaml`.

Once you have satisfied all of these requirements, copy and paste the `call_self` and `call_mutref` rules in `workflow/rules/call.smk` and rename them for your caller. The only thing you should need to change is the use of a `conda` or `container` directive.

#### Example

**Script**

```bash
#!/usr/bin/env bash
set -euo pipefail

exec 2>"${snakemake_log[0]}" # send all stderr from this script to the log file

aln="${snakemake_input[alignment]}"
ref="${snakemake_input[reference]}"
outvcf="${snakemake_output[vcf]}"

bcftools mpileup -f "$ref" \
    -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    -Ou --threads "${snakemake[threads]}" -Q 10 -x -M 10000 -h 100 "$aln" |
    bcftools call --ploidy 1 --threads "${snakemake[threads]}" -m --prior 0.005 \
        -o "$outvcf" --variants-only
```

**Conda environment**

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - bcftools=1.19
```

## [`pep/project_config.yaml`](./pep/project_config.yaml)

If you are unfamiliar with PEP (Portable Encapsulated Projects), then have a read through the documentation [here][pep].

### `sample_table`

This is the path to the sample table, relative to this `project_config.yaml` file. The specifications for a sample table can be read [here][pepsample]. For the purposes of this pipeline, here is a list of the columns and their meaning:

- `sample_name` which is a unique name for this sample (row).
- `species` the species of the sample (please use underscores instead of spaces)
- `taxid` the (NCBI) taxonomy ID for the species
- `reads_dir` see [here](#reads_dir) for more detail
- `reference_path` see [here](#reference_path) for more detail
- `illumina_1` see [here](#illumina_1-and-illumina_2) for more detail
- `illumina_2` see [here](#illumina_1-and-illumina_2) for more detail

### `reads_dir`

This column indicates the directory in which the reads exist for this sample. In our sample table, you will see values like `source1`. This is replaced during PEP valdation with the value indicated in the `project_config.yaml` section `sample_modifiers` > `derive` > `sources` > `source1`. See [this guide][peppathguide] for examples of how you can incorporate `sample_name` etc. into the path dynamically (if needed). If different samples occur in different locations, create a new source variable (e.g. `source6`) and add it to the list of `sources` as in our example.

There is an assumption about how this `reads_dir` is organised. The structure under this directory should follow the convention `{mode}/{model_version}/{model_name}/{sample_name}.fq.gz`, where `{mode}`, `{model_version}` and `{model_name}` are one of the sequencing modes, model versions, and names listed in [the pipeline config file](#configyaml) and `{sample_name}` is one of the sample names listed in the [sample table](#sample_table). Here is an example directory tree

```text
$ tree reads_dir/
└── simplex
    └── v4.3.0
        ├── dna_r10.4.1_e8.2_400bps_fast@v4.3.0
        │   ├── ATCC_10708__202309.fq.gz
        │   ├── ATCC_14035__202309.fq.gz
        │   ├── ATCC_17802__202309.fq.gz
        │   ├── ATCC_19119__202309.fq.gz
        │   ├── ATCC_25922__202309.fq.gz
        │   ├── ATCC_33560__202309.fq.gz
        │   ├── ATCC_35221__202309.fq.gz
        │   ├── ATCC_35897__202309.fq.gz
        │   └── ATCC_BAA-679__202309.fq.gz
        └── dna_r10.4.1_e8.2_400bps_sup@v4.3.0
            ├── ATCC_10708__202309.fq.gz
            ├── ATCC_14035__202309.fq.gz
            ├── ATCC_17802__202309.fq.gz
            ├── ATCC_19119__202309.fq.gz
            ├── ATCC_25922__202309.fq.gz
            ├── ATCC_33560__202309.fq.gz
            ├── ATCC_35221__202309.fq.gz
            ├── ATCC_35897__202309.fq.gz
            └── ATCC_BAA-679__202309.fq.gz
```

### `reference_path`

Similar to [`reads_dir`](#reads_dir), this column is a placeholder for the path to the reference genome of this sample. However, rather than being to a directory, this is a concrete path to a FASTA file. In our sample table, you will see values like `source2`. This is replaced during PEP valdation with the value indicated in the `project_config.yaml` section `sample_modifiers` > `derive` > `sources` > `source2`. An example from our config is `/data/ont/references/{sample_name}.fa`, where `{sample_name}` will be replaced with the sample name for that row, dynamically, on validation (when the pipeline is run). If you have references in various places that don't all follow the same naming convention then create a new source variable (e.g. `source8`) and add it to the list of `sources` as in our example.

### `illumina_1` and `illumina_2`

These columns are a placeholder for the path to the Illumina reads for each sample. In our sample table, you will see values like `source3` and `source4`. This is replaced during PEP valdation with the value indicated in the `project_config.yaml` section `sample_modifiers` > `derive` > `sources` > `source3`. An example from our config is `/data/illumina/{sample_name}_1.fq.gz`, where `{sample_name}` will be replaced with the sample name for that row, dynamically, on validation (when the pipeline is run). If you have files in various places that don't all follow the same naming convention then create a new source variable (e.g. `source8`) and add it to the list of `sources` as in our example.

[pepschema]: http://eido.databio.org/en/latest/writing-a-schema/
[pepconfig]: http://pep.databio.org/en/latest/specification/#project-config-file-specification
[pep]: http://pep.databio.org/en/latest/
[pepsample]: http://pep.databio.org/en/latest/specification/#sample-table-specification
[peppathguide]: http://pep.databio.org/en/latest/howto_eliminate_paths/
[mash]: https://github.com/marbl/Mash
[runtime]: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#snakefiles-standard-resources
[smk-script]: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#external-scripts