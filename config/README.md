# Configuring the pipeline

- [Configuring the pipeline](#configuring-the-pipeline)
  - [`config.yaml`](#configyaml)
  - [`pep/project_config.yaml`](#pepproject_configyaml)
    - [`sample_table`](#sample_table)
    - [`reads_dir`](#reads_dir)
    - [`reference_path`](#reference_path)

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

## [`pep/project_config.yaml`](./pep/project_config.yaml)

If you are unfamiliar with PEP (Portable Encapsulated Projects), then have a read through the documentation [here][pep].

### `sample_table`

This is the path to the sample table, relative to this `project_config.yaml` file. The specifications for a sample table can be read [here][pepsample]. For the purposes of this pipeline, here is a list of the columns and their meaning:

- `sample_name` which is a unique name for this sample (row).
- `organism` the organism of the sample (please use underscores instead of spaces)
- `reads_dir` see [here](#reads_dir) for more detail
- `reference_path` see [here](#reference_path) for more detail

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

[pepschema]: http://eido.databio.org/en/latest/writing-a-schema/
[pepconfig]: http://pep.databio.org/en/latest/specification/#project-config-file-specification
[pep]: http://pep.databio.org/en/latest/
[pepsample]: http://pep.databio.org/en/latest/specification/#sample-table-specification
[peppathguide]: http://pep.databio.org/en/latest/howto_eliminate_paths/
