# NanoVarBench

### Evaluating Nanopore-based bacterial variant calling

This repository holds the code for our paper (PREPRINT COMING SOON) which performs comprehensive benchmarking of SNP and indel variant calling accuracy across 14 diverse bacterial species using Oxford Nanopore Technologies (ONT) and Illumina sequencing.

You can find the results in that paper. Future updates after publication based on new tools, version, experiments etc. will be reported and shown here.

## Data

Accessions for all data can be found in [`config/accessions.csv`](./config/accessions.txt).

## Usage

See [the config docs](./config/README.md) for instructions on how to configure this pipeline for your data.

You will need the following packages to run the pipeline:
    - `snakemake`
    - `pandas`
    - `apptainer` or `singularity`
    - `conda`