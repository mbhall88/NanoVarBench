# NanoVarBench

### Evaluating Nanopore-based bacterial variant calling

This repository holds the code for our paper (PREPRINT COMING SOON) which performs comprehensive benchmarking of SNP and indel variant calling accuracy across 14 diverse bacterial species using Oxford Nanopore Technologies (ONT) and Illumina sequencing.

You can find the results in that paper. Future updates after publication based on new tools, version, experiments etc. will be reported and shown here.

## Citation

Coming soon...

## Data

Accessions for all data can be found in [`config/accessions.csv`](./config/accessions.txt).

*(Note 15/03/2024) we are in the process of uploading the data to SRA so please be patient as this is taking some time. We will also be uploading the pod5 data to SRA too. BioProject and BioSample accessions have been generated and uploads for the appropriate samples will appear there.*

## Usage

See [the config docs](./config/README.md) for instructions on how to configure this pipeline for your data.

You will need the following packages to run the pipeline:
    - `snakemake`
    - `pandas`
    - `apptainer` or `singularity`
    - `conda`

A script for submitting the master Snakemake job on a Slurm cluster can be found at [`scripts/submit_slurm.sh`](./scripts/submit_slurm.sh).