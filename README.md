# NanoVarBench

## Evaluating Nanopore-based bacterial variant calling

This repository holds the code for [our paper][doi] which performs comprehensive benchmarking of SNP and indel variant calling accuracy across 14 diverse bacterial species using Oxford Nanopore Technologies (ONT) and Illumina sequencing.

You can find the results in that paper. Future updates after publication based on new tools, versions, experiments etc. will be reported and shown here.

- [Citation](#citation)
- [Data](#data)
- [Usage](#usage)

## Citation

>  Benchmarking reveals superiority of deep learning variant callers on bacterial nanopore sequence data
Michael B. Hall, Ryan R. Wick, Louise M. Judd, An N. T. Nguyen, Eike J. Steinig, Ouli Xie, Mark R. Davies, Torsten Seemann, Timothy P. Stinear, Lachlan J. M. Coin
bioRxiv 2024.03.15.585313; doi: [10.1101/2024.03.15.585313][doi]

```bibtex
@article{hall_benchmarking_2024,
    title = {Benchmarking reveals superiority of deep learning variant callers on bacterial nanopore sequence data},
    url = {https://www.biorxiv.org/content/early/2024/03/16/2024.03.15.585313},
    doi = {10.1101/2024.03.15.585313},
    journal = {bioRxiv},
    author = {Hall, Michael B. and Wick, Ryan R. and Judd, Louise M. and Nguyen, An N. T. and Steinig, Eike J. and Xie, Ouli and Davies, Mark R. and Seemann, Torsten and Stinear, Timothy P. and Coin, Lachlan J. M.},
    year = {2024},
    pages = {2024.03.15.585313}
}
```

## Data

Accessions and DOIs for all data can be found in [`config/accessions.csv`](./config/accessions.csv).

The variant truthsets and associated data for making these is [archived on Zenodo][truth].

## Usage

See [the config docs](./config/README.md) for instructions on how to configure this pipeline for your data.

You will need the following packages to run the pipeline:

- `snakemake`
- `pandas`
- `apptainer` or `singularity`
- `conda`

A script for submitting the master Snakemake job on a Slurm cluster can be found at [`scripts/submit_slurm.sh`](./scripts/submit_slurm.sh).


[doi]: https://doi.org/10.1101/2024.03.15.585313 
[truth]: https://zenodo.org/doi/10.5281/zenodo.10867170