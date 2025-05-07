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
_eLife_ (2024) 13:RP98300; doi: [10.7554/eLife.98300][doi]

```bibtex
@article{hall_benchmarking_2024,
	title = {Benchmarking reveals superiority of deep learning variant callers on bacterial nanopore sequence data},
	volume = {13},
	copyright = {Creative Commons Attribution-ShareAlike 4.0 International License (CC-BY-SA)},
	issn = {2050-084X},
	url = {https://doi.org/10.7554/eLife.98300},
	doi = {10.7554/eLife.98300},
	urldate = {2024-10-17},
	journal = {eLife},
	author = {Hall, Michael B and Wick, Ryan R and Judd, Louise M and Nguyen, An N and Steinig, Eike J and Xie, Ouli and Davies, Mark and Seemann, Torsten and Stinear, Timothy P and Coin, Lachlan},
	editor = {Weigel, Detlef},
	month = oct,
	year = {2024},
	pages = {RP98300},
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


[doi]: https://doi.org/10.7554/eLife.98300 
[truth]: https://zenodo.org/doi/10.5281/zenodo.10867170
