# Scripts and Code for HyperGen Experiment

This repository contains the code and scripts to reproduce experimental results of HyperGen paper.


## Requirements

1. Put the benchmarked except for HyperGen in the `tools` folder
2. Install `HyperGen` from: [https://github.com/wh-xu/Hyper-Gen](https://github.com/wh-xu/Hyper-Gen)


## Pairwise ANI Analysis

1. `bash generate_pw_ani_dataset.sh` to create data for pairwise ANI analysis
2. `bash benchmark_pw_ani.sh` to benchmark different tools

## Database Search Analysis

1. `bash generate_search_dataset_link.sh` to create data for database search analysis
2. `bash benchmark_db_search.sh` to benchmark different tools


## BUSCO Analysis

1. Install `BUSCO` software from: [https://busco.ezlab.org](https://busco.ezlab.org)
2. Run `run_copy_fna_busco.py` to copy fna files in `BUSCO_fna_list.txt`

## Publication
1. Weihong Xu, Po-kai Hsu, Niema Moshiri, Shimeng Yu, and Tajana Rosing. "[HyperGen: Compact and Efficient Genome Sketching using Hyperdimensional Vectors](https://www.biorxiv.org/content/10.1101/2024.03.05.583605)." _Under review_.


## Contact
For more information, post an issue or send an email to <wexu@ucsd.edu>.