# scTensor-experiments

# About data preparation
- Synthetic datasets: synthetic_datasets/data/README.md
- Real datasets: real_datasets/data/README.md

# Requirements
- Bash: GNU bash, version 4.2.46(1)-release (x86_64-redhat-linux-gnu)
- Snakemake: 5.3.0
- Singularity: 3.5.3

# How to reproduce this workflow

```
snakemake -j 4 --use-singularity # Local Machine
snakemake -j 32 --cluster qsub --latency-wait 600 --use-singularity # Open Grid Engine
snakemake -j 32 --cluster sbatch --latency-wait 600 --use-singularity # Slurm
```
# License
Copyright (c) 2020 Koki Tsuyuzaki and RIKEN Bioinformatics Research Unit Released under the [Artistic License 2.0](http://www.perlfoundation.org/artistic_license_2_0).

# Authors
- Koki Tsuyuzaki
- Manabu Ishii
- Itoshi Nikaido