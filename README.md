# scTensor-experiments

## About data preparation
- Synthetic datasets: synthetic_datasets/data/README.md
- Real datasets: real_datasets/data/README.md

## Requirements
- Snakemake

## How to perform all experiments

```
snakemake -j 4 --use-singularity # Local Machine
snakemake -j 32 --cluster qsub --latency-wait 20 --use-singularity # Open Grid Engine
snakemake -j 32 --cluster sbatch --latency-wait 20 --use-singularity # Slurm
```