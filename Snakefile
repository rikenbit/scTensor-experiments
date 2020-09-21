subworkflow synthetic_datasets:
    workdir:
        "synthetic_datasets"
    snakefile:
        "synthetic_datasets/Snakefile"
    configfile:
        "synthetic_datasets/config.yaml"

subworkflow real_datasets:
    workdir:
        "real_datasets"
    snakefile:
        "real_datasets/Snakefile"
    configfile:
        "real_datasets/config.yaml"

rule all:
	input:
		synthetic_datasets("data/input.csv"),
		real_datasets("data/input.csv")