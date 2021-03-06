source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Data loading
load(infile)

# Label Permutation
out <- scTensor:::.celltypemergedtensor.hl(
    input=input,
    LR=LR,
    celltypes=celltypes,
    mergeas="mean",
    outerfunc="*")

# Save
save(out, file=outfile)
